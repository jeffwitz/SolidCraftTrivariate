//! Build exact tensor-product NURBS covers of the tall disk/cylinder model.
//!
//! Usage: cargo run --release --example cylinder_tensor_covers -- MODEL OUTDIR

use solidcraft_trivariate::json::{object, parse, Json};
use std::{collections::BTreeMap, env, fs};

type H = [i64; 4];
type Grid = Vec<Vec<Vec<Option<H>>>>;
type ControlNet = Vec<Vec<Vec<H>>>;
type SourceBlock = (Vec<usize>, ControlNet);
type SourceBlocks = BTreeMap<String, SourceBlock>;

fn n(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}
fn s(value: impl Into<String>) -> Json {
    Json::String(value.into())
}
fn b(value: bool) -> Json {
    Json::Bool(value)
}
fn arr(values: impl IntoIterator<Item = Json>) -> Json {
    Json::Array(values.into_iter().collect())
}
fn nums(values: &[usize]) -> Json {
    arr(values.iter().map(|x| n(*x)))
}
fn bools(values: &[bool]) -> Json {
    arr(values.iter().copied().map(b))
}
fn jstrs(values: &[String]) -> Json {
    arr(values.iter().cloned().map(s))
}
fn int_array(value: &Json) -> Result<Vec<i64>, String> {
    value
        .items()
        .map_err(str::to_owned)?
        .iter()
        .map(|x| x.number().map_err(str::to_owned))
        .collect()
}

fn read_grid(block: &Json) -> Result<(String, Vec<usize>, ControlNet), String> {
    let id = block
        .field("id")
        .map_err(str::to_owned)?
        .text()
        .map_err(str::to_owned)?
        .to_owned();
    let degrees: Vec<usize> = block
        .field("degrees")
        .map_err(str::to_owned)?
        .items()
        .map_err(str::to_owned)?
        .iter()
        .map(|x| x.number().map_err(str::to_owned))
        .collect::<Result<_, _>>()?;
    let controls = block
        .field("controls")
        .map_err(str::to_owned)?
        .items()
        .map_err(str::to_owned)?;
    if degrees.len() != 3
        || controls.len() != (degrees[0] + 1) * (degrees[1] + 1) * (degrees[2] + 1)
    {
        return Err(format!("{id}: invalid control dimensions"));
    }
    let mut grid = vec![vec![vec![[0; 4]; degrees[2] + 1]; degrees[1] + 1]; degrees[0] + 1];
    for i in 0..=degrees[0] {
        for j in 0..=degrees[1] {
            for k in 0..=degrees[2] {
                let p = int_array(&controls[(i * (degrees[1] + 1) + j) * (degrees[2] + 1) + k])?;
                if p.len() != 4 || p[3] <= 0 {
                    return Err(format!("{id}: invalid homogeneous control"));
                }
                grid[i][j][k] = [p[0], p[1], p[2], p[3]];
            }
        }
    }
    Ok((id, degrees, grid))
}

fn put(grid: &mut Grid, i: usize, j: usize, k: usize, value: H, label: &str) -> Result<(), String> {
    if let Some(old) = grid[i][j][k] {
        if old != value {
            return Err(format!("{label}: incompatible shared control {i},{j},{k}"));
        }
    }
    grid[i][j][k] = Some(value);
    Ok(())
}

fn put_point(slot: &mut Option<H>, value: H, label: &str) -> Result<(), String> {
    if let Some(old) = slot {
        if *old != value {
            return Err(format!("{label}: incompatible shared control"));
        }
    }
    *slot = Some(value);
    Ok(())
}

fn verify_source_ring(blocks: &SourceBlocks) -> Result<(), String> {
    for level in ["lower", "middle", "upper"] {
        for (a, au, c, cu) in [
            ("east", 0, "south", 2),
            ("east", 2, "north", 0),
            ("north", 2, "west", 0),
            ("south", 0, "west", 2),
        ] {
            let a_name = format!("disk_{level}_sector_{a}");
            let c_name = format!("disk_{level}_sector_{c}");
            let a_grid = &blocks[&a_name].1;
            let c_grid = &blocks[&c_name].1;
            for v in 0..2 {
                for w in 0..2 {
                    if a_grid[au][v][w] != c_grid[cu][v][w] {
                        return Err(format!("source ring seam mismatch: {a_name} to {c_name}"));
                    }
                }
            }
        }
    }
    Ok(())
}

fn knots(degree: usize, spans: usize) -> Json {
    let mut values = vec![0; degree + 1];
    for k in 1..spans {
        values.extend(std::iter::repeat_n(k, degree));
    }
    values.extend(std::iter::repeat_n(spans, degree + 1));
    arr(values.into_iter().map(n))
}

#[allow(clippy::needless_range_loop)] // Tensor indexing must follow the declared control-net axis order.
fn patch(
    id: &str,
    cells: &[String],
    degrees: [usize; 3],
    spans: [usize; 3],
    grid: Grid,
    span_info: Vec<Json>,
) -> Json {
    let shape = [
        spans[0] * degrees[0] + 1,
        spans[1] * degrees[1] + 1,
        spans[2] * degrees[2] + 1,
    ];
    let mut controls = Vec::new();
    for i in 0..shape[0] {
        for j in 0..shape[1] {
            for k in 0..shape[2] {
                let p = grid[i][j][k].expect("complete tensor net");
                controls.push(arr(p.into_iter().map(n)));
            }
        }
    }
    object([
        ("id", s(id)),
        ("patch_kind", s("multi_span_nurbs")),
        ("degrees", nums(&degrees)),
        ("knots", arr((0..3).map(|a| knots(degrees[a], spans[a])))),
        ("control_shape", nums(&shape)),
        ("controls_homogeneous", arr(controls)),
        ("source_blocks", jstrs(cells)),
        ("cells", jstrs(cells)),
        ("exact_piece_extraction_expected", b(true)),
        ("homogeneous_common_scale_denominator", n(1)),
        ("spans", arr(span_info)),
    ])
}

fn span(source: &str, intervals: [[usize; 2]; 3], angle: [usize; 2]) -> Json {
    object([
        ("source_block", s(source)),
        ("axis_map", nums(&[0, 1, 2])),
        ("flips", bools(&[false, false, false])),
        ("parameter_interval", arr([n(angle[0]), n(angle[1])])),
        (
            "parameter_intervals",
            arr(intervals.into_iter().map(|r| arr([n(r[0]), n(r[1])]))),
        ),
        ("degree_elevation", nums(&[0, 0, 0])),
    ])
}

fn source_span(block_id: &str, angle: [usize; 2], z: usize) -> Json {
    span(block_id, [[angle[0], angle[1]], [0, 1], [z, z + 1]], angle)
}

#[allow(clippy::needless_range_loop)] // Tensor indexing must follow source control-net axis order.
fn main_result(model: &Json, blocks: &SourceBlocks, sectors: usize) -> Result<Json, String> {
    verify_source_ring(blocks)?;
    for level in ["lower", "middle", "upper"] {
        for name in [
            "core",
            "sector_east",
            "sector_north",
            "sector_west",
            "sector_south",
        ] {
            let id = format!("disk_{level}_{name}");
            if !blocks.contains_key(&id) {
                return Err(format!("missing expected source block {id}"));
            }
        }
    }
    let ring_groups: Vec<Vec<&str>> = if sectors == 2 {
        vec![vec!["east", "north", "west", "south"]]
    } else {
        vec![vec!["east", "north"], vec!["west", "south"]]
    };
    let mut patches = Vec::new();
    // The square disk core spans the three axial source blocks.
    let mut core_grid = vec![vec![vec![None; 4]; 3]; 3];
    let mut core_cells = Vec::new();
    let mut core_spans = Vec::new();
    for (z, level) in ["lower", "middle", "upper"].iter().enumerate() {
        let name = format!("disk_{level}_core");
        let (_, g) = &blocks[&name];
        core_cells.push(name.clone());
        core_spans.push(source_span(&name, [0, 1], z));
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..2 {
                    put(&mut core_grid, i, j, z + k, g[i][j][k], &name)?;
                }
            }
        }
    }
    patches.push(patch(
        "disk_core",
        &core_cells,
        [2, 2, 1],
        [1, 1, 3],
        core_grid,
        core_spans,
    ));
    for (group_i, group) in ring_groups.iter().enumerate() {
        let id = if sectors == 2 {
            "disk_ring".to_owned()
        } else {
            format!("disk_ring_{}", group_i + 1)
        };
        let angular_spans = group.len();
        let mut grid = vec![vec![vec![None; 4]; 2]; 2 * angular_spans + 1];
        let mut cells = Vec::new();
        let mut span_info = Vec::new();
        for (local_angle, sector) in group.iter().enumerate() {
            for (z, level) in ["lower", "middle", "upper"].iter().enumerate() {
                let name = format!("disk_{level}_sector_{sector}");
                let (_, g) = &blocks[&name];
                cells.push(name.clone());
                span_info.push(source_span(&name, [local_angle, local_angle + 1], z));
                for u in 0..3 {
                    for v in 0..2 {
                        for w in 0..2 {
                            put_point(&mut grid[2 * local_angle + u][v][z + w], g[u][v][w], &name)?;
                        }
                    }
                }
            }
        }
        patches.push(patch(
            &id,
            &cells,
            [2, 1, 1],
            [angular_spans, 1, 3],
            grid,
            span_info,
        ));
    }
    let source_interfaces = model
        .field("topology")
        .map_err(str::to_owned)?
        .field("interfaces")
        .map_err(str::to_owned)?
        .items()
        .map_err(str::to_owned)?;
    let mut owners = BTreeMap::new();
    for p in &patches {
        for c in p
            .field("cells")
            .map_err(str::to_owned)?
            .items()
            .map_err(str::to_owned)?
        {
            owners.insert(
                c.text().unwrap().to_owned(),
                p.field("id").unwrap().text().unwrap().to_owned(),
            );
        }
    }
    let mut interfaces = Vec::new();
    let mut self_seams = Vec::new();
    let mut exact = 0;
    for iface in source_interfaces {
        exact += 1;
        let first = iface
            .field("first")
            .map_err(str::to_owned)?
            .items()
            .map_err(str::to_owned)?;
        let second = iface
            .field("second")
            .map_err(str::to_owned)?
            .items()
            .map_err(str::to_owned)?;
        let a = first[0].text().map_err(str::to_owned)?;
        let sa = first[1].text().map_err(str::to_owned)?;
        let c = second[0].text().map_err(str::to_owned)?;
        let sc = second[1].text().map_err(str::to_owned)?;
        let pa = &owners[a];
        let pc = &owners[c];
        let ma = mapped_side(a, sa, pa, sectors)?;
        let mb = mapped_side(c, sc, pc, sectors)?;
        let record = object([
            ("first", ma),
            ("second", mb),
            (
                "orientation",
                iface.field("orientation").map_err(str::to_owned)?.clone(),
            ),
            ("status", s("exact_match")),
        ]);
        if pa == pc {
            if sectors == 2
                && pa == "disk_ring"
                && ((a.ends_with("sector_east")
                    && sa == "u0"
                    && c.ends_with("sector_south")
                    && sc == "u1")
                    || (c.ends_with("sector_east")
                        && sc == "u0"
                        && a.ends_with("sector_south")
                        && sa == "u1"))
            {
                self_seams.push(record);
            }
        } else {
            interfaces.push(record);
        }
    }
    let exterior = model
        .field("topology")
        .map_err(str::to_owned)?
        .field("exterior")
        .map_err(str::to_owned)?
        .items()
        .map_err(str::to_owned)?;
    let mut boundary_faces = Vec::new();
    for face in exterior {
        let f = face.items().map_err(str::to_owned)?;
        let id = f[0].text().map_err(str::to_owned)?;
        let side = f[1].text().map_err(str::to_owned)?;
        let intervals = block_intervals(id, &owners[id], sectors)?;
        boundary_faces.push(object([
            ("patch_id", s(&owners[id])),
            ("source_block", s(id)),
            ("side", s(side)),
            ("patch_side", s(side)),
            ("axis_map", nums(&[0, 1, 2])),
            ("flips", bools(&[false, false, false])),
            ("parameter_intervals", interval_json(intervals)),
            ("status", s("extracted_exactly_from_model")),
        ]));
    }
    let model_sha = model
        .field("__sha256")
        .map(|x| x.text().unwrap().to_owned())
        .unwrap_or_default();
    let recipe_sha = model
        .field("__recipe_sha256")
        .map(|x| x.text().unwrap().to_owned())
        .unwrap_or_default();
    let all_cells = blocks.keys().cloned().collect::<Vec<_>>();
    let patch_count = patches.len();
    Ok(object([
        ("schema",s("gallery_exact_patch_map_v1")),("case_id",s(format!("box-cylinder-tall-cylinder-{sectors}"))),
        ("representation_kind",s("exact_disk_cylinder_tensor_nurbs_cover")),("method",s("exact_homogeneous_tensor_product_bezier_spans")),
        ("scope_note",s("Exact homogeneous tensor-product patch assembly from every source block. Angular and axial spans retain their source-cell parameter intervals; the periodic seam is recorded explicitly. This verifies representation and traces, not global injectivity or material coverage.")),
        ("source_model",s("cylinder.model.json")),("source_model_sha256",s(model_sha)),("source_recipe",s("cylinder.recipe.json")),("source_recipe_sha256",s(recipe_sha)),
        ("source_schema",model.field("schema").map_err(str::to_owned)?.clone()),("source_block_count",n(blocks.len())),("native_nurbs_patch_count",n(blocks.len())),
        ("patch_count",n(patch_count)),("exact_extraction_count",n(blocks.len())),("exact_interface_count",n(exact)),("source_interface_count",n(exact)),
        ("independent_sample_evaluations",n(0)),("axial_chain_patch_count",n(0)),("axial_chain_source_block_count",n(0)),("additional_exact_pair_patch_count",n(0)),
        ("patches",arr(patches)),("interfaces",arr(interfaces)),("boundary_faces",arr(boundary_faces)),("self_seams",arr(self_seams)),
        ("merged_source_blocks",jstrs(&all_cells)),("errors",arr(Vec::new())),
    ]))
}

fn mapped_side(block: &str, side: &str, patch: &str, sectors: usize) -> Result<Json, String> {
    let intervals = block_intervals(block, patch, sectors)?;
    Ok(object([
        ("patch_id", s(patch)),
        ("source_block", s(block)),
        ("side", s(side)),
        ("patch_side", s(side)),
        ("axis_map", nums(&[0, 1, 2])),
        ("flips", bools(&[false, false, false])),
        (
            "parameter_intervals",
            arr(intervals.into_iter().map(|r| arr([n(r[0]), n(r[1])]))),
        ),
    ]))
}

fn block_intervals(block: &str, patch: &str, sectors: usize) -> Result<[[usize; 2]; 3], String> {
    let z = if block.contains("lower") {
        0
    } else if block.contains("middle") {
        1
    } else if block.contains("upper") {
        2
    } else {
        return Err(format!("unknown source layer {block}"));
    };
    if block.ends_with("_core") {
        return Ok([[0, 1], [0, 1], [z, z + 1]]);
    }
    let a = ["east", "north", "west", "south"]
        .iter()
        .position(|x| block.ends_with(x))
        .ok_or_else(|| format!("unknown source sector {block}"))?;
    let local = if sectors == 2 || patch == "disk_ring_1" {
        a
    } else if patch == "disk_ring_2" {
        a.checked_sub(2).ok_or("sector assigned to wrong patch")?
    } else {
        return Err("unknown ring patch".into());
    };
    Ok([[local, local + 1], [0, 1], [z, z + 1]])
}

fn interval_json(intervals: [[usize; 2]; 3]) -> Json {
    arr(intervals.into_iter().map(|r| arr([n(r[0]), n(r[1])])))
}

fn sha256(path: &str) -> Result<String, String> {
    let output = std::process::Command::new("sha256sum")
        .arg(path)
        .output()
        .map_err(|e| format!("sha256sum: {e}"))?;
    if !output.status.success() {
        return Err("sha256sum failed".into());
    }
    String::from_utf8(output.stdout)
        .map_err(|e| e.to_string())
        .map(|x| x.split_whitespace().next().unwrap().to_owned())
}

fn run() -> Result<(), String> {
    let args: Vec<_> = env::args().collect();
    if args.len() != 3 {
        return Err("usage: cylinder_tensor_covers MODEL OUTDIR".into());
    }
    let model_text = fs::read_to_string(&args[1]).map_err(|e| e.to_string())?;
    let mut model = parse(&model_text).map_err(str::to_owned)?;
    validate_source(&model)?;
    solidcraft_trivariate::interchange::audit(&model)
        .map_err(|e| format!("source model audit: {e}"))?;
    let recipe_path = std::path::Path::new(&args[1]).with_file_name("cylinder.recipe.json");
    if let Json::Object(fields) = &mut model {
        fields.insert("__sha256".into(), s(sha256(&args[1])?));
        fields.insert(
            "__recipe_sha256".into(),
            s(sha256(recipe_path.to_str().ok_or("invalid recipe path")?)?),
        );
    }
    let mut blocks = BTreeMap::new();
    for block in model
        .field("blocks")
        .map_err(str::to_owned)?
        .items()
        .map_err(str::to_owned)?
    {
        let (id, d, g) = read_grid(block)?;
        blocks.insert(id, (d, g));
    }
    if blocks.len() != 15 {
        return Err(format!("expected 15 source blocks, got {}", blocks.len()));
    }
    for count in [2, 3] {
        let result = main_result(&model, &blocks, count)?;
        let file = format!("{}/cylinder-{count}.patches.json", args[2]);
        fs::write(&file, format!("{result}\n")).map_err(|e| e.to_string())?;
        println!("wrote {file}");
    }
    Ok(())
}

fn validate_source(model: &Json) -> Result<(), String> {
    if model
        .field("schema")
        .map_err(str::to_owned)?
        .text()
        .map_err(str::to_owned)?
        != "box_cylinder_boolean_suite_v1"
    {
        return Err("source model must use box_cylinder_boolean_suite_v1".into());
    }
    if model
        .field("operation")
        .map_err(str::to_owned)?
        .text()
        .map_err(str::to_owned)?
        != "cylinder"
    {
        return Err("source model operation must be cylinder".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn source() -> Json {
        parse(include_str!(
            "../docs/site/box-cylinder/tall/cylinder.model.json"
        ))
        .unwrap()
    }

    #[test]
    fn exact_covers_have_expected_topology_and_periodic_seam() {
        let model = source();
        validate_source(&model).unwrap();
        solidcraft_trivariate::interchange::audit(&model).unwrap();
        let mut blocks = BTreeMap::new();
        for block in model.field("blocks").unwrap().items().unwrap() {
            let (id, degrees, grid) = read_grid(block).unwrap();
            blocks.insert(id, (degrees, grid));
        }
        let two = main_result(&model, &blocks, 2).unwrap();
        let three = main_result(&model, &blocks, 3).unwrap();
        assert_eq!(two.field("patch_count").unwrap().number::<usize>(), Ok(2));
        assert_eq!(three.field("patch_count").unwrap().number::<usize>(), Ok(3));
        assert_eq!(two.field("self_seams").unwrap().items().unwrap().len(), 3);
        assert_eq!(three.field("self_seams").unwrap().items().unwrap().len(), 0);
    }

    #[test]
    fn rejects_non_cylinder_operations() {
        let mut model = source();
        let Json::Object(fields) = &mut model else {
            panic!("expected model object")
        };
        fields.insert("operation".into(), s("box"));
        assert_eq!(
            validate_source(&model),
            Err("source model operation must be cylinder".into())
        );
    }
}
fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        std::process::exit(1)
    }
}
