//! Deterministic X2 analysis-package export for dyadically refined brackets.

use crate::{
    adaptive_complex,
    adaptive_partition::{
        AdaptivePartition, ExteriorPatch, FacePatch, FaceRect, LeafKey, LeafSide,
    },
    adaptive_space::{self, AdaptiveSpace},
    analysis_space::{self, AnalysisSpace},
    boundary_traction::{integrate_traction, LoadPatch},
    bracket_audit::{self, BracketCertificate},
    bracket_json::decode_blocks,
    bracket_parameters_audit,
    complex_skeleton::encode_skeleton,
    constraint_system::{resolve, Lift, Ratio},
    incidence::BlockSide,
    json::{object, Json},
    solver_package,
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RefinementPlan {
    pub max_depth: u8,
    pub max_cells: usize,
    pub steps: Vec<Vec<LeafKey>>,
}

fn text(s: &str) -> Json {
    Json::String(s.to_owned())
}
fn number(n: impl ToString) -> Json {
    Json::Number(n.to_string())
}
fn numbers<T: ToString>(values: impl IntoIterator<Item = T>) -> Json {
    Json::Array(values.into_iter().map(number).collect())
}
fn arr<T>(values: impl IntoIterator<Item = T>, mut f: impl FnMut(T) -> Json) -> Json {
    Json::Array(values.into_iter().map(&mut f).collect())
}

/// Parse the strict, versioned list of successive active-leaf refinements.
pub fn parse_refinement_plan(plan: &Json) -> Result<RefinementPlan, &'static str> {
    plan.keys(&["schema", "max_depth", "max_cells", "steps"])?;
    if plan.field("schema")?.text()? != "adaptive_refinement_plan_v1" {
        return Err("unsupported adaptive refinement plan schema");
    }
    let max_depth = plan.field("max_depth")?.number::<u8>()?;
    let max_cells = plan.field("max_cells")?.number::<usize>()?;
    if max_depth > 20 || max_cells == 0 {
        return Err("invalid adaptive refinement budget");
    }
    let mut steps = Vec::new();
    for step in plan.field("steps")?.items()? {
        step.keys(&["refine"])?;
        let mut leaves = Vec::new();
        let mut unique = BTreeSet::new();
        for leaf in step.field("refine")?.items()? {
            leaf.keys(&["root", "level", "index"])?;
            let root = leaf.field("root")?.text()?.to_owned();
            let level = leaf.field("level")?.number::<u8>()?;
            let index: [u32; 3] = leaf
                .field("index")?
                .items()?
                .iter()
                .map(Json::number)
                .collect::<Result<Vec<_>, _>>()?
                .try_into()
                .map_err(|_| "expected three dyadic indices")?;
            if root.is_empty() || level > max_depth || index.iter().any(|&i| i >= (1_u32 << level))
            {
                return Err("invalid adaptive leaf request");
            }
            let key = LeafKey { root, level, index };
            if !unique.insert(key.clone()) {
                return Err("duplicate leaf in one refinement step");
            }
            leaves.push(key);
        }
        if leaves.is_empty() {
            return Err("empty adaptive refinement step");
        }
        leaves.sort();
        steps.push(leaves);
    }
    Ok(RefinementPlan {
        max_depth,
        max_cells,
        steps,
    })
}

fn side_name(side: &BlockSide) -> String {
    format!("{}{}", ["u", "v", "w"][side.axis], u8::from(side.upper))
}
fn leaf_json(leaf: &LeafKey) -> Json {
    object([
        ("root", text(&leaf.root)),
        ("level", number(leaf.level)),
        ("index", numbers(leaf.index)),
    ])
}
fn rect_json(rect: FaceRect) -> Json {
    arr(rect.ranges, numbers)
}
fn side_ref_json(adaptive: &AdaptiveSpace, side: &LeafSide) -> Result<Json, &'static str> {
    Ok(object([
        (
            "block",
            text(
                adaptive
                    .leaf_blocks
                    .get(&side.leaf)
                    .ok_or("leaf block is missing")?,
            ),
        ),
        ("side", text(&side_name(&side.root_side))),
        ("root_block", text(&side.root_side.block)),
        ("root_side", text(&side_name(&side.root_side))),
        ("leaf", leaf_json(&side.leaf)),
    ]))
}
fn face_json(
    adaptive: &AdaptiveSpace,
    patch: &FacePatch,
    denominator: u32,
) -> Result<Json, &'static str> {
    let o = patch.orientation;
    Ok(object([
        ("first", side_ref_json(adaptive, &patch.first)?),
        ("second", side_ref_json(adaptive, &patch.second)?),
        ("first_level", number(patch.first.leaf.level)),
        ("second_level", number(patch.second.leaf.level)),
        (
            "level_difference_first_minus_second",
            number(i16::from(patch.first.leaf.level) - i16::from(patch.second.leaf.level)),
        ),
        (
            "two_to_one_balanced",
            Json::Bool(patch.first.leaf.level.abs_diff(patch.second.leaf.level) <= 1),
        ),
        ("master_rect", rect_json(patch.master_rect)),
        ("first_subrect", rect_json(patch.first_subrect)),
        ("second_subrect", rect_json(patch.second_subrect)),
        ("grid_denominator", number(denominator)),
        (
            "orientation",
            object([
                ("swap", Json::Bool(o.swap)),
                ("flip_first", Json::Bool(o.flip_first)),
                ("flip_second", Json::Bool(o.flip_second)),
                ("maps", text("second_trace_into_first")),
            ]),
        ),
    ]))
}
fn exterior_json(
    adaptive: &AdaptiveSpace,
    patch: &ExteriorPatch,
    denominator: u32,
    feature: Option<&str>,
) -> Result<Json, &'static str> {
    let mut fields = match side_ref_json(adaptive, &patch.side)? {
        Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    fields.insert("root_rect".into(), rect_json(patch.root_rect));
    fields.insert("grid_denominator".into(), number(denominator));
    if let Some(tag) = feature {
        fields.insert("feature".into(), text(tag));
    }
    Ok(Json::Object(fields))
}
pub(crate) fn lift_json(lift: &Lift) -> Json {
    object([
        ("free_roots", numbers(&lift.free_roots)),
        (
            "master_index_convention",
            text("zero_based_position_in_free_roots"),
        ),
        (
            "lines",
            arr(&lift.lines, |line| {
                object([
                    ("dof", number(line.dof)),
                    (
                        "terms",
                        arr(&line.terms, |(master, ratio)| {
                            numbers([
                                master.to_string(),
                                ratio.numerator().to_string(),
                                ratio.denominator().to_string(),
                            ])
                        }),
                    ),
                    (
                        "offset",
                        numbers([line.offset.numerator(), line.offset.denominator()]),
                    ),
                ])
            }),
        ),
    ])
}
fn boundary_faces(
    adaptive: &AdaptiveSpace,
    partition: &AdaptivePartition,
    tags: &BTreeMap<BlockSide, String>,
    wanted: &str,
) -> Result<Vec<Json>, &'static str> {
    let denom = partition.grid_denominator();
    partition
        .exterior_patches()
        .iter()
        .filter(|p| tags.get(&p.side.root_side).is_some_and(|f| f == wanted))
        .map(|p| exterior_json(adaptive, p, denom, Some(wanted)))
        .collect()
}
fn adaptive_reduction(
    space: &AnalysisSpace,
    adaptive: &AdaptiveSpace,
) -> Result<(Json, Json, Json), &'static str> {
    let y50 = space
        .boundary_dofs
        .get("y50")
        .ok_or("adaptive space has no y50 selection")?;
    if y50.is_empty() {
        return Err("empty adaptive y50 selection");
    }
    let prescribed = y50.iter().map(|&d| (d, Ratio::zero())).collect::<Vec<_>>();
    let mounted = resolve(space.ndof, &adaptive.constraints, &prescribed)?;
    let nfull = 3 * space.ndof;
    let nfree = mounted.free_roots.len();
    let mut fixed = Vec::with_capacity(3 * y50.len());
    for &d in y50 {
        for c in 0..3 {
            fixed.push(3 * d + c);
        }
    }
    let free_vector = mounted
        .free_roots
        .iter()
        .flat_map(|&d| (0..3).map(move |c| 3 * d + c))
        .collect::<Vec<_>>();
    let mut entries = Vec::<[String; 4]>::new();
    for line in &mounted.lines {
        if line.offset != Ratio::zero() {
            return Err("adaptive y50 reduction has nonzero prescribed offset");
        }
        for &(master, r) in &line.terms {
            for c in 0..3 {
                entries.push([
                    (3 * line.dof + c).to_string(),
                    (c * nfree + master).to_string(),
                    r.numerator().to_string(),
                    r.denominator().to_string(),
                ]);
            }
        }
    }
    let clamping = object([
        ("feature", text("y50")),
        ("scalar_dofs", numbers(y50)),
        ("vector_dofs", numbers(&fixed)),
        (
            "component_order",
            Json::Array(["ux", "uy", "uz"].into_iter().map(text).collect()),
        ),
        ("value_mm", numbers([0, 0, 0])),
    ]);
    let reduction = object([
        ("full_vector_dofs", number(nfull)),
        ("reduced_vector_dofs", number(3 * nfree)),
        ("ordering", text("scalar_node_major_xyz")),
        ("reduced_ordering", text("component_major_free_roots")),
        ("free_master_scalar_dofs", numbers(&mounted.free_roots)),
        ("fixed_vector_dofs", numbers(&fixed)),
        ("free_vector_dofs", numbers(&free_vector)),
        (
            "C",
            object([
                ("shape", numbers([nfull, 3 * nfree])),
                (
                    "representation",
                    text("sparse_exact_scalar_lift_xyz_interleaved_rows_component_major_columns"),
                ),
                (
                    "unit_entries_row_column_numerator_denominator",
                    arr(entries, numbers),
                ),
                ("g_is_zero", Json::Bool(true)),
            ]),
        ),
        ("clamped_nodes", clamping.clone()),
    ]);
    Ok((clamping, reduction, lift_json(&adaptive.lift)))
}

fn leaf_domain(leaf: &LeafKey) -> Json {
    let scale = 1_u32 << leaf.level;
    object([
        ("numerator", arr(leaf.index, |i| numbers([i, i + 1]))),
        ("denominator", number(scale)),
    ])
}
fn design_face_json(
    adaptive: &AdaptiveSpace,
    patch: &ExteriorPatch,
    feature: &str,
    denom: u32,
) -> Result<Json, &'static str> {
    let side = side_ref_json(adaptive, &patch.side)?;
    let Json::Object(mut fields) = side else {
        unreachable!()
    };
    fields.insert("feature".into(), text(feature));
    fields.insert("root_rect".into(), rect_json(patch.root_rect));
    fields.insert("grid_denominator".into(), number(denom));
    Ok(Json::Object(fields))
}
fn adaptive_mesh_json(
    adaptive: &AdaptiveSpace,
    partition: &AdaptivePartition,
    certificate: &BracketCertificate,
    tags: &BTreeMap<BlockSide, String>,
    plan: &RefinementPlan,
) -> Result<(Json, Vec<Json>), &'static str> {
    let oriented_complex = adaptive_complex::build_oriented_complex(adaptive, partition)?;
    let report = partition.snapshot()?;
    let denom = partition.grid_denominator();
    let leaves = report
        .leaves
        .iter()
        .map(|leaf| {
            let id = adaptive
                .leaf_blocks
                .get(leaf)
                .ok_or("adaptive leaf has no block")?;
            Ok(object([
                ("root", text(&leaf.root)),
                ("level", number(leaf.level)),
                ("index", numbers(leaf.index)),
                ("block_id", text(id)),
                ("domain", leaf_domain(leaf)),
                (
                    "determinant_lower_bound_mm3",
                    number(
                        certificate.determinant_lower_bound_mm3
                            * 2.0_f64.powi(-3 * i32::from(leaf.level)),
                    ),
                ),
            ]))
        })
        .collect::<Result<Vec<_>, &'static str>>()?;
    let interfaces = report
        .interfaces
        .iter()
        .map(|face| face_json(adaptive, face, denom))
        .collect::<Result<Vec<_>, _>>()?;
    let exterior = report
        .exterior
        .iter()
        .map(|face| {
            exterior_json(
                adaptive,
                face,
                denom,
                tags.get(&face.side.root_side).map(String::as_str),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let design = report
        .exterior
        .iter()
        .filter_map(|face| {
            tags.get(&face.side.root_side)
                .map(|tag| design_face_json(adaptive, face, tag, denom))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let steps = plan
        .steps
        .iter()
        .map(|step| object([("refine", arr(step, leaf_json))]))
        .collect::<Vec<_>>();
    let mesh = object([
        ("schema", text("adaptive_q2_mesh_v1")),
        ("max_depth", number(partition.max_depth())),
        ("max_cells", number(partition.max_cells())),
        ("grid_denominator", number(denom)),
        ("refinement_step_count", number(plan.steps.len())),
        (
            "refinement_plan",
            object([
                ("schema", text("adaptive_refinement_plan_v1")),
                ("max_depth", number(plan.max_depth)),
                ("max_cells", number(plan.max_cells)),
                ("steps", Json::Array(steps)),
            ]),
        ),
        ("leaves", Json::Array(leaves)),
        ("interfaces", Json::Array(interfaces)),
        ("exterior_patches", Json::Array(exterior)),
        ("oriented_complex", oriented_complex.to_json()),
        ("exact_lift", lift_json(&adaptive.lift)),
        (
            "geometry_relation",
            text(
                "each leaf is the exact homogeneous dyadic restriction of its certified root chart",
            ),
        ),
        ("design_faces", Json::Array(design.clone())),
        (
            "root_determinant_lower_bound_mm3",
            number(certificate.determinant_lower_bound_mm3),
        ),
    ]);
    Ok((mesh, design))
}

fn replace_object_field(json: &mut Json, key: &str, value: Json) -> Result<(), &'static str> {
    let Json::Object(fields) = json else {
        return Err("expected package object");
    };
    fields.insert(key.to_owned(), value);
    Ok(())
}
fn get_object_mut<'a>(
    json: &'a mut Json,
    key: &str,
) -> Result<&'a mut BTreeMap<String, Json>, &'static str> {
    let Json::Object(fields) = json else {
        return Err("expected package object");
    };
    match fields.get_mut(key) {
        Some(Json::Object(value)) => Ok(value),
        _ => Err("expected nested package object"),
    }
}

/// Export X1's package schema with an exact A1/A2 adaptive mesh extension.
pub fn export_adaptive_package(
    model: &Json,
    plan_json: &Json,
    order: usize,
) -> Result<Json, &'static str> {
    let plan = parse_refinement_plan(plan_json)?;
    let mut package = solver_package::export_analysis(model, order)?;
    let schema = model.field("schema")?.text()?;
    if schema != "filleted_bracket_v1" && schema != "parameterized_bracket_v1" {
        return Err("adaptive export supports fixed and parameterized brackets");
    }
    let blocks = decode_blocks(model)?;
    let cert = if schema == "filleted_bracket_v1" {
        bracket_audit::audit_bracket(&blocks)?
    } else {
        let [rx, ry] = model.field("radii_sixteenths_mm")?.items()? else {
            return Err("invalid bracket radii");
        };
        bracket_parameters_audit::audit_bracket_parameters(rx.number()?, ry.number()?, &blocks)?
            .certificate
    };
    let parent = if schema == "filleted_bracket_v1" {
        analysis_space::build_analysis_space(&blocks)?
    } else {
        let [rx, ry] = model.field("radii_sixteenths_mm")?.items()? else {
            return Err("invalid bracket radii");
        };
        analysis_space::parameterized_bracket_space(rx.number()?, ry.number()?, &blocks)?
    };
    let roots = parent
        .blocks
        .iter()
        .map(|b| b.id.clone())
        .collect::<Vec<_>>();
    let mut partition =
        AdaptivePartition::new(&roots, &cert.incidence, plan.max_depth, plan.max_cells)?;
    for step in &plan.steps {
        partition.refine(step)?;
    }
    let adaptive = adaptive_space::build_adaptive_space(&parent, &partition)?;
    let space = &adaptive.space;
    let tags = cert
        .design_faces
        .iter()
        .cloned()
        .collect::<BTreeMap<_, _>>();
    let (mesh, design) = adaptive_mesh_json(&adaptive, &partition, &cert, &tags, &plan)?;
    let exterior = partition
        .exterior_patches()
        .iter()
        .map(|p| BlockSide {
            block: adaptive.leaf_blocks[&p.side.leaf].clone(),
            axis: p.side.root_side.axis,
            upper: p.side.root_side.upper,
        })
        .collect::<BTreeSet<_>>();
    let patches = partition
        .exterior_patches()
        .iter()
        .filter(|p| tags.get(&p.side.root_side).is_some_and(|f| f == "x60"))
        .map(|p| {
            LoadPatch::whole(BlockSide {
                block: adaptive.leaf_blocks[&p.side.leaf].clone(),
                axis: p.side.root_side.axis,
                upper: p.side.root_side.upper,
            })
        })
        .collect::<Vec<_>>();
    let load = integrate_traction(space, &exterior, &patches, order, |_, _| {
        [0., -1000. / 320., 0.]
    })?;
    if (load.area_mm2 - 320.).abs() > 1e-8
        || (load.resultant_n[1] + 1000.).abs() > 1e-8
        || load.resultant_n[0].abs() > 1e-10
        || load.resultant_n[2].abs() > 1e-10
    {
        return Err("adaptive x60 load does not preserve the certified area and resultant");
    }
    let volume = solver_package::integrate_volume(space, order)?;
    let (clamping, reduction, _) = adaptive_reduction(space, &adaptive)?;
    let skeleton = encode_skeleton(space)?;
    let y_faces = boundary_faces(&adaptive, &partition, &tags, "y50")?;
    let x_faces = boundary_faces(&adaptive, &partition, &tags, "x60")?;
    let Json::Object(fields) = &mut package else {
        return Err("solver package is not an object");
    };
    fields.insert(
        "schema".into(),
        text("solidcraft_adaptive_solver_package_v1"),
    );
    fields.insert("package_extension".into(), text("adaptive_q2_v1"));
    for (key, source_key) in [
        ("derivatives", "source_derivatives"),
        ("topology", "source_topology"),
        ("certificate", "source_certificate"),
    ] {
        let value = fields.remove(key).ok_or("missing source package data")?;
        fields.insert(source_key.into(), value);
    }
    fields.insert("analysis_derivatives".into(), object([
        ("status",text("unavailable")),
        ("scope",text("source_derivatives refers to root charts and the conforming root space only. Adaptive leaf and constraint derivatives have not been exported.")),
    ]));
    fields.insert("skeleton_scope".into(),text("Unsplit edges of individual leaves, identified by shared Q2 node indices. For hanging-face subrectangles use adaptive_mesh.interfaces."));
    fields.insert("adaptive_mesh".into(), mesh);
    fields.insert("design_faces".into(), Json::Array(design));
    fields.insert("skeleton".into(), skeleton);
    let analysis = get_object_mut(&mut package, "analysis_space")?;
    analysis.insert("scalar_ndof".into(), number(space.ndof));
    analysis.insert("vector_ndof".into(), number(3 * space.ndof));
    analysis.insert(
        "blocks".into(),
        arr(&space.blocks, solver_package::analysis_block_json),
    );
    analysis.insert(
        "boundary_dofs".into(),
        Json::Object(
            space
                .boundary_dofs
                .iter()
                .map(|(k, v)| (k.clone(), numbers(v)))
                .collect(),
        ),
    );
    let selections = get_object_mut(&mut package, "selections")?;
    if let Some(Json::Object(y50)) = selections.get_mut("y50") {
        y50.insert("faces".into(), Json::Array(y_faces));
        y50.insert("scalar_dofs".into(), numbers(&space.boundary_dofs["y50"]));
    }
    if let Some(Json::Object(x60)) = selections.get_mut("x60") {
        x60.insert("faces".into(), Json::Array(x_faces));
        x60.insert("area_mm2".into(), number(load.area_mm2));
        x60.insert("resultant_n".into(), numbers(load.resultant_n));
    }
    let load_object = get_object_mut(&mut package, "load")?;
    load_object.insert("area_mm2".into(), number(load.area_mm2));
    load_object.insert("resultant_n".into(), numbers(load.resultant_n));
    load_object.insert("nodal_rhs_n".into(), numbers(&load.nodal_force_n));
    replace_object_field(&mut package, "clamping", clamping)?;
    replace_object_field(&mut package, "reduction", reduction)?;
    let measures = get_object_mut(&mut package, "measures")?;
    measures.insert("volume_quadrature_mm3".into(), number(volume));
    measures.insert("x60_area_mm2".into(), number(load.area_mm2));
    Ok(package)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json::parse;

    fn model(rx: u8) -> Json {
        let recipe = parse(&format!(
            "{{\"family\":\"parameterized_bracket_v1\",\"radii_sixteenths_mm\":[{rx},64]}}"
        ))
        .unwrap();
        crate::interchange::build(&recipe).unwrap()
    }

    fn plan(steps: &str) -> Json {
        parse(&format!(
            "{{\"schema\":\"adaptive_refinement_plan_v1\",\"max_depth\":4,\"max_cells\":1500,\"steps\":[{steps}]}}"
        ))
        .unwrap()
    }

    fn leaf(root: &str, level: u8, index: [u32; 3]) -> Json {
        object([
            ("root", text(root)),
            ("level", number(level)),
            ("index", numbers(index)),
        ])
    }

    fn refinement_plan(levels: u8) -> Json {
        let focus = "hole_x_sector_east_z10_30";
        let mut steps = Vec::new();
        if levels >= 1 {
            steps.push(object([(
                "refine",
                Json::Array(vec![leaf(focus, 0, [0; 3])]),
            )]));
        }
        if levels >= 2 {
            steps.push(object([(
                "refine",
                Json::Array(vec![leaf(focus, 1, [1, 0, 0])]),
            )]));
        }
        object([
            ("schema", text("adaptive_refinement_plan_v1")),
            ("max_depth", number(4)),
            ("max_cells", number(1500)),
            ("steps", Json::Array(steps)),
        ])
    }

    #[test]
    fn exports_radius_four_and_three_point_five_h0_local1_local2() {
        for rx in [64, 56] {
            let model = model(rx);
            for level in 0..=2 {
                let package = export_adaptive_package(&model, &refinement_plan(level), 2).unwrap();
                assert_eq!(
                    package.field("schema").unwrap().text().unwrap(),
                    "solidcraft_adaptive_solver_package_v1"
                );
                assert_eq!(
                    package.field("package_extension").unwrap().text().unwrap(),
                    "adaptive_q2_v1"
                );
                let mesh = package.field("adaptive_mesh").unwrap();
                let leaves = mesh.field("leaves").unwrap().items().unwrap();
                assert_eq!(
                    leaves.len(),
                    package
                        .field("analysis_space")
                        .unwrap()
                        .field("blocks")
                        .unwrap()
                        .items()
                        .unwrap()
                        .len()
                );
                assert_eq!(
                    mesh.field("refinement_step_count")
                        .unwrap()
                        .number::<usize>()
                        .unwrap(),
                    usize::from(level)
                );
                assert_eq!(
                    package
                        .field("load")
                        .unwrap()
                        .field("area_mm2")
                        .unwrap()
                        .number::<f64>()
                        .unwrap(),
                    320.
                );
                assert_eq!(
                    package
                        .field("design_faces")
                        .unwrap()
                        .items()
                        .unwrap()
                        .len(),
                    mesh.field("exterior_patches")
                        .unwrap()
                        .items()
                        .unwrap()
                        .iter()
                        .filter(|p| p.field("feature").is_ok())
                        .count()
                );
                let exact_lift = mesh.field("exact_lift").unwrap();
                assert_eq!(
                    exact_lift
                        .field("master_index_convention")
                        .unwrap()
                        .text()
                        .unwrap(),
                    "zero_based_position_in_free_roots"
                );
                assert_eq!(
                    exact_lift.field("lines").unwrap().items().unwrap().len(),
                    package
                        .field("analysis_space")
                        .unwrap()
                        .field("scalar_ndof")
                        .unwrap()
                        .number::<usize>()
                        .unwrap()
                );
                let reduction = package.field("reduction").unwrap();
                assert!(reduction
                    .field("C")
                    .unwrap()
                    .field("unit_entries_row_column_numerator_denominator")
                    .is_ok());
                assert_eq!(
                    reduction
                        .field("C")
                        .unwrap()
                        .field("shape")
                        .unwrap()
                        .items()
                        .unwrap()[0]
                        .number::<usize>()
                        .unwrap(),
                    3 * package
                        .field("analysis_space")
                        .unwrap()
                        .field("scalar_ndof")
                        .unwrap()
                        .number::<usize>()
                        .unwrap()
                );
            }
        }
    }

    #[test]
    fn rejects_exhausted_budget_and_invalid_active_leaf_atomically() {
        let model = model(64);
        let over_budget=parse(r#"{"schema":"adaptive_refinement_plan_v1","max_depth":4,"max_cells":36,"steps":[{"refine":[{"root":"hole_x_sector_east_z10_30","level":0,"index":[0,0,0]}]}]}"#).unwrap();
        assert!(export_adaptive_package(&model, &over_budget, 2).is_err());
        let bad = plan(r#"{"refine":[{"root":"missing","level":0,"index":[0,0,0]}]}"#);
        assert!(export_adaptive_package(&model, &bad, 2).is_err());
    }
}
