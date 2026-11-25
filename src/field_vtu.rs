//! Sample a solved rational displacement field on separate linear hexahedral charts.
use crate::{
    analysis_space::AnalysisBlock,
    elasticity::{basis_geometry, Matrix},
    json::Json,
    projected_bulk::ProjectedElasticity,
    stress::von_mises,
};
use std::{collections::BTreeSet, fmt::Write};

fn number(v: &Json) -> Result<f64, &'static str> {
    let x = v.number::<f64>()?;
    if x.is_finite() {
        Ok(x)
    } else {
        Err("nonfinite result coefficient")
    }
}
fn block(v: &Json, ndof: usize) -> Result<AnalysisBlock, &'static str> {
    v.keys(&["id", "degrees", "control_order", "controls", "dofs"])?;
    if v.field("control_order")?.text()? != "u-major,v-major,w-fast"
        || v.field("degrees")?
            .items()?
            .iter()
            .map(Json::number::<usize>)
            .collect::<Result<Vec<_>, _>>()?
            != [2, 2, 2]
    {
        return Err("unsupported field control ordering or degree");
    }
    let controls = v
        .field("controls")?
        .items()?
        .iter()
        .map(|p| {
            let [a, b, c, w] = p.items()? else {
                return Err("expected homogeneous field point");
            };
            Ok([a.number()?, b.number()?, c.number()?, w.number()?])
        })
        .collect::<Result<Vec<_>, _>>()?;
    let dofs = v
        .field("dofs")?
        .items()?
        .iter()
        .map(Json::number)
        .collect::<Result<Vec<usize>, _>>()?;
    if controls.len() != 27 || dofs.len() != 27 || dofs.iter().any(|&id| id >= ndof) {
        return Err("invalid field block connectivity");
    }
    Ok(AnalysisBlock {
        id: v.field("id")?.text()?.to_owned(),
        controls,
        dofs,
    })
}

struct Material {
    mu: f64,
    bulk: f64,
}
struct ChartData {
    points: String,
    u: String,
    stress: String,
    mises: String,
    params: String,
}
fn row<const N: usize>(storage: &mut String, values: [f64; N]) -> Result<(), &'static str> {
    for x in values {
        if !x.is_finite() {
            return Err("nonfinite sampled field");
        }
        write!(storage, "{x:.16e} ").unwrap();
    }
    storage.push('\n');
    Ok(())
}
fn sample_chart(
    b: &AnalysisBlock,
    u: &[f64],
    material: &Material,
    pressure: Option<[f64; 8]>,
    n: usize,
) -> Result<ChartData, &'static str> {
    let mut data = ChartData {
        points: String::new(),
        u: String::new(),
        stress: String::new(),
        mises: String::new(),
        params: String::new(),
    };
    for i in 0..=n {
        for j in 0..=n {
            for k in 0..=n {
                let uvw = [i, j, k].map(|v| v as f64 / n as f64);
                let sample = basis_geometry(b, uvw)?;
                let mut displacement = [0.; 3];
                let mut grad: Matrix = [[0.; 3]; 3];
                for (local, &id) in b.dofs.iter().enumerate() {
                    for a in 0..3 {
                        let ua = u[3 * id + a];
                        displacement[a] += ua * sample.values[local];
                        for (component, derivative) in grad[a].iter_mut().zip(sample.gradients[local]) {
                            *component += ua * derivative;
                        }
                    }
                }
                let div = grad[0][0] + grad[1][1] + grad[2][2];
                let mean = match pressure {
                    Some(p) => ProjectedElasticity::mean_stress_at(p, uvw)?,
                    None => material.bulk * div,
                };
                let stress: Matrix = std::array::from_fn(|a| {
                    std::array::from_fn(|d| {
                        let shear = material.mu * (grad[a][d] + grad[d][a]);
                        if a == d {
                            shear + mean - (2. / 3.) * material.mu * div
                        } else {
                            shear
                        }
                    })
                });
                row(&mut data.points, sample.position)?;
                row(&mut data.u, displacement)?;
                row(&mut data.params, uvw)?;
                row(
                    &mut data.stress,
                    std::array::from_fn::<_, 9, _>(|v| stress[v / 3][v % 3]),
                )?;
                row(&mut data.mises, [von_mises(stress)])?;
            }
        }
    }
    Ok(data)
}
fn array(out: &mut String, name: &str, ty: &str, n: usize, content: &str) {
    writeln!(
        out,
        r#"<DataArray type="{ty}" Name="{name}" NumberOfComponents="{n}" format="ascii">"#
    )
    .unwrap();
    out.push_str(content);
    out.push_str("</DataArray>\n");
}
fn write_piece(out: &mut String, data: ChartData, n: usize, chart: usize) {
    let stride = n + 1;
    writeln!(
        out,
        r#"<Piece NumberOfPoints="{}" NumberOfCells="{}">"#,
        stride.pow(3),
        n.pow(3)
    )
    .unwrap();
    out.push_str("<Points>\n");
    array(out, "position_mm", "Float64", 3, &data.points);
    out.push_str("</Points>\n<PointData Vectors=\"displacement_mm\" Tensors=\"stress_mpa\" Scalars=\"von_mises_mpa\">\n");
    array(out, "displacement_mm", "Float64", 3, &data.u);
    array(out, "stress_mpa", "Float64", 9, &data.stress);
    array(out, "von_mises_mpa", "Float64", 1, &data.mises);
    array(out, "uvw", "Float64", 3, &data.params);
    out.push_str("</PointData>\n<CellData>\n");
    let mut ids = String::new();
    for _ in 0..n.pow(3) {
        writeln!(ids, "{chart}").unwrap();
    }
    array(out, "analysis_block_index", "Int32", 1, &ids);
    out.push_str("</CellData>\n<Cells>\n");
    let mut conn = String::new();
    let offsets = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                for [a, b, c] in offsets {
                    write!(conn, "{} ", ((i + a) * stride + j + b) * stride + k + c).unwrap();
                }
                conn.push('\n');
            }
        }
    }
    array(out, "connectivity", "Int64", 1, &conn);
    ids.clear();
    for i in 1..=n.pow(3) {
        writeln!(ids, "{}", 8 * i).unwrap();
    }
    array(out, "offsets", "Int64", 1, &ids);
    array(out, "types", "UInt8", 1, &"12\n".repeat(n.pow(3)));
    out.push_str("</Cells>\n</Piece>\n");
}

/// Sample the serialized solution. Chart indices refer to the ordered analysis_blocks array.
pub fn export_field_vtu(result: &Json, n: usize) -> Result<String, &'static str> {
    if result.field("schema")?.text()? != "solidcraft_mechanical_result_v1"
        || !(1..=32).contains(&n)
    {
        return Err("invalid field result schema or subdivision count");
    }
    let blocks = result.field("analysis_blocks")?.items()?;
    if blocks.is_empty() || blocks.len() > 1_000_000 / (n + 1).pow(3) {
        return Err("field sample budget exceeded");
    }
    let ndof = result.field("scalar_dofs")?.number::<usize>()?;
    let u = result
        .field("displacement_mm_xyz")?
        .items()?
        .iter()
        .map(number)
        .collect::<Result<Vec<_>, _>>()?;
    if ndof.checked_mul(3) != Some(u.len()) {
        return Err("invalid field displacement count");
    }
    let job = result.field("job")?;
    let material = job.field("material")?;
    let young = number(material.field("young_mpa")?)?;
    let nu = number(material.field("poisson")?)?;
    if young <= 0. || !(-1. < nu && nu < 0.5) {
        return Err("invalid field material");
    }
    let material = Material {
        mu: young / (2. * (1. + nu)),
        bulk: young / (3. * (1. - 2. * nu)),
    };
    let projected = match job.field("formulation")?.text()? {
        "displacement" => false,
        "projected_bulk_q1" => true,
        _ => return Err("unknown field formulation"),
    };
    let mut seen = BTreeSet::new();
    let mut out=String::from("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n");
    for (index, item) in blocks.iter().enumerate() {
        let b = block(item, ndof)?;
        if !seen.insert(b.id.clone()) {
            return Err("duplicate field chart");
        }
        let p = if projected {
            let values = result
                .field("mean_stress_q1_mpa_by_block")?
                .field(&b.id)?
                .items()?
                .iter()
                .map(number)
                .collect::<Result<Vec<_>, _>>()?;
            Some(
                values
                    .try_into()
                    .map_err(|_| "invalid recovered Q1 pressure length")?,
            )
        } else {
            None
        };
        write_piece(&mut out, sample_chart(&b, &u, &material, p, n)?, n, index);
    }
    out.push_str("</UnstructuredGrid>\n</VTKFile>\n");
    Ok(out)
}
