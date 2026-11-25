//! Display samples of audited geometry and its one-sided chart distortion.
use crate::{
    field_vtu::{array, row, write_cells},
    interchange,
    json::Json,
    quality, RationalBezierMap,
};
use std::fmt::Write;

fn map(block: &Json, affine_pair: bool) -> Result<RationalBezierMap, &'static str> {
    let degrees = if affine_pair {
        [1, 1, 1]
    } else {
        let [u, v, w] = block.field("degrees")?.items()? else {
            return Err("invalid geometry chart degrees");
        };
        [u.number()?, v.number()?, w.number()?]
    };
    let mut controls = Vec::new();
    for point in block.field("controls")?.items()? {
        let [x, y, z, w] = point.items()? else {
            return Err("invalid homogeneous chart point");
        };
        controls.push([x.number()?, y.number()?, z.number()?, w.number()?]);
    }
    RationalBezierMap::new(degrees, controls).map_err(|_| "invalid or unsupported geometry chart")
}

/// Each VTK piece retains the native block index. Point diagnostics are samples
/// of the original rational chart, not quality bounds for the linear display cells.
pub fn export_geometry_vtu(model: &Json, n: usize) -> Result<String, &'static str> {
    if !(1..=32).contains(&n) {
        return Err("invalid geometry subdivision count");
    }
    let blocks = model.field("blocks")?.items()?;
    if blocks.is_empty() || blocks.len() > 1_000_000 / (n + 1).pow(3) {
        return Err("geometry sample budget exceeded");
    }
    let audit = interchange::audit(model)?;
    if audit.field("status")?.text()? != "valid" {
        return Err("geometry export requires a valid audit");
    }
    let affine_pair = model.field("schema")?.text()? == "affine_pair_v1";
    let mut out=String::from("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n");
    for (id, block) in blocks.iter().enumerate() {
        let chart = map(block, affine_pair)?;
        let [mut points, mut params, mut det, mut stretches, mut condition, mut scaled] =
            std::array::from_fn(|_| String::new());
        for i in 0..=n {
            for j in 0..=n {
                for k in 0..=n {
                    let t = [i, j, k].map(|i| i as f64 / n as f64);
                    let sample = chart
                        .evaluate(t)
                        .map_err(|_| "geometry sample outside numerical range")?;
                    let jac = sample.jacobian_columns;
                    let q = quality::at_jacobian(jac)?;
                    let [a, b, c] = jac;
                    let determinant = a[0] * (b[1] * c[2] - b[2] * c[1])
                        + a[1] * (b[2] * c[0] - b[0] * c[2])
                        + a[2] * (b[0] * c[1] - b[1] * c[0]);
                    if determinant <= 0. {
                        return Err("nonpositive sampled geometry determinant");
                    }
                    row(&mut points, sample.position)?;
                    row(&mut params, t)?;
                    row(&mut det, [determinant])?;
                    row(&mut stretches, q.singular_values)?;
                    row(&mut condition, [q.condition_number])?;
                    row(&mut scaled, [q.scaled_determinant])?;
                }
            }
        }
        writeln!(
            out,
            r#"<Piece NumberOfPoints="{}" NumberOfCells="{}">"#,
            (n + 1).pow(3),
            n.pow(3)
        )
        .unwrap();
        out.push_str("<Points>\n");
        array(&mut out, "position_mm", "Float64", 3, &points);
        out.push_str("</Points>\n<PointData Scalars=\"condition_number\">\n");
        for (name, components, values) in [
            ("uvw", 3, &params),
            ("jacobian_determinant_mm3", 1, &det),
            ("singular_values_mm", 3, &stretches),
            ("condition_number", 1, &condition),
            ("scaled_jacobian", 1, &scaled),
        ] {
            array(&mut out, name, "Float64", components, values);
        }
        out.push_str("</PointData>\n<CellData Scalars=\"geometry_block_index\">\n");
        array(
            &mut out,
            "geometry_block_index",
            "Int32",
            1,
            &format!("{id}\n").repeat(n.pow(3)),
        );
        out.push_str("</CellData>\n");
        write_cells(&mut out, n);
    }
    out.push_str("</UnstructuredGrid>\n</VTKFile>\n");
    Ok(out)
}
