//! A conforming rational tensor-product field on audited whole-face NURBS charts.
use crate::{
    analysis_space::{face_axes, UnionFind},
    bracket::BracketBlock,
    bracket_json::decode_blocks,
    incidence::{enumerate_faces, BlockSide},
    interchange::audit,
    json::{object, Json},
};
use std::collections::BTreeMap;

fn index(q: [usize; 3], degree: [usize; 3]) -> usize {
    (q[0] * (degree[1] + 1) + q[1]) * (degree[2] + 1) + q[2]
}

pub(crate) fn elevate_to(block: &BracketBlock, target: [usize; 3]) -> Vec<[i128; 4]> {
    let mut d = block.degrees;
    let mut net: Vec<_> = block.controls.iter().map(|p| p.map(i128::from)).collect();
    // At most two exact elevations per axis, or a scale factor of 216.
    for axis in 0..3 {
        while d[axis] < target[axis] {
            let old = d;
            d[axis] += 1;
            let n = d[axis];
            let mut next = Vec::new();
            for i in 0..=d[0] {
                for j in 0..=d[1] {
                    for k in 0..=d[2] {
                        let mut a = [i, j, k];
                        let mut b = a;
                        let m = a[axis];
                        a[axis] = m.saturating_sub(1);
                        b[axis] = m.min(n - 1);
                        let left = net[index(a, old)];
                        let right = net[index(b, old)];
                        next.push(std::array::from_fn(|c| {
                            m as i128 * left[c] + (n - m) as i128 * right[c]
                        }));
                    }
                }
            }
            net = next;
        }
    }
    net
}

fn face_index(side: &BlockSide, a: usize, b: usize, p: usize) -> usize {
    let mut q = [0; 3];
    let [s, t] = face_axes(side.axis);
    q[side.axis] = usize::from(side.upper) * p;
    q[s] = a;
    q[t] = b;
    index(q, [p; 3])
}
fn n(v: impl ToString) -> Json {
    Json::Number(v.to_string())
}

/// Return a geometry-degree C0 field on any audited whole-face block model.
///
/// Each block is elevated exactly to the largest geometry degree, then shared
/// face Bernstein coefficients are identified using the audited incidence
/// orientation. The model family itself is not used to infer parent support.
pub fn export_audited_multiblock(model: &Json) -> Result<Json, &'static str> {
    let report = audit(model)?;
    if report.field("status")?.text()? != "valid" {
        return Err("uncertified NURBS geometry");
    }
    let blocks = decode_blocks(model)?;
    let links = enumerate_faces(&blocks)?;
    let p = blocks
        .iter()
        .flat_map(|b| b.degrees)
        .max()
        .ok_or("empty NURBS model")?;
    let count = (p + 1).pow(3);
    let by_id: BTreeMap<_, _> = blocks
        .iter()
        .enumerate()
        .map(|(i, b)| (b.id.as_str(), i))
        .collect();
    let mut union = UnionFind::new(
        blocks
            .len()
            .checked_mul(count)
            .ok_or("NURBS DOF overflow")?,
    );
    for link in &links.interfaces {
        let o = link.orientation;
        let left = by_id[link.first.block.as_str()] * count;
        let right = by_id[link.second.block.as_str()] * count;
        for a in 0..=p {
            for b in 0..=p {
                let s = if o.flip_first { p - a } else { a };
                let t = if o.flip_second { p - b } else { b };
                let (i, j) = if o.swap { (t, s) } else { (s, t) };
                union.join(
                    left + face_index(&link.first, a, b, p),
                    right + face_index(&link.second, i, j, p),
                );
            }
        }
    }
    let mut ids = BTreeMap::new();
    let mut elements = Vec::new();
    for (i, block) in blocks.iter().enumerate() {
        let controls = elevate_to(block, [p; 3]);
        let dofs = (0..count)
            .map(|j| {
                let next = ids.len();
                n(*ids.entry(union.find(i * count + j)).or_insert(next))
            })
            .collect();
        elements.push(object([
            ("id", Json::String(block.id.clone())),
            ("degrees", Json::Array(vec![n(p); 3])),
            (
                "controls",
                Json::Array(
                    controls
                        .iter()
                        .map(|v| Json::Array(v.iter().map(n).collect()))
                        .collect(),
                ),
            ),
            ("dofs", Json::Array(dofs)),
        ]));
    }
    Ok(object([
        (
            "schema",
            Json::String("solidcraft_nurbs_h1_space_v1".into()),
        ),
        ("source_model", model.clone()),
        ("scalar_ndof", n(ids.len())),
        ("tensor_degree", n(p)),
        ("trace_continuity", Json::String("C0".into())),
        ("blocks", Json::Array(elements)),
    ]))
}

/// Keep the original sketch-family allowlist for the stable public command.
pub fn export(model: &Json) -> Result<Json, &'static str> {
    if !matches!(
        model.field("schema")?.text()?,
        "nurbs_multiloop_v1"
            | "g2_cross_plate_v1"
            | "s3_face_pocket_v1"
            | "primitive_sweep_v1"
            | "nurbs_edge_loops_v1"
            | "nurbs_sketch_strip_v1"
            | "nurbs_cross_bores_v1"
            | "nurbs_cross_bores_bounded_v1"
            | "nurbs_cross_oval_v1"
            | "affine_image_v1"
    ) {
        return Err("expected a supported NURBS sketch model");
    }
    export_audited_multiblock(model)
}

/// Attach the physical volume integrated at the requested tensor Gauss order.
pub fn analysis_package(model: &Json, order: usize) -> Result<Json, &'static str> {
    let quadrature = crate::elasticity::gauss(order)?;
    let mut package = export(model)?;
    let mut volume = 0.;
    for b in decode_blocks(model)? {
        let map = crate::RationalBezierMap::new(
            b.degrees,
            b.controls.iter().map(|p| p.map(|x| x as f64)).collect(),
        )
        .map_err(|_| "unsupported NURBS analysis map")?;
        for &(u, a) in &quadrature {
            for &(v, b) in &quadrature {
                for &(w, c) in &quadrature {
                    let [x, y, z] = map
                        .evaluate([u, v, w])
                        .map_err(|_| "invalid NURBS quadrature evaluation")?
                        .jacobian_columns;
                    let det = x[0] * (y[1] * z[2] - y[2] * z[1])
                        + x[1] * (y[2] * z[0] - y[0] * z[2])
                        + x[2] * (y[0] * z[1] - y[1] * z[0]);
                    if !det.is_finite() || det <= 0. {
                        return Err("NURBS quadrature determinant out of range");
                    }
                    volume += a * b * c * det;
                }
            }
        }
    }
    if let Json::Object(fields) = &mut package {
        fields.insert("quadrature_order".into(), n(order));
        fields.insert("integrated_volume_mm3".into(), n(volume));
    }
    Ok(package)
}
