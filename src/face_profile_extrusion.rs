//! Native curved charts for a rational planar profile swept from a stock face.
//!
//! This module keeps the source Bézier curves in the tensor-product maps. It
//! never replaces a curved input edge by a chord. The rational controls are
//! reconstructed on a fixed declared lattice and that reconstruction bound
//! is checked against the caller's tolerance.

use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    bracket_json::encode,
    incidence::enumerate_faces,
    json::{Json, object},
};

const FAMILY: &str = "face_profile_extrusion_v1";
const SCALE: f64 = 1_000_000_000.0;
// A common gauge for halves through degree five and line elevation through degree five.
const SPAN_SCALE: i128 = 480;
const MAX_PATCHES: usize = 1024;

/// One planar tensor patch, degree `degree` in u and degree one in v.
/// Controls are homogeneous XYW, flattened as `u * 2 + v`; v=0 is the
/// source curve and v=1 is the inner chord.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SectionPatch {
    pub degree: usize,
    pub controls: Vec<[i64; 3]>,
}

/// Extrusion request in a right-handed rigid sketch frame.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceProfileExtrusion {
    /// Ordered rational Bézier spans, homogeneous (xw,yw,w) controls.
    pub spans: Vec<Vec<[f64; 3]>>,
    /// `[[xmin,xmax],[ymin,ymax]]`, in local sketch coordinates.
    pub bounds: [[f64; 2]; 2],
    pub thickness: f64,
    pub depth: f64,
    /// `add`, `cut`, or `tool`.
    pub operation: String,
    /// `[origin, x-axis, y-axis, outward normal]`.
    pub frame: [[f64; 3]; 4],
    /// Maximum permitted homogeneous reconstruction error in mm.
    pub tolerance: f64,
}

fn point(h: [f64; 3]) -> Result<[f64; 2], &'static str> {
    if !h.iter().all(|v| v.is_finite()) || h[2] <= 0.0 {
        return Err("profile_controls_require_finite_positive_weights");
    }
    Ok([h[0] / h[2], h[1] / h[2]])
}
fn q(v: f64) -> Result<i64, &'static str> {
    if !v.is_finite() || v.abs() > i64::MAX as f64 / SCALE {
        return Err("profile_rational_reconstruction_overflow");
    }
    Ok((v * SCALE).round() as i64)
}
fn normalize_row(row: &[[f64; 3]]) -> Result<Vec<[i64; 3]>, &'static str> {
    let raw = row
        .iter()
        .map(|p| Ok([q(p[0])?, q(p[1])?, q(p[2])?]))
        .collect::<Result<Vec<_>, &'static str>>()?;
    if raw.iter().any(|p| p[2] <= 0) {
        return Err("profile_reconstruction_nonpositive_weight");
    }
    Ok(raw)
}
type HomogeneousSpan = Vec<[i64; 3]>;
fn split_half(row: &[[i64; 3]]) -> Result<(HomogeneousSpan, HomogeneousSpan), &'static str> {
    let mut a = row.iter().map(|p| p.map(i128::from)).collect::<Vec<_>>();
    let scale = SPAN_SCALE;
    let mut left = vec![a[0].map(|x| x * scale)];
    let mut right = vec![a[a.len() - 1].map(|x| x * scale)];
    let mut level = 0_i128;
    while a.len() > 1 {
        level += 1;
        a = a
            .windows(2)
            .map(|p| std::array::from_fn(|k| p[0][k] + p[1][k]))
            .collect();
        let factor = scale >> level;
        left.push(a[0].map(|x| x * factor));
        right.push(a[a.len() - 1].map(|x| x * factor));
    }
    right.reverse();
    let cvt = |v: Vec<[i128; 3]>| -> Result<Vec<[i64; 3]>, &'static str> {
        v.into_iter()
            .map(|p| {
                Ok([
                    i64::try_from(p[0]).map_err(|_| "profile_split_overflow")?,
                    i64::try_from(p[1]).map_err(|_| "profile_split_overflow")?,
                    i64::try_from(p[2]).map_err(|_| "profile_split_overflow")?,
                ])
            })
            .collect()
    };
    Ok((cvt(left)?, cvt(right)?))
}
fn elevate_chord(a: [i64; 3], b: [i64; 3], degree: usize) -> Result<Vec<[i64; 3]>, &'static str> {
    (0..=degree)
        .map(|i| {
            let k = i as i128;
            let x = degree as i128;
            let factor = SPAN_SCALE / x;
            let mix = |d: usize| -> Result<i64, &'static str> {
                let n = (i128::from(a[d]) * (x - k) + i128::from(b[d]) * k) * factor;
                i64::try_from(n).map_err(|_| "profile_chord_overflow")
            };
            Ok([mix(0)?, mix(1)?, mix(2)?])
        })
        .collect()
}
fn patch(outer: Vec<[i64; 3]>, inner: Vec<[i64; 3]>) -> Result<SectionPatch, &'static str> {
    if outer.len() != inner.len() || !(2..=6).contains(&outer.len()) {
        return Err("profile_span_degree_out_of_range");
    }
    let degree = outer.len() - 1;
    let mut controls = Vec::with_capacity(outer.len() * 2);
    for i in 0..outer.len() {
        controls.push(outer[i]);
        controls.push(inner[i]);
    }
    if controls.iter().any(|p| p[2] <= 0) {
        return Err("profile_patch_nonpositive_weight");
    }
    Ok(SectionPatch { degree, controls })
}
fn verify_curve_row(row: &[[i64; 3]]) -> Result<(), &'static str> {
    if row.len() < 2 || row.iter().any(|p| p[2] <= 0) {
        return Err("profile_patch_nonpositive_weight");
    }
    Ok(())
}

/// Make exact curved annular patches and nonsingular center-fan quadrilaterals.
/// Each source span is split at one half so its inner chord shares a midpoint
/// vertex with the core mesh. The kernel center is the mean of unique boundary
/// vertices and must lie strictly to the left of every oriented polygon edge.
pub fn curved_ring_section(
    spans: &[Vec<[f64; 3]>],
    tolerance: f64,
) -> Result<(Vec<SectionPatch>, Vec<SectionPatch>), &'static str> {
    if spans.is_empty()
        || spans.len() > MAX_PATCHES / 2
        || !tolerance.is_finite()
        || tolerance <= 0.0
    {
        return Err("face_profile_span_or_tolerance_range");
    }
    let rows = spans
        .iter()
        .map(|s| {
            if !(2..=6).contains(&s.len()) {
                return Err("profile_degree_must_be_one_to_five");
            }
            normalize_row(s)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut reconstruction_bound = 0.0_f64;
    let mut min_weight = f64::INFINITY;
    let mut projected = Vec::new();
    for (source, rounded) in spans.iter().zip(&rows) {
        for (a, b) in source.iter().zip(rounded) {
            let pa = point(*a)?;
            let pb = point(b.map(|v| v as f64))?;
            reconstruction_bound = reconstruction_bound.max((pa[0] - pb[0]).hypot(pa[1] - pb[1]));
            min_weight = min_weight.min(a[2]);
            projected.push(pa);
        }
    }
    let diameter = projected
        .iter()
        .flat_map(|a| {
            projected
                .iter()
                .map(move |b| (a[0] - b[0]).hypot(a[1] - b[1]))
        })
        .fold(0.0_f64, f64::max);
    // Each homogeneous coefficient is rounded by at most 0.5/SCALE. The
    // first term bounds control point motion; the second bounds the change in
    // normalized rational Bernstein weights over the whole parameter span.
    let eta = 0.5 / (SCALE * min_weight);
    let rational_bound = reconstruction_bound + 2.0 * diameter * eta / (1.0 - eta);
    if !eta.is_finite() || eta >= 1.0 || rational_bound > tolerance * 0.25 {
        return Err("profile_rational_reconstruction_exceeds_declared_tolerance");
    }
    for row in &rows {
        verify_curve_row(row)?;
    }
    for i in 0..rows.len() {
        if *rows[i].last().ok_or("empty_profile_span")? != rows[(i + 1) % rows.len()][0] {
            return Err("profile_spans_must_join_at_shared_exact_points");
        }
    }
    let verts = rows
        .iter()
        .map(|h| point(h[0].map(|v| v as f64)))
        .collect::<Result<Vec<_>, _>>()?;
    let center = [
        verts.iter().map(|p| p[0]).sum::<f64>() / verts.len() as f64,
        verts.iter().map(|p| p[1]).sum::<f64>() / verts.len() as f64,
    ];
    let mut area = 0.0;
    for i in 0..verts.len() {
        let a = verts[i];
        let b = verts[(i + 1) % verts.len()];
        area += a[0] * b[1] - a[1] * b[0];
    }
    if area < 0.0 {
        let reversed = spans
            .iter()
            .rev()
            .map(|s| s.iter().rev().copied().collect())
            .collect::<Vec<_>>();
        return curved_ring_section(&reversed, tolerance);
    }
    if area == 0.0 {
        return Err("profile_has_zero_or_ambiguous_orientation");
    }
    for i in 0..verts.len() {
        let a = verts[i];
        let b = verts[(i + 1) % verts.len()];
        let va = [a[0] - center[0], a[1] - center[1]];
        let vb = [b[0] - center[0], b[1] - center[1]];
        let turn = va[0] * vb[1] - va[1] * vb[0];
        if turn <= 0.0 {
            return Err("profile_kernel_sectors_must_have_positive_turn");
        }
        let side = (b[0] - a[0]) * (center[1] - a[1]) - (b[1] - a[1]) * (center[0] - a[0]);
        if side <= tolerance * 1e-9 {
            return Err("profile_vertices_have_no_strict_common_kernel_center");
        }
        let span = &rows[i];
        for h in span
            .iter()
            .map(|p| point(p.map(|v| v as f64)))
            .collect::<Result<Vec<_>, _>>()?
        {
            let vh = [h[0] - center[0], h[1] - center[1]];
            let left = va[0] * vh[1] - va[1] * vh[0];
            let right = vh[0] * vb[1] - vh[1] * vb[0];
            if left < -tolerance * 1e-8 || right < -tolerance * 1e-8 {
                return Err("profile_span_control_hull_crosses_its_kernel_sector");
            }
        }
    }
    let winding = verts
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let b = verts[(i + 1) % verts.len()];
            let va = [a[0] - center[0], a[1] - center[1]];
            let vb = [b[0] - center[0], b[1] - center[1]];
            (va[0] * vb[1] - va[1] * vb[0]).atan2(va[0] * vb[0] + va[1] * vb[1])
        })
        .sum::<f64>();
    if (winding - std::f64::consts::TAU).abs() > 1e-8 {
        return Err("profile_boundary_must_wind_once_around_kernel");
    }
    for i in 0..verts.len() {
        for j in i + 1..verts.len() {
            if j == i + 1 || (i == 0 && j + 1 == verts.len()) {
                continue;
            }
            let a = verts[i];
            let b = verts[(i + 1) % verts.len()];
            let c = verts[j];
            let d = verts[(j + 1) % verts.len()];
            let orient = |p: [f64; 2], q: [f64; 2], r: [f64; 2]| {
                (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
            };
            if orient(a, b, c) * orient(a, b, d) <= 0.0 && orient(c, d, a) * orient(c, d, b) <= 0.0
            {
                return Err("profile_boundary_self_intersects");
            }
        }
    }
    let center_h = [q(center[0])?, q(center[1])?, q(1.0)?];
    let mut inner = Vec::with_capacity(rows.len());
    for p in rows.iter().map(|r| r[0]) {
        let xy = point(p.map(|v| v as f64))?;
        inner.push([
            q((xy[0] + center[0]) * 0.5)?,
            q((xy[1] + center[1]) * 0.5)?,
            q(1.0)?,
        ]);
    }
    let mut ring = Vec::with_capacity(rows.len() * 2);
    let mut midpoints = Vec::with_capacity(rows.len());
    for i in 0..rows.len() {
        let (left, right) = split_half(&rows[i])?;
        let a = inner[i];
        let b = inner[(i + 1) % inner.len()];
        let mid = [
            q((a[0] as f64 + b[0] as f64) / (2.0 * SCALE))?,
            q((a[1] as f64 + b[1] as f64) / (2.0 * SCALE))?,
            q(1.0)?,
        ];
        midpoints.push(mid);
        let (l, rr) = (left, right);
        ring.push(patch(l, elevate_chord(a, mid, rows[i].len() - 1)?)?);
        ring.push(patch(rr, elevate_chord(mid, b, rows[i].len() - 1)?)?);
    }
    let mut core = Vec::with_capacity(rows.len());
    for i in 0..rows.len() {
        let scale = |mut p: [i64; 3]| -> Result<[i64; 3], &'static str> {
            for v in &mut p {
                *v = i64::try_from(i128::from(*v) * SPAN_SCALE)
                    .map_err(|_| "profile_core_overflow")?;
            }
            Ok(p)
        };
        let vi = scale(inner[i])?;
        let prev = scale(midpoints[(i + inner.len() - 1) % inner.len()])?;
        let next = scale(midpoints[i])?;
        let c = scale(center_h)?;
        core.push(SectionPatch {
            degree: 1,
            controls: vec![vi, prev, next, c],
        });
    }
    Ok((ring, core))
}

/// Lift section patches into native homogeneous trivariate blocks.
pub fn lift_sections(
    patches: &[SectionPatch],
    depth: f64,
    frame: [[f64; 3]; 4],
    reverse_depth: bool,
) -> Result<Vec<BracketBlock>, &'static str> {
    if !depth.is_finite() || depth <= 0.0 || patches.is_empty() || patches.len() > MAX_PATCHES {
        return Err("invalid_profile_lift_parameters");
    }
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    for &axis in &frame[1..] {
        if !axis.iter().all(|v| v.is_finite()) || (dot(axis, axis) - 1.0).abs() > 1e-10 {
            return Err("profile_frame_axes_must_be_unit");
        }
    }
    if !frame[0].iter().all(|v| v.is_finite())
        || [
            dot(frame[1], frame[2]),
            dot(frame[1], frame[3]),
            dot(frame[2], frame[3]),
        ]
        .iter()
        .any(|v| v.abs() > 1e-10)
        || cross(frame[1], frame[2])
            .iter()
            .zip(frame[3])
            .any(|(a, b)| (*a - b).abs() > 1e-10)
    {
        return Err("profile_frame_must_be_right_handed_orthonormal");
    }
    patches
        .iter()
        .enumerate()
        .map(|(idx, p)| {
            if p.controls.len() != (p.degree + 1) * 2
                || !(1..=5).contains(&p.degree)
                || p.controls.iter().any(|x| x[2] <= 0)
            {
                return Err("invalid_profile_section_patch");
            }
            let controls = crate::face_profile_frame::lift(
                p,
                if reverse_depth { -depth } else { depth },
                frame,
            )?;
            let block = BracketBlock {
                id: format!("profile_{idx}"),
                degrees: [p.degree, 1, 1],
                controls,
            };
            Ok(block)
        })
        .collect()
}

/// Build a native profile tool or delegate the enclosing-stock partition.
pub fn build(r: &FaceProfileExtrusion) -> Result<Json, &'static str> {
    if !r.thickness.is_finite()
        || r.thickness <= 0.0
        || !r.depth.is_finite()
        || r.depth <= 0.0
        || !matches!(r.operation.as_str(), "tool" | "add" | "cut")
        || r.bounds.iter().any(|b| b[0] >= b[1])
    {
        return Err("invalid_face_profile_extrusion_parameters");
    }
    let frame_bound = frame_rounding_bound(r)?;
    if frame_bound + 2.0 / SCALE > r.tolerance * 0.25 {
        return Err("face_frame_rationalization_exceeds_declared_tolerance");
    }
    let (ring, core) = curved_ring_section(&r.spans, r.tolerance)?;
    for span in &r.spans {
        for h in span {
            let p = point(*h)?;
            if (0..2).any(|k| p[k] <= r.bounds[k][0] || p[k] >= r.bounds[k][1]) {
                return Err("profile_control_hull_must_be_strictly_inside_stock_face");
            }
        }
    }
    let mut blocks = if r.operation == "tool" {
        let mut all = [
            lift_sections(&ring, r.depth, r.frame, false)?,
            lift_sections(&core, r.depth, r.frame, false)?,
        ]
        .concat();
        for (i, b) in all.iter_mut().enumerate() {
            b.id = format!("profile_tool_{i}");
        }
        all
    } else {
        crate::face_profile_stock::layers(r, &core, &ring)?
    };
    let incidence = enumerate_faces(&blocks)?;
    let mut lower = f64::INFINITY;
    for b in &blocks {
        let c =
            crate::spatial_certificate::certify_spatial(b.degrees, &b.controls, Default::default());
        let Some(v) = c.determinant_lower_bound_mm3 else {
            return Err("face_profile_jacobian_unresolved");
        };
        lower = lower.min(v);
    }
    if !lower.is_finite() || lower <= 0.0 {
        return Err("face_profile_jacobian_nonpositive");
    }
    let mut cert = BracketCertificate {
        incidence,
        design_faces: Vec::new(),
        determinant_lower_bound_mm3: lower,
        radial_error_upper_mm: r.tolerance,
    };
    for side in &cert.incidence.exterior {
        cert.design_faces
            .push((side.clone(), "profile_tool".into()));
    }
    let mut model = encode(FAMILY, std::mem::take(&mut blocks), &cert);
    if let Json::Object(fields) = &mut model {
        fields.insert("source_recipe".into(), recipe_json(r));
        fields.insert(
            "determinant_lower_bound_mm3".into(),
            Json::Number(lower.to_string()),
        );
        fields.insert(
            "hausdorff_error_upper_bound_mm".into(),
            Json::Number(r.tolerance.to_string()),
        );
    }
    Ok(model)
}
fn frame_rounding_bound(r: &FaceProfileExtrusion) -> Result<f64, &'static str> {
    let mut extent = 0.0_f64;
    for span in &r.spans {
        for h in span {
            let p = point(*h)?;
            extent = extent.max(p[0].abs()).max(p[1].abs());
        }
    }
    if r.operation != "tool" {
        for b in r.bounds {
            extent = extent.max(b[0].abs()).max(b[1].abs());
        }
    }
    // Each frame-axis coefficient and each already-translated layer origin
    // is rounded once, in world coordinates. Include the stock perimeter grid.
    let grid = 0.5 / SCALE;
    let world = 3.0_f64.sqrt() * grid * (2.0 * extent + 1.0);
    let stock = if r.operation == "tool" {
        0.0
    } else {
        2.0_f64.sqrt() * grid
    };
    Ok(world + stock)
}
fn recipe_json(r: &FaceProfileExtrusion) -> Json {
    let arr = |v: &[f64]| Json::Array(v.iter().map(|x| Json::Number(x.to_string())).collect());
    object([
        ("family", Json::String(FAMILY.into())),
        ("operation", Json::String(r.operation.clone())),
        (
            "bounds_mm",
            Json::Array(r.bounds.iter().map(|v| arr(v)).collect()),
        ),
        ("thickness_mm", Json::Number(r.thickness.to_string())),
        ("depth_mm", Json::Number(r.depth.to_string())),
        ("tolerance_mm", Json::Number(r.tolerance.to_string())),
        (
            "frame",
            Json::Array(r.frame.iter().map(|v| arr(v)).collect()),
        ),
        (
            "spans",
            Json::Array(
                r.spans
                    .iter()
                    .map(|s| Json::Array(s.iter().map(|p| arr(p)).collect()))
                    .collect(),
            ),
        ),
    ])
}
/// Reconstruct the native maps from the stored recipe and reject altered nets.
pub fn audit(model: &Json) -> Result<Json, &'static str> {
    crate::face_profile_audit::audit(model)
}

#[cfg(test)]
mod tests {
    use super::*;
    const FRAME: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    #[test]
    fn exact_convex_polygon_tool_has_positive_curved_ring_and_core_maps() {
        let spans = vec![
            vec![[0., 0., 1.], [2., 0., 1.]],
            vec![[2., 0., 1.], [2., 2., 1.]],
            vec![[2., 2., 1.], [0., 2., 1.]],
            vec![[0., 2., 1.], [0., 0., 1.]],
        ];
        let (ring, core) = curved_ring_section(&spans, 1e-7).unwrap();
        assert_eq!(ring.len(), 8);
        assert_eq!(core.len(), 4);
        let blocks = [
            lift_sections(&ring, 3., FRAME, false).unwrap(),
            lift_sections(&core, 3., FRAME, false).unwrap(),
        ]
        .concat();
        assert_eq!(blocks.len(), 12);
        assert!(blocks.iter().all(|b| {
            crate::spatial_certificate::certify_spatial(b.degrees, &b.controls, Default::default())
                .determinant_lower_bound_mm3
                .is_some_and(|x| x > 0.)
        }));
    }
    #[test]
    fn rational_quarter_circle_remains_curved_in_the_native_controls() {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let spans = vec![
            vec![[1., 0., 1.], [1., 1., w], [0., 1., 1.]],
            vec![[0., 1., 1.], [-1., 1., w], [-1., 0., 1.]],
            vec![[-1., 0., 1.], [-1., -1., w], [0., -1., 1.]],
            vec![[0., -1., 1.], [1., -1., w], [1., 0., 1.]],
        ];
        let (ring, core) = curved_ring_section(&spans, 1e-6).unwrap();
        assert_eq!(ring.len(), 8);
        assert_eq!(core.len(), 4);
        // A source middle control is off the endpoint chord: the boundary row
        // is still quadratic after exact half-span extraction.
        assert!(ring.iter().any(|p| p.degree == 2));
        let blocks = [
            lift_sections(&ring, 1., FRAME, false).unwrap(),
            lift_sections(&core, 1., FRAME, false).unwrap(),
        ]
        .concat();
        assert!(blocks.iter().all(|b| {
            crate::spatial_certificate::certify_spatial(b.degrees, &b.controls, Default::default())
                .determinant_lower_bound_mm3
                .is_some_and(|x| x > 0.)
        }));
    }
}
