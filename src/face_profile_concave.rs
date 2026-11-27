//! Independent inward-ring fallback for simple non-star-shaped profiles.
use crate::face_profile_extrusion::{SectionPatch, lift_sections};

const SCALE: f64 = 1_000_000_000.0;
const GAUGE: i128 = 480;
type H = [i64; 3];
type P = [i64; 2];

fn q(x: f64) -> Result<i64, &'static str> {
    if !x.is_finite() || x.abs() > i64::MAX as f64 / SCALE {
        return Err("concave_profile_range");
    }
    Ok((x * SCALE).round() as i64)
}
fn xy(h: H) -> Result<[f64; 2], &'static str> {
    if h[2] <= 0 {
        return Err("concave_profile_nonpositive_weight");
    }
    Ok([h[0] as f64 / h[2] as f64, h[1] as f64 / h[2] as f64])
}
fn orient(a: P, b: P, c: P) -> i128 {
    (i128::from(b[0]) - i128::from(a[0])) * (i128::from(c[1]) - i128::from(a[1]))
        - (i128::from(b[1]) - i128::from(a[1])) * (i128::from(c[0]) - i128::from(a[0]))
}
fn area2(p: &[P]) -> i128 {
    (0..p.len())
        .map(|i| {
            i128::from(p[i][0]) * i128::from(p[(i + 1) % p.len()][1])
                - i128::from(p[i][1]) * i128::from(p[(i + 1) % p.len()][0])
        })
        .sum()
}
fn on_seg(a: P, b: P, p: P) -> bool {
    orient(a, b, p) == 0 && (0..2).all(|k| p[k] >= a[k].min(b[k]) && p[k] <= a[k].max(b[k]))
}
fn touch(a: P, b: P, c: P, d: P) -> bool {
    let (x, y, z, w) = (
        orient(a, b, c),
        orient(a, b, d),
        orient(c, d, a),
        orient(c, d, b),
    );
    (x.signum() != y.signum() && z.signum() != w.signum())
        || on_seg(a, b, c)
        || on_seg(a, b, d)
        || on_seg(c, d, a)
        || on_seg(c, d, b)
}
fn simple(p: &[P]) -> bool {
    if p.len() < 3 {
        return false;
    }
    for i in 0..p.len() {
        if p[i] == p[(i + 1) % p.len()]
            || orient(p[(i + p.len() - 1) % p.len()], p[i], p[(i + 1) % p.len()]) == 0
        {
            return false;
        }
        for j in i + 1..p.len() {
            if j == i + 1 || (i == 0 && j + 1 == p.len()) {
                continue;
            }
            if touch(p[i], p[(i + 1) % p.len()], p[j], p[(j + 1) % p.len()]) {
                return false;
            }
        }
    }
    true
}
fn inside(p: &[P], q: P) -> bool {
    let (mut inside, mut j) = (false, p.len() - 1);
    for i in 0..p.len() {
        let a = p[i];
        let b = p[j];
        if (a[1] > q[1]) != (b[1] > q[1]) {
            let cross = orient(a, b, q);
            if (cross > 0) == (b[1] > a[1]) {
                inside = !inside
            }
        }
        j = i
    }
    inside
}
fn split(row: &[H]) -> Result<(Vec<H>, Vec<H>), &'static str> {
    let mut a = row.iter().map(|p| p.map(i128::from)).collect::<Vec<_>>();
    let mut l = vec![a[0].map(|v| v * GAUGE)];
    let mut r = vec![a[a.len() - 1].map(|v| v * GAUGE)];
    let mut level = 0;
    while a.len() > 1 {
        level += 1;
        a = a
            .windows(2)
            .map(|p| std::array::from_fn(|k| p[0][k] + p[1][k]))
            .collect();
        let m = GAUGE >> level;
        l.push(a[0].map(|v| v * m));
        r.push(a[a.len() - 1].map(|v| v * m));
    }
    r.reverse();
    let cv = |v: Vec<[i128; 3]>| -> Result<Vec<H>, &'static str> {
        v.into_iter()
            .map(|p| {
                Ok([
                    i64::try_from(p[0]).map_err(|_| "concave_split_overflow")?,
                    i64::try_from(p[1]).map_err(|_| "concave_split_overflow")?,
                    i64::try_from(p[2]).map_err(|_| "concave_split_overflow")?,
                ])
            })
            .collect()
    };
    Ok((cv(l)?, cv(r)?))
}
fn chord(a: H, b: H, d: usize) -> Result<Vec<H>, &'static str> {
    (0..=d)
        .map(|i| {
            let x = d as i128;
            let k = i as i128;
            let f = GAUGE / x;
            let mut h = [0_i64; 3];
            for c in 0..3 {
                h[c] = i64::try_from((i128::from(a[c]) * (x - k) + i128::from(b[c]) * k) * f)
                    .map_err(|_| "concave_chord_overflow")?;
            }
            Ok(h)
        })
        .collect()
}
fn patch(a: Vec<H>, b: Vec<H>) -> Result<SectionPatch, &'static str> {
    if a.len() != b.len() || a.len() < 2 {
        return Err("concave_patch_degree_mismatch");
    }
    let mut controls = Vec::new();
    for i in 0..a.len() {
        controls.push(a[i]);
        controls.push(b[i]);
    }
    if controls.iter().any(|h| h[2] <= 0) {
        return Err("concave_patch_nonpositive_weight");
    }
    Ok(SectionPatch {
        degree: a.len() - 1,
        controls,
    })
}
fn midpoint(a: P, b: P) -> Result<P, &'static str> {
    Ok([
        i64::try_from((i128::from(a[0]) + i128::from(b[0])) / 2)
            .map_err(|_| "concave_midpoint_range")?,
        i64::try_from((i128::from(a[1]) + i128::from(b[1])) / 2)
            .map_err(|_| "concave_midpoint_range")?,
    ])
}
fn point_h(p: P) -> H {
    [p[0], p[1], SCALE as i64]
}
fn triangle_quads(t: [P; 3]) -> Result<[SectionPatch; 3], &'static str> {
    let ab = midpoint(t[0], t[1])?;
    let bc = midpoint(t[1], t[2])?;
    let ca = midpoint(t[2], t[0])?;
    let c = [
        i64::try_from((i128::from(t[0][0]) + i128::from(t[1][0]) + i128::from(t[2][0])) / 3)
            .map_err(|_| "concave_centroid_range")?,
        i64::try_from((i128::from(t[0][1]) + i128::from(t[1][1]) + i128::from(t[2][1])) / 3)
            .map_err(|_| "concave_centroid_range")?,
    ];
    let center = point_h(c);
    let make = |a: H, b: H, d: H| SectionPatch {
        degree: 1,
        controls: vec![a, d, b, center],
    };
    Ok([
        make(point_h(t[0]), point_h(ab), point_h(ca)),
        make(point_h(t[1]), point_h(bc), point_h(ab)),
        make(point_h(t[2]), point_h(ca), point_h(bc)),
    ])
}
fn ears(poly: &[P]) -> Result<Vec<[usize; 3]>, &'static str> {
    if !simple(poly) || area2(poly) <= 0 {
        return Err("concave_inner_polygon_invalid");
    }
    let mut ids = (0..poly.len()).collect::<Vec<_>>();
    let mut out = Vec::new();
    while ids.len() > 3 {
        let mut found = None;
        for k in 0..ids.len() {
            let a = ids[(k + ids.len() - 1) % ids.len()];
            let b = ids[k];
            let c = ids[(k + 1) % ids.len()];
            if orient(poly[a], poly[b], poly[c]) <= 0 {
                continue;
            }
            let has = ids
                .iter()
                .copied()
                .filter(|i| *i != a && *i != b && *i != c)
                .any(|i| {
                    orient(poly[a], poly[b], poly[i]) >= 0
                        && orient(poly[b], poly[c], poly[i]) >= 0
                        && orient(poly[c], poly[a], poly[i]) >= 0
                });
            if !has {
                found = Some((k, [a, b, c]));
                break;
            }
        }
        let (k, t) = found.ok_or("concave_ear_clipping_stalled")?;
        out.push(t);
        ids.remove(k);
    }
    out.push([ids[0], ids[1], ids[2]]);
    let sum: i128 = out
        .iter()
        .map(|t| orient(poly[t[0]], poly[t[1]], poly[t[2]]))
        .sum();
    if sum != area2(poly) {
        return Err("concave_triangulation_area_mismatch");
    }
    Ok(out)
}
fn offset(vertices: &[[f64; 2]], delta: f64) -> Result<Vec<P>, &'static str> {
    let n = vertices.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let p = vertices[i];
        let prev = vertices[(i + n - 1) % n];
        let next = vertices[(i + 1) % n];
        let u = [p[0] - prev[0], p[1] - prev[1]];
        let v = [next[0] - p[0], next[1] - p[1]];
        let lu = u[0].hypot(u[1]);
        let lv = v[0].hypot(v[1]);
        if lu == 0. || lv == 0. {
            return Err("concave_zero_edge");
        }
        let n1 = [-u[1] / lu, u[0] / lu];
        let n2 = [-v[1] / lv, v[0] / lv];
        let mut bis = [n1[0] + n2[0], n1[1] + n2[1]];
        let lb = bis[0].hypot(bis[1]);
        if lb < 1e-10 {
            return Err("concave_offset_bisector_undefined");
        }
        bis = [bis[0] / lb, bis[1] / lb];
        let den = bis[0] * n1[0] + bis[1] * n1[1];
        if den <= 1e-8 {
            return Err("concave_offset_bisector_reversed");
        }
        let d = delta / den;
        out.push([q(p[0] + bis[0] * d)?, q(p[1] + bis[1] * d)?]);
    }
    Ok(out)
}
fn chord_deviation(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let d = [b[0] - a[0], b[1] - a[1]];
    let n = d[0] * d[0] + d[1] * d[1];
    if n == 0.0 {
        return (p[0] - a[0]).hypot(p[1] - a[1]);
    }
    let t = (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / n).clamp(0., 1.);
    (p[0] - a[0] - t * d[0]).hypot(p[1] - a[1] - t * d[1])
}
fn segment_distance(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let d = [b[0] - a[0], b[1] - a[1]];
    let n = d[0] * d[0] + d[1] * d[1];
    if n == 0.0 {
        return (p[0] - a[0]).hypot(p[1] - a[1]);
    }
    let t = (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / n).clamp(0., 1.);
    (p[0] - a[0] - t * d[0]).hypot(p[1] - a[1] - t * d[1])
}
fn control_hulls_disjoint(rows: &[Vec<H>]) -> Result<bool, &'static str> {
    let hulls = rows
        .iter()
        .map(|r| r.iter().map(|h| xy(*h)).collect::<Result<Vec<_>, _>>())
        .collect::<Result<Vec<_>, _>>()?;
    for i in 0..hulls.len() {
        let end = hulls[i][hulls[i].len() - 1];
        let next = &hulls[(i + 1) % hulls.len()];
        let prev_direction = [
            end[0] - hulls[i][hulls[i].len() - 2][0],
            end[1] - hulls[i][hulls[i].len() - 2][1],
        ];
        let next_direction = [next[1][0] - end[0], next[1][1] - end[1]];
        let nd = [
            prev_direction[0] + next_direction[0],
            prev_direction[1] + next_direction[1],
        ];
        let nm = nd[0].hypot(nd[1]);
        if nm <= 1e-14 {
            return Ok(false);
        }
        let normal = [nd[0] / nm, nd[1] / nm];
        for p in &hulls[i][..hulls[i].len() - 1] {
            if (p[0] - end[0]) * normal[0] + (p[1] - end[1]) * normal[1] >= 0.0 {
                return Ok(false);
            }
        }
        for p in &next[1..] {
            if (p[0] - end[0]) * normal[0] + (p[1] - end[1]) * normal[1] <= 0.0 {
                return Ok(false);
            }
        }
    }
    for i in 0..hulls.len() {
        for j in i + 1..hulls.len() {
            if j == i + 1 || (i == 0 && j + 1 == hulls.len()) {
                continue;
            }
            let bounds = |h: &[[f64; 2]]| {
                [
                    [
                        h.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min),
                        h.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max),
                    ],
                    [
                        h.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min),
                        h.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max),
                    ],
                ]
            };
            let a = bounds(&hulls[i]);
            let b = bounds(&hulls[j]);
            if (0..2).all(|k| a[k][0] <= b[k][1] && b[k][0] <= a[k][1]) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn candidate(
    rows: &[Vec<H>],
    inner: &[P],
) -> Result<(Vec<SectionPatch>, Vec<SectionPatch>), &'static str> {
    let triangles = ears(inner)?;
    let mut ring = Vec::with_capacity(2 * rows.len());
    for i in 0..rows.len() {
        let (l, r) = split(&rows[i])?;
        let a = inner[i];
        let b = inner[(i + 1) % inner.len()];
        let m = midpoint(a, b)?;
        ring.push(patch(l, chord(point_h(a), point_h(m), rows[i].len() - 1)?)?);
        ring.push(patch(r, chord(point_h(m), point_h(b), rows[i].len() - 1)?)?);
    }
    let mut core = Vec::with_capacity(triangles.len() * 3);
    for t in triangles {
        let tq = [inner[t[0]], inner[t[1]], inner[t[2]]];
        if orient(tq[0], tq[1], tq[2]) <= 0 {
            return Err("concave_triangle_orientation_invalid");
        }
        core.extend(triangle_quads(tq)?);
    }
    for p in ring.iter().chain(core.iter()) {
        let b = lift_sections(
            std::slice::from_ref(p),
            1.0,
            [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            false,
        )?;
        if b.len() != 1 {
            return Err("concave_chart_lift_failed");
        }
        let cert = crate::spatial_certificate::certify_spatial(
            b[0].degrees,
            &b[0].controls,
            Default::default(),
        );
        if cert.determinant_lower_bound_mm3.is_none_or(|v| v <= 0.) {
            return Err("concave_chart_jacobian_unresolved");
        }
    }
    Ok((ring, core))
}

/// Build a curved inward ring and a midpoint-conforming triangulated core.
/// Returns the ring patches first, then regular quadrilateral core patches.
pub(crate) fn curved_ring_section(
    spans: &[Vec<[f64; 3]>],
    tolerance: f64,
) -> Result<(Vec<SectionPatch>, Vec<SectionPatch>), &'static str> {
    if !(4..=512).contains(&spans.len()) || !tolerance.is_finite() || tolerance <= 0. {
        return Err("concave_profile_input_range");
    }
    let mut rows = spans
        .iter()
        .map(|s| {
            if !(2..=6).contains(&s.len()) {
                return Err("concave_profile_degree_out_of_range");
            }
            let mut out = Vec::new();
            for h in s {
                if h.iter().any(|x| !x.is_finite()) || h[2] <= 0. {
                    return Err("concave_profile_invalid_weight");
                }
                out.push([q(h[0])?, q(h[1])?, q(h[2])?]);
            }
            Ok(out)
        })
        .collect::<Result<Vec<_>, _>>()?;
    for i in 0..rows.len() {
        if rows[i].last() != rows[(i + 1) % rows.len()].first() {
            return Err("concave_profile_spans_must_join_exactly");
        }
    }
    for row in &rows {
        for h in row {
            let p = xy(*h)?;
            if p.iter().any(|x| x.abs() > 1_000_000.) {
                return Err("concave_profile_coordinates_exceed_exact_predicate_range");
            }
        }
    }
    if !control_hulls_disjoint(&rows)? {
        return Err("concave_profile_nonadjacent_control_hulls_overlap");
    }
    let mut vertices = rows
        .iter()
        .map(|r| xy(r[0]))
        .collect::<Result<Vec<_>, _>>()?;
    for row in &rows {
        let controls = row.iter().map(|h| xy(*h)).collect::<Result<Vec<_>, _>>()?;
        let a = controls[0];
        let b = *controls.last().ok_or("concave_empty_span")?;
        let d = [b[0] - a[0], b[1] - a[1]];
        let n = d[0].hypot(d[1]);
        if n <= 1e-14 {
            return Err("concave_span_zero_chord");
        }
        let scalar = controls
            .iter()
            .map(|p| (p[0] - a[0]) * d[0] / n + (p[1] - a[1]) * d[1] / n)
            .collect::<Vec<_>>();
        if scalar.windows(2).any(|w| w[1] < w[0]) || scalar[0] >= scalar[scalar.len() - 1] {
            return Err("concave_span_not_certifiably_injective");
        }
    }
    let mut endpoints = vertices
        .iter()
        .map(|p| Ok([q(p[0])?, q(p[1])?]))
        .collect::<Result<Vec<_>, &'static str>>()?;
    let signed = area2(&endpoints);
    if signed == 0 {
        return Err("concave_profile_zero_area");
    }
    if signed < 0 {
        rows = rows
            .iter()
            .rev()
            .map(|s| s.iter().rev().copied().collect())
            .collect();
        vertices = rows
            .iter()
            .map(|r| xy(r[0]))
            .collect::<Result<Vec<_>, _>>()?;
        endpoints = vertices
            .iter()
            .map(|p| Ok([q(p[0])?, q(p[1])?]))
            .collect::<Result<Vec<_>, &'static str>>()?;
    }
    if !simple(&endpoints) {
        return Err("concave_profile_endpoint_loop_self_intersects");
    }
    let diameter = vertices
        .iter()
        .flat_map(|a| {
            vertices
                .iter()
                .map(move |b| (a[0] - b[0]).hypot(a[1] - b[1]))
        })
        .fold(0., f64::max);
    let mut deviation = 0.0_f64;
    for i in 0..rows.len() {
        let a = vertices[i];
        let b = vertices[(i + 1) % vertices.len()];
        for h in &rows[i] {
            deviation = deviation.max(chord_deviation(xy(*h)?, a, b));
        }
    }
    // The inner polygon must stay farther from the source polyline than the
    // rational control-hull deviation; this keeps the nested boundaries apart.
    let mut delta = (diameter * 1e-3).max((4.0 * deviation).max(tolerance * 16.0));
    for _ in 0..24 {
        if delta <= tolerance * 4.0 {
            break;
        }
        if let Ok(inner) = offset(&vertices, delta) {
            let separated = inner.iter().enumerate().all(|(i, p)| {
                let p = [p[0] as f64 / SCALE, p[1] as f64 / SCALE];
                let clearance = (0..vertices.len())
                    .map(|j| segment_distance(p, vertices[j], vertices[(j + 1) % vertices.len()]))
                    .fold(f64::INFINITY, f64::min);
                clearance > deviation + tolerance * 0.1 && inside(&endpoints, inner[i])
            });
            if simple(&inner) && area2(&inner) > 0 && separated {
                if let Ok(model) = candidate(&rows, &inner) {
                    return Ok(model);
                }
            }
        }
        delta *= 0.5;
    }
    Err("concave_profile_no_regular_inward_offset")
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn concave_c_profile_gets_curved_ring_and_exact_quad_core() {
        let p = [
            [0., 0.],
            [4., 0.],
            [4., 1.],
            [1., 1.],
            [1., 3.],
            [4., 3.],
            [4., 4.],
            [0., 4.],
        ];
        let spans = (0..p.len())
            .map(|i| {
                if i == 6 {
                    vec![
                        [p[i][0], p[i][1], 1.],
                        [
                            (p[i][0] + p[(i + 1) % p.len()][0]) * 0.5,
                            (p[i][1] + p[(i + 1) % p.len()][1]) * 0.5 - 0.01,
                            1.,
                        ],
                        [p[(i + 1) % p.len()][0], p[(i + 1) % p.len()][1], 1.],
                    ]
                } else {
                    vec![
                        [p[i][0], p[i][1], 1.],
                        [p[(i + 1) % p.len()][0], p[(i + 1) % p.len()][1], 1.],
                    ]
                }
            })
            .collect::<Vec<_>>();
        let (ring, core) = curved_ring_section(&spans, 1e-6).unwrap();
        assert_eq!(ring.len(), 16);
        assert!(!core.is_empty());
        assert!(ring.iter().any(|x| x.degree == 2));
        let recipe = crate::face_profile_extrusion::FaceProfileExtrusion {
            spans: spans.clone(),
            bounds: [[-1., 5.], [-1., 5.]],
            thickness: 2.,
            depth: 1.,
            operation: "tool".into(),
            frame: [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            tolerance: 1e-6,
        };
        let model = crate::face_profile_extrusion::build(&recipe).unwrap();
        crate::face_profile_extrusion::audit(&model).unwrap();
        let step = crate::step_boundary::export_step(&model).unwrap();
        assert!(step.contains("MANIFOLD_SOLID_BREP"));
    }
    #[test]
    fn self_crossing_loop_is_rejected() {
        let p = [[0., 0.], [3., 3.], [0., 3.], [3., 0.]];
        let spans = (0..4)
            .map(|i| {
                vec![
                    [p[i][0], p[i][1], 1.],
                    [p[(i + 1) % 4][0], p[(i + 1) % 4][1], 1.],
                ]
            })
            .collect::<Vec<_>>();
        assert!(curved_ring_section(&spans, 1e-6).is_err());
    }
}
