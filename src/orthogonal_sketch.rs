//! Plan disjoint local sweep charts directly from an orthogonal design sketch.
//! Cell occupancy and every seam are checked against the input feature supports.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    constraint_system::Ratio,
    incidence::enumerate_faces,
    sketch_layers::{area2, edges, inside, ratio, valid_name, validate, SketchLayers, SketchLoop},
    spatial_certificate::{certify_spatial, CertificateBudget, CertificateStatus},
};
use std::collections::{BTreeMap, BTreeSet};

pub struct Fillet {
    pub id: String,
    pub vertex: [i32; 2],
    pub radius: i32,
}
pub struct Bore {
    pub id: String,
    pub axis: usize,
    pub center: [i32; 3],
    pub radius: i32,
}
pub struct OrthogonalSketch {
    pub outer: SketchLoop,
    pub height: i32,
    pub fillets: Vec<Fillet>,
    pub bores: Vec<Bore>,
}
type Bounds = [[i32; 2]; 3];
struct Corner<'a> {
    f: &'a Fillet,
    sign: [i32; 2],
    origin: [i32; 2],
    size: i32,
    bounds: Bounds,
}
struct Tube<'a> {
    b: &'a Bore,
    bounds: Bounds,
    guard: i32,
}
pub struct SweepPartition {
    pub blocks: Vec<BracketBlock>,
    pub certificate: BracketCertificate,
    pub twice_constant: i64,
    pub four_pi: i64,
    pub junction: Option<([i32; 3], [usize; 2], i32, i32)>,
    pub planes: [Vec<i32>; 3],
    pub guards: Vec<(String, Bounds)>,
}
fn half(a: i32, b: i32) -> Ratio {
    Ratio::new(i128::from(a) + i128::from(b), 2).expect("dyadic coordinate")
}
fn in_base(j: &OrthogonalSketch, p: [Ratio; 3]) -> Result<bool, &'static str> {
    Ok(p[2].numerator() > 0
        && p[2].checked_sub(ratio(j.height))?.numerator() < 0
        && inside(&j.outer, [p[0], p[1]])?)
}
fn cell_at(j: &OrthogonalSketch, b: &Bounds) -> Result<bool, &'static str> {
    in_base(j, b.map(|t| half(t[0], t[1])))
}
fn enclosed(c: &Bounds, b: &Bounds) -> bool {
    (0..3).all(|k| b[k][0] <= c[k][0] && c[k][1] <= b[k][1])
}
fn overlap(a: &Bounds, b: &Bounds) -> bool {
    (0..3).all(|k| a[k][0].max(b[k][0]) < a[k][1].min(b[k][1]))
}
fn straight(id: String, c: Bounds) -> BracketBlock {
    let mut controls = Vec::with_capacity(8);
    for x in c[0] {
        for y in c[1] {
            for z in c[2] {
                controls.push([i64::from(x), i64::from(y), i64::from(z), 1]);
            }
        }
    }
    BracketBlock {
        id,
        degrees: [1; 3],
        controls,
    }
}
fn corner<'a>(j: &OrthogonalSketch, f: &'a Fillet) -> Result<Corner<'a>, &'static str> {
    if !j.outer.vertices.contains(&f.vertex) || !(1..=1024).contains(&f.radius) {
        return Err("invalid_sketch_fillet");
    }
    let mut notch = Vec::new();
    for x in [-1, 1] {
        for y in [-1, 1] {
            let p = [
                half(f.vertex[0], f.vertex[0] + x),
                half(f.vertex[1], f.vertex[1] + y),
            ];
            if !inside(&j.outer, p)? {
                notch.push([x, y]);
            }
        }
    }
    let &[sign] = notch.as_slice() else {
        return Err("fillet_requires_a_reentrant_right_angle");
    };
    let thickness = std::array::from_fn::<_, 2, _>(|k| {
        j.outer
            .vertices
            .iter()
            .map(|v| (f.vertex[k] - v[k]) * sign[k])
            .filter(|d| *d > 0)
            .min()
            .unwrap_or(0)
    });
    if thickness[0] <= 0 || thickness[0] != thickness[1] {
        return Err("fillet_corner_frame_requires_equal_clearance");
    }
    let origin = std::array::from_fn(|k| f.vertex[k] - thickness[k] * sign[k]);
    let size = thickness[0] + f.radius;
    let bounds = std::array::from_fn(|k| {
        if k == 2 {
            [0, j.height]
        } else {
            let h = origin[k] + sign[k] * size;
            [origin[k].min(h), origin[k].max(h)]
        }
    });
    Ok(Corner {
        f,
        origin,
        sign,
        size,
        bounds,
    })
}
fn tube<'a>(
    j: &OrthogonalSketch,
    b: &'a Bore,
    base: &[BTreeSet<i32>; 3],
) -> Result<Tube<'a>, &'static str> {
    if b.axis >= 3
        || !(1..=1024).contains(&b.radius)
        || b.center.iter().any(|v| !(-4096..=4096).contains(v))
    {
        return Err("invalid_sketch_bore");
    }
    if !in_base(j, b.center.map(ratio))? {
        return Err("bore_seed_outside_material");
    }
    let levels: Vec<_> = base[b.axis].iter().copied().collect();
    let mut runs: Vec<[i32; 2]> = Vec::new();
    for w in levels.windows(2) {
        let mut p = b.center.map(ratio);
        p[b.axis] = half(w[0], w[1]);
        if in_base(j, p)? {
            if let Some(last) = runs.last_mut() {
                if last[1] == w[0] {
                    last[1] = w[1];
                    continue;
                }
            }
            runs.push([w[0], w[1]]);
        }
    }
    let run = runs
        .into_iter()
        .find(|w| w[0] < b.center[b.axis] && b.center[b.axis] < w[1])
        .ok_or("ambiguous_bore_axis_interval")?;
    let guard = b.radius + 1;
    let bounds = std::array::from_fn(|k| {
        if k == b.axis {
            run
        } else {
            [b.center[k] - guard, b.center[k] + guard]
        }
    });
    Ok(Tube { b, bounds, guard })
}
fn transform_fillet(block: &mut BracketBlock, c: &Corner<'_>) -> Result<(), &'static str> {
    for p in &mut block.controls {
        for k in 0..2 {
            p[k] = i64::try_from(
                i128::from(c.origin[k]) * i128::from(p[3])
                    + i128::from(c.sign[k]) * i128::from(p[k]),
            )
            .map_err(|_| "fillet_transform_overflow")?;
        }
    }
    if c.sign[0] * c.sign[1] < 0 {
        let stride = (block.degrees[1] + 1) * (block.degrees[2] + 1);
        for i in 0..stride {
            block.controls.swap(i, block.degrees[0] * stride + i);
        }
    }
    Ok(())
}
fn side_points(block: &BracketBlock, axis: usize, upper: bool) -> Vec<[i64; 4]> {
    let d = block.degrees;
    block
        .controls
        .iter()
        .enumerate()
        .filter_map(|(i, &h)| {
            let uvw = [
                i / ((d[1] + 1) * (d[2] + 1)),
                i / (d[2] + 1) % (d[1] + 1),
                i % (d[2] + 1),
            ];
            (uvw[axis] == usize::from(upper) * d[axis]).then_some(h)
        })
        .collect()
}
fn planar_label(j: &OrthogonalSketch, points: &[[i64; 4]]) -> Option<String> {
    let coord =
        |p: &[i64; 4], k: usize, x: i32| i128::from(p[k]) - i128::from(x) * i128::from(p[3]);
    for (z, name) in [(0, "start"), (j.height, "end")] {
        if points.iter().all(|p| coord(p, 2, z) == 0) {
            return Some(name.into());
        }
    }
    for (a, b) in edges(&j.outer) {
        let k = usize::from(a[1] == b[1]);
        let along = 1 - k;
        if a[k] == b[k]
            && points.iter().all(|p| {
                coord(p, k, a[k]) == 0
                    && coord(p, along, a[along].min(b[along])) >= 0
                    && coord(p, along, a[along].max(b[along])) <= 0
            })
        {
            return Some(format!("outer/{}", j.outer.id));
        }
    }
    None
}
/// Five unisolvent evaluations prove the quartic circle identity for each
/// quadratic homogeneous curve. Identical transverse coefficients on the two
/// axial ends extend the identity through the entire straight sweep.
fn certify_circular_side(
    block: &BracketBlock,
    axis: usize,
    center: [i32; 3],
    radius: i32,
) -> Result<(), &'static str> {
    if block.degrees != [2, 1, 1] {
        return Err("unexpected_circular_feature_degrees");
    }
    let a = (axis + 1) % 3;
    let b = (axis + 2) % 3;
    let mut net = [[0_i128; 3]; 3];
    for (i, out) in net.iter_mut().enumerate() {
        let p = block.controls[i * 4 + 2];
        let q = block.controls[i * 4 + 3];
        if [p[a], p[b], p[3]] != [q[a], q[b], q[3]] {
            return Err("circular_feature_axial_mismatch");
        }
        *out = [
            i128::from(p[a]) - i128::from(center[a]) * i128::from(p[3]),
            i128::from(p[b]) - i128::from(center[b]) * i128::from(p[3]),
            i128::from(radius) * i128::from(p[3]),
        ];
    }
    for j in 0..=4i128 {
        let basis = [(4 - j).pow(2), 2 * j * (4 - j), j.pow(2)];
        let v =
            std::array::from_fn::<_, 3, _>(|k| (0..3).map(|i| basis[i] * net[i][k]).sum::<i128>());
        let sq = |x: i128| {
            x.checked_mul(x)
                .ok_or("circular_feature_certificate_overflow")
        };
        if sq(v[0])?.checked_add(sq(v[1])?) != Some(sq(v[2])?) {
            return Err("circular_feature_polynomial_identity_failed");
        }
    }
    Ok(())
}
/// Build and certify the local-chart partition. Interacting guard regions are explicitly refused.
pub fn build(j: &OrthogonalSketch) -> Result<SweepPartition, &'static str> {
    validate(&SketchLayers {
        outer: j.outer.clone(),
        height: j.height,
        shift: [0, 0],
        cuts: Vec::new(),
        guides: Vec::new(),
        force_triangulation: false,
        operations: None,
    })?;
    if edges(&j.outer).any(|(a, b)| a[0] != b[0] && a[1] != b[1])
        || j.bores.len() > 8
        || j.fillets.len() > 4
    {
        return Err("unsupported_orthogonal_sketch");
    }
    let mut names = BTreeSet::new();
    names.insert(&j.outer.id);
    for id in j
        .fillets
        .iter()
        .map(|v| &v.id)
        .chain(j.bores.iter().map(|v| &v.id))
    {
        if !valid_name(id) || !names.insert(id) {
            return Err("invalid_sketch_feature_identity");
        }
    }
    if let Some(p) = crate::steinmetz_partition::try_build(j)? {
        return Ok(p);
    }
    let mut levels: [BTreeSet<i32>; 3] = std::array::from_fn(|k| {
        if k == 2 {
            BTreeSet::from([0, j.height])
        } else {
            j.outer.vertices.iter().map(|p| p[k]).collect()
        }
    });
    let tubes = j
        .bores
        .iter()
        .map(|b| tube(j, b, &levels))
        .collect::<Result<Vec<_>, _>>()?;
    let corners = j
        .fillets
        .iter()
        .map(|f| corner(j, f))
        .collect::<Result<Vec<_>, _>>()?;
    let guards: Vec<_> = tubes
        .iter()
        .map(|t| (t.b.id.clone(), t.bounds))
        .chain(corners.iter().map(|c| (c.f.id.clone(), c.bounds)))
        .collect();
    for (idx, (_, g)) in guards.iter().enumerate() {
        if guards[..idx].iter().any(|(_, h)| overlap(g, h)) {
            return Err("interacting_feature_guards_require_a_transition_partition");
        }
        for k in 0..3 {
            levels[k].extend(g[k]);
        }
    }
    for t in &tubes {
        for (k, axis) in levels.iter_mut().enumerate() {
            if k != t.b.axis {
                axis.insert(t.b.center[k]);
            }
        }
    }
    let planes = levels.map(|s| s.into_iter().collect::<Vec<_>>());
    for t in &tubes {
        for (k, axis) in planes.iter().enumerate() {
            if k != t.b.axis
                && axis
                    .iter()
                    .any(|&x| x > t.bounds[k][0] && x < t.bounds[k][1] && x != t.b.center[k])
            {
                return Err("bore_guard_needs_additional_trace_subdivision");
            }
        }
    }
    let mut blocks = Vec::new();
    let mut physical = BTreeMap::new();
    for x in planes[0].windows(2) {
        for y in planes[1].windows(2) {
            for z in planes[2].windows(2) {
                let cell = [[x[0], x[1]], [y[0], y[1]], [z[0], z[1]]];
                let material = cell_at(j, &cell)?;
                if tubes.iter().any(|t| enclosed(&cell, &t.bounds)) {
                    if !material {
                        return Err("bore_guard_leaves_sketch_material");
                    }
                    continue;
                }
                if let Some(c) = corners.iter().find(|c| enclosed(&cell, &c.bounds)) {
                    let notch = (0..2)
                        .all(|k| (cell[k][0] + cell[k][1] - 2 * c.f.vertex[k]) * c.sign[k] > 0);
                    if material == notch {
                        return Err("fillet_guard_is_not_a_simple_corner");
                    }
                    continue;
                }
                if material {
                    blocks.push(straight(format!("cell_{}_{}_{}", x[0], y[0], z[0]), cell));
                }
            }
        }
    }
    let mut twice_constant = i128::from(area2(&j.outer)) * i128::from(j.height);
    let mut four_pi = 0i128;
    for t in tubes {
        for w in planes[t.b.axis]
            .windows(2)
            .filter(|w| w[0] >= t.bounds[t.b.axis][0] && w[1] <= t.bounds[t.b.axis][1])
        {
            let id = format!("bore_{}_{}_{}", t.b.id, w[0], w[1]);
            for block in crate::sketch_bore_chart::bore_charts(
                t.b.axis,
                t.b.center,
                [w[0], w[1]],
                t.guard,
                t.b.radius,
                &id,
            )? {
                certify_circular_side(&block, t.b.axis, t.b.center, t.b.radius)?;
                physical.insert((block.id.clone(), 1, true), format!("bore/{}/wall", t.b.id));
                blocks.push(block);
            }
        }
        four_pi -= 4
            * i128::from(t.b.radius).pow(2)
            * i128::from(t.bounds[t.b.axis][1] - t.bounds[t.b.axis][0]);
    }
    for c in corners {
        let vol = i128::from(c.f.radius).pow(2) * i128::from(j.height);
        twice_constant += 2 * vol;
        four_pi -= vol;
        for w in planes[2].windows(2) {
            let id = format!("fillet_{}_{}_{}", c.f.id, w[0], w[1]);
            for mut block in crate::sketch_fillet_chart::fillet_charts(
                [0, 0],
                c.size,
                c.f.radius,
                [w[0], w[1]],
                &id,
            )? {
                transform_fillet(&mut block, &c)?;
                let center = [
                    c.f.vertex[0] + c.sign[0] * c.f.radius,
                    c.f.vertex[1] + c.sign[1] * c.f.radius,
                    0,
                ];
                certify_circular_side(&block, 2, center, c.f.radius)?;
                physical.insert((block.id.clone(), 1, true), format!("fillet/{}", c.f.id));
                blocks.push(block);
            }
        }
    }
    blocks.sort_by(|a, b| a.id.cmp(&b.id));
    if blocks.is_empty() || blocks.len() > 10000 {
        return Err("orthogonal_sketch_block_budget");
    }
    let incidence = enumerate_faces(&blocks)?;
    let index: BTreeMap<_, _> = blocks.iter().map(|b| (&b.id, b)).collect();
    let mut tags = Vec::new();
    for side in &incidence.exterior {
        let b = index[&side.block];
        let tag = physical
            .remove(&(side.block.clone(), side.axis, side.upper))
            .or_else(|| planar_label(j, &side_points(b, side.axis, side.upper)))
            .ok_or("uncoupled_orthogonal_sketch_interface")?;
        tags.push((side.clone(), tag));
    }
    if !physical.is_empty() {
        return Err("physical_feature_face_buried");
    }
    let mut det = f64::INFINITY;
    for b in &blocks {
        let c = certify_spatial(b.degrees, &b.controls, CertificateBudget::default());
        if c.status != CertificateStatus::Valid {
            return Err("orthogonal_sketch_spatial_certificate_incomplete");
        }
        det = det.min(
            c.determinant_lower_bound_mm3
                .ok_or("orthogonal_sketch_missing_bound")?,
        );
    }
    Ok(SweepPartition {
        blocks,
        certificate: BracketCertificate {
            incidence,
            determinant_lower_bound_mm3: det,
            design_faces: tags,
            radial_error_upper_mm: 0.,
        },
        twice_constant: i64::try_from(twice_constant).map_err(|_| "sketch_volume_range")?,
        junction: None,
        four_pi: i64::try_from(four_pi).map_err(|_| "sketch_volume_range")?,
        planes,
        guards,
    })
}
