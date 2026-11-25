//! Independent layer-wise certificate for an exact translating bore.
//!
//! This audit consumes only serialized controls. It never invokes the bore
//! constructor. Each layer is certified as a positive, conforming planar
//! complex whose uncancelled, oriented edges form the expected plate/pocket/
//! circle boundary. Exact boundary coverage and positive Jacobians give degree
//! one, hence a partition without overlap or holes for this arrangement.

use crate::bracket::BracketBlock;
use crate::bracket_audit::BracketCertificate;
use crate::incidence::{enumerate_faces, BlockSide, Incidence};
use crate::planar::{certify_planar_refined, PlanarError};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::convert::TryFrom;

const CENTER_MIN_UM: i128 = -5_000;
const CENTER_MAX_UM: i128 = 14_000;
const RADIUS_UM: i128 = 2_000;
const SCALE_UM_PER_MM: i128 = 1_000;
const OUTER_HALF_UM: i128 = 10_000;
const POCKET_X_UM: [i128; 2] = [-6_000, -2_000];
const POCKET_Y_UM: [i128; 2] = [-4_000, 4_000];

type Point = [i128; 3];
type Trace = Vec<Point>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Layer {
    Lower,
    Upper,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum BoundaryFeature {
    XMin,
    XMax,
    YMin,
    YMax,
    Bore,
    PocketXMin,
    PocketXMax,
    PocketYMin,
    PocketYMax,
}

impl BoundaryFeature {
    fn tag(self) -> &'static str {
        match self {
            Self::XMin => "x_min",
            Self::XMax => "x_max",
            Self::YMin => "y_min",
            Self::YMax => "y_max",
            Self::Bore => "bore_wall",
            Self::PocketXMin => "pocket_x_min",
            Self::PocketXMax => "pocket_x_max",
            Self::PocketYMin => "pocket_y_min",
            Self::PocketYMax => "pocket_y_max",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Rat {
    n: i128,
    d: i128,
}

impl Rat {
    fn new(mut n: i128, mut d: i128) -> Result<Self, &'static str> {
        if d == 0 {
            return Err("travelling_bore_zero_denominator");
        }
        if d < 0 {
            n = n.checked_neg().ok_or("travelling_bore_integer_overflow")?;
            d = d.checked_neg().ok_or("travelling_bore_integer_overflow")?;
        }
        let g = gcd(n, d)?;
        Ok(Self { n: n / g, d: d / g })
    }

    fn integer(n: i128) -> Self {
        Self { n, d: 1 }
    }
}

#[derive(Clone)]
struct BoundaryEdge {
    block_id: String,
    axis: usize,
    upper: bool,
    trace: Trace,
    layer: Layer,
    feature: Option<BoundaryFeature>,
}

fn gcd(mut a: i128, mut b: i128) -> Result<i128, &'static str> {
    a = a.checked_abs().ok_or("travelling_bore_integer_overflow")?;
    b = b.checked_abs().ok_or("travelling_bore_integer_overflow")?;
    while b != 0 {
        (a, b) = (b, a % b);
    }
    Ok(a)
}

fn checked_mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("travelling_bore_integer_overflow")
}

fn checked_add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("travelling_bore_integer_overflow")
}

fn checked_sub(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_sub(b).ok_or("travelling_bore_integer_overflow")
}

fn normalize_trace(mut trace: Trace) -> Result<Trace, &'static str> {
    let divisor = trace
        .iter()
        .flatten()
        .try_fold(0_i128, |g, value| gcd(g, *value))?;
    if divisor <= 0 || trace.iter().any(|p| p[2] <= 0) {
        return Err("travelling_bore_invalid_trace_weight");
    }
    for point in &mut trace {
        for value in point {
            *value /= divisor;
        }
    }
    Ok(trace)
}

fn normalize_point(point: Point) -> Result<Point, &'static str> {
    let divisor = point.iter().try_fold(0_i128, |g, value| gcd(g, *value))?;
    if divisor <= 0 || point[2] <= 0 {
        return Err("travelling_bore_invalid_endpoint");
    }
    Ok(point.map(|value| value / divisor))
}

fn validate_extrusion(block: &BracketBlock) -> Result<Layer, &'static str> {
    if block.degrees != [2, 2, 1] || block.controls.len() != 18 {
        return Err("travelling_bore_control_layout");
    }
    let mut selected = None;
    for i in 0..3 {
        for j in 0..3 {
            let p = block.controls[(i * 3 + j) * 2];
            let q = block.controls[(i * 3 + j) * 2 + 1];
            let w = i128::from(p[3]);
            if w <= 0 || q[..2] != p[..2] || q[3] != p[3] {
                return Err("travelling_bore_extrusion_identity");
            }
            let lower = i128::from(p[2]);
            let upper = i128::from(q[2]);
            let current = if lower == 0 && upper == 4 * w {
                Layer::Lower
            } else if lower == 4 * w && upper == 8 * w {
                Layer::Upper
            } else {
                return Err("travelling_bore_layer_identity");
            };
            if selected.replace(current).is_some_and(|old| old != current) {
                return Err("travelling_bore_mixed_layer");
            }
        }
    }
    selected.ok_or("travelling_bore_empty_block")
}

fn profile(block: &BracketBlock) -> Vec<[i64; 3]> {
    (0..9)
        .map(|p| {
            let c = block.controls[2 * p];
            [c[0], c[1], c[3]]
        })
        .collect()
}

fn side_trace(block: &BracketBlock, axis: usize, upper: bool) -> Result<Trace, &'static str> {
    let mut points = Vec::with_capacity(3);
    match (axis, upper) {
        (0, false) => {
            for j in (0..3).rev() {
                let p = block.controls[j * 2];
                points.push([i128::from(p[0]), i128::from(p[1]), i128::from(p[3])]);
            }
        }
        (0, true) => {
            for j in 0..3 {
                let p = block.controls[(2 * 3 + j) * 2];
                points.push([i128::from(p[0]), i128::from(p[1]), i128::from(p[3])]);
            }
        }
        (1, false) => {
            for i in 0..3 {
                let p = block.controls[(i * 3) * 2];
                points.push([i128::from(p[0]), i128::from(p[1]), i128::from(p[3])]);
            }
        }
        (1, true) => {
            for i in (0..3).rev() {
                let p = block.controls[(i * 3 + 2) * 2];
                points.push([i128::from(p[0]), i128::from(p[1]), i128::from(p[3])]);
            }
        }
        _ => return Err("travelling_bore_invalid_side"),
    }
    normalize_trace(points)
}

fn paired_boundaries(blocks: &[BracketBlock]) -> Result<Vec<BoundaryEdge>, &'static str> {
    let mut by_layer = [Vec::<BoundaryEdge>::new(), Vec::<BoundaryEdge>::new()];
    for block in blocks {
        let l = validate_extrusion(block)?;
        let target = if l == Layer::Lower {
            &mut by_layer[0]
        } else {
            &mut by_layer[1]
        };
        for axis in 0..2 {
            for upper in [false, true] {
                target.push(BoundaryEdge {
                    block_id: block.id.clone(),
                    axis,
                    upper,
                    trace: side_trace(block, axis, upper)?,
                    layer: l,
                    feature: None,
                });
            }
        }
    }
    let mut boundary = Vec::new();
    for layer_edges in by_layer {
        let mut groups = BTreeMap::<Trace, Vec<usize>>::new();
        for (index, edge) in layer_edges.iter().enumerate() {
            groups.entry(edge.trace.clone()).or_default().push(index);
        }
        let mut consumed = BTreeSet::new();
        for (trace, indices) in &groups {
            if consumed.contains(trace) {
                continue;
            }
            let mut reverse = trace.clone();
            reverse.reverse();
            if &reverse == trace {
                return Err("travelling_bore_collapsed_edge");
            }
            match groups.get(&reverse) {
                Some(reverse_indices) => {
                    if indices.len() != 1 || reverse_indices.len() != 1 {
                        return Err("travelling_bore_ambiguous_edge_pair");
                    }
                    consumed.insert(trace.clone());
                    consumed.insert(reverse);
                }
                None => {
                    if indices.len() != 1 {
                        return Err("travelling_bore_duplicate_boundary_edge");
                    }
                    consumed.insert(trace.clone());
                    boundary.push(layer_edges[indices[0]].clone());
                }
            }
        }
    }
    Ok(boundary)
}

fn support_value(points: &[Point], axis: usize, um: i128) -> bool {
    points.iter().all(|p| {
        checked_mul(SCALE_UM_PER_MM, p[axis])
            .and_then(|scaled| checked_sub(scaled, checked_mul(um, p[2])?))
            == Ok(0)
    })
}

fn cross_delta(a: i128, b: i128, c: i128, d: i128) -> Result<i128, &'static str> {
    checked_sub(checked_mul(a, b)?, checked_mul(c, d)?)
}

fn monotonic_sign(
    trace: &[Point],
    component: usize,
    center_um: Option<i128>,
) -> Result<i8, &'static str> {
    let coord = |p: Point| -> Result<i128, &'static str> {
        let value =
            i128::from(i64::try_from(p[component]).map_err(|_| "travelling_bore_trace_range")?);
        if let Some(center) = center_um {
            checked_sub(
                checked_mul(SCALE_UM_PER_MM, value)?,
                checked_mul(center, p[2])?,
            )
        } else {
            checked_mul(SCALE_UM_PER_MM, value)
        }
    };
    let h = trace
        .iter()
        .copied()
        .map(coord)
        .collect::<Result<Vec<_>, _>>()?;
    let [h0, h1, h2] = [h[0], h[1], h[2]];
    let [w0, w1, w2] = [trace[0][2], trace[1][2], trace[2][2]];
    let n0 = checked_mul(2, cross_delta(h1, w0, h0, w1)?)?;
    let n1 = cross_delta(h2, w0, h0, w2)?;
    let n2 = checked_mul(2, cross_delta(h2, w1, h1, w2)?)?;
    let values = [n0, n1, n2];
    let positive = values.iter().all(|v| *v >= 0) && values.iter().any(|v| *v > 0);
    let negative = values.iter().all(|v| *v <= 0) && values.iter().any(|v| *v < 0);
    match (positive, negative) {
        (true, false) => Ok(1),
        (false, true) => Ok(-1),
        _ => Err("travelling_bore_nonmonotone_boundary"),
    }
}

fn compare_nonnegative(mut a: Rat, mut b: Rat) -> Result<Ordering, &'static str> {
    let mut reverse = false;
    loop {
        let aq = a.n / a.d;
        let bq = b.n / b.d;
        if aq != bq {
            let order = aq.cmp(&bq);
            return Ok(if reverse { order.reverse() } else { order });
        }
        let ar = a.n % a.d;
        let br = b.n % b.d;
        match (ar == 0, br == 0) {
            (true, true) => return Ok(Ordering::Equal),
            (true, false) => {
                let order = Ordering::Less;
                return Ok(if reverse { order.reverse() } else { order });
            }
            (false, true) => {
                let order = Ordering::Greater;
                return Ok(if reverse { order.reverse() } else { order });
            }
            (false, false) => {
                a = Rat::new(a.d, ar)?;
                b = Rat::new(b.d, br)?;
                reverse = !reverse;
            }
        }
    }
}

fn rat_cmp(a: Rat, b: Rat) -> Result<Ordering, &'static str> {
    match (a.n.signum(), b.n.signum()) {
        (x, y) if x != y => Ok(x.cmp(&y)),
        (0, 0) => Ok(Ordering::Equal),
        (-1, -1) => Ok(compare_nonnegative(
            Rat::new(
                a.n.checked_neg()
                    .ok_or("travelling_bore_integer_overflow")?,
                a.d,
            )?,
            Rat::new(
                b.n.checked_neg()
                    .ok_or("travelling_bore_integer_overflow")?,
                b.d,
            )?,
        )?
        .reverse()),
        _ => compare_nonnegative(a, b),
    }
}

fn sort_intervals(values: &mut [(Rat, Rat)]) -> Result<(), &'static str> {
    for i in 1..values.len() {
        let mut j = i;
        while j > 0 {
            let order =
                rat_cmp(values[j].0, values[j - 1].0)?.then(rat_cmp(values[j].1, values[j - 1].1)?);
            if order != Ordering::Less {
                break;
            }
            values.swap(j, j - 1);
            j -= 1;
        }
    }
    Ok(())
}

fn merge_intervals(mut values: Vec<(Rat, Rat)>) -> Result<Vec<(Rat, Rat)>, &'static str> {
    sort_intervals(&mut values)?;
    let mut merged = Vec::<(Rat, Rat)>::new();
    for interval in values {
        if rat_cmp(interval.0, interval.1)? != Ordering::Less {
            return Err("travelling_bore_zero_length_boundary");
        }
        if let Some(last) = merged.last_mut() {
            match rat_cmp(interval.0, last.1)? {
                Ordering::Less => return Err("travelling_bore_boundary_overlap"),
                Ordering::Equal => last.1 = interval.1,
                Ordering::Greater => merged.push(interval),
            }
        } else {
            merged.push(interval);
        }
    }
    Ok(merged)
}

fn sqrt_exact(n: i128) -> Result<Option<i128>, &'static str> {
    if n < 0 {
        return Ok(None);
    }
    let mut x = n;
    let mut y = x.checked_add(1).ok_or("travelling_bore_integer_overflow")? / 2;
    while y < x {
        x = y;
        y = x
            .checked_add(n / x)
            .ok_or("travelling_bore_integer_overflow")?
            / 2;
    }
    Ok((x.checked_mul(x).ok_or("travelling_bore_integer_overflow")? == n).then_some(x))
}

fn circle_identity(trace: &[Point], center_um: i128) -> Result<bool, &'static str> {
    for p in 0..=4_i128 {
        let q = 4 - p;
        let sigma = [
            checked_mul(q, q)?,
            checked_mul(2, checked_mul(p, q)?)?,
            checked_mul(p, p)?,
        ];
        let mut x = 0_i128;
        let mut y = 0_i128;
        let mut w = 0_i128;
        for i in 0..3 {
            x = checked_add(
                x,
                checked_mul(
                    sigma[i],
                    checked_sub(
                        checked_mul(SCALE_UM_PER_MM, trace[i][0])?,
                        checked_mul(center_um, trace[i][2])?,
                    )?,
                )?,
            )?;
            y = checked_add(
                y,
                checked_mul(sigma[i], checked_mul(SCALE_UM_PER_MM, trace[i][1])?)?,
            )?;
            w = checked_add(w, checked_mul(sigma[i], trace[i][2])?)?;
        }
        let lhs = checked_add(checked_mul(x, x)?, checked_mul(y, y)?)?;
        let rhs = checked_mul(checked_mul(2_000, w)?, checked_mul(2_000, w)?)?;
        if lhs != rhs {
            return Ok(false);
        }
    }
    Ok(true)
}

fn circle_signs(trace: &[Point], center_um: i128) -> Result<(i8, i8), &'static str> {
    let mut sx = 0_i8;
    let mut sy = 0_i8;
    for p in trace {
        let x = checked_sub(
            checked_mul(SCALE_UM_PER_MM, p[0])?,
            checked_mul(center_um, p[2])?,
        )?;
        let y = checked_mul(SCALE_UM_PER_MM, p[1])?;
        if x != 0 {
            let sign = x.signum() as i8;
            if sx != 0 && sx != sign {
                return Err("travelling_bore_arc_crosses_quadrant");
            }
            sx = sign;
        }
        if y != 0 {
            let sign = y.signum() as i8;
            if sy != 0 && sy != sign {
                return Err("travelling_bore_arc_crosses_quadrant");
            }
            sy = sign;
        }
    }
    if sx == 0 || sy == 0 {
        return Err("travelling_bore_arc_degenerate_quadrant");
    }
    Ok((sx, sy))
}

fn line_feature(trace: &[Point], layer: Layer) -> Result<Option<BoundaryFeature>, &'static str> {
    for x in [-OUTER_HALF_UM, OUTER_HALF_UM] {
        if support_value(trace, 0, x) {
            monotonic_sign(trace, 1, None)?;
            return Ok(Some(if x < 0 {
                BoundaryFeature::XMin
            } else {
                BoundaryFeature::XMax
            }));
        }
    }
    for y in [-OUTER_HALF_UM, OUTER_HALF_UM] {
        if support_value(trace, 1, y) {
            monotonic_sign(trace, 0, None)?;
            return Ok(Some(if y < 0 {
                BoundaryFeature::YMin
            } else {
                BoundaryFeature::YMax
            }));
        }
    }
    if layer == Layer::Upper {
        for x in POCKET_X_UM {
            if support_value(trace, 0, x) {
                monotonic_sign(trace, 1, None)?;
                return Ok(Some(if x == POCKET_X_UM[0] {
                    BoundaryFeature::PocketXMin
                } else {
                    BoundaryFeature::PocketXMax
                }));
            }
        }
        for y in POCKET_Y_UM {
            if support_value(trace, 1, y) {
                monotonic_sign(trace, 0, None)?;
                return Ok(Some(if y < 0 {
                    BoundaryFeature::PocketYMin
                } else {
                    BoundaryFeature::PocketYMax
                }));
            }
        }
    }
    Ok(None)
}

fn endpoint_coordinate(point: Point, axis: usize) -> Result<Rat, &'static str> {
    Rat::new(checked_mul(SCALE_UM_PER_MM, point[axis])?, point[2])
}

fn endpoint_circle_t(point: Point, center_um: i128) -> Result<Rat, &'static str> {
    Rat::new(
        checked_sub(
            checked_mul(SCALE_UM_PER_MM, point[0])?,
            checked_mul(center_um, point[2])?,
        )?
        .checked_abs()
        .ok_or("travelling_bore_integer_overflow")?,
        point[2],
    )
}

fn expected_circle_intervals(center_um: i128, sx: i8, layer: Layer) -> Vec<(i128, i128)> {
    let mut cuts = vec![0, RADIUS_UM];
    let outer_t = (OUTER_HALF_UM - center_um).abs();
    if outer_t > 0 && outer_t < RADIUS_UM {
        cuts.push(outer_t);
    }
    if layer == Layer::Upper {
        for wall in POCKET_X_UM {
            let t = (wall - center_um).abs();
            if t > 0 && t < RADIUS_UM {
                cuts.push(t);
            }
        }
    }
    cuts.sort_unstable();
    cuts.dedup();
    let pieces = cuts
        .windows(2)
        .filter_map(|pair| {
            let (a, b) = (pair[0], pair[1]);
            let world_twice = 2 * center_um + i128::from(sx) * (a + b);
            let inside_pocket = layer == Layer::Upper
                && world_twice > 2 * POCKET_X_UM[0]
                && world_twice < 2 * POCKET_X_UM[1];
            let inside_plate = world_twice < 2 * OUTER_HALF_UM;
            (!inside_pocket && inside_plate).then_some((a, b))
        })
        .collect::<Vec<_>>();
    let mut merged = Vec::<(i128, i128)>::new();
    for (lo, hi) in pieces {
        if let Some(last) = merged.last_mut() {
            if last.1 == lo {
                last.1 = hi;
                continue;
            }
        }
        merged.push((lo, hi));
    }
    merged
}

fn validate_arc_coverage(
    edges: &[BoundaryEdge],
    center_um: i128,
    layer: Layer,
) -> Result<(), &'static str> {
    for sx in [-1_i8, 1] {
        for sy in [-1_i8, 1] {
            let mut actual = Vec::<(Rat, Rat)>::new();
            for edge in edges {
                let trace = &edge.trace;
                if edge.feature != Some(BoundaryFeature::Bore) {
                    continue;
                }
                if !circle_identity(trace, center_um)? {
                    continue;
                }
                let (edge_sx, edge_sy) = circle_signs(trace, center_um)?;
                if (edge_sx, edge_sy) != (sx, sy) {
                    continue;
                }
                monotonic_sign(trace, 0, Some(center_um))?;
                let start = normalize_point(trace[0])?;
                let end = normalize_point(trace[2])?;
                let t0 = endpoint_circle_t(start, center_um)?;
                let t1 = endpoint_circle_t(end, center_um)?;
                if rat_cmp(t0, Rat::integer(RADIUS_UM))? == Ordering::Greater
                    || rat_cmp(t1, Rat::integer(RADIUS_UM))? == Ordering::Greater
                {
                    return Err("travelling_bore_arc_outside_radius");
                }
                let start_t = endpoint_circle_t(start, center_um)?;
                let end_t = endpoint_circle_t(end, center_um)?;
                let (lo, hi) = if rat_cmp(start_t, end_t)? == Ordering::Less {
                    (start_t, end_t)
                } else {
                    (end_t, start_t)
                };
                if rat_cmp(lo, hi)? == Ordering::Equal {
                    return Err("travelling_bore_zero_arc_span");
                }
                let cross = cross_delta(
                    checked_sub(
                        checked_mul(SCALE_UM_PER_MM, start[0])?,
                        checked_mul(center_um, start[2])?,
                    )?,
                    checked_mul(SCALE_UM_PER_MM, end[1])?,
                    checked_mul(SCALE_UM_PER_MM, start[1])?,
                    checked_sub(
                        checked_mul(SCALE_UM_PER_MM, end[0])?,
                        checked_mul(center_um, end[2])?,
                    )?,
                )?;
                if cross >= 0 {
                    return Err("travelling_bore_arc_orientation");
                }
                actual.push((lo, hi));
            }
            actual = merge_intervals(actual)?;
            let expected = expected_circle_intervals(center_um, sx, layer);
            if actual.len() != expected.len() {
                return Err(match (layer, sx, sy) {
                    (Layer::Lower, -1, -1) => "travelling_bore_arc_coverage_l0",
                    (Layer::Lower, -1, 1) => "travelling_bore_arc_coverage_l1",
                    (Layer::Lower, 1, -1) => "travelling_bore_arc_coverage_r0",
                    (Layer::Lower, 1, 1) => "travelling_bore_arc_coverage_r1",
                    (Layer::Upper, -1, -1) => "travelling_bore_arc_coverage_u_l0",
                    (Layer::Upper, -1, 1) => "travelling_bore_arc_coverage_u_l1",
                    (Layer::Upper, 1, -1) => "travelling_bore_arc_coverage_u_r0",
                    _ => "travelling_bore_arc_coverage_u_r1",
                });
            }
            for ((lo, hi), (expected_lo, expected_hi)) in actual.iter().zip(expected) {
                if rat_cmp(*lo, Rat::integer(expected_lo))? != Ordering::Equal
                    || rat_cmp(*hi, Rat::integer(expected_hi))? != Ordering::Equal
                {
                    return Err("travelling_bore_arc_endpoint_coverage");
                }
            }
        }
    }
    Ok(())
}

fn exact_y_cut(center_um: i128, wall_um: i128) -> Result<i128, &'static str> {
    let delta = (wall_um - center_um).abs();
    if delta >= RADIUS_UM {
        return Ok(0);
    }
    let radicand = checked_sub(
        checked_mul(RADIUS_UM, RADIUS_UM)?,
        checked_mul(delta, delta)?,
    )?;
    sqrt_exact(radicand)?.ok_or("travelling_bore_nonrepresentable_intersection")
}

fn expected_line_intervals(
    feature: BoundaryFeature,
    center_um: i128,
) -> Result<Vec<(Rat, Rat)>, &'static str> {
    let i = |x| Rat::integer(x);
    Ok(match feature {
        BoundaryFeature::XMax => {
            let cut = exact_y_cut(center_um, OUTER_HALF_UM)?;
            if cut == 0 {
                vec![(i(-OUTER_HALF_UM), i(OUTER_HALF_UM))]
            } else {
                vec![(i(-OUTER_HALF_UM), i(-cut)), (i(cut), i(OUTER_HALF_UM))]
            }
        }
        BoundaryFeature::XMin | BoundaryFeature::YMin | BoundaryFeature::YMax => {
            vec![(i(-OUTER_HALF_UM), i(OUTER_HALF_UM))]
        }
        BoundaryFeature::PocketYMin | BoundaryFeature::PocketYMax => {
            vec![(i(POCKET_X_UM[0]), i(POCKET_X_UM[1]))]
        }
        BoundaryFeature::PocketXMin | BoundaryFeature::PocketXMax => {
            let x = if feature == BoundaryFeature::PocketXMin {
                POCKET_X_UM[0]
            } else {
                POCKET_X_UM[1]
            };
            let cut = exact_y_cut(center_um, x)?;
            if cut == 0 {
                vec![(i(-4_000), i(4_000))]
            } else {
                vec![(i(-4_000), i(-cut)), (i(cut), i(4_000))]
            }
        }
        BoundaryFeature::Bore => return Err("travelling_bore_internal_feature_error"),
    })
}

fn edge_line_interval(edge: &BoundaryEdge) -> Result<(Rat, Rat), &'static str> {
    let (constant_axis, _) = if matches!(
        edge.feature,
        Some(
            BoundaryFeature::XMin
                | BoundaryFeature::XMax
                | BoundaryFeature::PocketXMin
                | BoundaryFeature::PocketXMax
        )
    ) {
        (0, 1)
    } else {
        (1, 0)
    };
    let variable_axis = 1 - constant_axis;
    let start = endpoint_coordinate(normalize_point(edge.trace[0])?, variable_axis)?;
    let end = endpoint_coordinate(normalize_point(edge.trace[2])?, variable_axis)?;
    if rat_cmp(start, end)? == Ordering::Less {
        Ok((start, end))
    } else {
        Ok((end, start))
    }
}

fn validate_line_coverage(
    edges: &[BoundaryEdge],
    center_um: i128,
    layer: Layer,
) -> Result<(), &'static str> {
    let features: &[BoundaryFeature] = if layer == Layer::Lower {
        &[
            BoundaryFeature::XMin,
            BoundaryFeature::XMax,
            BoundaryFeature::YMin,
            BoundaryFeature::YMax,
        ]
    } else {
        &[
            BoundaryFeature::XMin,
            BoundaryFeature::XMax,
            BoundaryFeature::YMin,
            BoundaryFeature::YMax,
            BoundaryFeature::PocketXMin,
            BoundaryFeature::PocketXMax,
            BoundaryFeature::PocketYMin,
            BoundaryFeature::PocketYMax,
        ]
    };
    for feature in features {
        let mut actual = edges
            .iter()
            .filter(|edge| edge.feature == Some(*feature))
            .map(edge_line_interval)
            .collect::<Result<Vec<_>, _>>()?;
        actual = merge_intervals(actual)?;
        let expected = expected_line_intervals(*feature, center_um)?;
        if actual.len() != expected.len() {
            return Err("travelling_bore_line_coverage");
        }
        for (actual, expected) in actual.iter().zip(expected) {
            if rat_cmp(actual.0, expected.0)? != Ordering::Equal
                || rat_cmp(actual.1, expected.1)? != Ordering::Equal
            {
                return Err("travelling_bore_line_coverage");
            }
        }
    }
    Ok(())
}

fn edge_orientation_is_valid(edge: &BoundaryEdge) -> Result<bool, &'static str> {
    let start = normalize_point(edge.trace[0])?;
    let end = normalize_point(edge.trace[2])?;
    let (axis, expected_sign) = match edge.feature.ok_or("travelling_bore_unclassified_edge")? {
        BoundaryFeature::XMin => (1, -1),
        BoundaryFeature::XMax => (1, 1),
        BoundaryFeature::YMin => (0, 1),
        BoundaryFeature::YMax => (0, -1),
        BoundaryFeature::PocketXMin => (1, 1),
        BoundaryFeature::PocketXMax => (1, -1),
        BoundaryFeature::PocketYMin => (0, -1),
        BoundaryFeature::PocketYMax => (0, 1),
        BoundaryFeature::Bore => return Ok(true),
    };
    let a = endpoint_coordinate(start, axis)?;
    let b = endpoint_coordinate(end, axis)?;
    let order = rat_cmp(a, b)?;
    Ok((expected_sign < 0 && order == Ordering::Greater)
        || (expected_sign > 0 && order == Ordering::Less))
}

fn boundary_cycles(edges: &[BoundaryEdge]) -> Result<usize, &'static str> {
    let mut degree = BTreeMap::<Point, (usize, usize)>::new();
    let mut next = BTreeMap::<Point, Point>::new();
    for edge in edges {
        let start = normalize_point(edge.trace[0])?;
        let end = normalize_point(edge.trace[2])?;
        degree.entry(start).or_default().1 += 1;
        degree.entry(end).or_default().0 += 1;
        if next.insert(start, end).is_some() {
            return Err("travelling_bore_boundary_branch");
        }
    }
    if degree
        .values()
        .any(|(incoming, outgoing)| *incoming != 1 || *outgoing != 1)
    {
        return Err("travelling_bore_open_boundary");
    }
    let mut visited = BTreeSet::new();
    let mut cycles = 0;
    for start in next.keys() {
        if visited.contains(start) {
            continue;
        }
        cycles += 1;
        let mut cursor = *start;
        loop {
            if !visited.insert(cursor) {
                if cursor != *start {
                    return Err("travelling_bore_boundary_cycle_overlap");
                }
                break;
            }
            cursor = *next.get(&cursor).ok_or("travelling_bore_open_boundary")?;
        }
    }
    Ok(cycles)
}

fn classify_boundaries(edges: &mut [BoundaryEdge], center_um: i128) -> Result<(), &'static str> {
    for edge in edges {
        if let Some(feature) = line_feature(&edge.trace, edge.layer)? {
            edge.feature = Some(feature);
            if !edge_orientation_is_valid(edge)? {
                return Err("travelling_bore_boundary_orientation");
            }
        } else if circle_identity(&edge.trace, center_um)? {
            let (sx, sy) = circle_signs(&edge.trace, center_um)?;
            let _ = (sx, sy);
            monotonic_sign(&edge.trace, 0, Some(center_um))?;
            let start = normalize_point(edge.trace[0])?;
            let end = normalize_point(edge.trace[2])?;
            let cross = cross_delta(
                checked_sub(
                    checked_mul(SCALE_UM_PER_MM, start[0])?,
                    checked_mul(center_um, start[2])?,
                )?,
                checked_mul(SCALE_UM_PER_MM, end[1])?,
                checked_mul(SCALE_UM_PER_MM, start[1])?,
                checked_sub(
                    checked_mul(SCALE_UM_PER_MM, end[0])?,
                    checked_mul(center_um, end[2])?,
                )?,
            )?;
            if cross >= 0 {
                return Err("travelling_bore_arc_orientation");
            }
            edge.feature = Some(BoundaryFeature::Bore);
        } else {
            return Err("travelling_bore_unexpected_boundary_curve");
        }
    }
    Ok(())
}

fn certify_layer(
    layer_blocks: &[&BracketBlock],
    center_um: i128,
    layer: Layer,
) -> Result<(Vec<BoundaryEdge>, f64, usize), &'static str> {
    let mut det = f64::INFINITY;
    for block in layer_blocks {
        let lower =
            certify_planar_refined([2, 2], &profile(block)).map_err(|error| match error {
                PlanarError::Invalid => "travelling_bore_nonpositive_jacobian",
                PlanarError::Indeterminate => "travelling_bore_indeterminate_jacobian",
            })?;
        if !lower.is_finite() || lower <= 0.0 {
            return Err("travelling_bore_nonpositive_jacobian");
        }
        det = det.min(lower * 4.0);
    }
    let all = layer_blocks
        .iter()
        .map(|b| (*b).clone())
        .collect::<Vec<_>>();
    let mut boundary = paired_boundaries(&all)?;
    boundary.retain(|edge| edge.layer == layer);
    classify_boundaries(&mut boundary, center_um)?;
    validate_line_coverage(&boundary, center_um, layer)?;
    validate_arc_coverage(&boundary, center_um, layer)?;
    let cycles = boundary_cycles(&boundary)?;
    let expected_cycles = match layer {
        Layer::Lower if center_um < 8_000 => 2,
        Layer::Lower => 1,
        Layer::Upper if !(0..=8_000).contains(&center_um) => 2,
        Layer::Upper => 3,
    };
    if cycles != expected_cycles {
        return Err("travelling_bore_boundary_component_count");
    }
    Ok((boundary, det, cycles))
}

fn face_controls(block: &BracketBlock, side: &BlockSide) -> Vec<[i64; 4]> {
    let fixed = if side.upper { 2 } else { 0 };
    match side.axis {
        0 => (0..3)
            .flat_map(|j| (0..2).map(move |k| block.controls[(fixed * 3 + j) * 2 + k]))
            .collect(),
        1 => (0..3)
            .flat_map(|i| (0..2).map(move |k| block.controls[(i * 3 + fixed) * 2 + k]))
            .collect(),
        2 => {
            let k = usize::from(side.upper);
            (0..3)
                .flat_map(|i| (0..3).map(move |j| block.controls[(i * 3 + j) * 2 + k]))
                .collect()
        }
        _ => Vec::new(),
    }
}

fn floor_face_is_inside_pocket(block: &BracketBlock, side: &BlockSide) -> bool {
    face_controls(block, side).iter().all(|p| {
        let w = i128::from(p[3]);
        let x = i128::from(p[0]);
        let y = i128::from(p[1]);
        x >= POCKET_X_UM[0] * w / SCALE_UM_PER_MM
            && x <= POCKET_X_UM[1] * w / SCALE_UM_PER_MM
            && y >= POCKET_Y_UM[0] * w / SCALE_UM_PER_MM
            && y <= POCKET_Y_UM[1] * w / SCALE_UM_PER_MM
    })
}

fn interface_is_outside_pocket(block: &BracketBlock, side: &BlockSide) -> bool {
    let controls = face_controls(block, side);
    [
        controls
            .iter()
            .all(|p| SCALE_UM_PER_MM * i128::from(p[0]) <= POCKET_X_UM[0] * i128::from(p[3])),
        controls
            .iter()
            .all(|p| SCALE_UM_PER_MM * i128::from(p[0]) >= POCKET_X_UM[1] * i128::from(p[3])),
        controls
            .iter()
            .all(|p| SCALE_UM_PER_MM * i128::from(p[1]) <= POCKET_Y_UM[0] * i128::from(p[3])),
        controls
            .iter()
            .all(|p| SCALE_UM_PER_MM * i128::from(p[1]) >= POCKET_Y_UM[1] * i128::from(p[3])),
    ]
    .into_iter()
    .any(|outside_halfspace| outside_halfspace)
}

fn validate_interfaces(incidence: &Incidence, blocks: &[BracketBlock]) -> Result<(), &'static str> {
    let by_id = blocks
        .iter()
        .map(|block| (block.id.as_str(), block))
        .collect::<BTreeMap<_, _>>();
    for interface in &incidence.interfaces {
        let first = by_id
            .get(interface.first.block.as_str())
            .ok_or("travelling_bore_unknown_interface_block")?;
        let second = by_id
            .get(interface.second.block.as_str())
            .ok_or("travelling_bore_unknown_interface_block")?;
        if interface.first.axis == 2 || interface.second.axis == 2 {
            if interface.first.axis != 2 || interface.second.axis != 2 {
                return Err("travelling_bore_interface_axis_mismatch");
            }
            let first_layer = validate_extrusion(first)?;
            let second_layer = validate_extrusion(second)?;
            let (lower_side, upper_side, upper_block) = match (first_layer, second_layer) {
                (Layer::Lower, Layer::Upper)
                    if interface.first.upper && !interface.second.upper =>
                {
                    (&interface.first, &interface.second, second)
                }
                (Layer::Upper, Layer::Lower)
                    if !interface.first.upper && interface.second.upper =>
                {
                    (&interface.second, &interface.first, first)
                }
                _ => return Err("travelling_bore_layer_cap_interface_mismatch"),
            };
            if lower_side.axis != 2
                || upper_side.axis != 2
                || !interface_is_outside_pocket(upper_block, upper_side)
            {
                return Err("travelling_bore_cap_interface_inside_pocket");
            }
        } else if validate_extrusion(first)? != validate_extrusion(second)? {
            return Err("travelling_bore_lateral_interface_layer_mismatch");
        }
    }
    Ok(())
}

fn side_tag(
    side: &BlockSide,
    boundaries: &[BoundaryEdge],
    by_id: &BTreeMap<&str, &BracketBlock>,
) -> Result<String, &'static str> {
    let block = by_id
        .get(side.block.as_str())
        .ok_or("travelling_bore_unknown_face")?;
    if side.axis == 2 {
        let l = validate_extrusion(block)?;
        return match (l, side.upper) {
            (Layer::Lower, false) => Ok("bottom".into()),
            (Layer::Upper, true) => Ok("top".into()),
            (Layer::Lower, true) => {
                if !floor_face_is_inside_pocket(block, side) {
                    return Err("travelling_bore_floor_outside_pocket");
                }
                Ok("pocket_floor".into())
            }
            (Layer::Upper, false) => Err("travelling_bore_missing_layer_interface"),
        };
    }
    boundaries
        .iter()
        .find(|edge| {
            edge.block_id == side.block && edge.axis == side.axis && edge.upper == side.upper
        })
        .and_then(|edge| edge.feature)
        .map(|feature| feature.tag().to_owned())
        .ok_or("travelling_bore_unlabelled_boundary_face")
}

/// Independently certify a serialized moving-bore partition at `cx_um`.
pub fn audit_travelling_bore(
    cx_um: i32,
    blocks: &[BracketBlock],
) -> Result<BracketCertificate, &'static str> {
    let center = i128::from(cx_um);
    if !(CENTER_MIN_UM..=CENTER_MAX_UM).contains(&center) {
        return Err("travelling_bore_center_out_of_range");
    }
    if blocks.is_empty() {
        return Err("travelling_bore_empty_complex");
    }
    for wall in [POCKET_X_UM[0], POCKET_X_UM[1], OUTER_HALF_UM] {
        if (wall - center).abs() == RADIUS_UM {
            return Err("travelling_bore_tangent_event");
        }
        let delta = (wall - center).abs();
        if delta < RADIUS_UM {
            let radicand = checked_sub(
                checked_mul(RADIUS_UM, RADIUS_UM)?,
                checked_mul(delta, delta)?,
            )?;
            if sqrt_exact(radicand)?.is_none() {
                return Err("travelling_bore_nonrepresentable_intersection");
            }
        }
    }
    let mut ids = BTreeSet::new();
    let mut lower = Vec::<&BracketBlock>::new();
    let mut upper = Vec::<&BracketBlock>::new();
    for block in blocks {
        if !ids.insert(block.id.as_str()) {
            return Err("travelling_bore_duplicate_block_id");
        }
        match validate_extrusion(block)? {
            Layer::Lower => lower.push(block),
            Layer::Upper => upper.push(block),
        }
    }
    if lower.is_empty() || upper.is_empty() {
        return Err("travelling_bore_missing_layer");
    }
    let (lower_boundary, lower_det, _) = certify_layer(&lower, center, Layer::Lower)?;
    let (upper_boundary, upper_det, _) = certify_layer(&upper, center, Layer::Upper)?;
    let mut boundaries = lower_boundary;
    boundaries.extend(upper_boundary);
    let incidence = enumerate_faces(blocks)?;
    validate_interfaces(&incidence, blocks)?;
    let by_id = blocks
        .iter()
        .map(|block| (block.id.as_str(), block))
        .collect::<BTreeMap<_, _>>();
    let mut design_faces = Vec::with_capacity(incidence.exterior.len());
    for side in &incidence.exterior {
        design_faces.push((side.clone(), side_tag(side, &boundaries, &by_id)?));
    }
    design_faces.sort_by(|a, b| a.0.cmp(&b.0));
    let determinant_lower_bound_mm3 = lower_det.min(upper_det);
    if !determinant_lower_bound_mm3.is_finite() || determinant_lower_bound_mm3 <= 0.0 {
        return Err("travelling_bore_determinant_certificate");
    }
    Ok(BracketCertificate {
        incidence,
        determinant_lower_bound_mm3,
        radial_error_upper_mm: 0.0,
        design_faces,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::travelling_bore::travelling_bore_blocks;

    #[test]
    fn certifies_separated_overlap_open_and_outside_regimes() {
        for center in [3_000, -800, -3_200, -4_800, 10_000, 11_920, 13_000] {
            let blocks = travelling_bore_blocks(center).unwrap();
            let certificate = audit_travelling_bore(center, &blocks)
                .unwrap_or_else(|error| panic!("center {center}: {error}"));
            assert!(certificate.determinant_lower_bound_mm3 > 0.0);
            assert_eq!(certificate.radial_error_upper_mm, 0.0);
            assert!(certificate
                .design_faces
                .iter()
                .any(|(_, tag)| tag == "x_max"));
            assert!(certificate
                .design_faces
                .iter()
                .any(|(_, tag)| tag == "pocket_floor"));
        }
    }

    #[test]
    fn accepts_independent_positive_projective_gauges() {
        let mut blocks = travelling_bore_blocks(-3_200).unwrap();
        for (index, block) in blocks.iter_mut().enumerate() {
            let factor = if index % 2 == 0 { 2_i64 } else { 3_i64 };
            for control in &mut block.controls {
                for coordinate in control {
                    *coordinate = coordinate.checked_mul(factor).unwrap();
                }
            }
        }
        let certificate = audit_travelling_bore(-3_200, &blocks).unwrap();
        assert!(certificate.determinant_lower_bound_mm3 > 0.0);
    }

    #[test]
    fn rejects_tangencies_irrational_cuts_mutations_and_missing_blocks() {
        for center in [-4_000, 0, 8_000, 12_000] {
            let built = travelling_bore_blocks(center);
            assert!(built.is_err());
            assert_eq!(
                audit_travelling_bore(center, &[]),
                Err("travelling_bore_empty_complex")
            );
        }
        let irrational = travelling_bore_blocks(-3_000);
        assert!(irrational.is_err());
        let blocks = travelling_bore_blocks(3_000).unwrap();
        assert_eq!(
            audit_travelling_bore(-3_000, &blocks),
            Err("travelling_bore_nonrepresentable_intersection")
        );
        let mut mutated = blocks.clone();
        mutated[0].controls[0][0] += 1;
        assert!(audit_travelling_bore(3_000, &mutated).is_err());
        let mut missing = blocks;
        missing.pop();
        assert!(audit_travelling_bore(3_000, &missing).is_err());
    }
}
