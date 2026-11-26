//! Conforming tensor-product charts for a periodic Z bore crossed by a
//! rectangular X bore. The outer cage is a 12-sector, positively oriented
//! partition of the rectangular stock boundary.

use crate::{
    bracket::BracketBlock,
    constraint_system::Ratio as R,
    incidence::{enumerate_faces, BlockSide},
    json::Json,
    nurbs_curve::BezierSpan,
    sketch_partition_plan::{crossings, SourceCurve},
    spatial_certificate::{certify_spatial, CertificateStatus},
};
use std::collections::{BTreeMap, BTreeSet};

type Error = &'static str;
type H = [R; 3];
type HNet = Vec<H>;
type Construction = (Vec<BracketBlock>, Vec<(BlockSide, String)>);

const ANCHOR_COUNT: usize = 12;
const LATERAL: [usize; 4] = [3, 4, 9, 10];

#[derive(Clone)]
struct ProfilePiece {
    controls: HNet,
    start: [R; 2],
    end: [R; 2],
    duration: R,
    source_span: usize,
}

#[derive(Clone, Copy)]
struct Anchor {
    point: [R; 2],
    outer: [R; 2],
    direction: Option<usize>,
}

#[derive(Clone)]
struct Tag {
    sector: usize,
    layer: usize,
    side_min: [R; 2],
    side_max: [R; 2],
}

fn int(x: i64) -> R {
    R::new(i128::from(x), 1).expect("integer ratio")
}

fn div(a: R, b: R) -> Result<R, Error> {
    a.checked_mul(R::new(b.denominator(), b.numerator())?)
}

fn point(h: H) -> Result<[R; 2], Error> {
    if h[2].numerator() <= 0 {
        return Err("periodic profile has a nonpositive denominator");
    }
    Ok([div(h[0], h[2])?, div(h[1], h[2])?])
}

fn compare(a: R, b: R) -> std::cmp::Ordering {
    a.cmp_value(&b)
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn lcm(a: i128, b: i128) -> Result<i128, Error> {
    let g = i128::try_from(gcd(a.unsigned_abs(), b.unsigned_abs()))
        .map_err(|_| "homogeneous scale exceeds i128")?;
    (a / g)
        .checked_mul(b)
        .ok_or("homogeneous scale exceeds i128")
}

/// Return exact Bézier spans over the active source domain. Each span keeps
/// the extraction gauge supplied by the shared blossom implementation.
pub(crate) fn profile_spans(curve: &Json) -> Result<Vec<BezierSpan>, Error> {
    let source = SourceCurve::parse(curve)?;
    crate::nurbs_curve::extract_active_bezier(source.degree, &source.knots, &source.controls)
}

fn ratio_net(span: &BezierSpan) -> HNet {
    span.controls.iter().map(|p| p.map(int)).collect()
}

fn scale_net(net: &mut HNet, scale: R) -> Result<(), Error> {
    for p in net {
        for v in p {
            *v = v.checked_mul(scale)?;
        }
    }
    Ok(())
}

fn reverse_net(net: &HNet) -> HNet {
    net.iter().copied().rev().collect()
}

fn split_net(net: &HNet, t: R) -> Result<(HNet, HNet), Error> {
    if t == R::zero() || t == R::one() {
        return Err("profile subdivision must be interior");
    }
    let mut row = net.clone();
    let mut left = vec![row[0]];
    let mut right = vec![*row.last().ok_or("empty profile span")?];
    while row.len() > 1 {
        let mut next_row = Vec::with_capacity(row.len() - 1);
        for pair in row.windows(2) {
            let mut q = [R::zero(); 3];
            for d in 0..3 {
                q[d] =
                    pair[0][d].checked_add(t.checked_mul(pair[1][d].checked_sub(pair[0][d])?)?)?;
            }
            next_row.push(q);
        }
        row = next_row;
        left.push(row[0]);
        right.push(*row.last().unwrap());
    }
    right.reverse();
    Ok((left, right))
}

fn direction(p: [R; 2], center: [R; 2]) -> Result<usize, Error> {
    let dx = p[0].checked_sub(center[0])?;
    let dy = p[1].checked_sub(center[1])?;
    if dx == R::zero() && dy == R::zero() {
        return Err("periodic profile anchor equals its center");
    }
    let abs = |v: R| -> Result<R, Error> {
        if v.numerator() < 0 {
            R::zero().checked_sub(v)
        } else {
            Ok(v)
        }
    };
    if abs(dx)? == abs(dy)? {
        Ok(match (dx.numerator().signum(), dy.numerator().signum()) {
            (1, 1) => 0,
            (-1, 1) => 2,
            (-1, -1) => 4,
            (1, -1) => 6,
            _ => return Err("invalid diagonal profile direction"),
        })
    } else if dx == R::zero() {
        Ok(if dy.numerator() > 0 { 1 } else { 5 })
    } else if dy == R::zero() {
        Ok(if dx.numerator() < 0 { 3 } else { 7 })
    } else {
        Err("periodic profile anchor is not on a cage ray")
    }
}

fn cross(a: [R; 2], b: [R; 2]) -> Result<R, Error> {
    a[0].checked_mul(b[1])?.checked_sub(a[1].checked_mul(b[0])?)
}

fn verify_polar_monotonicity(net: &[H], center: [R; 2]) -> Result<(), Error> {
    let points = net
        .iter()
        .map(|&h| point(h))
        .collect::<Result<Vec<_>, _>>()?;
    let vectors = points
        .iter()
        .map(|p| Ok([p[0].checked_sub(center[0])?, p[1].checked_sub(center[1])?]))
        .collect::<Result<Vec<[R; 2]>, Error>>()?;
    for i in 0..vectors.len() {
        for j in i + 1..vectors.len() {
            if cross(vectors[i], vectors[j])?.numerator() <= 0 {
                return Err("periodic profile span lacks a strict single-sector polar certificate");
            }
        }
    }
    Ok(())
}

struct Curve {
    bounds: [i64; 4],
    center: [R; 2],
    y_low: R,
    y_high: R,
    z_cut: [i64; 2],
    height: i64,
}

fn parse_curve(recipe: &Json) -> Result<(Curve, Vec<BezierSpan>), Error> {
    recipe.keys(&[
        "family",
        "bounds_xyz_mm",
        "profile",
        "rect_y_mm",
        "rect_z_mm",
    ])?;
    if recipe.field("family")?.text()? != "g2_cross_plate_v1" {
        return Err("unsupported cross-plate family");
    }
    let bounds = recipe.field("bounds_xyz_mm")?.items()?;
    let [xb, yb, zb] = bounds else {
        return Err("expected XYZ stock bounds");
    };
    let pair = |value: &Json| -> Result<[i64; 2], Error> {
        let [a, b] = value.items()? else {
            return Err("expected two stock bounds");
        };
        Ok([a.number()?, b.number()?])
    };
    let [xmin, xmax] = pair(xb)?;
    let [ymin, ymax] = pair(yb)?;
    let [zmin, zmax] = pair(zb)?;
    if [xmin, ymin, zmin] != [0; 3] || xmin >= xmax || ymin >= ymax || zmin >= zmax {
        return Err("cross-plate stock bounds must be positive and start at the origin");
    }
    let rect_y = recipe.field("rect_y_mm")?.items()?;
    let [yl, yh] = rect_y else {
        return Err("rectangle y bounds must be rational pairs");
    };
    let rational = |value: &Json| -> Result<R, Error> {
        let [n, d] = value.items()? else {
            return Err("rectangle ordinate must be a numerator/denominator pair");
        };
        let (n, d) = (n.number::<i64>()?, d.number::<i64>()?);
        if d <= 0 {
            return Err("rectangle denominator must be positive");
        }
        R::new(i128::from(n), i128::from(d))
    };
    let y_low = rational(yl)?;
    let y_high = rational(yh)?;
    if compare(y_low, int(ymin)) != std::cmp::Ordering::Greater
        || compare(y_high, int(ymax)) != std::cmp::Ordering::Less
        || compare(y_low, y_high) != std::cmp::Ordering::Less
    {
        return Err("rectangle y section must be strictly inside stock");
    }
    let z_cut = pair(recipe.field("rect_z_mm")?)?;
    if z_cut[0] <= zmin || z_cut[0] >= z_cut[1] || z_cut[1] >= zmax {
        return Err("rectangle z interval must be internal to the plate");
    }
    let source = SourceCurve::parse(recipe.field("profile")?)?;
    if source.degree != 3 {
        return Err("periodic profile must be cubic");
    }
    let spans = profile_spans(recipe.field("profile")?)?;
    if spans.len() < 8 {
        return Err("periodic profile needs eight active spans");
    }
    for p in source.controls {
        if p[2] <= 0
            || compare(R::new(i128::from(p[0]), i128::from(p[2]))?, int(xmin))
                != std::cmp::Ordering::Greater
            || compare(R::new(i128::from(p[0]), i128::from(p[2]))?, int(xmax))
                != std::cmp::Ordering::Less
            || compare(R::new(i128::from(p[1]), i128::from(p[2]))?, int(ymin))
                != std::cmp::Ordering::Greater
            || compare(R::new(i128::from(p[1]), i128::from(p[2]))?, int(ymax))
                != std::cmp::Ordering::Less
        {
            return Err("periodic profile control hull must be strictly inside stock");
        }
    }
    if spans[0].controls[0] != *spans.last().unwrap().controls.last().unwrap() {
        // Extraction spans may carry independent integer gauges, so compare
        // physical closure below instead of rejecting unequal homogeneous nets.
        let a = point(spans[0].controls[0].map(int))?;
        let b = point(spans.last().unwrap().controls.last().unwrap().map(int))?;
        if a != b {
            return Err("periodic profile does not close geometrically");
        }
    }
    let endpoints: Vec<_> = spans
        .iter()
        .flat_map(|s| [s.controls[0], *s.controls.last().unwrap()])
        .map(|p| point(p.map(int)))
        .collect::<Result<_, _>>()?;
    let extrema = |axis: usize| -> Result<(R, R), Error> {
        let mut values = endpoints.iter().map(|p| p[axis]);
        let first = values.next().ok_or("profile has no endpoints")?;
        Ok(values.fold((first, first), |(lo, hi), v| {
            (
                if compare(v, lo).is_lt() { v } else { lo },
                if compare(v, hi).is_gt() { v } else { hi },
            )
        }))
    };
    let (xmin_p, xmax_p) = extrema(0)?;
    let (ymin_p, ymax_p) = extrema(1)?;
    let center = [
        xmin_p.checked_add(xmax_p)?.checked_mul(R::new(1, 2)?)?,
        ymin_p.checked_add(ymax_p)?.checked_mul(R::new(1, 2)?)?,
    ];
    Ok((
        Curve {
            bounds: [xmin, xmax, ymin, ymax],
            center,
            y_low,
            y_high,
            z_cut,
            height: zmax,
        },
        spans,
    ))
}

fn gcd_point(net: &[H]) -> Result<(), Error> {
    if net.iter().any(|p| p[2].numerator() <= 0) {
        return Err("profile span has nonpositive weight");
    }
    Ok(())
}

fn canonical_spans(
    spans: &[BezierSpan],
    center: [R; 2],
) -> Result<Vec<(HNet, usize, bool)>, Error> {
    let mut nets: Vec<_> = spans.iter().map(ratio_net).collect();
    for i in 1..nets.len() {
        let p = *nets[i - 1].last().unwrap();
        let q = nets[i][0];
        scale_net(&mut nets[i], div(p[2], q[2])?)?;
        if nets[i][0] != p {
            return Err("profile knot traces disagree after gauge alignment");
        }
    }
    let first = nets[0][0];
    let last = *nets.last().unwrap().last().unwrap();
    if point(first)? != point(last)? || first != last {
        return Err("periodic profile endpoints are not coincident in one homogeneous gauge");
    }
    let unique: Vec<_> = nets
        .iter()
        .flat_map(|n| [n[0], *n.last().unwrap()])
        .map(point)
        .collect::<Result<_, _>>()?;
    let anchors: BTreeMap<usize, [R; 2]> = unique
        .into_iter()
        .filter_map(|p| direction(p, center).ok().map(|d| (d, p)))
        .collect();
    if anchors.len() != 8 {
        return Err("profile must have eight distinct cardinal/diagonal cage anchors");
    }

    let start = anchors[&0]; // northeast diagonal, independent of source index.
    let incident: Vec<_> = nets
        .iter()
        .enumerate()
        .filter_map(|(i, net)| {
            let a = point(net[0]).ok()?;
            let b = point(*net.last()?).ok()?;
            if a == start {
                Some((i, false, b))
            } else if b == start {
                Some((i, true, a))
            } else {
                None
            }
        })
        .collect();
    if incident.len() != 2 {
        return Err("periodic profile does not have a unique start anchor");
    }
    let mut chosen = None;
    for option in &incident {
        let end = option.2;
        let a = [
            start[0].checked_sub(center[0])?,
            start[1].checked_sub(center[1])?,
        ];
        let b = [
            end[0].checked_sub(center[0])?,
            end[1].checked_sub(center[1])?,
        ];
        if cross(a, b)?.numerator() > 0 {
            chosen = Some(option);
            break;
        }
    }
    let chosen = chosen.ok_or("profile orientation is not counterclockwise from northeast")?;
    let mut used = BTreeSet::new();
    let mut ordered = Vec::with_capacity(nets.len());
    let mut current = start;
    let mut first_piece = true;
    let mut candidate = (chosen.0, chosen.1);
    while ordered.len() < nets.len() {
        let (index, reverse) = candidate;
        if !used.insert(index) {
            return Err("profile span graph contains a cycle before closure");
        }
        let mut net = nets[index].clone();
        if reverse {
            net = reverse_net(&net);
        }
        if point(net[0])? != current {
            return Err("profile span graph has a disconnected endpoint");
        }
        gcd_point(&net)?;
        verify_polar_monotonicity(&net, center)?;
        current = point(*net.last().unwrap())?;
        ordered.push((net, index, reverse));
        first_piece = false;
        if current == start {
            break;
        }
        let next: Vec<_> = nets
            .iter()
            .enumerate()
            .filter(|(i, _)| !used.contains(i))
            .filter_map(|(i, net)| {
                let a = point(net[0]).ok()?;
                let b = point(*net.last()?).ok()?;
                if a == current {
                    Some((i, false))
                } else if b == current {
                    Some((i, true))
                } else {
                    None
                }
            })
            .collect();
        if next.len() != 1 {
            return Err("profile spans do not form one simple periodic chain");
        }
        candidate = next[0];
    }
    let _ = first_piece;
    if ordered.len() != nets.len() || current != start {
        return Err("profile chain did not close after one revolution");
    }
    Ok(ordered)
}

fn endpoint_directions(nets: &[HNet], center: [R; 2]) -> Result<BTreeMap<usize, [R; 2]>, Error> {
    let mut anchors = BTreeMap::new();
    for net in nets {
        for h in [net[0], *net.last().unwrap()] {
            let p = point(h)?;
            if let Ok(d) = direction(p, center) {
                if anchors.insert(d, p).is_some_and(|old| old != p) {
                    return Err("duplicate profile point on one cage ray");
                }
            }
        }
    }
    if anchors.len() != 8 {
        return Err("profile must expose one endpoint on each cage ray");
    }
    Ok(anchors)
}

fn anchor_outer(direction: usize, bounds: [i64; 4], center: [R; 2]) -> Result<[R; 2], Error> {
    let [xmin, xmax, ymin, ymax] = bounds;
    let [center_x, center_y] = center;
    Ok(match direction {
        0 => [int(xmax), int(ymax)],
        1 => [center_x, int(ymax)],
        2 => [int(xmin), int(ymax)],
        3 => [int(xmin), center_y],
        4 => [int(xmin), int(ymin)],
        5 => [center_x, int(ymin)],
        6 => [int(xmax), int(ymin)],
        7 => [int(xmax), center_y],
        _ => unreachable!(),
    })
}

fn add_ratio_point(
    map: &mut BTreeMap<[R; 2], Anchor>,
    point: [R; 2],
    outer: [R; 2],
    d: Option<usize>,
) -> Result<(), Error> {
    if map
        .insert(
            point,
            Anchor {
                point,
                outer,
                direction: d,
            },
        )
        .is_some()
    {
        return Err("duplicate geometric sector anchor");
    }
    Ok(())
}

fn split_piece_at(net: HNet, cuts: &[R]) -> Result<Vec<(HNet, R, R)>, Error> {
    let mut parts = vec![(net, R::zero(), R::one())];
    let mut previous = R::zero();
    let mut result = Vec::new();
    for cut in cuts {
        let local = div(cut.checked_sub(previous)?, R::one().checked_sub(previous)?)?;
        let (rest, lo, hi) = parts.pop().ok_or("missing profile subdivision")?;
        let (left, right) = split_net(&rest, local)?;
        let mid = lo.checked_add(hi.checked_sub(lo)?.checked_mul(local)?)?;
        result.push((left, lo, mid));
        parts.push((right, mid, hi));
        previous = *cut;
    }
    result.extend(parts);
    Ok(result)
}

fn line_controls(a: H, b: H, degree: usize) -> Result<HNet, Error> {
    if degree != 3 {
        return Err("outer cage curve requires cubic elevation");
    }
    let mut net = vec![a, b];
    while net.len() < degree + 1 {
        let n = net.len() - 1;
        let m = n + 1;
        let mut elevated = Vec::with_capacity(m + 1);
        elevated.push(net[0]);
        for i in 1..=n {
            let left = R::new(i as i128, m as i128)?;
            let right = R::one().checked_sub(left)?;
            let mut control = [R::zero(); 3];
            for (d, value) in control.iter_mut().enumerate() {
                *value = net[i - 1][d]
                    .checked_mul(left)?
                    .checked_add(net[i][d].checked_mul(right)?)?;
            }
            elevated.push(control);
        }
        elevated.push(*net.last().unwrap());
        net = elevated;
    }
    Ok(net)
}

fn outer_net(inner: &HNet, endpoints: [[R; 2]; 2], lateral_x: Option<i64>) -> Result<HNet, Error> {
    if let Some(x) = lateral_x {
        return inner
            .iter()
            .map(|p| {
                let scale = R::new(3, 1)?;
                Ok([
                    int(x).checked_mul(p[2])?.checked_mul(scale)?,
                    p[1].checked_mul(scale)?,
                    p[2].checked_mul(scale)?,
                ])
            })
            .collect();
    }
    let scale = R::new(3, 1)?;
    let a = [
        endpoints[0][0]
            .checked_mul(inner[0][2])?
            .checked_mul(scale)?,
        endpoints[0][1]
            .checked_mul(inner[0][2])?
            .checked_mul(scale)?,
        inner[0][2].checked_mul(scale)?,
    ];
    let b = [
        endpoints[1][0]
            .checked_mul(inner.last().unwrap()[2])?
            .checked_mul(scale)?,
        endpoints[1][1]
            .checked_mul(inner.last().unwrap()[2])?
            .checked_mul(scale)?,
        inner.last().unwrap()[2].checked_mul(scale)?,
    ];
    line_controls(a, b, inner.len() - 1)
}

fn block_controls(outer: &HNet, inner: &HNet, z0: i64, z1: i64) -> Result<Vec<[i64; 4]>, Error> {
    if outer.len() != 4 || inner.len() != 4 {
        return Err("cross-plate sector must have cubic profile degree");
    }
    let mut net = Vec::with_capacity(16);
    for i in 0..4 {
        for (v, p) in [outer[i], inner[i]].into_iter().enumerate() {
            let _ = v;
            for z in [z0, z1] {
                net.push([p[0], p[1], int(z).checked_mul(p[2])?, p[2]]);
            }
        }
    }
    let mut scale = 1_i128;
    for p in &net {
        for c in p {
            scale = lcm(scale, c.denominator())?;
        }
    }
    let mut integers = Vec::with_capacity(net.len());
    let mut common = 0_u128;
    for p in net {
        let mut q = [0_i128; 4];
        for d in 0..4 {
            q[d] = p[d]
                .numerator()
                .checked_mul(scale / p[d].denominator())
                .ok_or("cross-plate homogeneous coefficient exceeds i128")?;
            common = gcd(common, q[d].unsigned_abs());
        }
        integers.push(q);
    }
    let common =
        i128::try_from(common).map_err(|_| "cross-plate homogeneous scale exceeds i128")?;
    integers
        .into_iter()
        .map(|p| {
            p.map(|v| i64::try_from(v / common).map_err(|_| "cross-plate control exceeds i64"))
                .into_iter()
                .collect::<Result<Vec<_>, _>>()?
                .try_into()
                .map_err(|_| "invalid cross-plate control")
        })
        .collect()
}

fn boundary_side(points: [[R; 2]; 2], bounds: [i64; 4]) -> Result<usize, Error> {
    let [xmin, xmax, ymin, ymax] = bounds;
    if points[0][0] == int(xmin) && points[1][0] == int(xmin) {
        Ok(0)
    } else if points[0][0] == int(xmax) && points[1][0] == int(xmax) {
        Ok(1)
    } else if points[0][1] == int(ymin) && points[1][1] == int(ymin) {
        Ok(2)
    } else if points[0][1] == int(ymax) && points[1][1] == int(ymax) {
        Ok(3)
    } else {
        Err("outer anchor pair is not on a stock rectangle side")
    }
}

fn transverse_face_y(block: &BracketBlock, side: &BlockSide) -> Result<R, Error> {
    let u = if side.upper { 3 } else { 0 };
    let mut common = None;
    for v in 0..=1 {
        for z in 0..=1 {
            let control = block.controls[(u * 2 + v) * 2 + z];
            let y = R::new(i128::from(control[1]), i128::from(control[3]))?;
            if common.is_some_and(|value| value != y) {
                return Err("transverse exposed face is not a constant-y sketch edge");
            }
            common = Some(y);
        }
    }
    common.ok_or("transverse exposed face has no controls")
}

/// Construct the full stock minus the two intersecting through bores.
pub fn construct(recipe: &Json) -> Result<Construction, Error> {
    let (parameters, _) = crate::sketch_partition_plan::plan(recipe)?;
    construct_partition(&parameters)
}

/// Construct from a normalized plan emitted by the shared sketch planner.
pub(crate) fn construct_partition(parameters: &Json) -> Result<Construction, Error> {
    let (curve, spans) = parse_curve(parameters)?;
    let nets = canonical_spans(&spans, curve.center)?;
    let net_values: Vec<_> = nets.iter().map(|x| x.0.clone()).collect();
    let ray_anchors = endpoint_directions(&net_values, curve.center)?;
    let start = ray_anchors[&0];
    let end = ray_anchors[&0];
    let _ = (start, end);
    let mut anchor_map = BTreeMap::new();
    for (d, p) in &ray_anchors {
        add_ratio_point(
            &mut anchor_map,
            *p,
            anchor_outer(*d, curve.bounds, curve.center)?,
            Some(*d),
        )?;
    }
    let low_roots = crossings(&spans, curve.y_low)?;
    let high_roots = crossings(&spans, curve.y_high)?;
    let mut crossing_points = BTreeMap::new();
    for (roots, y) in [(&low_roots, curve.y_low), (&high_roots, curve.y_high)] {
        for root in roots {
            let side = if compare(root.position[0], curve.center[0]) == std::cmp::Ordering::Less {
                curve.bounds[0]
            } else {
                curve.bounds[1]
            };
            let outer = [int(side), y];
            add_ratio_point(&mut anchor_map, root.position, outer, None)?;
            crossing_points.insert((root.span, root.local), root.position);
        }
    }
    let mut anchor_list: Vec<_> = anchor_map.values().copied().collect();
    // The anchor order is the geometrically prescribed CCW outer cage. This
    // also makes source cyclic shifts and reversals produce the same plan.
    let mut ordered_anchors = Vec::with_capacity(ANCHOR_COUNT);
    let expected_anchors = [
        (Some(0), None, None),
        (Some(1), None, None),
        (Some(2), None, None),
        (None, Some(curve.bounds[0]), Some(curve.y_high)),
        (Some(3), None, None),
        (None, Some(curve.bounds[0]), Some(curve.y_low)),
        (Some(4), None, None),
        (Some(5), None, None),
        (Some(6), None, None),
        (None, Some(curve.bounds[1]), Some(curve.y_low)),
        (Some(7), None, None),
        (None, Some(curve.bounds[1]), Some(curve.y_high)),
    ];
    for (direction, side, y) in expected_anchors {
        if let Some(d) = direction {
            ordered_anchors.push(
                anchor_map
                    .values()
                    .find(|a| a.direction == Some(d))
                    .copied()
                    .ok_or("missing periodic cage ray")?,
            );
        } else {
            let target = [int(side.unwrap()), y.unwrap()];
            ordered_anchors.push(
                anchor_map
                    .values()
                    .find(|a| a.direction.is_none() && a.outer == target)
                    .copied()
                    .ok_or("missing transverse intersection anchor")?,
            );
        }
    }
    anchor_list.clear();
    let anchor_index: BTreeMap<_, _> = ordered_anchors
        .iter()
        .enumerate()
        .map(|(i, a)| (a.point, i))
        .collect();

    // Split source spans exactly at the four plan events and orient each piece
    // in canonical CCW order. The active parameter translation is immaterial:
    // geometric endpoint anchors, rather than knot indices, select the start.
    let mut source_cuts = BTreeMap::<usize, Vec<R>>::new();
    for root in low_roots.iter().chain(&high_roots) {
        source_cuts.entry(root.span).or_default().push(root.local);
    }
    for v in source_cuts.values_mut() {
        v.sort_by(|a, b| compare(*a, *b));
        v.dedup();
    }
    let mut split_by_span = Vec::<Vec<(HNet, R, R)>>::new();
    for (net, source_index, reversed) in &nets {
        let mut cuts = source_cuts.get(source_index).cloned().unwrap_or_default();
        if *reversed {
            cuts = cuts
                .into_iter()
                .map(|t| R::one().checked_sub(t))
                .collect::<Result<_, _>>()?;
        }
        cuts.sort_by(|a, b| compare(*a, *b));
        split_by_span.push(split_piece_at(net.clone(), &cuts)?);
    }
    // Recover the geometric chain from endpoint connectivity. This keeps
    // sector numbering invariant under cyclic reindexing and parameter reversal.
    let mut pieces = Vec::<ProfilePiece>::new();
    for (span_index, subnets) in split_by_span.iter().enumerate() {
        let source_index = nets[span_index].1;
        let knots = &spans[source_index].interval;
        for (net, lo, hi) in subnets {
            pieces.push(ProfilePiece {
                controls: net.clone(),
                start: point(net[0])?,
                end: point(*net.last().unwrap())?,
                duration: int(knots[1]
                    .checked_sub(knots[0])
                    .ok_or("source knot interval exceeds i64")?)
                .checked_mul(hi.checked_sub(*lo)?)?,
                source_span: source_index,
            });
        }
    }
    // `canonical_spans` already orders complete source spans counterclockwise
    // from the northeast anchor; de Casteljau children retain that order.
    let ordered_pieces = pieces;
    if ordered_pieces.first().map(|p| p.start) != Some(ordered_anchors[0].point) {
        return Err("profile does not start at the canonical northeast anchor");
    }

    // Assign every oriented piece to one of the twelve ordered anchor arcs.
    let mut sectors: [Vec<ProfilePiece>; ANCHOR_COUNT] = std::array::from_fn(|_| Vec::new());
    let mut current_sector = 0;
    for piece in ordered_pieces {
        if let Some(&index) = anchor_index.get(&piece.start) {
            current_sector = index;
        }
        if current_sector >= ANCHOR_COUNT {
            return Err("invalid cross-plate sector assignment");
        }
        sectors[current_sector].push(piece.clone());
        if anchor_index.contains_key(&piece.end) {
            current_sector = anchor_index[&piece.end];
        }
    }
    if sectors.iter().any(Vec::is_empty) {
        return Err("periodic profile omitted a required cage sector");
    }

    let mut blocks = Vec::new();
    let mut tags = BTreeMap::<String, Tag>::new();
    let z_events = [0, curve.z_cut[0], curve.z_cut[1], curve.height];
    for layer in 0..3 {
        let z0 = z_events[layer];
        let z1 = z_events[layer + 1];
        for sector in 0..ANCHOR_COUNT {
            if layer == 1 && LATERAL.contains(&sector) {
                continue;
            }
            let left_anchor = ordered_anchors[sector];
            let right_anchor = ordered_anchors[(sector + 1) % ANCHOR_COUNT];
            let is_lateral = LATERAL.contains(&sector);
            let line_side = if is_lateral {
                None
            } else {
                Some(boundary_side(
                    [left_anchor.outer, right_anchor.outer],
                    curve.bounds,
                )?)
            };
            let duration: R = sectors[sector]
                .iter()
                .try_fold(R::zero(), |sum, p| sum.checked_add(p.duration))?;
            let mut offset = R::zero();
            for (piece_index, piece) in sectors[sector].iter().enumerate() {
                let next_offset = offset.checked_add(piece.duration)?;
                let q0 = div(offset, duration)?;
                let q1 = div(next_offset, duration)?;
                let outer_points = if let Some(side) = line_side {
                    let interpolate = |q: R| -> Result<[R; 2], Error> {
                        Ok([
                            left_anchor.outer[0].checked_add(q.checked_mul(
                                right_anchor.outer[0].checked_sub(left_anchor.outer[0])?,
                            )?)?,
                            left_anchor.outer[1].checked_add(q.checked_mul(
                                right_anchor.outer[1].checked_sub(left_anchor.outer[1])?,
                            )?)?,
                        ])
                    };
                    let pts = [interpolate(q0)?, interpolate(q1)?];
                    let _ = side;
                    pts
                } else {
                    [
                        [
                            int(if sector < 6 {
                                curve.bounds[0]
                            } else {
                                curve.bounds[1]
                            }),
                            piece.start[1],
                        ],
                        [
                            int(if sector < 6 {
                                curve.bounds[0]
                            } else {
                                curve.bounds[1]
                            }),
                            piece.end[1],
                        ],
                    ]
                };
                let outer = outer_net(
                    &piece.controls,
                    outer_points,
                    if is_lateral {
                        Some(if sector < 6 {
                            curve.bounds[0]
                        } else {
                            curve.bounds[1]
                        })
                    } else {
                        None
                    },
                )?;
                let controls = block_controls(&outer, &piece.controls, z0, z1)?;
                let id = format!(
                    "g2_s{sector}_p{piece_index}_z{layer}_src{}",
                    piece.source_span
                );
                let block = BracketBlock {
                    id: id.clone(),
                    degrees: [3, 1, 1],
                    controls,
                };
                let certificate =
                    certify_spatial(block.degrees, &block.controls, Default::default());
                if certificate.status != CertificateStatus::Valid {
                    return Err(
                        "cross-plate sector lacks a positive Bernstein determinant certificate",
                    );
                }
                tags.insert(
                    id,
                    Tag {
                        sector,
                        layer,
                        side_min: outer_points[0],
                        side_max: outer_points[1],
                    },
                );
                blocks.push(block);
                offset = next_offset;
            }
        }
    }

    let incidence = enumerate_faces(&blocks)?;
    let blocks_by_id: BTreeMap<_, _> = blocks
        .iter()
        .map(|block| (block.id.as_str(), block))
        .collect();
    let mut faces = Vec::new();
    for side in incidence.exterior {
        let tag = tags
            .get(&side.block)
            .ok_or("missing cross-plate block provenance")?;
        let z0 = z_events[tag.layer];
        let z1 = z_events[tag.layer + 1];
        let feature = match side.axis {
            0 => {
                let neighbor = if side.upper {
                    (tag.sector + 1) % ANCHOR_COUNT
                } else {
                    (tag.sector + ANCHOR_COUNT - 1) % ANCHOR_COUNT
                };
                if tag.layer == 1 && LATERAL.contains(&neighbor) {
                    let block = blocks_by_id
                        .get(side.block.as_str())
                        .ok_or("missing exposed cross-plate block")?;
                    let face_y = transverse_face_y(block, &side)?;
                    if face_y == curve.y_low {
                        "rect_x/ymin".to_owned()
                    } else if face_y == curve.y_high {
                        "rect_x/ymax".to_owned()
                    } else {
                        return Err(
                            "exposed transverse side is not on a registered rectangle edge",
                        );
                    }
                } else {
                    return Err("unexpected exposed periodic sector seam");
                }
            }
            1 if side.upper => "nurbs_z_hole/wall".to_owned(),
            1 => {
                let p = tag.side_min;
                let q = tag.side_max;
                if p[0] == int(curve.bounds[0]) && q[0] == int(curve.bounds[0]) {
                    "stock/xmin".to_owned()
                } else if p[0] == int(curve.bounds[1]) && q[0] == int(curve.bounds[1]) {
                    "stock/xmax".to_owned()
                } else if p[1] == int(curve.bounds[2]) && q[1] == int(curve.bounds[2]) {
                    "stock/ymin".to_owned()
                } else if p[1] == int(curve.bounds[3]) && q[1] == int(curve.bounds[3]) {
                    "stock/ymax".to_owned()
                } else {
                    return Err("outer guard face is not on the stock boundary");
                }
            }
            2 if !side.upper && z0 == 0 => "stock/zmin".to_owned(),
            2 if side.upper && z1 == curve.height => "stock/zmax".to_owned(),
            2 if side.upper && z1 == curve.z_cut[0] && LATERAL.contains(&tag.sector) => {
                "rect_x/zmin".to_owned()
            }
            2 if !side.upper && z0 == curve.z_cut[1] && LATERAL.contains(&tag.sector) => {
                "rect_x/zmax".to_owned()
            }
            2 => return Err("unexpected exposed cross-plate layer face"),
            _ => return Err("unexpected cross-plate boundary side"),
        };
        faces.push((side, feature));
        let _ = z0;
        let _ = z1;
    }
    Ok((blocks, faces))
}

#[path = "crossed_tool_boolean.rs"]
pub mod boolean;
