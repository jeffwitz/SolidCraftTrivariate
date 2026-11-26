//! Exact recursive separation of disjoint convex polygonal cuts.
//!
//! A result is a list of alternative partitions. Each partition contains one
//! counter-clockwise host polygon per input hole, in the same order as `holes`.
//! Chords are inserted as shared edges, and all chord endpoints are propagated
//! to every host edge on which they lie so neighboring hosts remain conforming.

use crate::constraint_system::Ratio;
use std::collections::{BTreeMap, BTreeSet};

type Point = [Ratio; 2];

fn as_f64(x: Ratio) -> f64 {
    x.numerator() as f64 / x.denominator() as f64
}
fn r0() -> Ratio {
    Ratio::zero()
}
fn r1() -> Ratio {
    Ratio::new(1, 1).expect("one")
}
fn sub(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_sub(b)
}
fn add(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_add(b)
}
fn mul(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_mul(b)
}
fn div(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    mul(a, Ratio::new(b.denominator(), b.numerator())?)
}
fn cross(a: Point, b: Point, c: Point) -> Result<Ratio, &'static str> {
    sub(
        mul(sub(b[0], a[0])?, sub(c[1], a[1])?)?,
        mul(sub(b[1], a[1])?, sub(c[0], a[0])?)?,
    )
}
fn sign(x: Ratio) -> i8 {
    x.numerator().signum() as i8
}
fn between(a: Ratio, x: Ratio, b: Ratio) -> bool {
    (a.cmp_value(&x).is_le() && x.cmp_value(&b).is_le())
        || (b.cmp_value(&x).is_le() && x.cmp_value(&a).is_le())
}
fn on_segment(a: Point, b: Point, p: Point) -> Result<bool, &'static str> {
    Ok(sign(cross(a, b, p)?) == 0 && between(a[0], p[0], b[0]) && between(a[1], p[1], b[1]))
}
fn unique_push(points: &mut Vec<Point>, p: Point) {
    if !points.contains(&p) {
        points.push(p);
    }
}
fn cmp_point(a: Point, b: Point) -> std::cmp::Ordering {
    a[0].cmp_value(&b[0]).then_with(|| a[1].cmp_value(&b[1]))
}
fn line_key(a: Point, b: Point) -> Result<(u8, Ratio, Ratio), &'static str> {
    let dx = sub(b[0], a[0])?;
    let dy = sub(b[1], a[1])?;
    if sign(dx) == 0 {
        Ok((0, a[0], r0()))
    } else {
        let slope = div(dy, dx)?;
        Ok((1, slope, sub(a[1], mul(slope, a[0])?)?))
    }
}

/// Intersect the infinite line `a + t*(b-a)` with the boundary of a convex host.
fn line_chord(host: &[Point], a: Point, b: Point) -> Result<Option<(Point, Point)>, &'static str> {
    let d = [sub(b[0], a[0])?, sub(b[1], a[1])?];
    let mut hits = Vec::new();
    for i in 0..host.len() {
        let p = host[i];
        let q = host[(i + 1) % host.len()];
        let e = [sub(q[0], p[0])?, sub(q[1], p[1])?];
        let den = sub(mul(d[0], e[1])?, mul(d[1], e[0])?)?;
        if sign(den) == 0 {
            continue;
        }
        let pa = [sub(p[0], a[0])?, sub(p[1], a[1])?];
        let t = div(sub(mul(pa[0], e[1])?, mul(pa[1], e[0])?)?, den)?;
        let s = div(sub(mul(pa[0], d[1])?, mul(pa[1], d[0])?)?, den)?;
        if s.cmp_value(&r0()).is_lt() || s.cmp_value(&r1()).is_gt() {
            continue;
        }
        let x = [add(a[0], mul(t, d[0])?)?, add(a[1], mul(t, d[1])?)?];
        unique_push(&mut hits, x);
    }
    if hits.len() < 2 {
        return Ok(None);
    }
    // A convex polygon meets a line in at most one interval; order hits by t.
    let parameter = |p: Point| -> Result<Ratio, &'static str> {
        if sign(d[0]) != 0 {
            div(sub(p[0], a[0])?, d[0])
        } else {
            div(sub(p[1], a[1])?, d[1])
        }
    };
    let mut keyed = hits
        .into_iter()
        .map(|p| Ok((parameter(p)?, p)))
        .collect::<Result<Vec<_>, &'static str>>()?;
    keyed.sort_by(|x, y| x.0.cmp_value(&y.0));
    Ok(Some((keyed[0].1, keyed[keyed.len() - 1].1)))
}

fn insert_on_boundary(poly: &[Point], p: Point) -> Result<Vec<Point>, &'static str> {
    if poly.contains(&p) {
        return Ok(poly.to_vec());
    }
    let mut out = Vec::with_capacity(poly.len() + 1);
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        out.push(a);
        if on_segment(a, b, p)? {
            out.push(p);
        }
    }
    if out.len() == poly.len() {
        return Err("separator endpoint is not on host boundary");
    }
    Ok(out)
}

#[derive(Clone)]
struct Work {
    host: Vec<Point>,
    hole_ids: Vec<usize>,
}
type SplitCandidate = (Work, Work, (Point, Point));
type HostPlan = Vec<(usize, Vec<Point>)>;
type HostPlanAlternatives = Vec<HostPlan>;

fn candidate_lines(
    host: &[Point],
    holes: &[Vec<Point>],
    ids: &[usize],
    lookahead: bool,
) -> Result<Vec<(Point, Point)>, &'static str> {
    let mut lines = Vec::new();
    // Existing host vertices provide robust, exact separator directions.
    for i in 0..host.len() {
        for j in (i + 1)..host.len() {
            if j == i + 1 || (i == 0 && j + 1 == host.len()) {
                continue;
            }
            lines.push((host[i], host[j]));
        }
    }
    // Exact support separators for every nontrivial partition of the holes.
    // The normal is perpendicular to an edge of either group's convex hull;
    // a strict projection gap then yields a rational mid-gap separator.
    let n = ids.len();
    let mut masks = BTreeSet::new();
    if n <= 6 {
        for mask in 1usize..(1usize << n) - 1 {
            if mask & 1 == 0 {
                continue;
            } // one representative per complementary pair
            masks.insert(mask);
        }
    } else {
        // Singleton and pair groups cover local separations without enumerating
        // all 2^n subsets. Add spatial sweeps along independent axes so
        // diagonally separated clusters with overlapping x/y bounds are found.
        for i in 0..n {
            masks.insert(1usize << i);
            for j in (i + 1)..n {
                masks.insert((1usize << i) | (1usize << j));
            }
        }
        for (nx, ny) in [(1i128, 0i128), (0, 1), (1, 1), (1, -1)] {
            let mut ordered = ids
                .iter()
                .enumerate()
                .map(|(i, &id)| {
                    let mut sx = r0();
                    let mut sy = r0();
                    for p in &holes[id] {
                        sx = add(sx, p[0])?;
                        sy = add(sy, p[1])?;
                    }
                    let count = Ratio::new(holes[id].len() as i128, 1)?;
                    let cx = div(sx, count)?;
                    let cy = div(sy, count)?;
                    let projection =
                        add(mul(Ratio::new(nx, 1)?, cx)?, mul(Ratio::new(ny, 1)?, cy)?)?;
                    Ok::<_, &'static str>((projection, i))
                })
                .collect::<Result<Vec<_>, _>>()?;
            ordered.sort_by(|a, b| a.0.cmp_value(&b.0));
            for k in 1..n {
                let mask = ordered[..k].iter().fold(0usize, |m, x| m | (1usize << x.1));
                let complement = ((1usize << n) - 1) ^ mask;
                masks.insert(mask.min(complement));
            }
        }
    }
    for mask in masks {
        let a_pts = ids
            .iter()
            .enumerate()
            .filter(|(i, _)| mask & (1usize << i) != 0)
            .flat_map(|(_, &id)| holes[id].iter().copied())
            .collect::<Vec<_>>();
        let b_pts = ids
            .iter()
            .enumerate()
            .filter(|(i, _)| mask & (1usize << i) == 0)
            .flat_map(|(_, &id)| holes[id].iter().copied())
            .collect::<Vec<_>>();
        let ah = convex_hull(&a_pts);
        let bh = convex_hull(&b_pts);
        for hull in [&ah, &bh] {
            let stride = hull.len().div_ceil(32).max(1);
            for i in 0..hull.len() {
                if i % stride != 0 {
                    continue;
                }
                let p = hull[i];
                let q = hull[(i + 1) % hull.len()];
                let edge = [sub(q[0], p[0])?, sub(q[1], p[1])?];
                let normal = [edge[1], sub(r0(), edge[0])?];
                if sign(normal[0]) == 0 && sign(normal[1]) == 0 {
                    continue;
                }
                let project = |x: Point| -> Result<Ratio, &'static str> {
                    add(mul(normal[0], x[0])?, mul(normal[1], x[1])?)
                };
                let Ok(ap) = ah
                    .iter()
                    .copied()
                    .map(project)
                    .collect::<Result<Vec<_>, _>>()
                else {
                    continue;
                };
                let Ok(bp) = bh
                    .iter()
                    .copied()
                    .map(project)
                    .collect::<Result<Vec<_>, _>>()
                else {
                    continue;
                };
                let (amax, bmin) = (
                    *ap.iter()
                        .max_by(|x, y| x.cmp_value(y))
                        .expect("nonempty group"),
                    *bp.iter()
                        .min_by(|x, y| x.cmp_value(y))
                        .expect("nonempty group"),
                );
                let (bmax, amin) = (
                    *bp.iter()
                        .max_by(|x, y| x.cmp_value(y))
                        .expect("nonempty group"),
                    *ap.iter()
                        .min_by(|x, y| x.cmp_value(y))
                        .expect("nonempty group"),
                );
                let gap = if amax.cmp_value(&bmin).is_lt() {
                    Some((amax, bmin))
                } else if bmax.cmp_value(&amin).is_lt() {
                    Some((bmax, amin))
                } else {
                    None
                };
                let Some((lo, hi)) = gap else {
                    continue;
                };
                // Mid-gap and one-third offsets retain exact positive clearance
                // from both convex group hulls and create distinct viable chords.
                let mut weights = vec![Ratio::new(1, 3)?, Ratio::new(1, 2)?, Ratio::new(2, 3)?];
                if lookahead {
                    weights.extend([Ratio::new(1, 4)?, Ratio::new(3, 4)?]);
                }
                for weight in weights {
                    let offset = add(lo, mul(sub(hi, lo)?, weight)?)?;
                    let anchor = if sign(normal[0]) != 0 {
                        [div(offset, normal[0])?, Ratio::zero()]
                    } else {
                        [Ratio::zero(), div(offset, normal[1])?]
                    };
                    let direction = [sub(r0(), normal[1])?, normal[0]];
                    let endpoint = [add(anchor[0], direction[0])?, add(anchor[1], direction[1])?];
                    lines.push((anchor, endpoint));
                }
            }
        }
    }
    let mut unique = BTreeMap::new();
    for (a, b) in lines {
        if a == b {
            continue;
        }
        unique.entry(line_key(a, b)?).or_insert((a, b));
    }
    Ok(unique.into_values().collect())
}

fn convex_hull(points: &[Point]) -> Vec<Point> {
    let mut p = points.to_vec();
    p.sort_by(|a, b| cmp_point(*a, *b));
    p.dedup();
    if p.len() <= 2 {
        return p;
    }
    let half = |seq: &[Point]| {
        let mut h = Vec::new();
        for &x in seq {
            while h.len() >= 2
                && cross(h[h.len() - 2], h[h.len() - 1], x).is_ok_and(|v| sign(v) <= 0)
            {
                h.pop();
            }
            h.push(x);
        }
        h
    };
    let lo = half(&p);
    let rev = p.iter().rev().copied().collect::<Vec<_>>();
    let hi = half(&rev);
    lo[..lo.len() - 1]
        .iter()
        .chain(hi[..hi.len() - 1].iter())
        .copied()
        .collect()
}

// Chord endpoints are distinct bounded rational coordinates. Carry the
// rounding residuals through the square root so clearance ties do not
// depend on the system libm implementation of hypot.
fn chord_length(a: Point, b: Point) -> f64 {
    let dx = (as_f64(b[0]) - as_f64(a[0])).abs();
    let dy = (as_f64(b[1]) - as_f64(a[1])).abs();
    let (long, short) = (dx.max(dy), dx.min(dy));
    let long_sq = long * long;
    let short_sq = short * short;
    let sum = long_sq + short_sq;
    let root = sum.sqrt();
    if root == 0.0 {
        return root;
    }
    let sum_error = (long_sq - sum) + short_sq;
    let error = sum_error + long.mul_add(long, -long_sq) + short.mul_add(short, -short_sq);
    let residual = (-root).mul_add(root, sum) + error;
    root + residual / (2.0 * root)
}

fn split(
    work: &Work,
    holes: &[Vec<Point>],
    rank_clearance: bool,
    lookahead: bool,
) -> Result<Vec<SplitCandidate>, &'static str> {
    let mut candidates = Vec::new();
    for (la, lb) in candidate_lines(&work.host, holes, &work.hole_ids, lookahead)? {
        let Some((a, b)) = line_chord(&work.host, la, lb)? else {
            continue;
        };
        if a == b {
            continue;
        }
        let mut left = Vec::new();
        let mut right = Vec::new();
        let length = chord_length(a, b);
        let mut clearance = f64::INFINITY;
        let mut valid = true;
        for &hid in &work.hole_ids {
            let mut side = 0;
            for &p in &holes[hid] {
                let turn = cross(a, b, p)?;
                clearance = clearance.min(as_f64(turn).abs() / length);
                let s = sign(turn);
                if s == 0 || (side != 0 && side != s) {
                    valid = false;
                    break;
                }
                side = s;
            }
            if !valid {
                break;
            }
            if side > 0 {
                left.push(hid);
            } else {
                right.push(hid);
            }
        }
        if !valid || left.is_empty() || right.is_empty() {
            continue;
        }
        // Both searches balance the hole counts. The second favors wide gaps
        // and separatrices whose endpoints already belong to the host.
        let score = left.len().abs_diff(right.len());
        let new_vertices =
            usize::from(!work.host.contains(&a)) + usize::from(!work.host.contains(&b));
        let margin = if rank_clearance {
            clearance * (1.0 - 0.15 * new_vertices as f64)
        } else {
            0.0
        };
        candidates.push((score, a, b, left, right, margin));
    }
    candidates.sort_by(|x, y| {
        x.0.cmp(&y.0)
            .then_with(|| y.5.total_cmp(&x.5))
            .then_with(|| x.1.cmp(&y.1))
            .then_with(|| x.2.cmp(&y.2))
    });
    // Greedily retain geometrically diverse splits: same group with a parallel
    // chord is a near-duplicate unless no nonparallel candidate exists.
    let mut result = Vec::new();
    let mut rejected_parallel = Vec::new();
    for (_, a, b, left_ids, right_ids, _) in candidates {
        let mut boundary = insert_on_boundary(&work.host, a)?;
        boundary = insert_on_boundary(&boundary, b)?;
        let Some(i) = boundary.iter().position(|p| *p == a) else {
            continue;
        };
        let Some(j) = boundary.iter().position(|p| *p == b) else {
            continue;
        };
        if i == j {
            continue;
        }
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        let first = boundary[lo..=hi].to_vec();
        let second = boundary[hi..]
            .iter()
            .chain(boundary[..=lo].iter())
            .copied()
            .collect::<Vec<_>>();
        // Both slices traverse the shared chord oppositely. Choose hole side by exact orientation.
        let ref_point = first[1];
        let ref_side = sign(cross(a, b, ref_point)?);
        let (left_host, right_host) = if ref_side > 0 {
            (first, second)
        } else {
            (second, first)
        };
        let candidate = (
            Work {
                host: left_host,
                hole_ids: left_ids,
            },
            Work {
                host: right_host,
                hole_ids: right_ids,
            },
            (a, b),
        );
        let same_partition = |x: &[usize], y: &[usize]| {
            x == y
                || x.iter().copied().collect::<BTreeSet<_>>()
                    == y.iter().copied().collect::<BTreeSet<_>>()
        };
        let parallel_to_selected =
            result
                .iter()
                .any(|(x, _y, ch): &(Work, Work, (Point, Point))| {
                    let (ux, uy) = (sub(ch.1[0], ch.0[0]), sub(ch.1[1], ch.0[1]));
                    let (vx, vy) = (
                        sub(candidate.2 .1[0], candidate.2 .0[0]),
                        sub(candidate.2 .1[1], candidate.2 .0[1]),
                    );
                    same_partition(&x.hole_ids, &candidate.0.hole_ids)
                        && ux
                            .ok()
                            .zip(uy.ok())
                            .zip(vx.ok())
                            .zip(vy.ok())
                            .and_then(|(((ux, uy), vx), vy)| {
                                cross([r0(), r0()], [ux, uy], [vx, vy]).ok()
                            })
                            .is_some_and(|c| sign(c) == 0)
                });
        if result
            .iter()
            .all(|(x, y, _)| x.host != candidate.0.host || y.host != candidate.1.host)
        {
            if parallel_to_selected {
                rejected_parallel.push(candidate);
            } else {
                result.push(candidate);
            }
        }
        if result.len() >= if lookahead { 16 } else { 4 } {
            break;
        }
    }
    for candidate in rejected_parallel {
        if result.len() >= if lookahead { 16 } else { 4 } {
            break;
        }
        result.push(candidate);
    }
    Ok(result)
}

fn recurse(
    work: Work,
    holes: &[Vec<Point>],
    rank_clearance: bool,
    lookahead: bool,
    budget: &mut usize,
) -> Result<HostPlanAlternatives, &'static str> {
    if work.hole_ids.len() == 1 {
        return Ok(vec![vec![(work.hole_ids[0], work.host)]]);
    }
    if *budget == 0 {
        return Ok(Vec::new());
    }
    *budget -= 1;
    let splits = split(&work, holes, rank_clearance, lookahead)?;
    if splits.is_empty() {
        return Err("no separating chord for convex hosts");
    }
    let mut plans = Vec::new();
    for (a, b, _chord) in splits {
        let Ok(left) = recurse(a, holes, rank_clearance, lookahead, budget) else {
            continue;
        };
        if left.is_empty() {
            continue;
        }
        let Ok(right) = recurse(b, holes, rank_clearance, lookahead, budget) else {
            continue;
        };
        for l in &left {
            for r in &right {
                let mut plan = l.clone();
                plan.extend(r.clone());
                if plans.iter().all(|p: &HostPlan| p != &plan) {
                    plans.push(plan);
                }
                if plans.len() >= if lookahead { 8 } else { 4 } {
                    break;
                }
            }
            if plans.len() >= if lookahead { 8 } else { 4 } {
                break;
            }
        }
        if plans.len() >= 4 {
            break;
        }
    }
    Ok(plans)
}

fn close_host(poly: &[Point], special: &BTreeSet<Point>) -> Result<Vec<Point>, &'static str> {
    let mut out = Vec::new();
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        out.push(a);
        let mut points = Vec::new();
        for &p in special {
            if p != a && p != b && on_segment(a, b, p)? {
                points.push(p);
            }
        }
        let axis = if a[0] != b[0] { 0 } else { 1 };
        let direction = a[axis].cmp_value(&b[axis]);
        points.sort_by(|x, y| {
            let order = x[axis].cmp_value(&y[axis]);
            if direction.is_gt() {
                order.reverse()
            } else {
                order
            }
        });
        out.extend(points);
    }
    Ok(out)
}

/// Return bounded alternative recursive convex-host partitions.
/// Each alternative is indexed by input hole and contains its CCW host polygon.
pub(super) fn host_candidates(stock: &[Point], holes: &[Vec<Point>]) -> Vec<Vec<Vec<Point>>> {
    host_candidates_mode(stock, holes, false)
}

/// Wider deterministic search used only by the opt-in planar ring lookahead.
/// The added exact support offsets produce previously unseen chord endpoints.
pub(super) fn host_candidates_lookahead(
    stock: &[Point],
    holes: &[Vec<Point>],
) -> Vec<Vec<Vec<Point>>> {
    host_candidates_mode(stock, holes, true)
}

fn host_candidates_mode(
    stock: &[Point],
    holes: &[Vec<Point>],
    lookahead: bool,
) -> Vec<Vec<Vec<Point>>> {
    if stock.len() < 3 || holes.is_empty() || holes.iter().any(|h| h.len() < 3) {
        return Vec::new();
    }
    let start = Work {
        host: stock.to_vec(),
        hole_ids: (0..holes.len()).collect(),
    };
    let plan_limit = if lookahead { 8 } else { 4 };
    let mut budget = if lookahead { 128 } else { 64 };
    let mut plans =
        recurse(start.clone(), holes, false, lookahead, &mut budget).unwrap_or_default();
    if (2..=6).contains(&holes.len()) {
        let mut budget = if lookahead { 128 } else { 64 };
        for plan in recurse(start, holes, true, lookahead, &mut budget).unwrap_or_default() {
            if !plans.contains(&plan) {
                plans.push(plan)
            }
            if lookahead && plans.len() >= plan_limit {
                break;
            }
        }
    }
    plans
        .into_iter()
        .filter_map(|plan| {
            let special = plan
                .iter()
                .flat_map(|(_, p)| p.iter().copied())
                .collect::<BTreeSet<_>>();
            let mut result = vec![Vec::new(); holes.len()];
            for (id, poly) in plan {
                result[id] = close_host(&poly, &special).ok()?;
            }
            if result.iter().any(Vec::is_empty) {
                None
            } else {
                Some(result)
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: i128, y: i128) -> Point {
        [Ratio::new(x, 1).unwrap(), Ratio::new(y, 1).unwrap()]
    }
    fn diamond(cx: i128, cy: i128, rx: i128, ry: i128) -> Vec<Point> {
        vec![
            p(cx, cy - ry),
            p(cx + rx, cy),
            p(cx, cy + ry),
            p(cx - rx, cy),
        ]
    }
    fn in_ccw(poly: &[Point], point: Point) -> bool {
        (0..poly.len())
            .all(|i| cross(poly[i], poly[(i + 1) % poly.len()], point).is_ok_and(|v| sign(v) >= 0))
    }
    fn check_multiple_plans(holes: Vec<Vec<Point>>) {
        let stock = vec![p(-12, -10), p(12, -10), p(12, 10), p(-12, 10)];
        let plans = host_candidates(&stock, &holes);
        assert!(
            plans.len() >= 4,
            "expected four distinct plans, found {}",
            plans.len()
        );
        for plan in &plans {
            assert_eq!(plan.len(), holes.len());
            for (host, hole) in plan.iter().zip(&holes) {
                assert!(
                    hole.iter().all(|&vertex| in_ccw(host, vertex)),
                    "hole is outside its convex host"
                );
            }
        }
        for i in 0..plans.len() {
            for j in (i + 1)..plans.len() {
                assert_ne!(plans[i], plans[j], "plans must be geometrically distinct");
            }
        }
    }

    #[test]
    fn exact_support_lines_separate_diamonds_with_overlapping_projection() {
        // Their y projections coincide while their x projections are disjoint.
        check_multiple_plans(vec![diamond(-3, 0, 2, 3), diamond(3, 0, 2, 3)]);
    }

    #[test]
    fn exact_support_lines_produce_multiple_plans_for_three_holes() {
        check_multiple_plans(vec![
            diamond(-6, 0, 1, 2),
            diamond(0, 3, 2, 1),
            diamond(6, 0, 1, 2),
        ]);
    }

    #[test]
    fn exact_support_lines_separate_diamonds_overlapping_both_axis_bounds() {
        // Long diagonal parallelograms are separated along x-y while both their
        // x and y bounding intervals overlap.
        let ten = |n: i128| Ratio::new(n, 10).unwrap();
        let tilted = |cx: i128, cy: i128| {
            let c = [cx, cy];
            let offsets = [[-25, 15], [15, 25], [25, -15], [-15, -25]];
            offsets
                .iter()
                .map(|o| [ten(c[0] + o[0]), ten(c[1] + o[1])])
                .collect::<Vec<_>>()
        };
        let stock = vec![p(-12, -10), p(12, -10), p(12, 10), p(-12, 10)];
        let holes = vec![tilted(-23, -23), tilted(23, 23)];
        assert!(holes[0]
            .iter()
            .map(|p| p[0])
            .max_by(|a, b| a.cmp_value(b))
            .unwrap()
            .cmp_value(
                &holes[1]
                    .iter()
                    .map(|p| p[0])
                    .min_by(|a, b| a.cmp_value(b))
                    .unwrap()
            )
            .is_gt());
        let plans = host_candidates(&stock, &holes);
        assert!(
            !plans.is_empty(),
            "diagonal separating support must be found"
        );
        for plan in plans {
            for (host, hole) in plan.iter().zip(&holes) {
                assert!(hole.iter().all(|&vertex| in_ccw(host, vertex)));
            }
        }
    }

    #[test]
    fn seeded_polygon_cases_10_and_19_find_convex_hosts() {
        let as_points = |vertices: &[[i32; 2]]| {
            vertices
                .iter()
                .map(|xy| p(i128::from(xy[0]), i128::from(xy[1])))
                .collect::<Vec<_>>()
        };
        let cases: [(Vec<Point>, Vec<Vec<Point>>); 2] = [
            (
                as_points(&[
                    [-893, 52],
                    [-1263, -448],
                    [-1263, -818],
                    [-893, -948],
                    [-263, -818],
                    [477, -448],
                    [1107, 52],
                    [1477, 552],
                    [1477, 922],
                    [1107, 1052],
                    [477, 922],
                    [-263, 552],
                ]),
                vec![
                    as_points(&[
                        [187, 52],
                        [113, -48],
                        [113, -122],
                        [187, -148],
                        [313, -122],
                        [461, -48],
                        [587, 52],
                        [661, 152],
                        [661, 226],
                        [587, 252],
                        [461, 226],
                        [313, 152],
                    ]),
                    as_points(&[
                        [-373, 52],
                        [-447, -48],
                        [-447, -122],
                        [-373, -148],
                        [-247, -122],
                        [-99, -48],
                        [27, 52],
                        [101, 152],
                        [101, 226],
                        [27, 252],
                        [-99, 226],
                        [-247, 152],
                    ]),
                ],
            ),
            (
                as_points(&[
                    [167, -754],
                    [757, -1744],
                    [1117, -1964],
                    [1117, -1344],
                    [757, -124],
                    [167, 1246],
                    [-423, 2236],
                    [-783, 2456],
                    [-783, 1836],
                    [-423, 616],
                ]),
                vec![
                    as_points(&[
                        [167, 326],
                        [285, 128],
                        [357, 84],
                        [357, 208],
                        [285, 452],
                        [167, 726],
                        [49, 924],
                        [-23, 968],
                        [-23, 844],
                        [49, 600],
                    ]),
                    as_points(&[
                        [167, -234],
                        [285, -432],
                        [357, -476],
                        [357, -352],
                        [285, -108],
                        [167, 166],
                        [49, 364],
                        [-23, 408],
                        [-23, 284],
                        [49, 40],
                    ]),
                ],
            ),
        ];
        for (stock, holes) in cases {
            let plans = host_candidates(&stock, &holes);
            assert!(!plans.is_empty());
            for plan in plans {
                for (host, hole) in plan.iter().zip(&holes) {
                    assert!(hole.iter().all(|&vertex| in_ccw(host, vertex)));
                }
            }
        }
    }

    #[test]
    fn bounded_search_handles_twelve_sixty_four_vertex_loops() {
        let subdivided_rectangle = |cx: i32, cy: i32, half: i32| {
            let corners = [
                [cx - half, cy - half],
                [cx + half, cy - half],
                [cx + half, cy + half],
                [cx - half, cy + half],
            ];
            let mut loop_points = Vec::new();
            for edge in 0..4 {
                let a = corners[edge];
                let b = corners[(edge + 1) % 4];
                for k in 0..16 {
                    let x = a[0] + (b[0] - a[0]) * k / 16;
                    let y = a[1] + (b[1] - a[1]) * k / 16;
                    loop_points.push(p(i128::from(x), i128::from(y)));
                }
            }
            loop_points
        };
        let stock = subdivided_rectangle(0, 0, 100);
        let centers = [
            (-70, -60),
            (-25, -60),
            (25, -60),
            (70, -60),
            (-70, 0),
            (-25, 0),
            (25, 0),
            (70, 0),
            (-70, 60),
            (-25, 60),
            (25, 60),
            (70, 60),
        ];
        let holes = centers
            .iter()
            .map(|&(x, y)| subdivided_rectangle(x, y, 5))
            .collect::<Vec<_>>();
        let plans = host_candidates(&stock, &holes);
        assert!(
            !plans.is_empty(),
            "bounded candidate search should separate the grid"
        );
        assert!(plans.len() <= 4);
    }
}
