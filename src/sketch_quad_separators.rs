//! Propose local diagonals in a simple polygon, then join each quadrangular
//! cut to the corners of its enclosing cell. Every proposal still goes
//! through the arrangement, coverage, trace and Jacobian checks.
use crate::sketch_arrangement::{Role, Segment};
use crate::sketch_layers::{edges, orient as turn, signed_area2, SketchLayers, SketchLoop};
use std::collections::BTreeSet;

type Point = [i32; 2];
type Quad = [Point; 4];
type HostSplit = (Quad, Vec<usize>, Quad, Vec<usize>, [Point; 2]);
type ScoredHostSplit = (
    (usize, Point, Point),
    Quad,
    Vec<usize>,
    Quad,
    Vec<usize>,
    [Point; 2],
);

fn oriented_loop(contour: &SketchLoop) -> Vec<Point> {
    let mut points = contour.vertices.clone();
    if signed_area2(contour) < 0 {
        points.reverse();
    }
    let start = points.iter().enumerate().min_by_key(|(_, p)| *p).unwrap().0;
    points.rotate_left(start);
    points
}

fn convex(q: &Quad) -> bool {
    (0..4).all(|i| turn(q[i], q[(i + 1) % 4], q[(i + 2) % 4]) > 0)
}

fn inside_quad(q: &Quad, point: Point, strict: bool) -> bool {
    (0..4).all(|i| turn(q[i], q[(i + 1) % 4], point) >= i64::from(strict))
}

fn segment(id: String, a: Point, b: Point) -> Segment {
    Segment {
        id,
        a,
        b,
        role: Role::Construction,
    }
}

/// A deterministic candidate, not an exhaustive search over quadrangulations.
/// An ear contains four consecutive outer vertices. Its diagonal cannot cut
/// a hole or touch another vertex. Failure leaves the baseline available.
pub(crate) fn local_quad_guides(job: &SketchLayers) -> Option<Vec<Segment>> {
    if !job.guides.is_empty()
        || job.outer.vertices.len() < 4
        || !job.outer.vertices.len().is_multiple_of(2)
        || job.cuts.iter().any(|(l, _)| l.vertices.len() != 4)
    {
        return None;
    }
    let mut ring = oriented_loop(&job.outer);
    let mut cells = Vec::new();
    let mut guides = Vec::new();
    while ring.len() > 4 {
        let size = ring.len();
        let (start, ear) = (0..size).find_map(|i| {
            let q = std::array::from_fn(|j| ring[(i + j) % size]);
            if !convex(&q)
                || (4..size).any(|j| inside_quad(&q, ring[(i + j) % size], false))
                || job.cuts.iter().any(|(l, _)| {
                    edges(l).any(|(a, b)| crate::sketch_layers::hits(q[0], q[3], a, b))
                })
            {
                return None;
            }
            Some((i, q))
        })?;
        guides.push(segment(
            format!(
                "local-diagonal-{}-{}-{}-{}",
                ear[0][0], ear[0][1], ear[3][0], ear[3][1]
            ),
            ear[0],
            ear[3],
        ));
        ring.rotate_left(start);
        ring.drain(1..3);
        cells.push(ear);
    }
    let last: Quad = ring.try_into().ok()?;
    if !convex(&last) {
        return None;
    }
    cells.push(last);
    let mut occupied = vec![false; cells.len()];
    for (hole, _) in &job.cuts {
        let inner: Quad = oriented_loop(hole).try_into().ok()?;
        if !convex(&inner) {
            return None;
        }
        let (index, outer) = cells
            .iter()
            .enumerate()
            .find(|(_, q)| inner.iter().all(|&p| inside_quad(q, p, true)))?;
        if occupied[index] {
            return None;
        }
        occupied[index] = true;
        // With fixed winding there are four cyclic pairings. Among admissible
        // pairings choose the one with the least total squared guide length.
        join_rings(outer, &inner, &hole.id, &mut guides)?;
    }
    Some(guides)
}

fn lower_midpoint(a: i32, b: i32) -> i32 {
    (i64::from(a) + i64::from(b)).div_euclid(2) as i32
}

/// Split rectangular stock into disjoint strips, one per convex four-sided
/// cut. For a six-sided cut, add two opposite boundary vertices and a
/// diagonal of the cut so its footprint also conforms below a pocket floor.
pub(crate) fn rectangular_hosts(job: &SketchLayers, axis: usize) -> Option<Vec<Segment>> {
    if !job.guides.is_empty() || job.operations.is_some() || job.force_triangulation {
        return None;
    }
    let outer = oriented_loop(&job.outer);
    if outer.len() != 4 || !edges(&job.outer).all(|(a, b)| a[0] == b[0] || a[1] == b[1]) {
        return None;
    }
    let lo = outer[0];
    let hi = outer[2];
    let mut guides = Vec::new();
    if let [(hole, _)] = job.cuts.as_slice() {
        let inner = oriented_loop(hole);
        if inner.len() != 6
            || !(0..6).all(|i| turn(inner[i], inner[(i + 1) % 6], inner[(i + 2) % 6]) > 0)
        {
            return None;
        }
        let mut perimeter = Vec::new();
        for i in 0..4 {
            let a = outer[i];
            let b = outer[(i + 1) % 4];
            perimeter.push(a);
            if a[axis] != b[axis] {
                let mid = std::array::from_fn(|d| lower_midpoint(a[d], b[d]));
                if mid == a || mid == b {
                    return None;
                }
                perimeter.push(mid);
            }
        }
        join_rings(&perimeter, &inner, &hole.id, &mut guides)?;
        guides.push(segment(format!("cap-{}", hole.id), inner[0], inner[3]));
        return Some(guides);
    }
    if job.cuts.len() < 2 {
        return None;
    }
    let mut holes = Vec::new();
    for (hole, _) in &job.cuts {
        let inner: Quad = oriented_loop(hole).try_into().ok()?;
        if !convex(&inner) {
            return None;
        }
        let min = inner.iter().map(|p| p[axis]).min()?;
        let max = inner.iter().map(|p| p[axis]).max()?;
        holes.push((min, max, &hole.id, inner));
    }
    holes.sort_by_key(|h| (h.0, h.1));
    let mut left = lo[axis];
    for (i, (_, end, id, inner)) in holes.iter().enumerate() {
        let mut right = hi[axis];
        if let Some(next) = holes.get(i + 1) {
            right = lower_midpoint(*end, next.0);
            if right <= *end || right >= next.0 {
                return None;
            }
            let mut a = lo;
            let mut b = hi;
            a[axis] = right;
            b[axis] = right;
            guides.push(segment(format!("strip-{axis}-{right}"), a, b));
        }
        let host: Quad = std::array::from_fn(|j| {
            let mut p = outer[j];
            p[axis] = if p[axis] == lo[axis] { left } else { right };
            p
        });
        join_rings(&host, inner, id, &mut guides)?;
        left = right;
    }
    Some(guides)
}

fn join_rings(outer: &[Point], inner: &[Point], id: &str, guides: &mut Vec<Segment>) -> Option<()> {
    let n = inner.len();
    let shift = (0..n)
        .filter(|s| {
            (0..n).all(|i| {
                let j = (i + 1) % n;
                convex(&[outer[i], outer[j], inner[(j + s) % n], inner[(i + s) % n]])
            })
        })
        .min_by_key(|s| {
            (0..n)
                .map(|i| {
                    let j = (i + s) % n;
                    (0..2)
                        .map(|a| i64::from(outer[i][a] - inner[j][a]).pow(2))
                        .sum::<i64>()
                })
                .sum::<i64>()
        })?;
    for (i, &point) in outer.iter().enumerate() {
        let j = (i + shift) % n;
        guides.push(segment(format!("local-hole-{id}-{j}"), point, inner[j]));
    }
    Some(())
}

/// Insert integral points on a convex host boundary until it has one vertex
/// per vertex of an even convex cut. A host edge without an interior lattice
/// point is skipped; no rounded point is ever introduced.
fn split_ring_to_count(host: &[Point], count: usize) -> Option<Vec<Point>> {
    if count < 4 || !count.is_multiple_of(2) {
        return None;
    }
    let mut ring = host.to_vec();
    while ring.len() < count {
        let candidate = (0..ring.len())
            .filter_map(|i| {
                let a = ring[i];
                let b = ring[(i + 1) % ring.len()];
                let dx = i64::from(b[0] - a[0]);
                let dy = i64::from(b[1] - a[1]);
                let g = gcd(dx.abs(), dy.abs());
                (g > 1).then_some((
                    g,
                    i,
                    [
                        a[0] + (dx / g) as i32 * (g / 2) as i32,
                        a[1] + (dy / g) as i32 * (g / 2) as i32,
                    ],
                ))
            })
            .max_by_key(|(g, i, _)| (*g, std::cmp::Reverse(*i)));
        let (_, edge, point) = candidate?;
        if point == ring[edge] || point == ring[(edge + 1) % ring.len()] {
            return None;
        }
        ring.insert(edge + 1, point);
    }
    Some(ring)
}

fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}

fn convex_even_loop(loop_: &SketchLoop) -> Option<Vec<Point>> {
    if loop_.vertices.len() < 4 {
        return None;
    }
    let mut points = oriented_loop(loop_);
    if !(0..points.len()).all(|i| {
        turn(
            points[i],
            points[(i + 1) % points.len()],
            points[(i + 2) % points.len()],
        ) > 0
    }) {
        return None;
    }
    if !points.len().is_multiple_of(2) {
        let (_, edge, midpoint) = (0..points.len())
            .filter_map(|i| {
                let a = points[i];
                let b = points[(i + 1) % points.len()];
                let x = i64::from(a[0]) + i64::from(b[0]);
                let y = i64::from(a[1]) + i64::from(b[1]);
                (x % 2 == 0 && y % 2 == 0).then_some((
                    i64::from(b[0] - a[0]).pow(2) + i64::from(b[1] - a[1]).pow(2),
                    i,
                    [(x / 2) as i32, (y / 2) as i32],
                ))
            })
            .max_by_key(|(length, edge, _)| (*length, std::cmp::Reverse(*edge)))?;
        points.insert(edge + 1, midpoint);
    }
    Some(points)
}

fn ring_guides(host: &Quad, loop_: &SketchLoop, guides: &mut Vec<Segment>) -> Option<()> {
    let inner = convex_even_loop(loop_)?;
    let max_count = inner.len().saturating_add(16).min(32);
    let mut target = inner.len();
    while target <= max_count {
        let outer = split_ring_to_count(host, target)?;
        let expanded_inner = split_ring_to_count(&inner, target)?;
        let mut trial = Vec::new();
        if join_rings(&outer, &expanded_inner, &loop_.id, &mut trial).is_some() {
            cap_ring(&expanded_inner, &loop_.id, &mut trial);
            guides.extend(trial);
            return Some(());
        }
        target += 2;
    }
    None
}

fn cap_ring(ring: &[Point], id: &str, guides: &mut Vec<Segment>) {
    // Pair consecutive convex boundary edges into quads, sharing a fan of
    // exact diagonals from vertex zero. The arrangement performs final checks.
    for i in (3..ring.len() - 1).step_by(2) {
        guides.push(segment(format!("cap-{id}-{i}"), ring[0], ring[i]));
    }
}

fn line_value(line: [i128; 3], p: Point) -> i128 {
    line[0] * i128::from(p[0]) + line[1] * i128::from(p[1]) + line[2]
}

fn line_edge_intersection(line: [i128; 3], a: Point, b: Point) -> Option<Point> {
    let va = line_value(line, a);
    let vb = line_value(line, b);
    if va == vb || va * vb > 0 {
        return None;
    }
    let den = va - vb;
    let x = i128::from(a[0]) * den + i128::from(b[0] - a[0]) * va;
    let y = i128::from(a[1]) * den + i128::from(b[1] - a[1]) * va;
    if x % den != 0 || y % den != 0 {
        return None;
    }
    Some([i32::try_from(x / den).ok()?, i32::try_from(y / den).ok()?])
}

fn split_quad(host: &Quad, line: [i128; 3]) -> Option<(Quad, Quad)> {
    let mut points = Vec::new();
    for i in 0..4 {
        let a = host[i];
        let b = host[(i + 1) % 4];
        let va = line_value(line, a);
        let vb = line_value(line, b);
        if va == 0 {
            points.push((i, a));
        }
        if va * vb < 0 {
            points.push((i, line_edge_intersection(line, a, b)?));
        }
    }
    points.sort_by_key(|(edge, p)| (*edge, *p));
    points.dedup_by_key(|(_, p)| *p);
    if points.len() != 2 {
        return None;
    }
    let (i, a) = points[0];
    let (j, b) = points[1];
    if (i + 2) % 4 != j && (j + 2) % 4 != i {
        return None;
    }
    let collect = |start: usize, end: usize, first: Point, last: Point| {
        let mut p = vec![first];
        let mut k = (start + 1) % 4;
        while k != (end + 1) % 4 {
            p.push(host[k]);
            k = (k + 1) % 4;
        }
        p.push(last);
        p.try_into().ok()
    };
    let first = collect(i, j, a, b)?;
    let second = collect(j, i, b, a)?;
    (convex(&first) && convex(&second)).then_some((first, second))
}

fn centroid_bisector(first: &[Point], second: &[Point]) -> Option<[i128; 3]> {
    let sum = |ring: &[Point]| -> [i128; 2] {
        ring.iter().fold([0, 0], |mut acc, p| {
            acc[0] += i128::from(p[0]);
            acc[1] += i128::from(p[1]);
            acc
        })
    };
    let (a, b) = (sum(first), sum(second));
    let (na, nb) = (first.len() as i128, second.len() as i128);
    let d = [b[0] * na - a[0] * nb, b[1] * na - a[1] * nb];
    if d == [0, 0] {
        return None;
    }
    // d is parallel to the centroid separation; this is the exact
    // perpendicular bisector after clearing the centroid denominators.
    let c = -(d[0] * (a[0] * nb + b[0] * na) + d[1] * (a[1] * nb + b[1] * na));
    Some([2 * na * nb * d[0], 2 * na * nb * d[1], c])
}

fn chord_line(a: Point, b: Point) -> Option<[i128; 3]> {
    if a == b {
        return None;
    }
    let dx = i128::from(b[0] - a[0]);
    let dy = i128::from(b[1] - a[1]);
    Some([dy, -dx, -dy * i128::from(a[0]) + dx * i128::from(a[1])])
}

fn boundary_samples(host: &Quad, edge: usize) -> Vec<Point> {
    let a = host[edge];
    let b = host[(edge + 1) % 4];
    let mut points = BTreeSet::from([a, b]);
    let dx = i64::from(b[0] - a[0]);
    let dy = i64::from(b[1] - a[1]);
    let g = gcd(dx.abs(), dy.abs());
    for k in 1..16 {
        let step = g * k / 16;
        points.insert([
            (i64::from(a[0]) + dx * k / 16) as i32,
            (i64::from(a[1]) + dy * k / 16) as i32,
        ]);
        if step > 0 {
            points.insert([
                (i64::from(a[0]) + dx * step / g) as i32,
                (i64::from(a[1]) + dy * step / g) as i32,
            ]);
        }
    }
    points.into_iter().collect()
}

fn separated_groups(
    host: &Quad,
    indices: &[usize],
    cuts: &[(SketchLoop, [i32; 2])],
) -> Vec<HostSplit> {
    if indices.len() < 2 {
        return Vec::new();
    }
    let mut proposals = Vec::new();
    for ai in 0..indices.len() {
        for bi in ai + 1..indices.len() {
            let (a, b) = (indices[ai], indices[bi]);
            let (Some(ca), Some(cb)) = (convex_even_loop(&cuts[a].0), convex_even_loop(&cuts[b].0))
            else {
                return Vec::new();
            };
            if let Some(line) = centroid_bisector(&ca, &cb) {
                proposals.push(line);
            }
        }
    }
    // Search opposite-side chords with independently selected boundary
    // parameters. The fixed sampling budget is 16 intervals per host edge;
    // exact arrangement checks decide whether a proposal is usable.
    for first_edge in [0usize, 1] {
        let second_edge = (first_edge + 2) % 4;
        for a in boundary_samples(host, first_edge) {
            for b in boundary_samples(host, second_edge) {
                if let Some(line) = chord_line(a, b) {
                    proposals.push(line);
                }
            }
        }
    }
    let mut candidates: Vec<ScoredHostSplit> = Vec::new();
    for line in proposals {
        let Some((left, right)) = split_quad(host, line) else {
            continue;
        };
        let mut groups = [Vec::new(), Vec::new()];
        let mut valid = true;
        for &index in indices {
            let Some(ring) = convex_even_loop(&cuts[index].0) else {
                valid = false;
                break;
            };
            let signs: Vec<_> = ring.iter().map(|&p| line_value(line, p).signum()).collect();
            let side = if signs.iter().all(|&s| s < 0) {
                0
            } else if signs.iter().all(|&s| s > 0) {
                1
            } else {
                valid = false;
                break;
            };
            groups[side].push(index);
        }
        if !valid || groups.iter().any(Vec::is_empty) {
            continue;
        }
        let left_side = left
            .iter()
            .map(|&p| line_value(line, p))
            .sum::<i128>()
            .signum();
        let right_side = right
            .iter()
            .map(|&p| line_value(line, p))
            .sum::<i128>()
            .signum();
        if left_side == 0 || right_side == 0 || left_side == right_side {
            continue;
        }
        let (left_group, right_group) = if left_side < 0 {
            (groups[0].clone(), groups[1].clone())
        } else {
            (groups[1].clone(), groups[0].clone())
        };
        let divider = host
            .iter()
            .enumerate()
            .find_map(|(i, &p)| line_edge_intersection(line, p, host[(i + 1) % 4]))
            .and_then(|first| {
                host.iter()
                    .enumerate()
                    .find_map(|(i, &p)| {
                        let q = line_edge_intersection(line, p, host[(i + 1) % 4])?;
                        (q != first).then_some(q)
                    })
                    .map(|second| [first, second])
            });
        let Some(divider) = divider else { continue };
        let connects = |host: &Quad, group: &[usize]| {
            group.len() != 1 || {
                if convex_even_loop(&cuts[group[0]].0).is_none() {
                    return false;
                }
                let mut ignored = Vec::new();
                ring_guides(host, &cuts[group[0]].0, &mut ignored).is_some()
            }
        };
        if !connects(&left, &left_group) || !connects(&right, &right_group) {
            continue;
        }
        let balance = left_group.len().abs_diff(right_group.len());
        let score = (balance, divider[0], divider[1]);
        candidates.push((score, left, left_group, right, right_group, divider));
    }
    candidates.sort_by_key(|(score, ..)| *score);
    candidates
        .into_iter()
        .map(|(_, left, il, right, ir, divider)| (left, il, right, ir, divider))
        .collect()
}

fn partition_hosts(
    host: Quad,
    indices: &[usize],
    cuts: &[(SketchLoop, [i32; 2])],
    out: &mut Vec<(Quad, usize)>,
    guides: &mut Vec<Segment>,
    budget: &mut usize,
) -> Option<()> {
    if *budget == 0 {
        return None;
    }
    *budget -= 1;
    if indices.len() == 1 {
        out.push((host, indices[0]));
        return Some(());
    }
    for (a, ia, b, ib, divider) in separated_groups(&host, indices, cuts) {
        let mut next_out = out.clone();
        let mut next_guides = guides.clone();
        let name = format!(
            "oblique-group-{}-{}-{}-{}",
            divider[0][0], divider[0][1], divider[1][0], divider[1][1]
        );
        next_guides.push(segment(name, divider[0], divider[1]));
        if partition_hosts(a, &ia, cuts, &mut next_out, &mut next_guides, budget).is_some()
            && partition_hosts(b, &ib, cuts, &mut next_out, &mut next_guides, budget).is_some()
        {
            *out = next_out;
            *guides = next_guides;
            return Some(());
        }
    }
    None
}

/// Propose exact oblique separators between convex cuts in convex quadrilateral
/// hosts. The centroid bisector is retained only when it intersects opposite
/// host edges at integral points and strictly separates complete cut contours.
pub(crate) fn convex_host_guides(job: &SketchLayers) -> Option<Vec<Segment>> {
    if !job.guides.is_empty()
        || job.force_triangulation
        || job.outer.vertices.len() != 4
        || job
            .operations
            .as_ref()
            .is_some_and(|ops| ops.iter().any(|&add| add))
    {
        return None;
    }
    let outer: Quad = oriented_loop(&job.outer).try_into().ok()?;
    if !convex(&outer) {
        return None;
    }
    let mut guides = Vec::new();
    let mut hosts = Vec::new();
    if job.cuts.len() == 1 {
        hosts.push((outer, 0));
    } else {
        let mut budget = 64;
        partition_hosts(
            outer,
            &(0..job.cuts.len()).collect::<Vec<_>>(),
            &job.cuts,
            &mut hosts,
            &mut guides,
            &mut budget,
        )?;
    }
    for (host, index) in hosts {
        ring_guides(&host, &job.cuts[index].0, &mut guides)?;
    }
    Some(guides)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch_quads::guided_quad_arrangement;

    fn contour(id: &str, vertices: &[[i32; 2]]) -> SketchLoop {
        SketchLoop {
            id: id.into(),
            vertices: vertices.to_vec(),
        }
    }

    #[test]
    fn c_section_keeps_separators_local_to_each_hole() {
        let mut job = SketchLayers {
            outer: contour(
                "body",
                &[
                    [0, 0],
                    [20, 0],
                    [20, 4],
                    [6, 4],
                    [6, 16],
                    [20, 16],
                    [20, 20],
                    [0, 20],
                ],
            ),
            cuts: vec![
                (
                    contour("lower", &[[12, 1], [14, 1], [14, 3], [12, 3]]),
                    [0, 7],
                ),
                (
                    contour("upper", &[[15, 17], [17, 17], [17, 19], [15, 19]]),
                    [0, 7],
                ),
            ],
            guides: vec![],
            height: 7,
            shift: [3, -2],
            force_triangulation: false,
            operations: None,
        };
        assert!(guided_quad_arrangement(&job).unwrap().is_none());
        let guides = local_quad_guides(&job).unwrap();
        assert_eq!(guides.len(), 10);
        job.guides = guides.clone();
        let graph = guided_quad_arrangement(&job).unwrap().unwrap();
        assert_eq!(graph.faces.len(), 11); // nine stock cells and two through holes
        job.guides.clear();
        job.outer.vertices.reverse();
        job.cuts[0].0.vertices.rotate_left(2);
        assert_eq!(local_quad_guides(&job).unwrap(), guides);
        job.cuts[0].0.vertices.iter_mut().for_each(|p| p[0] -= 1);
        let stationary = |list: Vec<Segment>| {
            list.into_iter()
                .filter(|g| !g.id.starts_with("local-hole-lower-"))
                .collect::<Vec<_>>()
        };
        assert_eq!(
            stationary(local_quad_guides(&job).unwrap()),
            stationary(guides)
        );
    }

    #[test]
    fn diagonal_user_guides_preserve_only_conforming_convex_quads() {
        let body = contour("body", &[[0, 0], [60, 0], [60, 40], [0, 40]]);
        let pocket = contour("pocket", &[[20, 15], [40, 15], [40, 25], [20, 25]]);
        let guides = body
            .vertices
            .iter()
            .zip(&pocket.vertices)
            .enumerate()
            .map(|(i, (&a, &b))| segment(format!("g{i}"), a, b))
            .collect();
        let mut job = SketchLayers {
            outer: body,
            cuts: vec![(pocket, [14, 20])],
            guides,
            height: 20,
            shift: [0, 0],
            force_triangulation: false,
            operations: None,
        };
        assert_eq!(
            guided_quad_arrangement(&job).unwrap().unwrap().faces.len(),
            5
        );
        job.guides.pop();
        assert!(guided_quad_arrangement(&job).unwrap().is_none());
    }

    #[test]
    fn arbitrary_even_convex_ring_is_a_complete_quad_host() {
        let job = SketchLayers {
            outer: contour("stock", &[[0, 0], [100, 0], [100, 80], [0, 80]]),
            cuts: vec![(
                contour(
                    "hex",
                    &[[30, 30], [40, 25], [50, 30], [50, 40], [40, 45], [30, 40]],
                ),
                [0, 10],
            )],
            guides: vec![],
            height: 10,
            shift: [0, 0],
            force_triangulation: false,
            operations: None,
        };
        let guides = convex_host_guides(&job).unwrap();
        assert_eq!(
            guides
                .iter()
                .filter(|g| g.id.starts_with("local-hole-hex-"))
                .count(),
            6
        );
        let mut proposed = job.clone();
        proposed.guides = guides;
        assert_eq!(
            guided_quad_arrangement(&proposed)
                .unwrap()
                .unwrap()
                .faces
                .len(),
            8
        );
    }

    #[test]
    fn odd_convex_ring_is_evenized_by_an_exact_edge_split() {
        let mut job = SketchLayers {
            outer: contour("stock", &[[0, 0], [100, 0], [100, 100], [0, 100]]),
            cuts: vec![(
                contour("pent", &[[30, 30], [40, 25], [50, 30], [45, 45], [32, 43]]),
                [0, 10],
            )],
            guides: vec![],
            height: 10,
            shift: [0, 0],
            force_triangulation: false,
            operations: None,
        };
        for p in &mut job.outer.vertices {
            p.iter_mut().for_each(|v| *v *= 2);
        }
        for (loop_, _) in &mut job.cuts {
            for p in &mut loop_.vertices {
                p.iter_mut().for_each(|v| *v *= 2);
            }
        }
        let guides = convex_host_guides(&job).unwrap();
        assert_eq!(
            guides
                .iter()
                .filter(|g| g.id.starts_with("local-hole-pent-"))
                .count(),
            6
        );
        job.guides = guides;
        assert_eq!(
            guided_quad_arrangement(&job).unwrap().unwrap().faces.len(),
            8
        );
    }

    #[test]
    fn independent_opposite_edge_parameters_find_an_oblique_hole_separator() {
        let job = SketchLayers {
            outer: contour("stock", &[[0, 0], [100, 0], [100, 100], [0, 100]]),
            cuts: vec![
                (
                    contour("lower_left", &[[21, 36], [29, 36], [29, 44], [21, 44]]),
                    [0, 10],
                ),
                (
                    contour("upper_right", &[[71, 56], [79, 56], [79, 64], [71, 64]]),
                    [0, 10],
                ),
            ],
            guides: vec![],
            height: 10,
            shift: [0, 0],
            force_triangulation: false,
            operations: None,
        };
        let host: Quad = oriented_loop(&job.outer).try_into().unwrap();
        assert!(!separated_groups(&host, &[0, 1], &job.cuts).is_empty());
        let guides = convex_host_guides(&job).unwrap();
        let diagonal = guides
            .iter()
            .find(|g| g.id.starts_with("oblique-group-"))
            .unwrap();
        assert_ne!(diagonal.a[0], diagonal.b[0]);
        assert_ne!(diagonal.a[1], diagonal.b[1]);
        let mut proposed = job.clone();
        proposed.guides = guides;
        assert_eq!(
            guided_quad_arrangement(&proposed)
                .unwrap()
                .unwrap()
                .faces
                .len(),
            10
        );
    }
}
