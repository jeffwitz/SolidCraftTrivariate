//! From sketch contours to a conforming rational quadrilateral footprint.
use crate::{
    constraint_system::Ratio,
    sketch_arrangement::{Arrangement, FaceClass},
    sketch_convex_hosts::host_candidates,
    sketch_convex_rings::plan_convex_rings,
    sketch_layers::{orient, ratio, SketchLayers, SketchLoop},
    sketch_quads::{centroid, cross, midpoint},
};
type Point = [Ratio; 2];
type Quad = [Point; 4];
type HostPlan = Vec<Vec<Point>>;
type HostPlanner = fn(&[Point], &[Vec<Point>]) -> Vec<HostPlan>;

fn ccw(contour: &SketchLoop) -> Option<Vec<Point>> {
    let v = &contour.vertices;
    let n = v.len();
    let sign = orient(v[0], v[1], v[2]).signum();
    if sign == 0 || !(0..n).all(|i| orient(v[i], v[(i + 1) % n], v[(i + 2) % n]).signum() == sign) {
        return None;
    }
    let mut points: Vec<_> = v.iter().map(|p| p.map(ratio)).collect();
    if sign < 0 {
        points.reverse()
    }
    let first = points.iter().enumerate().min_by_key(|(_, p)| *p)?.0;
    points.rotate_left(first);
    Some(points)
}

// Choose a strictly convex ladder across the footprint of a blind cut.
// Rank rotations by the largest corner Jacobian distortion.
fn cap(poly: &[Point]) -> Option<Vec<Quad>> {
    let n = poly.len();
    let mut best: Option<(f64, Vec<Quad>)> = None;
    for start in 0..n / 2 {
        let p = |i: usize| poly[(start + i) % n];
        let quads: Vec<Quad> = (0..n / 2 - 1)
            .map(|i| [p(i), p(i + 1), p(n - 2 - i), p(n - 1 - i)])
            .collect();
        let mut worst = 0.0_f64;
        for q in &quads {
            for i in 0..4 {
                let [a, b, d] = [q[i], q[(i + 1) % 4], q[(i + 3) % 4]];
                let det = cross(a, b, d).ok()?;
                if det.numerator() <= 0 {
                    worst = f64::INFINITY;
                    break;
                }
                let [a, b, d] =
                    [a, b, d].map(|p| p.map(|v| v.numerator() as f64 / v.denominator() as f64));
                let norm: f64 = (0..2)
                    .map(|k| (b[k] - a[k]).powi(2) + (d[k] - a[k]).powi(2))
                    .sum();
                worst = worst.max(norm / (det.numerator() as f64 / det.denominator() as f64));
            }
        }
        if worst.is_finite() && best.as_ref().is_none_or(|(score, _)| worst < *score) {
            best = Some((worst, quads))
        }
    }
    best.map(|(_, quads)| quads)
}

pub(super) fn visit(job: &SketchLayers, mut consider: impl FnMut(Arrangement)) {
    visit_with(job, host_candidates, &mut consider)
}

pub(super) fn visit_lookahead(job: &SketchLayers, mut consider: impl FnMut(Arrangement)) {
    visit_with(
        job,
        crate::sketch_convex_hosts::host_candidates_lookahead,
        &mut consider,
    )
}

fn visit_with(
    job: &SketchLayers,
    host_candidates: HostPlanner,
    consider: &mut impl FnMut(Arrangement),
) {
    if job.force_triangulation
        || !job.guides.is_empty()
        || job.operations.is_some()
        || job.cuts.is_empty()
    {
        return;
    }
    let Some(stock) = ccw(&job.outer) else { return };
    let Some(mut holes) = job
        .cuts
        .iter()
        .map(|(l, _)| ccw(l))
        .collect::<Option<Vec<_>>>()
    else {
        return;
    };
    let mut caps = Vec::new();
    for (hole, (_, interval)) in holes.iter_mut().zip(&job.cuts) {
        if *interval == [0, job.height] {
            continue;
        }
        // A triangular floor needs three genuine corners per boundary cycle
        // and cannot be one regular quadrilateral. Use the same edge midpoints
        // for the three floor cells and for the adjacent ring.
        if hole.len() == 3 {
            let [a, b, c] = [hole[0], hole[1], hole[2]];
            let (Ok(ab), Ok(bc), Ok(ca), Ok(g)) = (
                midpoint(a, b),
                midpoint(b, c),
                midpoint(c, a),
                centroid([a, b, c]),
            ) else {
                return;
            };
            caps.extend(
                [[a, ab, g, ca], [b, bc, g, ab], [c, ca, g, bc]].map(|q| (q, FaceClass::Void)),
            );
            *hole = vec![a, ab, b, bc, c, ca];
            continue;
        }
        if hole.len() % 2 == 1 {
            let n = hole.len();
            let edge = (0..n)
                .max_by_key(|&i| {
                    let [a, b] = [hole[i], hole[(i + 1) % n]];
                    (0..2)
                        .map(|j| (a[j].numerator() - b[j].numerator()).pow(2))
                        .sum::<i128>()
                })
                .unwrap();
            let Ok(m) = midpoint(hole[edge], hole[(edge + 1) % n]) else {
                return;
            };
            hole.insert(edge + 1, m);
        }
        let Some(quads) = cap(hole) else { return };
        caps.extend(quads.into_iter().map(|q| (q, FaceClass::Void)));
    }
    for hosts in host_candidates(&stock, &holes) {
        let mut submit = |rings: Vec<Quad>| {
            let mut cells: Vec<_> = rings
                .into_iter()
                .map(|q| (q, FaceClass::Material))
                .collect();
            cells.extend_from_slice(&caps);
            if let Ok(graph) = crate::sketch_radial::graph(job, cells, "convex_ring") {
                consider(graph)
            }
        };
        if let Some(rings) = plan_convex_rings(&hosts, &holes) {
            submit(rings)
        }
        let mut order: Vec<_> = (0..holes.len()).collect();
        if let Some(rings) = crate::sketch_optional_ring::rings(&hosts, &holes, &order) {
            submit(rings)
        }
        if order.len() > 1 {
            order.reverse();
            if let Some(rings) = crate::sketch_optional_ring::rings(&hosts, &holes, &order) {
                submit(rings)
            }
        }
    }
}
