//! Propose local diagonals in a simple polygon, then join each quadrangular
//! cut to the corners of its enclosing cell. Every proposal still goes
//! through the arrangement, coverage, trace and Jacobian checks.
use crate::sketch_arrangement::{Role, Segment};
use crate::sketch_layers::{edges, orient as turn, signed_area2, SketchLayers, SketchLoop};

type Point = [i32; 2];
type Quad = [Point; 4];

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
        let shift = (0..4)
            .filter(|shift| {
                (0..4).all(|i| {
                    let j = (i + 1) % 4;
                    convex(&[
                        outer[i],
                        outer[j],
                        inner[(j + shift) % 4],
                        inner[(i + shift) % 4],
                    ])
                })
            })
            .min_by_key(|shift| {
                (0..4)
                    .map(|i| {
                        (0..2)
                            .map(|axis| {
                                i64::from(outer[i][axis] - inner[(i + shift) % 4][axis]).pow(2)
                            })
                            .sum::<i64>()
                    })
                    .sum::<i64>()
            })?;
        for (i, &p) in outer.iter().enumerate() {
            let j = (i + shift) % 4;
            guides.push(segment(format!("local-hole-{}-{j}", hole.id), p, inner[j]));
        }
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
}
