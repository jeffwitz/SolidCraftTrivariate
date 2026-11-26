//! A rational fan for one convex footprint inside convex stock.
//! Rays through both sets of corners keep every prescribed segment on the mesh.
use crate::constraint_system::Ratio;
use crate::sketch_arrangement::{
    Arrangement, ArrangementEdge, Face, FaceClass, FaceEdge, Provenance, Role,
};
use crate::sketch_layers::{edges, orient, SketchLayers, SketchLoop};
use crate::sketch_quads::{centroid, cross, midpoint};
use std::collections::BTreeMap;

type Point = [Ratio; 2];
type Vector = [i128; 2];
type Quad = [Point; 4];
type Cell = (Quad, FaceClass);

fn det(a: Vector, b: Vector) -> i128 {
    a[0] * b[1] - a[1] * b[0]
}
fn turn_sign(l: &SketchLoop) -> Option<i128> {
    let v = &l.vertices;
    let n = v.len();
    let s = orient(v[0], v[1], v[2]).signum();
    (s != 0 && (0..n).all(|i| orient(v[i], v[(i + 1) % n], v[(i + 2) % n]).signum() == s))
        .then_some(i128::from(s))
}
fn exit(l: &SketchLoop, sign: i128, c: Vector, n: i128, d: Vector) -> Result<Point, &'static str> {
    let mut limit: Option<Ratio> = None;
    for (a, b) in edges(l) {
        let e = [
            i128::from(b[0]) - i128::from(a[0]),
            i128::from(b[1]) - i128::from(a[1]),
        ];
        let dist = sign
            * det(
                e,
                [c[0] - n * i128::from(a[0]), c[1] - n * i128::from(a[1])],
            );
        if dist <= 0 {
            return Err("radial_origin_outside");
        }
        let velocity = -sign * det(e, d);
        if velocity > 0 {
            let t = Ratio::new(dist, velocity)?;
            if limit.is_none_or(|x| t.cmp_value(&x).is_lt()) {
                limit = Some(t);
            }
        }
    }
    let t = limit.ok_or("open_radial_loop")?;
    let coordinate = |k: usize| {
        Ratio::new(c[k], 1)?
            .checked_add(t.checked_mul(Ratio::new(d[k], 1)?)?)?
            .checked_mul(Ratio::new(1, n)?)
    };
    Ok([coordinate(0)?, coordinate(1)?])
}
fn convex(q: Quad) -> Result<bool, &'static str> {
    for i in 0..4 {
        if cross(q[i], q[(i + 1) % 4], q[(i + 2) % 4])?.numerator() <= 0 {
            return Ok(false);
        }
    }
    Ok(true)
}
fn cells(
    mut rays: Vec<(Point, Point)>,
    center: Point,
    through: bool,
) -> Result<Vec<Cell>, &'static str> {
    use FaceClass::{Material, Void};
    if !through && rays.len() % 2 == 1 {
        let (a, b) = rays[0];
        let (c, d) = rays[1];
        rays.insert(1, (midpoint(a, c)?, midpoint(b, d)?));
    }
    let n = rays.len();
    let ring: Vec<_> = (0..n)
        .map(|i| {
            let (a, b) = rays[i];
            let (d, c) = rays[(i + 1) % n];
            ([a, b, c, d], Material)
        })
        .collect();
    if through {
        return Ok(ring);
    }
    for offset in 0..2 {
        let mut out = ring.clone();
        for i in (offset..n).step_by(2) {
            let q = [center, rays[i].0, rays[(i + 1) % n].0, rays[(i + 2) % n].0];
            if !convex(q)? {
                break;
            }
            out.push((q, Void));
        }
        if out.len() == n + n / 2 {
            return Ok(out);
        }
    }
    let mut out = Vec::new();
    for ([a, b, c, d], _) in ring {
        let m = midpoint(a, d)?;
        let k = midpoint(b, c)?;
        let ca = midpoint(center, a)?;
        let cd = midpoint(center, d)?;
        let g = centroid([center, a, d])?;
        out.extend([
            ([a, b, k, m], Material),
            ([m, k, c, d], Material),
            ([center, ca, g, cd], Void),
            ([a, m, g, ca], Void),
            ([d, cd, g, m], Void),
        ]);
    }
    Ok(out)
}
pub(crate) fn graph(
    job: &SketchLayers,
    cells: Vec<Cell>,
    guide: &str,
) -> Result<Arrangement, &'static str> {
    let mut result = Arrangement {
        vertices: vec![],
        edges: vec![],
        faces: vec![],
    };
    let mut points = BTreeMap::new();
    let mut edge_ids = BTreeMap::new();
    let mut uses: Vec<[usize; 2]> = vec![];
    for (quad, class) in cells {
        if !convex(quad)? {
            return Err("nonconvex_radial_cell");
        }
        let mut vertices = Vec::new();
        for p in quad {
            let id = *points.entry(p).or_insert_with(|| {
                result.vertices.push(p);
                result.vertices.len() - 1
            });
            vertices.push(id);
        }
        let mut sides = Vec::new();
        for i in 0..4 {
            let (a, b) = (vertices[i], vertices[(i + 1) % 4]);
            let key = [a.min(b), a.max(b)];
            let id = if let Some(id) = edge_ids.get(&key) {
                *id
            } else {
                let [p, q] = key.map(|v| result.vertices[v]);
                let mut source = None;
                for contour in std::iter::once(&job.outer).chain(job.cuts.iter().map(|(l, _)| l)) {
                    for (j, (u, v)) in edges(contour).enumerate() {
                        let [u, v] = [u, v].map(|p| p.map(|x| Ratio::new(x.into(), 1).unwrap()));
                        if cross(u, v, p)?.numerator() == 0 && cross(u, v, q)?.numerator() == 0 {
                            let k = usize::from(u[0] == v[0]);
                            let lo = if u[k].cmp_value(&v[k]).is_lt() {
                                u[k]
                            } else {
                                v[k]
                            };
                            let hi = if lo == u[k] { v[k] } else { u[k] };
                            if [p[k], q[k]]
                                .into_iter()
                                .any(|x| x.cmp_value(&lo).is_lt() || x.cmp_value(&hi).is_gt())
                            {
                                continue;
                            }
                            source = Some(Provenance {
                                segment_id: format!("boundary/{}/{j}", contour.id),
                                role: Role::Boundary,
                                forward: p[k].cmp_value(&q[k]) == u[k].cmp_value(&v[k]),
                            });
                            break;
                        }
                    }
                    if source.is_some() {
                        break;
                    }
                }
                let id = result.edges.len();
                result.edges.push(ArrangementEdge {
                    vertices: key,
                    provenance: vec![source.unwrap_or_else(|| Provenance {
                        segment_id: format!("guide/{guide}/{id}"),
                        role: Role::Construction,
                        forward: true,
                    })],
                });
                uses.push([0; 2]);
                edge_ids.insert(key, id);
                id
            };
            let forward = a < b;
            uses[id][usize::from(forward)] += 1;
            sides.push(FaceEdge { edge: id, forward });
        }
        result.faces.push(Face {
            vertices,
            edges: sides,
            class,
        });
    }
    let stock = format!("boundary/{}/", job.outer.id);
    let through: Vec<_> = job
        .cuts
        .iter()
        .filter(|(_, z)| *z == [0, job.height])
        .map(|(l, _)| format!("boundary/{}/", l.id))
        .collect();
    for (edge, directions) in result.edges.iter().zip(uses) {
        let exterior = edge.provenance[0].segment_id.starts_with(&stock)
            || through
                .iter()
                .any(|prefix| edge.provenance[0].segment_id.starts_with(prefix));
        if (exterior && directions.iter().sum::<usize>() != 1)
            || (!exterior && directions != [1, 1])
        {
            return Err("radial_edge_incidence");
        }
    }
    let two = Ratio::new(2, 1)?;
    for p in &mut result.vertices {
        for x in p {
            *x = x.checked_mul(two)?
        }
    }
    Ok(result)
}
pub(crate) fn footprint(job: &SketchLayers) -> Result<Option<Arrangement>, &'static str> {
    if job.force_triangulation
        || !job.guides.is_empty()
        || job.operations.as_ref().is_some_and(|ops| ops != &[false])
        || job.cuts.len() != 1
    {
        return Ok(None);
    }
    let inner = &job.cuts[0].0;
    if inner.vertices.len() < 3
        || job.outer.vertices.len() < 3
        || inner.vertices.len() + job.outer.vertices.len() > 48
    {
        return Ok(None);
    }
    let (Some(si), Some(so)) = (turn_sign(inner), turn_sign(&job.outer)) else {
        return Ok(None);
    };
    let n = inner.vertices.len() as i128;
    let c = inner.vertices.iter().fold([0; 2], |v, p| {
        [v[0] + i128::from(p[0]), v[1] + i128::from(p[1])]
    });
    let mut dirs: Vec<_> = inner
        .vertices
        .iter()
        .chain(&job.outer.vertices)
        .map(|p| [n * i128::from(p[0]) - c[0], n * i128::from(p[1]) - c[1]])
        .collect();
    let half = |v: Vector| v[1] < 0 || (v[1] == 0 && v[0] < 0);
    dirs.sort_by(|&a, &b| half(a).cmp(&half(b)).then_with(|| det(b, a).cmp(&0)));
    dirs.dedup_by(|a, b| half(*a) == half(*b) && det(*a, *b) == 0);
    let mut rays = vec![];
    for d in dirs {
        rays.push((exit(inner, si, c, n, d)?, exit(&job.outer, so, c, n, d)?));
    }
    let center = [Ratio::new(c[0], n)?, Ratio::new(c[1], n)?];
    graph(
        job,
        cells(rays, center, job.cuts[0].1 == [0, job.height])?,
        "radial",
    )
    .map(Some)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rational_exits_keep_both_contours_and_cap_a_blind_odd_loop() {
        let mut job = SketchLayers {
            outer: SketchLoop {
                id: "stock".into(),
                vertices: vec![[0, 0], [42, 1], [46, 24], [30, 43], [1, 31]],
            },
            cuts: vec![(
                SketchLoop {
                    id: "cut".into(),
                    vertices: vec![[10, 12], [14, 22], [25, 15]],
                },
                [0, 10],
            )],
            height: 10,
            guides: vec![],
            operations: None,
            shift: [0, 0],
            strategy: crate::sketch_layers::PlannerStrategy::Automatic,
            force_triangulation: false,
        };
        let ring = footprint(&job).unwrap().unwrap();
        assert!(ring.vertices.iter().flatten().any(|v| v.denominator() > 1));
        assert_eq!(ring.faces.len(), 8);
        for (contour, n) in [("stock", 5), ("cut", 3)] {
            for i in 0..n {
                assert!(ring
                    .edges
                    .iter()
                    .flat_map(|e| &e.provenance)
                    .any(|p| p.segment_id == format!("boundary/{contour}/{i}")));
            }
        }
        job.cuts[0].1 = [3, 10];
        let cap = footprint(&job).unwrap().unwrap();
        assert!(cap.faces.iter().any(|f| f.class == FaceClass::Void));
        job.outer.vertices.reverse();
        assert_eq!(footprint(&job).unwrap().unwrap().vertices, cap.vertices);
    }
}
