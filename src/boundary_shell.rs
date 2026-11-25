//! Exact boundary incidence for a single closed shell of certified whole-face blocks.
use crate::{bracket::BracketBlock, bracket_json::decode_blocks, interchange::audit, json::Json};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug)]
pub struct ShellEdge {
    pub vertices: [usize; 2],
    pub controls: Vec<[i64; 4]>,
}
#[derive(Debug)]
pub struct ShellFace {
    pub feature: String,
    pub block: String,
    pub side: String,
    pub degrees: [usize; 2],
    pub controls: Vec<[i64; 4]>,
    /// In exterior winding order; the boolean denotes traversal in edge direction.
    pub loop_edges: Vec<(usize, bool)>,
    pub same_sense: bool,
}
#[derive(Debug)]
pub struct BoundaryShell {
    pub vertices: Vec<[i128; 4]>,
    pub edges: Vec<ShellEdge>,
    pub faces: Vec<ShellFace>,
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}
fn vertex(p: [i64; 4]) -> [i128; 4] {
    let p = p.map(i128::from);
    let g = p.into_iter().fold(0, gcd);
    p.map(|x| x / g)
}
fn canonical_edge(points: &[[i64; 4]]) -> Vec<[i128; 4]> {
    let mut net: Vec<_> = points.iter().map(|p| p.map(i128::from)).collect();
    if let [a, b] = net.as_slice() {
        net = vec![
            a.map(|x| 2 * x),
            std::array::from_fn(|i| a[i] + b[i]),
            b.map(|x| 2 * x),
        ];
    }
    let scale = net.iter().flatten().copied().fold(0, gcd);
    for p in &mut net {
        for h in p {
            *h /= scale;
        }
    }
    net
}

fn extract(block: &BracketBlock, axis: usize, upper: bool) -> ([usize; 2], Vec<[i64; 4]>) {
    let other = match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    };
    let degrees = other.map(|a| block.degrees[a]);
    let mut points = Vec::new();
    for i in 0..=degrees[0] {
        for j in 0..=degrees[1] {
            let mut q = [0; 3];
            q[axis] = usize::from(upper) * block.degrees[axis];
            q[other[0]] = i;
            q[other[1]] = j;
            points.push(
                block.controls
                    [(q[0] * (block.degrees[1] + 1) + q[1]) * (block.degrees[2] + 1) + q[2]],
            );
        }
    }
    (degrees, points)
}

fn reachable(start: usize, adjacent: &BTreeMap<usize, Vec<usize>>) -> usize {
    let mut visited = BTreeSet::new();
    let mut pending = vec![start];
    while let Some(v) = pending.pop() {
        if visited.insert(v) {
            pending.extend(&adjacent[&v]);
        }
    }
    visited.len()
}

/// Load a canonically audited model and construct its exterior shell, proving
/// exact trace identity before identifying two curves with common endpoints.
pub fn boundary_shell(model: &Json) -> Result<BoundaryShell, &'static str> {
    let report = audit(model)?;
    if report.field("method")?.text()? != "certified_fixed_extrusion_partition" {
        return Err("unsupported_boundary_model_family");
    }
    let blocks: BTreeMap<_, _> = decode_blocks(model)?
        .into_iter()
        .map(|b| (b.id.clone(), b))
        .collect();
    let mut vertices = Vec::new();
    let mut vertex_index = BTreeMap::new();
    let mut edges: Vec<ShellEdge> = Vec::new();
    let mut edge_index = BTreeMap::new();
    let mut uses: Vec<Vec<(usize, bool)>> = Vec::new();
    let mut faces = Vec::new();
    let mut links = BTreeMap::<usize, BTreeMap<usize, Vec<usize>>>::new();
    for face in model.field("design_faces")?.items()? {
        let id = face.field("block")?.text()?;
        let side = face.field("side")?.text()?;
        let (axis, upper) = match side {
            "u0" => (0, false),
            "u1" => (0, true),
            "v0" => (1, false),
            "v1" => (1, true),
            "w0" => (2, false),
            "w1" => (2, true),
            _ => return Err("unknown_boundary_side"),
        };
        let (degrees, controls) =
            extract(blocks.get(id).ok_or("missing_boundary_block")?, axis, upper);
        let [p, q] = degrees;
        let at = |i, j| controls[i * (q + 1) + j];
        let mut ring: Vec<Vec<_>> = vec![
            (0..=p).map(|i| at(i, 0)).collect(),
            (0..=q).map(|j| at(p, j)).collect(),
            (0..=p).rev().map(|i| at(i, q)).collect(),
            (0..=q).rev().map(|j| at(0, j)).collect(),
        ];
        let same_sense = upper ^ (axis == 1);
        if !same_sense {
            ring.reverse();
            for curve in &mut ring {
                curve.reverse();
            }
        }
        let mut loop_edges = Vec::new();
        let mut corners = Vec::new();
        for curve in ring {
            let endpoints = [curve[0], *curve.last().unwrap()].map(|h| {
                let key = vertex(h);
                *vertex_index.entry(key).or_insert_with(|| {
                    vertices.push(key);
                    vertices.len() - 1
                })
            });
            let [a, b] = endpoints;
            if a == b {
                return Err("collapsed_boundary_edge");
            }
            let key = [a.min(b), a.max(b)];
            let e = *edge_index.entry(key).or_insert_with(|| {
                edges.push(ShellEdge {
                    vertices: endpoints,
                    controls: curve.clone(),
                });
                uses.push(vec![]);
                edges.len() - 1
            });
            let forward = edges[e].vertices == endpoints;
            let mut net = canonical_edge(&curve);
            if !forward {
                net.reverse();
            }
            if net != canonical_edge(&edges[e].controls) {
                return Err("incompatible_boundary_edge_trace");
            }
            uses[e].push((faces.len(), forward));
            corners.push(a);
            loop_edges.push((e, forward));
        }
        for (j, v) in corners.into_iter().enumerate() {
            let a = loop_edges[(j + 3) % 4].0;
            let b = loop_edges[j].0;
            let graph = links.entry(v).or_default();
            graph.entry(a).or_default().push(b);
            graph.entry(b).or_default().push(a);
        }
        faces.push(ShellFace {
            feature: face.field("feature")?.text()?.into(),
            block: id.into(),
            side: side.into(),
            degrees,
            controls,
            loop_edges,
            same_sense,
        });
    }
    if faces.is_empty() {
        return Err("empty_boundary_shell");
    }
    let mut adjacent = BTreeMap::<usize, Vec<usize>>::new();
    for entries in uses {
        let [(a, first), (b, second)] = entries.as_slice() else {
            return Err("nonmanifold_boundary_edge");
        };
        if a == b || first == second {
            return Err("inconsistent_boundary_winding");
        }
        adjacent.entry(*a).or_default().push(*b);
        adjacent.entry(*b).or_default().push(*a);
    }
    if reachable(0, &adjacent) != faces.len() {
        return Err("multiple_boundary_shells");
    }
    for link in links.values() {
        if link.values().any(|v| v.len() != 2)
            || reachable(*link.keys().next().unwrap(), link) != link.len()
        {
            return Err("nonmanifold_boundary_vertex");
        }
    }
    Ok(BoundaryShell {
        vertices,
        edges,
        faces,
    })
}
