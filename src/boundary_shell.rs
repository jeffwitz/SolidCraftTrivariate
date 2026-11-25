//! Exact boundary incidence and closed-shell components of certified whole-face blocks.
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
    /// Indices into faces, with the material exterior first and each void after it.
    pub components: Vec<Vec<usize>>,
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
    while net.len() < 4 {
        let n = net.len();
        net = (0..=n)
            .map(|i| {
                std::array::from_fn(|c| {
                    i as i128 * net[i.saturating_sub(1)][c] + (n - i) as i128 * net[i.min(n - 1)][c]
                })
            })
            .collect();
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

fn reachable(start: usize, adjacent: &BTreeMap<usize, Vec<usize>>) -> BTreeSet<usize> {
    let mut visited = BTreeSet::new();
    let mut pending = vec![start];
    while let Some(v) = pending.pop() {
        if visited.insert(v) {
            pending.extend(&adjacent[&v]);
        }
    }
    visited
}

pub(crate) fn supports_source_method(method: &str) -> bool {
    matches!(
        method,
        "certified_fixed_extrusion_partition"
            | "certified_ruled_similarity_sweep"
            | "certified_coaxial_feature_program"
            | "certified_circular_bend"
            | "certified_orthogonal_feature_program"
            | "certified_polygon_prism"
            | "certified_perforated_polygon_prism"
            | crate::rational_sector_json::METHOD
            | crate::counterbores_json::METHOD
            | crate::countersinks_json::METHOD
            | crate::sketch_layers_json::METHOD
            | crate::nurbs_multiloop_json::METHOD
            | crate::nurbs_sketch_json::METHOD
            | crate::orthogonal_sketch_json::METHOD
            | crate::interval_bores_json::METHOD
            | crate::circular_polygon_json::METHOD
            | crate::nurbs_cross_json::METHOD
            | crate::nurbs_cross_json::BOUNDED_METHOD
            | crate::nurbs_cross_json::OVAL_METHOD
    )
}

/// Load a canonically audited model and construct its exterior shell, proving
/// exact trace identity before identifying two curves with common endpoints.
pub fn boundary_shell(model: &Json) -> Result<BoundaryShell, &'static str> {
    let report = audit(model)?;
    let method = report.field("method")?.text()?;
    if !supports_source_method(method) && method != crate::affine_image::METHOD {
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
    let expected_voids = if matches!(
        method,
        crate::affine_image::METHOD
            | crate::counterbores_json::METHOD
            | crate::countersinks_json::METHOD
            | crate::sketch_layers_json::METHOD
            | crate::nurbs_multiloop_json::METHOD
            | crate::nurbs_sketch_json::METHOD
            | crate::orthogonal_sketch_json::METHOD
            | crate::interval_bores_json::METHOD
            | "certified_coaxial_feature_program"
            | "certified_circular_bend"
            | "certified_orthogonal_feature_program"
    ) {
        report.field("cavities")?.number::<usize>()?
    } else {
        0
    };
    let mut components = Vec::new();
    let mut pending: BTreeSet<_> = (0..faces.len()).collect();
    while let Some(&seed) = pending.first() {
        let connected = reachable(seed, &adjacent);
        for i in &connected {
            pending.remove(i);
        }
        components.push(connected.into_iter().collect::<Vec<_>>());
    }
    if components.len() != expected_voids + 1 {
        return Err("unexpected_boundary_shell_count");
    }
    if expected_voids > 0 {
        // A coaxial plate retains its x=10 wall. For a periodic bend the
        // outer section curve is certified by its independent section audit.
        let bend_anchor = if method == "certified_circular_bend" {
            match model.field("section")?.text()? {
                "full_annulus_v1" => Some("section_circle_r4"),
                "pierced_square_v1" => Some("section_x10"),
                _ => return Err("unexpected_bend_void"),
            }
        } else {
            None
        };
        let anchor_x = if method == "certified_orthogonal_feature_program" {
            i128::from(report.field("minimum_x_mm")?.number::<i32>()?)
        } else {
            10
        };
        // Design face order is checked against the source during affine audit.
        let source_anchor = if method == crate::affine_image::METHOD {
            let original = boundary_shell(model.field("source_model")?)?;
            Some(original.components[0][0])
        } else {
            None
        };
        let mut outer = components.iter().enumerate().filter(|(_, component)| {
            component.iter().any(|&f| {
                if let Some(anchor) = source_anchor {
                    return f == anchor;
                }
                if matches!(
                    method,
                    crate::interval_bores_json::METHOD | crate::sketch_layers_json::METHOD
                ) {
                    return faces[f].feature == "start";
                }
                if method == crate::nurbs_multiloop_json::METHOD {
                    return faces[f].feature == "stock:xmin";
                }
                if let Some(tag) = bend_anchor {
                    return faces[f].feature == tag;
                }
                faces[f]
                    .controls
                    .iter()
                    .all(|p| i128::from(p[0]) == anchor_x * i128::from(p[3]))
            })
        });
        let index = outer.next().ok_or("outer_shell_missing")?.0;
        if outer.next().is_some() {
            return Err("ambiguous_outer_shell");
        }
        components.swap(0, index);
    }
    for link in links.values() {
        if link.values().any(|v| v.len() != 2)
            || reachable(*link.keys().next().unwrap(), link).len() != link.len()
        {
            return Err("nonmanifold_boundary_vertex");
        }
    }
    Ok(BoundaryShell {
        vertices,
        edges,
        faces,
        components,
    })
}
