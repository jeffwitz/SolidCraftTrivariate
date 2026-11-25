//! Edge and vertex incidence for a whole-face, conforming triquadratic complex.
//! Global analysis indices carry topology established by exact oriented faces.
//! Two coincident points with unrelated topological indices remain independent.
use crate::{
    analysis_space::AnalysisSpace,
    json::{object, Json},
};
use std::collections::BTreeMap;

fn number(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn numbers<T: ToString>(values: impl IntoIterator<Item = T>) -> Json {
    Json::Array(values.into_iter().map(number).collect())
}
fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}
fn normalized<const N: usize>(points: [[i128; 4]; N]) -> Result<[[i128; 4]; N], &'static str> {
    if points.iter().any(|p| p[3] <= 0 || p.contains(&i128::MIN)) {
        return Err("invalid skeleton homogeneous coefficients");
    }
    let common = points.iter().flatten().copied().fold(0, gcd);
    Ok(points.map(|p| p.map(|c| c / common)))
}
type EdgeUse = (String, usize, [usize; 2], bool);
type VertexUse = (String, [usize; 3]);
type Edge = ([[i128; 4]; 3], Vec<EdgeUse>);
type Vertex = ([i128; 4], Vec<VertexUse>);

/// Encode every edge and vertex, including embedded interior entities.
/// The local axis varies increasingly, and fixed-side flags follow the other
/// two axes in increasing order. Edge orientation uses increasing endpoint ID.
pub fn encode_skeleton(space: &AnalysisSpace) -> Result<Json, &'static str> {
    let mut edges: BTreeMap<[usize; 3], Edge> = BTreeMap::new();
    let mut vertices: BTreeMap<usize, Vertex> = BTreeMap::new();
    for block in &space.blocks {
        if block.dofs.len() != 27
            || block.controls.len() != 27
            || block.dofs.iter().any(|&v| v >= space.ndof)
        {
            return Err("invalid skeleton Q2 element");
        }
        let idx = |p: [usize; 3]| 9 * p[0] + 3 * p[1] + p[2];
        for corner in 0..8 {
            let side = [corner / 4, corner / 2 % 2, corner % 2];
            let n = idx(side.map(|v| 2 * v));
            let [point] = normalized([block.controls[n]])?;
            let entry = vertices
                .entry(block.dofs[n])
                .or_insert_with(|| (point, Vec::new()));
            if entry.0 != point {
                return Err("shared skeleton vertex differs");
            }
            entry.1.push((block.id.clone(), side));
        }
        for (axis, other) in [[1, 2], [0, 2], [0, 1]].into_iter().enumerate() {
            for signs in 0..4 {
                let side = [signs / 2, signs % 2];
                let indices: [usize; 3] = std::array::from_fn(|i| {
                    let mut p = [0; 3];
                    p[axis] = i;
                    p[other[0]] = 2 * side[0];
                    p[other[1]] = 2 * side[1];
                    idx(p)
                });
                let mut ids = indices.map(|i| block.dofs[i]);
                let mut controls = normalized(indices.map(|i| block.controls[i]))?;
                if ids[0] == ids[1] || ids[1] == ids[2] || ids[0] == ids[2] {
                    return Err("collapsed skeleton edge");
                }
                let forward = ids[0] < ids[2];
                if !forward {
                    ids.reverse();
                    controls.reverse();
                }
                let entry = edges.entry(ids).or_insert_with(|| (controls, Vec::new()));
                if entry.0 != controls {
                    return Err("shared skeleton edge trace differs");
                }
                entry.1.push((block.id.clone(), axis, side, forward));
            }
        }
    }
    let edges = edges
        .into_iter()
        .map(|(dofs, (controls, mut uses))| {
            uses.sort();
            let uses = uses
                .into_iter()
                .map(|(block, axis, side, forward)| {
                    object([
                        ("block", Json::String(block)),
                        ("axis", number(axis)),
                        ("fixed_sides", numbers(side)),
                        ("forward", Json::Bool(forward)),
                    ])
                })
                .collect();
            object([
                ("scalar_dofs", numbers(dofs)),
                ("vertices", numbers([dofs[0], dofs[2]])),
                (
                    "homogeneous_controls",
                    Json::Array(controls.into_iter().map(numbers).collect()),
                ),
                ("uses", Json::Array(uses)),
            ])
        })
        .collect();
    let vertices = vertices
        .into_iter()
        .map(|(dof, (control, mut uses))| {
            uses.sort();
            let uses = uses
                .into_iter()
                .map(|(block, side)| {
                    object([
                        ("block", Json::String(block)),
                        ("local_corner", numbers(side)),
                    ])
                })
                .collect();
            object([
                ("scalar_dof", number(dof)),
                ("homogeneous_control", numbers(control)),
                ("uses", Json::Array(uses)),
            ])
        })
        .collect();
    Ok(object([
        ("edges", Json::Array(edges)),
        ("vertices", Json::Array(vertices)),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{analysis_space::build_analysis_space, bracket::bracket_blocks};
    #[test]
    fn whole_complex_euler_characteristic_and_local_incidence_are_preserved() {
        let mut space = build_analysis_space(&bracket_blocks().unwrap()).unwrap();
        let skeleton = encode_skeleton(&space).unwrap();
        let edges = skeleton.field("edges").unwrap().items().unwrap();
        let vertices = skeleton.field("vertices").unwrap().items().unwrap();
        assert_eq!(vertices.len() as i32 - edges.len() as i32 + 161 - 36, -1);
        let use_count = |a: &[Json]| {
            a.iter()
                .map(|e| e.field("uses").unwrap().items().unwrap().len())
                .sum::<usize>()
        };
        assert_eq!(use_count(edges), 12 * 36);
        assert_eq!(use_count(vertices), 8 * 36);
        space.blocks.reverse();
        assert_eq!(skeleton, encode_skeleton(&space).unwrap());
    }
}
