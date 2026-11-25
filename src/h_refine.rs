//! Exact uniform subdivision and topological node numbering for Q2 blocks.

use crate::analysis_space::{AnalysisBlock, AnalysisSpace};
use std::collections::{BTreeMap, BTreeSet};

type ControlNet = Vec<[i128; 4]>;
type NodeKey = Vec<(usize, u8)>;
type FaceKey = (String, usize, bool);

#[derive(Clone)]
struct ChildNet {
    parent_id: String,
    bits: [usize; 3],
    controls: Vec<[i128; 4]>,
}

fn index(i: usize, j: usize, k: usize) -> usize {
    (i * 3 + j) * 3 + k
}

fn add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b)
        .ok_or("subdivision coefficient exceeds i128")
}

fn mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b)
        .ok_or("subdivision coefficient exceeds i128")
}

fn split_axis(
    controls: &[[i128; 4]],
    axis: usize,
) -> Result<(ControlNet, ControlNet), &'static str> {
    if controls.len() != 27 || axis >= 3 {
        return Err("invalid Q2 control grid for subdivision");
    }
    let mut lower = vec![[0i128; 4]; 27];
    let mut upper = vec![[0i128; 4]; 27];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let mut point = [i, j, k];
                if point[axis] != 0 {
                    continue;
                }
                let mut middle = point;
                let mut end = point;
                middle[axis] = 1;
                end[axis] = 2;
                let a = controls[index(point[0], point[1], point[2])];
                let b = controls[index(middle[0], middle[1], middle[2])];
                let c = controls[index(end[0], end[1], end[2])];
                let mut left_mid = [0i128; 4];
                let mut center = [0i128; 4];
                let mut right_mid = [0i128; 4];
                for component in 0..4 {
                    left_mid[component] = mul(add(a[component], b[component])?, 2)?;
                    center[component] =
                        add(add(a[component], mul(b[component], 2)?)?, c[component])?;
                    right_mid[component] = mul(add(b[component], c[component])?, 2)?;
                    let lower_values = [
                        mul(a[component], 4)?,
                        left_mid[component],
                        center[component],
                    ];
                    let upper_values = [
                        center[component],
                        right_mid[component],
                        mul(c[component], 4)?,
                    ];
                    for (coordinate, value) in lower_values.into_iter().enumerate() {
                        point[axis] = coordinate;
                        lower[index(point[0], point[1], point[2])][component] = value;
                    }
                    for (coordinate, value) in upper_values.into_iter().enumerate() {
                        point[axis] = coordinate;
                        upper[index(point[0], point[1], point[2])][component] = value;
                    }
                }
            }
        }
    }
    Ok((lower, upper))
}

fn children(block: &AnalysisBlock, split_w: bool) -> Result<Vec<ChildNet>, &'static str> {
    if block.controls.len() != 27
        || block.dofs.len() != 27
        || block.controls.iter().any(|p| p[3] <= 0)
    {
        return Err("invalid analysis block for Q2 subdivision");
    }
    let mut patches = vec![([0usize; 3], block.controls.clone())];
    for axis in 0..if split_w { 3 } else { 2 } {
        let mut next = Vec::with_capacity(patches.len() * 2);
        for (bits, controls) in patches {
            let (lower, upper) = split_axis(&controls, axis)?;
            for (bit, child_controls) in [(0usize, lower), (1usize, upper)] {
                let mut child_bits = bits;
                child_bits[axis] = bit;
                next.push((child_bits, child_controls));
            }
        }
        patches = next;
    }
    Ok(patches
        .into_iter()
        .map(|(bits, controls)| ChildNet {
            parent_id: block.id.clone(),
            bits,
            controls,
        })
        .collect())
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn corner_dofs(block: &AnalysisBlock) -> [usize; 8] {
    std::array::from_fn(|corner| {
        let i = 2 * ((corner >> 2) & 1);
        let j = 2 * ((corner >> 1) & 1);
        let k = 2 * (corner & 1);
        block.dofs[index(i, j, k)]
    })
}

fn topology_key(corners: &[usize; 8], grid: [usize; 3]) -> NodeKey {
    let mut weights = BTreeMap::<usize, u64>::new();
    for (corner, &dof) in corners.iter().enumerate() {
        let mut weight = 1u64;
        for (axis, &coordinate) in grid.iter().enumerate() {
            let coordinate = coordinate as u64;
            weight *= if (corner >> (2 - axis)) & 1 == 0 {
                4 - coordinate
            } else {
                coordinate
            };
        }
        if weight != 0 {
            *weights.entry(dof).or_default() += weight;
        }
    }
    let divisor = weights.values().copied().fold(0, gcd);
    debug_assert!(divisor > 0);
    weights
        .into_iter()
        .map(|(dof, weight)| (dof, (weight / divisor) as u8))
        .collect()
}

fn reduced_ratio(numerator: i128, denominator: i128) -> (i128, i128) {
    let mut a = numerator.unsigned_abs();
    let mut b = denominator as u128;
    while b != 0 {
        (a, b) = (b, a % b);
    }
    let divisor = a as i128;
    (numerator / divisor, denominator / divisor)
}

fn same_position(first: [i128; 4], second: [i128; 4]) -> Result<bool, &'static str> {
    if first[3] <= 0 || second[3] <= 0 {
        return Err("subdivided control has nonpositive weight");
    }
    Ok((0..3)
        .all(|axis| reduced_ratio(first[axis], first[3]) == reduced_ratio(second[axis], second[3])))
}

fn local_face_node(axis: usize, upper: bool, first: usize, second: usize) -> usize {
    let mut point = [0usize; 3];
    let free = match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    };
    point[axis] = if upper { 2 } else { 0 };
    point[free[0]] = first;
    point[free[1]] = second;
    index(point[0], point[1], point[2])
}

/// Subdivide each triquadratic block into eight exact children at half spans.
///
/// Node identity is derived from the trilinear rational weights on the eight
/// already identified parent corner DOFs, never from spatial proximity.
/// Homogeneous controls use exact integer midpoint subdivision with one factor
/// of four per split axis. The input is immutable; overflow or a mismatching
/// physical position rejects the complete operation without partial output.
pub fn refine_uniform(parent: &AnalysisSpace) -> Result<AnalysisSpace, &'static str> {
    refine(parent, true)
}

/// Split the two section directions, retaining the full thickness of every
/// block. All columns must be vertical with constant weights and ordered z.
pub fn refine_cross_section(parent: &AnalysisSpace) -> Result<AnalysisSpace, &'static str> {
    for b in &parent.blocks {
        if b.controls.len() != 27 || !b.controls.as_chunks::<3>().0.iter().all(|column| {
            column.windows(2).all(|p| p[0][..2] == p[1][..2] &&
                p[0][3] == p[1][3] && p[0][2] < p[1][2])
        }) {
            return Err("cross-section refinement requires vertical extrusion columns");
        }
    }
    refine(parent, false)
}

fn refine(parent: &AnalysisSpace, split_w: bool) -> Result<AnalysisSpace, &'static str> {
    if parent.blocks.is_empty() || parent.ndof == 0 {
        return Err("empty analysis space");
    }
    if parent
        .boundary_dofs
        .iter()
        .any(|(feature, dofs)| feature.is_empty() || dofs.iter().any(|dof| *dof >= parent.ndof))
    {
        return Err("invalid parent boundary selection");
    }
    let mut blocks: Vec<&AnalysisBlock> = parent.blocks.iter().collect();
    blocks.sort_by(|left, right| left.id.cmp(&right.id));
    let mut ids = BTreeSet::new();
    for block in &blocks {
        if block.id.is_empty()
            || !ids.insert(block.id.as_str())
            || block.controls.len() != 27
            || block.dofs.len() != 27
            || block.dofs.iter().any(|dof| *dof >= parent.ndof)
            || block.controls.iter().any(|control| control[3] <= 0)
        {
            return Err("invalid analysis block");
        }
    }
    let used_dofs: BTreeSet<usize> = blocks
        .iter()
        .flat_map(|block| block.dofs.iter().copied())
        .collect();
    if used_dofs.len() != parent.ndof {
        return Err("parent global degree-of-freedom count is inconsistent");
    }

    let mut parent_corners = BTreeMap::<String, [usize; 8]>::new();
    let mut face_features = BTreeMap::<FaceKey, String>::new();
    for block in &blocks {
        let corners = corner_dofs(block);
        parent_corners.insert(block.id.clone(), corners);
        for axis in 0..3 {
            for upper in [false, true] {
                let center_dof = block.dofs[local_face_node(axis, upper, 1, 1)];
                let matching: Vec<&String> = parent
                    .boundary_dofs
                    .iter()
                    .filter_map(|(feature, dofs)| dofs.contains(&center_dof).then_some(feature))
                    .collect();
                if matching.len() > 1 {
                    return Err("parent design face center has ambiguous feature ownership");
                }
                if let Some(feature) = matching.first() {
                    face_features.insert((block.id.clone(), axis, upper), (*feature).clone());
                }
            }
        }
    }

    let mut child_nets = Vec::<ChildNet>::with_capacity(blocks.len() * 8);
    for block in &blocks {
        child_nets.extend(children(block, split_w)?);
    }

    let mut keys = BTreeSet::<NodeKey>::new();
    let mut child_keys = Vec::with_capacity(child_nets.len());
    for child in &child_nets {
        let corners = parent_corners
            .get(&child.parent_id)
            .ok_or("missing parent corner DOFs")?;
        let grid_keys: Vec<NodeKey> = (0..3)
            .flat_map(|i| (0..3).flat_map(move |j| (0..3).map(move |k| [i, j, k])))
            .map(|local| {
                let grid = [
                    2 * child.bits[0] + local[0],
                    2 * child.bits[1] + local[1],
                    if split_w { 2 * child.bits[2] + local[2] } else { 2 * local[2] },
                ];
                let key = topology_key(corners, grid);
                keys.insert(key.clone());
                key
            })
            .collect();
        child_keys.push(grid_keys);
    }
    let global_ids: BTreeMap<NodeKey, usize> = keys
        .into_iter()
        .enumerate()
        .map(|(id, key)| (key, id))
        .collect();

    let mut position_by_key = BTreeMap::<NodeKey, [i128; 4]>::new();
    let mut refined_blocks = Vec::with_capacity(child_nets.len());
    let mut child_metadata = Vec::with_capacity(child_nets.len());
    for (child, keys) in child_nets.into_iter().zip(child_keys) {
        for (key, control) in keys.iter().cloned().zip(child.controls.iter().copied()) {
            if let Some(previous) = position_by_key.get(&key) {
                if !same_position(*previous, control)? {
                    return Err("topologically identified refined controls disagree");
                }
            } else {
                position_by_key.insert(key, control);
            }
        }
        let dofs = keys.iter().map(|key| global_ids[key]).collect::<Vec<_>>();
        let id = format!(
            "{}/{}{}{}",
            child.parent_id, child.bits[0], child.bits[1], child.bits[2]
        );
        child_metadata.push((child.parent_id, child.bits));
        refined_blocks.push(AnalysisBlock {
            id,
            controls: child.controls,
            dofs,
        });
    }

    let mut boundary_sets = BTreeMap::<String, BTreeSet<usize>>::new();
    for (child, (parent_id, bits)) in refined_blocks.iter().zip(child_metadata) {
        for (axis, &bit) in bits.iter().enumerate() {
            for upper in [false, true] {
                if (axis != 2 || split_w) && bit != usize::from(upper) {
                    continue;
                }
                let Some(feature) = face_features.get(&(parent_id.clone(), axis, upper)) else {
                    continue;
                };
                for a in 0..=2 {
                    for b in 0..=2 {
                        let local = local_face_node(axis, upper, a, b);
                        boundary_sets
                            .entry(feature.clone())
                            .or_default()
                            .insert(child.dofs[local]);
                    }
                }
            }
        }
    }
    if boundary_sets.len() != parent.boundary_dofs.len()
        || parent
            .boundary_dofs
            .keys()
            .any(|feature| !boundary_sets.contains_key(feature))
    {
        return Err("design-face subdivision did not preserve feature coverage");
    }

    Ok(AnalysisSpace {
        ndof: global_ids.len(),
        blocks: refined_blocks,
        boundary_dofs: boundary_sets
            .into_iter()
            .map(|(feature, dofs)| (feature, dofs.into_iter().collect()))
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn original_space() -> AnalysisSpace {
        let blocks = crate::bracket::bracket_blocks().unwrap();
        crate::analysis_space::build_analysis_space(&blocks).unwrap()
    }

    #[test]
    fn first_and_second_refinement_have_expected_topological_counts() {
        let coarse = original_space();
        let first = refine_uniform(&coarse).unwrap();
        assert_eq!(first.blocks.len(), 36 * 8);
        assert_eq!(first.ndof, 3_215);
        let second = refine_uniform(&first).unwrap();
        assert_eq!(second.blocks.len(), 36 * 64);
        assert_eq!(second.ndof, 21_951);
    }

    #[test]
    fn input_order_does_not_change_refined_ids_or_features() {
        let coarse = original_space();
        let forward = refine_uniform(&coarse).unwrap();
        let mut reversed = coarse.clone();
        reversed.blocks.reverse();
        let backward = refine_uniform(&reversed).unwrap();
        assert_eq!(forward, backward);
    }

    #[test]
    fn exact_affine_fields_are_preserved_by_each_child_net() {
        let coarse = original_space();
        let fine = refine_uniform(&coarse).unwrap();
        let parents: BTreeMap<&str, &AnalysisBlock> = coarse
            .blocks
            .iter()
            .map(|block| (block.id.as_str(), block))
            .collect();
        for block in &fine.blocks {
            assert_eq!(block.controls.len(), 27);
            assert_eq!(block.dofs.len(), 27);
            for point in &block.controls {
                assert!(point[3] > 0);
            }
            let (parent_id, child_bits) = block.id.rsplit_once('/').unwrap();
            let parent = parents[parent_id];
            let bits = [
                usize::from(child_bits.as_bytes()[0] == b'1'),
                usize::from(child_bits.as_bytes()[1] == b'1'),
                usize::from(child_bits.as_bytes()[2] == b'1'),
            ];
            let child_basis = [1i128, 2, 1];
            let mut child_h = [0i128; 4];
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        let factor = child_basis[i] * child_basis[j] * child_basis[k];
                        for (acc, component) in
                            child_h.iter_mut().zip(block.controls[index(i, j, k)])
                        {
                            *acc += factor * component;
                        }
                    }
                }
            }
            let quarter_basis = [[9i128, 6, 1], [1i128, 6, 9]];
            let mut parent_h = [0i128; 4];
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        let factor = quarter_basis[bits[0]][i]
                            * quarter_basis[bits[1]][j]
                            * quarter_basis[bits[2]][k];
                        for (acc, component) in
                            parent_h.iter_mut().zip(parent.controls[index(i, j, k)])
                        {
                            *acc += factor * component;
                        }
                    }
                }
            }
            assert_eq!(child_h, parent_h);

            let mut weight_sum = 0i128;
            let mut x_sum = 0i128;
            let mut y_sum = 0i128;
            let mut z_sum = 0i128;
            let mut affine_sum = 0i128;
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        let factor = child_basis[i] * child_basis[j] * child_basis[k];
                        let [x, y, z, weight] = block.controls[index(i, j, k)];
                        weight_sum += factor * weight;
                        x_sum += factor * x;
                        y_sum += factor * y;
                        z_sum += factor * z;
                        affine_sum += factor * (11 * weight + 2 * x - 3 * y + 5 * z);
                    }
                }
            }
            assert!(weight_sum > 0);
            assert_eq!(
                affine_sum,
                11 * weight_sum + 2 * x_sum - 3 * y_sum + 5 * z_sum
            );
        }
        let fine_ids: BTreeSet<usize> = fine
            .blocks
            .iter()
            .flat_map(|block| block.dofs.iter().copied())
            .collect();
        assert_eq!(fine_ids.len(), fine.ndof);
    }

    #[test]
    fn overflow_rejects_without_mutating_the_parent() {
        let mut coarse = original_space();
        coarse.blocks[0].controls[0][0] = i128::MAX;
        let before = coarse.clone();
        assert!(refine_uniform(&coarse).is_err());
        assert_eq!(coarse, before);
    }
}
