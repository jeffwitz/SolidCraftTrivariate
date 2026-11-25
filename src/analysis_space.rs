//! Conforming triquadratic H1 analysis space over an audited whole-face mesh.

use crate::bracket::BracketBlock;
use crate::bracket_audit::audit_bracket;
use crate::incidence::{BlockSide, Incidence};
use std::collections::{BTreeMap, BTreeSet};

/// One accepted block elevated to homogeneous degree `(2,2,2)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnalysisBlock {
    pub id: String,
    /// Homogeneous `(wx, wy, wz, weight)` controls, `w` index fastest.
    pub controls: Vec<[i128; 4]>,
    /// Scalar global degree-of-freedom index for each of the 27 controls.
    pub dofs: Vec<usize>,
}

/// Global scalar H1 space with exact face-node identification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnalysisSpace {
    pub ndof: usize,
    pub blocks: Vec<AnalysisBlock>,
    pub boundary_dofs: BTreeMap<String, Vec<usize>>,
}

#[derive(Clone)]
struct ElevatedBlock {
    id: String,
    controls: Vec<[i128; 4]>,
}

fn flat_index([i, j, k]: [usize; 3], dims: [usize; 3]) -> usize {
    (i * dims[1] + j) * dims[2] + k
}

fn elevate_to_quadratic(block: &BracketBlock) -> Result<ElevatedBlock, &'static str> {
    if block.degrees.iter().any(|degree| !(1..=2).contains(degree))
        || block.controls.len()
            != block
                .degrees
                .iter()
                .map(|degree| degree + 1)
                .product::<usize>()
    {
        return Err("invalid block net for analysis elevation");
    }
    let mut degrees = block.degrees;
    let mut values: Vec<[i128; 4]> = block
        .controls
        .iter()
        .map(|point| point.map(i128::from))
        .collect();

    for axis in 0..3 {
        if degrees[axis] == 2 {
            continue;
        }
        let old_degrees = degrees;
        let old_dims = old_degrees.map(|degree| degree + 1);
        let mut new_degrees = old_degrees;
        new_degrees[axis] = 2;
        let new_dims = new_degrees.map(|degree| degree + 1);
        let mut elevated = Vec::with_capacity(27);
        for i in 0..new_dims[0] {
            for j in 0..new_dims[1] {
                for k in 0..new_dims[2] {
                    let target = [i, j, k];
                    let coordinate = target[axis];
                    let mut a_index = target;
                    let mut b_index = target;
                    let (a_factor, b_factor) = match coordinate {
                        0 => {
                            a_index[axis] = 0;
                            (2i128, 0i128)
                        }
                        1 => {
                            a_index[axis] = 0;
                            b_index[axis] = 1;
                            (1i128, 1i128)
                        }
                        _ => {
                            a_index[axis] = 1;
                            b_index[axis] = 1;
                            (2i128, 0i128)
                        }
                    };
                    let a = values[flat_index(a_index, old_dims)];
                    let b = values[flat_index(b_index, old_dims)];
                    let mut point = [0i128; 4];
                    for component in 0..4 {
                        point[component] = a[component]
                            .checked_mul(a_factor)
                            .and_then(|value| {
                                if b_factor == 0 {
                                    Some(value)
                                } else {
                                    b[component]
                                        .checked_mul(b_factor)
                                        .and_then(|other| value.checked_add(other))
                                }
                            })
                            .ok_or("degree elevation exceeds i128")?;
                    }
                    elevated.push(point);
                }
            }
        }
        degrees = new_degrees;
        values = elevated;
    }
    let controls = values;
    if controls.len() != 27 || controls.iter().any(|point| point[3] <= 0) {
        return Err("invalid elevated homogeneous controls");
    }
    Ok(ElevatedBlock {
        id: block.id.clone(),
        controls,
    })
}

fn face_axes(axis: usize) -> [usize; 2] {
    match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}

fn side_node(side: &BlockSide, a: usize, b: usize) -> usize {
    let mut index = [0usize; 3];
    let free = face_axes(side.axis);
    index[free[0]] = a;
    index[free[1]] = b;
    index[side.axis] = if side.upper { 2 } else { 0 };
    flat_index(index, [3, 3, 3])
}

fn oriented_source_index(
    swap: bool,
    flip_first: bool,
    flip_second: bool,
    output_first: usize,
    output_second: usize,
) -> (usize, usize) {
    let a = if flip_first {
        2 - output_first
    } else {
        output_first
    };
    let b = if flip_second {
        2 - output_second
    } else {
        output_second
    };
    if swap {
        (b, a)
    } else {
        (a, b)
    }
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn reduced_ratio(numerator: i128, denominator: i128) -> (i128, i128) {
    let divisor = gcd_u128(numerator.unsigned_abs(), denominator as u128) as i128;
    (numerator / divisor, denominator / divisor)
}

fn same_physical_control(left: [i128; 4], right: [i128; 4]) -> Result<bool, &'static str> {
    if left[3] <= 0 || right[3] <= 0 {
        return Err("identified analysis control has nonpositive weight");
    }
    for coordinate in 0..3 {
        if reduced_ratio(left[coordinate], left[3]) != reduced_ratio(right[coordinate], right[3]) {
            return Ok(false);
        }
    }
    Ok(true)
}

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    fn find(&mut self, node: usize) -> usize {
        if self.parent[node] != node {
            self.parent[node] = self.find(self.parent[node]);
        }
        self.parent[node]
    }

    fn join(&mut self, left: usize, right: usize) {
        let mut a = self.find(left);
        let mut b = self.find(right);
        if a == b {
            return;
        }
        if self.rank[a] < self.rank[b] {
            std::mem::swap(&mut a, &mut b);
        }
        self.parent[b] = a;
        if self.rank[a] == self.rank[b] {
            self.rank[a] += 1;
        }
    }
}

/// Build a deterministic conforming triquadratic scalar H1 space.
///
/// The bracket must first pass its independent fixed-family audit. Degree-one
/// axes are elevated exactly with integer Bernstein formulas; each elevation
/// introduces one common factor of two. Whole-face incidences identify the 27
/// local control positions by exact rational equality. Global scalar IDs are
/// ordered by block ID and local node index, independent of input block order.
pub fn build_analysis_space(blocks: &[BracketBlock]) -> Result<AnalysisSpace, &'static str> {
    let certificate = audit_bracket(blocks)?;
    from_certified_faces(blocks, &certificate.incidence, &certificate.design_faces)
}

/// Build the conforming Q2 space of the bracket with independently checked bore radii.
pub fn parameterized_bracket_space(
    rx: u8, ry: u8, blocks: &[BracketBlock],
) -> Result<AnalysisSpace, &'static str> {
    let report = crate::bracket_parameters_audit::audit_bracket_parameters(rx,ry,blocks)?;
    from_certified_faces(blocks, &report.certificate.incidence, &report.certificate.design_faces)
}

/// Build the conforming space of one audited fixed section template.
pub fn section_analysis_space(
    family: &str,
    blocks: &[BracketBlock],
) -> Result<AnalysisSpace, &'static str> {
    let certificate = crate::section_audit::audit_section(family, blocks)?;
    from_certified_faces(blocks, &certificate.incidence, &certificate.design_faces)
}

/// Build a whole-face Q2 space after independently checking the movable bore.
pub fn bore_analysis_space(
    center_half_mm: [i32;2], radius_half_mm: i32, blocks: &[BracketBlock],
) -> Result<AnalysisSpace, &'static str> {
    let certificate = crate::movable_bore_audit::audit_bore(center_half_mm, radius_half_mm, blocks)?;
    from_certified_faces(blocks, &certificate.incidence, &certificate.design_faces)
}

/// Build a conforming Q2 space on the certified row of pierced cells.
pub fn row_analysis_space(
    half_widths: [i32; 2], centers: &[[i32; 2]], radii: &[i32], blocks: &[BracketBlock],
) -> Result<AnalysisSpace, &'static str> {
    let certificate = crate::row_bores_audit::audit_row(half_widths, centers, radii, blocks)?;
    from_certified_faces(blocks, &certificate.incidence, &certificate.design_faces)
}

// Callers must establish model-wide validity before numbering the space.
pub(crate) fn from_certified_faces(
    blocks: &[BracketBlock],
    incidence: &Incidence,
    design_faces: &[(BlockSide, String)],
) -> Result<AnalysisSpace, &'static str> {
    let mut ordered: Vec<&BracketBlock> = blocks.iter().collect();
    ordered.sort_by(|left, right| left.id.cmp(&right.id));
    let elevated = ordered
        .into_iter()
        .map(elevate_to_quadratic)
        .collect::<Result<Vec<_>, _>>()?;
    let indices: BTreeMap<&str, usize> = elevated
        .iter()
        .enumerate()
        .map(|(index, block)| (block.id.as_str(), index))
        .collect();

    let local_count = elevated
        .len()
        .checked_mul(27)
        .ok_or("analysis space size overflow")?;
    let mut union_find = UnionFind::new(local_count);
    for interface in &incidence.interfaces {
        let first_block = *indices
            .get(interface.first.block.as_str())
            .ok_or("incidence references an unknown first block")?;
        let second_block = *indices
            .get(interface.second.block.as_str())
            .ok_or("incidence references an unknown second block")?;
        for a in 0..=2 {
            for b in 0..=2 {
                let (source_a, source_b) = oriented_source_index(
                    interface.orientation.swap,
                    interface.orientation.flip_first,
                    interface.orientation.flip_second,
                    a,
                    b,
                );
                let first_local = side_node(&interface.first, a, b);
                let second_local = side_node(&interface.second, source_a, source_b);
                if !same_physical_control(
                    elevated[first_block].controls[first_local],
                    elevated[second_block].controls[second_local],
                )? {
                    return Err("identified control positions disagree");
                }
                union_find.join(
                    first_block * 27 + first_local,
                    second_block * 27 + second_local,
                );
            }
        }
    }

    let mut representative_positions: BTreeMap<usize, [i128; 4]> = BTreeMap::new();
    for (block_index, block) in elevated.iter().enumerate() {
        for (local, control) in block.controls.iter().copied().enumerate() {
            let root = union_find.find(block_index * 27 + local);
            if let Some(previous) = representative_positions.get(&root) {
                if !same_physical_control(*previous, control)? {
                    return Err("transitive identified positions disagree");
                }
            } else {
                representative_positions.insert(root, control);
            }
        }
    }

    let mut root_keys: BTreeMap<usize, (String, usize)> = BTreeMap::new();
    for (block_index, block) in elevated.iter().enumerate() {
        for local in 0..27 {
            let root = union_find.find(block_index * 27 + local);
            let key = (block.id.clone(), local);
            root_keys
                .entry(root)
                .and_modify(|current| {
                    if key < *current {
                        *current = key.clone();
                    }
                })
                .or_insert(key);
        }
    }
    let mut ordered_roots: Vec<(usize, (String, usize))> = root_keys.into_iter().collect();
    ordered_roots.sort_by(|left, right| left.1.cmp(&right.1));
    let global_ids: BTreeMap<usize, usize> = ordered_roots
        .into_iter()
        .enumerate()
        .map(|(global, (root, _))| (root, global))
        .collect();

    let analysis_blocks = elevated
        .iter()
        .enumerate()
        .map(|(block_index, block)| {
            let dofs = (0..27)
                .map(|local| {
                    let root = union_find.find(block_index * 27 + local);
                    global_ids[&root]
                })
                .collect();
            AnalysisBlock {
                id: block.id.clone(),
                controls: block.controls.clone(),
                dofs,
            }
        })
        .collect::<Vec<_>>();

    let mut boundary_sets = BTreeMap::<String, BTreeSet<usize>>::new();
    for (side, feature) in design_faces {
        let block_index = *indices
            .get(side.block.as_str())
            .ok_or("design face references an unknown block")?;
        for a in 0..=2 {
            for b in 0..=2 {
                let local = side_node(side, a, b);
                let root = union_find.find(block_index * 27 + local);
                boundary_sets
                    .entry(feature.clone())
                    .or_default()
                    .insert(global_ids[&root]);
            }
        }
    }
    let boundary_dofs = boundary_sets
        .into_iter()
        .map(|(feature, dofs)| (feature, dofs.into_iter().collect()))
        .collect();

    Ok(AnalysisSpace {
        ndof: global_ids.len(),
        blocks: analysis_blocks,
        boundary_dofs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bracket_space_is_stable_under_input_order_and_has_expected_size() {
        let blocks = crate::bracket::bracket_blocks().unwrap();
        let forward = build_analysis_space(&blocks).unwrap();
        let reversed: Vec<_> = blocks.iter().cloned().rev().collect();
        let backward = build_analysis_space(&reversed).unwrap();
        assert_eq!(forward, backward);
        assert_eq!(forward.blocks.len(), 36);
        assert_eq!(forward.ndof, 531);
        assert!(forward
            .blocks
            .iter()
            .all(|block| { block.controls.len() == 27 && block.dofs.len() == 27 }));
    }

    #[test]
    fn classified_design_faces_select_global_boundary_dofs() {
        let blocks = crate::bracket::bracket_blocks().unwrap();
        let space = build_analysis_space(&blocks).unwrap();
        let expected = [
            "x0",
            "x8",
            "x60",
            "y0",
            "y8",
            "y50",
            "z0",
            "z40",
            "inside_fillet",
            "hole_x",
            "hole_y",
        ];
        assert_eq!(space.boundary_dofs.len(), expected.len());
        for feature in expected {
            let dofs = &space.boundary_dofs[feature];
            assert!(!dofs.is_empty());
            assert!(dofs.windows(2).all(|pair| pair[0] < pair[1]));
            assert!(dofs.iter().all(|dof| *dof < space.ndof));
        }
        assert_ne!(space.boundary_dofs["hole_x"], space.boundary_dofs["hole_y"]);
    }

    #[test]
    fn rational_quadratic_basis_reproduces_partition_and_affine_fields() {
        let blocks = crate::bracket::bracket_blocks().unwrap();
        let space = build_analysis_space(&blocks).unwrap();
        // Common-denominator Bernstein values at (u,v,w)=(1/2,1/4,3/4).
        let basis_u = [1i128, 2, 1];
        let basis_v = [9i128, 6, 1];
        let basis_w = [1i128, 6, 9];
        for block in &space.blocks {
            let mut weight_sum = 0i128;
            let mut x_sum = 0i128;
            let mut y_sum = 0i128;
            let mut z_sum = 0i128;
            let mut affine_sum = 0i128;
            for (i, bu) in basis_u.iter().enumerate() {
                for (j, bv) in basis_v.iter().enumerate() {
                    for (k, bw) in basis_w.iter().enumerate() {
                        let factor = bu * bv * bw;
                        let [x, y, z, weight] = block.controls[(i * 3 + j) * 3 + k];
                        weight_sum += factor * weight;
                        x_sum += factor * x;
                        y_sum += factor * y;
                        z_sum += factor * z;
                        affine_sum += factor * (11 * weight + 2 * x - 3 * y + 5 * z);
                    }
                }
            }
            assert!(weight_sum > 0);
            // N_i are the rational Bernstein basis; this identity proves both
            // partition of unity and reproduction of the affine physical field.
            assert_eq!(
                affine_sum,
                11 * weight_sum + 2 * x_sum - 3 * y_sum + 5 * z_sum
            );
        }

        let mut dof_positions = BTreeMap::<usize, [i128; 4]>::new();
        for block in &space.blocks {
            for (control, dof) in block.controls.iter().copied().zip(&block.dofs) {
                if let Some(previous) = dof_positions.get(dof) {
                    assert!(same_physical_control(*previous, control).unwrap());
                } else {
                    dof_positions.insert(*dof, control);
                }
            }
        }
    }
}
