//! Exact prolongation between nested dyadic adaptive Q2 spaces.
//!
//! Each output row maps unconstrained source lift masters to one target lift
//! master. This transfer does not include essential boundary conditions; a
//! mechanical solver can transfer a complete xyz field and re-extract its
//! own boundary-constrained reduced vector.

use crate::{
    adaptive_partition::{AdaptivePartition, LeafKey},
    adaptive_space::AdaptiveSpace,
    analysis_space::AnalysisBlock,
    constraint_system::{Lift, Ratio},
    incidence::BlockSide,
};
use std::collections::{BTreeMap, BTreeSet};

type SparseRow = BTreeMap<usize, Ratio>;
type RestrictionMatrix = Vec<Vec<i128>>;
type RootInterfaceSignature = BTreeSet<(BlockSide, BlockSide, (bool, bool, bool))>;
type TopologySignature = (RootInterfaceSignature, BTreeSet<BlockSide>);

#[derive(Debug, Clone, PartialEq, Eq)]
struct SpaceSignature {
    ndof: usize,
    blocks: Vec<(String, Vec<[i128; 4]>, Vec<usize>)>,
    leaf_blocks: BTreeMap<LeafKey, String>,
    lift: Lift,
    boundary_patches: BTreeMap<String, Vec<crate::adaptive_partition::LeafSide>>,
}

/// Exact interpolation from source lift masters to target lift masters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdaptiveTransfer {
    /// Target-master rows, each mapping source-master indices to exact ratios.
    pub rows: Vec<Vec<(usize, Ratio)>>,
    source_lift: Lift,
    target_lift: Lift,
    source_signature: SpaceSignature,
    target_signature: SpaceSignature,
}

impl AdaptiveTransfer {
    pub fn source_master_count(&self) -> usize {
        self.source_lift.free_roots.len()
    }

    pub fn target_master_count(&self) -> usize {
        self.rows.len()
    }

    /// Apply P to component-major reduced values for the unconstrained lift.
    /// This vector does not include essential-boundary elimination.
    pub fn apply_reduced(&self, source: &[f64]) -> Result<Vec<f64>, &'static str> {
        let ns = self.source_master_count();
        if source.len() != 3 * ns || source.iter().any(|x| !x.is_finite()) {
            return Err("invalid source reduced field for adaptive transfer");
        }
        let nt = self.rows.len();
        let mut target = vec![0.; 3 * nt];
        for component in 0..3 {
            for (row_index, row) in self.rows.iter().enumerate() {
                let mut value = 0.;
                for &(column, ratio) in row {
                    let coefficient = ratio.numerator() as f64 / ratio.denominator() as f64;
                    value += coefficient * source[component * ns + column];
                }
                if !value.is_finite() {
                    return Err("nonfinite transferred reduced field");
                }
                target[component * nt + row_index] = value;
            }
        }
        Ok(target)
    }

    /// Transfer a full xyz-interleaved source field and reconstruct the full
    /// target field through its unconstrained exact lift.
    pub fn transfer_full_field(&self, source: &[f64]) -> Result<Vec<f64>, &'static str> {
        let source_ndof = self.source_lift.lines.len();
        if source.len() != 3 * source_ndof || source.iter().any(|x| !x.is_finite()) {
            return Err("invalid full source field for adaptive transfer");
        }
        let source_masters = self.source_lift.free_roots.len();
        let mut reduced = vec![0.; 3 * source_masters];
        for component in 0..3 {
            for (index, &dof) in self.source_lift.free_roots.iter().enumerate() {
                reduced[component * source_masters + index] = source[3 * dof + component];
            }
            for line in &self.source_lift.lines {
                let expected = evaluate_line(line, &reduced, component, source_masters)?;
                let actual = source[3 * line.dof + component];
                if (expected - actual).abs() > 1e-10 * (1. + expected.abs().max(actual.abs())) {
                    return Err("source full field violates its adaptive lift");
                }
            }
        }
        let target_reduced = self.apply_reduced(&reduced)?;
        let nt = self.rows.len();
        let mut full = vec![0.; 3 * self.target_lift.lines.len()];
        for component in 0..3 {
            for line in &self.target_lift.lines {
                full[3 * line.dof + component] =
                    evaluate_line(line, &target_reduced, component, nt)?;
            }
        }
        if full.iter().any(|x| !x.is_finite()) {
            return Err("nonfinite transferred full field");
        }
        Ok(full)
    }

    /// Compose two exact transfers whose intermediate adaptive spaces match.
    pub fn compose(&self, next: &Self) -> Result<Self, &'static str> {
        if self.target_signature != next.source_signature {
            return Err("adaptive transfer composition spaces do not match");
        }
        let mut rows = Vec::with_capacity(next.rows.len());
        for next_row in &next.rows {
            let mut composed = SparseRow::new();
            for &(middle, next_coefficient) in next_row {
                let source_row = self
                    .rows
                    .get(middle)
                    .ok_or("invalid intermediate transfer index")?;
                for &(source, source_coefficient) in source_row {
                    add_to_row(
                        &mut composed,
                        source,
                        next_coefficient.checked_mul(source_coefficient)?,
                    )?;
                }
            }
            validate_unit_row(&composed)?;
            rows.push(composed.into_iter().collect());
        }
        Ok(Self {
            rows,
            source_lift: self.source_lift.clone(),
            target_lift: next.target_lift.clone(),
            source_signature: self.source_signature.clone(),
            target_signature: next.target_signature.clone(),
        })
    }
}

fn evaluate_line(
    line: &crate::constraint_system::LiftLine,
    reduced: &[f64],
    component: usize,
    count: usize,
) -> Result<f64, &'static str> {
    if line.offset != Ratio::zero() {
        return Err("adaptive transfer does not support affine lift offsets");
    }
    let mut value = 0.;
    for &(index, coefficient) in &line.terms {
        if index >= count {
            return Err("invalid adaptive lift master index");
        }
        value += coefficient.numerator() as f64 / coefficient.denominator() as f64
            * reduced[component * count + index];
    }
    if !value.is_finite() {
        return Err("nonfinite adaptive lift evaluation");
    }
    Ok(value)
}

fn signature(space: &AdaptiveSpace) -> SpaceSignature {
    let mut blocks = space
        .space
        .blocks
        .iter()
        .map(|b| (b.id.clone(), b.controls.clone(), b.dofs.clone()))
        .collect::<Vec<_>>();
    blocks.sort_by(|a, b| a.0.cmp(&b.0));
    SpaceSignature {
        ndof: space.space.ndof,
        blocks,
        leaf_blocks: space.leaf_blocks.clone(),
        lift: space.lift.clone(),
        boundary_patches: space.boundary_patches.clone(),
    }
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn canonical_controls(block: &AnalysisBlock) -> Result<Vec<[i128; 4]>, &'static str> {
    if block.controls.len() != 27 || block.controls.iter().any(|p| p[3] <= 0) {
        return Err("invalid adaptive homogeneous control net");
    }
    let divisor = block
        .controls
        .iter()
        .flatten()
        .fold(0_u128, |g, &x| gcd_u128(g, x.unsigned_abs()));
    if divisor == 0 || divisor > i128::MAX as u128 {
        return Err("invalid homogeneous net scale");
    }
    let divisor = divisor as i128;
    Ok(block
        .controls
        .iter()
        .map(|p| p.map(|x| x / divisor))
        .collect())
}

fn split_matrix(bit: u8) -> [[i128; 3]; 3] {
    if bit == 0 {
        [[4, 0, 0], [2, 2, 0], [1, 2, 1]]
    } else {
        [[1, 2, 1], [0, 2, 2], [0, 0, 4]]
    }
}

fn flat([i, j, k]: [usize; 3]) -> usize {
    (i * 3 + j) * 3 + k
}

fn restriction_matrix(path: &[[u8; 3]]) -> Result<RestrictionMatrix, &'static str> {
    let mut matrix = vec![vec![0_i128; 27]; 27];
    for (i, row) in matrix.iter_mut().enumerate() {
        row[i] = 1;
    }
    for octant in path {
        if octant.iter().any(|&b| b > 1) {
            return Err("invalid dyadic transfer path");
        }
        for (axis, &bit) in octant.iter().enumerate() {
            let split = split_matrix(bit);
            let mut next = vec![vec![0_i128; 27]; 27];
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        let output = [i, j, k];
                        let row = &mut next[flat(output)];
                        for (source_axis, &factor) in split[output[axis]].iter().enumerate() {
                            let mut input = output;
                            input[axis] = source_axis;
                            if factor == 0 {
                                continue;
                            }
                            for (dst, &value) in row.iter_mut().zip(&matrix[flat(input)]) {
                                let product = factor
                                    .checked_mul(value)
                                    .ok_or("transfer restriction matrix exceeds i128")?;
                                *dst = dst
                                    .checked_add(product)
                                    .ok_or("transfer restriction matrix exceeds i128")?;
                            }
                        }
                    }
                }
            }
            matrix = next;
        }
        let divisor = matrix
            .iter()
            .flatten()
            .fold(0_u128, |g, &x| gcd_u128(g, x.unsigned_abs()));
        if divisor == 0 || divisor > i128::MAX as u128 {
            return Err("invalid transfer restriction normalization");
        }
        let divisor = divisor as i128;
        for value in matrix.iter_mut().flatten() {
            *value /= divisor;
        }
    }
    Ok(matrix)
}

fn descendant_path(ancestor: &LeafKey, descendant: &LeafKey) -> Result<Vec<[u8; 3]>, &'static str> {
    if ancestor.root != descendant.root || ancestor.level > descendant.level {
        return Err("adaptive transfer attempted a reduction or root change");
    }
    let delta = descendant.level - ancestor.level;
    for axis in 0..3 {
        if descendant.index[axis] >> delta != ancestor.index[axis] {
            return Err("target leaf is not nested in its source leaf");
        }
    }
    Ok(((ancestor.level + 1)..=descendant.level)
        .map(|level| {
            let shift = descendant.level - level;
            std::array::from_fn(|a| ((descendant.index[a] >> shift) & 1) as u8)
        })
        .collect())
}

fn ancestor_for<'a>(
    leaf: &LeafKey,
    source_leaves: &'a BTreeSet<LeafKey>,
) -> Result<&'a LeafKey, &'static str> {
    source_leaves
        .iter()
        .filter(|s| {
            s.root == leaf.root
                && s.level <= leaf.level
                && (0..3).all(|a| leaf.index[a] >> (leaf.level - s.level) == s.index[a])
        })
        .max_by_key(|s| s.level)
        .ok_or("target leaf has no source ancestor")
}

fn add_to_row(row: &mut SparseRow, index: usize, value: Ratio) -> Result<(), &'static str> {
    let sum = row
        .get(&index)
        .copied()
        .unwrap_or_else(Ratio::zero)
        .checked_add(value)?;
    if sum == Ratio::zero() {
        row.remove(&index);
    } else {
        row.insert(index, sum);
    }
    Ok(())
}

fn validate_unit_row(row: &SparseRow) -> Result<(), &'static str> {
    if row.is_empty() {
        return Err("empty adaptive transfer row");
    }
    let sum = row
        .values()
        .try_fold(Ratio::zero(), |a, &b| a.checked_add(b))?;
    if sum != Ratio::one() {
        return Err("adaptive transfer row does not preserve constants");
    }
    Ok(())
}

fn insert_scaled(
    target: &mut SparseRow,
    source: &[(usize, Ratio)],
    scale: Ratio,
) -> Result<(), &'static str> {
    for &(index, coefficient) in source {
        add_to_row(target, index, coefficient.checked_mul(scale)?)?;
    }
    Ok(())
}

fn topology_signature(partition: &AdaptivePartition) -> Result<TopologySignature, &'static str> {
    let report = partition.snapshot()?;
    let mut interfaces = BTreeSet::new();
    for p in report.interfaces {
        if p.first.leaf.root != p.second.leaf.root {
            interfaces.insert((
                p.first.root_side,
                p.second.root_side,
                (
                    p.orientation.swap,
                    p.orientation.flip_first,
                    p.orientation.flip_second,
                ),
            ));
        }
    }
    let exterior = report
        .exterior
        .into_iter()
        .map(|p| p.side.root_side)
        .collect();
    Ok((interfaces, exterior))
}

fn feature_signature(space: &AdaptiveSpace) -> BTreeMap<String, BTreeSet<BlockSide>> {
    space
        .boundary_patches
        .iter()
        .map(|(feature, sides)| {
            (
                feature.clone(),
                sides.iter().map(|s| s.root_side.clone()).collect(),
            )
        })
        .collect()
}

fn validate_lift(space: &AdaptiveSpace) -> Result<(), &'static str> {
    if space.lift.lines.len() != space.space.ndof
        || space.lift.free_roots.iter().any(|&i| i >= space.space.ndof)
        || space.lift.free_roots.windows(2).any(|w| w[0] >= w[1])
    {
        return Err("invalid or affine adaptive lift for transfer");
    }
    for (i, line) in space.lift.lines.iter().enumerate() {
        if line.dof != i || line.offset != Ratio::zero() {
            return Err("invalid or affine adaptive lift for transfer");
        }
        let mut previous = None;
        let mut sum = Ratio::zero();
        for &(index, coefficient) in &line.terms {
            if index >= space.lift.free_roots.len() || previous.is_some_and(|p| p >= index) {
                return Err("invalid adaptive lift terms");
            }
            previous = Some(index);
            sum = sum.checked_add(coefficient)?;
        }
        if sum != Ratio::one() {
            return Err("adaptive lift does not preserve constants");
        }
    }
    for (index, &dof) in space.lift.free_roots.iter().enumerate() {
        if space.lift.lines[dof].terms != [(index, Ratio::one())] {
            return Err("adaptive lift free root is not an identity row");
        }
    }
    Ok(())
}

fn validate_space_mesh(
    space: &AdaptiveSpace,
    partition: &AdaptivePartition,
) -> Result<(), &'static str> {
    let leaves = partition
        .snapshot()?
        .leaves
        .into_iter()
        .collect::<BTreeSet<_>>();
    if space.leaf_blocks.len() != leaves.len()
        || space.leaf_blocks.keys().cloned().collect::<BTreeSet<_>>() != leaves
        || space.space.blocks.len() != leaves.len()
    {
        return Err("adaptive transfer space does not match its partition");
    }
    let blocks = space
        .space
        .blocks
        .iter()
        .map(|b| (b.id.as_str(), b))
        .collect::<BTreeMap<_, _>>();
    if blocks.len() != space.space.blocks.len()
        || space.space.blocks.iter().any(|b| {
            b.controls.len() != 27
                || b.dofs.len() != 27
                || b.controls.iter().any(|h| h[3] <= 0)
                || b.dofs.iter().any(|&d| d >= space.space.ndof)
        })
        || space
            .leaf_blocks
            .values()
            .any(|id| !blocks.contains_key(id.as_str()))
        || space.leaf_blocks.values().collect::<BTreeSet<_>>().len() != leaves.len()
    {
        return Err("adaptive transfer has incomplete or malformed blocks");
    }
    Ok(())
}

/// Build the exact prolongation from a source adaptive mesh to its nested
/// target refinement. Geometry and all shared/hanging occurrences are checked.
pub fn build_transfer(
    source: &AdaptiveSpace,
    source_partition: &AdaptivePartition,
    target: &AdaptiveSpace,
    target_partition: &AdaptivePartition,
) -> Result<AdaptiveTransfer, &'static str> {
    validate_lift(source)?;
    validate_lift(target)?;
    validate_space_mesh(source, source_partition)?;
    validate_space_mesh(target, target_partition)?;
    if source_partition.roots() != target_partition.roots() {
        return Err("adaptive transfer root sets differ");
    }
    if topology_signature(source_partition)? != topology_signature(target_partition)?
        || feature_signature(source) != feature_signature(target)
    {
        return Err("adaptive transfer incidence or boundary selections differ");
    }
    let source_report = source_partition.snapshot()?;
    let target_report = target_partition.snapshot()?;
    if source.leaf_blocks.len() != source_report.leaves.len()
        || target.leaf_blocks.len() != target_report.leaves.len()
    {
        return Err("adaptive transfer spaces are incomplete");
    }
    let source_leaves = source_report
        .leaves
        .iter()
        .cloned()
        .collect::<BTreeSet<_>>();
    let source_blocks = source
        .space
        .blocks
        .iter()
        .map(|b| (b.id.as_str(), b))
        .collect::<BTreeMap<_, _>>();
    let target_blocks = target
        .space
        .blocks
        .iter()
        .map(|b| (b.id.as_str(), b))
        .collect::<BTreeMap<_, _>>();
    let mut full_rows = vec![None::<SparseRow>; target.space.ndof];
    for leaf in &target_report.leaves {
        let ancestor = ancestor_for(leaf, &source_leaves)?;
        let path = descendant_path(ancestor, leaf)?;
        let matrix = restriction_matrix(&path)?;
        let source_id = source
            .leaf_blocks
            .get(ancestor)
            .ok_or("source leaf block is missing")?;
        let target_id = target
            .leaf_blocks
            .get(leaf)
            .ok_or("target leaf block is missing")?;
        let sb = source_blocks
            .get(source_id.as_str())
            .ok_or("source block identity is missing")?;
        let tb = target_blocks
            .get(target_id.as_str())
            .ok_or("target block identity is missing")?;
        let restricted = crate::h_refine::restrict_q2_child(&sb.controls, &path)?;
        if canonical_controls(tb)?
            != canonical_controls(&AnalysisBlock {
                id: sb.id.clone(),
                controls: restricted,
                dofs: sb.dofs.clone(),
            })?
        {
            return Err("target geometry differs from exact dyadic restriction");
        }
        for (output, transform_row) in matrix.iter().enumerate() {
            let mut denominator = 0_i128;
            let mut weighted = Vec::new();
            for (input, &transform_factor) in transform_row.iter().enumerate() {
                let factor = transform_factor
                    .checked_mul(sb.controls[input][3])
                    .ok_or("transfer field weight exceeds i128")?;
                if factor != 0 {
                    denominator = denominator
                        .checked_add(factor)
                        .ok_or("transfer denominator exceeds i128")?;
                    weighted.push((input, factor));
                }
            }
            if denominator <= 0 {
                return Err("invalid restricted transfer weight");
            }
            let mut row = SparseRow::new();
            for (input, factor) in weighted {
                let ratio = Ratio::new(factor, denominator)?;
                let source_dof = sb.dofs[input];
                let line = source
                    .lift
                    .lines
                    .get(source_dof)
                    .ok_or("source lift line is missing")?;
                insert_scaled(&mut row, &line.terms, ratio)?;
            }
            validate_unit_row(&row)?;
            let target_dof = tb.dofs[output];
            if target_dof >= full_rows.len() {
                return Err("target local DOF index is invalid");
            }
            if let Some(previous) = &full_rows[target_dof] {
                if previous != &row {
                    return Err("target DOF occurrences have inconsistent transfer rows");
                }
            } else {
                full_rows[target_dof] = Some(row);
            }
        }
    }
    if full_rows.iter().any(Option::is_none) {
        return Err("target transfer does not cover every full DOF");
    }
    let full_rows = full_rows
        .into_iter()
        .map(Option::unwrap)
        .collect::<Vec<_>>();
    let mut rows = Vec::with_capacity(target.lift.free_roots.len());
    for &dof in &target.lift.free_roots {
        rows.push(
            full_rows
                .get(dof)
                .ok_or("target master DOF is invalid")?
                .iter()
                .map(|(&i, &r)| (i, r))
                .collect::<Vec<_>>(),
        );
    }
    for line in &target.lift.lines {
        let mut expected = SparseRow::new();
        for &(master, coefficient) in &line.terms {
            let row = rows.get(master).ok_or("target lift master is missing")?;
            for &(source, entry) in row {
                add_to_row(&mut expected, source, coefficient.checked_mul(entry)?)?;
            }
        }
        if expected != full_rows[line.dof] {
            return Err("target hanging equation conflicts with transfer");
        }
    }
    Ok(AdaptiveTransfer {
        rows,
        source_lift: source.lift.clone(),
        target_lift: target.lift.clone(),
        source_signature: signature(source),
        target_signature: signature(target),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{adaptive_partition::LeafKey, incidence::Incidence};

    fn cube_parent() -> AdaptiveSpace {
        let controls = (0..27)
            .map(|n| {
                let [i, j, k] = [n / 9, n / 3 % 3, n % 3];
                [
                    (i as i128) * 500,
                    (j as i128) * 500,
                    (k as i128) * 500,
                    1000,
                ]
            })
            .collect();
        let block = AnalysisBlock {
            id: "cube".into(),
            controls,
            dofs: (0..27).collect(),
        };
        let space = crate::analysis_space::AnalysisSpace {
            ndof: 27,
            blocks: vec![block],
            boundary_dofs: BTreeMap::new(),
        };
        AdaptiveSpace {
            space,
            lift: crate::constraint_system::resolve(27, &[], &[]).unwrap(),
            leaf_blocks: BTreeMap::from([(LeafKey::root("cube"), "cube".into())]),
            constraints: Vec::new(),
            boundary_patches: BTreeMap::new(),
            root_ndof: 27,
        }
    }

    fn partition(depth: u8) -> AdaptivePartition {
        let exterior = (0..3)
            .flat_map(|axis| {
                [false, true].map(move |upper| BlockSide {
                    block: "cube".into(),
                    axis,
                    upper,
                })
            })
            .collect();
        let incidence = Incidence {
            interfaces: Vec::new(),
            exterior,
        };
        let mut p = AdaptivePartition::new(&["cube".into()], &incidence, 4, 1000).unwrap();
        for _ in 0..depth {
            let leaves = p.leaves().iter().cloned().collect::<Vec<_>>();
            p.refine(&leaves).unwrap();
        }
        p
    }

    fn build_space(parent: &AdaptiveSpace, p: &AdaptivePartition) -> AdaptiveSpace {
        crate::adaptive_space::build_adaptive_space(&parent.space, p).unwrap()
    }

    fn d4_space() -> (AdaptiveSpace, AdaptivePartition) {
        use crate::incidence::{Interface, Orientation};
        let mut a_controls = Vec::new();
        let mut b_controls = Vec::new();
        let mut a_dofs = vec![usize::MAX; 27];
        let mut b_dofs = vec![usize::MAX; 27];
        let mut next = 0;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let n = (i * 3 + j) * 3 + k;
                    a_controls.push([
                        (i as i128) * 500,
                        (j as i128) * 500,
                        (k as i128) * 500,
                        1000,
                    ]);
                    b_controls.push([
                        1000 + (i as i128) * 500,
                        1000 - (k as i128) * 500,
                        (j as i128) * 500,
                        1000,
                    ]);
                    a_dofs[n] = next;
                    next += 1;
                }
            }
        }
        for j in 0..3 {
            for k in 0..3 {
                let a_node = (2 * 3 + j) * 3 + k;
                let b_node = k * 3 + (2 - j);
                b_dofs[b_node] = a_dofs[a_node];
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let n = (i * 3 + j) * 3 + k;
                    if b_dofs[n] == usize::MAX {
                        b_dofs[n] = next;
                        next += 1;
                    }
                }
            }
        }
        let blocks = vec![
            AnalysisBlock {
                id: "a".into(),
                controls: a_controls,
                dofs: a_dofs,
            },
            AnalysisBlock {
                id: "b".into(),
                controls: b_controls,
                dofs: b_dofs,
            },
        ];
        let side_a = BlockSide {
            block: "a".into(),
            axis: 0,
            upper: true,
        };
        let side_b = BlockSide {
            block: "b".into(),
            axis: 0,
            upper: false,
        };
        let orientation = Orientation {
            swap: true,
            flip_first: true,
            flip_second: false,
        };
        let mut exterior = Vec::new();
        for id in ["a", "b"] {
            for axis in 0..3 {
                for upper in [false, true] {
                    let side = BlockSide {
                        block: id.into(),
                        axis,
                        upper,
                    };
                    if side != side_a && side != side_b {
                        exterior.push(side);
                    }
                }
            }
        }
        let incidence = Incidence {
            interfaces: vec![Interface {
                first: side_a,
                second: side_b,
                orientation,
            }],
            exterior,
        };
        let roots = ["a".into(), "b".into()];
        let p = AdaptivePartition::new(&roots, &incidence, 3, 500).unwrap();
        let space = AdaptiveSpace {
            space: crate::analysis_space::AnalysisSpace {
                ndof: next,
                blocks,
                boundary_dofs: BTreeMap::new(),
            },
            lift: crate::constraint_system::resolve(next, &[], &[]).unwrap(),
            leaf_blocks: BTreeMap::from([
                (LeafKey::root("a"), "a".into()),
                (LeafKey::root("b"), "b".into()),
            ]),
            constraints: Vec::new(),
            boundary_patches: BTreeMap::new(),
            root_ndof: next,
        };
        (space, p)
    }

    #[test]
    fn identity_composition_and_affine_geometry_transfer_are_exact() {
        let parent = cube_parent();
        let p0 = partition(0);
        let s0 = build_space(&parent, &p0);
        let identity = build_transfer(&s0, &p0, &s0, &p0).unwrap();
        assert_eq!(identity.rows.len(), 27);
        for i in 0..27 {
            assert_eq!(identity.rows[i], vec![(i, Ratio::one())]);
        }
        let p1 = partition(1);
        let s1 = build_space(&parent, &p1);
        let p2 = partition(2);
        let s2 = build_space(&parent, &p2);
        let t01 = build_transfer(&s0, &p0, &s1, &p1).unwrap();
        let t12 = build_transfer(&s1, &p1, &s2, &p2).unwrap();
        let t02 = build_transfer(&s0, &p0, &s2, &p2).unwrap();
        assert_eq!(t01.compose(&t12).unwrap().rows, t02.rows);
        let source = (0..3 * 27)
            .map(|i| ((i * 37 % 101) as f64 - 50.) / 13.)
            .collect::<Vec<_>>();
        let direct = t01.apply_reduced(&source).unwrap();
        let transferred_full = t01
            .transfer_full_field(&source_to_full(&s0, &source))
            .unwrap();
        let recovered = full_to_reduced(&s1, &transferred_full);
        for (a, b) in direct.iter().zip(recovered) {
            assert!((a - b).abs() < 1e-12);
        }
        let affine = (0..27)
            .flat_map(|n| {
                let [i, j, k] = [n / 9, n / 3 % 3, n % 3];
                [
                    0.1 + 0.2 * i as f64 / 2.,
                    -0.3 + 0.7 * j as f64 / 2.,
                    0.4 - 0.1 * k as f64 / 2.,
                ]
            })
            .collect::<Vec<_>>();
        let transferred = t02.transfer_full_field(&affine).unwrap();
        for block in &s2.space.blocks {
            for (local, &dof) in block.dofs.iter().enumerate() {
                let h = block.controls[local];
                let xyz = std::array::from_fn::<_, 3, _>(|a| h[a] as f64 / h[3] as f64);
                let expected = [0.1 + 0.2 * xyz[0], -0.3 + 0.7 * xyz[1], 0.4 - 0.1 * xyz[2]];
                for c in 0..3 {
                    assert!((transferred[3 * dof + c] - expected[c]).abs() < 1e-12);
                }
            }
        }
    }

    fn source_to_full(space: &AdaptiveSpace, reduced: &[f64]) -> Vec<f64> {
        let count = space.lift.free_roots.len();
        let mut full = vec![0.; space.space.ndof * 3];
        for c in 0..3 {
            for line in &space.lift.lines {
                full[3 * line.dof + c] = line
                    .terms
                    .iter()
                    .map(|&(i, r)| {
                        r.numerator() as f64 / r.denominator() as f64 * reduced[c * count + i]
                    })
                    .sum();
            }
        }
        full
    }

    fn full_to_reduced(space: &AdaptiveSpace, full: &[f64]) -> Vec<f64> {
        let count = space.lift.free_roots.len();
        let mut reduced = vec![0.; 3 * count];
        for c in 0..3 {
            for (i, &dof) in space.lift.free_roots.iter().enumerate() {
                reduced[c * count + i] = full[3 * dof + c];
            }
        }
        reduced
    }

    #[test]
    fn rejects_reduction_and_changed_geometry_or_topology() {
        let parent = cube_parent();
        let p1 = partition(1);
        let s1 = build_space(&parent, &p1);
        let p0 = partition(0);
        let s0 = build_space(&parent, &p0);
        assert!(build_transfer(&s1, &p1, &s0, &p0).is_err());
        let mut changed = s1.clone();
        changed.space.blocks[0].controls[0][0] += 1;
        assert!(build_transfer(&s0, &p0, &changed, &p1).is_err());
    }

    #[test]
    fn arbitrary_field_transfer_preserves_d4_hanging_constraints() {
        let (parent, p0) = d4_space();
        let s0 = build_space(&parent, &p0);
        let mut p1 = p0.clone();
        p1.refine(&[LeafKey::root("a")]).unwrap();
        let s1 = build_space(&parent, &p1);
        let transfer = build_transfer(&s0, &p0, &s1, &p1).unwrap();
        let source = (0..3 * s0.lift.free_roots.len())
            .map(|i| ((i * 43 % 97) as f64 - 31.) / 19.)
            .collect::<Vec<_>>();
        let mapped = transfer.apply_reduced(&source).unwrap();
        let full = transfer
            .transfer_full_field(&source_to_full(&s0, &source))
            .unwrap();
        assert_eq!(mapped.len(), 3 * s1.lift.free_roots.len());
        let extracted = full_to_reduced(&s1, &full);
        for (a, b) in mapped.iter().zip(extracted) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn affine_field_transfers_to_local_level_two_bracket_refinement() {
        let blocks = crate::bracket::bracket_blocks().unwrap();
        let certificate = crate::bracket_audit::audit_bracket(&blocks).unwrap();
        let parent_space = crate::analysis_space::from_certified_faces(
            &blocks,
            &certificate.incidence,
            &certificate.design_faces,
        )
        .unwrap();
        let ids = blocks.iter().map(|b| b.id.clone()).collect::<Vec<_>>();
        let mut p0 = AdaptivePartition::new(&ids, &certificate.incidence, 4, 1500).unwrap();
        let focus = LeafKey::root("hole_x_sector_east_z10_30");
        p0.refine(std::slice::from_ref(&focus)).unwrap();
        let s0 = crate::adaptive_space::build_adaptive_space(&parent_space, &p0).unwrap();
        let mut p1 = p0.clone();
        p1.refine(&[LeafKey {
            root: focus.root,
            level: 1,
            index: [1, 0, 0],
        }])
        .unwrap();
        let s1 = crate::adaptive_space::build_adaptive_space(&parent_space, &p1).unwrap();
        let transfer = build_transfer(&s0, &p0, &s1, &p1).unwrap();
        let mut source = vec![0.; 3 * s0.space.ndof];
        for block in &s0.space.blocks {
            for (local, &dof) in block.dofs.iter().enumerate() {
                let h = block.controls[local];
                let xyz = std::array::from_fn::<_, 3, _>(|a| h[a] as f64 / h[3] as f64);
                let value = [
                    0.001 + 0.02 * xyz[0] - 0.03 * xyz[1],
                    -0.002 + 0.01 * xyz[1] + 0.015 * xyz[2],
                    0.003 - 0.012 * xyz[0] + 0.01 * xyz[2],
                ];
                for c in 0..3 {
                    source[3 * dof + c] = value[c];
                }
            }
        }
        let target = transfer.transfer_full_field(&source).unwrap();
        for block in &s1.space.blocks {
            for (local, &dof) in block.dofs.iter().enumerate() {
                let h = block.controls[local];
                let xyz = std::array::from_fn::<_, 3, _>(|a| h[a] as f64 / h[3] as f64);
                let expected = [
                    0.001 + 0.02 * xyz[0] - 0.03 * xyz[1],
                    -0.002 + 0.01 * xyz[1] + 0.015 * xyz[2],
                    0.003 - 0.012 * xyz[0] + 0.01 * xyz[2],
                ];
                for c in 0..3 {
                    assert!((target[3 * dof + c] - expected[c]).abs() < 2e-12);
                }
            }
        }
    }
}
