//! Exact Q2 analysis spaces for locally refined multipatch complexes.
//!
//! Same-level traces identify Bernstein coefficients by certified D4 face
//! incidence. Nonmatching traces contribute exact rational restriction rows;
//! they are resolved together so duplicate edge and corner equations agree.

use crate::{
    adaptive_partition::{AdaptivePartition, ExteriorPatch, FacePatch, LeafKey, LeafSide},
    analysis_space::{AnalysisBlock, AnalysisSpace},
    constraint_system::{self, Lift},
    face_constraints::{nested_constraints_wide, FaceNetWide, Row, Span},
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdaptiveSpace {
    /// Unconstrained Q2 finite-element space with only exact trace identities
    /// merged. Hanging coefficients remain present in this space.
    pub space: AnalysisSpace,
    /// Exact affine representation of every DOF using the independent masters.
    pub lift: Lift,
    /// Stable child-block identity by dyadic leaf.
    pub leaf_blocks: BTreeMap<LeafKey, String>,
    /// Canonical exact hanging-node equations after same-level/identity merges.
    pub constraints: Vec<(usize, Row)>,
    /// Semantic boundary feature to the originating root sides and leaves.
    pub boundary_patches: BTreeMap<String, Vec<LeafSide>>,
    pub root_ndof: usize,
}

#[derive(Debug, Clone)]
struct Child {
    leaf: LeafKey,
    block: AnalysisBlock,
}

#[derive(Debug, Clone)]
struct RawConstraint {
    dependent: usize,
    terms: Vec<(usize, i128)>,
    denominator: i128,
}

type BoundaryDofs = BTreeMap<String, Vec<usize>>;
type BoundaryPatches = BTreeMap<String, Vec<LeafSide>>;

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
    fn find(&mut self, n: usize) -> usize {
        if self.parent[n] != n {
            self.parent[n] = self.find(self.parent[n]);
        }
        self.parent[n]
    }
    fn join(&mut self, a: usize, b: usize) -> bool {
        let (mut a, mut b) = (self.find(a), self.find(b));
        if a == b {
            return false;
        }
        if self.rank[a] < self.rank[b] {
            std::mem::swap(&mut a, &mut b);
        }
        self.parent[b] = a;
        if self.rank[a] == self.rank[b] {
            self.rank[a] += 1;
        }
        true
    }
}

fn flat(i: usize, j: usize, k: usize) -> usize {
    (i * 3 + j) * 3 + k
}

fn face_axes(axis: usize) -> [usize; 2] {
    match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}

fn face_local(side_axis: usize, upper: bool, a: usize, b: usize) -> usize {
    let mut point = [0; 3];
    let tangent = face_axes(side_axis);
    point[side_axis] = if upper { 2 } else { 0 };
    point[tangent[0]] = a;
    point[tangent[1]] = b;
    flat(point[0], point[1], point[2])
}

fn face_net(block: &AnalysisBlock, axis: usize, upper: bool) -> Result<FaceNetWide, &'static str> {
    if block.controls.len() != 27 {
        return Err("invalid Q2 control net on an interface face");
    }
    let mut controls = Vec::with_capacity(9);
    for a in 0..3 {
        for b in 0..3 {
            let control = block.controls[face_local(axis, upper, a, b)];
            if control[3] <= 0 {
                return Err("nonpositive face weight");
            }
            controls.push(control);
        }
    }
    let divisor = controls
        .iter()
        .flatten()
        .copied()
        .fold(0_u128, |a, b| gcd_u128(a, b.unsigned_abs()));
    if divisor == 0 || divisor > i128::MAX as u128 {
        return Err("invalid common face scale");
    }
    let divisor = divisor as i128;
    let controls = controls
        .into_iter()
        .map(|point| point.map(|x| x / divisor))
        .collect();
    Ok(FaceNetWide {
        degrees: [2, 2],
        controls,
    })
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn same_position(a: [i128; 4], b: [i128; 4]) -> Result<bool, &'static str> {
    if a[3] <= 0 || b[3] <= 0 {
        return Err("nonpositive analysis weight");
    }
    for axis in 0..3 {
        if reduced_ratio(a[axis], a[3]) != reduced_ratio(b[axis], b[3]) {
            return Ok(false);
        }
    }
    Ok(true)
}

fn reduced_ratio(n: i128, d: i128) -> (i128, i128) {
    let mut a = n.unsigned_abs();
    let mut b = d as u128;
    while b != 0 {
        (a, b) = (b, a % b);
    }
    let divisor = a.max(1) as i128;
    (n / divisor, d / divisor)
}

fn oriented_source_index(
    swap: bool,
    flip_first: bool,
    flip_second: bool,
    a: usize,
    b: usize,
) -> (usize, usize) {
    let a = if flip_first { 2 - a } else { a };
    let b = if flip_second { 2 - b } else { b };
    if swap {
        (b, a)
    } else {
        (a, b)
    }
}

fn inverse_orientation(o: crate::incidence::Orientation) -> crate::incidence::Orientation {
    if o.swap {
        crate::incidence::Orientation {
            swap: true,
            flip_first: o.flip_second,
            flip_second: o.flip_first,
        }
    } else {
        o
    }
}

fn spans_for_patch(
    master: &LeafKey,
    axis: usize,
    patch: crate::adaptive_partition::FaceRect,
    max_depth: u8,
) -> Result<[Span; 2], &'static str> {
    let scale = 1_u32 << (max_depth - master.level);
    let tangent = face_axes(axis);
    let mut result = [Span::Full; 2];
    for d in 0..2 {
        let lo = master.index[tangent[d]] * scale;
        let hi = (master.index[tangent[d]] + 1) * scale;
        let mid = lo + (hi - lo) / 2;
        result[d] = if patch.ranges[d] == [lo, hi] {
            Span::Full
        } else if patch.ranges[d] == [lo, mid] {
            Span::Lower
        } else if patch.ranges[d] == [mid, hi] {
            Span::Upper
        } else {
            return Err("interface patch is not a dyadic full or half span");
        };
    }
    Ok(result)
}

fn root_block_face_feature(
    root: &AnalysisBlock,
    side_axis: usize,
    upper: bool,
    boundary_dofs: &BTreeMap<String, Vec<usize>>,
) -> Result<Option<String>, &'static str> {
    let center = root.dofs[face_local(side_axis, upper, 1, 1)];
    let matches = boundary_dofs
        .iter()
        .filter_map(|(feature, dofs)| dofs.contains(&center).then_some(feature.clone()))
        .collect::<Vec<_>>();
    if matches.len() > 1 {
        return Err("root face center has ambiguous design-feature ownership");
    }
    Ok(matches.into_iter().next())
}

fn validate_parent(
    parent: &AnalysisSpace,
    partition: &AdaptivePartition,
) -> Result<(), &'static str> {
    if parent.ndof == 0 || parent.blocks.is_empty() {
        return Err("empty parent analysis space");
    }
    let mut ids = BTreeSet::new();
    let mut used = BTreeSet::new();
    for block in &parent.blocks {
        if block.id.is_empty()
            || !ids.insert(block.id.as_str())
            || block.controls.len() != 27
            || block.dofs.len() != 27
            || block.controls.iter().any(|p| p[3] <= 0)
            || block.dofs.iter().any(|&d| d >= parent.ndof)
        {
            return Err("invalid parent analysis block");
        }
        used.extend(block.dofs.iter().copied());
    }
    if used.len() != parent.ndof
        || used.iter().copied().ne(0..parent.ndof)
        || ids.iter().map(|s| (*s).to_owned()).collect::<BTreeSet<_>>()
            != partition.roots().iter().cloned().collect()
    {
        return Err("parent analysis roots or scalar numbering disagree with partition");
    }
    let mut features = BTreeSet::new();
    for (name, dofs) in &parent.boundary_dofs {
        if name.is_empty()
            || dofs.is_empty()
            || dofs.iter().any(|&d| d >= parent.ndof)
            || dofs.iter().copied().collect::<BTreeSet<_>>().len() != dofs.len()
            || !features.insert(name)
        {
            return Err("invalid parent physical boundary selection");
        }
    }
    Ok(())
}

fn child_path(leaf: &LeafKey) -> Result<Vec<[u8; 3]>, &'static str> {
    if leaf.level > 20 || leaf.index.iter().any(|&i| i >= (1_u32 << leaf.level)) {
        return Err("invalid dyadic leaf address");
    }
    Ok((0..leaf.level)
        .map(|depth| {
            std::array::from_fn(|axis| ((leaf.index[axis] >> (leaf.level - depth - 1)) & 1) as u8)
        })
        .collect())
}

fn child_id(leaf: &LeafKey) -> String {
    let path = child_path(leaf).expect("validated leaf path");
    if path.is_empty() {
        return leaf.root.clone();
    }
    let suffix = path
        .iter()
        .map(|bits| format!("{}{}{}", bits[0], bits[1], bits[2]))
        .collect::<Vec<_>>()
        .join("/");
    format!("{}/{suffix}", leaf.root)
}

fn append_face_equations(
    patch: &FacePatch,
    first: &AnalysisBlock,
    second: &AnalysisBlock,
    max_depth: u8,
    raw: &mut Vec<RawConstraint>,
    uf: &mut UnionFind,
    local_block_indices: &BTreeMap<LeafKey, usize>,
) -> Result<(), &'static str> {
    let first_level = patch.first.leaf.level;
    let second_level = patch.second.leaf.level;
    if first_level == second_level {
        let master = face_net(
            first,
            patch.first.root_side.axis,
            patch.first.root_side.upper,
        )?;
        let slave = face_net(
            second,
            patch.second.root_side.axis,
            patch.second.root_side.upper,
        )?;
        let rows =
            nested_constraints_wide(&master, &slave, [Span::Full, Span::Full], patch.orientation)?;
        if rows.len() != 9 {
            return Err("whole-face trace validation returned an incomplete grid");
        }
        let first_block = local_block_indices[&patch.first.leaf];
        let second_block = local_block_indices[&patch.second.leaf];
        for a in 0..3 {
            for b in 0..3 {
                let (source_a, source_b) = oriented_source_index(
                    patch.orientation.swap,
                    patch.orientation.flip_first,
                    patch.orientation.flip_second,
                    a,
                    b,
                );
                let first_local = first_block * 27
                    + face_local(
                        patch.first.root_side.axis,
                        patch.first.root_side.upper,
                        a,
                        b,
                    );
                let second_local = second_block * 27
                    + face_local(
                        patch.second.root_side.axis,
                        patch.second.root_side.upper,
                        source_a,
                        source_b,
                    );
                if !same_position(
                    first.controls[face_local(
                        patch.first.root_side.axis,
                        patch.first.root_side.upper,
                        a,
                        b,
                    )],
                    second.controls[face_local(
                        patch.second.root_side.axis,
                        patch.second.root_side.upper,
                        source_a,
                        source_b,
                    )],
                )? {
                    return Err("whole-face identified controls disagree physically");
                }
                uf.join(first_local, second_local);
            }
        }
        return Ok(());
    }

    let (
        master_leaf,
        slave_leaf,
        master_block,
        slave_block,
        master_side,
        slave_side,
        subrect,
        orientation,
    ) = if first_level < second_level {
        (
            &patch.first.leaf,
            &patch.second.leaf,
            first,
            second,
            &patch.first.root_side,
            &patch.second.root_side,
            patch.first_subrect,
            patch.orientation,
        )
    } else {
        (
            &patch.second.leaf,
            &patch.first.leaf,
            second,
            first,
            &patch.second.root_side,
            &patch.first.root_side,
            patch.second_subrect,
            inverse_orientation(patch.orientation),
        )
    };
    if slave_leaf.level != master_leaf.level + 1 {
        return Err("adaptive face violates the one-level balance invariant");
    }
    let master_face = face_net(master_block, master_side.axis, master_side.upper)?;
    let slave_face = face_net(slave_block, slave_side.axis, slave_side.upper)?;
    let spans = spans_for_patch(master_leaf, master_side.axis, subrect, max_depth)?;
    let rows = nested_constraints_wide(&master_face, &slave_face, spans, orientation)?;
    if rows.len() != 9 {
        return Err("nested face trace returned an incomplete grid");
    }
    let master_block_index = local_block_indices[master_leaf];
    let slave_block_index = local_block_indices[slave_leaf];
    for (slave_index, row) in rows.into_iter().enumerate() {
        let a = slave_index / 3;
        let b = slave_index % 3;
        let dependent =
            slave_block_index * 27 + face_local(slave_side.axis, slave_side.upper, a, b);
        let terms = row
            .terms
            .into_iter()
            .map(|(master_index, numerator)| {
                let ma = master_index / 3;
                let mb = master_index % 3;
                (
                    master_block_index * 27
                        + face_local(master_side.axis, master_side.upper, ma, mb),
                    numerator,
                )
            })
            .collect();
        raw.push(RawConstraint {
            dependent,
            terms,
            denominator: row.denominator,
        });
    }
    Ok(())
}

fn collapse_identity_rows(raw: &[RawConstraint], uf: &mut UnionFind) -> Result<bool, &'static str> {
    let mut changed = false;
    for row in raw {
        let mut terms = BTreeMap::<usize, i128>::new();
        for &(local, numerator) in &row.terms {
            let root = uf.find(local);
            let value = terms
                .get(&root)
                .copied()
                .unwrap_or(0)
                .checked_add(numerator)
                .ok_or("hanging-node identity sum overflow")?;
            terms.insert(root, value);
        }
        if terms.len() == 1 {
            if let Some((&master, &numerator)) = terms.iter().next() {
                if numerator == row.denominator {
                    changed |= uf.join(row.dependent, master);
                }
            }
        }
    }
    Ok(changed)
}

/// If two face restrictions describe the same coefficient and differ in just
/// two master entries, their difference is c(u_i-u_j)=0. This identifies
/// coarse edge coefficients even when their cells meet only across an edge.
/// Every other relation is retained for the full constraint resolver.
fn collapse_implied_edge_identities(
    raw: &[RawConstraint],
    uf: &mut UnionFind,
) -> Result<bool, &'static str> {
    use crate::constraint_system::Ratio;
    type Weights = BTreeMap<usize, Ratio>;
    let mut grouped: BTreeMap<usize, Vec<Weights>> = BTreeMap::new();
    for row in raw {
        let mut weights = Weights::new();
        for &(id, numerator) in &row.terms {
            let slot = weights.entry(uf.find(id)).or_insert(Ratio::zero());
            *slot = slot.checked_add(Ratio::new(numerator, row.denominator)?)?;
        }
        let rows = grouped.entry(uf.find(row.dependent)).or_default();
        if !rows.contains(&weights) {
            rows.push(weights);
        }
    }
    let mut changed = false;
    for rows in grouped.values() {
        for (j, a) in rows.iter().enumerate() {
            for b in &rows[..j] {
                let mut difference = BTreeMap::new();
                for (&k, &v) in a {
                    let delta = v.checked_sub(b.get(&k).copied().unwrap_or(Ratio::zero()))?;
                    if delta != Ratio::zero() {
                        difference.insert(k, delta);
                    }
                }
                for (&k, &v) in b {
                    if !a.contains_key(&k) {
                        difference.insert(k, Ratio::zero().checked_sub(v)?);
                    }
                }
                if difference.len() == 2 {
                    let mut terms = difference.into_iter();
                    let (i, ci) = terms.next().unwrap();
                    let (k, ck) = terms.next().unwrap();
                    if ci.checked_add(ck)? == Ratio::zero() {
                        changed |= uf.join(i, k);
                    }
                }
            }
        }
    }
    Ok(changed)
}

fn local_to_global(uf: &mut UnionFind, children: &[Child]) -> (Vec<usize>, usize) {
    let mut representative = BTreeMap::<usize, (String, usize)>::new();
    for (block_index, child) in children.iter().enumerate() {
        for local in 0..27 {
            let root = uf.find(block_index * 27 + local);
            let key = (child.block.id.clone(), local);
            representative
                .entry(root)
                .and_modify(|old| *old = old.clone().min(key.clone()))
                .or_insert(key);
        }
    }
    let mut ordered = representative.into_iter().collect::<Vec<_>>();
    ordered.sort_by(|a, b| a.1.cmp(&b.1));
    let ids = ordered
        .iter()
        .enumerate()
        .map(|(id, (root, _))| (*root, id))
        .collect::<BTreeMap<_, _>>();
    let local = (0..children.len() * 27)
        .map(|index| ids[&uf.find(index)])
        .collect();
    (local, ids.len())
}

fn canonical_constraints(
    raw: &[RawConstraint],
    local_global: &[usize],
) -> Result<Vec<(usize, Row)>, &'static str> {
    let mut constraints = Vec::new();
    for raw_row in raw {
        let dependent = local_global[raw_row.dependent];
        let mut terms = BTreeMap::<usize, i128>::new();
        for &(local, numerator) in &raw_row.terms {
            let global = local_global[local];
            let sum = terms
                .get(&global)
                .copied()
                .unwrap_or(0)
                .checked_add(numerator)
                .ok_or("collapsed constraint numerator overflow")?;
            terms.insert(global, sum);
        }
        let sum = terms.values().try_fold(0_i128, |a, &b| {
            a.checked_add(b).ok_or("constraint denominator overflow")
        })?;
        if sum != raw_row.denominator || terms.values().any(|&n| n <= 0) {
            return Err("invalid collapsed hanging-node equation");
        }
        if terms.len() == 1 && terms.get(&dependent) == Some(&raw_row.denominator) {
            continue;
        }
        constraints.push((
            dependent,
            Row {
                terms: terms.into_iter().collect(),
                denominator: raw_row.denominator,
            },
        ));
    }
    constraints.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| a.1.denominator.cmp(&b.1.denominator))
            .then_with(|| a.1.terms.cmp(&b.1.terms))
    });
    constraints.dedup();
    Ok(constraints)
}

fn check_unified_positions(
    children: &[Child],
    local_global: &[usize],
    ndof: usize,
) -> Result<(), &'static str> {
    let mut positions = vec![None; ndof];
    for (block_index, child) in children.iter().enumerate() {
        for (local, &control) in child.block.controls.iter().enumerate() {
            let global = local_global[block_index * 27 + local];
            if let Some(previous) = positions[global] {
                if !same_position(previous, control)? {
                    return Err("identified adaptive controls have different physical positions");
                }
            } else {
                positions[global] = Some(control);
            }
        }
    }
    if positions.iter().any(Option::is_none) {
        return Err("adaptive numbering contains an unused scalar DOF");
    }
    Ok(())
}

fn inherited_boundary_data(
    parent: &AnalysisSpace,
    children: &[Child],
    local_block_indices: &BTreeMap<LeafKey, usize>,
    local_global: &[usize],
    exterior: &[ExteriorPatch],
) -> Result<(BoundaryDofs, BoundaryPatches), &'static str> {
    let root_blocks = parent
        .blocks
        .iter()
        .map(|b| (b.id.as_str(), b))
        .collect::<BTreeMap<_, _>>();
    let mut feature_by_side = BTreeMap::new();
    for root in &parent.blocks {
        for axis in 0..3 {
            for upper in [false, true] {
                if let Some(feature) =
                    root_block_face_feature(root, axis, upper, &parent.boundary_dofs)?
                {
                    feature_by_side.insert((root.id.clone(), axis, upper), feature);
                }
            }
        }
    }
    let mut boundary = parent
        .boundary_dofs
        .keys()
        .map(|feature| (feature.clone(), BTreeSet::new()))
        .collect::<BTreeMap<_, _>>();
    let mut patches = parent
        .boundary_dofs
        .keys()
        .map(|feature| (feature.clone(), Vec::new()))
        .collect::<BTreeMap<_, _>>();
    for exterior_patch in exterior {
        let root_side = &exterior_patch.side.root_side;
        let Some(feature) =
            feature_by_side.get(&(root_side.block.clone(), root_side.axis, root_side.upper))
        else {
            continue;
        };
        if !root_blocks.contains_key(root_side.block.as_str()) {
            return Err("exterior patch references an unknown root selection");
        }
        let block_index = *local_block_indices
            .get(&exterior_patch.side.leaf)
            .ok_or("exterior patch references an unknown leaf")?;
        let child = &children[block_index];
        for a in 0..3 {
            for b in 0..3 {
                let local = block_index * 27 + face_local(root_side.axis, root_side.upper, a, b);
                boundary
                    .get_mut(feature)
                    .ok_or("unknown inherited feature")?
                    .insert(local_global[local]);
            }
        }
        patches
            .get_mut(feature)
            .ok_or("unknown inherited feature")?
            .push(LeafSide {
                leaf: child.leaf.clone(),
                root_side: root_side.clone(),
            });
    }
    if boundary.values().any(BTreeSet::is_empty) {
        return Err("adaptive subdivision lost a physical boundary selection");
    }
    for sides in patches.values_mut() {
        sides.sort();
        sides.dedup();
    }
    Ok((
        boundary
            .into_iter()
            .map(|(name, dofs)| (name, dofs.into_iter().collect()))
            .collect(),
        patches,
    ))
}

/// Build the exactly restricted, conforming-plus-hanging Q2 space described
/// by an audited root `AnalysisSpace` and an A1 dyadic partition.
pub fn build_adaptive_space(
    parent: &AnalysisSpace,
    partition: &AdaptivePartition,
) -> Result<AdaptiveSpace, &'static str> {
    validate_parent(parent, partition)?;
    let report = partition.snapshot()?;
    let max_depth = partition.max_depth();
    let roots = parent
        .blocks
        .iter()
        .map(|block| (block.id.as_str(), block))
        .collect::<BTreeMap<_, _>>();
    let mut children = Vec::with_capacity(report.leaves.len());
    let mut ids = BTreeSet::new();
    for leaf in &report.leaves {
        let root = roots
            .get(leaf.root.as_str())
            .ok_or("adaptive leaf references unknown root block")?;
        let path = child_path(leaf)?;
        let controls = crate::h_refine::restrict_q2_child(&root.controls, &path)?;
        let id = child_id(leaf);
        if !ids.insert(id.clone()) {
            return Err("duplicate adaptive child block identity");
        }
        children.push(Child {
            leaf: leaf.clone(),
            block: AnalysisBlock {
                id,
                controls,
                dofs: Vec::new(),
            },
        });
    }
    children.sort_by(|a, b| a.leaf.cmp(&b.leaf));
    let local_block_indices = children
        .iter()
        .enumerate()
        .map(|(index, child)| (child.leaf.clone(), index))
        .collect::<BTreeMap<_, _>>();
    let local_count = children
        .len()
        .checked_mul(27)
        .ok_or("adaptive scalar local size overflow")?;
    let mut uf = UnionFind::new(local_count);
    let mut raw_constraints = Vec::new();
    for patch in &report.interfaces {
        let first = &children[*local_block_indices
            .get(&patch.first.leaf)
            .ok_or("interface references unknown first leaf")?]
        .block;
        let second = &children[*local_block_indices
            .get(&patch.second.leaf)
            .ok_or("interface references unknown second leaf")?]
        .block;
        if patch.first.root_side.block != patch.first.leaf.root
            || patch.second.root_side.block != patch.second.leaf.root
        {
            return Err("interface semantic root side does not match its leaf");
        }
        append_face_equations(
            patch,
            first,
            second,
            max_depth,
            &mut raw_constraints,
            &mut uf,
            &local_block_indices,
        )?;
    }
    loop {
        while collapse_identity_rows(&raw_constraints, &mut uf)? {}
        if !collapse_implied_edge_identities(&raw_constraints, &mut uf)? {
            break;
        }
    }
    let (local_global, ndof) = local_to_global(&mut uf, &children);
    check_unified_positions(&children, &local_global, ndof)?;
    let constraints = canonical_constraints(&raw_constraints, &local_global)?;
    let equations = constraints.clone();
    let lift = constraint_system::resolve(ndof, &equations, &[])?;
    let mut blocks = Vec::with_capacity(children.len());
    let mut leaf_blocks = BTreeMap::new();
    for (block_index, mut child) in children.iter().cloned().enumerate() {
        child.block.dofs = (0..27)
            .map(|local| local_global[block_index * 27 + local])
            .collect();
        leaf_blocks.insert(child.leaf.clone(), child.block.id.clone());
        blocks.push(child.block);
    }
    let (boundary_dofs, boundary_patches) = inherited_boundary_data(
        parent,
        &children,
        &local_block_indices,
        &local_global,
        &report.exterior,
    )?;
    let space = AnalysisSpace {
        ndof,
        blocks,
        boundary_dofs,
    };
    Ok(AdaptiveSpace {
        space,
        lift,
        leaf_blocks,
        constraints,
        boundary_patches,
        root_ndof: parent.ndof,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        adaptive_partition::LeafKey,
        bracket::{bracket_blocks, BracketBlock},
        bracket_audit::audit_bracket,
        incidence::{enumerate_faces, Incidence},
    };

    fn d4_roots(orientation: crate::incidence::Orientation) -> (AnalysisSpace, AdaptivePartition) {
        use crate::incidence::{BlockSide, Interface};
        let mut controls_a = Vec::new();
        let mut controls_b = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    controls_a.push([i as i128, j as i128, k as i128, 2]);
                    let (a, b) = oriented_source_index(
                        inverse_orientation(orientation).swap,
                        inverse_orientation(orientation).flip_first,
                        inverse_orientation(orientation).flip_second,
                        j,
                        k,
                    );
                    controls_b.push([2 + i as i128, a as i128, b as i128, 2]);
                }
            }
        }
        let root_block = |id: &str, controls: Vec<[i128; 4]>, base: usize| AnalysisBlock {
            id: id.into(),
            controls,
            dofs: (base..base + 27).collect(),
        };
        let roots = vec![
            root_block("a", controls_a, 0),
            root_block("b", controls_b, 27),
        ];
        let incidence = Incidence {
            interfaces: vec![Interface {
                first: BlockSide {
                    block: "a".into(),
                    axis: 0,
                    upper: true,
                },
                second: BlockSide {
                    block: "b".into(),
                    axis: 0,
                    upper: false,
                },
                orientation,
            }],
            exterior: [
                BlockSide {
                    block: "a".into(),
                    axis: 0,
                    upper: false,
                },
                BlockSide {
                    block: "a".into(),
                    axis: 1,
                    upper: false,
                },
                BlockSide {
                    block: "a".into(),
                    axis: 1,
                    upper: true,
                },
                BlockSide {
                    block: "a".into(),
                    axis: 2,
                    upper: false,
                },
                BlockSide {
                    block: "a".into(),
                    axis: 2,
                    upper: true,
                },
                BlockSide {
                    block: "b".into(),
                    axis: 0,
                    upper: true,
                },
                BlockSide {
                    block: "b".into(),
                    axis: 1,
                    upper: false,
                },
                BlockSide {
                    block: "b".into(),
                    axis: 1,
                    upper: true,
                },
                BlockSide {
                    block: "b".into(),
                    axis: 2,
                    upper: false,
                },
                BlockSide {
                    block: "b".into(),
                    axis: 2,
                    upper: true,
                },
            ]
            .into(),
        };
        let partition =
            AdaptivePartition::new(&["a".into(), "b".into()], &incidence, 2, 100).unwrap();
        let parent = AnalysisSpace {
            ndof: 54,
            blocks: roots,
            boundary_dofs: BTreeMap::new(),
        };
        (parent, partition)
    }

    fn root_space(blocks: &[BracketBlock]) -> (AnalysisSpace, Incidence, AdaptivePartition) {
        let certificate = audit_bracket(blocks).unwrap();
        let space = crate::analysis_space::from_certified_faces(
            blocks,
            &certificate.incidence,
            &certificate.design_faces,
        )
        .unwrap();
        let incidence = enumerate_faces(blocks).unwrap();
        let roots = blocks.iter().map(|b| b.id.clone()).collect::<Vec<_>>();
        let partition = AdaptivePartition::new(&roots, &incidence, 3, 50_000).unwrap();
        (space, incidence, partition)
    }

    fn check_affine_lift(adaptive: &AdaptiveSpace, axis: usize) {
        let mut full: Vec<Option<f64>> = vec![None; adaptive.space.ndof];
        for block in &adaptive.space.blocks {
            for (local, control) in block.controls.iter().enumerate() {
                let dof = block.dofs[local];
                let value = control[axis] as f64 / control[3] as f64;
                if let Some(previous) = full[dof] {
                    assert!((previous - value).abs() < 1e-12);
                } else {
                    full[dof] = Some(value);
                }
            }
        }
        let full = full.into_iter().map(Option::unwrap).collect::<Vec<_>>();
        let reduced = adaptive
            .lift
            .free_roots
            .iter()
            .map(|&dof| full[dof])
            .collect::<Vec<_>>();
        for line in &adaptive.lift.lines {
            let value = line
                .terms
                .iter()
                .map(|(i, coefficient)| {
                    reduced[*i] * coefficient.numerator() as f64 / coefficient.denominator() as f64
                })
                .sum::<f64>()
                + line.offset.numerator() as f64 / line.offset.denominator() as f64;
            assert!((value - full[line.dof]).abs() < 1e-10);
        }
    }

    #[test]
    fn h0_matches_certified_space_and_has_identity_lift() {
        let blocks = bracket_blocks().unwrap();
        let (root, _, partition) = root_space(&blocks);
        let adaptive = build_adaptive_space(&root, &partition).unwrap();
        assert_eq!(adaptive.space.ndof, 531);
        assert_eq!(adaptive.space.blocks.len(), 36);
        assert_eq!(adaptive.space, root);
        assert_eq!(adaptive.space.boundary_dofs, root.boundary_dofs);
        assert!(adaptive.constraints.is_empty());
        assert_eq!(adaptive.lift.free_roots.len(), root.ndof);
        assert!(adaptive.space.blocks.iter().all(|b| b.controls.len() == 27));
    }

    #[test]
    fn uniform_h1_matches_uniform_refinement_and_preserves_affines() {
        let blocks = bracket_blocks().unwrap();
        let (root, _, mut partition) = root_space(&blocks);
        let requests = partition.leaves().iter().cloned().collect::<Vec<_>>();
        let report = partition.refine(&requests).unwrap();
        let adaptive = build_adaptive_space(&root, &partition).unwrap();
        let reference = crate::h_refine::refine_uniform(&root).unwrap();
        assert_eq!(report.leaves.len(), 36 * 8);
        assert_eq!(adaptive.space.ndof, 3_215);
        assert_eq!(adaptive.space.ndof, reference.ndof);
        assert_eq!(adaptive.space.blocks.len(), reference.blocks.len());
        assert!(adaptive.constraints.is_empty());
        for axis in 0..3 {
            check_affine_lift(&adaptive, axis);
        }
    }

    #[test]
    fn edited_bracket_h1_preserves_certified_selections_and_affines() {
        let edited = crate::bracket_parameters::build_bracket_parameters(56, 64).unwrap();
        let parent =
            crate::analysis_space::parameterized_bracket_space(56, 64, &edited.blocks).unwrap();
        let incidence = enumerate_faces(&edited.blocks).unwrap();
        let roots = edited
            .blocks
            .iter()
            .map(|b| b.id.clone())
            .collect::<Vec<_>>();
        let mut partition = AdaptivePartition::new(&roots, &incidence, 2, 10_000).unwrap();
        let requests = partition.leaves().iter().cloned().collect::<Vec<_>>();
        partition.refine(&requests).unwrap();
        let adaptive = build_adaptive_space(&parent, &partition).unwrap();
        assert_eq!(adaptive.space.ndof, 3_215);
        assert_eq!(
            adaptive.space.boundary_dofs.keys().collect::<Vec<_>>(),
            parent.boundary_dofs.keys().collect::<Vec<_>>()
        );
        for axis in 0..3 {
            check_affine_lift(&adaptive, axis);
        }
    }

    #[test]
    fn local_interface_constraints_are_exact_and_consistent_at_edges() {
        let blocks = bracket_blocks().unwrap();
        let (root, _, mut partition) = root_space(&blocks);
        let target = root
            .blocks
            .iter()
            .find(|b| b.id.contains("hole_x"))
            .unwrap();
        let selected = LeafKey::root(target.id.clone());
        let report = partition.refine(&[selected]).unwrap();
        let adaptive = build_adaptive_space(&root, &partition).unwrap();
        assert_eq!(report.leaves.len(), 43);
        assert!(!adaptive.constraints.is_empty());
        assert!(adaptive.lift.free_roots.len() < adaptive.space.ndof);
        assert!(adaptive.boundary_patches.contains_key("x60"));
        assert!(adaptive.boundary_patches.contains_key("y50"));
        let hanging_patch = report
            .interfaces
            .iter()
            .find(|face| face.first.leaf.level != face.second.leaf.level)
            .unwrap();
        let (slave_leaf, slave_side) =
            if hanging_patch.first.leaf.level > hanging_patch.second.leaf.level {
                (&hanging_patch.first.leaf, &hanging_patch.first.root_side)
            } else {
                (&hanging_patch.second.leaf, &hanging_patch.second.root_side)
            };
        let slave_id = &adaptive.leaf_blocks[slave_leaf];
        let slave = adaptive
            .space
            .blocks
            .iter()
            .find(|b| &b.id == slave_id)
            .unwrap();
        let mut endpoint_rows = Vec::new();
        for a in [0, 2] {
            for b in [0, 2] {
                let dof = slave.dofs[face_local(slave_side.axis, slave_side.upper, a, b)];
                if let Some(row) = adaptive
                    .constraints
                    .iter()
                    .find(|(dependent, _)| *dependent == dof)
                {
                    endpoint_rows.push(row);
                }
            }
        }
        assert!(
            endpoint_rows.iter().any(|(_, row)| row.terms.len() > 1),
            "a fine face endpoint on the coarse midpoint must remain a weighted restriction"
        );
        for axis in 0..3 {
            check_affine_lift(&adaptive, axis);
        }
    }

    #[test]
    fn parameterized_hole_local_level_two_keeps_wide_exact_face_restrictions() {
        let geometry = crate::bracket_parameters::build_bracket_parameters(56, 64).unwrap();
        let parent =
            crate::analysis_space::parameterized_bracket_space(56, 64, &geometry.blocks).unwrap();
        let incidence = enumerate_faces(&geometry.blocks).unwrap();
        let roots = geometry
            .blocks
            .iter()
            .map(|b| b.id.clone())
            .collect::<Vec<_>>();
        let mut partition = AdaptivePartition::new(&roots, &incidence, 4, 1500).unwrap();
        let focus = LeafKey::root("hole_x_sector_east_z10_30");
        partition.refine(std::slice::from_ref(&focus)).unwrap();
        partition
            .refine(&[LeafKey {
                root: focus.root,
                level: 1,
                index: [1, 0, 0],
            }])
            .unwrap();
        let adaptive = build_adaptive_space(&parent, &partition).unwrap();
        assert!(!adaptive.constraints.is_empty());
        assert!(adaptive
            .leaf_blocks
            .keys()
            .any(|leaf| leaf.level == 2 && leaf.index == [2, 0, 0]));
        for axis in 0..3 {
            check_affine_lift(&adaptive, axis);
        }
    }

    #[test]
    fn edge_identification_requires_an_exact_two_term_difference() {
        let trace = |terms| RawConstraint {
            dependent: 4,
            terms,
            denominator: 2,
        };
        let mut groups = UnionFind::new(5);
        let a = trace(vec![(0, 1), (1, 1)]);
        let b = trace(vec![(2, 1), (3, 1)]);
        assert!(!collapse_implied_edge_identities(&[a, b], &mut groups).unwrap());
        let a = trace(vec![(0, 1), (1, 1)]);
        let b = trace(vec![(0, 1), (3, 1)]);
        assert!(collapse_implied_edge_identities(&[a, b], &mut groups).unwrap());
        assert_eq!(groups.find(1), groups.find(3));
        assert_ne!(groups.find(0), groups.find(1));
        assert_ne!(groups.find(2), groups.find(1));
    }

    #[test]
    fn altered_or_malformed_parent_state_is_refused() {
        let blocks = bracket_blocks().unwrap();
        let (root, _, partition) = root_space(&blocks);
        let mut altered = root.clone();
        altered.blocks[0].controls[0][3] = 0;
        assert!(build_adaptive_space(&altered, &partition).is_err());
        let mut wrong_ids = root;
        wrong_ids.blocks[0].dofs[0] = wrong_ids.ndof;
        assert!(build_adaptive_space(&wrong_ids, &partition).is_err());
    }

    #[test]
    fn hanging_restrictions_are_valid_for_all_eight_d4_orientations() {
        for swap in [false, true] {
            for flip_first in [false, true] {
                for flip_second in [false, true] {
                    let orientation = crate::incidence::Orientation {
                        swap,
                        flip_first,
                        flip_second,
                    };
                    let (parent, mut partition) = d4_roots(orientation);
                    let before = build_adaptive_space(&parent, &partition).unwrap();
                    assert_eq!(before.space.ndof, 45);
                    assert!(before.constraints.is_empty());
                    partition.refine(&[LeafKey::root("a")]).unwrap();
                    let adaptive = build_adaptive_space(&parent, &partition).unwrap();
                    assert_eq!(adaptive.space.blocks.len(), 9);
                    assert!(!adaptive.constraints.is_empty());
                    for axis in 0..3 {
                        check_affine_lift(&adaptive, axis);
                    }
                }
            }
        }
    }
}
