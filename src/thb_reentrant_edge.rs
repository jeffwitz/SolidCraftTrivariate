//! Exact multilevel THB spaces for a scalar L-section corner benchmark.
//!
//! Three positively oriented affine charts cover the reentrant section. The
//! module is deliberately limited to this chart arrangement and to cubic
//! tensor B-splines; it is not a multipatch adaptation API.

use crate::{
    constraint_system::Ratio as Q,
    parent_field::axis_rows,
    spline_constraints::{add, Row},
};
use std::collections::{BTreeMap, BTreeSet};

type Result<T> = std::result::Result<T, &'static str>;
type Sparse = Vec<Row>;
const DEGREE: usize = 3;
const CHARTS: [Chart; 3] = [
    Chart {
        origin: [0., 0.],
        u: [-1., 0.],
        v: [0., -1.],
    },
    Chart {
        origin: [0., 0.],
        u: [0., -1.],
        v: [1., 0.],
    },
    Chart {
        origin: [0., 0.],
        u: [0., 1.],
        v: [-1., 0.],
    },
];

#[derive(Clone, Copy)]
struct Chart {
    origin: [f64; 2],
    u: [f64; 2],
    v: [f64; 2],
}

#[derive(Clone)]
struct BasisFunction {
    source_level: usize,
    source_id: usize,
    by_level: Vec<BTreeMap<usize, Q>>,
    global: usize,
}

#[derive(Clone)]
struct Patch {
    functions: Vec<BasisFunction>,
    shapes: Vec<[usize; 2]>,
}

#[derive(Clone, Copy, Debug)]
struct LeafCell {
    patch: usize,
    level: usize,
    interval: [[Q; 2]; 2],
}

/// One integration cell and its exact-extraction THB basis.
#[derive(Clone, Debug)]
pub struct IntegrationLeaf {
    pub patch: usize,
    pub level: usize,
    /// Parameter bounds `[[umin,umax],[vmin,vmax]]`.
    pub uv: [[f64; 2]; 2],
    pub dofs: Vec<usize>,
    /// `bernstein[4*i+j][local_dof]` is one extracted basis coefficient.
    pub bernstein: Vec<Vec<f64>>,
}

/// A locally refined, exactly nested scalar basis on the L section.
#[derive(Clone)]
pub struct ReentrantHierarchy {
    pub maximum_level: usize,
    pub active_scalar_dofs: usize,
    pub uniform_fine_scalar_dofs: usize,
    pub leaves: usize,
    pub exact_seam_trace_pairs: usize,
    is_uniform: bool,
    widths: Vec<usize>,
    patches: Vec<Patch>,
    knots: Vec<[Vec<Q>; 2]>,
    inclusions: Vec<[Sparse; 2]>,
    leaves_data: Vec<LeafCell>,
}

fn q(n: i64, d: i64) -> Q {
    Q::new(i128::from(n), i128::from(d)).expect("small rational")
}

fn rational_add(row: &mut Row, index: usize, value: Q) -> Result<()> {
    let sum = row
        .get(&index)
        .copied()
        .unwrap_or(Q::zero())
        .checked_add(value)?;
    if sum == Q::zero() {
        row.remove(&index);
    } else {
        row.insert(index, sum);
    }
    Ok(())
}

fn index(shape: [usize; 2], i: usize, j: usize) -> usize {
    i * shape[1] + j
}

fn intervals(knots: &[Q]) -> Vec<[Q; 2]> {
    knots
        .windows(2)
        .filter_map(|p| p[0].cmp_value(&p[1]).is_lt().then_some([p[0], p[1]]))
        .collect()
}

fn refine_knots(knots: &[Q]) -> Result<Vec<Q>> {
    let mut result = knots.to_vec();
    for [a, b] in intervals(knots) {
        result.push(a.checked_add(b)?.checked_div(q(2, 1))?);
    }
    result.sort_by(Q::cmp_value);
    Ok(result)
}

fn base_knots() -> Vec<Q> {
    let mut knots = vec![Q::zero(); 4];
    knots.extend((1..8).map(|i| q(i, 8)));
    knots.extend(vec![Q::one(); 4]);
    knots
}

fn shape(knots: &[Vec<Q>; 2]) -> [usize; 2] {
    [knots[0].len() - 4, knots[1].len() - 4]
}

fn inverse(matrix: &[Vec<Q>]) -> Result<Vec<Vec<Q>>> {
    let n = matrix.len();
    if n != DEGREE + 1 || matrix.iter().any(|r| r.len() != n) {
        return Err("THB local extraction matrix has wrong shape");
    }
    let mut a = vec![vec![Q::zero(); 2 * n]; n];
    for i in 0..n {
        a[i][..n].copy_from_slice(&matrix[i]);
        a[i][n + i] = Q::one();
    }
    for c in 0..n {
        let p = (c..n)
            .find(|&r| a[r][c] != Q::zero())
            .ok_or("THB local extraction is singular")?;
        a.swap(c, p);
        let pivot = Q::new(a[c][c].denominator(), a[c][c].numerator())?;
        for v in &mut a[c] {
            *v = v.checked_mul(pivot)?;
        }
        for r in 0..n {
            if r == c || a[r][c] == Q::zero() {
                continue;
            }
            let factor = Q::zero().checked_sub(a[r][c])?;
            let row = a[c].clone();
            for (dst, src) in a[r].iter_mut().zip(row) {
                *dst = dst.checked_add(factor.checked_mul(src)?)?;
            }
        }
    }
    Ok(a.into_iter().map(|r| r[n..].to_vec()).collect())
}

fn axis_inclusion(old: &[Q], new: &[Q]) -> Result<Sparse> {
    let old_count = old.len() - 4;
    let new_count = new.len() - 4;
    let mut rows = vec![Row::new(); new_count];
    for cell in intervals(new) {
        let old_rows = axis_rows(old, DEGREE, cell)?;
        let new_rows = axis_rows(new, DEGREE, cell)?;
        let active: Vec<_> = new_rows
            .iter()
            .flat_map(|r| r.keys().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        if active.len() != DEGREE + 1 {
            return Err("THB fine interval has an invalid active count");
        }
        let matrix: Vec<Vec<Q>> = new_rows
            .iter()
            .map(|r| active.iter().map(|i| r.get(i).copied().unwrap_or(Q::zero())).collect())
            .collect();
        let inv = inverse(&matrix)?;
        for (j, &fine_id) in active.iter().enumerate() {
            let mut coeffs = Row::new();
            for (r, old_row) in old_rows.iter().enumerate() {
                let scale = inv[j][r];
                for (&id, &value) in old_row {
                    add(&mut coeffs, id, value.checked_mul(scale)?)?;
                }
            }
            if rows[fine_id].is_empty() {
                rows[fine_id] = coeffs;
            } else if rows[fine_id] != coeffs {
                return Err("THB inclusion differs across knot spans");
            }
        }
    }
    if rows.len() != new_count || rows.iter().any(|r| r.keys().any(|i| *i >= old_count)) {
        return Err("invalid exact THB axis inclusion");
    }
    Ok(rows)
}


fn tensor_support(knots: &[Vec<Q>; 2], id: usize, shape: [usize; 2]) -> [[Q; 2]; 2] {
    let i = id / shape[1];
    let j = id % shape[1];
    [[knots[0][i], knots[0][i + 4]], [knots[1][j], knots[1][j + 4]]]
}

fn inside_region(support: [[Q; 2]; 2], end: Q) -> bool {
    (0..2).all(|a| {
        !support[a][0].cmp_value(&Q::zero()).is_lt()
            && !support[a][1].cmp_value(&end).is_gt()
    })
}

fn all_inside(knots: &[Vec<Q>; 2], id: usize, end: Q) -> bool {
    inside_region(tensor_support(knots, id, shape(knots)), end)
}

fn region_end(level: usize, widths: &[usize]) -> Q {
    if level == 0 {
        Q::one()
    } else {
        q(widths[level - 1] as i64, (8_i64) << (level - 1))
    }
}

fn make_patch(
    knots: &[[Vec<Q>; 2]],
    inclusions: &[[Sparse; 2]],
    max_level: usize,
    widths: &[usize],
) -> Result<Patch> {
    let shapes: Vec<_> = knots.iter().map(shape).collect();
    let mut functions = Vec::new();
    for level in 0..=max_level {
        let count = shapes[level][0] * shapes[level][1];
        for id in 0..count {
            let active = if level == 0 {
                max_level == 0 || !all_inside(&knots[0], id, region_end(1, widths))
            } else if level < max_level {
                all_inside(&knots[level], id, region_end(level, widths))
                    && !all_inside(&knots[level], id, region_end(level + 1, widths))
            } else {
                all_inside(&knots[level], id, region_end(level, widths))
            };
            if !active {
                continue;
            }
            let mut by_level = vec![BTreeMap::new(); max_level + 1];
            by_level[level].insert(id, Q::one());
            for next in level + 1..=max_level {
                let mut expanded_2d = BTreeMap::new();
                for fine_i in 0..inclusions[next - 1][0].len() {
                    for fine_j in 0..inclusions[next - 1][1].len() {
                        let fine_id = index(shapes[next], fine_i, fine_j);
                        let mut coefficient = Q::zero();
                        for (&old_i, &ci) in &inclusions[next - 1][0][fine_i] {
                            for (&old_j, &cj) in &inclusions[next - 1][1][fine_j] {
                                if let Some(&v) = by_level[next - 1]
                                    .get(&index(shapes[next - 1], old_i, old_j))
                                {
                                    coefficient = coefficient.checked_add(ci.checked_mul(cj)?.checked_mul(v)?)?;
                                }
                            }
                        }
                        if coefficient != Q::zero()
                            && !all_inside(&knots[next], fine_id, region_end(next, widths))
                        {
                            expanded_2d.insert(fine_id, coefficient);
                        }
                    }
                }
                by_level[next] = expanded_2d;
            }
            functions.push(BasisFunction {
                source_level: level,
                source_id: id,
                by_level,
                global: 0,
            });
        }
    }
    Ok(Patch { functions, shapes })
}

/// Build the fixed, dyadically nested THB hierarchy. Levels are 0 through 3.
pub fn build_hierarchy(maximum_level: usize) -> Result<ReentrantHierarchy> {
    if maximum_level > 3 {
        return Err("reentrant-edge THB pilot supports maximum level 3");
    }
    let widths = [4; 3];
    build_hierarchy_with_widths(maximum_level, &widths[..maximum_level])
}

/// Build the same hierarchy with 4–8 prepared spans along each nested box.
pub fn build_hierarchy_with_widths(
    maximum_level: usize,
    widths: &[usize],
) -> Result<ReentrantHierarchy> {
    if maximum_level > 3 {
        return Err("reentrant-edge THB pilot supports maximum level 3");
    }
    if widths.len() != maximum_level || widths.iter().any(|w| !(4..=8).contains(w)) {
        return Err("THB box widths must provide 4..=8 spans at every active level");
    }
    let mut axes = vec![[base_knots(), base_knots()]];
    let mut inclusions = Vec::new();
    for _ in 0..maximum_level {
        let previous = axes.last().ok_or("missing base knots")?;
        let next = [refine_knots(&previous[0])?, refine_knots(&previous[1])?];
        let axis0 = axis_inclusion(&previous[0], &next[0])?;
        let axis1 = axis_inclusion(&previous[1], &next[1])?;
        let patch_inclusion = [axis0, axis1];
        inclusions.push(patch_inclusion);
        axes.push(next);
    }
    let mut patches = Vec::new();
    for _ in 0..3 {
        patches.push(make_patch(&axes, &inclusions, maximum_level, widths)?);
    }
    let (global_count, seam_pairs) = identify_traces(&mut patches, maximum_level)?;
    let leaves_data = make_leaves(&axes, maximum_level, widths)?;
    let uniform_count = build_uniform(maximum_level)?.active_scalar_dofs;
    Ok(ReentrantHierarchy {
        maximum_level,
        active_scalar_dofs: global_count,
        uniform_fine_scalar_dofs: uniform_count,
        leaves: leaves_data.len(),
        exact_seam_trace_pairs: seam_pairs,
        is_uniform: false,
        widths: widths.to_vec(),
        patches,
        knots: axes,
        inclusions,
        leaves_data,
    })
}

fn trace_vector(function: &BasisFunction, shape: [usize; 2], edge: usize) -> Vec<Q> {
    let n = shape[0].max(shape[1]);
    let mut result = vec![Q::zero(); n];
    let coeffs = function.by_level.last().expect("level data");
    for (&id, &value) in coeffs {
        let i = id / shape[1];
        let j = id % shape[1];
        let (normal, tangent) = match edge {
            0 => (i == 0, j),
            1 => (i + 1 == shape[0], j),
            2 => (j == 0, i),
            _ => (j + 1 == shape[1], i),
        };
        if normal {
            result[tangent] = result[tangent].checked_add(value).expect("small trace sum");
        }
    }
    while result.last() == Some(&Q::zero()) {
        result.pop();
    }
    result
}

fn identify_traces(patches: &mut [Patch], level: usize) -> Result<(usize, usize)> {
    let local_offsets: Vec<usize> = (0..patches.len())
        .scan(0, |sum, p| {
            let start = *sum;
            *sum += patches[p].functions.len();
            Some(start)
        })
        .collect();
    let total: usize = patches.iter().map(|p| p.functions.len()).sum();
    let mut parent: Vec<usize> = (0..total).collect();
    fn find(parent: &mut [usize], i: usize) -> usize {
        if parent[i] != i {
            parent[i] = find(parent, parent[i]);
        }
        parent[i]
    }
    fn union(parent: &mut [usize], a: usize, b: usize) {
        let (ra, rb) = (find(parent, a), find(parent, b));
        if ra != rb {
            parent[rb] = ra;
        }
    }
    let mut pairs = 0;
    for (a, edge_a, b, edge_b) in [(0, 0, 1, 2), (0, 2, 2, 0)] {
        let mut left = BTreeMap::<Vec<Q>, usize>::new();
        let mut right = BTreeMap::<Vec<Q>, usize>::new();
        for (i, function) in patches[a].functions.iter().enumerate() {
            let trace = trace_vector(function, patches[a].shapes[level], edge_a);
            if trace.iter().any(|x| *x != Q::zero()) && left.insert(trace, i).is_some() {
                return Err("duplicate nonzero THB trace on a chart edge");
            }
        }
        for (i, function) in patches[b].functions.iter().enumerate() {
            let trace = trace_vector(function, patches[b].shapes[level], edge_b);
            if trace.iter().any(|x| *x != Q::zero()) && right.insert(trace, i).is_some() {
                return Err("duplicate nonzero THB trace on a chart edge");
            }
        }
        if left.len() != right.len() || left.keys().ne(right.keys()) {
            return Err("THB chart interface trace spaces differ");
        }
        for (trace, &ia) in &left {
            let ib = right[trace];
            union(&mut parent, local_offsets[a] + ia, local_offsets[b] + ib);
            pairs += 1;
        }
    }
    let roots: BTreeSet<_> = (0..total).map(|i| find(&mut parent, i)).collect();
    let map: BTreeMap<_, _> = roots.iter().enumerate().map(|(i, &r)| (r, i)).collect();
    for p in 0..patches.len() {
        for i in 0..patches[p].functions.len() {
            let root = find(&mut parent, local_offsets[p] + i);
            patches[p].functions[i].global = map[&root];
        }
    }
    Ok((roots.len(), pairs))
}

fn make_leaves(axes: &[[Vec<Q>; 2]], max_level: usize, widths: &[usize]) -> Result<Vec<LeafCell>> {
    let mut leaves = Vec::new();
    for patch in 0..3 {
        for (level, axis) in axes.iter().enumerate().take(max_level + 1) {
            let u = intervals(&axis[0]);
            let v = intervals(&axis[1]);
            for a in &u {
                for b in &v {
                    let cell = [*a, *b];
                    let inside = inside_region(cell, region_end(level, widths));
                    let include = if max_level == 0 {
                        true
                    } else if level == 0 {
                        !inside_region(cell, region_end(1, widths))
                    } else if level == max_level {
                        inside
                    } else {
                        inside && !inside_region(cell, region_end(level + 1, widths))
                    };
                    if include {
                        leaves.push(LeafCell { patch, level, interval: cell });
                    }
                }
            }
        }
    }
    if leaves.is_empty() || leaves.len() > 8192 {
        return Err("THB integration leaf count exceeds its checked budget");
    }
    Ok(leaves)
}

/// Build the corresponding globally uniform cubic space at a selected level.
pub fn build_uniform(level: usize) -> Result<ReentrantHierarchy> {
    if level > 3 {
        return Err("reentrant-edge reference supports maximum level 3");
    }
    let mut axes = vec![[base_knots(), base_knots()]];
    let mut inclusions = Vec::new();
    for _ in 0..level {
        let previous = axes.last().ok_or("missing base knots")?;
        let next = [refine_knots(&previous[0])?, refine_knots(&previous[1])?];
        inclusions.push([
            axis_inclusion(&previous[0], &next[0])?,
            axis_inclusion(&previous[1], &next[1])?,
        ]);
        axes.push(next);
    }
    let final_shape = shape(&axes[level]);
    let mut patches = Vec::new();
    for _ in 0..3 {
        let mut functions = Vec::new();
        for id in 0..final_shape[0] * final_shape[1] {
            let mut by_level = vec![BTreeMap::new(); level + 1];
            by_level[level].insert(id, Q::one());
            functions.push(BasisFunction {
                source_level: level,
                source_id: id,
                by_level,
                global: 0,
            });
        }
        patches.push(Patch {
            functions,
            shapes: axes.iter().map(shape).collect(),
        });
    }
    let (global_count, seam_pairs) = identify_traces(&mut patches, level)?;
    let mut leaves_data = Vec::new();
    for patch in 0..3 {
        for u in intervals(&axes[level][0]) {
            for v in intervals(&axes[level][1]) {
                leaves_data.push(LeafCell { patch, level, interval: [u, v] });
            }
        }
    }
    Ok(ReentrantHierarchy {
        maximum_level: level,
        active_scalar_dofs: global_count,
        uniform_fine_scalar_dofs: global_count,
        leaves: leaves_data.len(),
        exact_seam_trace_pairs: seam_pairs,
        is_uniform: true,
        widths: Vec::new(),
        patches,
        knots: axes,
        inclusions,
        leaves_data,
    })
}

impl ReentrantHierarchy {
    /// Return compact cellwise Bernstein extractions for independent integration and solving.
    pub fn integration_leaves(&self) -> Result<Vec<IntegrationLeaf>> {
        let mut tensor_index = Vec::with_capacity(self.patches.len());
        for patch in &self.patches {
            let mut levels = Vec::with_capacity(self.maximum_level + 1);
            for level in 0..=self.maximum_level {
                let shape = patch.shapes[level];
                let mut index = vec![Vec::<(usize, Q)>::new(); shape[0] * shape[1]];
                for (function_id, function) in patch.functions.iter().enumerate() {
                    for (&tensor_id, &coefficient) in &function.by_level[level] {
                        index[tensor_id].push((function_id, coefficient));
                    }
                }
                levels.push(index);
            }
            tensor_index.push(levels);
        }
        let mut axis_cache = BTreeMap::<(usize, usize, Q, Q), Vec<Row>>::new();
        let mut result = Vec::with_capacity(self.leaves_data.len());
        for cell in &self.leaves_data {
            let patch = self.patches.get(cell.patch).ok_or("unknown THB chart")?;
            let u_key = (cell.level, 0, cell.interval[0][0], cell.interval[0][1]);
            let v_key = (cell.level, 1, cell.interval[1][0], cell.interval[1][1]);
            if let std::collections::btree_map::Entry::Vacant(entry) = axis_cache.entry(u_key) {
                entry.insert(axis_rows(&self.knots[cell.level][0], DEGREE, cell.interval[0])?);
            }
            if let std::collections::btree_map::Entry::Vacant(entry) = axis_cache.entry(v_key) {
                entry.insert(axis_rows(&self.knots[cell.level][1], DEGREE, cell.interval[1])?);
            }
            let u_rows = &axis_cache[&u_key];
            let v_rows = &axis_cache[&v_key];
            let active_u: BTreeSet<_> = u_rows.iter().flat_map(|r| r.keys().copied()).collect();
            let active_v: BTreeSet<_> = v_rows.iter().flat_map(|r| r.keys().copied()).collect();
            let mut local = BTreeMap::<usize, BTreeMap<usize, Q>>::new();
            for &iu in &active_u {
                for &iv in &active_v {
                    let tensor_id = index(patch.shapes[cell.level], iu, iv);
                    for &(function_id, coefficient) in &tensor_index[cell.patch][cell.level][tensor_id] {
                        rational_add(local.entry(function_id).or_default(), tensor_id, coefficient)?;
                    }
                }
            }
            let mut dofs = Vec::new();
            let mut bernstein = vec![Vec::new(); 16];
            for (function_id, coefficients) in local {
                let mut extraction = [Q::zero(); 16];
                for i in 0..=DEGREE {
                    for j in 0..=DEGREE {
                        let mut value = Q::zero();
                        for (&tensor_id, &coefficient) in &coefficients {
                            let iu = tensor_id / patch.shapes[cell.level][1];
                            let iv = tensor_id % patch.shapes[cell.level][1];
                            let a = u_rows[i].get(&iu).copied().unwrap_or(Q::zero());
                            let b = v_rows[j].get(&iv).copied().unwrap_or(Q::zero());
                            if a != Q::zero() && b != Q::zero() {
                                value = value.checked_add(coefficient.checked_mul(a)?.checked_mul(b)?)?;
                            }
                        }
                        extraction[i * 4 + j] = value;
                    }
                }
                if extraction.iter().all(|x| *x == Q::zero()) {
                    continue;
                }
                dofs.push(patch.functions[function_id].global);
                for (row, value) in bernstein.iter_mut().zip(extraction) {
                    row.push(q_value(value));
                }
            }
            if dofs.is_empty() || bernstein.iter().any(|row| row.len() != dofs.len()) {
                return Err("THB integration leaf has an invalid active extraction");
            }
            result.push(IntegrationLeaf {
                patch: cell.patch,
                level: cell.level,
                uv: cell.interval.map(|x| x.map(q_value)),
                dofs,
                bernstein,
            });
        }
        Ok(result)
    }

    /// Map local chart parameters to physical section coordinates.
    pub fn map_chart_point(&self, patch: usize, uv: [f64; 2]) -> Result<[f64; 2]> {
        let chart = *CHARTS.get(patch).ok_or("unknown L-section chart")?;
        if uv.iter().any(|x| !x.is_finite()) {
            return Err("non-finite chart parameter");
        }
        Ok([
            chart.origin[0] + chart.u[0] * uv[0] + chart.v[0] * uv[1],
            chart.origin[1] + chart.u[1] * uv[0] + chart.v[1] * uv[1],
        ])
    }

    /// Source levels that contribute nonzero basis functions on each integration leaf.
    pub fn active_levels_by_leaf(&self) -> Result<Vec<Vec<usize>>> {
        let integration = self.integration_leaves()?;
        let mut result = Vec::with_capacity(integration.len());
        for (cell, extracted) in self.leaves_data.iter().zip(integration) {
            let mut levels = BTreeSet::new();
            let dofs: BTreeSet<_> = extracted.dofs.into_iter().collect();
            for function in &self.patches[cell.patch].functions {
                if dofs.contains(&function.global)
                    && !function.by_level[cell.level].is_empty()
                {
                    levels.insert(function.source_level);
                }
            }
            result.push(levels.into_iter().collect());
        }
        Ok(result)
    }

    /// Verify the complete, nonoverlapping partition of each affine chart on the finest grid.
    pub fn audit_leaf_cover(&self) -> Result<bool> {
        let side = 8usize
            .checked_shl(self.maximum_level as u32)
            .ok_or("THB audit grid overflow")?;
        let mut occupancy = vec![vec![false; side * side]; 3];
        for leaf in &self.leaves_data {
            let mut bounds = [[0usize; 2]; 2];
            for (axis, axis_bounds) in bounds.iter_mut().enumerate() {
                for (edge, bound) in axis_bounds.iter_mut().enumerate() {
                    let value = leaf.interval[axis][edge];
                    let scaled = value
                        .numerator()
                        .checked_mul(side as i128)
                        .ok_or("THB grid index overflow")?;
                    if scaled % value.denominator() != 0 {
                        return Err("THB leaf endpoint is off the finest grid");
                    }
                    *bound = usize::try_from(scaled / value.denominator())
                        .map_err(|_| "THB grid index outside allocation")?;
                }
            }
            for i in bounds[0][0]..bounds[0][1] {
                for j in bounds[1][0]..bounds[1][1] {
                    let offset = i * side + j;
                    if occupancy[leaf.patch][offset] {
                        return Ok(false);
                    }
                    occupancy[leaf.patch][offset] = true;
                }
            }
        }
        Ok(occupancy.iter().all(|chart| chart.iter().all(|covered| *covered)))
    }

    /// Check the exact partition of unity in the common finest tensor basis.
    pub fn audit_exact_partition(&self) -> Result<bool> {
        let final_shape = self.patches[0].shapes[self.maximum_level];
        let count = final_shape[0] * final_shape[1];
        for patch in &self.patches {
            let mut sums = BTreeMap::<usize, Q>::new();
            for function in &patch.functions {
                let coefficients = function
                    .by_level
                    .last()
                    .ok_or("missing THB finest-level coefficients")?;
                for (&id, &coefficient) in coefficients {
                    rational_add(&mut sums, id, coefficient)?;
                }
            }
            if (0..count).any(|id| sums.get(&id).copied().unwrap_or(Q::zero()) != Q::one()) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Count active hierarchical functions by their source level, before seam identification.
    pub fn active_by_level(&self) -> Vec<usize> {
        let mut counts = vec![0; self.maximum_level + 1];
        for patch in &self.patches {
            for function in &patch.functions {
                counts[function.source_level] += 1;
            }
        }
        counts
    }

    /// Sample a patch basis and return physical section gradients.
    pub fn sample_patch(
        &self,
        patch_id: usize,
        level: usize,
        uv: [f64; 2],
    ) -> Result<Vec<(usize, f64, [f64; 2])>> {
        let patch = self.patches.get(patch_id).ok_or("unknown L-section chart")?;
        if level > self.maximum_level
            || uv.iter().any(|x| !x.is_finite() || !(0.0..=1.0).contains(x))
        {
            return Err("invalid THB sample coordinates");
        }
        let knots_u: Vec<_> = self.knots[level][0].iter().map(|&x| q_value(x)).collect();
        let knots_v: Vec<_> = self.knots[level][1].iter().map(|&x| q_value(x)).collect();
        let (nu, du) = bspline_values(&knots_u, DEGREE, uv[0]);
        let (nv, dv) = bspline_values(&knots_v, DEGREE, uv[1]);
        let chart = CHARTS[patch_id];
        let mut values = Vec::new();
        for function in &patch.functions {
            let coefficients = &function.by_level[level];
            if coefficients.is_empty() {
                continue;
            }
            let mut value = 0.0;
            let mut derivative_u = 0.0;
            let mut derivative_v = 0.0;
            for (&id, &coefficient) in coefficients {
                let i = id / patch.shapes[level][1];
                let j = id % patch.shapes[level][1];
                let c = q_value(coefficient);
                value += c * nu[i] * nv[j];
                derivative_u += c * du[i] * nv[j];
                derivative_v += c * nu[i] * dv[j];
            }
            if value != 0.0 || derivative_u != 0.0 || derivative_v != 0.0 {
                values.push((
                    function.global,
                    value,
                    [
                        chart.u[0] * derivative_u + chart.v[0] * derivative_v,
                        chart.u[1] * derivative_u + chart.v[1] * derivative_v,
                    ],
                ));
            }
        }
        Ok(values)
    }

    /// Exact coarse-space inclusion into the active THB basis.
    fn transfer_from(&self, coarse: &Self) -> Result<Sparse> {
        if self.is_uniform || coarse.is_uniform || self.maximum_level != coarse.maximum_level + 1 {
            return Err("THB transfer requires consecutive nested hierarchy levels");
        }
        if coarse.widths != self.widths[..coarse.maximum_level] {
            return Err("THB transfer box schedule differs on its common levels");
        }
        let added_level = self.maximum_level;
        let mut result = vec![Row::new(); self.active_scalar_dofs];
        for patch_id in 0..3 {
            let coarse_patch = &coarse.patches[patch_id];
            let fine_patch = &self.patches[patch_id];
            let fine_active: BTreeMap<_, _> = fine_patch
                .functions
                .iter()
                .map(|function| ((function.source_level, function.source_id), function.global))
                .collect();
            let top_level_active: BTreeMap<_, _> = fine_patch
                .functions
                .iter()
                .filter(|f| f.source_level == added_level)
                .map(|f| (f.source_id, f.global))
                .collect();
            for coarse_function in &coarse_patch.functions {
                let source = coarse_function.global;
                let mut coefficients = BTreeMap::<usize, Q>::new();
                if let Some(&same_level) = fine_active
                    .get(&(coarse_function.source_level, coarse_function.source_id))
                {
                    coefficients.insert(same_level, Q::one());
                }
                let expanded = tensor_apply(
                    &self.inclusions[coarse.maximum_level],
                    coarse_patch.shapes[coarse.maximum_level],
                    &coarse_function.by_level[coarse.maximum_level],
                )?;
                for (fine_id, coefficient) in expanded {
                    if !all_inside(
                        &self.knots[added_level],
                        fine_id,
                        region_end(added_level, &self.widths),
                    ) {
                        continue;
                    }
                    let global = top_level_active
                        .get(&fine_id)
                        .copied()
                        .ok_or("missing active THB function in the added level")?;
                    let entry = coefficients.entry(global).or_insert(Q::zero());
                    *entry = entry.checked_add(coefficient)?;
                }
                if let Some(&same_level) = fine_patch
                    .functions
                    .iter()
                    .find(|f| (f.source_level, f.source_id) == (coarse_function.source_level, coarse_function.source_id))
                    .map(|f| &f.global)
                {
                    coefficients.entry(same_level).or_insert(Q::one());
                }
                for (global, coefficient) in coefficients {
                    if let Some(previous) = result[global].insert(source, coefficient) {
                        if previous != coefficient {
                            return Err("exact THB transfer disagrees across a seam");
                        }
                    }
                }
            }
        }
        Ok(result)
    }

    fn transfer_from_prepared(&self) -> Result<Sparse> {
        let base = build_hierarchy(0)?;
        let mut result = vec![Row::new(); self.active_scalar_dofs];
        for source in 0..base.active_scalar_dofs {
            for patch_id in 0..3 {
                let base_patch = &base.patches[patch_id];
                let mut coefficients = BTreeMap::new();
                for function in &base_patch.functions {
                    if function.global == source {
                        coefficients.insert(function.source_id, Q::one());
                    }
                }
                if coefficients.is_empty() {
                    continue;
                }
                let patch = &self.patches[patch_id];
                let mut hierarchical = BTreeMap::<usize, Q>::new();
                for level in 0..=self.maximum_level {
                    let active: BTreeMap<_, _> = patch
                        .functions
                        .iter()
                        .filter(|f| f.source_level == level)
                        .map(|f| (f.source_id, f.global))
                        .collect();
                    if level == self.maximum_level {
                        for (&id, &coefficient) in &coefficients {
                            let global = active
                                .get(&id)
                                .copied()
                                .ok_or("coarse transfer reached an inactive finest function")?;
                            hierarchical.insert(global, coefficient);
                        }
                        break;
                    }
                    for (&id, &coefficient) in &coefficients {
                        if let Some(&global) = active.get(&id) {
                            hierarchical.insert(global, coefficient);
                        }
                    }
                    let fine_shape = patch.shapes[level + 1];
                    let mut residual = BTreeMap::new();
                    for i in 0..fine_shape[0] {
                        for j in 0..fine_shape[1] {
                            let fine_id = index(fine_shape, i, j);
                            if !all_inside(
                                &self.knots[level + 1],
                                fine_id,
                                region_end(level + 1, &self.widths),
                            ) {
                                continue;
                            }
                            let mut value = Q::zero();
                            for (&old_i, &a) in &self.inclusions[level][0][i] {
                                for (&old_j, &b) in &self.inclusions[level][1][j] {
                                    if let Some(&old) = coefficients.get(&index(patch.shapes[level], old_i, old_j)) {
                                        value = value.checked_add(a.checked_mul(b)?.checked_mul(old)?)?;
                                    }
                                }
                            }
                            if value != Q::zero() {
                                residual.insert(fine_id, value);
                            }
                        }
                    }
                    coefficients = residual;
                }
                for (global, coefficient) in hierarchical {
                    if let Some(previous) = result[global].insert(source, coefficient) {
                        if previous != coefficient {
                            return Err("prepared THB transfer disagrees across a seam");
                        }
                    }
                }
            }
        }
        Ok(result)
    }

    /// Check exact representation of every prepared spline in the THB space.
    pub fn audit_exact_transfer(&self) -> Result<bool> {
        if self.is_uniform {
            if self.maximum_level == 0 {
                return Ok(true);
            }
            return Err("exact transfer audit is defined for THB spaces, not uniform controls");
        }
        if self.maximum_level == 0 {
            return Ok(true);
        }
        let base = build_hierarchy(0)?;
        let transfer = self.transfer_from_prepared()?;
        if sparse_rank(&transfer)? != base.active_scalar_dofs {
            return Ok(false);
        }
        for patch_id in 0..3 {
            let base_patch = &base.patches[patch_id];
            let fine_patch = &self.patches[patch_id];
            for base_function in &base_patch.functions {
                let column = base_function.global;
                let mut target = BTreeMap::from([(base_function.source_id, Q::one())]);
                for level in 0..self.maximum_level {
                    target = tensor_apply(
                        &self.inclusions[level],
                        fine_patch.shapes[level],
                        &target,
                    )?;
                }
                let mut reconstructed = BTreeMap::<usize, Q>::new();
                for function in &fine_patch.functions {
                    let coefficient = transfer[function.global]
                        .get(&column)
                        .copied()
                        .unwrap_or(Q::zero());
                    for (&id, &factor) in &function.by_level[self.maximum_level] {
                        if coefficient != Q::zero() && factor != Q::zero() {
                            rational_add(&mut reconstructed, id, coefficient.checked_mul(factor)?)?;
                        }
                    }
                }
                if reconstructed != target {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    /// Check exact nesting and rank for every consecutive THB level pair.
    pub fn audit_nested_transfers(&self) -> Result<bool> {
        if self.is_uniform {
            if self.maximum_level == 0 {
                return Ok(true);
            }
            return Err("nested transfer audit is defined for THB spaces, not uniform controls");
        }
        if self.maximum_level == 0 {
            return Ok(true);
        }
        let mut coarse = build_hierarchy(0)?;
        for level in 1..=self.maximum_level {
            let fine = build_hierarchy_with_widths(level, &self.widths[..level])?;
            if !audit_transfer_pair(&coarse, &fine)? {
                return Ok(false);
            }
            coarse = fine;
        }
        Ok(true)
    }
}

fn audit_transfer_pair(coarse: &ReentrantHierarchy, fine: &ReentrantHierarchy) -> Result<bool> {
        let transfer = fine.transfer_from(coarse)?;
        let base_ndof = coarse.active_scalar_dofs;
        let mut rank_rows = Vec::<Row>::new();
        for row in &transfer {
            if !row.is_empty() {
                rank_rows.push(row.clone());
            }
        }
        if sparse_rank(&rank_rows)? != base_ndof {
            return Ok(false);
        }
        for patch_id in 0..3 {
            let coarse_patch = &coarse.patches[patch_id];
            let fine_patch = &fine.patches[patch_id];
            for coarse_function in &coarse_patch.functions {
                let column = coarse_function.global;
                let mut target = coarse_function.by_level[coarse.maximum_level].clone();
                for level in coarse.maximum_level..fine.maximum_level {
                    target = tensor_apply(
                        &fine.inclusions[level],
                        fine_patch.shapes[level],
                        &target,
                    )?;
                }
                let mut reconstructed = BTreeMap::<usize, Q>::new();
                for function in &fine_patch.functions {
                    let coefficient = transfer[function.global]
                        .get(&column)
                        .copied()
                        .unwrap_or(Q::zero());
                    if coefficient == Q::zero() {
                        continue;
                    }
                    for (&id, &factor) in &function.by_level[fine.maximum_level] {
                        rational_add(&mut reconstructed, id, coefficient.checked_mul(factor)?)?;
                    }
                }
                if reconstructed != target {
                    return Ok(false);
                }
            }
        }
        Ok(true)
}

fn tensor_apply(
    axes: &[Sparse; 2],
    old_shape: [usize; 2],
    coefficients: &BTreeMap<usize, Q>,
) -> Result<BTreeMap<usize, Q>> {
    let new_shape = [axes[0].len(), axes[1].len()];
    let mut result = BTreeMap::new();
    for i in 0..new_shape[0] {
        for j in 0..new_shape[1] {
            let mut value = Q::zero();
            for (&oi, &a) in &axes[0][i] {
                for (&oj, &b) in &axes[1][j] {
                    if let Some(&c) = coefficients.get(&index(old_shape, oi, oj)) {
                        value = value.checked_add(a.checked_mul(b)?.checked_mul(c)?)?;
                    }
                }
            }
            if value != Q::zero() {
                result.insert(index(new_shape, i, j), value);
            }
        }
    }
    Ok(result)
}

fn sparse_rank(rows: &[Row]) -> Result<usize> {
    let mut pivots = BTreeMap::<usize, Row>::new();
    for original in rows {
        let mut row = original.clone();
        while let Some((&column, &coefficient)) = row.first_key_value() {
            if let Some(pivot) = pivots.get(&column) {
                let factor = Q::zero().checked_sub(coefficient)?;
                for (&j, &value) in pivot {
                    add(&mut row, j, factor.checked_mul(value)?)?;
                }
            } else {
                let scale = Q::new(coefficient.denominator(), coefficient.numerator())?;
                for value in row.values_mut() {
                    *value = value.checked_mul(scale)?;
                }
                pivots.insert(column, row);
                break;
            }
        }
    }
    Ok(pivots.len())
}

fn q_value(value: Q) -> f64 {
    value.numerator() as f64 / value.denominator() as f64
}

fn bspline_values(knots: &[f64], degree: usize, x: f64) -> (Vec<f64>, Vec<f64>) {
    let mut current: Vec<f64> = knots
        .windows(2)
        .map(|w| f64::from(w[0] <= x && x < w[1]))
        .collect();
    for d in 1..=degree {
        let mut next = vec![0.0; current.len() - 1];
        for i in 0..next.len() {
            let left = knots[i + d] - knots[i];
            let right = knots[i + d + 1] - knots[i + 1];
            if left > 0.0 {
                next[i] += (x - knots[i]) / left * current[i];
            }
            if right > 0.0 {
                next[i] += (knots[i + d + 1] - x) / right * current[i + 1];
            }
        }
        current = next;
    }
    let mut lower = knots
        .windows(2)
        .map(|w| f64::from(w[0] <= x && x < w[1]))
        .collect::<Vec<_>>();
    for d in 1..degree {
        let mut next = vec![0.0; lower.len() - 1];
        for i in 0..next.len() {
            let left = knots[i + d] - knots[i];
            let right = knots[i + d + 1] - knots[i + 1];
            if left > 0.0 {
                next[i] += (x - knots[i]) / left * lower[i];
            }
            if right > 0.0 {
                next[i] += (knots[i + d + 1] - x) / right * lower[i + 1];
            }
        }
        lower = next;
    }
    let n = knots.len() - degree - 1;
    let mut derivative = vec![0.0; n];
    for i in 0..n {
        let left = knots[i + degree] - knots[i];
        let right = knots[i + degree + 1] - knots[i + 1];
        if left > 0.0 {
            derivative[i] += degree as f64 / left * lower[i];
        }
        if right > 0.0 {
            derivative[i] -= degree as f64 / right * lower[i + 1];
        }
    }
    (current, derivative)
}
