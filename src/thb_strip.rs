//! Exact one-level truncated hierarchical B-splines on the recovered strip.
//!
//! The pilot keeps the geometry parent separate from the field basis. It uses
//! a uniformly prepared tricubic tensor space, then one strictly interior
//! dyadic refinement region. It is a THB basis, not a T-spline construction.

use crate::{
    constraint_system::{Lift, LiftLine, Ratio as Q},
    extracted_space::{ElementBasisSample, ElementGeometrySample},
    inherited_space::val,
    json::{object, Json},
    parent_field::axis_rows,
    spline_constraints::{add, Row},
    spline_support::{recover_nurbs_strip_parent, TensorSupport},
    RationalBezierMap,
};
use std::collections::{BTreeMap, BTreeSet};

type Result<T> = std::result::Result<T, &'static str>;
type Sparse = Vec<Row>;
type V = [f64; 3];

fn index(shape: [usize; 3], [i, j, k]: [usize; 3]) -> usize {
    (i * shape[1] + j) * shape[2] + k
}

fn tensor_shape(knots: [&[Q]; 3], degrees: [usize; 3]) -> [usize; 3] {
    std::array::from_fn(|axis| knots[axis].len() - degrees[axis] - 1)
}

fn q(value: i64) -> Q {
    Q::new(i128::from(value), 1).expect("integer rational")
}

fn add_scaled(into: &mut Row, row: &Row, scale: Q) -> Result<()> {
    for (&index, &value) in row {
        add(into, index, value.checked_mul(scale)?)?;
    }
    Ok(())
}

fn unique_intervals(knots: &[Q]) -> Vec<[Q; 2]> {
    knots
        .windows(2)
        .filter_map(|pair| {
            pair[0]
                .cmp_value(&pair[1])
                .is_lt()
                .then_some([pair[0], pair[1]])
        })
        .collect()
}

fn midpoint(a: Q, b: Q) -> Result<Q> {
    a.checked_add(b)?.checked_div(Q::new(2, 1)?)
}

fn refine_knots_once(knots: &[Q]) -> Result<Vec<Q>> {
    let mut result = knots.to_vec();
    for [lo, hi] in unique_intervals(knots) {
        result.push(midpoint(lo, hi)?);
    }
    result.sort_by(Q::cmp_value);
    Ok(result)
}

fn elevation_matrix(from: usize, to: usize) -> Result<Vec<Vec<Q>>> {
    if from > to || to > 3 {
        return Err("unsupported Bernstein elevation");
    }
    fn choose(n: usize, k: usize) -> i128 {
        if k > n {
            return 0;
        }
        let k = k.min(n - k);
        (1..=k).fold(1_i128, |value, i| value * (n - k + i) as i128 / i as i128)
    }
    let mut matrix = vec![vec![Q::zero(); from + 1]; to + 1];
    for (i, row) in matrix.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            let other = i.saturating_sub(j);
            if j <= i && other <= to - from {
                *value = Q::new(choose(from, j) * choose(to - from, other), choose(to, i))?;
            }
        }
    }
    Ok(matrix)
}

fn inverse(matrix: &[Vec<Q>]) -> Result<Vec<Vec<Q>>> {
    let n = matrix.len();
    if n == 0 || matrix.iter().any(|row| row.len() != n) {
        return Err("local spline extraction is not square");
    }
    let mut augmented = vec![vec![Q::zero(); n * 2]; n];
    for i in 0..n {
        augmented[i][..n].copy_from_slice(&matrix[i]);
        augmented[i][n + i] = Q::one();
    }
    for col in 0..n {
        let pivot = (col..n)
            .find(|&row| augmented[row][col] != Q::zero())
            .ok_or("singular local spline extraction")?;
        augmented.swap(col, pivot);
        let inverse_pivot = Q::new(
            augmented[col][col].denominator(),
            augmented[col][col].numerator(),
        )?;
        for value in &mut augmented[col] {
            *value = value.checked_mul(inverse_pivot)?;
        }
        for row in 0..n {
            if row == col || augmented[row][col] == Q::zero() {
                continue;
            }
            let factor = Q::zero().checked_sub(augmented[row][col])?;
            let pivot_row = augmented[col].clone();
            for (target, source) in augmented[row].iter_mut().zip(&pivot_row) {
                *target = target.checked_add(factor.checked_mul(*source)?)?;
            }
        }
    }
    Ok(augmented.into_iter().map(|row| row[n..].to_vec()).collect())
}

/// Exact inclusion matrix: rows are new spline functions, columns are old.
fn axis_inclusion(
    old_knots: &[Q],
    old_degree: usize,
    new_knots: &[Q],
    new_degree: usize,
) -> Result<Sparse> {
    if old_degree > new_degree || old_degree > 3 || new_degree > 3 {
        return Err("unsupported spline degree inclusion");
    }
    let old_count = old_knots.len() - old_degree - 1;
    let new_count = new_knots.len() - new_degree - 1;
    let mut rows = vec![Row::new(); new_count];
    let elevate = elevation_matrix(old_degree, new_degree)?;
    for interval in unique_intervals(new_knots) {
        let old = axis_rows(old_knots, old_degree, interval)?;
        let new = axis_rows(new_knots, new_degree, interval)?;
        let active: Vec<_> = new
            .iter()
            .flat_map(|row| row.keys().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        if active.len() != new_degree + 1 {
            return Err("local fine spline basis has wrong support count");
        }
        let matrix: Vec<Vec<Q>> = new
            .iter()
            .map(|row| {
                active
                    .iter()
                    .map(|index| row.get(index).copied().unwrap_or(Q::zero()))
                    .collect()
            })
            .collect();
        let inv = inverse(&matrix)?;
        let mut rhs = vec![Row::new(); new_degree + 1];
        for bernstein in 0..=new_degree {
            for old_bernstein in 0..=old_degree {
                add_scaled(
                    &mut rhs[bernstein],
                    &old[old_bernstein],
                    elevate[bernstein][old_bernstein],
                )?;
            }
        }
        for (local_fine, &fine_id) in active.iter().enumerate() {
            let mut transfer = Row::new();
            for old_bernstein in 0..=new_degree {
                add_scaled(
                    &mut transfer,
                    &rhs[old_bernstein],
                    inv[local_fine][old_bernstein],
                )?;
            }
            for (&old_id, &coefficient) in &transfer {
                if coefficient == Q::zero() {
                    continue;
                }
                if let Some(existing) = rows[fine_id].get(&old_id) {
                    if *existing != coefficient {
                        return Err("inconsistent exact spline inclusion across spans");
                    }
                } else {
                    rows[fine_id].insert(old_id, coefficient);
                }
            }
        }
    }
    if rows.len() != new_count || rows.iter().any(|row| row.keys().any(|&i| i >= old_count)) {
        return Err("invalid spline inclusion matrix");
    }
    for row in &rows {
        let sum = row
            .values()
            .try_fold(Q::zero(), |sum, &value| sum.checked_add(value))?;
        if sum != Q::one() {
            return Err("spline inclusion row does not preserve constants");
        }
    }
    Ok(rows)
}

fn tensor_inclusion(axes: [&Sparse; 3], old_counts: [usize; 3]) -> Result<Sparse> {
    let fine_shape = [axes[0].len(), axes[1].len(), axes[2].len()];
    let total = fine_shape.iter().product();
    let mut rows = vec![Row::new(); total];
    for i in 0..fine_shape[0] {
        for j in 0..fine_shape[1] {
            for k in 0..fine_shape[2] {
                let row = &mut rows[index(fine_shape, [i, j, k])];
                for (&a, &ca) in &axes[0][i] {
                    for (&b, &cb) in &axes[1][j] {
                        for (&c, &cc) in &axes[2][k] {
                            let old_id = index(old_counts, [a, b, c]);
                            add(row, old_id, ca.checked_mul(cb)?.checked_mul(cc)?)?;
                        }
                    }
                }
            }
        }
    }
    Ok(rows)
}

fn compose(a: &[Row], b: &[Row]) -> Result<Sparse> {
    let mut result = vec![Row::new(); a.len()];
    for (row_id, row) in a.iter().enumerate() {
        for (&middle, &coefficient) in row {
            let source = b
                .get(middle)
                .ok_or("spline transfer composition index out of range")?;
            add_scaled(&mut result[row_id], source, coefficient)?;
        }
    }
    Ok(result)
}

fn apply(rows: &[Row], values: &[Q]) -> Result<Vec<Q>> {
    rows.iter()
        .map(|row| {
            row.iter().try_fold(Q::zero(), |sum, (&i, &coefficient)| {
                sum.checked_add(
                    coefficient.checked_mul(
                        *values
                            .get(i)
                            .ok_or("spline coefficient index out of range")?,
                    )?,
                )
            })
        })
        .collect()
}

fn apply_controls(rows: &[Row], controls: &[[Q; 4]]) -> Result<Vec<[Q; 4]>> {
    let mut result = Vec::with_capacity(rows.len());
    for row in rows {
        let mut control = [Q::zero(); 4];
        for (&id, &coefficient) in row {
            let source = controls
                .get(id)
                .ok_or("homogeneous transfer control index out of range")?;
            for c in 0..4 {
                control[c] = control[c].checked_add(coefficient.checked_mul(source[c])?)?;
            }
        }
        result.push(control);
    }
    Ok(result)
}

fn exact_sparse_rank(rows: &[Row]) -> Result<usize> {
    let mut pivots: BTreeMap<usize, Row> = BTreeMap::new();
    for source in rows {
        let mut row = source.clone();
        while let Some((&column, &leading)) = row.first_key_value() {
            if let Some(pivot) = pivots.get(&column) {
                let factor = leading;
                for (&j, &value) in pivot {
                    let updated = row
                        .get(&j)
                        .copied()
                        .unwrap_or(Q::zero())
                        .checked_sub(factor.checked_mul(value)?)?;
                    if updated == Q::zero() {
                        row.remove(&j);
                    } else {
                        row.insert(j, updated);
                    }
                }
            } else {
                let inverse = Q::new(leading.denominator(), leading.numerator())?;
                for value in row.values_mut() {
                    *value = value.checked_mul(inverse)?;
                }
                pivots.insert(column, row);
                break;
            }
        }
    }
    Ok(pivots.len())
}

fn choose(n: usize, k: usize) -> i128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    (1..=k).fold(1_i128, |value, i| value * (n - k + i) as i128 / i as i128)
}

fn bernstein_product(a: &[Q], b: &[Q]) -> Result<Vec<Q>> {
    let (m, n) = (a.len() - 1, b.len() - 1);
    let mut result = vec![Q::zero(); m + n + 1];
    for (i, &a_value) in a.iter().enumerate() {
        for (j, &b_value) in b.iter().enumerate() {
            let k = i + j;
            let factor = Q::new(choose(m, i) * choose(n, j), choose(m + n, k))?;
            result[k] =
                result[k].checked_add(a_value.checked_mul(b_value)?.checked_mul(factor)?)?;
        }
    }
    Ok(result)
}

fn split_bernstein(coefficients: &[Q]) -> Result<(Vec<Q>, Vec<Q>)> {
    let n = coefficients.len() - 1;
    let half = Q::new(1, 2)?;
    let mut work = coefficients.to_vec();
    let mut left = vec![Q::zero(); n + 1];
    let mut right = vec![Q::zero(); n + 1];
    left[0] = work[0];
    right[n] = work[n];
    for level in 1..=n {
        for i in 0..=n - level {
            work[i] = work[i].checked_add(work[i + 1])?.checked_mul(half)?;
        }
        left[level] = work[0];
        right[n - level] = work[n - level];
    }
    Ok((left, right))
}

fn certify_positive_bernstein(coefficients: &[Q], depth: usize) -> Result<bool> {
    if coefficients.iter().all(|x| *x > Q::zero()) {
        return Ok(true);
    }
    if coefficients.iter().all(|x| *x < Q::zero()) || depth == 0 {
        return Ok(false);
    }
    let (left, right) = split_bernstein(coefficients)?;
    Ok(certify_positive_bernstein(&left, depth - 1)?
        && certify_positive_bernstein(&right, depth - 1)?)
}

/// Certify the signed Jacobian of the source extrusion on every original U span.
fn certify_parent_orientation(support: &TensorSupport) -> Result<()> {
    let knots: Vec<Q> = support.knots[0].iter().copied().map(q).collect();
    let shape = support.control_shape;
    let source_spans = unique_intervals(&knots);
    for interval in source_spans {
        let extraction = axis_rows(&knots, 3, interval)?;
        let mut x = [Q::zero(); 4];
        let mut lower = [Q::zero(); 4];
        let mut upper = [Q::zero(); 4];
        let mut weight = [Q::zero(); 4];
        for (r, row) in extraction.iter().enumerate() {
            for (&u, &coefficient) in row {
                let low = support.homogeneous_controls[index(shape, [u, 0, 0])];
                let high = support.homogeneous_controls[index(shape, [u, 1, 0])];
                x[r] = x[r].checked_add(coefficient.checked_mul(Q::new(low[0], 1)?)?)?;
                lower[r] = lower[r].checked_add(coefficient.checked_mul(Q::new(low[1], 1)?)?)?;
                upper[r] = upper[r].checked_add(coefficient.checked_mul(Q::new(high[1], 1)?)?)?;
                weight[r] = weight[r].checked_add(coefficient.checked_mul(Q::new(low[3], 1)?)?)?;
            }
        }
        if weight.iter().any(|w| *w <= Q::zero()) {
            return Err("source chart has a nonpositive rational denominator");
        }
        let dx: Vec<_> = (0..3)
            .map(|i| x[i + 1].checked_sub(x[i]).and_then(|v| v.checked_mul(q(3))))
            .collect::<Result<_>>()?;
        let dw: Vec<_> = (0..3)
            .map(|i| {
                weight[i + 1]
                    .checked_sub(weight[i])
                    .and_then(|v| v.checked_mul(q(3)))
            })
            .collect::<Result<_>>()?;
        let mut numerator = bernstein_product(&dx, &weight)?;
        let second = bernstein_product(&x, &dw)?;
        for (a, b) in numerator.iter_mut().zip(second) {
            *a = a.checked_sub(b)?;
        }
        let gap: Vec<_> = upper
            .iter()
            .zip(&lower)
            .map(|(&a, &b)| a.checked_sub(b))
            .collect::<Result<_>>()?;
        let determinant_numerator = bernstein_product(&numerator, &gap)?;
        if !certify_positive_bernstein(&determinant_numerator, 10)? {
            return Err("source strip orientation is reversed, singular, or uncertified");
        }
    }
    Ok(())
}

fn parent_knots(support: &TensorSupport) -> [Vec<Q>; 3] {
    std::array::from_fn(|axis| support.knots[axis].iter().copied().map(q).collect())
}

fn cubic_axes(support: &TensorSupport) -> [Vec<Q>; 3] {
    let source = parent_knots(support);
    [
        source[0].clone(),
        vec![
            Q::zero(),
            Q::zero(),
            Q::zero(),
            Q::zero(),
            Q::one(),
            Q::one(),
            Q::one(),
            Q::one(),
        ],
        vec![
            Q::zero(),
            Q::zero(),
            Q::zero(),
            Q::zero(),
            Q::one(),
            Q::one(),
            Q::one(),
            Q::one(),
        ],
    ]
}

fn uniform_levels(initial: &[Vec<Q>; 3], count: usize) -> Result<Vec<[Vec<Q>; 3]>> {
    let mut levels = vec![initial.clone()];
    for _ in 0..count {
        let previous = levels.last().ok_or("missing initial knot level")?;
        let next = [
            refine_knots_once(&previous[0])?,
            refine_knots_once(&previous[1])?,
            refine_knots_once(&previous[2])?,
        ];
        levels.push(next);
    }
    Ok(levels)
}

fn parent_to_initial(support: &TensorSupport, initial_knots: &[Vec<Q>; 3]) -> Result<Sparse> {
    let parent_knots = parent_knots(support);
    let parent_degrees = support.degrees;
    let initial_degrees = [3; 3];
    let axes: [Result<Sparse>; 3] = std::array::from_fn(|axis| {
        axis_inclusion(
            &parent_knots[axis],
            parent_degrees[axis],
            &initial_knots[axis],
            initial_degrees[axis],
        )
    });
    let old_counts =
        std::array::from_fn(|axis| parent_knots[axis].len() - parent_degrees[axis] - 1);
    let [a, b, c] = axes;
    tensor_inclusion([&a?, &b?, &c?], old_counts)
}

fn stage_inclusion(old_knots: &[Vec<Q>; 3], new_knots: &[Vec<Q>; 3]) -> Result<Sparse> {
    let axes: [Result<Sparse>; 3] =
        std::array::from_fn(|axis| axis_inclusion(&old_knots[axis], 3, &new_knots[axis], 3));
    let old_counts = std::array::from_fn(|axis| old_knots[axis].len() - 4);
    let [a, b, c] = axes;
    tensor_inclusion([&a?, &b?, &c?], old_counts)
}

fn box_contains(region: &[[Q; 2]; 3], cell: &[[Q; 2]; 3]) -> bool {
    (0..3).all(|axis| {
        !cell[axis][0].cmp_value(&region[axis][0]).is_lt()
            && !cell[axis][1].cmp_value(&region[axis][1]).is_gt()
    })
}

fn support_inside(
    knots: &[Vec<Q>; 3],
    degrees: [usize; 3],
    shape: [usize; 3],
    id: usize,
    region: &[[Q; 2]; 3],
) -> bool {
    let mut remainder = id;
    let k = remainder % shape[2];
    remainder /= shape[2];
    let j = remainder % shape[1];
    let i = remainder / shape[1];
    [i, j, k].iter().enumerate().all(|(axis, &basis)| {
        let lo = knots[axis][basis];
        let hi = knots[axis][basis + degrees[axis] + 1];
        !lo.cmp_value(&region[axis][0]).is_lt() && !hi.cmp_value(&region[axis][1]).is_gt()
    })
}

fn tensor_rows(axes: [&[Row]; 3], basis_shape: [usize; 3]) -> Result<Vec<Row>> {
    let bernstein_shape = [axes[0].len(), axes[1].len(), axes[2].len()];
    let mut rows = vec![Row::new(); bernstein_shape.iter().product()];
    for i in 0..bernstein_shape[0] {
        for j in 0..bernstein_shape[1] {
            for k in 0..bernstein_shape[2] {
                let row = &mut rows[index(bernstein_shape, [i, j, k])];
                for (&a, &ca) in &axes[0][i] {
                    for (&b, &cb) in &axes[1][j] {
                        for (&c, &cc) in &axes[2][k] {
                            add(
                                row,
                                index(basis_shape, [a, b, c]),
                                ca.checked_mul(cb)?.checked_mul(cc)?,
                            )?;
                        }
                    }
                }
            }
        }
    }
    Ok(rows)
}

fn exact_cell_rows(knots: &[Vec<Q>; 3], cell: &[[Q; 2]; 3], shape: [usize; 3]) -> Result<Vec<Row>> {
    let axes = [
        axis_rows(&knots[0], 3, cell[0])?,
        axis_rows(&knots[1], 3, cell[1])?,
        axis_rows(&knots[2], 3, cell[2])?,
    ];
    tensor_rows([&axes[0], &axes[1], &axes[2]], shape)
}

fn compact_rows(
    rows: &[Row],
    global_to_dof: &BTreeMap<usize, usize>,
) -> Result<(Vec<usize>, Vec<Row>)> {
    let ids: Vec<_> = rows
        .iter()
        .flat_map(|row| row.keys().copied())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let dofs = ids
        .iter()
        .map(|id| {
            global_to_dof
                .get(id)
                .copied()
                .ok_or("local extraction references an inactive basis function")
        })
        .collect::<Result<Vec<_>>>()?;
    let local_columns: BTreeMap<_, _> = ids
        .iter()
        .enumerate()
        .map(|(local, &id)| (id, local))
        .collect();
    let compact: Vec<_> = rows
        .iter()
        .map(|row| {
            let mut result = Row::new();
            for (&id, &coefficient) in row {
                if coefficient < Q::zero() {
                    return Err("THB extraction has a negative Bernstein coefficient");
                }
                result.insert(local_columns[&id], coefficient);
            }
            Ok(result)
        })
        .collect::<Result<Vec<_>>>()?;
    Ok((dofs, compact))
}

#[allow(clippy::too_many_arguments)]
fn make_thb_rows(
    cell: &[[Q; 2]; 3],
    in_region: bool,
    prepared_knots: &[Vec<Q>; 3],
    fine_knots: &[Vec<Q>; 3],
    prepared_shape: [usize; 3],
    fine_shape: [usize; 3],
    coarse_active: &[usize],
    fine_active: &[usize],
    fine_to_prepared: &[Row],
    coarse_to_dof: &BTreeMap<usize, usize>,
    fine_to_dof: &BTreeMap<usize, usize>,
) -> Result<(Vec<usize>, Vec<Row>)> {
    let coarse = exact_cell_rows(prepared_knots, cell, prepared_shape)?;
    let fine = if in_region {
        exact_cell_rows(fine_knots, cell, fine_shape)?
    } else {
        vec![Row::new(); 64]
    };
    let mut rows = vec![Row::new(); 64];
    for (bernstein_row, target) in rows.iter_mut().enumerate() {
        for &id in coarse_active {
            if let Some(&value) = coarse[bernstein_row].get(&id) {
                target.insert(id, value);
            }
        }
        if in_region {
            for &fine_id in fine_active {
                let coefficient = fine[bernstein_row]
                    .get(&fine_id)
                    .copied()
                    .unwrap_or(Q::zero());
                if coefficient == Q::zero() {
                    continue;
                }
                for &coarse_id in coarse_active {
                    let p = fine_to_prepared[fine_id]
                        .get(&coarse_id)
                        .copied()
                        .unwrap_or(Q::zero());
                    if p != Q::zero() {
                        let old = target.get(&coarse_id).copied().unwrap_or(Q::zero());
                        let updated = old.checked_sub(p.checked_mul(coefficient)?)?;
                        if updated == Q::zero() {
                            target.remove(&coarse_id);
                        } else {
                            target.insert(coarse_id, updated);
                        }
                    }
                }
                target.insert(
                    fine_id + prepared_shape.iter().product::<usize>(),
                    coefficient,
                );
            }
        }
        if target
            .values()
            .try_fold(Q::zero(), |sum, &x| sum.checked_add(x))?
            != Q::one()
        {
            return Err("THB polynomial extraction does not preserve partition of unity");
        }
    }
    let mut global_map = coarse_to_dof.clone();
    for (&fine_id, &dof) in fine_to_dof {
        global_map.insert(fine_id + prepared_shape.iter().product::<usize>(), dof);
    }
    compact_rows(&rows, &global_map)
}

fn rational_transfer(polynomial: &[Row], old_weights: &[Q], new_weights: &[Q]) -> Result<Sparse> {
    if polynomial
        .iter()
        .any(|row| row.keys().any(|&i| i >= old_weights.len()))
        || polynomial.len() != new_weights.len()
    {
        return Err("rational spline transfer dimensions do not match");
    }
    let result: Sparse = polynomial
        .iter()
        .zip(new_weights)
        .map(|(row, &weight)| {
            row.iter()
                .try_fold(Row::new(), |mut result, (&i, &coefficient)| {
                    result.insert(
                        i,
                        coefficient
                            .checked_mul(old_weights[i])?
                            .checked_div(weight)?,
                    );
                    Ok(result)
                })
        })
        .collect::<Result<_>>()?;
    if result.iter().any(|row| {
        row.values()
            .try_fold(Q::zero(), |sum, &value| sum.checked_add(value))
            .ok()
            != Some(Q::one())
    }) {
        return Err("rational transfer does not preserve constant fields");
    }
    Ok(result)
}

/// Construct the fixed one-level interior-refinement pilot from an audited strip model.
pub fn build_thb_strip(model: &Json) -> Result<ThbStrip> {
    build_thb_strip_for_cells(model, [2, 1, 1], [6, 3, 3])
}

/// Build the strip model from its retained NURBS sketch recipe and refine it.
pub fn build_thb_strip_from_recipe(recipe: &Json) -> Result<ThbStrip> {
    let model = crate::nurbs_sketch_json::build(recipe)?;
    build_thb_strip(&model)
}

/// Build the strip model from a recipe with an explicit half-open prepared-cell mark.
pub fn build_thb_strip_from_recipe_cells(
    recipe: &Json,
    start: [usize; 3],
    end: [usize; 3],
) -> Result<ThbStrip> {
    let model = crate::nurbs_sketch_json::build(recipe)?;
    build_thb_strip_for_cells(&model, start, end)
}

/// Construct a one-level refinement marked by the half-open prepared-cell box
/// `start..end`. The box must be strictly interior and expose at least one
/// fine basis function.
pub fn build_thb_strip_for_cells(
    model: &Json,
    start: [usize; 3],
    end: [usize; 3],
) -> Result<ThbStrip> {
    let support = recover_nurbs_strip_parent(model)?;
    crate::spline_support::initial_mesh_admission(&support)?;
    if support.degrees[0] != 3 {
        return Err("THB strip pilot currently requires a cubic source curve");
    }
    let source_span_count = unique_intervals(&parent_knots(&support)[0]).len();
    let prepared_u = support.control_shape[0]
        .checked_add(
            3usize
                .checked_mul(source_span_count)
                .ok_or("THB size estimate overflow")?,
        )
        .ok_or("THB size estimate overflow")?;
    let prepared_estimate = prepared_u
        .checked_mul(49)
        .ok_or("THB size estimate overflow")?;
    let fine_estimate = prepared_u
        .checked_add(
            4usize
                .checked_mul(source_span_count)
                .ok_or("THB size estimate overflow")?,
        )
        .ok_or("THB size estimate overflow")?
        .checked_mul(121)
        .ok_or("THB size estimate overflow")?;
    let leaf_estimate = source_span_count
        .checked_mul(512)
        .ok_or("THB size estimate overflow")?;
    if source_span_count == 0
        || source_span_count > 8
        || prepared_estimate > 20_000
        || fine_estimate > 100_000
        || leaf_estimate > 4096
    {
        return Err("source strip exceeds the checked THB pilot allocation budget");
    }
    certify_parent_orientation(&support)?;
    let initial_knots = cubic_axes(&support);
    let levels = uniform_levels(&initial_knots, 2)?;
    let prepared_knots = &levels[2];
    let initial_shape = tensor_shape(
        [&initial_knots[0], &initial_knots[1], &initial_knots[2]],
        [3; 3],
    );
    let prepared_shape = tensor_shape(
        [&prepared_knots[0], &prepared_knots[1], &prepared_knots[2]],
        [3; 3],
    );
    let parent_transfer = parent_to_initial(&support, &initial_knots)?;
    let p01 = stage_inclusion(&levels[0], &levels[1])?;
    let p12 = stage_inclusion(&levels[1], &levels[2])?;
    let initial_to_prepared = compose(&p12, &p01)?;
    let fine_knots = [
        refine_knots_once(&prepared_knots[0])?,
        refine_knots_once(&prepared_knots[1])?,
        refine_knots_once(&prepared_knots[2])?,
    ];
    let fine_shape = tensor_shape([&fine_knots[0], &fine_knots[1], &fine_knots[2]], [3; 3]);
    let prepared_to_fine = stage_inclusion(prepared_knots, &fine_knots)?;
    let spans = [
        unique_intervals(&prepared_knots[0]),
        unique_intervals(&prepared_knots[1]),
        unique_intervals(&prepared_knots[2]),
    ];
    if (0..3)
        .any(|axis| start[axis] == 0 || start[axis] >= end[axis] || end[axis] >= spans[axis].len())
    {
        return Err("THB marked cell box must be strictly interior and nonempty");
    }
    let region =
        std::array::from_fn(|axis| [spans[axis][start[axis]][0], spans[axis][end[axis] - 1][1]]);
    let coarse_active: Vec<_> = (0..prepared_shape.iter().product())
        .filter(|&i| !support_inside(prepared_knots, [3; 3], prepared_shape, i, &region))
        .collect();
    let fine_active: Vec<_> = (0..fine_shape.iter().product())
        .filter(|&i| support_inside(&fine_knots, [3; 3], fine_shape, i, &region))
        .collect();
    if fine_active.is_empty()
        || coarse_active.is_empty()
        || coarse_active.len() + fine_active.len() > 100_000
    {
        return Err("THB marked box has no usable fine basis or exceeds the active basis budget");
    }
    let coarse_to_dof: BTreeMap<_, _> = coarse_active
        .iter()
        .enumerate()
        .map(|(dof, &id)| (id, dof))
        .collect();
    let fine_to_dof: BTreeMap<_, _> = fine_active
        .iter()
        .enumerate()
        .map(|(offset, &id)| (id, coarse_active.len() + offset))
        .collect();
    let mut prepared_to_thb = vec![Row::new(); coarse_active.len() + fine_active.len()];
    for (&id, &dof) in &coarse_to_dof {
        prepared_to_thb[dof].insert(id, Q::one());
    }
    for (&fine_id, &dof) in &fine_to_dof {
        prepared_to_thb[dof] = prepared_to_fine[fine_id].clone();
    }
    let parent_weights: Vec<Q> = support
        .homogeneous_controls
        .iter()
        .map(|control| Q::new(control[3], 1))
        .collect::<std::result::Result<_, _>>()?;
    let initial_weights = apply(&parent_transfer, &parent_weights)?;
    let prepared_weights = apply(&initial_to_prepared, &initial_weights)?;
    let thb_weights = apply(&prepared_to_thb, &prepared_weights)?;
    if thb_weights.iter().any(|weight| *weight <= Q::zero()) {
        return Err("THB rational weights must remain positive");
    }
    let transfer_rank = exact_sparse_rank(&prepared_to_thb)?;
    if transfer_rank != prepared_weights.len() {
        return Err("prepared-to-THB transfer lost prepared-space rank");
    }
    let fine_active_set: BTreeSet<_> = fine_active.iter().copied().collect();
    let mut thb_in_fine = vec![Row::new(); coarse_active.len() + fine_active.len()];
    for (coarse_local, &coarse_id) in coarse_active.iter().enumerate() {
        for (fine_id, insertion) in prepared_to_fine.iter().enumerate() {
            if !fine_active_set.contains(&fine_id) {
                if let Some(&coefficient) = insertion.get(&coarse_id) {
                    if coefficient != Q::zero() {
                        thb_in_fine[coarse_local].insert(fine_id, coefficient);
                    }
                }
            }
        }
    }
    for (fine_local, &fine_id) in fine_active.iter().enumerate() {
        thb_in_fine[coarse_active.len() + fine_local].insert(fine_id, Q::one());
    }
    let uniform_fine_rank = exact_sparse_rank(&thb_in_fine)?;
    if uniform_fine_rank != thb_in_fine.len() {
        return Err("THB basis is not independent in the uniform fine tensor space");
    }
    let p_initial = rational_transfer(&parent_transfer, &parent_weights, &initial_weights)?;
    let p_prepared = rational_transfer(&initial_to_prepared, &initial_weights, &prepared_weights)?;
    let p_thb = rational_transfer(&prepared_to_thb, &prepared_weights, &thb_weights)?;
    let parent_geometry: Vec<[Q; 4]> = support
        .homogeneous_controls
        .iter()
        .map(|control| -> Result<[Q; 4]> {
            Ok([
                Q::new(control[0], 1)?,
                Q::new(control[1], 1)?,
                Q::new(control[2], 1)?,
                Q::new(control[3], 1)?,
            ])
        })
        .collect::<Result<_>>()?;
    let initial_geometry = apply_controls(&parent_transfer, &parent_geometry)?;
    let prepared_geometry = apply_controls(&initial_to_prepared, &initial_geometry)?;
    let thb_geometry = apply_controls(&prepared_to_thb, &prepared_geometry)?;
    let mut leaves = Vec::new();
    let mut boundary = Vec::new();
    let coarse_intervals = [
        unique_intervals(&prepared_knots[0]),
        unique_intervals(&prepared_knots[1]),
        unique_intervals(&prepared_knots[2]),
    ];
    for (iu, u) in coarse_intervals[0].iter().enumerate() {
        for (iv, v) in coarse_intervals[1].iter().enumerate() {
            for (iw, w) in coarse_intervals[2].iter().enumerate() {
                let coarse_cell = [*u, *v, *w];
                let refined = box_contains(&region, &coarse_cell);
                let subintervals = if refined {
                    let mut result = [[[Q::zero(); 2]; 2]; 3];
                    for axis in 0..3 {
                        let mid = midpoint(coarse_cell[axis][0], coarse_cell[axis][1])?;
                        result[axis] = [[coarse_cell[axis][0], mid], [mid, coarse_cell[axis][1]]];
                    }
                    result
                } else {
                    coarse_cell.map(|interval| [interval, interval])
                };
                let children = if refined { 8 } else { 1 };
                for child in 0..children {
                    let cell = if refined {
                        std::array::from_fn(|axis| subintervals[axis][(child >> (2 - axis)) & 1])
                    } else {
                        coarse_cell
                    };
                    let (dofs, extraction) = make_thb_rows(
                        &cell,
                        refined,
                        prepared_knots,
                        &fine_knots,
                        prepared_shape,
                        fine_shape,
                        &coarse_active,
                        &fine_active,
                        &prepared_to_fine,
                        &coarse_to_dof,
                        &fine_to_dof,
                    )?;
                    let prepared_global = exact_cell_rows(prepared_knots, &cell, prepared_shape)?;
                    let prepared_map: BTreeMap<_, _> = prepared_global
                        .iter()
                        .flat_map(|r| r.keys().copied())
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .map(|id| (id, id))
                        .collect();
                    let (prepared_dofs, prepared_extraction) =
                        compact_rows(&prepared_global, &prepared_map)?;
                    let prepared_weights_local = prepared_dofs
                        .iter()
                        .map(|&id| prepared_weights[id])
                        .collect::<Vec<_>>();
                    let field_weights =
                        dofs.iter().map(|&dof| thb_weights[dof]).collect::<Vec<_>>();
                    let parent_cell = cell;
                    let geom_axes = [
                        axis_rows(
                            &support.knots[0].iter().copied().map(q).collect::<Vec<_>>(),
                            support.degrees[0],
                            parent_cell[0],
                        )?,
                        axis_rows(
                            &support.knots[1].iter().copied().map(q).collect::<Vec<_>>(),
                            support.degrees[1],
                            parent_cell[1],
                        )?,
                        axis_rows(
                            &support.knots[2].iter().copied().map(q).collect::<Vec<_>>(),
                            support.degrees[2],
                            parent_cell[2],
                        )?,
                    ];
                    let geom_extraction = tensor_rows(
                        [&geom_axes[0], &geom_axes[1], &geom_axes[2]],
                        support.control_shape,
                    )?;
                    let geometry_controls_exact: Vec<[Q; 4]> = geom_extraction
                        .iter()
                        .map(|row| -> Result<[Q; 4]> {
                            let mut control = [Q::zero(); 4];
                            for (&id, &coefficient) in row {
                                for (c, value) in control.iter_mut().enumerate() {
                                    *value = value.checked_add(coefficient.checked_mul(
                                        Q::new(support.homogeneous_controls[id][c], 1)?,
                                    )?)?;
                                }
                            }
                            Ok(control)
                        })
                        .collect::<Result<_>>()?;
                    let geometry_controls: Vec<[f64; 4]> = geometry_controls_exact
                        .iter()
                        .map(|control| control.map(val))
                        .collect();
                    let geometry = RationalBezierMap::new(support.degrees, geometry_controls)
                        .map_err(|_| "invalid restricted parent geometry")?;
                    let center = geometry
                        .evaluate([0.5; 3])
                        .map_err(|_| "invalid restricted parent geometry")?;
                    if dot(
                        center.jacobian_columns[0],
                        cross(center.jacobian_columns[1], center.jacobian_columns[2]),
                    ) <= 0.
                    {
                        return Err(
                            "reversed or singular source chart is unsupported by this THB pilot",
                        );
                    }
                    let id = format!(
                        "l{}/{}/{}/{}",
                        if refined { 1 } else { 0 },
                        iu,
                        iv,
                        iw * children + child
                    );
                    let leaf_id = leaves.len();
                    for (axis, interval) in cell.iter().enumerate() {
                        if interval[0] == q(support.domain[axis][0]) {
                            boundary.push((leaf_id, axis, false));
                        }
                        if interval[1] == q(support.domain[axis][1]) {
                            boundary.push((leaf_id, axis, true));
                        }
                    }
                    let numeric_extraction = numeric_rows(&extraction);
                    let prepared_numeric_extraction = numeric_rows(&prepared_extraction);
                    leaves.push(ThbLeaf {
                        id,
                        level: if refined { 1 } else { 0 },
                        intervals: cell,
                        dofs,
                        extraction: extraction.clone(),
                        field_weights,
                        prepared_dofs,
                        prepared_extraction,
                        prepared_weights: prepared_weights_local,
                        geometry_controls: geometry_controls_exact,
                        numeric_extraction,
                        numeric_weights: Vec::new(),
                        prepared_numeric_extraction,
                        prepared_numeric_weights: Vec::new(),
                        geometry,
                    });
                    let leaf = leaves.last_mut().unwrap();
                    leaf.numeric_weights = leaf.field_weights.iter().copied().map(val).collect();
                    leaf.prepared_numeric_weights =
                        leaf.prepared_weights.iter().copied().map(val).collect();
                }
            }
        }
    }
    if leaves.is_empty() || leaves.len() > 4096 {
        return Err("THB refined leaf count exceeds the pilot budget");
    }
    // Verify exact field-space nesting on every extracted Bernstein row.
    let elevation = elevation_matrix(1, 3)?;
    for leaf in &leaves {
        let direct_elevated: Vec<[Q; 4]> = (0..64)
            .map(|row| -> Result<[Q; 4]> {
                let i = row / 16;
                let j = (row / 4) % 4;
                let k = row % 4;
                let mut control = [Q::zero(); 4];
                for v in 0..2 {
                    for w in 0..2 {
                        let source = leaf.geometry_controls[(i * 2 + v) * 2 + w];
                        let factor = elevation[j][v].checked_mul(elevation[k][w])?;
                        for c in 0..4 {
                            control[c] = control[c].checked_add(factor.checked_mul(source[c])?)?;
                        }
                    }
                }
                Ok(control)
            })
            .collect::<Result<_>>()?;
        let initial_cell = exact_cell_rows(&initial_knots, &leaf.intervals, initial_shape)?;
        for row in 0..64 {
            let mut initial_control = [Q::zero(); 4];
            for (&id, &coefficient) in &initial_cell[row] {
                for c in 0..4 {
                    initial_control[c] = initial_control[c]
                        .checked_add(coefficient.checked_mul(initial_geometry[id][c])?)?;
                }
            }
            if initial_control != direct_elevated[row] {
                return Err("exact homogeneous parent-to-initial geometry identity failed");
            }
            let mut thb_control = [Q::zero(); 4];
            for (&local, &coefficient) in &leaf.extraction[row] {
                let dof = leaf.dofs[local];
                for c in 0..4 {
                    thb_control[c] = thb_control[c]
                        .checked_add(coefficient.checked_mul(thb_geometry[dof][c])?)?;
                }
            }
            if thb_control != direct_elevated[row] {
                return Err("exact homogeneous prepared-to-THB geometry identity failed");
            }
        }
        for (row, thb_extraction) in leaf.extraction.iter().enumerate() {
            let mut reconstructed = Row::new();
            for (&local, &basis_value) in thb_extraction {
                let dof = leaf.dofs[local];
                for (&prepared_id, &transfer) in &prepared_to_thb[dof] {
                    add(
                        &mut reconstructed,
                        prepared_id,
                        basis_value.checked_mul(transfer)?,
                    )?;
                }
            }
            let old_local = &leaf.prepared_extraction[row];
            let old: Row = old_local
                .iter()
                .map(|(&local, &value)| (leaf.prepared_dofs[local], value))
                .collect();
            if reconstructed != old {
                return Err("exact THB nesting identity failed on a leaf");
            }
        }
    }
    let initial_ndof = initial_shape.iter().product();
    let prepared_ndof = prepared_shape.iter().product();
    let uniform_fine_ndof = fine_shape.iter().product();
    let z_bounds = support
        .homogeneous_controls
        .iter()
        .map(|c| c[2] as f64 / c[3] as f64)
        .fold([f64::INFINITY, f64::NEG_INFINITY], |mut b, z| {
            b[0] = b[0].min(z);
            b[1] = b[1].max(z);
            b
        });
    Ok(ThbStrip {
        leaves,
        source_parent_ndof: support.homogeneous_controls.len(),
        ndof: thb_weights.len(),
        initial_ndof,
        prepared_ndof,
        uniform_fine_ndof,
        transfer_rank,
        uniform_fine_rank,
        coarse_active_count: coarse_active.len(),
        fine_active_count: fine_active.len(),
        parent_to_initial: parent_transfer,
        initial_to_prepared,
        prepared_to_thb,
        rational_parent_to_initial: p_initial,
        rational_initial_to_prepared: p_prepared,
        rational_prepared_to_thb: p_thb,
        weights: thb_weights,
        prepared_weights,
        boundary,
        z_bounds,
        common_partition_verified: true,
        prepared_shape,
        region,
        initial_knots: initial_knots.clone(),
        prepared_knots: prepared_knots.clone(),
        fine_knots,
        coarse_active_ids: coarse_active,
        fine_active_ids: fine_active,
        refined_cell_range: std::array::from_fn(|axis| [start[axis], end[axis]]),
        source_degrees: support.degrees,
        source_knots: support.knots.clone(),
        source_control_shape: support.control_shape,
        source_geometry_controls: support.homogeneous_controls.clone(),
        source_recipe: support.generating_recipe,
    })
}

#[derive(Debug, Clone)]
pub struct ThbLeaf {
    pub id: String,
    pub level: usize,
    pub intervals: [[Q; 2]; 3],
    pub dofs: Vec<usize>,
    /// Bernstein rows by local active THB function.
    pub extraction: Vec<Row>,
    pub field_weights: Vec<Q>,
    /// Prepared-space extraction on the same geometric leaf, before truncation.
    pub prepared_dofs: Vec<usize>,
    pub prepared_extraction: Vec<Row>,
    pub prepared_weights: Vec<Q>,
    pub geometry_controls: Vec<[Q; 4]>,
    numeric_extraction: Vec<Vec<(usize, f64)>>,
    numeric_weights: Vec<f64>,
    prepared_numeric_extraction: Vec<Vec<(usize, f64)>>,
    prepared_numeric_weights: Vec<f64>,
    geometry: RationalBezierMap,
}

impl ThbLeaf {
    pub fn sample(&self, uvw: V) -> Result<(ElementGeometrySample, ElementBasisSample)> {
        if uvw
            .iter()
            .any(|x| !x.is_finite() || !(0. ..=1.).contains(x))
        {
            return Err("leaf parameter outside the closed unit cube");
        }
        let geometry = self
            .geometry
            .evaluate(uvw)
            .map_err(|_| "invalid THB leaf geometry")?;
        let j = geometry.jacobian_columns;
        let cof = [cross(j[1], j[2]), cross(j[2], j[0]), cross(j[0], j[1])];
        let det = dot(j[0], cof[0]);
        if !det.is_finite() || det <= 0. {
            return Err("nonpositive THB geometry Jacobian");
        }
        let geometry_sample = ElementGeometrySample {
            position: geometry.position,
            jacobian_columns: geometry.jacobian_columns,
            determinant: det,
        };
        let field_sample = evaluate_field(
            uvw,
            &self.numeric_extraction,
            &self.numeric_weights,
            &self.dofs,
            cof,
            det,
        )?;
        Ok((geometry_sample, field_sample))
    }

    pub fn sample_prepared(&self, uvw: V) -> Result<(ElementGeometrySample, ElementBasisSample)> {
        if uvw
            .iter()
            .any(|x| !x.is_finite() || !(0. ..=1.).contains(x))
        {
            return Err("leaf parameter outside the closed unit cube");
        }
        let geometry = self
            .geometry
            .evaluate(uvw)
            .map_err(|_| "invalid THB leaf geometry")?;
        let j = geometry.jacobian_columns;
        let cof = [cross(j[1], j[2]), cross(j[2], j[0]), cross(j[0], j[1])];
        let det = dot(j[0], cof[0]);
        if !det.is_finite() || det <= 0. {
            return Err("nonpositive THB geometry Jacobian");
        }
        Ok((
            ElementGeometrySample {
                position: geometry.position,
                jacobian_columns: geometry.jacobian_columns,
                determinant: det,
            },
            evaluate_field(
                uvw,
                &self.prepared_numeric_extraction,
                &self.prepared_numeric_weights,
                &self.prepared_dofs,
                cof,
                det,
            )?,
        ))
    }
}

fn numeric_rows(rows: &[Row]) -> Vec<Vec<(usize, f64)>> {
    rows.iter()
        .map(|row| row.iter().map(|(&i, &c)| (i, val(c))).collect())
        .collect()
}

fn evaluate_field(
    uvw: V,
    extraction: &[Vec<(usize, f64)>],
    weights: &[f64],
    dofs: &[usize],
    cof: [V; 3],
    det: f64,
) -> Result<ElementBasisSample> {
    let basis = [
        bernstein(3, uvw[0]),
        bernstein(3, uvw[1]),
        bernstein(3, uvw[2]),
    ];
    let mut polynomial = vec![0.; dofs.len()];
    let mut derivatives = vec![[0.; 3]; dofs.len()];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                let row_id = (i * 4 + j) * 4 + k;
                let row = &extraction[row_id];
                let b = basis[0].0[i] * basis[1].0[j] * basis[2].0[k];
                let db = [
                    basis[0].1[i] * basis[1].0[j] * basis[2].0[k],
                    basis[0].0[i] * basis[1].1[j] * basis[2].0[k],
                    basis[0].0[i] * basis[1].0[j] * basis[2].1[k],
                ];
                for &(local, c) in row {
                    polynomial[local] += c * b;
                    for axis in 0..3 {
                        derivatives[local][axis] += c * db[axis];
                    }
                }
            }
        }
    }
    let mut denominator = 0.;
    let mut ddenominator = [0.; 3];
    let mut field_values = vec![0.; dofs.len()];
    let mut parametric_gradients = vec![[0.; 3]; dofs.len()];
    for local in 0..dofs.len() {
        let weight = weights[local];
        field_values[local] = weight * polynomial[local];
        denominator += field_values[local];
        for axis in 0..3 {
            parametric_gradients[local][axis] = weight * derivatives[local][axis];
            ddenominator[axis] += parametric_gradients[local][axis];
        }
    }
    if !denominator.is_finite() || denominator <= 0. {
        return Err("nonpositive THB rational field denominator");
    }
    let values: Vec<_> = field_values.iter().map(|x| x / denominator).collect();
    let physical_gradients = parametric_gradients
        .iter()
        .zip(&values)
        .map(|(gradient, value)| {
            let local_gradient = std::array::from_fn::<_, 3, _>(|axis| {
                (gradient[axis] - value * ddenominator[axis]) / denominator
            });
            std::array::from_fn(|component| {
                (local_gradient[0] * cof[0][component]
                    + local_gradient[1] * cof[1][component]
                    + local_gradient[2] * cof[2][component])
                    / det
            })
        })
        .collect();
    Ok(ElementBasisSample {
        dofs: dofs.to_vec(),
        values,
        physical_gradients,
    })
}

fn bernstein(degree: usize, t: f64) -> (Vec<f64>, Vec<f64>) {
    let values = match degree {
        1 => vec![1. - t, t],
        3 => vec![
            (1. - t).powi(3),
            3. * t * (1. - t).powi(2),
            3. * t.powi(2) * (1. - t),
            t.powi(3),
        ],
        _ => unreachable!("THB pilot degrees are cubic"),
    };
    let derivatives = match degree {
        1 => vec![-1., 1.],
        3 => vec![
            -3. * (1. - t).powi(2),
            3. * (1. - t) * (1. - 3. * t),
            3. * t * (2. - 3. * t),
            3. * t.powi(2),
        ],
        _ => unreachable!("THB pilot degrees are cubic"),
    };
    (values, derivatives)
}

fn dot(a: V, b: V) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// One strict interior THB refinement of the exact NURBS strip support.
pub struct ThbStrip {
    pub leaves: Vec<ThbLeaf>,
    pub source_parent_ndof: usize,
    pub ndof: usize,
    pub initial_ndof: usize,
    pub prepared_ndof: usize,
    pub uniform_fine_ndof: usize,
    pub transfer_rank: usize,
    pub uniform_fine_rank: usize,
    pub coarse_active_count: usize,
    pub fine_active_count: usize,
    pub parent_to_initial: Sparse,
    pub initial_to_prepared: Sparse,
    pub prepared_to_thb: Sparse,
    pub rational_parent_to_initial: Sparse,
    pub rational_initial_to_prepared: Sparse,
    pub rational_prepared_to_thb: Sparse,
    pub weights: Vec<Q>,
    pub prepared_weights: Vec<Q>,
    pub boundary: Vec<(usize, usize, bool)>,
    pub z_bounds: [f64; 2],
    pub common_partition_verified: bool,
    pub prepared_shape: [usize; 3],
    pub region: [[Q; 2]; 3],
    pub initial_knots: [Vec<Q>; 3],
    pub prepared_knots: [Vec<Q>; 3],
    pub fine_knots: [Vec<Q>; 3],
    pub coarse_active_ids: Vec<usize>,
    pub fine_active_ids: Vec<usize>,
    pub refined_cell_range: [[usize; 2]; 3],
    pub source_degrees: [usize; 3],
    pub source_knots: [Vec<i64>; 3],
    pub source_control_shape: [usize; 3],
    pub source_geometry_controls: Vec<[i128; 4]>,
    pub source_recipe: Json,
}

impl ThbStrip {
    pub fn sample(
        &self,
        leaf: usize,
        uvw: V,
    ) -> Result<(ElementGeometrySample, ElementBasisSample)> {
        self.leaves
            .get(leaf)
            .ok_or("THB leaf index out of range")?
            .sample(uvw)
    }

    pub fn sample_prepared(
        &self,
        leaf: usize,
        uvw: V,
    ) -> Result<(ElementGeometrySample, ElementBasisSample)> {
        self.leaves
            .get(leaf)
            .ok_or("THB leaf index out of range")?
            .sample_prepared(uvw)
    }

    /// Homogeneous zero Dirichlet trace on the physical bottom (w=0) in the
    /// uniformly prepared tensor basis. Tensor indexing is w-fast.
    pub fn prepared_zero_bottom_lift(&self) -> Result<Lift> {
        let mut constrained = BTreeSet::new();
        for i in 0..self.prepared_shape[0] {
            for j in 0..self.prepared_shape[1] {
                constrained.insert(index(self.prepared_shape, [i, j, 0]));
            }
        }
        let free_roots: Vec<_> = (0..self.prepared_ndof)
            .filter(|i| !constrained.contains(i))
            .collect();
        let lines = (0..self.prepared_ndof)
            .map(|dof| LiftLine {
                dof,
                terms: free_roots
                    .binary_search(&dof)
                    .map(|root| vec![(root, Q::one())])
                    .unwrap_or_default(),
                offset: Q::zero(),
            })
            .collect();
        Ok(Lift { free_roots, lines })
    }

    pub fn zero_face_lift(&self, face: &[(usize, usize, bool)]) -> Result<Lift> {
        let mut constrained = BTreeSet::new();
        for &(leaf_id, axis, upper) in face {
            let leaf = self
                .leaves
                .get(leaf_id)
                .ok_or("THB boundary leaf out of range")?;
            if axis >= 3 {
                return Err("THB boundary axis outside 0..3");
            }
            let bernstein_index = if upper { 3 } else { 0 };
            let tangential: Vec<_> = (0..3).filter(|&candidate| candidate != axis).collect();
            for a in 0..4 {
                for b in 0..4 {
                    let mut at = [0; 3];
                    at[axis] = bernstein_index;
                    at[tangential[0]] = a;
                    at[tangential[1]] = b;
                    let row = (at[0] * 4 + at[1]) * 4 + at[2];
                    constrained.extend(leaf.extraction[row].keys().map(|&local| leaf.dofs[local]));
                }
            }
        }
        let free_roots: Vec<_> = (0..self.ndof)
            .filter(|i| !constrained.contains(i))
            .collect();
        let mut lines = Vec::with_capacity(self.ndof);
        for dof in 0..self.ndof {
            let terms = free_roots
                .binary_search(&dof)
                .map(|root| vec![(root, Q::one())])
                .unwrap_or_default();
            lines.push(LiftLine {
                dof,
                terms,
                offset: Q::zero(),
            });
        }
        Ok(Lift { free_roots, lines })
    }

    pub fn transfer_field(&self, coarse: &[Q]) -> Result<Vec<Q>> {
        apply(&self.rational_prepared_to_thb, coarse)
    }

    pub fn document(&self) -> Json {
        let exact = |value: Q| {
            Json::Array(vec![
                Json::String(value.numerator().to_string()),
                Json::String(value.denominator().to_string()),
            ])
        };
        let sparse = |matrix: &[Row]| {
            Json::Array(
                matrix
                    .iter()
                    .map(|row| {
                        Json::Array(
                            row.iter()
                                .map(|(&column, &value)| {
                                    Json::Array(vec![
                                        Json::Number(column.to_string()),
                                        exact(value),
                                    ])
                                })
                                .collect(),
                        )
                    })
                    .collect(),
            )
        };
        let exact_knots = |axes: &[Vec<Q>; 3]| {
            Json::Array(
                axes.iter()
                    .map(|axis| Json::Array(axis.iter().map(|&x| exact(x)).collect()))
                    .collect(),
            )
        };
        object([
            ("schema", Json::String("thb_strip_exact_v1".into())),
            (
                "profile",
                Json::String("interior_dyadic_thb_cubic_v1".into()),
            ),
            (
                "basis_name",
                Json::String("strictly truncated hierarchical B-spline".into()),
            ),
            (
                "degrees",
                Json::Array(
                    vec![3, 3, 3]
                        .into_iter()
                        .map(|x| Json::Number(x.to_string()))
                        .collect(),
                ),
            ),
            (
                "storage_order",
                Json::String("w-fast, then v, then u".into()),
            ),
            ("source_recipe", self.source_recipe.clone()),
            (
                "source_degrees",
                Json::Array(
                    self.source_degrees
                        .iter()
                        .map(|x| Json::Number(x.to_string()))
                        .collect(),
                ),
            ),
            (
                "source_knots",
                Json::Array(
                    self.source_knots
                        .iter()
                        .map(|axis| {
                            Json::Array(axis.iter().map(|x| Json::Number(x.to_string())).collect())
                        })
                        .collect(),
                ),
            ),
            (
                "source_control_shape",
                Json::Array(
                    self.source_control_shape
                        .iter()
                        .map(|x| Json::Number(x.to_string()))
                        .collect(),
                ),
            ),
            (
                "source_geometry_controls_homogeneous",
                Json::Array(
                    self.source_geometry_controls
                        .iter()
                        .map(|control| {
                            Json::Array(
                                control
                                    .iter()
                                    .map(|x| Json::String(x.to_string()))
                                    .collect(),
                            )
                        })
                        .collect(),
                ),
            ),
            (
                "initial_scalar_dofs",
                Json::Number(self.initial_ndof.to_string()),
            ),
            (
                "source_parent_scalar_dofs",
                Json::Number(self.source_parent_ndof.to_string()),
            ),
            (
                "prepared_scalar_dofs",
                Json::Number(self.prepared_ndof.to_string()),
            ),
            ("local_scalar_dofs", Json::Number(self.ndof.to_string())),
            (
                "uniform_fine_scalar_dofs",
                Json::Number(self.uniform_fine_ndof.to_string()),
            ),
            (
                "coarse_active",
                Json::Number(self.coarse_active_count.to_string()),
            ),
            (
                "fine_active",
                Json::Number(self.fine_active_count.to_string()),
            ),
            (
                "prepared_to_thb_rank",
                Json::Number(self.transfer_rank.to_string()),
            ),
            (
                "uniform_fine_basis_rank",
                Json::Number(self.uniform_fine_rank.to_string()),
            ),
            (
                "field_weights",
                Json::Array(self.weights.iter().copied().map(exact).collect()),
            ),
            (
                "prepared_shape",
                Json::Array(
                    self.prepared_shape
                        .iter()
                        .map(|x| Json::Number(x.to_string()))
                        .collect(),
                ),
            ),
            (
                "coarse_active_ids",
                Json::Array(
                    self.coarse_active_ids
                        .iter()
                        .map(|x| Json::Number(x.to_string()))
                        .collect(),
                ),
            ),
            (
                "fine_active_ids",
                Json::Array(
                    self.fine_active_ids
                        .iter()
                        .map(|x| Json::Number(x.to_string()))
                        .collect(),
                ),
            ),
            (
                "refined_cell_range",
                Json::Array(
                    self.refined_cell_range
                        .iter()
                        .map(|range| {
                            Json::Array(range.iter().map(|x| Json::Number(x.to_string())).collect())
                        })
                        .collect(),
                ),
            ),
            ("initial_knots", exact_knots(&self.initial_knots)),
            ("prepared_knots", exact_knots(&self.prepared_knots)),
            ("fine_knots", exact_knots(&self.fine_knots)),
            (
                "prepared_weights",
                Json::Array(self.prepared_weights.iter().copied().map(exact).collect()),
            ),
            (
                "leaves",
                Json::Array(
                    self.leaves
                        .iter()
                        .map(|leaf| {
                            object([
                                ("id", Json::String(leaf.id.clone())),
                                ("level", Json::Number(leaf.level.to_string())),
                                (
                                    "intervals",
                                    Json::Array(
                                        leaf.intervals
                                            .iter()
                                            .map(|interval| {
                                                Json::Array(
                                                    interval.iter().copied().map(exact).collect(),
                                                )
                                            })
                                            .collect(),
                                    ),
                                ),
                                (
                                    "dofs",
                                    Json::Array(
                                        leaf.dofs
                                            .iter()
                                            .map(|x| Json::Number(x.to_string()))
                                            .collect(),
                                    ),
                                ),
                                ("extraction", sparse(&leaf.extraction)),
                                (
                                    "field_weights",
                                    Json::Array(
                                        leaf.field_weights.iter().copied().map(exact).collect(),
                                    ),
                                ),
                                (
                                    "prepared_dofs",
                                    Json::Array(
                                        leaf.prepared_dofs
                                            .iter()
                                            .map(|x| Json::Number(x.to_string()))
                                            .collect(),
                                    ),
                                ),
                                ("prepared_extraction", sparse(&leaf.prepared_extraction)),
                                (
                                    "prepared_field_weights",
                                    Json::Array(
                                        leaf.prepared_weights.iter().copied().map(exact).collect(),
                                    ),
                                ),
                                (
                                    "geometry_controls_homogeneous",
                                    Json::Array(
                                        leaf.geometry_controls
                                            .iter()
                                            .map(|control| {
                                                Json::Array(
                                                    control.iter().copied().map(exact).collect(),
                                                )
                                            })
                                            .collect(),
                                    ),
                                ),
                            ])
                        })
                        .collect(),
                ),
            ),
            (
                "rational_transfer_prepared_to_thb",
                sparse(&self.rational_prepared_to_thb),
            ),
            (
                "rational_transfer_parent_to_initial",
                sparse(&self.rational_parent_to_initial),
            ),
            (
                "rational_transfer_initial_to_prepared",
                sparse(&self.rational_initial_to_prepared),
            ),
            ("exact_parent_to_initial", sparse(&self.parent_to_initial)),
            (
                "exact_initial_to_prepared",
                sparse(&self.initial_to_prepared),
            ),
            ("exact_prepared_to_thb", sparse(&self.prepared_to_thb)),
            (
                "common_partition_verified",
                Json::Bool(self.common_partition_verified),
            ),
            (
                "interior_region",
                Json::Array(
                    self.region
                        .iter()
                        .map(|interval| Json::Array(interval.iter().copied().map(exact).collect()))
                        .collect(),
                ),
            ),
        ])
    }
}
