//! Linear-elastic solve for the audited fixed bracket.

use crate::analysis_space::{AnalysisBlock, AnalysisSpace};
use crate::elasticity::{assemble, basis_geometry, gauss, Matrix, Stiffness};
use std::collections::BTreeSet;

const YOUNG_MPA: f64 = 210_000.0;
const POISSON: f64 = 0.3;
const RESULTANT_N: f64 = 1_000.0;
const LOADED_AREA_MM2: f64 = 320.0;
const PCG_TOLERANCE: f64 = 1e-10;
const MAX_PCG_ITERATIONS: usize = 20_000;

/// Measurable outputs of the fixed bracket elasticity solve.
#[derive(Debug, Clone, PartialEq)]
pub struct BracketResult {
    pub volume_mm3: f64,
    pub load_area_x60_mm2: f64,
    pub ndof: usize,
    pub pcg_iterations: usize,
    pub relative_residual: f64,
    pub applied_force_n: [f64; 3],
    pub reaction_force_n: [f64; 3],
    pub strain_energy_n_mm: f64,
    pub external_work_n_mm: f64,
    pub mean_y_displacement_mm: f64,
    pub point_displacement_mm: [f64; 3],
    pub point_y_displacement_mm: f64,
    /// Global rational Q2 displacement coefficients, scalar DOF then x/y/z.
    pub displacement_mm: Vec<f64>,
}

#[derive(Clone, Copy)]
struct Face {
    block: usize,
    axis: usize,
    upper: bool,
}

fn flat_index(i: usize, j: usize, k: usize) -> usize {
    (i * 3 + j) * 3 + k
}

fn coordinate_on_face(
    block: &AnalysisBlock,
    face: Face,
    coordinate: usize,
    value: i128,
) -> Result<bool, &'static str> {
    let fixed = if face.upper { 2 } else { 0 };
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let index = [i, j, k];
                if index[face.axis] == fixed {
                    let point = block.controls[flat_index(i, j, k)];
                    let expected = value
                        .checked_mul(point[3])
                        .ok_or("loaded-face plane comparison overflow")?;
                    if point[coordinate] != expected {
                        return Ok(false);
                    }
                }
            }
        }
    }
    Ok(true)
}

fn loaded_faces(space: &AnalysisSpace) -> Result<Vec<Face>, &'static str> {
    let mut result = Vec::new();
    for (block_index, block) in space.blocks.iter().enumerate() {
        if block.controls.len() != 27 || block.dofs.len() != 27 {
            return Err("invalid analysis block in solve");
        }
        for axis in 0..3 {
            for upper in [false, true] {
                let face = Face {
                    block: block_index,
                    axis,
                    upper,
                };
                if coordinate_on_face(block, face, 0, 60)? {
                    result.push(face);
                }
            }
        }
    }
    if result.is_empty() {
        return Err("no boundary faces lie on x=60");
    }
    Ok(result)
}

fn free_axes(fixed: usize) -> [usize; 2] {
    match fixed {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(vector: [f64; 3]) -> f64 {
    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt()
}

fn integrate_load(
    space: &AnalysisSpace,
    order: usize,
    faces: &[Face],
) -> Result<(Vec<f64>, f64), &'static str> {
    let points = gauss(order)?;
    let mut force = vec![0.0; 3 * space.ndof];
    let mut area = 0.0;
    for &face in faces {
        let block = &space.blocks[face.block];
        let [axis_a, axis_b] = free_axes(face.axis);
        for &(qa, wa) in &points {
            for &(qb, wb) in &points {
                let mut uvw = [0.0; 3];
                uvw[face.axis] = if face.upper { 1.0 } else { 0.0 };
                uvw[axis_a] = qa;
                uvw[axis_b] = qb;
                let sample = basis_geometry(block, uvw)?;
                let surface_jacobian = norm(cross(
                    sample.jacobian_columns[axis_a],
                    sample.jacobian_columns[axis_b],
                ));
                let measure = surface_jacobian * wa * wb;
                if !measure.is_finite() || measure <= 0.0 {
                    return Err("invalid loaded-face area element");
                }
                area += measure;
                let traction = -RESULTANT_N / LOADED_AREA_MM2;
                for local in 0..27 {
                    let dof = block.dofs[local];
                    if dof >= space.ndof {
                        return Err("analysis degree of freedom is out of range");
                    }
                    force[3 * dof + 1] += sample.values[local] * traction * measure;
                }
            }
        }
    }
    if !area.is_finite() || (area - LOADED_AREA_MM2).abs() > 1e-8 {
        return Err("x60 loaded-face area does not equal 320 mm^2");
    }
    let total_y: f64 = force.iter().skip(1).step_by(3).sum();
    let other_force: f64 = force
        .iter()
        .enumerate()
        .filter(|(index, _)| index % 3 != 1)
        .map(|(_, value)| value.abs())
        .sum();
    if (total_y + RESULTANT_N).abs() > 1e-8 || other_force > 1e-10 {
        return Err("x60 traction resultant does not equal -1000 N in y");
    }
    Ok((force, area))
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn inverse_3x3(matrix: Matrix) -> Result<Matrix, &'static str> {
    let mut augmented = [[0.0; 6]; 3];
    for row in 0..3 {
        for column in 0..3 {
            augmented[row][column] = matrix[row][column];
        }
        augmented[row][row + 3] = 1.0;
    }
    for column in 0..3 {
        let pivot = (column..3)
            .max_by(|&a, &b| {
                augmented[a][column]
                    .abs()
                    .total_cmp(&augmented[b][column].abs())
            })
            .ok_or("empty diagonal block")?;
        if !augmented[pivot][column].is_finite() || augmented[pivot][column].abs() <= 1e-30 {
            return Err("singular block-Jacobi diagonal");
        }
        augmented.swap(column, pivot);
        let scale = augmented[column][column];
        for value in &mut augmented[column] {
            *value /= scale;
        }
        let normalized = augmented[column];
        for (row, values) in augmented.iter_mut().enumerate() {
            if row == column {
                continue;
            }
            let factor = values[column];
            for (value, pivot_value) in values.iter_mut().zip(normalized) {
                *value -= factor * pivot_value;
            }
        }
    }
    let inverse =
        std::array::from_fn(|row| std::array::from_fn(|column| augmented[row][column + 3]));
    if inverse.iter().flatten().all(|value| value.is_finite()) {
        Ok(inverse)
    } else {
        Err("nonfinite block-Jacobi inverse")
    }
}

fn apply_block_preconditioner(residual: &[f64], inverse: &[Matrix]) -> Vec<f64> {
    let mut result = vec![0.0; residual.len()];
    for (node, matrix) in inverse.iter().enumerate() {
        let source = [
            residual[3 * node],
            residual[3 * node + 1],
            residual[3 * node + 2],
        ];
        for component in 0..3 {
            result[3 * node + component] = matrix[component][0] * source[0]
                + matrix[component][1] * source[1]
                + matrix[component][2] * source[2];
        }
    }
    result
}

fn apply_restricted(
    stiffness: &Stiffness,
    ndof: usize,
    free: &[usize],
    vector: &[f64],
) -> Result<Vec<f64>, &'static str> {
    let mut full = vec![0.0; 3 * ndof];
    for (local, &global) in free.iter().enumerate() {
        full[3 * global..3 * global + 3].copy_from_slice(&vector[3 * local..3 * local + 3]);
    }
    let product = stiffness.apply(&full)?;
    Ok(free
        .iter()
        .flat_map(|&global| product[3 * global..3 * global + 3].iter().copied())
        .collect())
}

fn diagonal_blocks(stiffness: &Stiffness, ndof: usize) -> Result<Vec<Matrix>, &'static str> {
    if stiffness.offsets.len() != ndof + 1 {
        return Err("invalid stiffness row offsets");
    }
    let mut diagonal = vec![[[0.0; 3]; 3]; ndof];
    for (row, block) in diagonal.iter_mut().enumerate() {
        for entry in stiffness.offsets[row]..stiffness.offsets[row + 1] {
            if stiffness.columns[entry] == row {
                *block = stiffness.values[entry];
                break;
            }
        }
    }
    Ok(diagonal)
}

fn pcg(
    stiffness: &Stiffness,
    ndof: usize,
    free: &[usize],
    rhs: &[f64],
) -> Result<(Vec<f64>, usize, f64), &'static str> {
    if rhs.len() != 3 * free.len() || rhs.is_empty() {
        return Err("invalid free load vector");
    }
    let diagonal = diagonal_blocks(stiffness, ndof)?;
    let inverse = free
        .iter()
        .map(|&dof| inverse_3x3(diagonal[dof]))
        .collect::<Result<Vec<_>, _>>()?;

    let rhs_norm = dot(rhs, rhs).sqrt();
    if !rhs_norm.is_finite() || rhs_norm <= 0.0 {
        return Err("free load has zero or nonfinite norm");
    }
    let mut solution = vec![0.0; rhs.len()];
    let mut residual = rhs.to_vec();
    let mut z = apply_block_preconditioner(&residual, &inverse);
    let mut direction = z.clone();
    let mut rz = dot(&residual, &z);
    if !rz.is_finite() || rz <= 0.0 {
        return Err("block-Jacobi preconditioner is not positive");
    }

    for iteration in 1..=MAX_PCG_ITERATIONS {
        let product = apply_restricted(stiffness, ndof, free, &direction)?;
        let denominator = dot(&direction, &product);
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err("PCG encountered a nonpositive curvature direction");
        }
        let alpha = rz / denominator;
        for index in 0..solution.len() {
            solution[index] += alpha * direction[index];
            residual[index] -= alpha * product[index];
        }
        let relative = dot(&residual, &residual).sqrt() / rhs_norm;
        if !relative.is_finite() {
            return Err("PCG residual became nonfinite");
        }
        if relative <= PCG_TOLERANCE {
            let actual_product = apply_restricted(stiffness, ndof, free, &solution)?;
            for index in 0..residual.len() {
                residual[index] = rhs[index] - actual_product[index];
            }
            let actual_relative = dot(&residual, &residual).sqrt() / rhs_norm;
            if actual_relative <= PCG_TOLERANCE {
                return Ok((solution, iteration, actual_relative));
            }
            z = apply_block_preconditioner(&residual, &inverse);
            direction.clone_from(&z);
            rz = dot(&residual, &z);
            if !rz.is_finite() || rz <= 0.0 {
                return Err("PCG true residual is not positive");
            }
            continue;
        }
        z = apply_block_preconditioner(&residual, &inverse);
        let next_rz = dot(&residual, &z);
        if !next_rz.is_finite() || next_rz <= 0.0 {
            return Err("PCG preconditioned residual is not positive");
        }
        let beta = next_rz / rz;
        for index in 0..direction.len() {
            direction[index] = z[index] + beta * direction[index];
        }
        rz = next_rz;
    }
    Err("PCG failed to reach relative residual tolerance")
}

fn face_grid_control(block: &AnalysisBlock, face: Face, a: usize, b: usize) -> [i128; 4] {
    let free = free_axes(face.axis);
    let mut index = [0usize; 3];
    index[face.axis] = if face.upper { 2 } else { 0 };
    index[free[0]] = a;
    index[free[1]] = b;
    block.controls[flat_index(index[0], index[1], index[2])]
}

fn physical(control: [i128; 4]) -> [f64; 3] {
    let weight = control[3] as f64;
    [
        control[0] as f64 / weight,
        control[1] as f64 / weight,
        control[2] as f64 / weight,
    ]
}

fn point_on_axis_aligned_box_face(
    block: &AnalysisBlock,
    face: Face,
    target: [f64; 3],
) -> Result<Option<[f64; 3]>, &'static str> {
    let weight = block.controls[0][3];
    if weight <= 0 || block.controls.iter().any(|point| point[3] != weight) {
        return Ok(None);
    }
    let free = free_axes(face.axis);
    let mut parameters = [0.0; 3];
    parameters[face.axis] = if face.upper { 1.0 } else { 0.0 };

    let p00 = physical(face_grid_control(block, face, 0, 0));
    let p20 = physical(face_grid_control(block, face, 2, 0));
    let p02 = physical(face_grid_control(block, face, 0, 2));
    let mut changed_dimensions = BTreeSet::new();
    for (parameter_axis, p0, p2) in [(free[0], p00, p20), (free[1], p00, p02)] {
        let deltas = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let changed: Vec<usize> = (0..3).filter(|&axis| deltas[axis].abs() > 1e-12).collect();
        if changed.len() != 1 || !changed_dimensions.insert(changed[0]) {
            return Ok(None);
        }
        let component = changed[0];
        parameters[parameter_axis] = (target[component] - p0[component]) / deltas[component];
        if !(-1e-12..=1.0 + 1e-12).contains(&parameters[parameter_axis]) {
            return Ok(None);
        }
    }

    // This inverse is deliberately restricted to the known affine end boxes;
    // it does not attempt Newton inversion on a curved/refined geometry.
    for a in 0..=2 {
        for b in 0..=2 {
            let point = physical(face_grid_control(block, face, a, b));
            let expected = [
                p00[0]
                    + (a as f64 / 2.0) * (p20[0] - p00[0])
                    + (b as f64 / 2.0) * (p02[0] - p00[0]),
                p00[1]
                    + (a as f64 / 2.0) * (p20[1] - p00[1])
                    + (b as f64 / 2.0) * (p02[1] - p00[1]),
                p00[2]
                    + (a as f64 / 2.0) * (p20[2] - p00[2])
                    + (b as f64 / 2.0) * (p02[2] - p00[2]),
            ];
            if (0..3).any(|axis| (point[axis] - expected[axis]).abs() > 1e-9) {
                return Ok(None);
            }
        }
    }

    let sample = basis_geometry(block, parameters)?;
    if (0..3).any(|axis| (sample.position[axis] - target[axis]).abs() > 1e-8) {
        return Ok(None);
    }
    Ok(Some(parameters))
}

fn point_displacement(
    space: &AnalysisSpace,
    faces: &[Face],
    full_solution: &[f64],
    target: [f64; 3],
) -> Result<[f64; 3], &'static str> {
    let mut values = Vec::new();
    for &face in faces {
        let block = &space.blocks[face.block];
        if let Some(parameters) = point_on_axis_aligned_box_face(block, face, target)? {
            let sample = basis_geometry(block, parameters)?;
            let displacement = std::array::from_fn(|component| {
                (0..27)
                    .map(|local| {
                        sample.values[local] * full_solution[3 * block.dofs[local] + component]
                    })
                    .sum::<f64>()
            });
            values.push(displacement);
        }
    }
    let first = *values
        .first()
        .ok_or("point is not on a loaded affine end face")?;
    if values
        .iter()
        .skip(1)
        .any(|other| (0..3).any(|component| (other[component] - first[component]).abs() > 1e-8))
    {
        return Err("point interpolation disagrees across loaded faces");
    }
    Ok(first)
}

/// Solve the audited bracket with a uniform total load of 1000 N in `-y`.
pub fn solve_bracket(space: &AnalysisSpace, order: usize) -> Result<BracketResult, &'static str> {
    if space.ndof == 0 || space.blocks.is_empty() {
        return Err("analysis space is empty");
    }
    let fixed_dofs: BTreeSet<usize> = space
        .boundary_dofs
        .get("y50")
        .ok_or("analysis space has no y50 Dirichlet boundary")?
        .iter()
        .copied()
        .collect();
    if fixed_dofs.is_empty() || fixed_dofs.iter().any(|dof| *dof >= space.ndof) {
        return Err("invalid y50 Dirichlet set");
    }

    let faces = loaded_faces(space)?;
    let stiffness = assemble(space, YOUNG_MPA, POISSON, order)?;
    let (load, area) = integrate_load(space, order, &faces)?;
    let applied_force =
        std::array::from_fn(|component| load.iter().skip(component).step_by(3).sum::<f64>());

    let free: Vec<usize> = (0..space.ndof)
        .filter(|dof| !fixed_dofs.contains(dof))
        .collect();
    if free.is_empty() {
        return Err("Dirichlet boundary fixes every degree of freedom");
    }
    let rhs: Vec<f64> = free
        .iter()
        .flat_map(|&dof| load[3 * dof..3 * dof + 3].iter().copied())
        .collect();
    let (free_solution, iterations, relative_residual) = pcg(&stiffness, space.ndof, &free, &rhs)?;
    let mut solution = vec![0.0; 3 * space.ndof];
    for (local, &global) in free.iter().enumerate() {
        solution[3 * global..3 * global + 3]
            .copy_from_slice(&free_solution[3 * local..3 * local + 3]);
    }

    let product = stiffness.apply(&solution)?;
    let mut reaction_force = [0.0; 3];
    for dof in fixed_dofs.iter().copied() {
        for component in 0..3 {
            reaction_force[component] += product[3 * dof + component] - load[3 * dof + component];
        }
    }
    let external_work = dot(&load, &solution);
    let strain_energy = 0.5 * dot(&solution, &product);
    let mean_y = external_work / -RESULTANT_N;
    let displacement_at_point = point_displacement(space, &faces, &solution, [60.0, 4.0, 20.0])?;
    if !external_work.is_finite()
        || !strain_energy.is_finite()
        || !mean_y.is_finite()
        || displacement_at_point.iter().any(|value| !value.is_finite())
    {
        return Err("nonfinite elasticity result");
    }
    let imbalance: f64 = (0..3)
        .map(|a| (reaction_force[a] + applied_force[a]).powi(2))
        .sum();
    if imbalance.sqrt() > 1e-8 * RESULTANT_N
        || external_work <= 0.0
        || (external_work - 2.0 * strain_energy).abs() > 1e-8 * external_work
    {
        return Err("force or work balance outside the fixed solve tolerance");
    }

    Ok(BracketResult {
        volume_mm3: stiffness.integrated_volume,
        load_area_x60_mm2: area,
        ndof: space.ndof,
        pcg_iterations: iterations,
        relative_residual,
        applied_force_n: applied_force,
        reaction_force_n: reaction_force,
        strain_energy_n_mm: strain_energy,
        external_work_n_mm: external_work,
        mean_y_displacement_mm: mean_y,
        point_displacement_mm: displacement_at_point,
        point_y_displacement_mm: displacement_at_point[1],
        displacement_mm: solution,
    })
}
