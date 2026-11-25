//! Fixed-regime radius sensitivity for the audited bracket.
//!
//! The x-bore radius is differentiated in physical millimetres at 4 mm.
//! Independent reconstructions at symmetric radius offsets check the
//! integrated stiffness derivative, compliance and loaded-face mean motion.

use solidcraft_trivariate::{
    analysis_space::{parameterized_bracket_space, AnalysisBlock, AnalysisSpace},
    bracket_parameters::build_bracket_parameters,
    bracket_solve::bracket_load,
    constrained_elasticity::solve_constrained,
    constraint_system::{resolve, Lift, Ratio},
    elasticity::{assemble, basis_geometry, Stiffness},
    elasticity_sensitivity::{assemble_sensitivity, MaterialJet},
    mechanical_sensitivity::{differentiate, LiftRate, MechanicalDerivative, SensitivityResult},
    stress::fillet_limits,
};
use std::collections::BTreeMap;
use std::time::Instant;

const YOUNG_MPA: f64 = 210_000.0;
const POISSON: f64 = 0.3;
const LOAD_RESULTANT_N: f64 = -1_000.0;
const BASE_RADIUS_TICKS: u8 = 64;
const RADIUS_OFFSETS_TICKS: [u8; 4] = [8, 4, 2, 1];

fn elevate_axis(
    values: Vec<[f64; 4]>,
    mut degrees: [usize; 3],
    axis: usize,
) -> Result<(Vec<[f64; 4]>, [usize; 3]), &'static str> {
    if degrees[axis] != 1 {
        return Ok((values, degrees));
    }
    let old_dims = degrees.map(|degree| degree + 1);
    let mut new_degrees = degrees;
    new_degrees[axis] = 2;
    let new_dims = new_degrees.map(|degree| degree + 1);
    let mut output = Vec::with_capacity(27);
    for i in 0..new_dims[0] {
        for j in 0..new_dims[1] {
            for k in 0..new_dims[2] {
                let index = [i, j, k];
                let coordinate = index[axis];
                let (mut left, mut right) = (index, index);
                let value = match coordinate {
                    0 => {
                        left[axis] = 0;
                        values[(left[0] * old_dims[1] + left[1]) * old_dims[2] + left[2]]
                            .map(|entry| 2.0 * entry)
                    }
                    1 => {
                        left[axis] = 0;
                        right[axis] = 1;
                        let a = values[(left[0] * old_dims[1] + left[1]) * old_dims[2] + left[2]];
                        let b =
                            values[(right[0] * old_dims[1] + right[1]) * old_dims[2] + right[2]];
                        std::array::from_fn(|component| a[component] + b[component])
                    }
                    _ => {
                        left[axis] = 1;
                        values[(left[0] * old_dims[1] + left[1]) * old_dims[2] + left[2]]
                            .map(|entry| 2.0 * entry)
                    }
                };
                output.push(value);
            }
        }
    }
    degrees = new_degrees;
    if output.len() != new_dims.iter().product::<usize>() {
        return Err("invalid elevated sensitivity net");
    }
    Ok((output, degrees))
}

fn elevate_jet(
    controls: Vec<[f64; 4]>,
    degrees: [usize; 3],
) -> Result<Vec<[f64; 4]>, &'static str> {
    let mut values = controls;
    let mut shape = degrees;
    // Match the analysis-space axis order and its unnormalised linear maps.
    for axis in 0..3 {
        (values, shape) = elevate_axis(values, shape, axis)?;
    }
    if shape != [2, 2, 2] {
        return Err("sensitivity net did not reach Q2");
    }
    Ok(values)
}

fn elevate_integer_controls(
    block: &solidcraft_trivariate::bracket::BracketBlock,
) -> Result<Vec<[i128; 4]>, &'static str> {
    let mut values = block
        .controls
        .iter()
        .map(|p| p.map(i128::from))
        .collect::<Vec<_>>();
    let mut degrees = block.degrees;
    for axis in 0..3 {
        if degrees[axis] == 2 {
            continue;
        }
        let old_dims = degrees.map(|degree| degree + 1);
        let mut new_degrees = degrees;
        new_degrees[axis] = 2;
        let new_dims = new_degrees.map(|degree| degree + 1);
        let mut output = Vec::with_capacity(new_dims.iter().product());
        for i in 0..new_dims[0] {
            for j in 0..new_dims[1] {
                for k in 0..new_dims[2] {
                    let mut index = [i, j, k];
                    let coordinate = index[axis];
                    let value = match coordinate {
                        0 => {
                            index[axis] = 0;
                            values[(index[0] * old_dims[1] + index[1]) * old_dims[2] + index[2]]
                                .map(|value| 2 * value)
                        }
                        1 => {
                            let mut right = index;
                            index[axis] = 0;
                            right[axis] = 1;
                            let left_value = values
                                [(index[0] * old_dims[1] + index[1]) * old_dims[2] + index[2]];
                            let right_value = values
                                [(right[0] * old_dims[1] + right[1]) * old_dims[2] + right[2]];
                            std::array::from_fn(|component| {
                                left_value[component] + right_value[component]
                            })
                        }
                        _ => {
                            index[axis] = 1;
                            values[(index[0] * old_dims[1] + index[1]) * old_dims[2] + index[2]]
                                .map(|value| 2 * value)
                        }
                    };
                    output.push(value);
                }
            }
        }
        degrees = new_degrees;
        values = output;
    }
    Ok(values)
}

struct Model {
    space: AnalysisSpace,
    radius_x_jet: Vec<Vec<[f64; 4]>>,
}

fn model(radius_x_ticks: u8) -> Result<Model, &'static str> {
    let parameterized = build_bracket_parameters(radius_x_ticks, BASE_RADIUS_TICKS)?;
    let space =
        parameterized_bracket_space(radius_x_ticks, BASE_RADIUS_TICKS, &parameterized.blocks)?;
    let source_blocks: BTreeMap<_, _> = parameterized
        .blocks
        .iter()
        .map(|block| (block.id.as_str(), block))
        .collect();
    let radius_x_jet = space
        .blocks
        .iter()
        .map(|block| {
            let source = *source_blocks
                .get(block.id.as_str())
                .ok_or("analysis space block has no source geometry")?;
            let raw = parameterized.derivatives[0]
                .get(&block.id)
                .cloned()
                .unwrap_or_else(|| vec![[0_i64; 4]; block.controls.len()]);
            elevate_jet(
                raw.into_iter().map(|p| p.map(|v| v as f64)).collect(),
                source.degrees,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    for analysis in &space.blocks {
        let source = *source_blocks
            .get(analysis.id.as_str())
            .ok_or("analysis-space block has no source geometry")?;
        let values = elevate_integer_controls(source)?;
        if values != analysis.controls {
            return Err("analysis-space elevation differs from sensitivity elevation");
        }
    }
    Ok(Model {
        space,
        radius_x_jet,
    })
}

fn fixed_lifts(space: &AnalysisSpace) -> Result<[Lift; 3], &'static str> {
    let fixed = space
        .boundary_dofs
        .get("y50")
        .ok_or("analysis space has no y50 boundary")?
        .iter()
        .map(|&dof| (dof, Ratio::zero()))
        .collect::<Vec<_>>();
    let component = || resolve(space.ndof, &[], &fixed);
    Ok([component()?, component()?, component()?])
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn face_axes(axis: usize) -> [usize; 2] {
    match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}

fn face_control(block: &AnalysisBlock, axis: usize, upper: bool, a: usize, b: usize) -> [i128; 4] {
    let free = face_axes(axis);
    let mut index = [0usize; 3];
    index[axis] = if upper { 2 } else { 0 };
    index[free[0]] = a;
    index[free[1]] = b;
    block.controls[(index[0] * 3 + index[1]) * 3 + index[2]]
}

fn physical(control: [i128; 4]) -> [f64; 3] {
    std::array::from_fn(|axis| control[axis] as f64 / control[3] as f64)
}

fn tip_y_displacement(space: &AnalysisSpace, displacement: &[f64]) -> Result<f64, &'static str> {
    let target = [60.0, 4.0, 20.0];
    let mut samples = Vec::new();
    for block in &space.blocks {
        for axis in 0..3 {
            for upper in [false, true] {
                let face: Vec<_> = (0..=2)
                    .flat_map(|a| (0..=2).map(move |b| face_control(block, axis, upper, a, b)))
                    .collect();
                if face.iter().any(|h| h[0] != 60 * h[3]) {
                    continue;
                }
                let p00 = physical(face_control(block, axis, upper, 0, 0));
                let p20 = physical(face_control(block, axis, upper, 2, 0));
                let p02 = physical(face_control(block, axis, upper, 0, 2));
                let free = face_axes(axis);
                let a = [p20[1] - p00[1], p02[1] - p00[1]];
                let b = [p20[2] - p00[2], p02[2] - p00[2]];
                let determinant = a[0] * b[1] - a[1] * b[0];
                if determinant.abs() <= 1e-20 {
                    continue;
                }
                let dy = target[1] - p00[1];
                let dz = target[2] - p00[2];
                let mut uvw = [0.0; 3];
                uvw[axis] = f64::from(u8::from(upper));
                uvw[free[0]] = (dy * b[1] - a[1] * dz) / determinant;
                uvw[free[1]] = (a[0] * dz - dy * b[0]) / determinant;
                if !(0.0..=1.0).contains(&uvw[free[0]]) || !(0.0..=1.0).contains(&uvw[free[1]]) {
                    continue;
                }
                let sample = basis_geometry(block, uvw)?;
                if (0..3)
                    .any(|component| (sample.position[component] - target[component]).abs() > 1e-8)
                {
                    continue;
                }
                let value = block
                    .dofs
                    .iter()
                    .enumerate()
                    .map(|(local, &dof)| sample.values[local] * displacement[3 * dof + 1])
                    .sum::<f64>();
                samples.push(value);
            }
        }
    }
    let first = *samples.first().ok_or("tip point is not on x60 boundary")?;
    if samples.iter().any(|value| (value - first).abs() > 1e-8) {
        return Err("tip interpolation disagrees between x60 faces");
    }
    Ok(first)
}

struct SolvedModel {
    compliance: f64,
    mean_y: f64,
    residual: f64,
    displacement: Vec<f64>,
    stiffness: Stiffness,
}

fn solve_model(candidate: &Model, order: usize) -> Result<SolvedModel, &'static str> {
    let stiffness = assemble(&candidate.space, YOUNG_MPA, POISSON, order)?;
    let (force, area) = bracket_load(&candidate.space, order)?;
    if (area - 320.0).abs() > 1e-8 {
        return Err("x60 load area changed under the radius edit");
    }
    let lifts = fixed_lifts(&candidate.space)?;
    let solution = solve_constrained(&stiffness, lifts.each_ref(), &force, 1e-10, 20_000)?;
    let compliance = dot(&force, &solution.displacement_mm);
    Ok(SolvedModel {
        compliance,
        mean_y: compliance / LOAD_RESULTANT_N,
        residual: solution.relative_residual,
        displacement: solution.displacement_mm,
        stiffness,
    })
}

fn apply_and_dot(stiffness: &Stiffness, left: &[f64], right: &[f64]) -> Result<f64, &'static str> {
    Ok(dot(left, &stiffness.apply(right)?))
}

fn stiffness_relative_error(
    plus: &Stiffness,
    minus: &Stiffness,
    rate: &Stiffness,
    h: f64,
) -> Result<f64, &'static str> {
    if plus.offsets != minus.offsets
        || plus.offsets != rate.offsets
        || plus.columns != minus.columns
        || plus.columns != rate.columns
    {
        return Err("stiffness derivative sparsity changed");
    }
    let mut error_sq = 0.0;
    let mut reference_sq = 0.0;
    for ((p, m), d) in plus.values.iter().zip(&minus.values).zip(&rate.values) {
        for row in 0..3 {
            for column in 0..3 {
                let fd = (p[row][column] - m[row][column]) / (2.0 * h);
                error_sq += (fd - d[row][column]).powi(2);
                reference_sq += d[row][column].powi(2);
            }
        }
    }
    Ok(error_sq.sqrt() / reference_sq.sqrt().max(1e-30))
}

fn objective_sensitivity(
    model: &Model,
    force: &[f64],
    stiffness: &Stiffness,
    stiffness_rate: &Stiffness,
    lifts: &[Lift; 3],
    gradient: &[f64],
) -> Result<SensitivityResult, &'static str> {
    let lift_rates: [Vec<LiftRate>; 3] =
        std::array::from_fn(|_| vec![LiftRate::default(); model.space.ndof]);
    let zero_force_rate = vec![0.0; force.len()];
    differentiate(
        MechanicalDerivative {
            stiffness,
            stiffness_rate,
            force,
            force_rate: &zero_force_rate,
            lifts: lifts.each_ref(),
            lift_rates: lift_rates.each_ref().map(Vec::as_slice),
            objective_gradient: gradient,
            objective_partial: 0.0,
            fixed_regime: true,
        },
        1e-10,
        20_000,
    )
}

fn run_order(order: usize) -> Result<String, &'static str> {
    let total = Instant::now();
    let center = model(BASE_RADIUS_TICKS)?;
    let sensitivity = assemble_sensitivity(
        &center.space,
        &center.radius_x_jet,
        MaterialJet {
            young_mpa: YOUNG_MPA,
            poisson: POISSON,
            d_young_mpa: 0.0,
            d_poisson: 0.0,
        },
        order,
    )?;
    let (force, area) = bracket_load(&center.space, order)?;
    let lifts = fixed_lifts(&center.space)?;
    let compliance_result = objective_sensitivity(
        &center,
        &force,
        &sensitivity.stiffness,
        &sensitivity.derivative,
        &lifts,
        &force,
    )?;
    let mean_gradient: Vec<_> = force.iter().map(|value| value / LOAD_RESULTANT_N).collect();
    let mean_result = objective_sensitivity(
        &center,
        &force,
        &sensitivity.stiffness,
        &sensitivity.derivative,
        &lifts,
        &mean_gradient,
    )?;
    let compliance = dot(&force, &compliance_result.primal.displacement_mm);
    let mean_y = compliance / LOAD_RESULTANT_N;
    let compliance_explicit = -apply_and_dot(
        &sensitivity.derivative,
        &compliance_result.primal.displacement_mm,
        &compliance_result.primal.displacement_mm,
    )?;
    let tangent_adjoint_error =
        (compliance_result.objective_tangent - compliance_result.objective_adjoint).abs()
            / compliance_result.objective_tangent.abs().max(1e-30);
    let mean_tangent_adjoint_error =
        (mean_result.objective_tangent - mean_result.objective_adjoint).abs()
            / mean_result.objective_tangent.abs().max(1e-30);
    if tangent_adjoint_error > 1e-7
        || mean_tangent_adjoint_error > 1e-7
        || (area - 320.0).abs() > 1e-8
    {
        return Err("compliance tangent/adjoint or load-area check failed");
    }

    let fillet = fillet_limits(&center.space, &compliance_result.primal.displacement_mm)?;
    let edited = model(56)?;
    let edited_solution = solve_model(&edited, order)?;
    let edited_fillet = fillet_limits(&edited.space, &edited_solution.displacement)?;
    let edited_tip_y = tip_y_displacement(&edited.space, &edited_solution.displacement)?;
    let mut fd = Vec::new();
    for offset in RADIUS_OFFSETS_TICKS {
        let h = f64::from(offset) / 16.0;
        let plus = model(BASE_RADIUS_TICKS + offset)?;
        let minus = model(BASE_RADIUS_TICKS - offset)?;
        let plus_solution = solve_model(&plus, order)?;
        let minus_solution = solve_model(&minus, order)?;
        let plus_compliance = plus_solution.compliance;
        let minus_compliance = minus_solution.compliance;
        let plus_mean = plus_solution.mean_y;
        let minus_mean = minus_solution.mean_y;
        let centered_compliance = (plus_compliance - minus_compliance) / (2.0 * h);
        let centered_mean = (plus_mean - minus_mean) / (2.0 * h);
        let kprime_error = stiffness_relative_error(
            &plus_solution.stiffness,
            &minus_solution.stiffness,
            &sensitivity.derivative,
            h,
        )?;
        fd.push((
            h,
            plus_compliance,
            minus_compliance,
            centered_compliance,
            centered_mean,
            plus_solution.residual.max(minus_solution.residual),
            kprime_error,
        ));
    }
    let taylor_remainders: Vec<_> = fd
        .iter()
        .map(|record| {
            let h = record.0;
            let rp = (record.1 - compliance - h * compliance_result.objective_adjoint).abs();
            let rm = (record.2 - compliance + h * compliance_result.objective_adjoint).abs();
            rp.max(rm)
        })
        .collect();
    let taylor_ratios: Vec<_> = taylor_remainders
        .windows(2)
        .map(|pair| pair[0] / pair[1].max(1e-300))
        .collect();
    let order_two_divisions = taylor_ratios.iter().filter(|ratio| **ratio >= 3.5).count();
    let records = fd.iter().enumerate().map(|(index, record)| {
        format!(
            "{{\"h_mm\":{:.17e},\"plus_compliance_n_mm\":{:.17e},\"minus_compliance_n_mm\":{:.17e},\"centered_compliance_derivative\":{:.17e},\"centered_mean_y_derivative_mm_per_mm\":{:.17e},\"compliance_derivative_error\":{:.17e},\"mean_y_derivative_error\":{:.17e},\"taylor_remainder_n_mm\":{:.17e},\"taylor_ratio_to_next\":{},\"max_fd_residual\":{:.17e},\"kprime_relative_error\":{:.17e}}}",
            record.0, record.1, record.2, record.3, record.4,
            (record.3 - compliance_result.objective_adjoint).abs(),
            (record.4 - mean_result.objective_adjoint).abs(),
            taylor_remainders[index],
            taylor_ratios.get(index).map(|ratio| format!("{ratio:.17e}")).unwrap_or_else(|| "null".to_owned()),
            record.5, record.6,
        )
    }).collect::<Vec<_>>().join(",");
    let fillet_json = fillet
        .iter()
        .map(|(id, value)| format!("\"{}\":{:.17e}", id, value))
        .collect::<Vec<_>>()
        .join(",");
    let edited_fillet_json = edited_fillet
        .iter()
        .map(|(id, value)| format!("\"{}\":{:.17e}", id, value))
        .collect::<Vec<_>>()
        .join(",");
    Ok(format!(
        "{{\"quadrature_order\":{order},\"radius_x_mm\":4.0,\"radius_y_mm\":4.0,\"scalar_dofs\":{},\"blocks\":{},\"load_area_mm2\":{area:.17e},\"compliance_n_mm\":{compliance:.17e},\"compliance_derivative_n_mm_per_mm\":{:.17e},\"compliance_derivative_explicit_minus_uKprimeu\":{compliance_explicit:.17e},\"mean_y_displacement_mm\":{mean_y:.17e},\"mean_y_derivative_mm_per_mm\":{:.17e},\"tangent_adjoint_relative_error\":{tangent_adjoint_error:.17e},\"mean_tangent_adjoint_relative_error\":{mean_tangent_adjoint_error:.17e},\"tangent_residual\":{:.17e},\"adjoint_residual\":{:.17e},\"fillet_limits_mpa\":{{{fillet_json}}},\"edited_radius_x_mm\":3.5,\"edited_mean_y_displacement_mm\":{:.17e},\"edited_tip_y_displacement_mm\":{edited_tip_y:.17e},\"edited_fillet_limits_mpa\":{{{edited_fillet_json}}},\"edited_solver_residual\":{:.17e},\"fd_cases\":[{records}],\"quadratic_taylor_divisions\":{order_two_divisions},\"total_wall_seconds\":{:.9},\"qualification\":{{\"at_least_two_order_two_divisions\":{qualified}}}}}",
        center.space.ndof, center.space.blocks.len(), compliance_result.objective_adjoint,
        mean_result.objective_adjoint, compliance_result.tangent_relative_residual,
        compliance_result.adjoint_relative_residual, edited_solution.mean_y,
        edited_solution.residual, total.elapsed().as_secs_f64(),
        qualified = order_two_divisions >= 2,
    ))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let results = [run_order(6)?, run_order(8)?];
    println!("{{\"family\":\"parameterized_bracket_v1\",\"radius_parameter\":\"x bore radius in millimetres\",\"base_radius_mm\":4.0,\"reference_y_radius_mm\":4.0,\"radius_steps_mm\":[0.5,0.25,0.125,0.0625],\"material\":{{\"young_mpa\":{YOUNG_MPA},\"poisson\":{POISSON}}},\"results\":[{}]}}", results.join(","));
    Ok(())
}
