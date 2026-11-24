use solidcraft_trivariate::analysis_space::{build_analysis_space, AnalysisBlock, AnalysisSpace};
use solidcraft_trivariate::bracket::bracket_blocks;
use solidcraft_trivariate::elasticity::{assemble, basis_geometry, Stiffness};
use std::collections::BTreeMap;

const E: f64 = 210_000.0;
const NU: f64 = 0.3;

fn cube_space() -> AnalysisSpace {
    let mut controls = Vec::with_capacity(27);
    for i in 0..3_i128 {
        for j in 0..3_i128 {
            for k in 0..3_i128 {
                // Physical control points are (i/2,j/2,k/2).
                controls.push([i, j, k, 2]);
            }
        }
    }
    AnalysisSpace {
        ndof: 27,
        blocks: vec![AnalysisBlock {
            id: "unit_cube".into(),
            controls,
            dofs: (0..27).collect(),
        }],
        boundary_dofs: BTreeMap::new(),
    }
}

fn control_positions(space: &AnalysisSpace) -> Vec<Option<[f64; 3]>> {
    let mut positions: Vec<Option<[f64; 3]>> = vec![None; space.ndof];
    for block in &space.blocks {
        assert_eq!(block.controls.len(), 27);
        assert_eq!(block.dofs.len(), 27);
        for (&dof, control) in block.dofs.iter().zip(&block.controls) {
            let point = std::array::from_fn(|axis| control[axis] as f64 / control[3] as f64);
            if let Some(previous) = positions[dof] {
                for axis in 0..3 {
                    assert!(
                        (previous[axis] - point[axis]).abs() < 1e-12,
                        "global control {dof} is geometrically inconsistent"
                    );
                }
            } else {
                positions[dof] = Some(point);
            }
        }
    }
    assert!(positions.iter().all(Option::is_some));
    positions
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(values: &[f64]) -> f64 {
    values.iter().map(|value| value * value).sum::<f64>().sqrt()
}

fn apply_and_energy(stiffness: &Stiffness, displacement: &[f64]) -> (Vec<f64>, f64) {
    let residual = stiffness.apply(displacement).unwrap();
    let energy = 0.5
        * displacement
            .iter()
            .zip(&residual)
            .map(|(a, b)| a * b)
            .sum::<f64>();
    (residual, energy)
}

fn affine_displacement(space: &AnalysisSpace, gradient: [[f64; 3]; 3]) -> Vec<f64> {
    let positions = control_positions(space);
    let mut displacement = vec![0.0; 3 * space.ndof];
    for (dof, point) in positions.into_iter().enumerate() {
        let point = point.unwrap();
        for component in 0..3 {
            displacement[3 * dof + component] = dot(gradient[component], point);
        }
    }
    displacement
}

fn assert_symmetric(stiffness: &Stiffness, relative_tolerance: f64) -> f64 {
    let n = stiffness.offsets.len() - 1;
    let mut max_entry: f64 = 0.0;
    let mut max_error: f64 = 0.0;
    for row in 0..n {
        for entry in stiffness.offsets[row]..stiffness.offsets[row + 1] {
            let column = stiffness.columns[entry];
            let reciprocal = (stiffness.offsets[column]..stiffness.offsets[column + 1])
                .find(|&candidate| stiffness.columns[candidate] == row)
                .expect("symmetric sparse pattern");
            for a in 0..3 {
                for b in 0..3 {
                    let left = stiffness.values[entry][a][b];
                    let right = stiffness.values[reciprocal][b][a];
                    max_entry = max_entry.max(left.abs()).max(right.abs());
                    max_error = max_error.max((left - right).abs());
                }
            }
        }
    }
    let relative = max_error / max_entry.max(f64::MIN_POSITIVE);
    assert!(
        relative <= relative_tolerance,
        "relative stiffness asymmetry {relative:e}"
    );
    relative
}

#[test]
fn cube_q2_space_has_partition_gradients_and_six_rigid_modes() {
    let space = cube_space();
    let block = &space.blocks[0];
    let coordinates: Vec<[f64; 3]> = block
        .controls
        .iter()
        .map(|p| {
            [
                p[0] as f64 / p[3] as f64,
                p[1] as f64 / p[3] as f64,
                p[2] as f64 / p[3] as f64,
            ]
        })
        .collect();

    for point in [
        [0.0, 0.0, 0.0],
        [0.17, 0.43, 0.91],
        [0.5, 0.5, 0.5],
        [1.0, 1.0, 1.0],
    ] {
        let sample = basis_geometry(block, point).unwrap();
        let sum_values: f64 = sample.values.iter().sum();
        assert!((sum_values - 1.0).abs() < 1e-12);
        for axis in 0..3 {
            let gradient_sum: f64 = sample.gradients.iter().map(|g| g[axis]).sum();
            assert!(gradient_sum.abs() < 1e-12);
        }
        for physical_axis in 0..3 {
            let reproduced: f64 = sample
                .values
                .iter()
                .zip(&coordinates)
                .map(|(value, node)| value * node[physical_axis])
                .sum();
            assert!((reproduced - sample.position[physical_axis]).abs() < 1e-12);
            for derivative_axis in 0..3 {
                let gradient: f64 = sample
                    .gradients
                    .iter()
                    .zip(&coordinates)
                    .map(|(g, node)| g[derivative_axis] * node[physical_axis])
                    .sum();
                let expected = if physical_axis == derivative_axis {
                    1.0
                } else {
                    0.0
                };
                assert!((gradient - expected).abs() < 1e-12);
            }
        }
    }

    let stiffness = assemble(&space, E, NU, 3).unwrap();
    assert!((stiffness.integrated_volume - 1.0).abs() < 1e-12);
    let positions = control_positions(&space);
    let mut max_relative_residual: f64 = 0.0;
    for mode in 0..6 {
        let mut displacement = vec![0.0; 3 * space.ndof];
        for (dof, point) in positions.iter().enumerate() {
            let [x, y, z] = point.unwrap();
            let vector = match mode {
                0 => [1.0, 0.0, 0.0],
                1 => [0.0, 1.0, 0.0],
                2 => [0.0, 0.0, 1.0],
                3 => [0.0, -z, y],
                4 => [z, 0.0, -x],
                _ => [-y, x, 0.0],
            };
            displacement[3 * dof..3 * dof + 3].copy_from_slice(&vector);
        }
        let residual = stiffness.apply(&displacement).unwrap();
        let relative = norm(&residual) / (E * stiffness.integrated_volume * norm(&displacement));
        max_relative_residual = max_relative_residual.max(relative);
        assert!(
            relative <= 1e-12,
            "rigid mode {mode} relative residual {relative:e}"
        );
    }
    eprintln!("cube six rigid modes: max relative residual {max_relative_residual:.3e}");
}

#[test]
#[ignore = "heavy bracket Q2 stiffness; run cargo test --release --test elasticity -- --ignored"]
fn bracket_affine_field_has_analytic_energy_quadrature_convergence_and_symmetric_stiffness() {
    let blocks = bracket_blocks().unwrap();
    let space = build_analysis_space(&blocks).unwrap();
    let gradient = [
        [0.01, 0.002, -0.001],
        [0.002, -0.003, 0.004],
        [-0.001, 0.004, 0.005],
    ];
    let displacement = affine_displacement(&space, gradient);
    let strain_norm_squared: f64 = gradient.iter().flatten().map(|x| x * x).sum();
    let trace = gradient[0][0] + gradient[1][1] + gradient[2][2];
    let mu = E / (2.0 * (1.0 + NU));
    let lambda = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
    let energy_density = mu * strain_norm_squared + 0.5 * lambda * trace * trace;
    let exact_volume = 34_080.0 - 616.0 * std::f64::consts::PI;
    let exact_energy = energy_density * exact_volume;

    let coarse = assemble(&space, E, NU, 8).unwrap();
    let fine = assemble(&space, E, NU, 16).unwrap();
    let (_, energy8) = apply_and_energy(&coarse, &displacement);
    let (_, energy16) = apply_and_energy(&fine, &displacement);
    let relative_energy8 = (energy8 - exact_energy).abs() / exact_energy;
    let relative_energy16 = (energy16 - exact_energy).abs() / exact_energy;
    let relative_quadrature_delta = (energy16 - energy8).abs() / exact_energy;
    let volume_error16 = (fine.integrated_volume - exact_volume).abs();
    let symmetry = assert_symmetric(&fine, 1e-12);
    eprintln!(
        "bracket Q2: volume={:.12}, volume_error={:.3e}, energy8_rel={:.3e}, energy16_rel={:.3e}, quadrature_delta_rel={:.3e}, symmetry_rel={:.3e}",
        fine.integrated_volume, volume_error16, relative_energy8, relative_energy16,
        relative_quadrature_delta, symmetry
    );
    assert!(
        volume_error16 <= 1e-8,
        "volume error {volume_error16:e} mm^3"
    );
    assert!(
        relative_energy8 <= 1e-8,
        "relative energy error {relative_energy8:e}"
    );
    assert!(
        relative_energy16 <= 1e-8,
        "relative energy error {relative_energy16:e}"
    );
    assert!(
        relative_quadrature_delta <= 1e-8,
        "quadrature variation {relative_quadrature_delta:e}"
    );
}
