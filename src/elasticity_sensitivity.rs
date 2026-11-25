//! Material derivatives of the assembled isotropic elasticity operator.
//!
//! Geometry jets are derivatives of homogeneous Q2 controls in the same scale
//! as each block's controls. The differentiation is valid while degrees,
//! connectivity, material law, and quadrature order remain fixed.

use crate::analysis_space::{AnalysisBlock, AnalysisSpace};
use crate::elasticity::{gauss, Matrix, Stiffness, Vector};
use std::collections::{BTreeMap, BTreeSet};

/// Material parameters and their derivatives per unit design parameter.
#[derive(Debug, Clone, Copy)]
pub struct MaterialJet {
    pub young_mpa: f64,
    pub poisson: f64,
    pub d_young_mpa: f64,
    pub d_poisson: f64,
}

/// Geometry and rational basis data, including first design derivatives.
#[derive(Debug, Clone)]
pub struct GeometryJetSample {
    pub position_mm: Vector,
    pub position_derivative: Vector,
    pub jacobian_columns: Matrix,
    pub jacobian_column_rates: Matrix,
    pub basis: [f64; 27],
    pub basis_derivative: [f64; 27],
    pub gradients_per_mm: [Vector; 27],
    pub gradient_derivatives: [Vector; 27],
    pub determinant_mm3: f64,
    pub determinant_derivative: f64,
}

/// Primal and differentiated stiffness assembled with identical CSR structure.
pub struct ElasticitySensitivity {
    pub stiffness: Stiffness,
    pub derivative: Stiffness,
    pub integrated_volume_derivative_mm3: f64,
}

/// Evaluate the rational geometry, field basis, and their material derivatives.
pub fn evaluate_geometry_jet(
    block: &AnalysisBlock,
    control_jets: &[[f64; 4]],
    uvw: Vector,
) -> Result<GeometryJetSample, &'static str> {
    if block.controls.len() != 27
        || control_jets.len() != 27
        || uvw
            .iter()
            .any(|&x| !x.is_finite() || !(0.0..=1.0).contains(&x))
    {
        return Err("invalid Q2 geometry jet dimensions or parameters");
    }
    if control_jets.iter().flatten().any(|x| !x.is_finite()) {
        return Err("nonfinite homogeneous control jet");
    }
    if block.controls.iter().any(|p| p[3] <= 0) {
        return Err("nonpositive homogeneous weight");
    }

    let pairs = uvw.map(|x| {
        (
            [(1. - x) * (1. - x), 2. * x * (1. - x), x * x],
            [-2. * (1. - x), 2. - 4. * x, 2. * x],
        )
    });
    let basis_1d = pairs.map(|pair| pair.0);
    let dbasis_1d = pairs.map(|pair| pair.1);
    let mut h = [0.; 4];
    let mut hp = [0.; 4];
    let mut ha = [[0.; 4]; 3];
    let mut hap = [[0.; 4]; 3];
    let mut b = [0.; 27];
    let mut ba = [[0.; 27]; 3];
    for n in 0..27 {
        let i = n / 9;
        let j = n / 3 % 3;
        let k = n % 3;
        let val = basis_1d[0][i] * basis_1d[1][j] * basis_1d[2][k];
        let deriv = [
            dbasis_1d[0][i] * basis_1d[1][j] * basis_1d[2][k],
            basis_1d[0][i] * dbasis_1d[1][j] * basis_1d[2][k],
            basis_1d[0][i] * basis_1d[1][j] * dbasis_1d[2][k],
        ];
        b[n] = val;
        for a in 0..3 {
            ba[a][n] = deriv[a];
        }
        let control = block.controls[n];
        let jet = control_jets[n];
        for c in 0..4 {
            let value = control[c] as f64;
            h[c] += val * value;
            hp[c] += val * jet[c];
            for a in 0..3 {
                ha[a][c] += deriv[a] * value;
                hap[a][c] += deriv[a] * jet[c];
            }
        }
    }
    if !h[3].is_finite() || h[3] <= 0. || h.iter().chain(&hp).any(|x| !x.is_finite()) {
        return Err("invalid homogeneous denominator");
    }
    let position: Vector = std::array::from_fn(|c| h[c] / h[3]);
    let position_p: Vector = std::array::from_fn(|c| (hp[c] - position[c] * hp[3]) / h[3]);
    let mut rational = [0.; 27];
    let mut rational_p = [0.; 27];
    let mut rational_a = [[0.; 27]; 3];
    let mut rational_ap = [[0.; 27]; 3];
    let mut jac = [[0.; 3]; 3];
    let mut jac_p_numerator = [[0.; 3]; 3];
    for n in 0..27 {
        let weight = block.controls[n][3] as f64;
        let weight_p = control_jets[n][3];
        rational[n] = b[n] * weight / h[3];
        rational_p[n] = (b[n] * weight_p - rational[n] * hp[3]) / h[3];
        for a in 0..3 {
            rational_a[a][n] = (ba[a][n] * weight - rational[n] * ha[a][3]) / h[3];
            rational_ap[a][n] = (ba[a][n] * weight_p
                - rational_p[n] * ha[a][3]
                - rational[n] * hap[a][3]
                - rational_a[a][n] * hp[3])
                / h[3];
        }
    }
    for a in 0..3 {
        for c in 0..3 {
            jac[a][c] = (ha[a][c] - position[c] * ha[a][3]) / h[3];
            jac_p_numerator[a][c] =
                (hap[a][c] - position_p[c] * ha[a][3] - position[c] * hap[a][3]) / h[3];
        }
    }
    let jac_p: Matrix = std::array::from_fn(|a| {
        std::array::from_fn(|c| jac_p_numerator[a][c] - jac[a][c] * hp[3] / h[3])
    });
    if rational
        .iter()
        .chain(&rational_p)
        .chain(jac.iter().flatten())
        .chain(jac_p.iter().flatten())
        .any(|x| !x.is_finite())
    {
        return Err("nonfinite rational geometry jet");
    }
    let (det, inverse) = inverse_and_det(jac)?;
    let det_p = det
        * (0..3)
            .map(|i| (0..3).map(|a| inverse[i][a] * jac_p[a][i]).sum::<f64>())
            .sum::<f64>();
    if !det_p.is_finite() {
        return Err("nonfinite determinant derivative");
    }
    let inverse_p: Matrix = std::array::from_fn(|c| {
        std::array::from_fn(|a| {
            -(0..3)
                .flat_map(|i| (0..3).map(move |j| inverse[c][i] * jac_p[i][j] * inverse[j][a]))
                .sum::<f64>()
        })
    });
    let gradients: [Vector; 27] = std::array::from_fn(|n| {
        std::array::from_fn(|c| (0..3).map(|a| inverse[c][a] * rational_a[a][n]).sum())
    });
    let gradients_p: [Vector; 27] = std::array::from_fn(|n| {
        std::array::from_fn(|c| {
            (0..3)
                .map(|a| inverse_p[c][a] * rational_a[a][n] + inverse[c][a] * rational_ap[a][n])
                .sum()
        })
    });
    if det <= 0.
        || !det.is_finite()
        || gradients
            .iter()
            .flatten()
            .chain(gradients_p.iter().flatten())
            .any(|x| !x.is_finite())
    {
        return Err("singular or nonfinite physical geometry jet");
    }
    Ok(GeometryJetSample {
        position_mm: position,
        position_derivative: position_p,
        jacobian_columns: jac,
        jacobian_column_rates: jac_p,
        basis: rational,
        basis_derivative: rational_p,
        gradients_per_mm: gradients,
        gradient_derivatives: gradients_p,
        determinant_mm3: det,
        determinant_derivative: det_p,
    })
}

/// Assemble primal `K`, its material derivative `K′`, and the volume jet.
///
/// `control_jets` is ordered exactly like `space.blocks`; each entry contains
/// 27 `(dXw,dYw,dZw,dw)` quadruples in the same homogeneous block scale.
pub fn assemble_sensitivity(
    space: &AnalysisSpace,
    control_jets: &[Vec<[f64; 4]>],
    material: MaterialJet,
    order: usize,
) -> Result<ElasticitySensitivity, &'static str> {
    if space.blocks.is_empty() || space.ndof == 0 || control_jets.len() != space.blocks.len() {
        return Err("incomplete sensitivity space or jets");
    }
    if !material.young_mpa.is_finite()
        || material.young_mpa <= 0.
        || !material.poisson.is_finite()
        || !(-1. < material.poisson && material.poisson < 0.5)
        || !material.d_young_mpa.is_finite()
        || !material.d_poisson.is_finite()
    {
        return Err("invalid material jet");
    }
    let points = gauss(order)?;
    let mu = material.young_mpa / (2. * (1. + material.poisson));
    let lambda = material.young_mpa * material.poisson
        / ((1. + material.poisson) * (1. - 2. * material.poisson));
    let mu_p = material.d_young_mpa / (2. * (1. + material.poisson))
        - material.young_mpa * material.d_poisson / (2. * (1. + material.poisson).powi(2));
    let lambda_den = (1. + material.poisson) * (1. - 2. * material.poisson);
    let lambda_den_p = (-1. - 4. * material.poisson) * material.d_poisson;
    let lambda_p = (material.d_young_mpa * material.poisson
        + material.young_mpa * material.d_poisson)
        / lambda_den
        - material.young_mpa * material.poisson * lambda_den_p / lambda_den.powi(2);
    if [mu, lambda, mu_p, lambda_p].iter().any(|x| !x.is_finite()) {
        return Err("nonfinite Lamé parameter jet");
    }

    let mut rows = vec![BTreeMap::<usize, Matrix>::new(); space.ndof];
    let mut derivative_rows = vec![BTreeMap::<usize, Matrix>::new(); space.ndof];
    let mut volume = 0.;
    let mut volume_p = 0.;
    let mut block_ids = BTreeSet::new();
    for (block_index, block) in space.blocks.iter().enumerate() {
        if block.controls.len() != 27
            || block.dofs.len() != 27
            || block.dofs.iter().any(|&dof| dof >= space.ndof)
            || control_jets[block_index].len() != 27
            || block.dofs.iter().copied().collect::<BTreeSet<_>>().len() != 27
            || !block_ids.insert(block.id.as_str())
        {
            return Err("invalid block, numbering, or incomplete control jet");
        }
        let mut element = vec![[[0.; 3]; 3]; 27 * 27];
        let mut element_p = vec![[[0.; 3]; 3]; 27 * 27];
        for &(u, wu) in &points {
            for &(v, wv) in &points {
                for &(w, ww) in &points {
                    let sample =
                        evaluate_geometry_jet(block, &control_jets[block_index], [u, v, w])?;
                    let quadrature_weight = wu * wv * ww;
                    let measure = sample.determinant_mm3 * quadrature_weight;
                    let measure_p = sample.determinant_derivative * quadrature_weight;
                    volume += measure;
                    volume_p += measure_p;
                    for i in 0..27 {
                        for j in 0..=i {
                            let gi = sample.gradients_per_mm[i];
                            let gj = sample.gradients_per_mm[j];
                            let gip = sample.gradient_derivatives[i];
                            let gjp = sample.gradient_derivatives[j];
                            let dot_ij = dot(gi, gj);
                            let dot_ijp = dot(gip, gj) + dot(gi, gjp);
                            for a in 0..3 {
                                for b in 0..3 {
                                    let constitutive = lambda * gi[a] * gj[b]
                                        + mu * gi[b] * gj[a]
                                        + if a == b { mu * dot_ij } else { 0. };
                                    let constitutive_p = lambda_p * gi[a] * gj[b]
                                        + lambda * (gip[a] * gj[b] + gi[a] * gjp[b])
                                        + mu_p * gi[b] * gj[a]
                                        + mu * (gip[b] * gj[a] + gi[b] * gjp[a])
                                        + if a == b {
                                            mu_p * dot_ij + mu * dot_ijp
                                        } else {
                                            0.
                                        };
                                    element[i * 27 + j][a][b] += measure * constitutive;
                                    element_p[i * 27 + j][a][b] +=
                                        measure_p * constitutive + measure * constitutive_p;
                                }
                            }
                        }
                    }
                }
            }
        }
        for i in 0..27 {
            let ii = block.dofs[i];
            for j in 0..=i {
                let jj = block.dofs[j];
                accumulate(&mut rows, ii, jj, element[i * 27 + j]);
                accumulate(&mut derivative_rows, ii, jj, element_p[i * 27 + j]);
                if i != j {
                    accumulate(&mut rows, jj, ii, transpose(element[i * 27 + j]));
                    accumulate(
                        &mut derivative_rows,
                        jj,
                        ii,
                        transpose(element_p[i * 27 + j]),
                    );
                }
            }
        }
    }
    if !volume.is_finite() || !volume_p.is_finite() {
        return Err("nonfinite integrated volume jet");
    }
    let stiffness = to_stiffness(rows, volume)?;
    let derivative = to_stiffness(derivative_rows, volume_p)?;
    if stiffness.offsets != derivative.offsets || stiffness.columns != derivative.columns {
        return Err("primal and derivative CSR patterns differ");
    }
    Ok(ElasticitySensitivity {
        stiffness,
        derivative,
        integrated_volume_derivative_mm3: volume_p,
    })
}

fn dot(a: Vector, b: Vector) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn transpose(a: Matrix) -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i]))
}

fn accumulate(rows: &mut [BTreeMap<usize, Matrix>], row: usize, col: usize, value: Matrix) {
    let target = rows[row].entry(col).or_default();
    for i in 0..3 {
        for j in 0..3 {
            target[i][j] += value[i][j];
        }
    }
}

fn to_stiffness(
    rows: Vec<BTreeMap<usize, Matrix>>,
    integrated_volume: f64,
) -> Result<Stiffness, &'static str> {
    let mut offsets = vec![0];
    let mut columns = Vec::new();
    let mut values = Vec::new();
    for row in rows {
        for (column, value) in row {
            if value.iter().flatten().any(|x| !x.is_finite()) {
                return Err("nonfinite stiffness coefficient");
            }
            columns.push(column);
            values.push(value);
        }
        offsets.push(columns.len());
    }
    let stiffness = Stiffness {
        offsets,
        columns,
        values,
        integrated_volume,
    };
    stiffness.validate()?;
    Ok(stiffness)
}

fn inverse_and_det(a: Matrix) -> Result<(f64, Matrix), &'static str> {
    let cofactors: Matrix = [
        [
            a[1][1] * a[2][2] - a[1][2] * a[2][1],
            a[1][2] * a[2][0] - a[1][0] * a[2][2],
            a[1][0] * a[2][1] - a[1][1] * a[2][0],
        ],
        [
            a[0][2] * a[2][1] - a[0][1] * a[2][2],
            a[0][0] * a[2][2] - a[0][2] * a[2][0],
            a[0][1] * a[2][0] - a[0][0] * a[2][1],
        ],
        [
            a[0][1] * a[1][2] - a[0][2] * a[1][1],
            a[0][2] * a[1][0] - a[0][0] * a[1][2],
            a[0][0] * a[1][1] - a[0][1] * a[1][0],
        ],
    ];
    let det = a[0][0] * cofactors[0][0] + a[0][1] * cofactors[0][1] + a[0][2] * cofactors[0][2];
    if !det.is_finite() || det <= 0. {
        return Err("singular or inverted geometry map");
    }
    let inverse = transpose(cofactors).map(|row| row.map(|value| value / det));
    if inverse.iter().flatten().any(|x| !x.is_finite()) {
        return Err("nonfinite inverse geometry map");
    }
    Ok((det, inverse))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis_space::AnalysisBlock;

    fn cube_space(scale: i128) -> AnalysisSpace {
        let mut controls = Vec::with_capacity(27);
        for i in 0..3_i128 {
            for j in 0..3_i128 {
                for k in 0..3_i128 {
                    controls.push([i * scale, j * scale, k * scale, 2 * scale]);
                }
            }
        }
        AnalysisSpace {
            ndof: 27,
            blocks: vec![AnalysisBlock {
                id: "cube".into(),
                controls,
                dofs: (0..27).collect(),
            }],
            boundary_dofs: BTreeMap::new(),
        }
    }

    fn constant_jets(space: &AnalysisSpace, value: [f64; 4]) -> Vec<Vec<[f64; 4]>> {
        vec![vec![value; 27]; space.blocks.len()]
    }

    fn close_matrix(a: &Stiffness, b: &Stiffness, tolerance: f64) {
        assert_eq!(a.offsets, b.offsets);
        assert_eq!(a.columns, b.columns);
        let error = a
            .values
            .iter()
            .zip(&b.values)
            .flat_map(|(x, y)| {
                x.iter()
                    .flatten()
                    .zip(y.iter().flatten())
                    .map(|(u, v)| (u - v).abs())
            })
            .fold(0.0_f64, f64::max);
        let scale = a
            .values
            .iter()
            .flatten()
            .flatten()
            .map(|x| x.abs())
            .fold(0.0_f64, f64::max);
        assert!(
            error <= tolerance * scale.max(1.0),
            "matrix difference {error}, scale {scale}"
        );
    }

    #[test]
    fn physical_dilation_derivative_is_k_for_unit_rate() {
        let space = cube_space(1);
        let jets = vec![space.blocks[0]
            .controls
            .iter()
            .map(|p| [p[0] as f64, p[1] as f64, p[2] as f64, 0.])
            .collect()];
        let result = assemble_sensitivity(
            &space,
            &jets,
            MaterialJet {
                young_mpa: 210000.,
                poisson: 0.3,
                d_young_mpa: 0.,
                d_poisson: 0.,
            },
            3,
        )
        .unwrap();
        close_matrix(&result.stiffness, &result.derivative, 1e-11);
        assert!((result.integrated_volume_derivative_mm3 - 3.).abs() < 1e-12);
    }

    #[test]
    fn rigid_translation_and_global_homogeneous_scale_have_zero_derivative() {
        let space = cube_space(1);
        let translation = vec![vec![[2., -3., 4., 0.]; 27]];
        let result = assemble_sensitivity(
            &space,
            &translation,
            MaterialJet {
                young_mpa: 10.,
                poisson: 0.25,
                d_young_mpa: 0.,
                d_poisson: 0.,
            },
            3,
        )
        .unwrap();
        close_matrix(&result.derivative, &zero_like(&result.derivative), 1e-12);
        let scale = vec![space.blocks[0]
            .controls
            .iter()
            .map(|p| p.map(|x| x as f64))
            .collect()];
        let result = assemble_sensitivity(
            &space,
            &scale,
            MaterialJet {
                young_mpa: 10.,
                poisson: 0.25,
                d_young_mpa: 0.,
                d_poisson: 0.,
            },
            3,
        )
        .unwrap();
        close_matrix(&result.derivative, &zero_like(&result.derivative), 1e-12);
    }

    #[test]
    fn material_only_derivative_and_invalid_jets() {
        let space = cube_space(1);
        let zero = constant_jets(&space, [0.; 4]);
        let result = assemble_sensitivity(
            &space,
            &zero,
            MaterialJet {
                young_mpa: 100.,
                poisson: 0.3,
                d_young_mpa: 100.,
                d_poisson: 0.,
            },
            3,
        )
        .unwrap();
        close_matrix(&result.stiffness, &result.derivative, 1e-11);
        let poisson_jet = MaterialJet {
            young_mpa: 100.,
            poisson: 0.3,
            d_young_mpa: 0.,
            d_poisson: 0.2,
        };
        let result = assemble_sensitivity(&space, &zero, poisson_jet, 3).unwrap();
        let delta = 1e-5;
        let plus = crate::elasticity::assemble(
            &space,
            poisson_jet.young_mpa,
            poisson_jet.poisson + delta * poisson_jet.d_poisson,
            3,
        )
        .unwrap();
        let minus = crate::elasticity::assemble(
            &space,
            poisson_jet.young_mpa,
            poisson_jet.poisson - delta * poisson_jet.d_poisson,
            3,
        )
        .unwrap();
        let error = plus
            .values
            .iter()
            .zip(&minus.values)
            .zip(&result.derivative.values)
            .flat_map(|((p, m), d)| {
                p.iter()
                    .flatten()
                    .zip(m.iter().flatten())
                    .zip(d.iter().flatten())
                    .map(|((p, m), d)| ((p - m) / (2. * delta) - d).abs())
            })
            .fold(0.0_f64, f64::max);
        let scale = result
            .derivative
            .values
            .iter()
            .flatten()
            .flatten()
            .map(|x| x.abs())
            .fold(0.0_f64, f64::max);
        assert!(error <= 1e-7 * scale.max(1.));
        let mut incomplete = zero.clone();
        incomplete[0].pop();
        assert!(assemble_sensitivity(
            &space,
            &incomplete,
            MaterialJet {
                young_mpa: 1.,
                poisson: 0.3,
                d_young_mpa: 0.,
                d_poisson: 0.
            },
            2
        )
        .is_err());
        let mut nonfinite = zero;
        nonfinite[0][0][0] = f64::NAN;
        assert!(assemble_sensitivity(
            &space,
            &nonfinite,
            MaterialJet {
                young_mpa: 1.,
                poisson: 0.3,
                d_young_mpa: 0.,
                d_poisson: 0.
            },
            2
        )
        .is_err());
        let mut singular = space.clone();
        for control in &mut singular.blocks[0].controls {
            control[2] = 0;
        }
        assert!(assemble_sensitivity(
            &singular,
            &constant_jets(&singular, [0.; 4]),
            MaterialJet {
                young_mpa: 1.,
                poisson: 0.3,
                d_young_mpa: 0.,
                d_poisson: 0.
            },
            2
        )
        .is_err());
        let mut bad_numbering = space.clone();
        bad_numbering.blocks[0].dofs[1] = bad_numbering.blocks[0].dofs[0];
        assert!(assemble_sensitivity(
            &bad_numbering,
            &constant_jets(&bad_numbering, [0.; 4]),
            MaterialJet {
                young_mpa: 1.,
                poisson: 0.3,
                d_young_mpa: 0.,
                d_poisson: 0.
            },
            2
        )
        .is_err());
    }

    #[test]
    fn curved_variable_weight_and_height_jets_match_centered_reassembly() {
        let scale = 300_000_000_000_i128;
        let mut block_controls = Vec::with_capacity(27);
        let mut block_jets = Vec::with_capacity(27);
        for i in 0..3_i128 {
            for j in 0..3_i128 {
                for k in 0..3_i128 {
                    let sum = i + j + k;
                    let weight = 2 * scale * (100 + sum);
                    block_controls.push([
                        i * scale * (100 + sum),
                        j * scale * (100 + sum),
                        k * scale * (100 + sum),
                        weight,
                    ]);
                    block_jets.push([
                        (scale * i * j / 50) as f64,
                        (-scale * j * k / 40) as f64,
                        (scale * k * (i + j) / 80) as f64,
                        (scale * (sum - 3) / 100) as f64,
                    ]);
                }
            }
        }
        let space = AnalysisSpace {
            ndof: 27,
            blocks: vec![AnalysisBlock {
                id: "curved_weighted_cube".into(),
                controls: block_controls,
                dofs: (0..27).collect(),
            }],
            boundary_dofs: BTreeMap::new(),
        };
        let jets = vec![block_jets];
        let material = MaterialJet {
            young_mpa: 1200.,
            poisson: 0.29,
            d_young_mpa: 0.,
            d_poisson: 0.,
        };
        let order = 3;
        let derivative = assemble_sensitivity(&space, &jets, material, order).unwrap();
        let primal =
            crate::elasticity::assemble(&space, material.young_mpa, material.poisson, order)
                .unwrap();
        close_matrix(&primal, &derivative.stiffness, 1e-11);

        let delta = 1e-5;
        let shifted = |parameter: f64| {
            let mut perturbed = space.clone();
            for (control, jet) in perturbed.blocks[0].controls.iter_mut().zip(&jets[0]) {
                for axis in 0..4 {
                    let step = (jet[axis] * parameter).round() as i128;
                    control[axis] += step;
                }
            }
            perturbed
        };
        let plus = crate::elasticity::assemble(
            &shifted(delta),
            material.young_mpa,
            material.poisson,
            order,
        )
        .unwrap();
        let minus = crate::elasticity::assemble(
            &shifted(-delta),
            material.young_mpa,
            material.poisson,
            order,
        )
        .unwrap();
        assert_eq!(plus.offsets, derivative.derivative.offsets);
        assert_eq!(plus.columns, derivative.derivative.columns);
        let max_error = plus
            .values
            .iter()
            .zip(&minus.values)
            .zip(&derivative.derivative.values)
            .flat_map(|((p, m), d)| {
                p.iter()
                    .flatten()
                    .zip(m.iter().flatten())
                    .zip(d.iter().flatten())
                    .map(|((p, m), d)| ((p - m) / (2. * delta) - d).abs())
            })
            .fold(0.0_f64, f64::max);
        let max_scale = derivative
            .derivative
            .values
            .iter()
            .flatten()
            .flatten()
            .map(|x| x.abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_error <= 1e-7 * max_scale.max(1.),
            "centered K′ error {max_error}, scale {max_scale}"
        );
        let taylor_error = |assembled: &Stiffness, step: f64| {
            assembled
                .values
                .iter()
                .zip(&derivative.stiffness.values)
                .zip(&derivative.derivative.values)
                .flat_map(|((actual, base), tangent)| {
                    actual
                        .iter()
                        .flatten()
                        .zip(base.iter().flatten())
                        .zip(tangent.iter().flatten())
                        .map(move |((actual, base), tangent)| {
                            (actual - base - step * tangent).abs()
                        })
                })
                .fold(0.0_f64, f64::max)
        };
        let taylor_delta = 1e-2;
        let taylor_plus = crate::elasticity::assemble(
            &shifted(taylor_delta),
            material.young_mpa,
            material.poisson,
            order,
        )
        .unwrap();
        let taylor_half = crate::elasticity::assemble(
            &shifted(taylor_delta * 0.5),
            material.young_mpa,
            material.poisson,
            order,
        )
        .unwrap();
        let full_taylor_error = taylor_error(&taylor_plus, taylor_delta);
        let half_taylor_error = taylor_error(&taylor_half, taylor_delta * 0.5);
        let ratio = full_taylor_error / half_taylor_error;
        assert!(
            (3.5..=4.5).contains(&ratio),
            "Taylor remainder ratio {ratio}"
        );
    }

    fn zero_like(matrix: &Stiffness) -> Stiffness {
        Stiffness {
            offsets: matrix.offsets.clone(),
            columns: matrix.columns.clone(),
            values: vec![[[0.; 3]; 3]; matrix.values.len()],
            integrated_volume: 0.,
        }
    }
}
