//! One-sided stress recovery from physical gradients of rational basis functions.

use crate::analysis_space::{AnalysisBlock, AnalysisSpace};
use crate::elasticity::{basis_geometry, Matrix, Vector};

pub fn cauchy_stress(
    block: &AnalysisBlock,
    coefficients: &[f64],
    uvw: Vector,
    young: f64,
    poisson: f64,
) -> Result<Matrix, &'static str> {
    if !young.is_finite() || young <= 0.0 || !(-1.0 < poisson && poisson < 0.5) {
        return Err("invalid elastic material");
    }
    if block.dofs.len() != 27 || block.dofs.iter().any(|&i| i >= coefficients.len() / 3) {
        return Err("invalid displacement field");
    }
    let sample = basis_geometry(block, uvw)?;
    let mut grad = [[0.; 3]; 3];
    for (&i, g) in block.dofs.iter().zip(sample.gradients) {
        for (a, row) in grad.iter_mut().enumerate() {
            for (b, v) in row.iter_mut().enumerate() {
                *v += coefficients[3 * i + a] * g[b];
            }
        }
    }
    let trace = grad[0][0] + grad[1][1] + grad[2][2];
    let mu = young / (2. * (1. + poisson));
    let lambda = young * poisson / ((1. + poisson) * (1. - 2. * poisson));
    Ok(std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            mu * (grad[i][j] + grad[j][i]) + if i == j { lambda * trace } else { 0. }
        })
    }))
}

pub fn von_mises(s: Matrix) -> f64 {
    (0.5 * ((s[0][0] - s[1][1]).powi(2)
        + (s[1][1] - s[2][2]).powi(2)
        + (s[2][2] - s[0][0]).powi(2))
        + 3. * (s[0][1].powi(2) + s[1][2].powi(2) + s[0][2].powi(2)))
    .sqrt()
}

/// Evaluate all incident block limits at the mid-height fillet bisector.
/// C0 fields share a displacement; their stresses can differ.
pub fn fillet_limits(
    space: &AnalysisSpace,
    coefficients: &[f64],
) -> Result<Vec<(String, f64)>, &'static str> {
    let q = 14. - 6. / 2f64.sqrt();
    let mut result = Vec::new();
    for b in space.blocks.iter().filter(|b| b.id.starts_with("elbow_")) {
        for i in [0, 2] {
            for j in [0, 2] {
                let corner = |k: usize| {
                    let p = b.controls[(i * 3 + j) * 3 + k];
                    [
                        p[0] as f64 / p[3] as f64,
                        p[1] as f64 / p[3] as f64,
                        p[2] as f64 / p[3] as f64,
                    ]
                };
                let a = corner(0);
                let c = corner(2);
                if (a[0] - q).abs() > 1e-11 || (a[1] - q).abs() > 1e-11 || 20. < a[2] || 20. > c[2]
                {
                    continue;
                }
                let w = (20. - a[2]) / (c[2] - a[2]);
                let s = cauchy_stress(
                    b,
                    coefficients,
                    [i as f64 / 2., j as f64 / 2., w],
                    210000.,
                    0.3,
                )?;
                result.push((b.id.clone(), von_mises(s)));
            }
        }
    }
    if result.is_empty() {
        return Err("fillet probe outside the analysis space");
    }
    Ok(result)
}
