//! Tangent and adjoint derivatives of a fixed-regime affine-constrained elastic system.

use crate::constrained_elasticity::{solve_constrained, ConstrainedSolution};
use crate::constraint_system::{Lift, Ratio};
use crate::elasticity::Stiffness;

/// Time/parameter derivative of one lift row, indexed like the primal component lift.
#[derive(Clone, Debug, Default)]
pub struct LiftRate {
    pub terms: Vec<(usize, f64)>,
    pub offset: f64,
}

pub struct MechanicalDerivative<'a> {
    pub stiffness: &'a Stiffness,
    pub stiffness_rate: &'a Stiffness,
    pub force: &'a [f64],
    pub force_rate: &'a [f64],
    pub lifts: [&'a Lift; 3],
    pub lift_rates: [&'a [LiftRate]; 3],
    /// Partial derivative of the scalar objective with respect to physical DOFs.
    pub objective_gradient: &'a [f64],
    /// Explicit partial derivative of the objective at fixed physical displacement.
    pub objective_partial: f64,
    pub fixed_regime: bool,
}

pub struct SensitivityResult {
    pub primal: ConstrainedSolution,
    pub physical_tangent: Vec<f64>,
    pub reduced_adjoint: Vec<f64>,
    pub objective_tangent: f64,
    pub objective_adjoint: f64,
    pub tangent_relative_residual: f64,
    pub adjoint_relative_residual: f64,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Solve and differentiate using component-specific constraints. Derivatives
/// are supplied at one parameter value. Sparse rate rows may use signed terms.
/// Root identities must be those returned by the exact constraint resolver.
pub fn differentiate(
    model: MechanicalDerivative<'_>,
    tolerance: f64,
    max_iterations: usize,
) -> Result<SensitivityResult, &'static str> {
    if !model.fixed_regime {
        return Err("mechanical derivative unavailable at regime event");
    }
    let size = model.stiffness.size();
    model.stiffness_rate.validate()?;
    if model.stiffness_rate.size() != size
        || !model.objective_partial.is_finite()
        || [model.force_rate, model.objective_gradient]
            .iter()
            .any(|v| v.len() != size || v.iter().any(|x| !x.is_finite()))
    {
        return Err("invalid mechanical derivative data");
    }
    let n = size / 3;
    let mut homogeneous = model.lifts.map(Clone::clone);
    let mut component_start = [0; 3];
    let mut reduced_count = 0;
    for (c, lift) in homogeneous.iter_mut().enumerate() {
        if lift.lines.len() != n || model.lift_rates[c].len() != n {
            return Err("invalid lift derivative dimension");
        }
        component_start[c] = reduced_count;
        reduced_count += lift.free_roots.len();
        for (j, &root) in lift.free_roots.iter().enumerate() {
            let line = lift.lines.get(root).ok_or("invalid physical root index")?;
            if line.terms.as_slice() != [(j, Ratio::one())] || line.offset != Ratio::zero() {
                return Err("free root does not have an identity lift row");
            }
        }
        for (line, rate) in lift.lines.iter_mut().zip(model.lift_rates[c]) {
            if !rate.offset.is_finite()
                || rate
                    .terms
                    .iter()
                    .any(|(j, v)| *j >= lift.free_roots.len() || !v.is_finite())
            {
                return Err("invalid lift derivative row");
            }
            line.offset = Ratio::zero();
        }
    }
    let refs = homogeneous.each_ref();
    let primal = solve_constrained(
        model.stiffness,
        model.lifts,
        model.force,
        tolerance,
        max_iterations,
    )?;
    let ku = model.stiffness.apply(&primal.displacement_mm)?;
    let dku = model.stiffness_rate.apply(&primal.displacement_mm)?;
    let mut motion = vec![0.; size];
    let mut transposed = vec![0.; reduced_count];
    for (c, &base) in component_start.iter().enumerate() {
        for (i, row) in model.lift_rates[c].iter().enumerate() {
            let id = 3 * i + c;
            motion[id] = row.offset;
            for &(j, dc) in &row.terms {
                motion[id] += dc * primal.reduced_values[base + j];
                transposed[base + j] += dc * (model.force[id] - ku[id]);
            }
        }
    }
    let kv = model.stiffness.apply(&motion)?;
    let mut rhs: Vec<_> = (0..size)
        .map(|i| model.force_rate[i] - dku[i] - kv[i])
        .collect();
    // C has an identity submatrix at its free roots, so applying C^T to
    // this physical injection adds exactly C'^T(f-Ku) to the reduced RHS.
    for (c, lift) in model.lifts.iter().enumerate() {
        for (j, &id) in lift.free_roots.iter().enumerate() {
            rhs[3 * id + c] += transposed[component_start[c] + j];
        }
    }
    let tangent = solve_constrained(model.stiffness, refs, &rhs, tolerance, max_iterations)?;
    let adjoint = solve_constrained(
        model.stiffness,
        refs,
        model.objective_gradient,
        tolerance,
        max_iterations,
    )?;
    let physical_tangent: Vec<_> = tangent
        .displacement_mm
        .iter()
        .zip(&motion)
        .map(|(a, b)| a + b)
        .collect();
    let objective_tangent =
        model.objective_partial + dot(model.objective_gradient, &physical_tangent);
    let objective_adjoint = model.objective_partial
        + dot(model.objective_gradient, &motion)
        + dot(&adjoint.displacement_mm, &rhs);
    if !objective_tangent.is_finite()
        || !objective_adjoint.is_finite()
        || physical_tangent.iter().any(|v| !v.is_finite())
    {
        return Err("nonfinite mechanical sensitivity");
    }
    Ok(SensitivityResult {
        primal,
        physical_tangent,
        reduced_adjoint: adjoint.reduced_values,
        objective_tangent,
        objective_adjoint,
        tangent_relative_residual: tangent.relative_residual,
        adjoint_relative_residual: adjoint.relative_residual,
    })
}
