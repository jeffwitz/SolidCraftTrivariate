//! Constrained C1 reparameterization for two adjacent unit-width bricks.
//!
//! The map on brick `i` is `x_i(u)=i+u+a_i*u*(1-u)`. The equality
//! `a_0+a_1=0` matches the interface speed, while `|a_i|<=3/4` keeps the
//! derivative bounded below by `1/4`. For fixed data weights and regularizer,
//! the exact L2 objective reduces to a scalar strictly convex quadratic except
//! when both data and regularization vanish.

const BOUND: f64 = 0.75;

/// Selection regime of the constrained minimizer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActiveRegime {
    /// Unique minimizer strictly inside both speed bounds.
    Interior,
    /// Unique minimizer at the positive bound.
    UpperBound,
    /// Unique minimizer at the negative bound.
    LowerBound,
    /// The unconstrained solution meets a bound; its two-sided derivative is unavailable.
    ActivationThreshold,
    /// Floating-point uncertainty overlaps a bound, so sensitivity is not certified.
    NumericalUncertainty,
    /// The objective is constant on the complete feasible interval.
    NonUnique,
}

/// KKT diagnostics for the equality-reduced quadratic program.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KktDiagnostics {
    /// Equality multiplier for `a0+a1=0`.
    pub equality_multiplier: f64,
    /// Multipliers of `a0<=3/4` and `-a0<=3/4`.
    pub upper_bound_multiplier: f64,
    pub lower_bound_multiplier: f64,
    /// Maximum absolute stationarity, feasibility, duality, and complementarity residual.
    pub residual: f64,
    pub equality_residual: f64,
    pub primal_bound_violation: f64,
    pub dual_violation: f64,
    pub complementarity_residual: f64,
    /// Rank of the single equality row.
    pub equality_rank: usize,
    /// Second derivative of the objective along parameter `a0`, `a1=-a0`.
    pub tangent_curvature: f64,
}

/// Solution family and sensitivity for the two-brick map.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReparameterizationResult {
    /// The selected coefficients when the minimizer is unique; otherwise `None`.
    pub coefficients: Option<[f64; 2]>,
    /// Projection of the complete minimizer set onto `a0`.
    pub a0_range: [f64; 2],
    /// Projection of the complete minimizer set onto `a1`; values remain correlated by `a1=-a0`.
    pub a1_range: [f64; 2],
    pub unique: bool,
    pub regime: ActiveRegime,
    pub objective: f64,
    pub kkt: KktDiagnostics,
    /// Derivative of `[a0,a1]` with respect to `theta`, if unique and two-sided.
    pub coefficient_derivative_wrt_theta: Option<[f64; 2]>,
    /// Positive activation threshold in `theta`, when the data quadratic is nonzero.
    pub activation_threshold: Option<f64>,
}

impl ReparameterizationResult {
    /// Evaluate the physical map of one brick at a closed-domain parameter.
    pub fn map(&self, brick: usize, u: f64) -> Result<f64, &'static str> {
        if brick > 1 || !u.is_finite() || !(0.0..=1.0).contains(&u) {
            return Err("invalid brick or parameter");
        }
        let coefficients = self.coefficients.ok_or("map is set-valued")?;
        Ok(brick as f64 + u + coefficients[brick] * u * (1.0 - u))
    }

    /// Evaluate the map derivative with respect to `u`.
    pub fn map_derivative(&self, brick: usize, u: f64) -> Result<f64, &'static str> {
        if brick > 1 || !u.is_finite() || !(0.0..=1.0).contains(&u) {
            return Err("invalid brick or parameter");
        }
        let coefficients = self.coefficients.ok_or("map is set-valued")?;
        Ok(1.0 + coefficients[brick] * (1.0 - 2.0 * u))
    }

    /// Evaluate `∂x_i/∂theta`; unavailable at a nonunique point or activation threshold.
    pub fn map_derivative_wrt_theta(&self, brick: usize, u: f64) -> Result<f64, &'static str> {
        if brick > 1 || !u.is_finite() || !(0.0..=1.0).contains(&u) {
            return Err("invalid brick or parameter");
        }
        let derivative = self
            .coefficient_derivative_wrt_theta
            .ok_or("map derivative is unavailable in this regime")?;
        Ok(derivative[brick] * u * (1.0 - u))
    }
}

/// Minimize the weighted L2 mismatch to `1±theta(1−2u)` plus quadratic
/// regularization, with C1 interface matching and the fixed speed bound.
///
/// The objective convention is
/// `Σ_i data_weights[i] * ∫(x_i'(u)-target_i(u))² du + regularization*(a0²+a1²)`.
/// The profiles use `target_0=1+theta*(1-2u)` and
/// `target_1=1-theta*(1-2u)`. No clipping is applied to `theta`; the map
/// constraint is enforced exactly through the active scalar bound.
pub fn optimize_reparameterization(
    theta: f64,
    data_weights: [f64; 2],
    regularization: f64,
) -> Result<ReparameterizationResult, &'static str> {
    if !theta.is_finite()
        || data_weights.iter().any(|w| !w.is_finite() || *w < 0.)
        || !regularization.is_finite()
        || regularization < 0.
    {
        return Err("invalid reparameterization objective data");
    }

    let alpha = [data_weights[0] / 3.0, data_weights[1] / 3.0];
    if (0..2).any(|i| data_weights[i] > 0.0 && alpha[i] == 0.0) {
        return Err("data weight is below the resolvable objective range");
    }
    let data_curvature = alpha[0] + alpha[1];
    if data_weights.iter().any(|w| *w > 0.0) && data_curvature == 0.0 {
        return Err("data curvature is below the resolvable objective range");
    }
    let denominator = data_curvature + 2.0 * regularization;
    if !denominator.is_finite() {
        return Err("objective coefficient overflow");
    }

    if denominator == 0. {
        let coefficients = None;
        let residual = 0.;
        return Ok(ReparameterizationResult {
            coefficients,
            a0_range: [-BOUND, BOUND],
            a1_range: [-BOUND, BOUND],
            unique: false,
            regime: ActiveRegime::NonUnique,
            objective: 0.,
            kkt: KktDiagnostics {
                equality_multiplier: 0.,
                upper_bound_multiplier: 0.,
                lower_bound_multiplier: 0.,
                residual,
                equality_residual: 0.,
                primal_bound_violation: 0.,
                dual_violation: 0.,
                complementarity_residual: 0.,
                equality_rank: 1,
                tangent_curvature: 0.,
            },
            coefficient_derivative_wrt_theta: None,
            activation_threshold: None,
        });
    }

    let response = data_curvature / denominator;
    if data_curvature > 0.0 && response == 0.0 {
        return Err("data response is below the resolvable objective range");
    }
    let unconstrained = response * theta;
    if !response.is_finite() || !unconstrained.is_finite() {
        return Err("unconstrained coefficient overflow");
    }
    let activation_threshold = if data_curvature > 0. {
        let threshold = BOUND / response;
        if !threshold.is_finite() {
            return Err("activation threshold exceeds finite floating-point range");
        }
        Some(threshold)
    } else {
        None
    };
    let threshold_tolerance = 32. * f64::EPSILON * BOUND.max(unconstrained.abs()).max(1.);
    let distance_to_bound = (unconstrained.abs() - BOUND).abs();
    let at_threshold = unconstrained.abs() == BOUND;
    let uncertain = !at_threshold && distance_to_bound <= threshold_tolerance;
    let (a0, regime) = if at_threshold {
        (unconstrained, ActiveRegime::ActivationThreshold)
    } else if unconstrained > BOUND {
        (
            BOUND,
            if uncertain {
                ActiveRegime::NumericalUncertainty
            } else {
                ActiveRegime::UpperBound
            },
        )
    } else if unconstrained < -BOUND {
        (
            -BOUND,
            if uncertain {
                ActiveRegime::NumericalUncertainty
            } else {
                ActiveRegime::LowerBound
            },
        )
    } else {
        (
            unconstrained,
            if uncertain {
                ActiveRegime::NumericalUncertainty
            } else {
                ActiveRegime::Interior
            },
        )
    };
    let coefficients = [a0, -a0];
    let objective = alpha[0] * (a0 - theta).powi(2)
        + alpha[1] * (-a0 + theta).powi(2)
        + regularization * (a0 * a0 + coefficients[1] * coefficients[1]);
    if !objective.is_finite() {
        return Err("objective value overflow");
    }

    let grad0 = 2. * alpha[0] * (a0 - theta) + 2. * regularization * a0;
    let grad1 = 2. * alpha[1] * (-a0 + theta) - 2. * regularization * a0;
    let equality_multiplier = -grad1;
    let reduced_gradient = grad0 + equality_multiplier;
    let upper_multiplier = if a0 == BOUND && unconstrained > BOUND {
        (-reduced_gradient).max(0.)
    } else {
        0.
    };
    let lower_multiplier = if a0 == -BOUND && unconstrained < -BOUND {
        reduced_gradient.max(0.)
    } else {
        0.
    };
    let stationarity0 = grad0 + equality_multiplier + upper_multiplier - lower_multiplier;
    let stationarity1 = grad1 + equality_multiplier;
    let equality_residual = (coefficients[0] + coefficients[1]).abs();
    let primal_bound_violation = (coefficients[0].abs() - BOUND).max(0.);
    let dual_violation = (-upper_multiplier).max(0.).max((-lower_multiplier).max(0.));
    let complementarity_residual = (upper_multiplier * (BOUND - a0))
        .abs()
        .max((lower_multiplier * (BOUND + a0)).abs());
    let residual = stationarity0
        .abs()
        .max(stationarity1.abs())
        .max(equality_residual)
        .max(primal_bound_violation)
        .max(dual_violation)
        .max(complementarity_residual);
    let tangent_curvature = 2. * denominator;
    if [
        grad0,
        grad1,
        equality_multiplier,
        upper_multiplier,
        lower_multiplier,
        residual,
        tangent_curvature,
    ]
    .iter()
    .any(|x| !x.is_finite())
    {
        return Err("nonfinite KKT quantities");
    }

    let unique = true;
    let coefficient_derivative_wrt_theta = match regime {
        ActiveRegime::Interior => Some([response, -response]),
        ActiveRegime::UpperBound | ActiveRegime::LowerBound => Some([0., 0.]),
        ActiveRegime::ActivationThreshold
        | ActiveRegime::NumericalUncertainty
        | ActiveRegime::NonUnique => None,
    };
    Ok(ReparameterizationResult {
        coefficients: Some(coefficients),
        a0_range: [a0, a0],
        a1_range: [-a0, -a0],
        unique,
        regime,
        objective,
        kkt: KktDiagnostics {
            equality_multiplier,
            upper_bound_multiplier: upper_multiplier,
            lower_bound_multiplier: lower_multiplier,
            residual,
            equality_residual,
            primal_bound_violation,
            dual_violation,
            complementarity_residual,
            equality_rank: 1,
            tangent_curvature,
        },
        coefficient_derivative_wrt_theta,
        activation_threshold,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interior_solution_satisfies_c1_kkt_and_domain_bounds() {
        let result = optimize_reparameterization(0.2, [1., 2.], 0.1).unwrap();
        assert_eq!(result.regime, ActiveRegime::Interior);
        assert!(result.unique);
        assert!((result.coefficients.unwrap()[0] - 1. / 6.).abs() < 1e-15);
        assert!((result.coefficients.unwrap()[0] + result.coefficients.unwrap()[1]).abs() < 1e-15);
        assert!(result.kkt.residual <= 1e-12);
        assert_eq!(result.kkt.equality_rank, 1);
        assert!(result.kkt.tangent_curvature > 0.);
        assert_eq!(result.map(0, 0.).unwrap(), 0.);
        assert_eq!(result.map(1, 1.).unwrap(), 2.);
        for brick in 0..2 {
            for u in [0., 0.25, 0.5, 0.75, 1.] {
                assert!(result.map_derivative(brick, u).unwrap() >= 0.25 - 1e-15);
            }
        }
        assert!((result.map_derivative_wrt_theta(0, 0.5).unwrap() - 5. / 24.).abs() < 1e-15);
    }

    #[test]
    fn bound_activation_is_continuous_but_has_no_two_sided_map_derivative() {
        let below = optimize_reparameterization(0.75 - 1e-6, [1., 1.], 0.).unwrap();
        let at = optimize_reparameterization(0.75, [1., 1.], 0.).unwrap();
        let above = optimize_reparameterization(0.75 + 1e-6, [1., 1.], 0.).unwrap();
        let negative_at = optimize_reparameterization(-0.75, [1., 1.], 0.).unwrap();
        assert_eq!(below.regime, ActiveRegime::Interior);
        assert_eq!(at.regime, ActiveRegime::ActivationThreshold);
        assert_eq!(above.regime, ActiveRegime::UpperBound);
        assert_eq!(negative_at.regime, ActiveRegime::ActivationThreshold);
        assert_eq!(negative_at.coefficients.unwrap()[0], -0.75);
        assert!(negative_at.coefficient_derivative_wrt_theta.is_none());
        assert!(below.coefficients.unwrap()[0] < 0.75);
        assert_eq!(at.coefficients.unwrap()[0], 0.75);
        assert_eq!(above.coefficients.unwrap()[0], 0.75);
        assert!(at.coefficient_derivative_wrt_theta.is_none());
        assert!(above.coefficient_derivative_wrt_theta.is_some());
        assert!(at.kkt.residual <= 1e-12);
        assert!(above.kkt.residual <= 1e-12);
        assert!((below.map(0, 0.5).unwrap() - at.map(0, 0.5).unwrap()).abs() < 1e-6);
        assert!((at.map(0, 0.5).unwrap() - above.map(0, 0.5).unwrap()).abs() < 1e-15);
        assert!((below.map_derivative_wrt_theta(0, 0.5).unwrap() - 0.25).abs() < 1e-15);
        assert_eq!(above.map_derivative_wrt_theta(0, 0.5).unwrap(), 0.);
    }

    #[test]
    fn no_data_and_no_regularization_returns_the_full_feasible_minimizer_set() {
        let result = optimize_reparameterization(0.4, [0., 0.], 0.).unwrap();
        assert_eq!(result.regime, ActiveRegime::NonUnique);
        assert!(!result.unique);
        assert!(result.coefficients.is_none());
        assert_eq!(result.a0_range, [-0.75, 0.75]);
        assert_eq!(result.a1_range, [-0.75, 0.75]);
        assert!(result.coefficient_derivative_wrt_theta.is_none());
        assert_eq!(result.kkt.tangent_curvature, 0.);
        assert_eq!(result.kkt.residual, 0.);
        assert!(result.map(0, 0.5).is_err());
    }

    #[test]
    fn regularization_alone_selects_zero_and_bad_inputs_are_refused() {
        let regularized = optimize_reparameterization(100., [0., 0.], 1.).unwrap();
        assert_eq!(regularized.regime, ActiveRegime::Interior);
        assert_eq!(regularized.coefficients, Some([0., 0.]));
        assert_eq!(regularized.coefficient_derivative_wrt_theta, Some([0., 0.]));
        assert!(optimize_reparameterization(f64::NAN, [1., 1.], 0.).is_err());
        assert!(optimize_reparameterization(0., [-1., 1.], 0.).is_err());
        assert!(optimize_reparameterization(0., [1., 1.], f64::INFINITY).is_err());
    }

    #[test]
    fn threshold_ulp_band_preserves_interior_and_reports_uncertainty() {
        let below_theta = f64::from_bits(0.75_f64.to_bits() - 1);
        let above_theta = f64::from_bits(0.75_f64.to_bits() + 1);
        let below = optimize_reparameterization(below_theta, [1., 1.], 0.).unwrap();
        let at = optimize_reparameterization(0.75, [1., 1.], 0.).unwrap();
        let above = optimize_reparameterization(above_theta, [1., 1.], 0.).unwrap();

        assert_eq!(below.regime, ActiveRegime::NumericalUncertainty);
        assert_eq!(below.coefficients.unwrap()[0], below_theta);
        assert_eq!(below.kkt.upper_bound_multiplier, 0.);
        assert_eq!(at.regime, ActiveRegime::ActivationThreshold);
        assert_eq!(at.kkt.upper_bound_multiplier, 0.);
        assert_eq!(above.regime, ActiveRegime::NumericalUncertainty);
        assert_eq!(above.coefficients.unwrap()[0], 0.75);
        assert!(above.kkt.upper_bound_multiplier > 0.);
        for result in [below, at, above] {
            assert!(result.coefficient_derivative_wrt_theta.is_none());
            assert!(result.kkt.residual <= 1e-12);
        }
    }

    #[test]
    fn underflowing_weights_and_unrepresentable_threshold_are_refused() {
        let smallest_subnormal = f64::from_bits(1);
        assert!(optimize_reparameterization(0., [smallest_subnormal, 0.], 0.).is_err());
        assert!(optimize_reparameterization(0., [1.5e-323, 0.], 0.1).is_err());
    }
}
