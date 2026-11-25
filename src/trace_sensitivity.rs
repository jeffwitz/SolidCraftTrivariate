//! First derivatives of normalized rational field transfer on an audited subface.

use crate::face_constraints::{nested_constraints, FaceNet, Span};
use crate::incidence::Orientation;

#[derive(Debug, Clone, PartialEq)]
pub struct CoefficientJet {
    pub master: usize,
    pub value: f64,
    pub derivative: f64,
}

/// Differentiate the scalar transfer for a fixed orientation and degree pair.
///
/// Weight derivatives use the same homogeneous scale as their respective net.
/// An independent time-dependent global rescaling of either net leaves field
/// transfer invariant. Position-jet compatibility and differentiability of the
/// selected charts are separate preconditions for a moving geometry problem.
pub fn nested_weight_jets(
    master: &FaceNet,
    master_dw: &[f64],
    slave: &FaceNet,
    slave_dw: &[f64],
    spans: [Span; 2],
    orientation: Orientation,
) -> Result<Vec<Vec<CoefficientJet>>, &'static str> {
    let rows = nested_constraints(master, slave, spans, orientation)?;
    if master_dw.len() != master.controls.len() || slave_dw.len() != slave.controls.len() {
        return Err("weight derivative count mismatch");
    }
    if master_dw.iter().chain(slave_dw).any(|v| !v.is_finite()) {
        return Err("nonfinite weight derivative");
    }
    let log_derivative: Vec<_> = master_dw
        .iter()
        .zip(&master.controls)
        .map(|(dw, p)| dw / p[3] as f64)
        .collect();
    let mut gauge = None;
    let mut result = Vec::with_capacity(rows.len());
    for (index, row) in rows.into_iter().enumerate() {
        let value = |numerator: i128| numerator as f64 / row.denominator as f64;
        let mean: f64 = row
            .terms
            .iter()
            .map(|&(j, n)| value(n) * log_derivative[j])
            .sum();
        let observed = slave_dw[index] / slave.controls[index][3] as f64;
        let scale_derivative = observed - mean;
        // The only admissible difference is a single logarithmic scale rate.
        // This numerical check bounds roundoff of at most nine weighted terms.
        let expected = *gauge.get_or_insert(scale_derivative);
        let max_rate = row
            .terms
            .iter()
            .map(|&(j, _)| log_derivative[j].abs())
            .fold(observed.abs().max(expected.abs()).max(1.), f64::max);
        if !scale_derivative.is_finite()
            || (scale_derivative - expected).abs() > 64. * f64::EPSILON * max_rate
        {
            return Err("incompatible restricted weight derivatives");
        }
        let mut terms = Vec::with_capacity(row.terms.len());
        for (master, numerator) in row.terms {
            let coefficient = value(numerator);
            let derivative = coefficient * (log_derivative[master] - mean);
            if !derivative.is_finite() {
                return Err("nonfinite field coefficient derivative");
            }
            terms.push(CoefficientJet {
                master,
                value: coefficient,
                derivative,
            });
        }
        result.push(terms);
    }
    Ok(result)
}
