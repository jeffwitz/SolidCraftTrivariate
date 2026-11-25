//! Independent bounded-fidelity certificate for quantized travelling-bore arcs.
//!
//! The serialized rational maps retain exact planar cuts and incidences, while
//! each circular arc may differ slightly from the ideal radius-2 mm circle.
//! Outward Bernstein intervals bound the full-arc radial residual; exact
//! projection and angular monotonicity plus endpoint enclosures bound the
//! Hausdorff distance, including clipped endpoints. These fidelity bounds are
//! kept separate from the positive-Jacobian and degree-one partition proof.

use crate::bracket::BracketBlock;
use crate::bracket_audit::BracketCertificate;
use crate::planar::{binomial_ratio, Interval};
use crate::travelling_bore_audit::{audit_travelling_bore_with_arc_verifier, ArcMetric, Point};
use std::convert::TryFrom;

const FIDELITY_BUDGET_MM: f64 = 1e-8;

/// Audited rational complex with explicit distances to the ideal circle family.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundedTravellingBoreCertificate {
    /// Local regularity, incidence, boundary tags and injectivity certificate.
    pub certificate: BracketCertificate,
    /// Conservative maximum deviation from the ideal radius-2 mm circle.
    pub radial_error_upper_mm: f64,
    /// Conservative two-sided Hausdorff distance to the ideal clipped arcs.
    pub hausdorff_error_upper_mm: f64,
    /// Required family-wide fidelity limit.
    pub fidelity_budget_mm: f64,
}

fn mul(a: Interval, b: Interval) -> Interval {
    a.mul(b)
}

fn bernstein_product_degree2(a: [Interval; 3], b: [Interval; 3]) -> [Interval; 5] {
    let mut result = [Interval::exact(0.0); 5];
    for i in 0..=2 {
        for j in 0..=2 {
            let factor = binomial_ratio(2, i, 2, j);
            let value = mul(mul(a[i], b[j]), factor);
            result[i + j] = result[i + j].add(value);
        }
    }
    result
}

fn next_up(value: f64) -> f64 {
    if value.is_nan() || value == f64::INFINITY {
        return value;
    }
    if value == 0.0 {
        return f64::from_bits(1);
    }
    let bits = value.to_bits();
    f64::from_bits(if value > 0.0 { bits + 1 } else { bits - 1 })
}

fn next_down(value: f64) -> f64 {
    if value.is_nan() || value == f64::NEG_INFINITY {
        return value;
    }
    if value == 0.0 {
        return -f64::from_bits(1);
    }
    let bits = value.to_bits();
    f64::from_bits(if value > 0.0 { bits - 1 } else { bits + 1 })
}

fn square_root_enclosure(value: Interval) -> Result<Interval, &'static str> {
    if !value.finite() || value.hi < 0.0 {
        return Err("bounded_bore_negative_endpoint_radicand");
    }
    let lower = value.lo.max(0.0).sqrt();
    let upper = value.hi.max(0.0).sqrt();
    let enclosure = Interval {
        lo: next_down(lower).max(0.0),
        hi: next_up(upper),
    };
    if enclosure.finite() {
        Ok(enclosure)
    } else {
        Err("bounded_bore_endpoint_sqrt_range")
    }
}

fn interval_abs_upper(value: Interval) -> Result<f64, &'static str> {
    if !value.finite() {
        return Err("bounded_bore_interval_overflow");
    }
    Ok(next_up(value.lo.abs().max(value.hi.abs())))
}

fn endpoint_error(point: Point, center_um: i128) -> Result<f64, &'static str> {
    let relative_x_numerator = point[0]
        .checked_mul(1_000)
        .and_then(|x| {
            point[2]
                .checked_mul(center_um)
                .and_then(|c| x.checked_sub(c))
        })
        .ok_or("bounded_bore_endpoint_overflow")?;
    let y_numerator = point[1]
        .checked_mul(1_000)
        .ok_or("bounded_bore_endpoint_overflow")?;
    let radius_numerator = point[2]
        .checked_mul(2_000)
        .ok_or("bounded_bore_endpoint_overflow")?;
    if y_numerator == 0
        && relative_x_numerator
            .checked_abs()
            .ok_or("bounded_bore_endpoint_overflow")?
            == radius_numerator
    {
        return Ok(0.0);
    }
    let values = point
        .iter()
        .map(|value| {
            i64::try_from(*value)
                .map(Interval::from_i64_scaled)
                .map_err(|_| "bounded_bore_control_range")
        })
        .collect::<Result<Vec<_>, _>>()?;
    let x = values[0];
    let y = values[1];
    let w = values[2];
    if w.lo <= 0.0 {
        return Err("bounded_bore_nonpositive_weight");
    }
    let center = Interval::from_i64_scaled(
        i64::try_from(center_um).map_err(|_| "bounded_bore_center_range")?,
    )
    .div(Interval::from_i64_scaled(1_000));
    let x_physical = x.div(w);
    let y_physical = y.div(w);
    let relative_x = x_physical.sub(center);
    let radius_squared = Interval::exact(4.0).sub(relative_x.mul(relative_x));
    let ideal_abs_y = square_root_enclosure(radius_squared)?;
    let actual_abs_y = if point[1] < 0 {
        y_physical.mul(Interval::exact(-1.0))
    } else {
        y_physical
    };
    interval_abs_upper(actual_abs_y.sub(ideal_abs_y))
}

/// Enclose `X²+Y²−4W²` in Bernstein form and convert it to a radial bound.
fn bounded_arc_metric(trace: &[Point], center_um: i128) -> Result<Option<ArcMetric>, &'static str> {
    if trace.len() != 3 {
        return Err("bounded_bore_arc_degree");
    }
    let mut x = [Interval::exact(0.0); 3];
    let mut y = [Interval::exact(0.0); 3];
    let mut w = [Interval::exact(0.0); 3];
    for i in 0..3 {
        let [hx, hy, hw] = trace[i];
        let hx = i64::try_from(hx).map_err(|_| "bounded_bore_control_range")?;
        let hy = i64::try_from(hy).map_err(|_| "bounded_bore_control_range")?;
        let hw = i64::try_from(hw).map_err(|_| "bounded_bore_control_range")?;
        x[i] = Interval::from_i64_scaled(hx);
        y[i] = Interval::from_i64_scaled(hy);
        w[i] = Interval::from_i64_scaled(hw);
        if hw <= 0 || w[i].lo <= 0.0 {
            return Err("bounded_bore_nonpositive_weight");
        }
    }
    let center = Interval::from_i64_scaled(
        i64::try_from(center_um).map_err(|_| "bounded_bore_center_range")?,
    )
    .div(Interval::from_i64_scaled(1_000));
    let relative_x = std::array::from_fn(|i| x[i].sub(center.mul(w[i])));
    let xx = bernstein_product_degree2(relative_x, relative_x);
    let yy = bernstein_product_degree2(y, y);
    let ww = bernstein_product_degree2(w, w);
    let residual: [Interval; 5] =
        std::array::from_fn(|i| xx[i].add(yy[i]).sub(ww[i].mul(Interval::exact(4.0))));
    let residual_abs = residual.iter().try_fold(0.0_f64, |bound, coefficient| {
        Ok::<_, &'static str>(bound.max(interval_abs_upper(*coefficient)?))
    })?;
    let minimum_weight = w
        .iter()
        .map(|weight| weight.lo)
        .fold(f64::INFINITY, f64::min);
    if minimum_weight <= 0.0 || !minimum_weight.is_finite() {
        return Err("bounded_bore_weight_enclosure");
    }
    let denominator = Interval::exact(2.0)
        .mul(Interval::exact(minimum_weight))
        .mul(Interval::exact(minimum_weight));
    if denominator.lo <= 0.0 || !denominator.finite() {
        return Err("bounded_bore_radial_denominator");
    }
    let radial = Interval {
        lo: 0.0,
        hi: residual_abs,
    }
    .div(denominator)
    .hi;
    if !radial.is_finite() || radial > FIDELITY_BUDGET_MM {
        return Err("bounded_bore_radial_fidelity_exceeded");
    }

    let endpoint = endpoint_error(trace[0], center_um)?.max(endpoint_error(trace[2], center_um)?);
    if endpoint > FIDELITY_BUDGET_MM {
        return Err("bounded_bore_endpoint_fidelity_exceeded");
    }
    // Radial projection moves every point by at most `radial`. The angle of a
    // projected endpoint differs from its ideal clipped endpoint by at most
    // the chord bound π/2 times (radial + vertical endpoint error). Combining
    // both directions and relaxing π/2 to 2 gives this conservative factor 3.
    let hausdorff = next_up(3.0 * next_up(radial + endpoint));
    if !hausdorff.is_finite() || hausdorff > FIDELITY_BUDGET_MM {
        return Err("bounded_bore_hausdorff_fidelity_exceeded");
    }
    Ok(Some(ArcMetric {
        radial_error_upper_mm: radial,
        hausdorff_error_upper_mm: hausdorff,
    }))
}

/// Independently certify a serialized quantized travelling-bore complex.
///
/// The returned local determinant, boundary partition, and injectivity proof
/// concern the represented rational complex. The two error bounds compare its
/// bore arcs to the ideal circular family and are separately limited to
/// `1e-8 mm`.
pub fn audit_travelling_bore_bounded(
    center_um: i32,
    blocks: &[BracketBlock],
) -> Result<BoundedTravellingBoreCertificate, &'static str> {
    let (certificate, radial_error_upper_mm, hausdorff_error_upper_mm) =
        audit_travelling_bore_with_arc_verifier(center_um, blocks, true, bounded_arc_metric)?;
    if radial_error_upper_mm > FIDELITY_BUDGET_MM || hausdorff_error_upper_mm > FIDELITY_BUDGET_MM {
        return Err("bounded_bore_family_fidelity_exceeded");
    }
    Ok(BoundedTravellingBoreCertificate {
        certificate,
        radial_error_upper_mm,
        hausdorff_error_upper_mm,
        fidelity_budget_mm: FIDELITY_BUDGET_MM,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::travelling_bore_bounded::travelling_bore_bounded_blocks;

    #[test]
    fn certifies_required_bounded_arcs_and_travel_regimes() {
        for center in [
            -5_000, -4_999, -4_001, -3_999, -3_600, -3_200, -2_001, -2_000, -1_999, -800, -400,
            -100, -80, -1, 1, 1_000, 3_000, 3_700, 7_999, 8_001, 8_080, 8_400, 8_800, 9_999,
            10_000, 10_001, 11_200, 11_600, 11_900, 11_920, 11_999, 12_001, 13_000, 14_000,
        ] {
            let blocks = travelling_bore_bounded_blocks(center).unwrap();
            let report = audit_travelling_bore_bounded(center, &blocks)
                .unwrap_or_else(|error| panic!("center {center}: {error}"));
            assert!(report.certificate.determinant_lower_bound_mm3 > 0.0);
            assert!(report.radial_error_upper_mm <= FIDELITY_BUDGET_MM);
            assert!(report.hausdorff_error_upper_mm <= FIDELITY_BUDGET_MM);
        }
    }

    #[test]
    fn preserves_per_block_projective_gauge() {
        let mut blocks = travelling_bore_bounded_blocks(-100).unwrap();
        for (index, block) in blocks.iter_mut().enumerate() {
            let factor = if index % 2 == 0 { 2_i64 } else { 3_i64 };
            for control in &mut block.controls {
                for coordinate in control {
                    *coordinate = coordinate.checked_mul(factor).unwrap();
                }
            }
        }
        let report = audit_travelling_bore_bounded(-100, &blocks).unwrap();
        assert!(report.hausdorff_error_upper_mm <= FIDELITY_BUDGET_MM);
    }

    #[test]
    fn rejects_events_out_of_domain_and_materially_perturbed_arcs() {
        for center in [-4_000, 0, 8_000, 12_000, -5_001, 14_001] {
            assert!(travelling_bore_bounded_blocks(center).is_err());
        }
        let mut blocks = travelling_bore_bounded_blocks(3_000).unwrap();
        let curved = blocks
            .iter_mut()
            .find(|block| {
                block.controls.iter().any(|control| {
                    control[0].abs() > 1_000_000_000_000 && control[1].abs() > 1_000_000_000_000
                })
            })
            .unwrap();
        curved.controls[5][1] += 1_000_000;
        assert!(audit_travelling_bore_bounded(3_000, &blocks).is_err());
    }
}
