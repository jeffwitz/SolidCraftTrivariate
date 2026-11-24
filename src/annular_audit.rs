//! Strict coefficient audit and continuous bounds for the fixed annular pair.
//!
//! The oracle reconstructs the stored binary64 controls independently of the
//! builder. Its continuous claims apply to the real rational map defined by
//! those stored values over the complete closed parameter cubes.
//!
//! Let `q=M/2^53` be the stored middle weight, `s=1-u`, and `t=u`. The planar
//! map has homogeneous denominator `D=s²+2qst+t²` and numerators
//! `X=s²+2qst`, `Y=t²+2qst`, scaled by radius `2-v`. Bernstein convexity gives
//! `q<=D<=1`. Direct expansion gives
//! `X²+Y²-D²=(4q²-2)s²t²` and
//! `XY'-YX'=2q(s²+t²)+2st>=q`. Thus its planar determinant is at least `q`.
//! The stored `q` is checked exactly from its binary mantissa: `q>7/10` and
//! `|q²-1/2|<2^-51`, the latter following from the u128 comparison
//! `|M²-2^105|<2^55`. Rationalizing the radius error then gives
//! `2|4q²-2|/(16q²)<2^-50` mm.
//!
//! At every integer axial level `n<=2000`, the exact product `q*n` and its
//! rounded binary64 value, expressed in units of `2^-53`, differ by at most
//! `2^10`; the audit checks this integer inequality for all 2001 levels using
//! u128 arithmetic. Therefore each rounded product error is at most `2^-43`.
//! Since `D>=q>1/2` and `2st<=1/2`, the plane elevation error is below
//! `2^-43` mm and the axial derivative is greater than `h-2^-42`. For every
//! allowed `h>=1`, this is greater than `3h/4`. The full determinant exceeds
//! `(7/10)*(3h/4)=21h/40>h/2`; `min_height/2` is consequently a conservative
//! lower bound in mm³. The angular coordinate is strictly increasing, the
//! radius decreases with `v`, and elevation increases with `w`; this module
//! reports local continuous bounds and does not certify other families.

use super::annular::{AnnularBlock, AnnularPair};

const Q_BITS: u64 = 0x3fe6_a09e_667f_3bcd;
const HEIGHT_LIMIT_MM: u32 = 1_000;

/// Conservative bounds for the certified stored annular pair, in millimetres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AnnularCertificate {
    /// Lower bound for the homogeneous denominator on either closed cube.
    pub min_weight: f64,
    /// Lower bound for det(DX), with parameter columns ordered (angle, radius, z).
    pub determinant_lower_bound_mm3: f64,
    /// Upper bound on radial deviation from the nominal quarter-annulus, in mm.
    pub radial_error_upper_mm: f64,
    /// Upper bound on axial deviation from the nominal planes, in mm.
    pub height_error_upper_mm: f64,
}

fn q_mantissa() -> u128 {
    (1u128 << 52) | ((Q_BITS & ((1u64 << 52) - 1)) as u128)
}

fn exact_product_error_within_bound(level: u32) -> bool {
    if level == 0 {
        return true;
    }
    let q = f64::from_bits(Q_BITS);
    let product = q * f64::from(level);
    let bits = product.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    if exponent == 0 || exponent == 0x7ff {
        return false;
    }
    // For these levels product lies in [0.7, 1415), so this shift is
    // nonnegative. Multiplication by 2^53 makes the stored value an integer.
    let shift = exponent - 1022;
    if shift < 0 {
        return false;
    }
    let product_mantissa = (1u128 << 52) | ((bits & ((1u64 << 52) - 1)) as u128);
    let stored_scaled = product_mantissa << (shift as u32);
    let exact_scaled = q_mantissa() * (level as u128);
    stored_scaled.abs_diff(exact_scaled) <= (1u128 << 10)
}

fn validate_fixed_dyadic_bounds() -> bool {
    let mantissa = q_mantissa();
    let q_squared_numerator = mantissa * mantissa;
    let half_scaled = 1u128 << 105;
    let q_square_error_numerator = q_squared_numerator.abs_diff(half_scaled);

    // q>7/10 and 2*|q^2-1/2| <= 2^-50 imply the radial error bound.
    let q_above_seven_tenths = 10 * mantissa > 7 * (1u128 << 53);
    let radial_error_below_dyadic = q_square_error_numerator < (1u128 << 55);
    let product_rounding_bounded = (0..=2_000).all(exact_product_error_within_bound);
    q_above_seven_tenths && radial_error_below_dyadic && product_rounding_bounded
}

fn same_bits(left: f64, right: f64) -> bool {
    left.to_bits() == right.to_bits()
}

fn expected_block(id: &'static str, z0: u32, height: u32) -> AnnularBlock {
    let q = f64::from_bits(Q_BITS);
    let angular = [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let angular_weights = [1.0, q, 1.0];
    let radii = [2.0, 1.0];
    let mut controls = Vec::with_capacity(12);
    for angle in 0..3 {
        for radius in radii {
            for local_z in [z0, z0 + height] {
                let weight = angular_weights[angle];
                let (cx, cy) = angular[angle];
                controls.push([
                    weight * radius * cx,
                    weight * radius * cy,
                    weight * f64::from(local_z),
                    weight,
                ]);
            }
        }
    }
    AnnularBlock {
        id,
        degrees: [2, 1, 1],
        controls,
    }
}

fn same_block(actual: &AnnularBlock, expected: &AnnularBlock) -> bool {
    actual.id == expected.id
        && actual.degrees == expected.degrees
        && actual.controls.len() == expected.controls.len()
        && actual
            .controls
            .iter()
            .zip(&expected.controls)
            .all(|(left, right)| {
                left.iter()
                    .zip(right)
                    .all(|(a, b)| a.is_finite() && same_bits(*a, *b))
            })
}

/// Audit exact stored controls and return conservative continuous bounds.
///
/// The accepted schema contains exactly two degree-(2,1,1) blocks with
/// heights in `[1,1000]` mm. Every coefficient is compared bit-for-bit with an
/// independent reconstruction, including rounded products `q*n` for all
/// possible axial levels. The determinant, radial and axial bounds apply over
/// both entire closed cubes, not a sample grid.
pub fn audit_annular(model: &AnnularPair) -> Result<AnnularCertificate, &'static str> {
    if !(1..=HEIGHT_LIMIT_MM).contains(&model.heights_mm[0])
        || !(1..=HEIGHT_LIMIT_MM).contains(&model.heights_mm[1])
    {
        return Err("invalid_heights");
    }
    if model.blocks.len() != 2 {
        return Err("invalid_block_count");
    }
    if model.interface != ["lower.w1", "upper.w0", "identity_uv"] {
        return Err("invalid_interface");
    }

    let lower = expected_block("lower", 0, model.heights_mm[0]);
    let upper = expected_block("upper", model.heights_mm[0], model.heights_mm[1]);
    if !same_block(&model.blocks[0], &lower) || !same_block(&model.blocks[1], &upper) {
        return Err("coefficient_mismatch");
    }
    if model.blocks[0].id == model.blocks[1].id {
        return Err("duplicate_block_identity");
    }
    if !validate_fixed_dyadic_bounds() {
        return Err("internal_bound_verification_failed");
    }

    let minimum_height = model.heights_mm[0].min(model.heights_mm[1]) as f64;
    Ok(AnnularCertificate {
        min_weight: 0.5,
        determinant_lower_bound_mm3: minimum_height / 2.0,
        radial_error_upper_mm: 2.0_f64.powi(-50),
        height_error_upper_mm: 2.0_f64.powi(-43),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dyadic_bound_checks_are_exact_and_hold() {
        assert!(validate_fixed_dyadic_bounds());
    }
}
