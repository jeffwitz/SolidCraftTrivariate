//! Outward-rounded local regularity certificate for small rational Bézier patches.
//!
//! The supported degree families are `(1,1)` and `(2,1)`. Homogeneous integer
//! controls are scaled by the exact power of two `2^-63` before conversion.
//! Every converted endpoint encloses the exact integer; every arithmetic
//! operation expands its result by one binary64 value on each side. This
//! includes Bernstein product factors and the final division by an upper bound
//! for `W^3`.

/// Why a planar regularity certificate was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlanarError {
    /// The degrees, control count, or positive-weight premise is invalid, or a
    /// closed-square corner proves that the determinant is nonpositive.
    Invalid,
    /// Valid input did not yield a strictly positive Bernstein lower bound.
    Indeterminate,
}

impl std::fmt::Display for PlanarError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(match self {
            Self::Invalid => "invalid planar patch or nonpositive corner determinant",
            Self::Indeterminate => "outward interval Bernstein bound is inconclusive",
        })
    }
}

impl std::error::Error for PlanarError {}

#[derive(Clone, Copy, Debug)]
struct Interval {
    lo: f64,
    hi: f64,
}

impl Interval {
    fn exact(value: f64) -> Self {
        Self {
            lo: value,
            hi: value,
        }
    }

    fn from_i64_scaled(value: i64) -> Self {
        // 2^-63 is exact, and all scaled i64 values remain normal binary64.
        let scale = f64::from_bits((960u64) << 52);
        let rounded = (value as f64) * scale;
        Self {
            lo: next_down(rounded),
            hi: next_up(rounded),
        }
    }

    fn add(self, rhs: Self) -> Self {
        Self {
            lo: next_down(self.lo + rhs.lo),
            hi: next_up(self.hi + rhs.hi),
        }
    }

    fn sub(self, rhs: Self) -> Self {
        Self {
            lo: next_down(self.lo - rhs.hi),
            hi: next_up(self.hi - rhs.lo),
        }
    }

    fn mul(self, rhs: Self) -> Self {
        let products = [
            self.lo * rhs.lo,
            self.lo * rhs.hi,
            self.hi * rhs.lo,
            self.hi * rhs.hi,
        ];
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for value in products {
            low = low.min(value);
            high = high.max(value);
        }
        Self {
            lo: next_down(low),
            hi: next_up(high),
        }
    }

    fn div(self, rhs: Self) -> Self {
        if rhs.lo <= 0.0 || !rhs.lo.is_finite() || !rhs.hi.is_finite() {
            return Self {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            };
        }
        let quotients = [
            self.lo / rhs.lo,
            self.lo / rhs.hi,
            self.hi / rhs.lo,
            self.hi / rhs.hi,
        ];
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for value in quotients {
            low = low.min(value);
            high = high.max(value);
        }
        Self {
            lo: next_down(low),
            hi: next_up(high),
        }
    }

    fn finite(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite() && self.lo <= self.hi
    }
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

#[derive(Clone)]
struct Bernstein {
    degrees: [usize; 2],
    values: Vec<Interval>,
}

impl Bernstein {
    fn at(&self, i: usize, j: usize) -> Interval {
        self.values[i * (self.degrees[1] + 1) + j]
    }

    fn constant_zero(degrees: [usize; 2]) -> Self {
        Self {
            degrees,
            values: vec![Interval::exact(0.0); (degrees[0] + 1) * (degrees[1] + 1)],
        }
    }

    fn derivative(&self, axis: usize) -> Self {
        let [du, dv] = self.degrees;
        if axis == 0 {
            let mut values = Vec::with_capacity(du * (dv + 1));
            for i in 0..du {
                for j in 0..=dv {
                    values.push(
                        self.at(i + 1, j)
                            .sub(self.at(i, j))
                            .mul(Interval::exact(du as f64)),
                    );
                }
            }
            Self {
                degrees: [du - 1, dv],
                values,
            }
        } else {
            let mut values = Vec::with_capacity((du + 1) * dv);
            for i in 0..=du {
                for j in 0..dv {
                    values.push(
                        self.at(i, j + 1)
                            .sub(self.at(i, j))
                            .mul(Interval::exact(dv as f64)),
                    );
                }
            }
            Self {
                degrees: [du, dv - 1],
                values,
            }
        }
    }

    fn add(&self, rhs: &Self) -> Self {
        debug_assert_eq!(self.degrees, rhs.degrees);
        Self {
            degrees: self.degrees,
            values: self
                .values
                .iter()
                .zip(&rhs.values)
                .map(|(a, b)| a.add(*b))
                .collect(),
        }
    }

    fn subtract(&self, rhs: &Self) -> Self {
        debug_assert_eq!(self.degrees, rhs.degrees);
        Self {
            degrees: self.degrees,
            values: self
                .values
                .iter()
                .zip(&rhs.values)
                .map(|(a, b)| a.sub(*b))
                .collect(),
        }
    }

    /// Restrict a scalar polynomial along one axis, preserving its values
    /// in parent coordinates (no derivative rescaling).
    fn bisect(&self, axis: usize) -> [Self; 2] {
        let mut children = [
            Self::constant_zero(self.degrees),
            Self::constant_zero(self.degrees),
        ];
        let degree = self.degrees[axis];
        for line in 0..=self.degrees[1 - axis] {
            let index = |k| {
                if axis == 0 {
                    k * (self.degrees[1] + 1) + line
                } else {
                    line * (self.degrees[1] + 1) + k
                }
            };
            let mut work: Vec<_> = (0..=degree).map(|k| self.values[index(k)]).collect();
            for depth in 0..=degree {
                children[0].values[index(depth)] = work[0];
                children[1].values[index(degree - depth)] = work[degree - depth];
                for k in 0..degree - depth {
                    work[k] = work[k].add(work[k + 1]).mul(Interval::exact(0.5));
                }
            }
        }
        children
    }

    fn positive_lower(&self, remaining: u8, axis: usize) -> Result<f64, PlanarError> {
        let [u, v] = self.degrees;
        if [(0, 0), (u, 0), (0, v), (u, v)]
            .iter()
            .any(|&(i, j)| self.at(i, j).hi <= 0.)
        {
            return Err(PlanarError::Invalid);
        }
        let lower = self
            .values
            .iter()
            .map(|p| p.lo)
            .fold(f64::INFINITY, f64::min);
        if lower > 0. {
            return Ok(lower);
        }
        if remaining == 0 {
            return Err(PlanarError::Indeterminate);
        }
        let [a, b] = self.bisect(axis);
        let left = a.positive_lower(remaining - 1, 1 - axis);
        if left == Err(PlanarError::Invalid) {
            return left;
        }
        match (left, b.positive_lower(remaining - 1, 1 - axis)) {
            (Ok(x), Ok(y)) => Ok(x.min(y)),
            (_, Err(PlanarError::Invalid)) => Err(PlanarError::Invalid),
            _ => Err(PlanarError::Indeterminate),
        }
    }

    fn product(&self, rhs: &Self) -> Self {
        let [am, an] = self.degrees;
        let [bm, bn] = rhs.degrees;
        let degrees = [am + bm, an + bn];
        let mut result = Self::constant_zero(degrees);
        for i in 0..=am {
            for j in 0..=an {
                for k in 0..=bm {
                    for ell in 0..=bn {
                        let factor_u = binomial_ratio(am, i, bm, k);
                        let factor_v = binomial_ratio(an, j, bn, ell);
                        let factor = factor_u.mul(factor_v);
                        let value = self.at(i, j).mul(rhs.at(k, ell)).mul(factor);
                        let index = (i + k) * (degrees[1] + 1) + j + ell;
                        result.values[index] = result.values[index].add(value);
                    }
                }
            }
        }
        result
    }
}

fn binomial(n: usize, k: usize) -> u64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result = 1u64;
    for i in 1..=k {
        result = result * ((n - k + i) as u64) / (i as u64);
    }
    result
}

fn binomial_ratio(m: usize, i: usize, n: usize, j: usize) -> Interval {
    let numerator = (binomial(m, i) * binomial(n, j)) as f64;
    let denominator = binomial(m + n, i + j) as f64;
    Interval::exact(numerator).div(Interval::exact(denominator))
}

fn determinant_numerator(degrees: [usize; 2], controls: &[[i64; 3]]) -> Bernstein {
    let [du, dv] = degrees;
    let components: [Bernstein; 3] = std::array::from_fn(|axis| Bernstein {
        degrees,
        values: controls
            .iter()
            .map(|point| Interval::from_i64_scaled(point[axis]))
            .collect(),
    });
    let hu: [Bernstein; 3] = std::array::from_fn(|axis| components[axis].derivative(0));
    let hv: [Bernstein; 3] = std::array::from_fn(|axis| components[axis].derivative(1));
    // det(Hu,Hv,H) = Hu_x(Hv_y W-Hv_w Y) + cyclic terms.
    let mut sum = Bernstein::constant_zero([3 * du - 1, 3 * dv - 1]);
    for (a, b, c) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
        let positive = hu[a].product(&hv[b]).product(&components[c]);
        let negative = hu[a]
            .product(&hv[(b + 1) % 3])
            .product(&components[(a + 1) % 3]);
        sum = sum.add(&positive.subtract(&negative));
    }
    sum
}

/// Certify strictly positive orientation over the complete closed square.
///
/// `controls` are flattened by the second parameter index, with each triple
/// `(wx, wy, weight)`. The return value is a conservative lower bound for the
/// physical determinant `det(H_u,H_v,H)/W^3`; errors deliberately distinguish
/// invalid input or a nonpositive corner from an inconclusive interior bound.
pub fn certify_planar(degrees: [usize; 2], controls: &[[i64; 3]]) -> Result<f64, PlanarError> {
    certify_planar_with_depth(degrees, controls, 0)
}

/// Refine inconclusive Bernstein envelopes on the closed parameter square.
/// A depth budget bounds the number of successive binary cuts on any path.
pub fn certify_planar_refined(
    degrees: [usize; 2],
    controls: &[[i64; 3]],
) -> Result<f64, PlanarError> {
    certify_planar_with_depth(degrees, controls, 12)
}

fn certify_planar_with_depth(
    degrees: [usize; 2],
    controls: &[[i64; 3]],
    depth: u8,
) -> Result<f64, PlanarError> {
    let [du, dv] = degrees;
    if !matches!(degrees, [1, 1] | [2, 1])
        || controls.len() != (du + 1) * (dv + 1)
        || controls.iter().any(|point| point[2] <= 0)
    {
        return Err(PlanarError::Invalid);
    }
    let numerator = determinant_numerator(degrees, controls);
    if numerator.values.iter().any(|value| !value.finite()) {
        return Err(PlanarError::Indeterminate);
    }

    let numerator_lower = numerator.positive_lower(depth, 0)?;
    let weight_upper = controls
        .iter()
        .map(|point| Interval::from_i64_scaled(point[2]).hi)
        .fold(0.0, f64::max);
    let weight_cubed_upper = Interval::exact(weight_upper)
        .mul(Interval::exact(weight_upper))
        .mul(Interval::exact(weight_upper))
        .hi;
    let lower_bound = next_down(numerator_lower / weight_cubed_upper);
    if !lower_bound.is_finite() || lower_bound <= 0.0 {
        return Err(PlanarError::Indeterminate);
    }
    Ok(lower_bound)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn subdivision_resolves_both_a_positive_minimum_and_a_negative_midpoint() {
        // Quadratic Bernstein controls (1,c,1) have exact extremum
        // (1+c)/2 at 1/2. These polynomials are constant along the second axis.
        let polynomial = |c| Bernstein {
            degrees: [2, 1],
            values: [1., 1., c, c, 1., 1.].map(Interval::exact).into(),
        };
        assert_eq!(
            polynomial(-0.75).positive_lower(0, 0),
            Err(PlanarError::Indeterminate)
        );
        let bound = polynomial(-0.75).positive_lower(1, 0).unwrap();
        assert!(bound > 0.124999 && bound <= 0.125);
        assert_eq!(
            polynomial(-2.).positive_lower(1, 0),
            Err(PlanarError::Invalid)
        );
        assert_eq!(
            polynomial(-1.).positive_lower(4, 0),
            Err(PlanarError::Indeterminate)
        );
    }

    #[test]
    fn affine_map_has_positive_exact_enclosure() {
        let controls = [[0, 0, 1], [0, 1, 1], [1, 0, 1], [1, 1, 1]];
        let lower = certify_planar([1, 1], &controls).unwrap();
        assert!(lower > 0.0 && lower <= 1.0);
    }

    #[test]
    fn reversed_affine_map_is_invalid_at_a_corner() {
        let controls = [[0, 0, 1], [0, 1, 1], [1, 0, 1], [-1, 1, 1]];
        assert_eq!(certify_planar([1, 1], &controls), Err(PlanarError::Invalid));
    }

    #[test]
    fn homogeneous_integer_scale_does_not_change_certificate_status() {
        let base = [[0, 0, 1], [0, 1, 1], [1, 0, 1], [1, 1, 1]];
        let scaled = base.map(|point| [point[0] * 7, point[1] * 7, point[2] * 7]);
        assert!(certify_planar([1, 1], &base).is_ok());
        assert!(certify_planar([1, 1], &scaled).is_ok());
    }

    #[test]
    fn east_pierced_plate_profile_bound_is_below_independent_center_value() {
        // E2 east-sector controls, integerized by the same positive common
        // factor used for native storage. Evaluate H and its two derivatives
        // directly at (1/2,1/2), independently of Bernstein products here.
        let controls: [[i64; 3]; 6] = [
            [90071992547409920, -90071992547409920, 9007199254740992],
            [25476206690103092, -25476206690103092, 9007199254740992],
            [90071992547409920, 0, 9007199254740992],
            [36028797018963976, 0, 6369051672525773],
            [90071992547409920, 90071992547409920, 9007199254740992],
            [25476206690103092, 25476206690103092, 9007199254740992],
        ];
        let certified = certify_planar([2, 1], &controls).unwrap();
        let scale = f64::from_bits(960u64 << 52);
        let bu = [0.25, 0.5, 0.25];
        let du = [-1.0, 0.0, 1.0];
        let bv = [0.5, 0.5];
        let dv = [-1.0, 1.0];
        let mut h = [0.0; 3];
        let mut hu = [0.0; 3];
        let mut hv = [0.0; 3];
        for i in 0..3 {
            for j in 0..2 {
                let point = controls[i * 2 + j];
                let factor = bu[i] * bv[j];
                for c in 0..3 {
                    let value = point[c] as f64 * scale;
                    h[c] += factor * value;
                    hu[c] += du[i] * bv[j] * value;
                    hv[c] += bu[i] * dv[j] * value;
                }
            }
        }
        let numerator = hu[0] * (hv[1] * h[2] - hv[2] * h[1])
            - hu[1] * (hv[0] * h[2] - hv[2] * h[0])
            + hu[2] * (hv[0] * h[1] - hv[1] * h[0]);
        let center_value = numerator / h[2].powi(3);
        assert!(center_value > 0.0);
        assert!(certified <= center_value);
    }
}
