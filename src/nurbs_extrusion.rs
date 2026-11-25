//! Extrusion of paired rational NURBS boundary curves into tensor-product strips.
//!
//! The admitted class has a common rational x-map and common weights on each
//! synchronized Bezier span. The two y-curves must remain strictly separated,
//! and x must be strictly monotone. These conditions give an exact positive
//! Jacobian certificate for the resulting `[p, 1, 1]` blocks.

use crate::{bracket::BracketBlock, constraint_system::Ratio, nurbs_curve::BezierSpan};
use std::convert::TryFrom;

type HomogeneousCurve = Vec<[i128; 3]>;
type AlignedPair = (HomogeneousCurve, HomogeneousCurve);

fn checked_add(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_add(b)
}

fn checked_sub(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_sub(b)
}

fn checked_mul(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_mul(b)
}

fn checked_div(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    let numerator = a
        .numerator()
        .checked_mul(b.denominator())
        .ok_or("rational division overflow")?;
    let denominator = a
        .denominator()
        .checked_mul(b.numerator())
        .ok_or("rational denominator overflow")?;
    Ratio::new(numerator, denominator)
}

fn choose(n: usize, k: usize) -> i128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result = 1_i128;
    for i in 1..=k {
        result = result * (n - i + 1) as i128 / i as i128;
    }
    result
}

fn ratio(value: i128) -> Result<Ratio, &'static str> {
    Ratio::new(value, 1)
}

fn bernstein_to_power(values: &[Ratio]) -> Result<Vec<Ratio>, &'static str> {
    let degree = values.len().checked_sub(1).ok_or("empty Bezier controls")?;
    let mut power = vec![Ratio::zero(); degree + 1];
    for (i, &value) in values.iter().enumerate() {
        for (k, slot) in power.iter_mut().enumerate().skip(i) {
            let coefficient = choose(degree, i)
                .checked_mul(choose(degree - i, k - i))
                .ok_or("polynomial coefficient overflow")?;
            let signed = if (k - i) % 2 == 0 {
                coefficient
            } else {
                -coefficient
            };
            *slot = checked_add(*slot, checked_mul(value, ratio(signed)?)?)?;
        }
    }
    Ok(power)
}

fn multiply_power(a: &[Ratio], b: &[Ratio]) -> Result<Vec<Ratio>, &'static str> {
    let length = a
        .len()
        .checked_add(b.len())
        .and_then(|n| n.checked_sub(1))
        .ok_or("polynomial degree overflow")?;
    let mut product = vec![Ratio::zero(); length];
    for (i, &left) in a.iter().enumerate() {
        for (j, &right) in b.iter().enumerate() {
            product[i + j] = checked_add(product[i + j], checked_mul(left, right)?)?;
        }
    }
    Ok(product)
}

fn derivative_power(values: &[Ratio]) -> Result<Vec<Ratio>, &'static str> {
    if values.len() <= 1 {
        return Ok(vec![Ratio::zero()]);
    }
    values
        .iter()
        .enumerate()
        .skip(1)
        .map(|(degree, &value)| checked_mul(value, ratio(degree as i128)?))
        .collect()
}

fn power_to_bernstein(values: &[Ratio], degree: usize) -> Result<Vec<Ratio>, &'static str> {
    if values.len() > degree + 1 {
        return Err("polynomial degree exceeds Bernstein degree");
    }
    let mut result = vec![Ratio::zero(); degree + 1];
    for (i, slot) in result.iter_mut().enumerate() {
        for (k, &coefficient) in values.iter().enumerate().take(i.min(values.len() - 1) + 1) {
            let denominator = choose(degree, k);
            if denominator == 0 {
                return Err("unsupported Bernstein degree");
            }
            let factor = checked_div(ratio(choose(i, k))?, ratio(denominator)?)?;
            *slot = checked_add(*slot, checked_mul(coefficient, factor)?)?;
        }
    }
    Ok(result)
}

fn split_bernstein(values: &[Ratio]) -> Result<(Vec<Ratio>, Vec<Ratio>), &'static str> {
    let n = values.len();
    let mut row = values.to_vec();
    let mut left = vec![Ratio::zero(); n];
    let mut right = vec![Ratio::zero(); n];
    left[0] = row[0];
    right[n - 1] = row[n - 1];
    for level in 1..n {
        for i in 0..n - level {
            row[i] = checked_div(checked_add(row[i], row[i + 1])?, ratio(2)?)?;
        }
        left[level] = row[0];
        right[n - level - 1] = row[n - level - 1];
    }
    Ok((left, right))
}

/// Prove a Bernstein polynomial has one strict sign; return its sign and a
/// positive rational lower bound on its absolute value.
fn strict_sign(values: &[Ratio]) -> Result<(i8, Ratio), &'static str> {
    fn recurse(values: &[Ratio], depth: u8) -> Result<(i8, Ratio), &'static str> {
        let positive = values.iter().all(|value| value.numerator() > 0);
        let negative = values.iter().all(|value| value.numerator() < 0);
        if positive || negative {
            let sign = if positive { 1 } else { -1 };
            let lower = values
                .iter()
                .map(|value| {
                    if sign > 0 {
                        Ok(*value)
                    } else {
                        Ratio::new(
                            value
                                .numerator()
                                .checked_neg()
                                .ok_or("rational negation overflow")?,
                            value.denominator(),
                        )
                    }
                })
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .min()
                .ok_or("empty polynomial")?;
            return Ok((sign, lower));
        }
        if depth == 0 {
            return Err("strict sign is not certified by Bernstein subdivision");
        }
        let (left, right) = split_bernstein(values)?;
        let (left_sign, left_bound) = recurse(&left, depth - 1)?;
        let (right_sign, right_bound) = recurse(&right, depth - 1)?;
        if left_sign != right_sign {
            return Err("polynomial changes sign");
        }
        Ok((left_sign, left_bound.min(right_bound)))
    }
    recurse(values, 20)
}

fn monotone_x_sign(xw: &[[i128; 2]]) -> Result<(i8, Ratio), &'static str> {
    let xs = xw
        .iter()
        .map(|point| ratio(point[0]))
        .collect::<Result<Vec<_>, _>>()?;
    let ws = xw
        .iter()
        .map(|point| ratio(point[1]))
        .collect::<Result<Vec<_>, _>>()?;
    let x = bernstein_to_power(&xs)?;
    let w = bernstein_to_power(&ws)?;
    let first = multiply_power(&derivative_power(&x)?, &w)?;
    let second = multiply_power(&x, &derivative_power(&w)?)?;
    let mut numerator = vec![Ratio::zero(); first.len().max(second.len())];
    for (i, &value) in first.iter().enumerate() {
        numerator[i] = checked_add(numerator[i], value)?;
    }
    for (i, &value) in second.iter().enumerate() {
        numerator[i] = checked_sub(numerator[i], value)?;
    }
    let degree = xw.len().saturating_mul(2).saturating_sub(3);
    let bernstein = power_to_bernstein(&numerator, degree)?;
    strict_sign(&bernstein)
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn gcd_values(values: impl IntoIterator<Item = i64>) -> Result<i128, &'static str> {
    let mut divisor = 0_u128;
    for value in values {
        divisor = gcd(divisor, i128::from(value).unsigned_abs());
    }
    if divisor == 0 || divisor > i128::MAX as u128 {
        return Err("invalid homogeneous gauge");
    }
    Ok(divisor as i128)
}

fn align_spans(lower: &BezierSpan, upper: &BezierSpan) -> Result<AlignedPair, &'static str> {
    if lower.interval != upper.interval
        || lower.interval[0] >= lower.interval[1]
        || lower.controls.len() != upper.controls.len()
        || !(2..=4).contains(&lower.controls.len())
    {
        return Err("NURBS spans must have matching nonempty intervals and degree 1..=3");
    }
    if lower
        .controls
        .iter()
        .chain(&upper.controls)
        .any(|point| point[2] <= 0)
    {
        return Err("NURBS span weights must be strictly positive");
    }
    let lower_gauge = gcd_values(lower.controls.iter().flat_map(|point| [point[0], point[2]]))?;
    let upper_gauge = gcd_values(upper.controls.iter().flat_map(|point| [point[0], point[2]]))?;
    let common = gcd(lower_gauge as u128, upper_gauge as u128) as i128;
    let lower_scale = upper_gauge / common;
    let upper_scale = lower_gauge / common;
    let lower_net = lower
        .controls
        .iter()
        .map(|point| {
            point
                .map(|value| {
                    i128::from(value)
                        .checked_mul(lower_scale)
                        .ok_or("aligned NURBS control overflow")
                })
                .into_iter()
                .collect::<Result<Vec<_>, _>>()
                .and_then(|values| values.try_into().map_err(|_| "invalid NURBS control"))
        })
        .collect::<Result<Vec<[i128; 3]>, _>>()?;
    let upper_net = upper
        .controls
        .iter()
        .map(|point| {
            point
                .map(|value| {
                    i128::from(value)
                        .checked_mul(upper_scale)
                        .ok_or("aligned NURBS control overflow")
                })
                .into_iter()
                .collect::<Result<Vec<_>, _>>()
                .and_then(|values| values.try_into().map_err(|_| "invalid NURBS control"))
        })
        .collect::<Result<Vec<[i128; 3]>, _>>()?;
    if lower_net
        .iter()
        .zip(&upper_net)
        .any(|(a, b)| a[0] != b[0] || a[2] != b[2])
    {
        return Err("paired NURBS spans do not share the same rational x-map");
    }
    let divisor = lower_net
        .iter()
        .chain(&upper_net)
        .flatten()
        .fold(0_u128, |value, item| gcd(value, item.unsigned_abs()));
    if divisor == 0 || divisor > i128::MAX as u128 {
        return Err("invalid aligned NURBS gauge");
    }
    let divisor = divisor as i128;
    let reduce = |net: &[[i128; 3]]| -> Result<Vec<[i128; 3]>, &'static str> {
        net.iter()
            .map(|point| {
                if point.iter().any(|value| value % divisor != 0) {
                    return Err("unaligned NURBS gauge reduction");
                }
                Ok(point.map(|value| value / divisor))
            })
            .collect()
    };
    Ok((reduce(&lower_net)?, reduce(&upper_net)?))
}

fn same_physical_endpoint(a: [i128; 3], b: [i128; 3]) -> Result<bool, &'static str> {
    if a[2] <= 0 || b[2] <= 0 {
        return Err("nonpositive NURBS endpoint weight");
    }
    for coordinate in 0..2 {
        let left = a[coordinate]
            .checked_mul(b[2])
            .ok_or("NURBS endpoint comparison overflow")?;
        let right = b[coordinate]
            .checked_mul(a[2])
            .ok_or("NURBS endpoint comparison overflow")?;
        if left != right {
            return Ok(false);
        }
    }
    Ok(true)
}

fn f64_lower_bound(value: Ratio) -> Result<f64, &'static str> {
    if value.numerator() <= 0 {
        return Err("determinant bound is not positive");
    }
    let approximate = value.numerator() as f64 / value.denominator() as f64;
    if !approximate.is_finite() {
        return Err("determinant bound is outside f64 range");
    }
    if approximate == 0.0 {
        return Ok(0.0);
    }
    // Both integer conversions and the division can round upward. Move down
    // several representable values; exact rational sign remains authoritative.
    let mut bound = approximate;
    for _ in 0..4 {
        bound = f64::from_bits(bound.to_bits() - 1);
    }
    Ok(bound.max(0.0))
}

/// Return a conservative floating lower bound for the determinant of one
/// admitted extrusion block. Exact sign proof is performed with rational
/// Bernstein bounds; an underflowed numeric bound is reported as zero.
pub fn determinant_bound(block: &BracketBlock) -> Result<f64, &'static str> {
    determinant_bound_at(block, 0)
}

/// Certify a strip extruded from the specified integer altitude.
pub(crate) fn determinant_bound_at(block: &BracketBlock, bottom: i64) -> Result<f64, &'static str> {
    if block.degrees[1..] != [1, 1] || !(1..=3).contains(&block.degrees[0]) {
        return Err("not a supported NURBS extrusion block");
    }
    let p = block.degrees[0];
    if block.controls.len() != (p + 1) * 4 {
        return Err("invalid NURBS extrusion control count");
    }
    let index = |u: usize, v: usize, w: usize| (u * 2 + v) * 2 + w;
    let mut common_height = None;
    for u in 0..=p {
        let base = block.controls[index(u, 0, 0)];
        if base[3] <= 0 || i128::from(base[2]) != i128::from(bottom) * i128::from(base[3]) {
            return Err("invalid bottom extrusion control");
        }
        for v in 0..=1 {
            for w in 0..=1 {
                let point = block.controls[index(u, v, w)];
                if point[3] != base[3]
                    || point[0] != base[0]
                    || point[1] != block.controls[index(u, v, 0)][1]
                {
                    return Err("extrusion x-map or weights vary across the strip");
                }
                if w == 0 {
                    if i128::from(point[2]) != i128::from(bottom) * i128::from(point[3]) {
                        return Err("inconsistent bottom extrusion altitude");
                    }
                } else {
                    let weight = i128::from(point[3]);
                    let z = i128::from(point[2]);
                    if weight <= 0 || z % weight != 0 {
                        return Err("nonintegral extrusion height");
                    }
                    let height = z / weight - i128::from(bottom);
                    if height <= 0 || common_height.is_some_and(|expected| expected != height) {
                        return Err("extrusion height is not positive and constant");
                    }
                    common_height = Some(height);
                }
            }
        }
    }
    let xw = (0..=p)
        .map(|u| {
            let point = block.controls[index(u, 0, 0)];
            [i128::from(point[0]), i128::from(point[3])]
        })
        .collect::<Vec<_>>();
    let (derivative_sign, derivative_bound) = monotone_x_sign(&xw)?;
    if derivative_sign <= 0 {
        return Err("extrusion x-map is not positively oriented");
    }
    let lower_y = (0..=p)
        .map(|u| i128::from(block.controls[index(u, 0, 0)][1]))
        .collect::<Vec<_>>();
    let upper_y = (0..=p)
        .map(|u| i128::from(block.controls[index(u, 1, 0)][1]))
        .collect::<Vec<_>>();
    let gap = lower_y
        .iter()
        .zip(upper_y)
        .map(|(&a, b)| ratio(b.checked_sub(a).ok_or("separation overflow")?))
        .collect::<Result<Vec<_>, _>>()?;
    let (gap_sign, gap_bound) = strict_sign(&gap)?;
    if gap_sign <= 0 {
        return Err("upper NURBS boundary is not strictly above lower boundary");
    }
    let weight_max = (0..=p)
        .map(|u| i128::from(block.controls[index(u, 0, 0)][3]))
        .max()
        .ok_or("empty weight net")?;
    if weight_max <= 0 {
        return Err("nonpositive extrusion weight");
    }
    let height = common_height.ok_or("missing extrusion height")?;
    if height <= 0 {
        return Err("extrusion height must be positive");
    }
    let denominator = checked_mul(
        ratio(weight_max)?,
        checked_mul(ratio(weight_max)?, ratio(weight_max)?)?,
    )?;
    let bound = checked_div(
        checked_mul(ratio(height)?, checked_mul(derivative_bound, gap_bound)?)?,
        denominator,
    )?;
    f64_lower_bound(bound)
}

/// Extrude matching rational boundary curves through a positive height.
///
/// Each input span is a degree-1..3 homogeneous Bezier net `[wx, wy, w]`.
/// Matching intervals are paired in order. Per-span homogeneous gauges may
/// differ; the routine aligns them exactly with checked integer arithmetic.
/// The paired x-map and weights must be identical after alignment, x must be
/// strictly monotone on each span with one global direction, and the upper y
/// curve must be strictly above the lower curve everywhere.
pub fn extrude_strips(
    lower: &[BezierSpan],
    upper: &[BezierSpan],
    height: i32,
) -> Result<Vec<BracketBlock>, &'static str> {
    if lower.is_empty() || lower.len() != upper.len() || height <= 0 {
        return Err("NURBS extrusion requires paired spans and positive height");
    }
    let mut aligned = lower
        .iter()
        .zip(upper)
        .map(|(lo, hi)| align_spans(lo, hi))
        .collect::<Result<Vec<_>, _>>()?;
    for index in 0..lower.len().saturating_sub(1) {
        if lower[index].interval[1] != lower[index + 1].interval[0]
            || upper[index].interval[1] != upper[index + 1].interval[0]
        {
            return Err("NURBS span intervals are not contiguous");
        }
        if !same_physical_endpoint(
            *aligned[index].0.last().ok_or("empty lower span")?,
            *aligned[index + 1].0.first().ok_or("empty lower span")?,
        )? || !same_physical_endpoint(
            *aligned[index].1.last().ok_or("empty upper span")?,
            *aligned[index + 1].1.first().ok_or("empty upper span")?,
        )? {
            return Err("NURBS boundary spans are not position-continuous");
        }
    }
    let mut direction = None;
    for (lo, _) in &aligned {
        let sign = monotone_x_sign(
            &lo.iter()
                .map(|point| [point[0], point[2]])
                .collect::<Vec<_>>(),
        )?
        .0;
        if direction.is_some_and(|old| old != sign) {
            return Err("NURBS x-direction changes between spans");
        }
        direction = Some(sign);
    }
    let reverse = direction == Some(-1);
    if reverse {
        aligned.reverse();
        for (lo, hi) in &mut aligned {
            lo.reverse();
            hi.reverse();
        }
    }

    let mut blocks = Vec::with_capacity(aligned.len());
    for (span_index, ((lo, hi), knot_span)) in aligned
        .into_iter()
        .zip(if reverse {
            lower
                .iter()
                .rev()
                .zip(upper.iter().rev())
                .collect::<Vec<_>>()
        } else {
            lower.iter().zip(upper.iter()).collect::<Vec<_>>()
        })
        .enumerate()
    {
        let degree = lo.len() - 1;
        if knot_span.0.interval != knot_span.1.interval {
            return Err("NURBS span intervals do not match");
        }
        let gap = lo
            .iter()
            .zip(&hi)
            .map(|(a, b)| ratio(b[1].checked_sub(a[1]).ok_or("separation overflow")?))
            .collect::<Result<Vec<_>, _>>()?;
        let (gap_sign, _) =
            strict_sign(&gap).map_err(|_| "upper NURBS curve is not strictly separated")?;
        if gap_sign <= 0 {
            return Err("upper NURBS curve is not strictly above lower curve");
        }

        let mut controls = Vec::with_capacity((degree + 1) * 4);
        for i in 0..=degree {
            for v in 0..=1 {
                for w in 0..=1 {
                    let source = if v == 0 { lo[i] } else { hi[i] };
                    let weight = source[2];
                    let z = if w == 0 {
                        0
                    } else {
                        i128::from(height)
                            .checked_mul(weight)
                            .ok_or("extrusion z control overflow")?
                    };
                    controls.push([
                        i64::try_from(source[0]).map_err(|_| "extrusion x control exceeds i64")?,
                        i64::try_from(source[1]).map_err(|_| "extrusion y control exceeds i64")?,
                        i64::try_from(z).map_err(|_| "extrusion z control exceeds i64")?,
                        i64::try_from(weight).map_err(|_| "extrusion weight exceeds i64")?,
                    ]);
                }
            }
        }
        let block = BracketBlock {
            id: format!("nurbs_strip_{span_index}_{}", knot_span.0.interval[0]),
            degrees: [degree, 1, 1],
            controls,
        };
        determinant_bound(&block)?;
        blocks.push(block);
    }
    Ok(blocks)
}
