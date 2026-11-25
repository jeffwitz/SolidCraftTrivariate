//! Certified piecewise-linear graph approximation for exact planar NURBS.
//!
//! Every accepted segment carries an exact Bernstein convex-hull bound on
//! its lateral deviation from the endpoint chord. Unsupported or exhausted
//! subdivisions are rejected instead of returning an uncertified polyline.

use crate::{
    constraint_system::Ratio,
    nurbs_cross::NurbsCurve2,
    nurbs_curve::{extract_bezier, BezierSpan},
};

type Error = &'static str;
type HPoint = [Ratio; 3];
type EPoint = [Ratio; 2];
const MAX_SEGMENTS: usize = 65_536;

fn int(value: i64) -> Ratio {
    Ratio::new(i128::from(value), 1).expect("unit denominator")
}

fn div(a: Ratio, b: Ratio) -> Result<Ratio, Error> {
    if b.numerator() == 0 {
        return Err("zero rational denominator");
    }
    a.checked_mul(Ratio::new(b.denominator(), b.numerator())?)
}

fn sign(value: Ratio) -> i8 {
    value.numerator().signum() as i8
}

fn euclidean(point: HPoint) -> Result<EPoint, Error> {
    if point[2].numerator() <= 0 {
        return Err("NURBS span has a nonpositive weight");
    }
    Ok([div(point[0], point[2])?, div(point[1], point[2])?])
}

fn blend(a: HPoint, b: HPoint) -> Result<HPoint, Error> {
    let half = Ratio::new(1, 2)?;
    let mut result = a;
    for c in 0..3 {
        result[c] = a[c].checked_add(half.checked_mul(b[c].checked_sub(a[c])?)?)?;
    }
    Ok(result)
}

fn split_half(net: &[HPoint]) -> Result<(Vec<HPoint>, Vec<HPoint>), Error> {
    let mut row = net.to_vec();
    let mut left = Vec::with_capacity(net.len());
    let mut right = Vec::with_capacity(net.len());
    while !row.is_empty() {
        left.push(row[0]);
        right.push(row[row.len() - 1]);
        row = row
            .windows(2)
            .map(|pair| blend(pair[0], pair[1]))
            .collect::<Result<_, _>>()?;
    }
    right.reverse();
    Ok((left, right))
}

fn monotone_y(net: &[HPoint]) -> Result<(EPoint, EPoint), Error> {
    let points = net
        .iter()
        .copied()
        .map(euclidean)
        .collect::<Result<Vec<_>, _>>()?;
    if points.len() < 2 {
        return Err("NURBS Bezier span needs at least two controls");
    }
    for pair in points.windows(2) {
        if sign(pair[1][1].checked_sub(pair[0][1])?) < 0 {
            return Err("NURBS ordinate controls do not prove graph monotonicity");
        }
    }
    let first = points[0];
    let last = *points.last().unwrap();
    if sign(last[1].checked_sub(first[1])?) <= 0 {
        return Err("NURBS ordinate must increase strictly over each span");
    }
    Ok((first, last))
}

/// Bound the lateral residual against the chord through the endpoint values.
/// The residual is affine in the Euclidean point, so its extrema over the
/// positive-weight rational Bezier curve are bounded by its control values.
fn chord_error(net: &[HPoint], endpoints: (EPoint, EPoint)) -> Result<Ratio, Error> {
    let [x0, y0] = endpoints.0;
    let [x1, y1] = endpoints.1;
    let dx = x1.checked_sub(x0)?;
    let dy = y1.checked_sub(y0)?;
    if dy.numerator() <= 0 {
        return Err("NURBS ordinate must increase strictly over each segment");
    }
    let mut bound = Ratio::zero();
    for &control in net {
        let [x, y] = euclidean(control)?;
        let chord_x = x0.checked_add(dx.checked_mul(div(y.checked_sub(y0)?, dy)?)?)?;
        let residual = x.checked_sub(chord_x)?;
        let magnitude = Ratio::new(
            residual
                .numerator()
                .checked_abs()
                .ok_or("NURBS residual magnitude overflow")?,
            residual.denominator(),
        )?;
        if sign(magnitude.checked_sub(bound)?) > 0 {
            bound = magnitude;
        }
    }
    Ok(bound)
}

fn append_certified(
    net: Vec<HPoint>,
    tolerance: Ratio,
    max_depth: usize,
    depth: usize,
    output: &mut Vec<EPoint>,
) -> Result<(), Error> {
    if output.len() > MAX_SEGMENTS {
        return Err("certified NURBS graph exceeds global segment budget");
    }
    let endpoints = monotone_y(&net)?;
    let error = chord_error(&net, endpoints)?;
    if sign(error.checked_sub(tolerance)?) <= 0 {
        if output.last().copied() != Some(endpoints.0) {
            if !output.is_empty() {
                return Err("NURBS subdivision endpoints do not join exactly");
            }
            output.push(endpoints.0);
        }
        output.push(endpoints.1);
        if output.len() > MAX_SEGMENTS + 1 {
            return Err("certified NURBS graph exceeds global segment budget");
        }
        return Ok(());
    }
    if depth >= max_depth {
        return Err("NURBS graph tolerance not certified within subdivision depth");
    }
    let (left, right) = split_half(&net)?;
    append_certified(left, tolerance, max_depth, depth + 1, output)?;
    append_certified(right, tolerance, max_depth, depth + 1, output)
}

fn approximate_span(
    span: &BezierSpan,
    tolerance: Ratio,
    max_depth: usize,
    output: &mut Vec<EPoint>,
) -> Result<(), Error> {
    if span.controls.len() < 2 || span.controls.iter().any(|point| point[2] <= 0) {
        return Err("NURBS span controls must have positive weights");
    }
    let net = span.controls.iter().map(|point| point.map(int)).collect();
    append_certified(net, tolerance, max_depth, 0, output)
}

/// Linearize a planar NURBS curve under an exact uniform lateral tolerance.
///
/// The returned points are Euclidean `(coordinate, ordinate)` pairs and are
/// strictly increasing in ordinate. `max_depth` is the maximum bisection
/// depth per extracted Bezier span and must not exceed 12.
pub fn certified_graph(
    curve: &NurbsCurve2,
    tolerance: Ratio,
    max_depth: usize,
) -> Result<Vec<[Ratio; 2]>, Error> {
    if tolerance.numerator() <= 0 {
        return Err("NURBS graph tolerance must be positive");
    }
    if max_depth > 12 {
        return Err("NURBS graph subdivision depth exceeds 12");
    }
    let spans = extract_bezier(curve.degree, &curve.knots, &curve.controls)?;
    if spans.is_empty() {
        return Err("NURBS graph has empty support");
    }
    let mut points = Vec::new();
    for span in &spans {
        approximate_span(span, tolerance, max_depth, &mut points)?;
    }
    if points.len() < 2 {
        return Err("NURBS graph has insufficient points");
    }
    for pair in points.windows(2) {
        if sign(pair[1][1].checked_sub(pair[0][1])?) <= 0 {
            return Err("NURBS graph ordinates are not strictly increasing");
        }
    }
    Ok(points)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quadratic_graph() -> NurbsCurve2 {
        NurbsCurve2 {
            degree: 2,
            knots: vec![0, 0, 0, 1, 1, 1],
            controls: vec![[0, 0, 1], [1, 1, 1], [0, 2, 1]],
        }
    }

    #[test]
    fn certifies_exact_midpoint_subdivision_and_rejects_exhausted_depth() {
        let curve = quadratic_graph();
        let tolerance = Ratio::new(1, 4).unwrap();
        assert_eq!(
            certified_graph(&curve, tolerance, 0),
            Err("NURBS graph tolerance not certified within subdivision depth")
        );
        let points = certified_graph(&curve, tolerance, 1).unwrap();
        assert_eq!(points.len(), 3);
        assert_eq!(points[0], [Ratio::zero(), Ratio::zero()]);
        assert_eq!(points[1], [Ratio::new(1, 2).unwrap(), Ratio::one()]);
        assert_eq!(points[2], [Ratio::zero(), Ratio::new(2, 1).unwrap()]);
    }

    #[test]
    fn refuses_nonmonotone_ordinate_controls_and_invalid_limits() {
        let mut curve = quadratic_graph();
        curve.controls[1][1] = 3;
        assert!(certified_graph(&curve, Ratio::one(), 2).is_err());
        assert_eq!(
            certified_graph(&quadratic_graph(), Ratio::zero(), 4),
            Err("NURBS graph tolerance must be positive")
        );
        assert_eq!(
            certified_graph(&quadratic_graph(), Ratio::one(), 13),
            Err("NURBS graph subdivision depth exceeds 12")
        );
    }
}
