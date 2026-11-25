//! Exact extraction of rational Bézier spans from clamped B-spline curves.
use crate::constraint_system::Ratio;

/// One polynomial rational Bézier span over an original knot interval.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BezierSpan {
    pub interval: [i64; 2],
    /// Homogeneous controls `(w*x, w*y, w)` in a positive integral gauge.
    pub controls: Vec<[i64; 3]>,
}

type HPoint = [Ratio; 3];

fn knot_difference(a: i64, b: i64) -> i128 {
    i128::from(a) - i128::from(b)
}

fn blend_checked(a: HPoint, b: HPoint, alpha: Ratio) -> Result<HPoint, &'static str> {
    let complement = Ratio::one().checked_sub(alpha)?;
    let mut result = [Ratio::zero(); 3];
    for k in 0..3 {
        result[k] = a[k]
            .checked_mul(alpha)?
            .checked_add(b[k].checked_mul(complement)?)?;
    }
    Ok(result)
}

fn insert_knot(
    degree: usize,
    knots: &[i64],
    controls: &[HPoint],
    value: i64,
) -> Result<(Vec<i64>, Vec<HPoint>), &'static str> {
    let last = controls
        .len()
        .checked_sub(1)
        .ok_or("empty NURBS control net")?;
    let multiplicity = knots.iter().filter(|&&k| k == value).count();
    if multiplicity > degree {
        return Err("cannot insert a knot beyond its Bezier multiplicity");
    }
    let span = (degree..=last)
        .find(|&i| knots[i] <= value && value < knots[i + 1])
        .ok_or("knot insertion value outside active domain")?;
    let mut next = vec![[Ratio::zero(); 3]; controls.len() + 1];
    let unchanged = span - degree + 1;
    next[..unchanged].copy_from_slice(&controls[..unchanged]);
    let blend_start = span - degree + 1;
    let blend_end = span - multiplicity;
    if blend_start <= blend_end {
        for i in blend_start..=blend_end {
            let denominator = knot_difference(knots[i + degree], knots[i]);
            let numerator = knot_difference(value, knots[i]);
            if denominator <= 0 || numerator < 0 || numerator > denominator {
                return Err("invalid knot insertion blend interval");
            }
            let alpha = Ratio::new(numerator, denominator)?;
            next[i] = blend_checked(controls[i], controls[i - 1], alpha)?;
        }
    }
    let shifted = span - multiplicity + 1;
    next[shifted..last + 2].copy_from_slice(&controls[shifted - 1..=last]);
    let mut next_knots = knots.to_vec();
    next_knots.insert(span + 1, value);
    Ok((next_knots, next))
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

pub(crate) fn normalize<const N: usize>(
    controls: &[[Ratio; N]],
) -> Result<Vec<[i64; N]>, &'static str> {
    let mut common = 1i128;
    for p in controls {
        for x in p {
            let denominator = x.denominator();
            let gcd = i128::try_from(gcd_u128(common as u128, denominator as u128))
                .map_err(|_| "NURBS denominator scale overflow")?;
            common = (common / gcd)
                .checked_mul(denominator)
                .ok_or("NURBS denominator scale overflow")?;
        }
    }
    let mut integers = Vec::with_capacity(controls.len());
    let mut divisor = 0u128;
    for p in controls {
        let mut point = [0i128; N];
        for k in 0..N {
            point[k] = p[k]
                .numerator()
                .checked_mul(common / p[k].denominator())
                .ok_or("NURBS homogeneous coefficient overflow")?;
            divisor = gcd_u128(divisor, point[k].unsigned_abs());
        }
        integers.push(point);
    }
    if divisor == 0 {
        return Err("zero NURBS homogeneous span");
    }
    let divisor = i128::try_from(divisor).map_err(|_| "NURBS normalization overflow")?;
    integers
        .into_iter()
        .map(|p| {
            let mut point = [0; N];
            for (dst, value) in point.iter_mut().zip(p) {
                *dst =
                    i64::try_from(value / divisor).map_err(|_| "NURBS coefficient exceeds i64")?;
            }
            if point[N - 1] <= 0 {
                return Err("NURBS span has nonpositive homogeneous weight");
            }
            Ok(point)
        })
        .collect()
}

/// Extract all nonzero polynomial spans of a clamped rational B-spline.
/// Knot insertion and homogeneous normalization use checked exact rationals.
pub fn extract_bezier(
    degree: usize,
    knots: &[i64],
    controls: &[[i64; 3]],
) -> Result<Vec<BezierSpan>, &'static str> {
    if !(1..=3).contains(&degree) {
        return Err("NURBS degree must be in 1..=3");
    }
    if controls.len() < degree + 1 || knots.len() != controls.len() + degree + 1 {
        return Err("NURBS knot/control count mismatch");
    }
    if knots.windows(2).any(|w| w[0] > w[1]) {
        return Err("NURBS knots must be nondecreasing");
    }
    let mut run_start = 0;
    while run_start < knots.len() {
        let mut run_end = run_start + 1;
        while run_end < knots.len() && knots[run_end] == knots[run_start] {
            run_end += 1;
        }
        if run_end - run_start > degree + 1 {
            return Err("NURBS knot multiplicity exceeds degree plus one");
        }
        run_start = run_end;
    }
    let last = controls.len() - 1;
    let start = knots[degree];
    let end = knots[last + 1];
    if start >= end
        || knots[..=degree].iter().any(|&k| k != start)
        || knots[knots.len() - degree - 1..].iter().any(|&k| k != end)
    {
        return Err("NURBS curve must be nonempty and clamped at both ends");
    }
    if controls.iter().any(|p| p[2] <= 0) {
        return Err("NURBS weights must be strictly positive");
    }
    let mut spans = Vec::new();
    let mut original_spans = Vec::new();
    for i in degree..=last {
        if knots[i] < knots[i + 1] {
            original_spans.push([knots[i], knots[i + 1]]);
        }
    }
    if original_spans.is_empty() {
        return Err("NURBS curve has no nonzero knot spans");
    }
    let mut refined_knots = knots.to_vec();
    let mut refined_controls: Vec<HPoint> = controls
        .iter()
        .map(|p| p.map(|x| Ratio::new(i128::from(x), 1).expect("integer rational")))
        .collect();
    let internal: Vec<i64> = knots[degree + 1..=last]
        .iter()
        .copied()
        .filter(|&u| start < u && u < end)
        .collect();
    for value in internal {
        loop {
            let multiplicity = refined_knots.iter().filter(|&&k| k == value).count();
            if multiplicity >= degree {
                break;
            }
            let inserted = insert_knot(degree, &refined_knots, &refined_controls, value)?;
            refined_knots = inserted.0;
            refined_controls = inserted.1;
        }
    }
    for interval in original_spans {
        let span = (degree..refined_controls.len())
            .find(|&i| refined_knots[i] == interval[0] && refined_knots[i + 1] == interval[1])
            .ok_or("failed to locate extracted NURBS span")?;
        let bezier = &refined_controls[span - degree..=span];
        if bezier.iter().any(|p| p[2].numerator() <= 0) {
            return Err("knot insertion produced a nonpositive homogeneous weight");
        }
        spans.push(BezierSpan {
            interval,
            controls: normalize(bezier)?,
        });
    }
    Ok(spans)
}

/// Extract each nonzero span on the active parameter interval [U[p], U[n+1]].
/// Open, clamped and periodically extended input use the same local blossom.
/// Each extracted net retains its original interval and an exact integral gauge.
pub fn extract_active_bezier(
    degree: usize,
    knots: &[i64],
    controls: &[[i64; 3]],
) -> Result<Vec<BezierSpan>, &'static str> {
    if !(1..=3).contains(&degree)
        || controls.len() <= degree
        || controls.len().checked_add(degree + 1) != Some(knots.len())
    {
        return Err("unsupported sketch degree or knot/control count");
    }
    if controls.iter().any(|p| p[2] <= 0)
        || knots.windows(2).any(|p| p[1] < p[0])
        || knots.windows(degree + 2).any(|p| p[0] == p[degree + 1])
    {
        return Err("invalid sketch knot sequence or nonpositive weight");
    }
    let mut spans = Vec::new();
    for k in degree..controls.len() {
        let (lo, hi) = (knots[k], knots[k + 1]);
        if lo == hi {
            continue;
        }
        let initial: Vec<HPoint> = controls[k - degree..=k]
            .iter()
            .map(|p| p.map(|c| Ratio::new(i128::from(c), 1).expect("integer rational")))
            .collect();
        let mut net = Vec::with_capacity(degree + 1);
        for upper_count in 0..=degree {
            let mut row = initial.clone();
            // B_j is the blossom at p-j lower and j upper knot arguments.
            for r in 1..=degree {
                let t = if r <= degree - upper_count { lo } else { hi };
                for j in (r..=degree).rev() {
                    let i = k - degree + j;
                    let alpha = Ratio::new(
                        knot_difference(t, knots[i]),
                        knot_difference(knots[i + degree + 1 - r], knots[i]),
                    )?;
                    row[j] = blend_checked(row[j], row[j - 1], alpha)?;
                }
            }
            net.push(row[degree]);
        }
        spans.push(BezierSpan {
            interval: [lo, hi],
            controls: normalize(&net)?,
        });
    }
    if spans.is_empty() {
        return Err("empty active sketch parameter interval");
    }
    Ok(spans)
}
