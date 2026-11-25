//! Exact local cuts for the left wing of an S3 radial sector.
//!
//! This module certifies only the supplied sector rows. In particular, it does
//! not establish that the pocket's far X limit remains inside the bore void;
//! that predicate must be checked globally by the caller.

use crate::{constraint_system::Ratio, sketch_rational_event::simple_zero};

type Error = &'static str;
type Point = [Ratio; 3];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SliceStatus {
    InsideBand,
    OutsideBand,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParameterSlice {
    pub interval: [Ratio; 2],
    pub status: SliceStatus,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SectorPlan {
    /// All exact parameter cuts, including the sector endpoints 0 and 1.
    pub u_cuts: Vec<Ratio>,
    /// The connected part of the source sector whose ordinate lies in the
    /// pocket band. `None` means the band does not meet this sector.
    pub remove_interval: Option<[Ratio; 2]>,
    /// The radial guide lies at this exact X coordinate.
    pub guide_x: Ratio,
    /// Whether the outer radial chart is removed by the pocket.
    pub omit_outer: bool,
    /// Whether the inner radial chart is removed by the pocket.
    pub omit_inner: bool,
    /// This sector is attached to the stock's xmin boundary.
    pub wall_xmin: bool,
    /// Parameter intervals and their exact ordinate-band classification.
    pub slices: Vec<ParameterSlice>,
}

fn divide(a: Ratio, b: Ratio) -> Result<Ratio, Error> {
    if b.numerator() == 0 {
        return Err("S3 sector row has a zero homogeneous weight");
    }
    a.checked_mul(Ratio::new(b.denominator(), b.numerator())?)
}

fn ordinate_coefficients(row: &[Point], level: Ratio) -> Result<Vec<Ratio>, Error> {
    row.iter()
        .map(|p| p[1].checked_sub(level.checked_mul(p[2])?))
        .collect()
}

fn mix(a: Ratio, b: Ratio, t: Ratio) -> Result<Ratio, Error> {
    a.checked_add(b.checked_sub(a)?.checked_mul(t)?)
}

fn split(row: &[Point], t: Ratio) -> Result<(Vec<Point>, Vec<Point>), Error> {
    if row.len() < 2 || row.len() > 4 {
        return Err("S3 sector rows must be linear through cubic");
    }
    let mut work = row.to_vec();
    let mut left = vec![work[0]];
    let mut right = vec![*work.last().ok_or("empty S3 row")?];
    while work.len() > 1 {
        work = work
            .windows(2)
            .map(|pair| {
                Ok([
                    mix(pair[0][0], pair[1][0], t)?,
                    mix(pair[0][1], pair[1][1], t)?,
                    mix(pair[0][2], pair[1][2], t)?,
                ])
            })
            .collect::<Result<_, Error>>()?;
        left.push(work[0]);
        right.push(*work.last().ok_or("empty S3 split row")?);
    }
    right.reverse();
    Ok((left, right))
}

fn restrict(row: &[Point], lo: Ratio, hi: Ratio) -> Result<Vec<Point>, Error> {
    if lo.cmp_value(&Ratio::zero()).is_lt()
        || hi.cmp_value(&Ratio::one()).is_gt()
        || !lo.cmp_value(&hi).is_lt()
    {
        return Err("S3 parameter interval must be a nonempty subset of [0,1]");
    }
    let left = if lo == Ratio::zero() {
        row.to_vec()
    } else {
        split(row, lo)?.1
    };
    if hi == Ratio::one() {
        return Ok(left);
    }
    let local_hi = divide(hi.checked_sub(lo)?, Ratio::one().checked_sub(lo)?)?;
    Ok(split(&left, local_hi)?.0)
}

fn exact_roots(coefficients: &[Ratio]) -> Result<Vec<Ratio>, Error> {
    if !(2..=4).contains(&coefficients.len()) {
        return Err("S3 sector rows must be linear through cubic");
    }
    let mut roots = Vec::new();
    let mut work = coefficients.to_vec();
    let endpoint_root = work[0] == Ratio::zero() || work.last() == Some(&Ratio::zero());
    if work[0] == Ratio::zero() && work.last() == Some(&Ratio::zero()) {
        return Err("S3 ordinate cut has coincident endpoint roots");
    }
    if work[0] == Ratio::zero() {
        roots.push(Ratio::zero());
        work.remove(0);
    }
    if work.last() == Some(&Ratio::zero()) {
        roots.push(Ratio::one());
        work.pop();
    }
    if work.len() >= 2 {
        let signs: Vec<_> = work
            .iter()
            .map(|c| c.numerator().signum())
            .filter(|&sign| sign != 0)
            .collect();
        let variations = signs.windows(2).filter(|pair| pair[0] != pair[1]).count();
        if variations > 1 || (endpoint_root && variations > 0) {
            return Err("S3 ordinate cut has multiple or uncertified roots in one sector span");
        }
        if variations == 1 {
            roots.push(
                simple_zero(&work)
                    .ok_or("S3 ordinate cut lacks an exact simple rational root certificate")?,
            );
        }
    }
    roots.sort_by(|a, b| a.cmp_value(b));
    roots.dedup();
    Ok(roots)
}

fn eval_scalar(coefficients: &[Ratio], t: Ratio) -> Result<Ratio, Error> {
    let mut row = coefficients.to_vec();
    while row.len() > 1 {
        row = row
            .windows(2)
            .map(|pair| mix(pair[0], pair[1], t))
            .collect::<Result<_, _>>()?;
    }
    Ok(row[0])
}

fn inside_band(row: &[Point], ymin: Ratio, ymax: Ratio, t: Ratio) -> Result<bool, Error> {
    let low = eval_scalar(&ordinate_coefficients(row, ymin)?, t)?;
    let high = eval_scalar(&ordinate_coefficients(row, ymax)?, t)?;
    Ok(low.cmp_value(&Ratio::zero()).is_gt() && high.cmp_value(&Ratio::zero()).is_lt())
}

fn binomial(n: usize, k: usize) -> i128 {
    let k = k.min(n - k);
    (0..k).fold(1_i128, |value, i| value * (n - i) as i128 / (i + 1) as i128)
}

fn bernstein_product(a: &[Ratio], b: &[Ratio]) -> Result<Vec<Ratio>, Error> {
    let degree_a = a.len() - 1;
    let degree_b = b.len() - 1;
    let degree = degree_a + degree_b;
    let mut product = vec![Ratio::zero(); degree + 1];
    for (i, &left) in a.iter().enumerate() {
        for (j, &right) in b.iter().enumerate() {
            let factor = Ratio::new(
                binomial(degree_a, i) * binomial(degree_b, j),
                binomial(degree, i + j),
            )?;
            product[i + j] =
                product[i + j].checked_add(left.checked_mul(right)?.checked_mul(factor)?)?;
        }
    }
    Ok(product)
}

fn same_rational_ordinate(a: &[Point], b: &[Point]) -> Result<bool, Error> {
    let ay: Vec<_> = a.iter().map(|p| p[1]).collect();
    let aw: Vec<_> = a.iter().map(|p| p[2]).collect();
    let by: Vec<_> = b.iter().map(|p| p[1]).collect();
    let bw: Vec<_> = b.iter().map(|p| p[2]).collect();
    Ok(bernstein_product(&ay, &bw)? == bernstein_product(&by, &aw)?)
}

/// Conservatively attest that `pocket_xmax` lies inside the simple S2 bore on
/// the requested ordinate band. The caller must supply every right-wing source
/// row pair for that band, after validating the simple-hole S2 topology. Each
/// outer row must lie on the same stock xmax side; each inner row must be
/// strictly to the right of the pocket limit. Monotone row intervals must
/// cover the full band, and every band endpoint crossing must have an exact
/// simple rational certificate.
pub fn certify_xmax_void(
    outer_rows: &[Vec<Point>],
    inner_rows: &[Vec<Point>],
    pocket_xmax: Ratio,
    ymin: Ratio,
    ymax: Ratio,
) -> Result<(), Error> {
    if outer_rows.is_empty()
        || outer_rows.len() != inner_rows.len()
        || !ymin.cmp_value(&ymax).is_lt()
    {
        return Err("S3 global bore check requires paired rows and a nonempty ordinate band");
    }
    let mut stock_xmax: Option<Ratio> = None;
    let mut y_intervals = Vec::new();
    for (outer, inner) in outer_rows.iter().zip(inner_rows) {
        if !(2..=4).contains(&outer.len()) || outer.len() != inner.len() {
            return Err("S3 global bore rows must be matching linear-through-cubic spans");
        }
        if outer.iter().chain(inner).any(|p| p[2].numerator() <= 0) {
            return Err("S3 global bore rows require positive homogeneous weights");
        }
        if !same_rational_ordinate(outer, inner)? {
            return Err("S3 global bore rows do not have identical rational ordinates");
        }
        let x0 = divide(outer[0][0], outer[0][2])?;
        if outer.iter().any(|p| divide(p[0], p[2]).ok() != Some(x0))
            || stock_xmax.is_some_and(|x| x != x0)
        {
            return Err(
                "S3 global bore input must contain only the right-wing rows on one stock side",
            );
        }
        stock_xmax = Some(x0);
        for p in inner {
            let margin = p[0].checked_sub(pocket_xmax.checked_mul(p[2])?)?;
            if !margin.cmp_value(&Ratio::zero()).is_gt() {
                return Err("pocket xmax is not strictly inside the global right bore flank");
            }
        }

        let ordinates = inner
            .iter()
            .map(|p| divide(p[1], p[2]))
            .collect::<Result<Vec<_>, _>>()?;
        let ascending = ordinates.windows(2).all(|p| !p[0].cmp_value(&p[1]).is_gt());
        let descending = ordinates.windows(2).all(|p| !p[0].cmp_value(&p[1]).is_lt());
        if !ascending && !descending {
            return Err("S3 global bore rows must be monotone in ordinate");
        }
        let (low, high) = if ascending {
            (
                *ordinates.first().ok_or("empty S3 ordinate row")?,
                *ordinates.last().ok_or("empty S3 ordinate row")?,
            )
        } else {
            (
                *ordinates.last().ok_or("empty S3 ordinate row")?,
                *ordinates.first().ok_or("empty S3 ordinate row")?,
            )
        };
        if low.cmp_value(&ymax).is_lt() && high.cmp_value(&ymin).is_gt() {
            for level in [ymin, ymax] {
                let coeff = ordinate_coefficients(inner, level)?;
                let _ = exact_roots(&coeff)?;
            }
            y_intervals.push([low, high]);
        }
    }
    y_intervals.sort_by(|a, b| a[0].cmp_value(&b[0]));
    let mut covered = ymin;
    for [lo, hi] in y_intervals {
        if hi.cmp_value(&covered).is_le() {
            continue;
        }
        if lo.cmp_value(&covered).is_gt() {
            return Err("right-wing source rows do not cover the complete pocket ordinate band");
        }
        if hi.cmp_value(&covered).is_gt() {
            covered = hi;
        }
        if !covered.cmp_value(&ymax).is_lt() {
            return Ok(());
        }
    }
    Err("right-wing source rows do not cover the complete pocket ordinate band")
}

/// Plan exact band cuts for one sector on the stock's left side.
///
/// The outer and inner rows retain the same source parameter direction and
/// must represent a left-facing horizontal ray. `x_max` containment in the
/// bore void is deliberately a global caller-side check.
pub fn plan_sector(
    outer: &[Point],
    inner: &[Point],
    stock_xmin: Ratio,
    pocket_xmin: Ratio,
    pocket_xmax: Ratio,
    ymin: Ratio,
    ymax: Ratio,
) -> Result<SectorPlan, Error> {
    if !(2..=4).contains(&outer.len()) || outer.len() != inner.len() {
        return Err("S3 sector rows must have matching linear-through-cubic degree");
    }
    if !ymin.cmp_value(&ymax).is_lt()
        || !pocket_xmin.cmp_value(&pocket_xmax).is_lt()
        || !stock_xmin.cmp_value(&pocket_xmin).is_le()
    {
        return Err("S3 pocket bounds are degenerate or outside the left stock side");
    }
    let wall_xmin = outer
        .iter()
        .all(|p| p[2].numerator() > 0 && divide(p[0], p[2]).ok() == Some(stock_xmin));
    if !wall_xmin {
        return Err("S3 sector outer row is not on stock xmin");
    }
    if outer.iter().chain(inner).any(|p| p[2].numerator() <= 0) {
        return Err("S3 sector is not a matching left-wing sector");
    }
    if !same_rational_ordinate(outer, inner)? {
        return Err("S3 sector rows do not have identical rational ordinates");
    }

    let low_roots = exact_roots(&ordinate_coefficients(inner, ymin)?)?;
    let high_roots = exact_roots(&ordinate_coefficients(inner, ymax)?)?;
    let mut cuts = vec![Ratio::zero(), Ratio::one()];
    cuts.extend(low_roots);
    cuts.extend(high_roots);
    cuts.sort_by(|a, b| a.cmp_value(b));
    cuts.dedup();

    let mut slices = Vec::new();
    for pair in cuts.windows(2) {
        if pair[0] == pair[1] {
            continue;
        }
        let midpoint = pair[0]
            .checked_add(pair[1])?
            .checked_mul(Ratio::new(1, 2)?)?;
        slices.push(ParameterSlice {
            interval: [pair[0], pair[1]],
            status: if inside_band(inner, ymin, ymax, midpoint)? {
                SliceStatus::InsideBand
            } else {
                SliceStatus::OutsideBand
            },
        });
    }
    let inside: Vec<_> = slices
        .iter()
        .filter(|s| s.status == SliceStatus::InsideBand)
        .collect();
    if inside.len() > 1 {
        return Err("S3 pocket band is disconnected on one left-wing sector span");
    }
    let remove_interval = inside.first().map(|s| s.interval);

    if let Some([lo, hi]) = remove_interval {
        let section = restrict(inner, lo, hi)?;
        let xmin_values: Vec<_> = section
            .iter()
            .map(|p| p[0].checked_sub(pocket_xmin.checked_mul(p[2])?))
            .collect::<Result<_, Error>>()?;
        if xmin_values
            .iter()
            .any(|v| v.cmp_value(&Ratio::zero()).is_le())
        {
            return Err("pocket xmin does not lie strictly inside the left-wing ray");
        }
        let xmax_values: Vec<_> = section
            .iter()
            .map(|p| pocket_xmax.checked_mul(p[2])?.checked_sub(p[0]))
            .collect::<Result<_, Error>>()?;
        if xmax_values
            .iter()
            .any(|v| v.cmp_value(&Ratio::zero()).is_le())
        {
            return Err("pocket xmax is not to the right of the certified left bore abscissa");
        }
    }

    let at_stock_edge = pocket_xmin == stock_xmin;
    Ok(SectorPlan {
        u_cuts: cuts,
        remove_interval,
        guide_x: pocket_xmin,
        omit_outer: at_stock_edge,
        omit_inner: true,
        wall_xmin: pocket_xmin != stock_xmin,
        slices,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn n(value: i128) -> Ratio {
        Ratio::new(value, 1).unwrap()
    }

    fn rows() -> (Vec<Point>, Vec<Point>) {
        (
            vec![[n(0), n(0), n(1)], [n(0), n(40), n(1)]],
            vec![[n(60), n(0), n(1)], [n(60), n(40), n(1)]],
        )
    }

    #[test]
    fn finds_exact_band_cuts_and_certifies_local_left_ray() {
        let (outer, inner) = rows();
        let plan = plan_sector(&outer, &inner, n(0), n(55), n(75), n(10), n(30)).unwrap();
        assert_eq!(
            plan.u_cuts,
            [
                n(0),
                Ratio::new(1, 4).unwrap(),
                Ratio::new(3, 4).unwrap(),
                n(1)
            ]
        );
        assert_eq!(
            plan.remove_interval,
            Some([Ratio::new(1, 4).unwrap(), Ratio::new(3, 4).unwrap()])
        );
        assert_eq!(plan.guide_x, n(55));
        assert!(!plan.omit_outer && plan.omit_inner && plan.wall_xmin);
        assert_eq!(plan.slices.len(), 3);
        assert_eq!(plan.slices[1].status, SliceStatus::InsideBand);
    }

    #[test]
    fn accepts_band_bounded_by_sector_endpoints_for_open_side_case() {
        let (outer, inner) = rows();
        let plan = plan_sector(&outer, &inner, n(0), n(0), n(75), n(0), n(40)).unwrap();
        assert_eq!(plan.remove_interval, Some([n(0), n(1)]));
        assert!(plan.omit_outer && plan.omit_inner);
        assert!(!plan.wall_xmin);
    }

    #[test]
    fn rejects_nonmatching_or_nonleft_sector_rows() {
        let (mut outer, inner) = rows();
        outer[1][1] = n(39);
        assert!(plan_sector(&outer, &inner, n(0), n(55), n(75), n(10), n(30)).is_err());
        let (mut outer, inner) = rows();
        outer[0][0] = n(1);
        assert!(plan_sector(&outer, &inner, n(0), n(55), n(75), n(10), n(30)).is_err());
    }

    #[test]
    fn isolates_simple_root_when_a_bernstein_control_is_zero() {
        assert_eq!(
            exact_roots(&[n(1), n(0), n(-1)]).unwrap(),
            [Ratio::new(1, 2).unwrap()]
        );
    }

    #[test]
    fn compares_rational_ordinate_functions_not_only_control_poles() {
        let a = vec![[n(0), n(0), n(1)], [n(1), n(2), n(2)]];
        // Same rational y trace under a constant homogeneous gauge factor.
        let b = vec![[n(0), n(0), n(2)], [n(2), n(4), n(4)]];
        assert!(same_rational_ordinate(&a, &b).unwrap());
        let c = vec![[n(0), n(0), n(2)], [n(2), n(5), n(4)]];
        assert!(!same_rational_ordinate(&a, &c).unwrap());
    }

    #[test]
    fn globally_certifies_right_flank_and_ordinate_band_coverage() {
        let outer = vec![[n(100), n(0), n(1)], [n(100), n(40), n(1)]];
        let inner = vec![[n(80), n(0), n(1)], [n(80), n(40), n(1)]];
        assert!(certify_xmax_void(
            std::slice::from_ref(&outer),
            std::slice::from_ref(&inner),
            n(75),
            n(10),
            n(30)
        )
        .is_ok());
        let too_far = vec![[n(70), n(0), n(1)], [n(70), n(40), n(1)]];
        assert!(certify_xmax_void(&[outer], &[too_far], n(75), n(10), n(30)).is_err());
    }
}
