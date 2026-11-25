//! Recover simple rational level crossings of low-degree Bernstein polynomials.
use crate::constraint_system::Ratio;

/// When the Bernstein sign sequence has one variation and both endpoints are
/// nonzero, a single simple root lies inside (0,1). Floating point only suggests
/// a parameter. Exact de Casteljau evaluation certifies each returned root.
pub(crate) fn simple_zero(c: &[Ratio]) -> Option<Ratio> {
    if !(2..=4).contains(&c.len()) || c[0].numerator() == 0 || c.last()?.numerator() == 0 {
        return None;
    }
    let mut sign = c[0].numerator().signum();
    let mut changes = 0;
    for coefficient in c {
        let next = coefficient.numerator().signum();
        if next != 0 && next != sign {
            changes += 1;
            sign = next;
        }
    }
    if changes != 1 {
        return None;
    }
    let values: Vec<_> = c
        .iter()
        .map(|v| v.numerator() as f64 / v.denominator() as f64)
        .collect();
    let mut lo = 0.;
    let mut hi = 1.;
    for _ in 0..60 {
        let mid = (lo + hi) / 2.;
        let mut row = [0.; 4];
        row[..c.len()].copy_from_slice(&values);
        for n in (1..c.len()).rev() {
            for i in 0..n {
                row[i] = (1. - mid) * row[i] + mid * row[i + 1];
            }
        }
        if row[0] == 0. {
            lo = mid;
            hi = mid;
            break;
        }
        if row[0].is_sign_positive() == values[0].is_sign_positive() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    // Continued-fraction convergents of the isolated numerical estimate.
    let mut x = (lo + hi) / 2.;
    let (mut n_old, mut n) = (0_i128, 1_i128);
    let (mut d_old, mut d) = (1_i128, 0_i128);
    for _ in 0..32 {
        let a = x.floor();
        if !a.is_finite() || a > 1_000_000. {
            break;
        }
        (n_old, n) = (n, a as i128 * n + n_old);
        (d_old, d) = (d, a as i128 * d + d_old);
        if d > 1_000_000 {
            break;
        }
        if n > 0 && n < d {
            let t = Ratio::new(n, d).ok()?;
            if exact_value(c, t) == Some(Ratio::zero()) {
                return Some(t);
            }
        }
        let fraction = x - a;
        if fraction == 0. {
            break;
        }
        x = 1. / fraction;
    }
    None
}

fn exact_value(c: &[Ratio], t: Ratio) -> Option<Ratio> {
    let complement = Ratio::one().checked_sub(t).ok()?;
    let mut row = c.to_vec();
    for n in (1..c.len()).rev() {
        for i in 0..n {
            row[i] = row[i]
                .checked_mul(complement)
                .ok()?
                .checked_add(row[i + 1].checked_mul(t).ok()?)
                .ok()?;
        }
    }
    Some(row[0])
}

#[cfg(test)]
mod tests {
    use super::*;
    fn integer(n: i128) -> Ratio {
        Ratio::new(n, 1).unwrap()
    }

    #[test]
    fn nondyadic_rational_event_survives_a_change_of_knot_interval() {
        let polynomial = [-3, -1, 4, 12].map(integer);
        // The degree-three Bernstein controls represent -3+6t+9t².
        assert_eq!(simple_zero(&polynomial), Some(Ratio::new(1, 3).unwrap()));
        assert_eq!(
            simple_zero(&[-123, 874].map(integer)),
            Some(Ratio::new(123, 997).unwrap())
        );
    }

    #[test]
    fn approximate_irrational_or_ambiguous_event_is_never_admitted() {
        // 2t²-1 on the unit interval has an irrational zero.
        assert_eq!(simple_zero(&[-1, -1, 1].map(integer)), None);
        assert_eq!(simple_zero(&[1, -2, 1].map(integer)), None);
        assert_eq!(simple_zero(&[0, 1].map(integer)), None);
    }
}
