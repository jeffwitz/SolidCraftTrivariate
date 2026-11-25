//! Exact scalar field constraints for rectangular restrictions of rational
//! Bézier faces. Coefficients are transformed only after a homogeneous trace
//! identity has been established for the whole subface.

use crate::incidence::Orientation;

#[derive(Debug, Clone)]
pub struct FaceNet {
    pub degrees: [usize; 2],
    pub controls: Vec<[i64; 4]>,
}

/// Wide exact face net used internally by repeatedly restricted adaptive
/// charts. The public `FaceNet` remains the compact i64 interchange type.
#[derive(Debug, Clone)]
pub(crate) struct FaceNetWide {
    pub degrees: [usize; 2],
    pub controls: Vec<[i128; 4]>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Span {
    Full,
    Lower,
    Upper,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    /// (master coefficient index, positive numerator) in the master's native layout.
    pub terms: Vec<(usize, i128)>,
    pub denominator: i128,
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}

fn validate(net: &FaceNet) -> Result<(), &'static str> {
    if net.degrees.iter().any(|&d| !(1..=2).contains(&d))
        || net.controls.len() != (net.degrees[0] + 1) * (net.degrees[1] + 1)
        || net.controls.iter().any(|p| p[3] <= 0)
    {
        return Err("invalid face net");
    }
    Ok(())
}

// The rows share a common power-of-two denominator. It cancels in both
// the projective trace comparison and the normalized field coefficients.
fn restriction(p: usize, out: usize, span: Span) -> Result<Vec<Vec<i128>>, &'static str> {
    if out < p {
        return Err("slave degree does not contain master degree");
    }
    let rows = match (p, span) {
        (1, Span::Full) => vec![vec![1, 0], vec![0, 1]],
        (1, Span::Lower) => vec![vec![2, 0], vec![1, 1]],
        (1, Span::Upper) => vec![vec![1, 1], vec![0, 2]],
        (2, Span::Full) => vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]],
        (2, Span::Lower) => vec![vec![4, 0, 0], vec![2, 2, 0], vec![1, 2, 1]],
        (2, Span::Upper) => vec![vec![1, 2, 1], vec![0, 2, 2], vec![0, 0, 4]],
        _ => return Err("unsupported trace degree"),
    };
    if p == out {
        return Ok(rows);
    }
    Ok(vec![
        rows[0].iter().map(|x| 2 * x).collect(),
        rows[0].iter().zip(&rows[1]).map(|(a, b)| a + b).collect(),
        rows[1].iter().map(|x| 2 * x).collect(),
    ])
}

/// The orientation maps the slave net to the master's two parameter axes.
/// Returned rows are indexed in the original, unpermuted slave storage.
pub fn nested_constraints(
    master: &FaceNet,
    slave: &FaceNet,
    spans: [Span; 2],
    o: Orientation,
) -> Result<Vec<Row>, &'static str> {
    validate(master)?;
    validate(slave)?;
    let master = FaceNetWide {
        degrees: master.degrees,
        controls: master.controls.iter().map(|p| p.map(i128::from)).collect(),
    };
    let slave = FaceNetWide {
        degrees: slave.degrees,
        controls: slave.controls.iter().map(|p| p.map(i128::from)).collect(),
    };
    nested_constraints_wide(&master, &slave, spans, o)
}

pub(crate) fn nested_constraints_wide(
    master: &FaceNetWide,
    slave: &FaceNetWide,
    spans: [Span; 2],
    o: Orientation,
) -> Result<Vec<Row>, &'static str> {
    validate_wide(master)?;
    validate_wide(slave)?;
    let [p, q] = master.degrees;
    let [r, s] = if o.swap {
        [slave.degrees[1], slave.degrees[0]]
    } else {
        slave.degrees
    };
    let u = restriction(p, r, spans[0])?;
    let v = restriction(q, s, spans[1])?;
    let mut expected = vec![[0_i128; 4]; slave.controls.len()];
    let mut constraints = Vec::with_capacity(expected.len());
    for (i, ur) in u.iter().enumerate() {
        for (j, vr) in v.iter().enumerate() {
            let a = if o.flip_first { r - i } else { i };
            let b = if o.flip_second { s - j } else { j };
            let slave_index = if o.swap {
                b * (slave.degrees[1] + 1) + a
            } else {
                a * (s + 1) + b
            };
            let mut terms = vec![];
            for (m, &um) in ur.iter().enumerate() {
                for (n, &vn) in vr.iter().enumerate() {
                    let factor = um * vn;
                    if factor == 0 {
                        continue;
                    }
                    let index = m * (q + 1) + n;
                    let point = master.controls[index];
                    for (dest, value) in expected[slave_index].iter_mut().zip(point) {
                        *dest = dest
                            .checked_add(
                                factor
                                    .checked_mul(value)
                                    .ok_or("restricted face coefficient exceeds i128")?,
                            )
                            .ok_or("restricted face coefficient exceeds i128")?;
                    }
                    terms.push((
                        index,
                        factor
                            .checked_mul(point[3])
                            .ok_or("constraint coefficient exceeds i128")?,
                    ));
                }
            }
            let common = terms.iter().map(|(_, weight)| *weight).fold(0, gcd);
            for (_, weight) in &mut terms {
                *weight /= common;
            }
            let denominator = terms.iter().map(|(_, w)| w).sum();
            constraints.push((slave_index, Row { terms, denominator }));
        }
    }
    // Net-wide normalization excludes independent pointwise rescaling.
    // Each axis contributes at most a factor four, so i64 inputs and their
    // exact restricted net fit i128 even at the limits of signed i64.
    let scale = expected.iter().flatten().copied().fold(0, gcd);
    let other = slave.controls.iter().flatten().copied().fold(0, gcd);
    if expected
        .iter()
        .flatten()
        .zip(slave.controls.iter().flatten())
        .any(|(e, &a)| e / scale != a / other)
    {
        return Err("incompatible restricted homogeneous trace");
    }
    constraints.sort_by_key(|(index, _)| *index);
    Ok(constraints.into_iter().map(|(_, row)| row).collect())
}

fn validate_wide(net: &FaceNetWide) -> Result<(), &'static str> {
    if net.degrees.iter().any(|&d| !(1..=2).contains(&d))
        || net.controls.len() != (net.degrees[0] + 1) * (net.degrees[1] + 1)
        || net.controls.iter().any(|p| p[3] <= 0)
    {
        return Err("invalid face net");
    }
    Ok(())
}
