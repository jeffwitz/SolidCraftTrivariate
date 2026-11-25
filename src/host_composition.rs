//! Exact composition of supported trivariate blocks with piecewise-affine hosts.
//!
//! The input block's geometric validity is deliberately outside this module. This
//! module checks only the representation and the stated monotone-x composition domain.

use crate::bracket::BracketBlock;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RationalParameter {
    pub numerator: i128,
    pub denominator: i128,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AffinePiece {
    pub matrix: [[i64; 3]; 3],
    pub translation: [i64; 3],
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PiecewiseAffineHost {
    /// Strictly increasing integer x coordinates separating the affine pieces.
    pub breaks: Vec<i64>,
    /// Must contain exactly `breaks.len() + 1` maps.
    pub pieces: Vec<AffinePiece>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HostNodeStatus {
    NotCrossed,
    AtParameterStart,
    AtParameterEnd,
    InternalRational(RationalParameter),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HostNodeEvent {
    pub plane_x: i64,
    pub status: HostNodeStatus,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HostRepresentation {
    Unsplit,
    SplitAtRationalNodes,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComposedBlock {
    pub block: BracketBlock,
    pub parameter_interval: [RationalParameter; 2],
    pub host_piece_index: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HostComposition {
    pub blocks: Vec<ComposedBlock>,
    pub node_events: Vec<HostNodeEvent>,
    pub representation: HostRepresentation,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HostCompositionError {
    InvalidHost,
    NonPositiveDeterminant,
    DiscontinuousHost,
    InvalidBlock,
    NonPositiveWeight,
    XVariesAcrossVOrW,
    NonMonotoneX,
    IrrationalPreimage,
    IntegerCapacityExceeded,
}

type H = [i128; 4];

fn add(a: i128, b: i128) -> Result<i128, HostCompositionError> {
    a.checked_add(b)
        .ok_or(HostCompositionError::IntegerCapacityExceeded)
}
fn sub(a: i128, b: i128) -> Result<i128, HostCompositionError> {
    a.checked_sub(b)
        .ok_or(HostCompositionError::IntegerCapacityExceeded)
}
fn mul(a: i128, b: i128) -> Result<i128, HostCompositionError> {
    a.checked_mul(b)
        .ok_or(HostCompositionError::IntegerCapacityExceeded)
}
fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

impl RationalParameter {
    fn new(n: i128, d: i128) -> Result<Self, HostCompositionError> {
        if d == 0 {
            return Err(HostCompositionError::InvalidBlock);
        }
        let (n, d) = if d < 0 {
            (
                n.checked_neg()
                    .ok_or(HostCompositionError::IntegerCapacityExceeded)?,
                d.checked_neg()
                    .ok_or(HostCompositionError::IntegerCapacityExceeded)?,
            )
        } else {
            (n, d)
        };
        let g = i128::try_from(gcd(n.unsigned_abs(), d.unsigned_abs()))
            .map_err(|_| HostCompositionError::IntegerCapacityExceeded)?;
        Ok(Self {
            numerator: n / g,
            denominator: d / g,
        })
    }
}

fn determinant(m: &[[i64; 3]; 3]) -> Result<i128, HostCompositionError> {
    let m = m.map(|r| r.map(i128::from));
    let a = sub(mul(m[1][1], m[2][2])?, mul(m[1][2], m[2][1])?)?;
    let b = sub(mul(m[1][0], m[2][2])?, mul(m[1][2], m[2][0])?)?;
    let c = sub(mul(m[1][0], m[2][1])?, mul(m[1][1], m[2][0])?)?;
    add(sub(mul(m[0][0], a)?, mul(m[0][1], b)?)?, mul(m[0][2], c)?)
}

fn validate_host(host: &PiecewiseAffineHost) -> Result<(), HostCompositionError> {
    if host.pieces.len() != host.breaks.len() + 1
        || host.pieces.is_empty()
        || host.breaks.windows(2).any(|w| w[0] >= w[1])
    {
        return Err(HostCompositionError::InvalidHost);
    }
    for p in &host.pieces {
        if determinant(&p.matrix)? <= 0 {
            return Err(HostCompositionError::NonPositiveDeterminant);
        }
    }
    for (i, &k) in host.breaks.iter().enumerate() {
        let l = &host.pieces[i];
        let r = &host.pieces[i + 1];
        for out in 0..3 {
            for axis in 1..3 {
                if l.matrix[out][axis] != r.matrix[out][axis] {
                    return Err(HostCompositionError::DiscontinuousHost);
                }
            }
            let lv = i128::from(l.matrix[out][0])
                .checked_mul(i128::from(k))
                .and_then(|v| v.checked_add(i128::from(l.translation[out])));
            let rv = i128::from(r.matrix[out][0])
                .checked_mul(i128::from(k))
                .and_then(|v| v.checked_add(i128::from(r.translation[out])));
            if lv.is_none() || lv != rv {
                return Err(HostCompositionError::DiscontinuousHost);
            }
        }
    }
    Ok(())
}

fn block_shape(block: &BracketBlock) -> Result<[usize; 3], HostCompositionError> {
    if block.degrees != [1, 1, 1] && block.degrees != [2, 1, 1] {
        return Err(HostCompositionError::InvalidBlock);
    }
    let shape = block.degrees.map(|d| d + 1);
    if block.controls.len() != shape.iter().product::<usize>() {
        return Err(HostCompositionError::InvalidBlock);
    }
    if block.controls.iter().any(|p| p[3] <= 0) {
        return Err(HostCompositionError::NonPositiveWeight);
    }
    Ok(shape)
}
fn idx(shape: [usize; 3], i: usize, j: usize, k: usize) -> usize {
    (i * shape[1] + j) * shape[2] + k
}

fn endpoint_x(
    block: &BracketBlock,
    shape: [usize; 3],
    i: usize,
) -> Result<(i128, i128), HostCompositionError> {
    let first = block.controls[idx(shape, i, 0, 0)];
    for j in 0..shape[1] {
        for k in 0..shape[2] {
            let p = block.controls[idx(shape, i, j, k)];
            if p[0] != first[0] || p[3] != first[3] {
                return Err(HostCompositionError::XVariesAcrossVOrW);
            }
        }
    }
    Ok((i128::from(first[0]), i128::from(first[3])))
}

fn isqrt(n: u128) -> u128 {
    if n < 2 {
        return n;
    }
    let mut x = n;
    let mut y = x.div_ceil(2);
    while y < x {
        x = y;
        y = (x + n / x) / 2;
    }
    x
}
fn preimage(q: &[i128]) -> Result<RationalParameter, HostCompositionError> {
    if q.len() == 2 {
        let den = sub(q[0], q[1])?;
        return RationalParameter::new(q[0], den);
    }
    let a = add(sub(q[0], mul(2, q[1])?)?, q[2])?;
    let b = mul(2, sub(q[1], q[0])?)?;
    let c = q[0];
    if a == 0 {
        if b == 0 {
            return Err(HostCompositionError::InvalidBlock);
        }
        return RationalParameter::new(
            c.checked_neg()
                .ok_or(HostCompositionError::IntegerCapacityExceeded)?,
            b,
        );
    }
    let disc = sub(mul(b, b)?, mul(4, mul(a, c)?)?)?;
    if disc < 0 {
        return Err(HostCompositionError::InvalidBlock);
    }
    let root = isqrt(disc as u128);
    if root.checked_mul(root) != Some(disc as u128) {
        return Err(HostCompositionError::IrrationalPreimage);
    }
    let root = i128::try_from(root).map_err(|_| HostCompositionError::IntegerCapacityExceeded)?;
    let den = mul(2, a)?;
    let neg_b = b
        .checked_neg()
        .ok_or(HostCompositionError::IntegerCapacityExceeded)?;
    let candidates = [
        RationalParameter::new(sub(neg_b, root)?, den)?,
        RationalParameter::new(add(neg_b, root)?, den)?,
    ];
    candidates
        .into_iter()
        .find(|p| p.numerator > 0 && p.numerator < p.denominator)
        .ok_or(HostCompositionError::InvalidBlock)
}

fn restrict_u(
    block: &BracketBlock,
    shape: [usize; 3],
    a: RationalParameter,
    b: RationalParameter,
) -> Result<Vec<H>, HostCompositionError> {
    let degree = block.degrees[0];
    let den = a
        .denominator
        .checked_mul(b.denominator)
        .ok_or(HostCompositionError::IntegerCapacityExceeded)?;
    let an = mul(a.numerator, b.denominator)?;
    let bn = mul(b.numerator, a.denominator)?;
    let mut out = vec![[0; 4]; block.controls.len()];
    for j in 0..shape[1] {
        for k in 0..shape[2] {
            let ps: Vec<H> = (0..=degree)
                .map(|i| block.controls[idx(shape, i, j, k)].map(i128::from))
                .collect();
            let restricted: Vec<H> = if degree == 1 {
                let make = |n: i128| -> Result<H, HostCompositionError> {
                    let mut v = [0; 4];
                    for c in 0..4 {
                        v[c] = add(mul(sub(den, n)?, ps[0][c])?, mul(n, ps[1][c])?)?;
                    }
                    Ok(v)
                };
                vec![make(an)?, make(bn)?]
            } else {
                let mut controls = Vec::with_capacity(3);
                for (x, y) in [(an, an), (an, bn), (bn, bn)] {
                    let mut v = [0; 4];
                    for c in 0..4 {
                        let c1 = mul(2, sub(ps[1][c], ps[0][c])?)?;
                        let c2 = add(sub(ps[0][c], mul(2, ps[1][c])?)?, ps[2][c])?;
                        v[c] = add(
                            add(
                                mul(2, mul(ps[0][c], mul(den, den)?)?)?,
                                mul(c1, mul(den, add(x, y)?)?)?,
                            )?,
                            mul(2, mul(c2, mul(x, y)?)?)?,
                        )?;
                    }
                    controls.push(v);
                }
                controls
            };
            for i in 0..=degree {
                out[idx(shape, i, j, k)] = restricted[i];
            }
        }
    }
    Ok(out)
}

fn apply_piece(
    mut net: Vec<H>,
    piece: &AffinePiece,
) -> Result<Vec<[i64; 4]>, HostCompositionError> {
    let mut g = 0u128;
    for p in &mut net {
        let xyz = [p[0], p[1], p[2]];
        for r in 0..3 {
            let mut v = mul(i128::from(piece.translation[r]), p[3])?;
            for (c, value) in xyz.iter().enumerate() {
                v = add(v, mul(i128::from(piece.matrix[r][c]), *value)?)?;
            }
            p[r] = v;
        }
        for x in p {
            g = gcd(g, x.unsigned_abs());
        }
    }
    let g = i128::try_from(g.max(1)).map_err(|_| HostCompositionError::IntegerCapacityExceeded)?;
    net.into_iter()
        .map(|p| {
            let mut q = [0i64; 4];
            for c in 0..4 {
                q[c] = i64::try_from(p[c] / g)
                    .map_err(|_| HostCompositionError::IntegerCapacityExceeded)?;
            }
            Ok(q)
        })
        .collect()
}

/// Compose a rational `[1,1,1]` or `[2,1,1]` block with a continuous host.
pub fn compose_with_host(
    block: &BracketBlock,
    host: &PiecewiseAffineHost,
) -> Result<HostComposition, HostCompositionError> {
    validate_host(host)?;
    let shape = block_shape(block)?;
    let mut xs = Vec::new();
    for i in 0..shape[0] {
        xs.push(endpoint_x(block, shape, i)?);
    }
    for i in 0..xs.len() - 1 {
        if mul(xs[i].0, xs[i + 1].1)? >= mul(xs[i + 1].0, xs[i].1)? {
            return Err(HostCompositionError::NonMonotoneX);
        }
    }
    let start = xs[0];
    let mut events = Vec::new();
    let mut cuts = vec![RationalParameter {
        numerator: 0,
        denominator: 1,
    }];
    for &k in &host.breaks {
        let mut q = Vec::with_capacity(xs.len());
        for (x, w) in &xs {
            q.push(sub(*x, mul(i128::from(k), *w)?)?);
        }
        let status = if q[0] == 0 {
            HostNodeStatus::AtParameterStart
        } else if *q.last().unwrap() == 0 {
            HostNodeStatus::AtParameterEnd
        } else if q[0] < 0 && *q.last().unwrap() > 0 {
            let p = preimage(&q)?;
            cuts.push(p);
            HostNodeStatus::InternalRational(p)
        } else {
            HostNodeStatus::NotCrossed
        };
        events.push(HostNodeEvent { plane_x: k, status });
    }
    cuts.push(RationalParameter {
        numerator: 1,
        denominator: 1,
    });
    // Host breaks and strict x monotonicity already put internal roots in order.
    let mut result = Vec::new();
    let start_piece = host
        .breaks
        .iter()
        .filter(|&&k| start.0 >= i128::from(k) * start.1)
        .count();
    for (segment, pair) in cuts.windows(2).enumerate() {
        let (a, b) = (pair[0], pair[1]);
        let piece_index = start_piece + segment;
        let net = restrict_u(block, shape, a, b)?;
        let controls = apply_piece(net, &host.pieces[piece_index])?;
        let id = format!("{}_host_{}", block.id, piece_index);
        result.push(ComposedBlock {
            block: BracketBlock {
                id,
                degrees: block.degrees,
                controls,
            },
            parameter_interval: [a, b],
            host_piece_index: piece_index,
        });
    }
    Ok(HostComposition {
        blocks: result,
        node_events: events,
        representation: if cuts.len() > 2 {
            HostRepresentation::SplitAtRationalNodes
        } else {
            HostRepresentation::Unsplit
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn identity() -> AffinePiece {
        AffinePiece {
            matrix: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            translation: [0, 0, 0],
        }
    }
    #[test]
    fn affine_host_identity_preserves_controls() {
        let b = BracketBlock {
            id: "unit".into(),
            degrees: [1, 1, 1],
            controls: vec![
                [0, 0, 0, 1],
                [0, 0, 1, 1],
                [0, 1, 0, 1],
                [0, 1, 1, 1],
                [1, 0, 0, 1],
                [1, 0, 1, 1],
                [1, 1, 0, 1],
                [1, 1, 1, 1],
            ],
        };
        let h = PiecewiseAffineHost {
            breaks: vec![],
            pieces: vec![identity()],
        };
        let out = compose_with_host(&b, &h).unwrap();
        assert_eq!(out.blocks.len(), 1);
        assert_eq!(out.blocks[0].block.controls, b.controls);
    }
    #[test]
    fn rejects_host_with_nonpositive_determinant() {
        let b = BracketBlock {
            id: "x".into(),
            degrees: [1, 1, 1],
            controls: vec![[0, 0, 0, 1]; 8],
        };
        let h = PiecewiseAffineHost {
            breaks: vec![],
            pieces: vec![AffinePiece {
                matrix: [[0; 3]; 3],
                translation: [0; 3],
            }],
        };
        assert_eq!(
            compose_with_host(&b, &h),
            Err(HostCompositionError::NonPositiveDeterminant)
        );
    }
    #[test]
    fn degenerate_quadratic_root_is_solved_linearly() {
        assert_eq!(
            preimage(&[-1, 0, 1]).unwrap(),
            RationalParameter {
                numerator: 1,
                denominator: 2
            }
        );
    }
    #[test]
    fn affine_output_capacity_overflow_is_rejected() {
        let mut controls = Vec::new();
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    controls.push([i as i64, j as i64, k as i64, 1]);
                }
            }
        }
        let b = BracketBlock {
            id: "overflow".into(),
            degrees: [1, 1, 1],
            controls,
        };
        let h = PiecewiseAffineHost {
            breaks: vec![],
            pieces: vec![AffinePiece {
                matrix: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                translation: [0, i64::MAX, 0],
            }],
        };
        assert_eq!(
            compose_with_host(&b, &h),
            Err(HostCompositionError::IntegerCapacityExceeded)
        );
    }
    #[test]
    fn determinant_arithmetic_overflow_is_rejected() {
        let b = BracketBlock {
            id: "overflow".into(),
            degrees: [1, 1, 1],
            controls: vec![[0, 0, 0, 1]; 8],
        };
        let m = i64::MAX;
        let h = PiecewiseAffineHost {
            breaks: vec![],
            pieces: vec![AffinePiece {
                matrix: [[m, 0, 0], [0, m, 0], [0, 0, m]],
                translation: [0; 3],
            }],
        };
        assert_eq!(
            compose_with_host(&b, &h),
            Err(HostCompositionError::IntegerCapacityExceeded)
        );
    }
}
