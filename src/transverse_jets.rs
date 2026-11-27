//! Strong physical gradient matching for independent straight extrusions.
//! On a vertical face, at fixed axial Bernstein index, all cleared gradient
//! equations are polynomials in one parameter. Exact unisolvent evaluation
//! generates their coefficient row space, including the endpoints.
use crate::{
    constraint_system::Ratio as Q,
    incidence::{BlockSide, Interface},
    inherited_space::{ratio_div, Section},
    spline_constraints::{add, Result, Row},
};

fn mul(a: Q, b: Q) -> Result<Q> {
    a.checked_mul(b)
}
fn diff_product(a: Q, b: Q, c: Q, d: Q) -> Result<Q> {
    mul(a, b)?.checked_sub(mul(c, d)?)
}
fn bernstein(p: usize, t: Q) -> Result<(Vec<Q>, Vec<Q>)> {
    let s = Q::one().checked_sub(t)?;
    let mut prev = vec![Q::one()];
    for level in 1..=p {
        let mut next = vec![Q::zero(); level + 1];
        for (i, &x) in prev.iter().enumerate() {
            next[i] = next[i].checked_add(mul(s, x)?)?;
            next[i + 1] = next[i + 1].checked_add(mul(t, x)?)?;
        }
        if level == p {
            let mut der = vec![Q::zero(); p + 1];
            let deg = Q::new(p as i128, 1)?;
            for (i, &x) in prev.iter().enumerate() {
                let v = mul(x, deg)?;
                der[i] = der[i].checked_sub(v)?;
                der[i + 1] = der[i + 1].checked_add(v)?;
            }
            return Ok((next, der));
        }
        prev = next;
    }
    Err("degree zero transverse direction")
}

fn section_gradients(s: &Section, side: &BlockSide, t: Q) -> Result<Vec<[Q; 2]>> {
    let t_fixed = Q::new(i128::from(side.upper), 1)?;
    let u = if side.axis == 0 { t_fixed } else { t };
    let v = if side.axis == 1 { t_fixed } else { t };
    let (bu, du) = bernstein(s.d[0], u)?;
    let (bv, dv) = bernstein(s.d[1], v)?;
    let mut h = [[Q::zero(); 3]; 3]; // point, u derivative, v derivative (X,Y,W)
    let mut funcs = Vec::new();
    for i in 0..=s.d[0] {
        for j in 0..=s.d[1] {
            let p = s.block.controls[2 * (i * (s.d[1] + 1) + j)];
            let b = [mul(bu[i], bv[j])?, mul(du[i], bv[j])?, mul(bu[i], dv[j])?];
            let mut weighted = [Q::zero(); 3];
            for (k, &f) in b.iter().enumerate() {
                for (c, &x) in [p[0], p[1], p[3]].iter().enumerate() {
                    h[k][c] = h[k][c].checked_add(mul(f, Q::new(x.into(), 1)?)?)?;
                }
                weighted[k] = mul(f, Q::new(p[3].into(), 1)?)?;
            }
            funcs.push(weighted);
        }
    }
    let w = h[0][2];
    let xu = diff_product(h[1][0], w, h[0][0], h[1][2])?;
    let xv = diff_product(h[2][0], w, h[0][0], h[2][2])?;
    let yu = diff_product(h[1][1], w, h[0][1], h[1][2])?;
    let yv = diff_product(h[2][1], w, h[0][1], h[2][2])?;
    let det = diff_product(xu, yv, xv, yu)?;
    if det.numerator() <= 0 {
        return Err("source section does not have positive orientation");
    }
    funcs
        .into_iter()
        .map(|b| {
            let nu = diff_product(b[1], w, b[0], h[1][2])?;
            let nv = diff_product(b[2], w, b[0], h[2][2])?;
            Ok([
                ratio_div(diff_product(nu, yv, nv, yu)?, det)?,
                ratio_div(diff_product(nv, xu, nu, xv)?, det)?,
            ])
        })
        .collect()
}

pub(crate) fn append(
    link: &Interface,
    a: &Section,
    b: &Section,
    ids: &[usize],
    rows: &mut Vec<Row>,
) -> Result<()> {
    let o = link.orientation;
    if o.swap || o.flip_second || a.shear != b.shear || a.z != b.z || a.d[2] != b.d[2] {
        return Err("unsupported transverse axial orientation or interval");
    }
    let pa = a.d[1 - link.first.axis];
    let pb = b.d[1 - link.second.axis];
    // Each rational map's cleared Jacobian entries have degree <= 2p.
    // Numerator and determinant of a physical gradient each have degree <= 4p.
    let bound = 4 * (pa + pb);
    // Probe endpoints and low-height rationals first. They define the same
    // polynomial row space and avoid gratuitous denominator growth in elimination.
    let mut probes = vec![Q::zero(), Q::one()];
    for d in 2..=bound + 1 {
        for k in 1..d {
            let t = Q::new(k as i128, d as i128)?;
            if t.denominator() == d as i128 {
                probes.push(t);
            }
        }
        if probes.len() > bound {
            break;
        }
    }
    let mut constraints = Vec::new();
    for t in probes.into_iter().take(bound + 1) {
        let a_grad = section_gradients(a, &link.first, t)?;
        let b_grad = section_gradients(
            b,
            &link.second,
            if o.flip_first {
                Q::one().checked_sub(t)?
            } else {
                t
            },
        )?;
        for w in 0..=a.d[2] {
            for c in 0..2 {
                let mut row = Row::new();
                for (s, grad, sign) in [(a, &a_grad, Q::one()), (b, &b_grad, Q::new(-1, 1)?)] {
                    for (i, g) in grad.iter().enumerate() {
                        let at = [i / (s.d[1] + 1), i % (s.d[1] + 1), w];
                        add(&mut row, ids[s.flat(at)], mul(sign, g[c])?)?;
                    }
                }
                if !row.is_empty() {
                    constraints.push(row);
                }
            }
        }
    }
    // Reduce this one interface before intersecting constraints at edges and
    // vertices. The local relations typically have far smaller numerators
    // and supports than the redundant physical point equations.
    let mut local_ids: Vec<_> = constraints.iter().flat_map(|r| r.keys().copied()).collect();
    local_ids.sort_unstable();
    local_ids.dedup();
    let local = constraints
        .iter()
        .map(|r| {
            r.iter()
                .map(|(&i, &c)| (local_ids.binary_search(&i).unwrap(), c))
                .collect()
        })
        .collect::<Vec<Row>>();
    let (lift, free) = crate::spline_constraints::conditioned_kernel(local_ids.len(), &local)?;
    for (i, r) in lift.iter().enumerate() {
        if free.contains(&i) {
            continue;
        }
        let mut equation = Row::from([(local_ids[i], Q::one())]);
        for (&j, &c) in r {
            add(&mut equation, local_ids[free[j]], Q::zero().checked_sub(c)?)?;
        }
        rows.push(equation);
    }
    Ok(())
}
