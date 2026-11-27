//! Scalar harmonic Ritz benchmark on the three extruded charts of an L section.
//! Displacements are independent of height. The outer trace is mass-projected;
//! stiffness integration is exact for the piecewise bicubic polynomial basis.
use crate::{
    bezier::basis,
    constrained_elasticity::solve_constrained,
    constraint_system::{resolve, Ratio},
    elasticity::{gauss, Stiffness},
    json::{object, Json},
    thb_reentrant_edge::{IntegrationLeaf, ReentrantHierarchy},
};
use std::collections::BTreeMap;
type Result<T> = std::result::Result<T, &'static str>;
type Rows = Vec<BTreeMap<usize, f64>>;
fn scalar(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn chart(p: usize, [u, v]: [f64; 2]) -> [f64; 2] {
    match p {
        0 => [-u, -v],
        1 => [v, -u],
        2 => [-v, u],
        _ => unreachable!(),
    }
}
fn exact([x, y]: [f64; 2]) -> (f64, [f64; 2]) {
    let theta = (-y).atan2(x).rem_euclid(std::f64::consts::TAU);
    let r = x.hypot(y);
    let k = 2. / 3.;
    let (st, ct) = (k * theta).sin_cos();
    let scale = k * r.powf(k - 1.);
    (
        r.powf(k) * st,
        [
            scale * (st * x / r + ct * y / r),
            scale * (st * y / r - ct * x / r),
        ],
    )
}
fn widths(e: &IntegrationLeaf) -> [f64; 2] {
    e.uv.map(|[a, b]| b - a)
}
fn xy(e: &IntegrationLeaf, t: [f64; 2]) -> [f64; 2] {
    chart(
        e.patch,
        std::array::from_fn(|k| e.uv[k][0] + widths(e)[k] * t[k]),
    )
}
/// Exterior section boundary. The two interior chart seams are eliminated in the basis.
fn boundary(e: &IntegrationLeaf) -> Vec<(usize, usize)> {
    let mut faces = Vec::new();
    for a in 0..2 {
        if e.uv[a][1] == 1. {
            faces.push((a, 1))
        }
        if e.uv[a][0] == 0. && ((e.patch == 1 && a == 0) || (e.patch == 2 && a == 1)) {
            faces.push((a, 0))
        }
    }
    faces
}
/// Apply the kernel PCG solver in one active scalar component. The independent
/// affine lift contributes the known nonzero Dirichlet values before reduction.
fn solve_scalar(rows: &Rows, rhs: &[f64], active: &[usize]) -> Result<(Vec<f64>, f64)> {
    let n = active.len();
    let mut index = vec![usize::MAX; rows.len()];
    for (j, &i) in active.iter().enumerate() {
        index[i] = j;
    }
    let mut k = Stiffness {
        offsets: vec![0],
        columns: Vec::new(),
        values: Vec::new(),
        integrated_volume: 0.,
    };
    let mut f = vec![0.; 3 * n];
    for (j, &i) in active.iter().enumerate() {
        f[3 * j] = rhs[i];
        for (&t, &v) in &rows[i] {
            if index[t] != usize::MAX {
                let mut m = [[0.; 3]; 3];
                m[0][0] = v;
                k.values.push(m);
                k.columns.push(index[t]);
            }
        }
        k.offsets.push(k.columns.len());
    }
    let free = resolve(n, &[], &[])?;
    let fixed = resolve(
        n,
        &[],
        &(0..n).map(|i| (i, Ratio::zero())).collect::<Vec<_>>(),
    )?;
    let solution = solve_constrained(&k, [&free, &fixed, &fixed], &f, 1e-12, 20000)?;
    Ok((
        solution
            .displacement_mm
            .as_chunks::<3>()
            .0
            .iter()
            .map(|r| r[0])
            .collect(),
        solution.relative_residual,
    ))
}
fn project_trace(leaves: &[IntegrationLeaf], count: usize) -> Result<(Vec<f64>, Vec<usize>)> {
    let mut mass = vec![BTreeMap::new(); count];
    let mut load = vec![0.; count];
    let q = gauss(16)?;
    for e in leaves {
        for (a, b) in boundary(e) {
            let t_axis = 1 - a;
            for &(t, w) in &q {
                let mut loc = [t; 2];
                loc[a] = b as f64;
                let (bern, _) = basis(3, t);
                let wt = w * widths(e)[t_axis];
                let value = exact(xy(e, loc)).0;
                let mut trace = vec![0.; e.dofs.len()];
                for (j, &nj) in bern.iter().enumerate() {
                    let ij = if a == 0 { 12 * b + j } else { 4 * j + 3 * b };
                    for (v, c) in trace.iter_mut().zip(&e.bernstein[ij]) {
                        *v += nj * c;
                    }
                }
                for (i, &id) in e.dofs.iter().enumerate().filter(|(i, _)| trace[*i] != 0.) {
                    load[id] += wt * trace[i] * value;
                    for (j, &jd) in e.dofs.iter().enumerate().filter(|(j, _)| trace[*j] != 0.) {
                        *mass[id].entry(jd).or_default() += wt * trace[i] * trace[j];
                    }
                }
            }
        }
    }
    let active: Vec<_> = (0..count).filter(|&i| !mass[i].is_empty()).collect();
    let (trace, residual) = solve_scalar(&mass, &load, &active)?;
    if residual > 1e-10 {
        return Err("boundary projection residual too large");
    }
    let mut field = vec![0.; count];
    for (&i, x) in active.iter().zip(trace) {
        field[i] = x;
    }
    Ok((field, (0..count).filter(|&i| mass[i].is_empty()).collect()))
}
/// Tensor Bernstein stiffness on the unit square. All cells are square in
/// physical coordinates, so the Jacobian cancels the two gradient scales.
fn stiffness(leaves: &[IntegrationLeaf], count: usize) -> Result<Rows> {
    let q = gauss(4)?;
    let mut k = [[0.; 16]; 16];
    for &(u, wu) in &q {
        for &(v, wv) in &q {
            let (bu, du) = basis(3, u);
            let (bv, dv) = basis(3, v);
            let grads: [[f64; 2]; 16] =
                std::array::from_fn(|i| [du[i / 4] * bv[i % 4], bu[i / 4] * dv[i % 4]]);
            for (i, a) in grads.iter().enumerate() {
                for (j, b) in grads.iter().enumerate() {
                    k[i][j] += wu * wv * (a[0] * b[0] + a[1] * b[1]);
                }
            }
        }
    }
    let mut rows = vec![BTreeMap::new(); count];
    for e in leaves {
        let [hx, hy] = widths(e);
        if hx <= 0. || hx != hy {
            return Err("expected positive dyadic square");
        }
        let mut temp = vec![vec![0.; e.dofs.len()]; 16];
        for (i, t) in temp.iter_mut().enumerate() {
            for (j, ej) in e.bernstein.iter().enumerate() {
                for (a, &c) in t.iter_mut().zip(ej) {
                    *a += k[i][j] * c;
                }
            }
        }
        for (j, ej) in e.bernstein.iter().enumerate() {
            for (&i, &c) in e.dofs.iter().zip(ej).filter(|(_, v)| **v != 0.) {
                for (&t, &x) in e.dofs.iter().zip(&temp[j]) {
                    *rows[i].entry(t).or_default() += c * x;
                }
            }
        }
    }
    Ok(rows)
}
fn value(c: &[f64; 16], loc: [f64; 2]) -> (f64, [f64; 2]) {
    let (bu, du) = basis(3, loc[0]);
    let (bv, dv) = basis(3, loc[1]);
    let mut f = 0.;
    let mut grad = [0.; 2];
    for (i, x) in c.iter().enumerate() {
        f += x * bu[i / 4] * bv[i % 4];
        grad[0] += x * du[i / 4] * bv[i % 4];
        grad[1] += x * bu[i / 4] * dv[i % 4];
    }
    (f, grad)
}
fn integrate(leaves: &[IntegrationLeaf], fields: &[[f64; 16]], order: usize) -> Result<[f64; 3]> {
    let q = gauss(order)?;
    let mut e2 = 0.;
    let mut g2 = 0.;
    let mut u2 = 0.;
    let mut du2 = 0.;
    let mut be = 0.;
    let mut bn = 0.;
    for (e, c) in leaves.iter().zip(fields) {
        let [hx, hy] = widths(e);
        let mut point = |t: [f64; 2], weight: f64| {
            let (u, g) = exact(xy(e, t));
            let (uh, d) = value(c, t);
            let gh = chart(e.patch, [d[0] / hx, d[1] / hy]);
            e2 += weight * ((g[0] - gh[0]).powi(2) + (g[1] - gh[1]).powi(2));
            g2 += weight * (g[0] * g[0] + g[1] * g[1]);
            u2 += weight * (u - uh).powi(2);
            du2 += weight * u * u;
        };
        let touches = e.uv.iter().all(|i| i[0] == 0.);
        for &(s, ws) in &q {
            for &(t, wt) in &q {
                let w = ws * wt * hx * hy;
                if touches {
                    // Split the corner square into two triangles. Duffy radius
                    // r=s³ resolves the fractional power in the energy integrand.
                    let r = s.powi(3);
                    let jac = 3. * s * s * r;
                    point([r, r * t], w * jac);
                    point([r * t, r], w * jac);
                } else {
                    point([s, t], w);
                }
            }
        }
        for (a, b) in boundary(e) {
            for &(s, w) in &q {
                let mut t = [s; 2];
                t[a] = b as f64;
                let u = exact(xy(e, t)).0;
                let uh = value(c, t).0;
                be += w * widths(e)[1 - a] * (u - uh).powi(2);
                bn += w * widths(e)[1 - a] * u * u;
            }
        }
    }
    Ok([(e2 / g2).sqrt(), (u2 / du2).sqrt(), (be / bn).sqrt()])
}
/// Solve the planar scalar Ritz problem and integrate the error directly on
/// its physical cells. The non-polynomial trace is projected independently at
/// each level. Area, dimensionality and the complete quadrature history are reported.
pub fn solve(space: &ReentrantHierarchy) -> Result<Json> {
    let leaves = space.integration_leaves()?;
    let count = space.active_scalar_dofs;
    let (mut field, free) = project_trace(&leaves, count)?;
    let k = stiffness(&leaves, count)?;
    let rhs = k
        .iter()
        .map(|r| -r.iter().map(|(&j, &c)| field[j] * c).sum::<f64>())
        .collect::<Vec<_>>();
    let (solution, residual) = solve_scalar(&k, &rhs, &free)?;
    for (&i, x) in free.iter().zip(solution) {
        field[i] = x
    }
    let fields = leaves
        .iter()
        .map(|e| {
            std::array::from_fn(|j| {
                e.dofs
                    .iter()
                    .zip(&e.bernstein[j])
                    .map(|(&i, v)| field[i] * v)
                    .sum()
            })
        })
        .collect::<Vec<_>>();
    let mut history = Vec::new();
    let mut previous: Option<[f64; 2]> = None;
    let mut stable = false;
    let mut norms = [0.; 3];
    for q in [8, 12, 16, 24, 32] {
        norms = integrate(&leaves, &fields, q)?;
        if let Some([h, l]) = previous {
            stable = (norms[0] / h - 1_f64).abs() <= 2e-3 && (norms[1] / l - 1_f64).abs() <= 5e-4;
        }
        previous = Some([norms[0], norms[1]]);
        history.push(object([
            ("order", scalar(q)),
            ("relative_H1_seminorm", scalar(norms[0])),
            ("relative_L2", scalar(norms[1])),
            ("relative_boundary_L2", scalar(norms[2])),
        ]));
    }
    let area = leaves
        .iter()
        .map(|e| widths(e).iter().product::<f64>())
        .sum::<f64>();
    if (area - 3.).abs() > 1e-13 || residual > 1e-10 || !stable {
        return Err("singular benchmark failed area, residual or quadrature acceptance");
    }
    Ok(object([
        ("active_scalar_dofs", scalar(count)),
        ("free_scalar_dofs", scalar(free.len())),
        ("maximum_level", scalar(space.maximum_level)),
        ("leaves", scalar(leaves.len())),
        ("area", scalar(area)),
        ("relative_free_residual", scalar(residual)),
        ("relative_H1_seminorm", scalar(norms[0])),
        ("relative_L2", scalar(norms[1])),
        ("relative_boundary_L2", scalar(norms[2])),
        ("quadrature_history", Json::Array(history)),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn polar_branch_matches_both_notch_sides_and_all_chart_gradients() {
        assert!(exact([0.4, 0.]).0.abs() < 1e-14);
        assert!(exact([0., 0.4]).0.abs() < 1e-14);
        for p in 0..3 {
            let x = chart(p, [0.23, 0.47]);
            let (_, g) = exact(x);
            for a in 0..2 {
                let mut lo = x;
                let mut hi = x;
                lo[a] -= 1e-6;
                hi[a] += 1e-6;
                assert!(((exact(hi).0 - exact(lo).0) / 2e-6 - g[a]).abs() < 1e-8);
            }
        }
    }
}
