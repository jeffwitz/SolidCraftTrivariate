//! Three-dimensional isotropic elasticity on an audited rational Q2 space.
//! Physical coordinates are mm, Young's modulus is MPa, and forces are N.

use super::analysis_space::{AnalysisBlock, AnalysisSpace};
use std::collections::BTreeMap;

pub type Vector = [f64; 3];
pub type Matrix = [Vector; 3];

fn dot(a: Vector, b: Vector) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: Vector, b: Vector) -> Vector {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Gauss-Legendre abscissas and weights on `[0,1]`.
pub fn gauss(order: usize) -> Result<Vec<(f64, f64)>, &'static str> {
    if !(1..=32).contains(&order) {
        return Err("quadrature order must be in 1..=32");
    }
    fn poly(n: usize, x: f64) -> (f64, f64) {
        let mut a = 1.;
        let mut b = x;
        for k in 2..=n {
            let fk = k as f64;
            (a, b) = (b, ((2. * fk - 1.) * x * b - (fk - 1.) * a) / fk);
        }
        (b, n as f64 * (a - x * b) / (1. - x * x))
    }
    let mut result = Vec::with_capacity(order);
    for k in 0..order {
        let mut x = (std::f64::consts::PI * (k as f64 + 0.75) / (order as f64 + 0.5)).cos();
        for _ in 0..32 {
            let (p, dp) = poly(order, x);
            let next = x - p / dp;
            if (next - x).abs() <= 2. * f64::EPSILON {
                x = next;
                break;
            }
            x = next;
        }
        let (_, dp) = poly(order, x);
        result.push(((1. - x) / 2., 1. / ((1. - x * x) * dp * dp)));
    }
    Ok(result)
}

pub struct BasisSample {
    pub position: Vector,
    pub values: [f64; 27],
    pub gradients: [Vector; 27],
    pub jacobian_columns: Matrix,
    pub determinant: f64,
}

pub fn basis_geometry(block: &AnalysisBlock, uvw: Vector) -> Result<BasisSample, &'static str> {
    if block.controls.len() != 27 || uvw.iter().any(|v| !(0.0..=1.0).contains(v)) {
        return Err("invalid Q2 block or parameters");
    }
    let basis = uvw.map(|x| [(1. - x) * (1. - x), 2. * x * (1. - x), x * x]);
    let derivative = uvw.map(|x| [-2. * (1. - x), 2. - 4. * x, 2. * x]);
    let mut values = [0.; 27];
    let mut grads = [[0.; 3]; 27];
    let mut denominator = 0.;
    let mut dden = [0.; 3];
    let mut position = [0.; 3];
    let mut coords = [[0.; 3]; 27];
    for (n, &p) in block.controls.iter().enumerate() {
        if p[3] <= 0 {
            return Err("nonpositive weight");
        }
        let [i, j, k] = [n / 9, n / 3 % 3, n % 3];
        let weight = p[3] as f64;
        let [bu, bv, bw] = [basis[0][i], basis[1][j], basis[2][k]];
        values[n] = weight * bu * bv * bw;
        grads[n] = [
            weight * derivative[0][i] * bv * bw,
            weight * bu * derivative[1][j] * bw,
            weight * bu * bv * derivative[2][k],
        ];
        denominator += values[n];
        for a in 0..3 {
            coords[n][a] = p[a] as f64 / weight;
            position[a] += values[n] * coords[n][a];
            dden[a] += grads[n][a];
        }
    }
    for x in &mut position {
        *x /= denominator;
    }
    let mut jac = [[0.; 3]; 3];
    for n in 0..27 {
        values[n] /= denominator;
        for a in 0..3 {
            grads[n][a] = (grads[n][a] - values[n] * dden[a]) / denominator;
            for c in 0..3 {
                jac[a][c] += grads[n][a] * (coords[n][c] - position[c]);
            }
        }
    }
    let cof = [
        cross(jac[1], jac[2]),
        cross(jac[2], jac[0]),
        cross(jac[0], jac[1]),
    ];
    let det = dot(jac[0], cof[0]);
    if !det.is_finite() || det <= 0. {
        return Err("nonpositive sampled Jacobian");
    }
    for g in &mut grads {
        let local = *g;
        *g = std::array::from_fn(|c| {
            (local[0] * cof[0][c] + local[1] * cof[1][c] + local[2] * cof[2][c]) / det
        });
    }
    Ok(BasisSample {
        position,
        values,
        gradients: grads,
        jacobian_columns: jac,
        determinant: det,
    })
}

pub struct Stiffness {
    pub offsets: Vec<usize>,
    pub columns: Vec<usize>,
    pub values: Vec<Matrix>,
    pub integrated_volume: f64,
}

impl Stiffness {
    pub fn size(&self) -> usize {
        3 * self.offsets.len().saturating_sub(1)
    }
    /// Validate block CSR before using a matrix provided by a caller.
    pub fn validate(&self) -> Result<(), &'static str> {
        let rows = self.offsets.len().checked_sub(1).ok_or("missing CSR offsets")?;
        if self.offsets.first() != Some(&0) || self.offsets.last() != Some(&self.columns.len())
            || self.columns.len() != self.values.len()
            || self.offsets.windows(2).any(|pair| pair[0] > pair[1])
            || self.columns.iter().any(|&col| col >= rows)
            || self.values.iter().flatten().flatten().any(|value| !value.is_finite()) {
            return Err("invalid block CSR matrix");
        }
        Ok(())
    }
    /// Evaluate f-Ku with compensated products and additions. This is intended
    /// for the final residual check, where internal forces nearly cancel.
    pub fn equilibrium_residual(&self, u:&[f64], f:&[f64]) -> Result<Vec<f64>, &'static str> {
        self.validate()?;
        if u.len()!=self.size() || f.len()!=u.len() || u.iter().chain(f).any(|x|!x.is_finite()) {
            return Err("invalid equilibrium vectors");
        }
        let mut residual=f.to_vec();
        for (i,components) in residual.as_chunks_mut::<3>().0.iter_mut().enumerate() {
            for (c,out) in components.iter_mut().enumerate() {
                let mut sum=*out;
                let mut correction=0.;
                for k in self.offsets[i]..self.offsets[i+1] {
                    let j=3*self.columns[k];
                    for a in 0..3 {
                        let coefficient=-self.values[k][c][a];
                        let term=coefficient*u[j+a];
                        let product_error=coefficient.mul_add(u[j+a],-term);
                        let next=sum+term;
                        let addition_error=if sum.abs()>=term.abs() {(sum-next)+term} else {(term-next)+sum};
                        correction+=addition_error+product_error;
                        sum=next;
                    }
                }
                *out=sum+correction;
                if !out.is_finite() {return Err("nonfinite equilibrium residual")}
            }
        }
        Ok(residual)
    }

    pub fn apply(&self, x: &[f64]) -> Result<Vec<f64>, &'static str> {
        if x.len() != self.size() {
            return Err("wrong displacement vector length");
        }
        let mut y = vec![0.; x.len()];
        for (i, out) in y.as_chunks_mut::<3>().0.iter_mut().enumerate() {
            for k in self.offsets[i]..self.offsets[i + 1] {
                let j = self.columns[k];
                let v = [x[3 * j], x[3 * j + 1], x[3 * j + 2]];
                for (c, row) in self.values[k].iter().enumerate() {
                    out[c] += dot(*row, v);
                }
            }
        }
        Ok(y)
    }
}

pub fn assemble(
    space: &AnalysisSpace,
    young: f64,
    poisson: f64,
    order: usize,
) -> Result<Stiffness, &'static str> {
    if !young.is_finite() || young <= 0. || !(-1.0 < poisson && poisson < 0.5) {
        return Err("invalid isotropic material");
    }
    let mu = young / (2. * (1. + poisson));
    let lambda = young * poisson / ((1. + poisson) * (1. - 2. * poisson));
    assemble_lame(space, lambda, mu, order)
}

/// Integrate the isotropic bilinear form from its two Lamé coefficients.
pub(crate) fn assemble_lame(
    space: &AnalysisSpace,
    lambda: f64,
    mu: f64,
    order: usize,
) -> Result<Stiffness, &'static str> {
    let points = gauss(order)?;
    let mut rows = vec![BTreeMap::<usize, Matrix>::new(); space.ndof];
    let mut volume = 0.;
    for block in &space.blocks {
        if block.dofs.len() != 27 || block.dofs.iter().any(|&d| d >= space.ndof) {
            return Err("invalid local DOF map");
        }
        let mut element = vec![[[0.; 3]; 3]; 27 * 27];
        for &(u, wu) in &points {
            for &(v, wv) in &points {
                for &(w, ww) in &points {
                    let sample = basis_geometry(block, [u, v, w])?;
                    let measure = sample.determinant * wu * wv * ww;
                    volume += measure;
                    for (i, &gi) in sample.gradients.iter().enumerate() {
                        for (j, &gj) in sample.gradients.iter().enumerate().take(i + 1) {
                            let diag = mu * dot(gi, gj);
                            let entry = &mut element[27 * i + j];
                            for a in 0..3 {
                                for b in 0..3 {
                                    let elastic = lambda * gi[a] * gj[b]
                                        + mu * gi[b] * gj[a]
                                        + if a == b { diag } else { 0. };
                                    entry[a][b] += measure * elastic;
                                }
                            }
                        }
                    }
                }
            }
        }
        for (i, &ii) in block.dofs.iter().enumerate() {
            for (j, &jj) in block.dofs.iter().enumerate().take(i + 1) {
                let entry = element[i * 27 + j];
                let target = rows[ii].entry(jj).or_default();
                for a in 0..3 {
                    for b in 0..3 {
                        target[a][b] += entry[a][b];
                    }
                }
                if i != j {
                    let target = rows[jj].entry(ii).or_default();
                    for a in 0..3 {
                        for b in 0..3 {
                            target[a][b] += entry[b][a];
                        }
                    }
                }
            }
        }
    }
    let mut offsets = vec![0];
    let mut columns = vec![];
    let mut values = vec![];
    for row in rows {
        for (j, e) in row {
            columns.push(j);
            values.push(e);
        }
        offsets.push(columns.len());
    }
    Ok(Stiffness {
        offsets,
        columns,
        values,
        integrated_volume: volume,
    })
}
