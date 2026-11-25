//! Q2 elasticity with the volumetric strain projected onto a discontinuous
//! trilinear pressure space. Local pressure variables are statically condensed.
//! Positive recovered mean stress denotes tension.
use crate::{
    analysis_space::AnalysisSpace,
    elasticity::{assemble_lame, basis_geometry, gauss, Stiffness},
};
use std::collections::BTreeMap;

const N: usize = 81;
struct ElementProjection {
    nodes: [usize; 27],
    /// (M^-1 B) in local component-interleaved displacement order.
    recovery: Box<[[f64; N]; 8]>,
}
pub struct ProjectedElasticity {
    pub stiffness: Stiffness,
    bulk_modulus: f64,
    elements: BTreeMap<String, ElementProjection>,
}

fn pressure_basis([u, v, w]: [f64; 3]) -> [f64; 8] {
    let basis = [[1. - u, u], [1. - v, v], [1. - w, w]];
    std::array::from_fn(|i| basis[0][i / 4] * basis[1][i / 2 % 2] * basis[2][i % 2])
}

impl ProjectedElasticity {
    /// Q1 coefficients of the projected mean stress κ div u on this element.
    pub fn mean_stress_coefficients(&self, id: &str, u: &[f64]) -> Result<[f64; 8], &'static str> {
        let elem = self
            .elements
            .get(id)
            .ok_or("unknown projected stress block")?;
        if u.len() != self.stiffness.size() {
            return Err("invalid projected stress displacement");
        }
        let local: [f64; N] = std::array::from_fn(|i| u[3 * elem.nodes[i / 3] + i % 3]);
        if local.iter().any(|x| !x.is_finite()) {
            return Err("nonfinite projected stress displacement");
        }
        let values = std::array::from_fn(|j| {
            self.bulk_modulus
                * elem.recovery[j]
                    .iter()
                    .zip(local)
                    .map(|(a, b)| a * b)
                    .sum::<f64>()
        });
        if values.iter().any(|v| !v.is_finite()) {
            return Err("nonfinite recovered mean stress");
        }
        Ok(values)
    }
    /// Evaluate reconstructed hydrostatic stress at a point in the closed chart.
    pub fn mean_stress_at(coefficients: [f64; 8], uvw: [f64; 3]) -> Result<f64, &'static str> {
        if uvw.iter().any(|v| !(0.0..=1.0).contains(v))
            || coefficients.iter().any(|p| !p.is_finite())
        {
            return Err("invalid pressure evaluation");
        }
        Ok(pressure_basis(uvw)
            .iter()
            .zip(coefficients)
            .map(|(a, b)| a * b)
            .sum())
    }
}

/// Assemble fully integrated deviatoric stiffness and the positive local Schur complement.
/// Quadrature order is shared by the pressure mass, divergence and deviatoric terms.
pub fn assemble_projected_bulk(
    space: &AnalysisSpace,
    young: f64,
    poisson: f64,
    order: usize,
) -> Result<ProjectedElasticity, &'static str> {
    if !young.is_finite() || young <= 0. || !(-1.0 < poisson && poisson < 0.5) {
        return Err("invalid projected isotropic material");
    }
    let mu = young / (2. * (1. + poisson));
    let bulk = young / (3. * (1. - 2. * poisson));
    if !bulk.is_finite() || !mu.is_finite() {
        return Err("nonfinite projected moduli");
    }
    let quad = gauss(order)?;
    if quad.len() < 2 {
        return Err("pressure mass requires at least two nodes per axis");
    }
    let mut stiffness = assemble_lame(space, -2. * mu / 3., mu, order)?;
    let mut elements = BTreeMap::new();
    for b in &space.blocks {
        let mut mass = [[0.; 8]; 8];
        let mut coupling = Box::new([[0.; N]; 8]);
        for &(u, wu) in &quad {
            for &(v, wv) in &quad {
                for &(w, ww) in &quad {
                    let t = [u, v, w];
                    let s = basis_geometry(b, t)?;
                    let p = pressure_basis(t);
                    let dv = wu * wv * ww * s.determinant;
                    for (i, pi) in p.iter().enumerate() {
                        let weight = dv * pi;
                        for (j, pj) in p.iter().enumerate().take(i + 1) {
                            mass[i][j] += weight * pj;
                        }
                        for (node, g) in s.gradients.iter().enumerate() {
                            for (c, gc) in g.iter().enumerate() {
                                coupling[i][3 * node + c] += weight * gc;
                            }
                        }
                    }
                }
            }
        }
        // Lower Cholesky M=LLᵀ and forward solution F=L^-1 B.
        for i in 0..8 {
            for j in 0..=i {
                let mut v = mass[i][j];
                for (&a, &b) in mass[i][..j].iter().zip(&mass[j][..j]) {
                    v -= a * b;
                }
                mass[i][j] = if i == j {
                    if !v.is_finite() || v <= 0. {
                        return Err("singular local Q1 pressure mass");
                    }
                    v.sqrt()
                } else {
                    v / mass[j][j]
                };
            }
            for col in 0..N {
                let mut v = coupling[i][col];
                for j in 0..i {
                    v -= mass[i][j] * coupling[j][col];
                }
                coupling[i][col] = v / mass[i][i];
            }
        }
        for (i, &ii) in b.dofs.iter().enumerate() {
            let row = stiffness.offsets[ii]..stiffness.offsets[ii + 1];
            for (j, &jj) in b.dofs.iter().enumerate() {
                let pos = row.start
                    + stiffness.columns[row.clone()]
                        .binary_search(&jj)
                        .map_err(|_| "missing deviatoric CSR block")?;
                for a in 0..3 {
                    for d in 0..3 {
                        let gram = coupling
                            .iter()
                            .map(|f| f[3 * i + a] * f[3 * j + d])
                            .sum::<f64>();
                        stiffness.values[pos][a][d] += bulk * gram;
                    }
                }
            }
        }
        // Backward solution Lᵀ T=F for pressure recovery after displacement solve.
        for i in (0..8).rev() {
            for col in 0..N {
                let mut v = coupling[i][col];
                for j in i + 1..8 {
                    v -= mass[j][i] * coupling[j][col];
                }
                coupling[i][col] = v / mass[i][i];
            }
        }
        let element = ElementProjection {
            nodes: b
                .dofs
                .as_slice()
                .try_into()
                .map_err(|_| "invalid projected DOF count")?,
            recovery: coupling,
        };
        if elements.insert(b.id.clone(), element).is_some() {
            return Err("duplicate projected elasticity block");
        }
    }
    stiffness.validate()?;
    Ok(ProjectedElasticity {
        stiffness,
        bulk_modulus: bulk,
        elements,
    })
}
