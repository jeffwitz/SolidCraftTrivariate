//! Variable-size rational tensor-extracted elements and native elasticity.
//!
//! Axis extraction matrices use Bernstein rows and local spline-function
//! columns. Tensor products use the local ordering `(i * n_v + j) * n_w + k`.
//! Geometry and field evaluation are kept separate at the public sample
//! boundary, so a compact parent geometry can be integrated with a distinct
//! extracted field space.

use crate::constraint_system::Ratio;
use crate::elasticity::{gauss, Matrix, Stiffness, Vector};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq)]
pub struct ElementGeometrySample {
    pub position: Vector,
    /// Rows are parameter directions; columns are physical coordinates.
    pub jacobian_columns: Matrix,
    pub determinant: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ElementBasisSample {
    pub dofs: Vec<usize>,
    pub values: Vec<f64>,
    pub physical_gradients: Vec<Vector>,
}

/// Sparse univariate extraction matrix stored by Bernstein row. Coefficients
/// remain exact rationals; each row has strictly increasing local-function
/// columns and zero entries are omitted.
#[derive(Debug, Clone, PartialEq)]
pub struct AxisExtraction {
    pub degree: usize,
    pub basis_count: usize,
    pub row_offsets: Vec<usize>,
    pub columns: Vec<usize>,
    pub coefficients: Vec<Ratio>,
}

impl AxisExtraction {
    pub fn validate(&self) -> Result<(), &'static str> {
        if !(1..=3).contains(&self.degree)
            || self.basis_count == 0
            || self.row_offsets.len() != self.degree + 2
            || self.row_offsets.first() != Some(&0)
            || self.row_offsets.last() != Some(&self.columns.len())
            || self.columns.len() != self.coefficients.len()
            || self.row_offsets.windows(2).any(|r| r[0] > r[1])
            || self.columns.iter().any(|&i| i >= self.basis_count)
            || self
                .coefficients
                .iter()
                .any(|x| x.denominator() <= 0 || x.numerator() == 0)
        {
            return Err("invalid sparse axis extraction");
        }
        for row in 0..=self.degree {
            let cols = &self.columns[self.row_offsets[row]..self.row_offsets[row + 1]];
            if cols.windows(2).any(|pair| pair[0] >= pair[1]) {
                return Err("axis extraction columns are not strictly ordered");
            }
            let coefficients = &self.coefficients[self.row_offsets[row]..self.row_offsets[row + 1]];
            let sum = coefficients
                .iter()
                .try_fold(Ratio::zero(), |acc, &value| acc.checked_add(value))?;
            if sum != Ratio::one() || coefficients.iter().any(|value| value.numerator() < 0) {
                return Err("axis extraction row does not preserve polynomial partition");
            }
        }
        let mut present = vec![false; self.basis_count];
        for &column in &self.columns {
            present[column] = true;
        }
        if present.iter().any(|seen| !*seen) {
            return Err("axis extraction contains an inactive local function");
        }
        Ok(())
    }

    fn apply(&self, bernstein: &[f64], derivative: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let mut values = vec![0.; self.basis_count];
        let mut derivs = vec![0.; self.basis_count];
        for row in 0..=self.degree {
            for entry in self.row_offsets[row]..self.row_offsets[row + 1] {
                let col = self.columns[entry];
                let coefficient = self.coefficients[entry];
                let c = coefficient.numerator() as f64 / coefficient.denominator() as f64;
                values[col] += bernstein[row] * c;
                derivs[col] += derivative[row] * c;
            }
        }
        (values, derivs)
    }
}

/// A tensor-product element with independent rational geometry and field
/// weights. `geometry_controls` are homogeneous `(wx, wy, wz, w)` values.
#[derive(Debug, Clone, PartialEq)]
pub struct ExtractedElement {
    pub axes: [AxisExtraction; 3],
    pub dofs: Vec<usize>,
    pub field_weights: Vec<f64>,
    pub geometry_controls: Vec<[f64; 4]>,
}

impl ExtractedElement {
    pub fn validate(&self) -> Result<(), &'static str> {
        for axis in &self.axes {
            axis.validate()?;
        }
        let count = self.axes.iter().try_fold(1usize, |n, a| {
            n.checked_mul(a.basis_count)
                .ok_or("element basis size overflow")
        })?;
        if self.dofs.len() != count
            || self.field_weights.len() != count
            || self.geometry_controls.len() != count
            || self
                .field_weights
                .iter()
                .any(|w| !w.is_finite() || *w <= 0.)
            || self
                .geometry_controls
                .iter()
                .any(|p| p.iter().any(|x| !x.is_finite()) || p[3] <= 0.)
        {
            return Err("invalid extracted element data");
        }
        let mut ids = self.dofs.clone();
        ids.sort_unstable();
        if ids.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err("duplicate global DOF in extracted element");
        }
        Ok(())
    }

    /// Evaluate rational geometry and the rational field basis at one local
    /// point. The physical gradients use the inverse Jacobian of the geometry.
    pub fn sample(
        &self,
        uvw: Vector,
    ) -> Result<(ElementGeometrySample, ElementBasisSample), &'static str> {
        self.validate()?;
        self.sample_valid(uvw)
    }

    fn sample_valid(
        &self,
        uvw: Vector,
    ) -> Result<(ElementGeometrySample, ElementBasisSample), &'static str> {
        if uvw
            .iter()
            .any(|x| !x.is_finite() || !(0. ..=1.).contains(x))
        {
            return Err("element parameter outside the closed unit cube");
        }
        let mut axis_values = Vec::with_capacity(3);
        let mut axis_derivs = Vec::with_capacity(3);
        for (axis, &t) in self.axes.iter().zip(&uvw) {
            let (bernstein, derivative) = bernstein(axis.degree, t)?;
            let (values, derivs) = axis.apply(&bernstein, &derivative);
            axis_values.push(values);
            axis_derivs.push(derivs);
        }

        let dims = self.axes.each_ref().map(|a| a.basis_count);
        let n = self.dofs.len();
        let mut polynomial = vec![0.; n];
        let mut dpoly = vec![[0.; 3]; n];
        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    let q = (i * dims[1] + j) * dims[2] + k;
                    let a = axis_values[0][i];
                    let b = axis_values[1][j];
                    let c = axis_values[2][k];
                    polynomial[q] = a * b * c;
                    dpoly[q] = [
                        axis_derivs[0][i] * b * c,
                        a * axis_derivs[1][j] * c,
                        a * b * axis_derivs[2][k],
                    ];
                }
            }
        }

        let mut field_num = vec![0.; n];
        let mut field_dnum = vec![[0.; 3]; n];
        let mut field_den = 0.;
        let mut field_dden = [0.; 3];
        let mut geom_den = 0.;
        let mut geom_dden = [0.; 3];
        let mut geom_num = [0.; 3];
        let mut geom_dnum = [[0.; 3]; 3];
        for q in 0..n {
            let fw = self.field_weights[q];
            let gw = self.geometry_controls[q][3];
            field_num[q] = fw * polynomial[q];
            field_den += field_num[q];
            geom_den += gw * polynomial[q];
            for a in 0..3 {
                field_dnum[q][a] = fw * dpoly[q][a];
                field_dden[a] += field_dnum[q][a];
                geom_dden[a] += gw * dpoly[q][a];
            }
            for c in 0..3 {
                let x = self.geometry_controls[q][c] / gw;
                geom_num[c] += gw * polynomial[q] * x;
                for a in 0..3 {
                    geom_dnum[a][c] += gw * dpoly[q][a] * x;
                }
            }
        }
        if !field_den.is_finite() || field_den <= 0. || !geom_den.is_finite() || geom_den <= 0. {
            return Err("nonpositive extracted rational denominator");
        }
        let position = geom_num.map(|x| x / geom_den);
        let mut jac = [[0.; 3]; 3];
        for a in 0..3 {
            for c in 0..3 {
                jac[a][c] = (geom_dnum[a][c] - position[c] * geom_dden[a]) / geom_den;
            }
        }
        let cof = [
            cross(jac[1], jac[2]),
            cross(jac[2], jac[0]),
            cross(jac[0], jac[1]),
        ];
        let det = dot(jac[0], cof[0]);
        if !det.is_finite() || det <= 0. {
            return Err("nonpositive extracted element Jacobian");
        }

        let mut values = vec![0.; n];
        let mut param_grads = vec![[0.; 3]; n];
        for q in 0..n {
            values[q] = field_num[q] / field_den;
            for a in 0..3 {
                param_grads[q][a] = (field_dnum[q][a] - values[q] * field_dden[a]) / field_den;
            }
        }
        let physical_gradients = param_grads
            .iter()
            .map(|g| {
                [0, 1, 2].map(|c| (g[0] * cof[0][c] + g[1] * cof[1][c] + g[2] * cof[2][c]) / det)
            })
            .collect();
        Ok((
            ElementGeometrySample {
                position,
                jacobian_columns: jac,
                determinant: det,
            },
            ElementBasisSample {
                dofs: self.dofs.clone(),
                values,
                physical_gradients,
            },
        ))
    }
}

fn bernstein(degree: usize, t: f64) -> Result<(Vec<f64>, Vec<f64>), &'static str> {
    let values = match degree {
        1 => vec![1. - t, t],
        2 => vec![(1. - t).powi(2), 2. * t * (1. - t), t.powi(2)],
        3 => vec![
            (1. - t).powi(3),
            3. * t * (1. - t).powi(2),
            3. * t.powi(2) * (1. - t),
            t.powi(3),
        ],
        _ => return Err("extraction degree must be in 1..=3"),
    };
    let derivatives = match degree {
        1 => vec![-1., 1.],
        2 => vec![-2. * (1. - t), 2. - 4. * t, 2. * t],
        3 => vec![
            -3. * (1. - t).powi(2),
            3. * (1. - t) * (1. - 3. * t),
            3. * t * (2. - 3. * t),
            3. * t.powi(2),
        ],
        _ => return Err("extraction degree must be in 1..=3"),
    };
    Ok((values, derivatives))
}

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

pub struct NativeAssembly {
    pub stiffness: Stiffness,
    pub force: Vec<f64>,
    pub quadrature_points: usize,
    pub cell_count: usize,
    pub quadrature_order: usize,
    /// Deterministic fingerprint of the ordered physical quadrature points.
    pub quadrature_signature: u64,
}

/// Assemble variable-size extracted elements with sparse global storage.
/// `sample(cell,uvw)` must return physical geometry and basis data at the same
/// point. The closure is called once per quadrature point and no global dense
/// matrix is formed.
pub fn assemble_sampled<F, B>(
    ndof: usize,
    cell_count: usize,
    lambda: f64,
    mu: f64,
    order: usize,
    mut sample: F,
    mut body_force: B,
) -> Result<NativeAssembly, &'static str>
where
    F: FnMut(usize, Vector) -> Result<(ElementGeometrySample, ElementBasisSample), &'static str>,
    B: FnMut(Vector) -> Vector,
{
    if ndof == 0
        || cell_count == 0
        || !lambda.is_finite()
        || !mu.is_finite()
        || mu <= 0.
        || lambda + (2. / 3.) * mu <= 0.
    {
        return Err("invalid extracted elasticity problem");
    }
    let points = gauss(order)?;
    let mut rows = vec![BTreeMap::<usize, Matrix>::new(); ndof];
    let mut force = vec![0.; 3 * ndof];
    let mut volume = 0.;
    let mut point_count = 0usize;
    let mut qsignature = 0xcbf29ce484222325u64;
    for cell in 0..cell_count {
        let mut local_ids: Option<Vec<usize>> = None;
        let mut local = Vec::<Matrix>::new();
        for &(u, wu) in &points {
            for &(v, wv) in &points {
                for &(w, ww) in &points {
                    let (geometry, basis) = sample(cell, [u, v, w])?;
                    if !geometry.determinant.is_finite()
                        || geometry.determinant <= 0.
                        || geometry.position.iter().any(|x| !x.is_finite())
                        || basis.dofs.is_empty()
                        || basis.values.len() != basis.dofs.len()
                        || basis.physical_gradients.len() != basis.dofs.len()
                        || basis.dofs.iter().any(|&d| d >= ndof)
                        || basis.values.iter().any(|x| !x.is_finite())
                        || basis
                            .physical_gradients
                            .iter()
                            .flatten()
                            .any(|x| !x.is_finite())
                    {
                        return Err("invalid geometry or basis sample");
                    }
                    if let Some(ids) = &local_ids {
                        if ids != &basis.dofs {
                            return Err("element DOF map changed across quadrature points");
                        }
                    } else {
                        if basis.dofs.iter().copied().collect::<BTreeSet<_>>().len()
                            != basis.dofs.len()
                        {
                            return Err("duplicate global DOF in sampled element");
                        }
                        local_ids = Some(basis.dofs.clone());
                        let n = basis.dofs.len();
                        local = vec![[[0.; 3]; 3]; n * n];
                    }
                    let measure = geometry.determinant * wu * wv * ww;
                    for word in [
                        cell as u64,
                        u.to_bits(),
                        v.to_bits(),
                        w.to_bits(),
                        geometry.determinant.to_bits(),
                        geometry.position[0].to_bits(),
                        geometry.position[1].to_bits(),
                        geometry.position[2].to_bits(),
                    ]
                    .into_iter()
                    .chain(
                        geometry
                            .jacobian_columns
                            .iter()
                            .flatten()
                            .map(|x| x.to_bits()),
                    ) {
                        qsignature ^= word;
                        qsignature = qsignature.wrapping_mul(0x100000001b3);
                    }
                    volume += measure;
                    point_count += 1;
                    let force_value = body_force(geometry.position);
                    if force_value.iter().any(|x| !x.is_finite()) {
                        return Err("nonfinite body force");
                    }
                    for (a, &dof) in basis.dofs.iter().enumerate() {
                        for component in 0..3 {
                            force[3 * dof + component] +=
                                measure * basis.values[a] * force_value[component];
                        }
                    }
                    for (i, &gi) in basis.physical_gradients.iter().enumerate() {
                        for (j, &gj) in basis.physical_gradients.iter().enumerate().take(i + 1) {
                            let diag = mu * dot(gi, gj);
                            let entry = &mut local[i * basis.dofs.len() + j];
                            for a in 0..3 {
                                for b in 0..3 {
                                    let value = lambda * gi[a] * gj[b]
                                        + mu * gi[b] * gj[a]
                                        + if a == b { diag } else { 0. };
                                    entry[a][b] += measure * value;
                                }
                            }
                        }
                    }
                }
            }
        }
        let ids = local_ids.ok_or("empty extracted element")?;
        let n = ids.len();
        for i in 0..n {
            for j in 0..=i {
                let entry = local[i * n + j];
                let target = rows[ids[i]].entry(ids[j]).or_default();
                for a in 0..3 {
                    for b in 0..3 {
                        target[a][b] += entry[a][b];
                    }
                }
                if i != j {
                    let target = rows[ids[j]].entry(ids[i]).or_default();
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
    let mut columns = Vec::new();
    let mut values = Vec::new();
    for row in rows {
        for (j, value) in row {
            columns.push(j);
            values.push(value);
        }
        offsets.push(columns.len());
    }
    let stiffness = Stiffness {
        offsets,
        columns,
        values,
        integrated_volume: volume,
    };
    stiffness.validate()?;
    Ok(NativeAssembly {
        stiffness,
        force,
        quadrature_points: point_count,
        cell_count,
        quadrature_order: order,
        quadrature_signature: qsignature,
    })
}

/// Convenience assembler for a uniform list of extracted tensor elements.
pub fn assemble_elements(
    elements: &[ExtractedElement],
    ndof: usize,
    young: f64,
    poisson: f64,
    order: usize,
) -> Result<NativeAssembly, &'static str> {
    if !young.is_finite() || young <= 0. || !(-1. < poisson && poisson < 0.5) {
        return Err("invalid isotropic material");
    }
    for element in elements {
        element.validate()?;
    }
    let lambda = young * poisson / ((1. + poisson) * (1. - 2. * poisson));
    let mu = young / (2. * (1. + poisson));
    assemble_sampled(
        ndof,
        elements.len(),
        lambda,
        mu,
        order,
        |cell, uvw| elements[cell].sample_valid(uvw),
        |_| [0.; 3],
    )
}

/// A manufactured axial extension solution on an extrusion-like solid:
/// `u=(0,0,sin(pi*(z-z0)/L))`. It is zero on both axial ends, is driven by the
/// analytic body force and receives exact Cauchy traction on all design faces.
/// `(lambda+2*mu)*(pi/L)^2*sin(pi*(z-z0)/L)` in z.
#[derive(Debug, Clone)]
pub struct ManufacturedResult {
    pub scalar_dofs: usize,
    pub free_scalar_dofs: usize,
    pub quadrature_points: usize,
    pub matrix_nonzeros: usize,
    pub volume_mm3: f64,
    pub l2_displacement_error_mm_sqrt_mm3: f64,
    pub energy_error_sqrt_n_mm: f64,
    pub relative_residual: f64,
    pub iterations: usize,
    pub strain_energy_n_mm: f64,
    pub applied_torque_n_mm: f64,
    pub reaction_torque_n_mm: f64,
    pub torque_balance_relative: f64,
    pub assembly_ms: f64,
    pub solve_ms: f64,
    pub error_evaluation_ms: f64,
    pub torque_relative_error: f64,
    pub energy_relative_error: f64,
    pub displacement_relative_error: f64,
}

#[derive(Debug, Clone)]
pub struct GalerkinComparison {
    pub same_quadrature_points: bool,
    pub stiffness_relative_frobenius_defect: f64,
    pub force_relative_l2_defect: f64,
    pub c0_scalar_dofs: usize,
    pub constrained_scalar_dofs: usize,
    pub matrix_nonzeros_c0: usize,
    pub matrix_nonzeros_constrained: usize,
}

/// Verify the exact subspace assembly identities on a common physical rule.
/// The pullback is applied in sparse coefficient form; no global dense matrix
/// or dense constraint operator is constructed.
pub fn compare_galerkin<F>(
    c0: &crate::inherited_space::InheritedSpace,
    constrained: &crate::inherited_space::InheritedSpace,
    young_mpa: f64,
    poisson: f64,
    order: usize,
    body_force: F,
) -> Result<GalerkinComparison, &'static str>
where
    F: Fn(Vector) -> Vector,
{
    if !young_mpa.is_finite()
        || young_mpa <= 0.
        || !(-1. < poisson && poisson < 0.5)
        || c0.cells.len() != constrained.cells.len()
        || c0.ndof == 0
    {
        return Err("invalid Galerkin comparison spaces or material");
    }
    let lambda = young_mpa * poisson / ((1. + poisson) * (1. - 2. * poisson));
    let mu = young_mpa / (2. * (1. + poisson));
    let assemble_space = |space: &crate::inherited_space::InheritedSpace| {
        assemble_sampled(
            space.ndof,
            space.cells.len(),
            lambda,
            mu,
            order,
            |cell, uvw| space.cells[cell].sample(uvw),
            &body_force,
        )
    };
    let k0 = assemble_space(c0)?;
    let k1 = assemble_space(constrained)?;
    compare_assembled(&k0, &k1, c0.ndof, constrained.ndof, &constrained.lift)
}

/// Check Galerkin stiffness/load identities for two assemblies on exactly the
/// same ordered physical quadrature rule. `lift` maps each large-space scalar
/// coefficient to the smaller-space coefficients using exact rational rows.
pub fn compare_assembled(
    large: &NativeAssembly,
    small: &NativeAssembly,
    large_ndof: usize,
    small_ndof: usize,
    lift: &[BTreeMap<usize, crate::constraint_system::Ratio>],
) -> Result<GalerkinComparison, &'static str> {
    large.stiffness.validate()?;
    small.stiffness.validate()?;
    if large_ndof == 0
        || small_ndof == 0
        || large.stiffness.offsets.len() != large_ndof + 1
        || small.stiffness.offsets.len() != small_ndof + 1
        || large.force.len() != 3 * large_ndof
        || small.force.len() != 3 * small_ndof
        || lift.len() != large_ndof
    {
        return Err("Galerkin assembly dimensions do not match coefficient spaces");
    }
    if large.quadrature_points != small.quadrature_points
        || large.cell_count != small.cell_count
        || large.quadrature_order != small.quadrature_order
        || large.quadrature_signature != small.quadrature_signature
    {
        return Err("Galerkin assemblies do not use the same physical quadrature points");
    }
    let prolong: Vec<Vec<(usize, f64)>> = lift
        .iter()
        .map(|row| {
            row.iter()
                .map(|(&id, &ratio)| (id, ratio.numerator() as f64 / ratio.denominator() as f64))
                .collect()
        })
        .collect();
    if prolong
        .iter()
        .flatten()
        .any(|(id, value)| *id >= small_ndof || !value.is_finite())
    {
        return Err("invalid exact Galerkin transfer");
    }
    let mut pulled = BTreeMap::<(usize, usize), Matrix>::new();
    for (i, terms_i) in prolong.iter().enumerate() {
        for entry in large.stiffness.offsets[i]..large.stiffness.offsets[i + 1] {
            let j = large.stiffness.columns[entry];
            let block = large.stiffness.values[entry];
            for &(a, ca) in terms_i {
                for &(b, cb) in &prolong[j] {
                    let target = pulled.entry((a, b)).or_default();
                    for r in 0..3 {
                        for s in 0..3 {
                            target[r][s] += ca * cb * block[r][s];
                        }
                    }
                }
            }
        }
    }
    let mut direct = BTreeMap::<(usize, usize), Matrix>::new();
    for i in 0..small_ndof {
        for e in small.stiffness.offsets[i]..small.stiffness.offsets[i + 1] {
            direct.insert((i, small.stiffness.columns[e]), small.stiffness.values[e]);
        }
    }
    let keys: BTreeSet<_> = pulled.keys().chain(direct.keys()).copied().collect();
    let (mut diff_sq, mut ref_sq) = (0., 0.);
    for key in keys {
        let a = pulled.get(&key).copied().unwrap_or_default();
        let b = direct.get(&key).copied().unwrap_or_default();
        for r in 0..3 {
            for s in 0..3 {
                diff_sq += (a[r][s] - b[r][s]).powi(2);
                ref_sq += b[r][s].powi(2);
            }
        }
    }
    let mut pulled_force = vec![0.; 3 * small_ndof];
    for (old, terms) in prolong.iter().enumerate() {
        for &(new, coef) in terms {
            for c in 0..3 {
                pulled_force[3 * new + c] += coef * large.force[3 * old + c];
            }
        }
    }
    let force_diff = pulled_force
        .iter()
        .zip(&small.force)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt();
    let force_ref = small.force.iter().map(|v| v * v).sum::<f64>().sqrt();
    Ok(GalerkinComparison {
        same_quadrature_points: true,
        stiffness_relative_frobenius_defect: diff_sq.sqrt() / ref_sq.sqrt().max(f64::MIN_POSITIVE),
        force_relative_l2_defect: force_diff / force_ref.max(f64::MIN_POSITIVE),
        c0_scalar_dofs: large_ndof,
        constrained_scalar_dofs: small_ndof,
        matrix_nonzeros_c0: large.stiffness.columns.len(),
        matrix_nonzeros_constrained: small.stiffness.columns.len(),
    })
}

/// Assemble, constrain, solve and measure the smooth manufactured axial
/// solution. The callback must return the same geometry/field point for both
/// assembly and error evaluation; only per-element matrices and global vectors
/// are stored, never a dense global operator.
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
pub fn solve_axial_manufactured<F>(
    ndof: usize,
    control_xyz: &[Vector],
    z_bounds: [f64; 2],
    cell_count: usize,
    young_mpa: f64,
    poisson: f64,
    order: usize,
    tolerance: f64,
    max_iterations: usize,
    boundary_faces: &[(usize, usize, bool)],
    mut sample: F,
) -> Result<ManufacturedResult, &'static str>
where
    F: FnMut(usize, Vector) -> Result<(ElementGeometrySample, ElementBasisSample), &'static str>,
{
    if ndof == 0
        || control_xyz.len() != ndof
        || control_xyz.iter().flatten().any(|x| !x.is_finite())
        || !z_bounds[0].is_finite()
        || !z_bounds[1].is_finite()
        || z_bounds[1] <= z_bounds[0]
        || !young_mpa.is_finite()
        || young_mpa <= 0.
        || !(-1. < poisson && poisson < 0.5)
    {
        return Err("invalid manufactured elasticity domain or material");
    }
    let lambda = young_mpa * poisson / ((1. + poisson) * (1. - 2. * poisson));
    let mu = young_mpa / (2. * (1. + poisson));
    let length = z_bounds[1] - z_bounds[0];
    let wave = std::f64::consts::PI / length;
    let assembly_start = std::time::Instant::now();
    let mut assembled = assemble_sampled(ndof, cell_count, lambda, mu, order, &mut sample, |x| {
        let value = (wave * (x[2] - z_bounds[0])).sin();
        [0., 0., (lambda + 2. * mu) * wave * wave * value]
    })?;
    let surface_q = gauss(order)?;
    for &(cell, axis, upper) in boundary_faces {
        if cell >= cell_count || axis > 2 {
            return Err("invalid manufactured boundary face reference");
        }
        let tangents = crate::analysis_space::face_axes(axis);
        for &(q0, w0) in &surface_q {
            for &(q1, w1) in &surface_q {
                let mut uvw = [0.; 3];
                uvw[tangents[0]] = q0;
                uvw[tangents[1]] = q1;
                uvw[axis] = if upper { 1. } else { 0. };
                let (geometry, basis) = sample(cell, uvw)?;
                let a = geometry.jacobian_columns[tangents[0]];
                let b = geometry.jacobian_columns[tangents[1]];
                let raw = cross(a, b);
                let area = norm(raw);
                if !area.is_finite() || area <= 0. {
                    return Err("invalid manufactured boundary surface Jacobian");
                }
                let orientation = if dot(raw, geometry.jacobian_columns[axis]) >= 0. {
                    1.
                } else {
                    -1.
                };
                let outward = if upper { orientation } else { -orientation };
                let normal = raw.map(|x| outward * x / area);
                let z = geometry.position[2];
                let wave = std::f64::consts::PI / length;
                let derivative = wave * (wave * (z - z_bounds[0])).cos();
                let stress = [
                    lambda * derivative,
                    lambda * derivative,
                    (lambda + 2. * mu) * derivative,
                ];
                let traction: Vector = std::array::from_fn(|c| stress[c] * normal[c]);
                let measure = area * w0 * w1;
                for i in 0..basis.dofs.len() {
                    for c in 0..3 {
                        assembled.force[3 * basis.dofs[i] + c] +=
                            measure * basis.values[i] * traction[c];
                    }
                }
            }
        }
    }
    let assembly_ms = assembly_start.elapsed().as_secs_f64() * 1000.;
    let boundary_tolerance = 1e-10 * length.max(1.);
    let prescribed: Vec<_> = control_xyz
        .iter()
        .enumerate()
        .filter_map(|(dof, p)| {
            ((p[2] - z_bounds[0]).abs() <= boundary_tolerance
                || (p[2] - z_bounds[1]).abs() <= boundary_tolerance)
                .then_some((dof, crate::constraint_system::Ratio::zero()))
        })
        .collect();
    if prescribed.is_empty() {
        return Err("manufactured axial solve found no end-face coefficients");
    }
    let lift = crate::constraint_system::resolve(ndof, &[], &prescribed)?;
    let solve_start = std::time::Instant::now();
    let solution = crate::constrained_elasticity::solve_constrained(
        &assembled.stiffness,
        [&lift, &lift, &lift],
        &assembled.force,
        tolerance,
        max_iterations,
    )?;
    let solve_ms = solve_start.elapsed().as_secs_f64() * 1000.;

    let error_start = std::time::Instant::now();
    let q = gauss(order)?;
    let mut l2 = 0.;
    let mut energy = 0.;
    for cell in 0..cell_count {
        for &(u, wu) in &q {
            for &(v, wv) in &q {
                for &(w, ww) in &q {
                    let (geometry, basis) = sample(cell, [u, v, w])?;
                    if basis.dofs.len() != basis.values.len()
                        || basis.dofs.len() != basis.physical_gradients.len()
                    {
                        return Err("invalid extracted error-evaluation sample");
                    }
                    let mut uh = [0.; 3];
                    let mut grad = [[0.; 3]; 3]; // component, physical derivative
                    for a in 0..basis.dofs.len() {
                        let dof = basis.dofs[a];
                        if dof >= ndof {
                            return Err("error sample DOF out of range");
                        }
                        for c in 0..3 {
                            let coefficient = solution.displacement_mm[3 * dof + c];
                            uh[c] += basis.values[a] * coefficient;
                            for d in 0..3 {
                                grad[c][d] += basis.physical_gradients[a][d] * coefficient;
                            }
                        }
                    }
                    let z = geometry.position[2];
                    let phase = wave * (z - z_bounds[0]);
                    let exact = [0., 0., phase.sin()];
                    let exact_grad = [0., 0., wave * phase.cos()];
                    let mut div_error = 0.;
                    let mut strain_sq = 0.;
                    for c in 0..3 {
                        let error = grad[c][c] - if c == 2 { exact_grad[2] } else { 0. };
                        div_error += error;
                        for d in 0..3 {
                            let e_cd = 0.5 * (grad[c][d] + grad[d][c]);
                            let exact_e = if c == 2 && d == 2 { exact_grad[2] } else { 0. };
                            let delta = e_cd - exact_e;
                            strain_sq += delta * delta;
                        }
                    }
                    let displacement_sq = (0..3).map(|c| (uh[c] - exact[c]).powi(2)).sum::<f64>();
                    let measure = geometry.determinant * wu * wv * ww;
                    l2 += measure * displacement_sq;
                    energy += measure * (lambda * div_error * div_error + 2. * mu * strain_sq);
                }
            }
        }
    }
    if energy < -1e-10 || l2 < -1e-10 || !energy.is_finite() || !l2.is_finite() {
        return Err("invalid manufactured error integral");
    }
    let error_evaluation_ms = error_start.elapsed().as_secs_f64() * 1000.;
    Ok(ManufacturedResult {
        scalar_dofs: ndof,
        free_scalar_dofs: lift.free_roots.len(),
        quadrature_points: assembled.quadrature_points,
        matrix_nonzeros: assembled.stiffness.columns.len(),
        volume_mm3: assembled.stiffness.integrated_volume,
        l2_displacement_error_mm_sqrt_mm3: l2.max(0.).sqrt(),
        energy_error_sqrt_n_mm: energy.max(0.).sqrt(),
        relative_residual: solution.relative_residual,
        iterations: solution.iterations,
        strain_energy_n_mm: solution.strain_energy_n_mm,
        applied_torque_n_mm: 0.,
        reaction_torque_n_mm: 0.,
        torque_balance_relative: 0.,
        assembly_ms,
        solve_ms,
        error_evaluation_ms,
        torque_relative_error: 0.,
        energy_relative_error: 0.,
        displacement_relative_error: 0.,
    })
}

/// Native Saint-Venant torsion solve on an audited extrusion space. The top
/// traction is `G*twist*(-dy,dx,0)` and the bottom trace is clamped. This
/// is admitted only for its exact circular-cylinder source family.
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
fn solve_extrusion_torsion(
    space: &crate::inherited_space::InheritedSpace,
    young_mpa: f64,
    poisson: f64,
    center_xy: [f64; 2],
    radius_mm: f64,
    twist_per_mm: f64,
    order: usize,
    tolerance: f64,
    max_iterations: usize,
) -> Result<ManufacturedResult, &'static str> {
    if !young_mpa.is_finite()
        || young_mpa <= 0.
        || !(-1. < poisson && poisson < 0.5)
        || center_xy.iter().any(|x| !x.is_finite())
        || !radius_mm.is_finite()
        || radius_mm <= 0.
        || !twist_per_mm.is_finite()
        || twist_per_mm == 0.
    {
        return Err("invalid extrusion torsion parameters");
    }
    let lambda = young_mpa * poisson / ((1. + poisson) * (1. - 2. * poisson));
    let mu = young_mpa / (2. * (1. + poisson));
    let assembly_start = std::time::Instant::now();
    let mut assembled = assemble_sampled(
        space.ndof,
        space.cells.len(),
        lambda,
        mu,
        order,
        |cell, uvw| space.cells[cell].sample(uvw),
        |_| [0.; 3],
    )?;
    let zmin = space
        .cells
        .iter()
        .map(|cell| cell.z_bounds[0])
        .fold(f64::INFINITY, f64::min);
    let zmax = space
        .cells
        .iter()
        .map(|cell| cell.z_bounds[1])
        .fold(f64::NEG_INFINITY, f64::max);
    if !zmin.is_finite() || !zmax.is_finite() || zmax <= zmin {
        return Err("invalid extrusion axial bounds");
    }
    let points = gauss(order)?;
    let mut applied_torque = 0.;
    let mut resultant = [0.; 3];
    let face_tol = 1e-10 * (zmax - zmin).max(1.);
    for (ci, cell) in space.cells.iter().enumerate() {
        if (cell.z_bounds[1] - zmax).abs() > face_tol {
            continue;
        }
        for &(u, wu) in &points {
            for &(v, wv) in &points {
                let (geometry, basis) = cell.sample([u, v, 1.])?;
                let a = geometry.jacobian_columns[0];
                let b = geometry.jacobian_columns[1];
                let area_jac = norm(cross(a, b));
                if !area_jac.is_finite() || area_jac <= 0. {
                    return Err("invalid torsion end-face Jacobian");
                }
                let dx = geometry.position[0] - center_xy[0];
                let dy = geometry.position[1] - center_xy[1];
                let traction = [-mu * twist_per_mm * dy, mu * twist_per_mm * dx, 0.];
                let measure = area_jac * wu * wv;
                resultant = std::array::from_fn(|c| resultant[c] + measure * traction[c]);
                applied_torque += measure * (dx * traction[1] - dy * traction[0]);
                for i in 0..basis.dofs.len() {
                    for c in 0..3 {
                        assembled.force[3 * basis.dofs[i] + c] +=
                            measure * basis.values[i] * traction[c];
                    }
                }
            }
        }
        let _ = ci;
    }
    let assembly_ms = assembly_start.elapsed().as_secs_f64() * 1000.;
    let prescribed: Vec<_> = space
        .control_xyz
        .iter()
        .enumerate()
        .filter_map(|(dof, p)| {
            ((p[2] - zmin).abs() <= face_tol)
                .then_some((dof, crate::constraint_system::Ratio::zero()))
        })
        .collect();
    if prescribed.is_empty() {
        return Err("torsion solve found no clamped end-face coefficients");
    }
    let lift = crate::constraint_system::resolve(space.ndof, &[], &prescribed)?;
    let solve_start = std::time::Instant::now();
    let solution = crate::constrained_elasticity::solve_constrained(
        &assembled.stiffness,
        [&lift, &lift, &lift],
        &assembled.force,
        tolerance,
        max_iterations,
    )?;
    let solve_ms = solve_start.elapsed().as_secs_f64() * 1000.;

    let error_start = std::time::Instant::now();
    let mut l2 = 0.;
    let mut energy_error = 0.;
    for (ci, cell) in space.cells.iter().enumerate() {
        for &(u, wu) in &points {
            for &(v, wv) in &points {
                for &(w, ww) in &points {
                    let (geometry, basis) = cell.sample([u, v, w])?;
                    let mut uh = [0.; 3];
                    let mut grad = [[0.; 3]; 3];
                    for a in 0..basis.dofs.len() {
                        let dof = basis.dofs[a];
                        if dof >= space.ndof {
                            return Err("torsion sample DOF out of range");
                        }
                        for c in 0..3 {
                            let coef = solution.displacement_mm[3 * dof + c];
                            uh[c] += basis.values[a] * coef;
                            for d in 0..3 {
                                grad[c][d] += basis.physical_gradients[a][d] * coef;
                            }
                        }
                    }
                    let dx = geometry.position[0] - center_xy[0];
                    let dy = geometry.position[1] - center_xy[1];
                    let phi = twist_per_mm * (geometry.position[2] - zmin);
                    let exact = [-phi * dy, phi * dx, 0.];
                    let exact_grad = [
                        [-0., 0., -twist_per_mm * dy],
                        [0., 0., twist_per_mm * dx],
                        [0.; 3],
                    ];
                    l2 += geometry.determinant
                        * wu
                        * wv
                        * ww
                        * (0..3).map(|c| (uh[c] - exact[c]).powi(2)).sum::<f64>();
                    let mut div = 0.;
                    let mut strain_sq = 0.;
                    for c in 0..3 {
                        div += grad[c][c] - exact_grad[c][c];
                        for d in 0..3 {
                            let e = 0.5 * (grad[c][d] + grad[d][c]);
                            let ee = 0.5 * (exact_grad[c][d] + exact_grad[d][c]);
                            strain_sq += (e - ee).powi(2);
                        }
                    }
                    energy_error += geometry.determinant
                        * wu
                        * wv
                        * ww
                        * (lambda * div * div + 2. * mu * strain_sq);
                }
            }
        }
        let _ = ci;
    }
    let mut reaction = [0.; 3];
    let mut reaction_torque = 0.;
    let internal = assembled.stiffness.apply(&solution.displacement_mm)?;
    for (dof, p) in space.control_xyz.iter().enumerate() {
        if (p[2] - zmin).abs() <= face_tol {
            let f = [
                internal[3 * dof] - assembled.force[3 * dof],
                internal[3 * dof + 1] - assembled.force[3 * dof + 1],
                internal[3 * dof + 2] - assembled.force[3 * dof + 2],
            ];
            for c in 0..3 {
                reaction[c] += f[c];
            }
            reaction_torque += (p[0] - center_xy[0]) * f[1] - (p[1] - center_xy[1]) * f[0];
        }
    }
    let torque_balance_relative =
        (reaction_torque + applied_torque).abs() / applied_torque.abs().max(f64::MIN_POSITIVE);
    let exact_torque = mu * twist_per_mm * std::f64::consts::PI * radius_mm.powi(4) / 2.;
    let height = zmax - zmin;
    let exact_energy = exact_torque * twist_per_mm * height / 2.;
    let exact_l2 =
        (std::f64::consts::PI * radius_mm.powi(4) / 2. * twist_per_mm.powi(2) * height.powi(3)
            / 3.)
            .sqrt();
    let error_evaluation_ms = error_start.elapsed().as_secs_f64() * 1000.;
    let _ = resultant;
    Ok(ManufacturedResult {
        scalar_dofs: space.ndof,
        free_scalar_dofs: lift.free_roots.len(),
        quadrature_points: assembled.quadrature_points,
        matrix_nonzeros: assembled.stiffness.columns.len(),
        volume_mm3: assembled.stiffness.integrated_volume,
        l2_displacement_error_mm_sqrt_mm3: l2.max(0.).sqrt(),
        energy_error_sqrt_n_mm: energy_error.max(0.).sqrt(),
        relative_residual: solution.relative_residual,
        iterations: solution.iterations,
        strain_energy_n_mm: solution.strain_energy_n_mm,
        applied_torque_n_mm: applied_torque,
        reaction_torque_n_mm: reaction_torque,
        torque_balance_relative,
        assembly_ms,
        solve_ms,
        error_evaluation_ms,
        torque_relative_error: (applied_torque / exact_torque - 1.).abs(),
        energy_relative_error: (solution.strain_energy_n_mm / exact_energy - 1.).abs(),
        displacement_relative_error: l2.max(0.).sqrt() / exact_l2,
    })
}

fn norm(v: Vector) -> f64 {
    dot(v, v).sqrt()
}

fn model_boundary_faces(
    model: &crate::json::Json,
    space: &crate::inherited_space::InheritedSpace,
) -> Result<Vec<(usize, usize, bool)>, &'static str> {
    let mut cell_by_id = BTreeMap::new();
    for (i, cell) in space.cells.iter().enumerate() {
        cell_by_id.insert(cell.id.as_str(), i);
    }
    let mut faces = BTreeSet::new();
    for face in model.field("design_faces")?.items()? {
        let id = face.field("block")?.text()?;
        let side = face.field("side")?.text()?;
        let (axis, upper) = match side {
            "u0" => (0, false),
            "u1" => (0, true),
            "v0" => (1, false),
            "v1" => (1, true),
            "w0" => (2, false),
            "w1" => (2, true),
            _ => return Err("unknown design face side"),
        };
        let cell = *cell_by_id
            .get(id)
            .ok_or("design face does not resolve to inherited cell")?;
        faces.insert((cell, axis, upper));
    }
    Ok(faces.into_iter().collect())
}

/// Solve an inherited spline field using a compact versioned JSON job.
/// Supported cases are `cylinder_torsion` and `axial_manufactured`.
pub fn solve_inherited(
    model: &crate::json::Json,
    job: &crate::json::Json,
) -> Result<crate::json::Json, &'static str> {
    use crate::json::{object, Json};
    fn number<T: std::str::FromStr>(job: &Json, key: &str, default: T) -> Result<T, &'static str> {
        match job.field(key) {
            Ok(value) => value.number(),
            Err(_) => Ok(default),
        }
    }
    fn bool_value(job: &Json, key: &str, default: bool) -> Result<bool, &'static str> {
        match job.field(key) {
            Ok(Json::Bool(v)) => Ok(*v),
            Ok(_) => Err("expected a JSON boolean"),
            Err(_) => Ok(default),
        }
    }
    fn as_number(value: impl ToString) -> Json {
        Json::Number(value.to_string())
    }
    let fields = match job {
        Json::Object(fields) => fields,
        _ => return Err("inherited elasticity job must be an object"),
    };
    const ALLOWED: &[&str] = &[
        "schema",
        "case",
        "degree",
        "continuity",
        "restore_strip",
        "transverse_interfaces",
        "young_mpa",
        "poisson",
        "quadrature_order",
        "tolerance",
        "max_iterations",
        "center_xy_mm",
        "twist_per_mm",
    ];
    if fields.keys().any(|key| !ALLOWED.contains(&key.as_str())) {
        return Err("unknown inherited elasticity job field");
    }
    for key in [
        "schema",
        "case",
        "degree",
        "continuity",
        "young_mpa",
        "poisson",
        "quadrature_order",
        "tolerance",
        "max_iterations",
    ] {
        if !fields.contains_key(key) {
            return Err("missing required inherited elasticity job field");
        }
    }
    if job.field("schema")?.text()? != "inherited_elasticity_job_v1" {
        return Err("unsupported inherited elasticity job schema");
    }
    let case = job.field("case")?.text()?;
    if case == "cylinder_torsion"
        && (!fields.contains_key("center_xy_mm") || !fields.contains_key("twist_per_mm"))
    {
        return Err("cylinder torsion job requires center_xy_mm and twist_per_mm");
    }
    if case != "cylinder_torsion" && case != "axial_manufactured" {
        return Err("unsupported inherited elasticity case");
    }
    let degree: usize = number(job, "degree", 0)?;
    let continuity: usize = number(job, "continuity", 0)?;
    let restore_strip = bool_value(job, "restore_strip", false)?;
    let mut transverse = Vec::<(String, String)>::new();
    if let Ok(records) = job.field("transverse_interfaces") {
        for record in records.items()? {
            let pair = record.items()?;
            if pair.len() != 2 {
                return Err("transverse interface needs a block pair");
            }
            transverse.push((pair[0].text()?.to_owned(), pair[1].text()?.to_owned()));
        }
    }
    let space =
        crate::inherited_space::build_space(model, degree, continuity, restore_strip, &transverse)?;
    let young: f64 = number(job, "young_mpa", 0.)?;
    let poisson: f64 = number(job, "poisson", 0.)?;
    let order: usize = number(job, "quadrature_order", 0)?;
    let tolerance: f64 = number(job, "tolerance", 0.)?;
    let max_iterations: usize = number(job, "max_iterations", 0)?;
    let result = match case {
        "cylinder_torsion" => {
            if model.field("schema")?.text()? != "box_cylinder_boolean_suite_v1"
                || model.field("operation")?.text()? != "cylinder"
                || model.field("tool")?.field("kind")?.text()? != "cylinder"
                || model.field("tool")?.field("axis")?.text()? != "z"
            {
                return Err(
                    "cylinder torsion case requires the audited z-axis cylinder source model",
                );
            }
            let tool = model.field("tool")?;
            let source_center = tool.field("center_mm")?.items()?;
            if source_center.len() != 2 {
                return Err("invalid cylinder center metadata");
            }
            let source_center = [
                source_center[0].number::<f64>()?,
                source_center[1].number::<f64>()?,
            ];
            let radius = tool.field("radius_mm")?.number::<f64>()?;
            let center = match job.field("center_xy_mm") {
                Ok(v) => {
                    let a = v.items()?;
                    if a.len() != 2 {
                        return Err("center_xy_mm needs two coordinates");
                    }
                    [a[0].number()?, a[1].number()?]
                }
                Err(_) => [0., 0.],
            };
            if center != source_center {
                return Err("torsion center does not match the source cylinder");
            }
            let twist: f64 = number(job, "twist_per_mm", 1e-3)?;
            solve_extrusion_torsion(
                &space,
                young,
                poisson,
                center,
                radius,
                twist,
                order,
                tolerance,
                max_iterations,
            )?
        }
        "axial_manufactured" => solve_axial_manufactured(
            space.ndof,
            &space.control_xyz,
            [
                space
                    .cells
                    .iter()
                    .map(|c| c.z_bounds[0])
                    .fold(f64::INFINITY, f64::min),
                space
                    .cells
                    .iter()
                    .map(|c| c.z_bounds[1])
                    .fold(f64::NEG_INFINITY, f64::max),
            ],
            space.cells.len(),
            young,
            poisson,
            order,
            tolerance,
            max_iterations,
            &model_boundary_faces(model, &space)?,
            |cell, uvw| space.cells[cell].sample(uvw),
        )?,
        _ => return Err("unsupported inherited elasticity case"),
    };
    Ok(object([
        (
            "schema",
            Json::String("inherited_elasticity_result_v1".into()),
        ),
        ("case", Json::String(case.into())),
        ("degree", as_number(degree)),
        ("continuity", as_number(continuity)),
        ("scalar_dofs", as_number(result.scalar_dofs)),
        ("vector_dofs", as_number(3 * result.scalar_dofs)),
        ("free_scalar_dofs", as_number(result.free_scalar_dofs)),
        ("quadrature_points", as_number(result.quadrature_points)),
        ("matrix_nonzeros", as_number(result.matrix_nonzeros)),
        ("volume_mm3", as_number(result.volume_mm3)),
        (
            "l2_displacement_error",
            as_number(result.l2_displacement_error_mm_sqrt_mm3),
        ),
        ("energy_error", as_number(result.energy_error_sqrt_n_mm)),
        ("relative_residual", as_number(result.relative_residual)),
        ("iterations", as_number(result.iterations)),
        ("strain_energy_n_mm", as_number(result.strain_energy_n_mm)),
        ("applied_torque_n_mm", as_number(result.applied_torque_n_mm)),
        (
            "reaction_torque_n_mm",
            as_number(result.reaction_torque_n_mm),
        ),
        (
            "torque_balance_relative",
            as_number(result.torque_balance_relative),
        ),
        ("assembly_ms", as_number(result.assembly_ms)),
        ("solve_ms", as_number(result.solve_ms)),
        ("error_evaluation_ms", as_number(result.error_evaluation_ms)),
        (
            "torque_relative_error",
            as_number(result.torque_relative_error),
        ),
        (
            "energy_relative_error",
            as_number(result.energy_relative_error),
        ),
        (
            "displacement_relative_error",
            as_number(result.displacement_relative_error),
        ),
        ("c0_scalar_dofs", as_number(space.c0_ndof)),
        ("inherited_join_count", as_number(space.common_parent_faces)),
        (
            "transverse_c1_join_count",
            as_number(space.transverse_faces),
        ),
        ("retained_c0_join_count", as_number(space.retained_faces)),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn identity_axis(degree: usize) -> AxisExtraction {
        let n = degree + 1;
        AxisExtraction {
            degree,
            basis_count: n,
            row_offsets: (0..=n).collect(),
            columns: (0..n).collect(),
            coefficients: vec![Ratio::one(); n],
        }
    }

    fn affine_cube(degree: usize) -> ExtractedElement {
        let axes = [
            identity_axis(degree),
            identity_axis(degree),
            identity_axis(degree),
        ];
        let n = degree + 1;
        let mut dofs = Vec::new();
        let mut controls = Vec::new();
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    dofs.push((i * n + j) * n + k);
                    let x = i as f64 / degree as f64;
                    let y = j as f64 / degree as f64;
                    let z = k as f64 / degree as f64;
                    controls.push([x, y, z, 1.]);
                }
            }
        }
        ExtractedElement {
            axes,
            field_weights: vec![1.; dofs.len()],
            dofs,
            geometry_controls: controls,
        }
    }

    #[test]
    fn cubic_sparse_axis_extraction_preserves_partition_and_derivative_sum() {
        let axis = AxisExtraction {
            degree: 3,
            basis_count: 4,
            row_offsets: vec![0, 1, 3, 5, 6],
            columns: vec![0, 0, 1, 1, 2, 3],
            coefficients: vec![
                Ratio::one(),
                Ratio::new(1, 2).unwrap(),
                Ratio::new(1, 2).unwrap(),
                Ratio::new(1, 2).unwrap(),
                Ratio::new(1, 2).unwrap(),
                Ratio::one(),
            ],
        };
        axis.validate().unwrap();
        for t in [0., 0.13, 0.5, 0.91, 1.] {
            let (bern, derivative) = bernstein(3, t).unwrap();
            let (values, gradients) = axis.apply(&bern, &derivative);
            assert!((values.iter().sum::<f64>() - 1.).abs() < 2e-15);
            assert!(gradients.iter().sum::<f64>().abs() < 2e-14);
        }
    }

    #[test]
    fn affine_cubic_element_reproduces_geometry_and_physical_gradients() {
        let element = affine_cube(3);
        let (geometry, basis) = element.sample([0.27, 0.41, 0.73]).unwrap();
        for (a, expected) in geometry.position.into_iter().zip([0.27, 0.41, 0.73]) {
            assert!((a - expected).abs() < 2e-14);
        }
        assert!((geometry.determinant - 1.).abs() < 2e-13);
        assert!((basis.values.iter().sum::<f64>() - 1.).abs() < 2e-14);
        for c in 0..3 {
            let sum = basis.physical_gradients.iter().map(|g| g[c]).sum::<f64>();
            assert!(sum.abs() < 2e-13);
        }
    }

    #[test]
    fn variable_size_assembly_is_sparse_and_integrates_cube_volume() {
        let q2 = affine_cube(2);
        let q3 = affine_cube(3);
        let a = assemble_elements(&[q2], 27, 1., 0.25, 3).unwrap();
        let b = assemble_elements(&[q3], 64, 1., 0.25, 4).unwrap();
        assert!((a.stiffness.integrated_volume - 1.).abs() < 1e-12);
        assert!((b.stiffness.integrated_volume - 1.).abs() < 1e-12);
        assert_eq!(a.force.len(), 81);
        assert_eq!(b.force.len(), 192);
        assert!(a.stiffness.validate().is_ok());
        assert!(b.stiffness.validate().is_ok());
    }
}
