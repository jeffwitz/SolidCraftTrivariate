//! Additive geometric preconditioning for nested, homogeneously mounted Q2 spaces.
//!
//! The coarse matrix is assembled as CᵀKC. Each finer level adds its positive
//! inverse reduced diagonal to the prolongation of the coarser solver.

use crate::{
    adaptive_partition::AdaptivePartition,
    adaptive_space::AdaptiveSpace,
    adaptive_transfer::build_transfer,
    constrained_elasticity::ReducedPreconditioner,
    constraint_system::{Lift, Ratio},
    elasticity::Stiffness,
};

/// Geometry, dyadic mesh and common homogeneous scalar mounting for xyz.
#[derive(Clone, Copy)]
pub struct MountedSpace<'a> {
    pub adaptive: &'a AdaptiveSpace,
    pub partition: &'a AdaptivePartition,
    pub lift: &'a Lift,
}

struct Level {
    rows: Vec<Vec<(usize, f64)>>,
    source_count: usize,
    inverse_diagonal: Vec<f64>,
}

/// A coarse Cholesky factor and exact nested-space transfer stencils.
pub struct MultilevelPreconditioner {
    coarse: Vec<f64>,
    first_column: Vec<usize>,
    levels: Vec<Level>,
    last_mount: Lift,
    last_mesh: Option<AdaptiveSpace>,
}

fn coefficient(r: Ratio) -> f64 {
    r.numerator() as f64 / r.denominator() as f64
}

fn validate(lift: &Lift, stiffness: &Stiffness) -> Result<(), &'static str> {
    stiffness.validate()?;
    for i in 0..stiffness.offsets.len() - 1 {
        let range = stiffness.offsets[i]..stiffness.offsets[i + 1];
        if stiffness.columns[range].windows(2).any(|w| w[0] >= w[1]) {
            return Err("multilevel CSR columns must be strictly ordered");
        }
    }
    for i in 0..stiffness.offsets.len() - 1 {
        for p in stiffness.offsets[i]..stiffness.offsets[i + 1] {
            let j = stiffness.columns[p];

            let row = stiffness.offsets[j]..stiffness.offsets[j + 1];
            let q = stiffness.columns[row.clone()]
                .binary_search(&i)
                .map_err(|_| "asymmetric multilevel sparsity")?
                + row.start;
            if j < i {
                continue;
            }
            let a = stiffness.values[p];
            let b = stiffness.values[q];
            let scale = a
                .iter()
                .chain(&b)
                .flatten()
                .fold(1_f64, |v, x| v.max(x.abs()));
            for (c, r) in a.iter().enumerate() {
                for (d, &x) in r.iter().enumerate() {
                    if (x - b[d][c]).abs() > 2e-12 * scale {
                        return Err("asymmetric multilevel stiffness");
                    }
                }
            }
        }
    }
    let n = lift.free_roots.len();
    if n == 0
        || stiffness.offsets.len() != lift.lines.len() + 1
        || lift.free_roots.windows(2).any(|p| p[0] >= p[1])
    {
        return Err("invalid multilevel mounting dimensions");
    }
    for (i, line) in lift.lines.iter().enumerate() {
        if line.dof != i
            || line.offset != Ratio::zero()
            || line.terms.windows(2).any(|p| p[0].0 >= p[1].0)
            || line
                .terms
                .iter()
                .any(|&(j, w)| j >= n || w == Ratio::zero())
        {
            return Err("invalid homogeneous multilevel mounting");
        }
    }
    for (k, &root) in lift.free_roots.iter().enumerate() {
        if lift
            .lines
            .get(root)
            .is_none_or(|line| line.terms != [(k, Ratio::one())])
        {
            return Err("inconsistent multilevel free roots");
        }
    }
    Ok(())
}

fn mounted_indices(space: MountedSpace<'_>) -> Result<Vec<Option<usize>>, &'static str> {
    let full = &space.adaptive.lift;
    let mount = space.lift;
    if full.lines.len() != mount.lines.len()
        || mount
            .free_roots
            .iter()
            .any(|r| full.free_roots.binary_search(r).is_err())
    {
        return Err("multilevel mounting is not a restriction of adaptive roots");
    }
    let indices: Vec<_> = full
        .free_roots
        .iter()
        .map(|root| mount.free_roots.binary_search(root).ok())
        .collect();
    for (adaptive, mounted) in full.lines.iter().zip(&mount.lines) {
        if adaptive.offset != Ratio::zero()
            || adaptive
                .terms
                .iter()
                .filter_map(|&(j, w)| indices.get(j).copied().flatten().map(|k| (k, w)))
                .ne(mounted.terms.iter().copied())
        {
            return Err("mounting changes an unconstrained interpolation stencil");
        }
    }
    Ok(indices)
}

impl MultilevelPreconditioner {
    /// Factor the complete mounted coarse matrix, subject to an explicit byte cap.
    pub fn new(k: &Stiffness, lift: &Lift, max_factor_bytes: usize) -> Result<Self, &'static str> {
        validate(lift, k)?;
        let n = 3 * lift.free_roots.len();
        let size = n
            .checked_mul(n)
            .filter(|s| *s <= max_factor_bytes / 8)
            .ok_or("coarse factor exceeds memory cap")?;
        let mut a = vec![0.; size];
        for (i, row) in lift.lines.iter().enumerate() {
            for p in k.offsets[i]..k.offsets[i + 1] {
                let col = &lift.lines[k.columns[p]];
                for &(ii, wi) in &row.terms {
                    for &(jj, wj) in &col.terms {
                        let weight = coefficient(wi) * coefficient(wj);
                        for c in 0..3 {
                            for d in 0..3 {
                                a[(3 * ii + c) * n + 3 * jj + d] += weight * k.values[p][c][d];
                            }
                        }
                    }
                }
            }
        }
        let mut first_column = vec![0; n];
        for i in 0..n {
            let row = i * n;
            first_column[i] = (0..i).find(|&j| a[row + j] != 0.).unwrap_or(i);
            for j in 0..i {
                let x = a[row + j];
                let y = a[j * n + i];
                let scale = a[row + i].abs().sqrt() * a[j * n + j].abs().sqrt();
                if !x.is_finite() || !y.is_finite() || (x - y).abs() > 2e-12 * scale.max(1.) {
                    return Err("asymmetric coarse stiffness");
                }
            }
        }
        // The lower skyline of a symmetric matrix contains its Cholesky fill.
        for i in 0..n {
            for j in first_column[i]..=i {
                let mut value = a[i * n + j];
                for q in first_column[i].max(first_column[j])..j {
                    value -= a[i * n + q] * a[j * n + q];
                }
                if i == j {
                    if !value.is_finite() || value <= 0. {
                        return Err("coarse stiffness is not positive definite");
                    }
                    a[i * n + j] = value.sqrt();
                } else {
                    a[i * n + j] = value / a[j * n + j];
                }
            }
        }
        Ok(Self {
            coarse: a,
            first_column,
            levels: Vec::new(),
            last_mount: lift.clone(),
            last_mesh: None,
        })
    }

    /// Append a nested refinement after checking the homogeneous mounting commutes with prolongation.
    pub fn add_level(
        &mut self,
        source: MountedSpace<'_>,
        target: MountedSpace<'_>,
        k: &Stiffness,
    ) -> Result<(), &'static str> {
        if &self.last_mount != source.lift
            || self
                .last_mesh
                .as_ref()
                .is_some_and(|previous| previous != source.adaptive)
        {
            return Err("multilevel hierarchy source mismatch");
        }
        validate(target.lift, k)?;
        let s = mounted_indices(source)?;
        let t = mounted_indices(target)?;
        let transfer = build_transfer(
            source.adaptive,
            source.partition,
            target.adaptive,
            target.partition,
        )?;
        let mut rows = Vec::new();
        for (row, index) in transfer.rows.iter().zip(t) {
            let entries: Vec<_> = row
                .iter()
                .filter_map(|&(j, w)| s[j].map(|k| (k, coefficient(w))))
                .collect();
            if let Some(i) = index {
                if rows.len() != i {
                    return Err("unordered multilevel target mounting");
                }
                rows.push(entries);
            } else if !entries.is_empty() {
                return Err("prolongation violates homogeneous mounting");
            }
        }
        let nt = rows.len();
        let mut diag = vec![0.; 3 * nt];
        for (i, row) in target.lift.lines.iter().enumerate() {
            for p in k.offsets[i]..k.offsets[i + 1] {
                let col = &target.lift.lines[k.columns[p]].terms;
                for &(master, w) in &row.terms {
                    if let Ok(j) = col.binary_search_by_key(&master, |x| x.0) {
                        let weight = coefficient(w) * coefficient(col[j].1);
                        for c in 0..3 {
                            diag[c * nt + master] += weight * k.values[p][c][c];
                        }
                    }
                }
            }
        }
        for x in &mut diag {
            if !x.is_finite() || *x <= 0. {
                return Err("nonpositive reduced multilevel diagonal");
            }
            *x = 1. / *x;
        }
        self.levels.push(Level {
            rows,
            source_count: source.lift.free_roots.len(),
            inverse_diagonal: diag,
        });
        self.last_mount = target.lift.clone();
        self.last_mesh = Some(target.adaptive.clone());
        Ok(())
    }

    pub fn level_count(&self) -> usize {
        1 + self.levels.len()
    }
    pub fn transfer_nonzeros(&self) -> usize {
        self.levels.iter().flat_map(|l| &l.rows).map(Vec::len).sum()
    }

    fn cycle(&self, depth: usize, r: &[f64]) -> Vec<f64> {
        if depth == 0 {
            let n = self.first_column.len();
            let m = n / 3;
            let mut x: Vec<_> = (0..n).map(|i| r[(i % 3) * m + i / 3]).collect();
            for i in 0..n {
                for j in self.first_column[i]..i {
                    x[i] -= self.coarse[i * n + j] * x[j];
                }
                x[i] /= self.coarse[i * n + i];
            }
            for i in (0..n).rev() {
                x[i] /= self.coarse[i * n + i];
                for j in self.first_column[i]..i {
                    x[j] -= self.coarse[i * n + j] * x[i];
                }
            }
            return (0..n).map(|i| x[3 * (i % m) + i / m]).collect();
        }
        let l = &self.levels[depth - 1];
        let ns = l.source_count;
        let nt = l.rows.len();
        let mut restricted = vec![0.; 3 * ns];
        for (i, row) in l.rows.iter().enumerate() {
            for &(j, w) in row {
                for c in 0..3 {
                    restricted[c * ns + j] += w * r[c * nt + i];
                }
            }
        }
        let coarse = self.cycle(depth - 1, &restricted);
        let mut z: Vec<_> = r
            .iter()
            .zip(&l.inverse_diagonal)
            .map(|(r, d)| r * d)
            .collect();
        for (i, row) in l.rows.iter().enumerate() {
            for &(j, w) in row {
                for c in 0..3 {
                    z[c * nt + i] += w * coarse[c * ns + j];
                }
            }
        }
        z
    }
}
impl ReducedPreconditioner for MultilevelPreconditioner {
    fn apply(&self, residual: &[f64]) -> Result<Vec<f64>, &'static str> {
        if residual.len() != 3 * self.last_mount.free_roots.len()
            || residual.iter().any(|x| !x.is_finite())
        {
            return Err("invalid multilevel residual");
        }
        let z = self.cycle(self.levels.len(), residual);
        if z.iter().any(|x| !x.is_finite()) {
            return Err("nonfinite multilevel correction");
        }
        Ok(z)
    }
}
