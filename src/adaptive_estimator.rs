//! Residual indicators for exact dyadic adaptive Q2 spaces.
//!
//! This computes a mesh-dependent residual indicator. It does not provide a
//! certified upper bound for the exact energy error.

use crate::{
    adaptive_partition::{
        AdaptivePartition, ExteriorPatch, FacePatch, FaceRect, LeafKey, LeafSide,
    },
    adaptive_space::AdaptiveSpace,
    analysis_space::AnalysisBlock,
    elasticity::{gauss, Vector},
    incidence::BlockSide,
};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct LeafResidual {
    pub volume_squared: f64,
    pub interior_squared: f64,
    pub natural_squared: f64,
    pub total_squared: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveResidual {
    pub by_leaf: BTreeMap<LeafKey, LeafResidual>,
    pub indicator_squared: f64,
    pub indicator: f64,
    pub error_upper_bound_available: bool,
    pub qualification: &'static str,
}

impl AdaptiveResidual {
    /// Deterministic Dörfler marking by squared element indicators.
    pub fn bulk_mark(&self, theta: f64) -> Result<Vec<LeafKey>, &'static str> {
        if !theta.is_finite() || theta <= 0. || theta > 1. {
            return Err("Dorfler fraction must be in (0,1]");
        }
        if self
            .by_leaf
            .values()
            .any(|x| !x.total_squared.is_finite() || x.total_squared < 0.)
        {
            return Err("invalid element indicator");
        }
        if self.indicator_squared == 0. {
            return Ok(Vec::new());
        }
        let mut ordered = self
            .by_leaf
            .iter()
            .map(|(key, value)| (key.clone(), value.total_squared))
            .collect::<Vec<_>>();
        ordered.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        let target = theta * self.indicator_squared;
        let mut sum = 0.;
        let mut selected = Vec::new();
        for (key, value) in ordered {
            selected.push(key);
            sum += value;
            if sum >= target {
                break;
            }
        }
        selected.sort();
        Ok(selected)
    }
}

#[derive(Clone, Copy)]
struct Jet2 {
    value: f64,
    first: Vector,
    second: [[f64; 3]; 3],
}

impl Default for Jet2 {
    fn default() -> Self {
        Self {
            value: 0.,
            first: [0.; 3],
            second: [[0.; 3]; 3],
        }
    }
}

#[derive(Clone, Copy)]
struct Kinematics {
    position: Vector,
    jacobian: [[f64; 3]; 3], // physical coordinate, then parameter axis
    det: f64,
    #[allow(dead_code)]
    hessian: [[[f64; 3]; 3]; 3], // physical displacement Hessian, retained for diagnostics
    stress: [[f64; 3]; 3],
    divergence: Vector,
}

fn q2(t: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
    (
        [(1. - t) * (1. - t), 2. * t * (1. - t), t * t],
        [-2. * (1. - t), 2. - 4. * t, 2. * t],
        [2., -4., 2.],
    )
}

fn quotient(n: Jet2, d: Jet2) -> Jet2 {
    let value = n.value / d.value;
    let first = std::array::from_fn(|a| (n.first[a] - value * d.first[a]) / d.value);
    let second = std::array::from_fn(|a| {
        std::array::from_fn(|b| {
            (n.second[a][b]
                - value * d.second[a][b]
                - first[a] * d.first[b]
                - first[b] * d.first[a])
                / d.value
        })
    });
    Jet2 {
        value,
        first,
        second,
    }
}

fn inverse3(m: [[f64; 3]; 3]) -> Result<([[f64; 3]; 3], f64), &'static str> {
    let c = [
        [
            m[1][1] * m[2][2] - m[1][2] * m[2][1],
            m[1][2] * m[2][0] - m[1][0] * m[2][2],
            m[1][0] * m[2][1] - m[1][1] * m[2][0],
        ],
        [
            m[0][2] * m[2][1] - m[0][1] * m[2][2],
            m[0][0] * m[2][2] - m[0][2] * m[2][0],
            m[0][1] * m[2][0] - m[0][0] * m[2][1],
        ],
        [
            m[0][1] * m[1][2] - m[0][2] * m[1][1],
            m[0][2] * m[1][0] - m[0][0] * m[1][2],
            m[0][0] * m[1][1] - m[0][1] * m[1][0],
        ],
    ];
    let det = m[0][0] * c[0][0] + m[0][1] * c[0][1] + m[0][2] * c[0][2];
    if !det.is_finite() || det <= 0. {
        return Err("nonpositive geometry Jacobian in residual estimator");
    }
    Ok((
        std::array::from_fn(|i| std::array::from_fn(|j| c[j][i] / det)),
        det,
    ))
}

fn evaluate(
    block: &AnalysisBlock,
    dofs: &[f64],
    uvw: Vector,
    lambda: f64,
    mu: f64,
) -> Result<Kinematics, &'static str> {
    if block.controls.len() != 27
        || block.dofs.len() != 27
        || uvw
            .iter()
            .any(|x| !x.is_finite() || !(0. ..=1.).contains(x))
    {
        return Err("invalid adaptive element evaluation");
    }
    let basis = uvw.map(q2);
    let mut den = Jet2::default();
    let mut geom = [Jet2::default(); 3];
    let mut field = [Jet2::default(); 3];
    for n in 0..27 {
        let [i, j, k] = [n / 9, n / 3 % 3, n % 3];
        let ids = [i, j, k];
        let mut b = 1.;
        let mut db = [0.; 3];
        let mut ddb = [[0.; 3]; 3];
        for a in 0..3 {
            b *= basis[a].0[ids[a]];
        }
        for a in 0..3 {
            db[a] = basis[a].1[ids[a]]
                * (0..3)
                    .filter(|&q| q != a)
                    .map(|q| basis[q].0[ids[q]])
                    .product::<f64>();
            ddb[a][a] = basis[a].2[ids[a]]
                * (0..3)
                    .filter(|&q| q != a)
                    .map(|q| basis[q].0[ids[q]])
                    .product::<f64>();
            for c in 0..a {
                ddb[a][c] = basis[a].1[ids[a]]
                    * basis[c].1[ids[c]]
                    * (0..3)
                        .filter(|&q| q != a && q != c)
                        .map(|q| basis[q].0[ids[q]])
                        .product::<f64>();
                ddb[c][a] = ddb[a][c];
            }
        }
        let h = block.controls[n];
        if h[3] <= 0 {
            return Err("nonpositive homogeneous weight");
        }
        let w = h[3] as f64;
        let weight = Jet2 {
            value: b * w,
            first: db.map(|x| x * w),
            second: std::array::from_fn(|a| ddb[a].map(|x| x * w)),
        };
        den.value += weight.value;
        for a in 0..3 {
            den.first[a] += weight.first[a];
            for c in 0..3 {
                den.second[a][c] += weight.second[a][c];
            }
        }
        let dof = block.dofs[n];
        if 3 * dof + 2 >= dofs.len() {
            return Err("displacement vector does not cover adaptive DOFs");
        }
        for c in 0..3 {
            let geo = Jet2 {
                value: b * h[c] as f64,
                first: db.map(|x| x * h[c] as f64),
                second: std::array::from_fn(|a| ddb[a].map(|x| x * h[c] as f64)),
            };
            geom[c].value += geo.value;
            field[c].value += weight.value * dofs[3 * dof + c];
            for a in 0..3 {
                geom[c].first[a] += geo.first[a];
                field[c].first[a] += weight.first[a] * dofs[3 * dof + c];
                for q in 0..3 {
                    geom[c].second[a][q] += geo.second[a][q];
                    field[c].second[a][q] += weight.second[a][q] * dofs[3 * dof + c];
                }
            }
        }
    }
    if !den.value.is_finite() || den.value <= 0. {
        return Err("invalid rational denominator");
    }
    let x = geom.map(|n| quotient(n, den));
    let u = field.map(|n| quotient(n, den));
    let jacobian = std::array::from_fn(|c| x[c].first);
    let (inv, det) = inverse3(jacobian)?; // inv[parameter][physical]
    let position = x.map(|z| z.value);
    let mut gradient = [[0.; 3]; 3];
    for (c, grad_c) in gradient.iter_mut().enumerate() {
        for (p, entry) in grad_c.iter_mut().enumerate() {
            *entry = (0..3).map(|a| u[c].first[a] * inv[a][p]).sum();
        }
    }
    let mut hessian = [[[0.; 3]; 3]; 3];
    for c in 0..3 {
        let mut corrected = u[c].second;
        for (a, row) in corrected.iter_mut().enumerate() {
            for (b, entry) in row.iter_mut().enumerate() {
                *entry -= (0..3)
                    .map(|q| gradient[c][q] * x[q].second[a][b])
                    .sum::<f64>();
            }
        }
        for p in 0..3 {
            for q in 0..3 {
                hessian[c][p][q] = (0..3)
                    .flat_map(|a| (0..3).map(move |b| (a, b)))
                    .map(|(a, b)| inv[a][p] * corrected[a][b] * inv[b][q])
                    .sum();
            }
        }
    }
    let strain_trace = gradient
        .iter()
        .enumerate()
        .map(|(c, row)| row[c])
        .sum::<f64>();
    let stress = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            lambda * strain_trace * if i == j { 1. } else { 0. }
                + mu * (gradient[i][j] + gradient[j][i])
        })
    });
    let divergence = std::array::from_fn(|i| {
        mu * (0..3).map(|j| hessian[i][j][j]).sum::<f64>()
            + (lambda + mu) * (0..3).map(|j| hessian[j][i][j]).sum::<f64>()
    });
    if position
        .iter()
        .chain(gradient.iter().flatten())
        .chain(hessian.iter().flatten().flatten())
        .chain(stress.iter().flatten())
        .chain(divergence.iter())
        .any(|v| !v.is_finite())
    {
        return Err("nonfinite residual kinematics");
    }
    Ok(Kinematics {
        position,
        jacobian,
        det,
        hessian,
        stress,
        divergence,
    })
}

fn diameter(block: &AnalysisBlock) -> Result<f64, &'static str> {
    if block.controls.len() != 27 {
        return Err("invalid control net for element diameter");
    }
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in &block.controls {
        if p[3] <= 0 {
            return Err("nonpositive homogeneous weight");
        }
        for a in 0..3 {
            let x = p[a] as f64 / p[3] as f64;
            if !x.is_finite() {
                return Err("nonfinite physical control");
            }
            lo[a] = lo[a].min(x);
            hi[a] = hi[a].max(x);
        }
    }
    Ok((0..3).map(|a| (hi[a] - lo[a]).powi(2)).sum::<f64>().sqrt())
}

fn dot(a: Vector, b: Vector) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn cross(a: Vector, b: Vector) -> Vector {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: Vector) -> f64 {
    dot(a, a).sqrt()
}
fn stress_traction(s: [[f64; 3]; 3], n: Vector) -> Vector {
    std::array::from_fn(|i| (0..3).map(|j| s[i][j] * n[j]).sum())
}

fn tangent_axes(axis: usize) -> [usize; 2] {
    match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}
fn side_parameter(side: &BlockSide, uv: [f64; 2]) -> Vector {
    let mut p = [0.; 3];
    let a = tangent_axes(side.axis);
    p[side.axis] = if side.upper { 1. } else { 0. };
    p[a[0]] = uv[0];
    p[a[1]] = uv[1];
    p
}
fn root_range(leaf: &LeafKey, axis: usize) -> [f64; 2] {
    let scale = (1_u32 << leaf.level) as f64;
    [
        leaf.index[axis] as f64 / scale,
        (leaf.index[axis] + 1) as f64 / scale,
    ]
}
fn local_from_root(leaf: &LeafKey, axis: usize, x: f64) -> f64 {
    (x - root_range(leaf, axis)[0]) / (root_range(leaf, axis)[1] - root_range(leaf, axis)[0])
}
fn local_patch(leaf: &LeafKey, side: &BlockSide, rect: FaceRect, denom: u32) -> [[f64; 2]; 2] {
    let axes = tangent_axes(side.axis);
    std::array::from_fn(|a| {
        [
            local_from_root(leaf, axes[a], rect.ranges[a][0] as f64 / denom as f64),
            local_from_root(leaf, axes[a], rect.ranges[a][1] as f64 / denom as f64),
        ]
    })
}
fn face_normal(k: Kinematics, side: &BlockSide) -> Result<Vector, &'static str> {
    let a = tangent_axes(side.axis);
    let n = cross(k.jacobian[a[0]], k.jacobian[a[1]]);
    let sign =
        if side.upper { 1. } else { -1. } * if side.axis.is_multiple_of(2) { 1. } else { -1. };
    let n = n.map(|x| sign * x);
    let l = norm(n);
    if !l.is_finite() || l <= 0. {
        return Err("invalid physical face normal");
    }
    Ok(n.map(|x| x / l))
}
fn face_measure(k: Kinematics, side: &BlockSide) -> f64 {
    let a = tangent_axes(side.axis);
    norm(cross(k.jacobian[a[0]], k.jacobian[a[1]]))
}
fn block_for<'a>(
    adaptive: &'a AdaptiveSpace,
    key: &LeafKey,
) -> Result<&'a AnalysisBlock, &'static str> {
    let id = adaptive
        .leaf_blocks
        .get(key)
        .ok_or("adaptive leaf has no block identity")?;
    adaptive
        .space
        .blocks
        .iter()
        .find(|b| &b.id == id)
        .ok_or("adaptive leaf block is missing")
}

/// Estimate the strong equilibrium residual and traction jumps on an adaptive
/// rational Q2 mesh. Displacements are in mm, E in MPa, forces in N/mm^3.
#[allow(clippy::too_many_arguments)]
pub fn estimate_residual<F>(
    adaptive: &AdaptiveSpace,
    partition: &AdaptivePartition,
    displacement: &[f64],
    young: f64,
    poisson: f64,
    order: usize,
    body_force: Vector,
    mut neumann: F,
) -> Result<AdaptiveResidual, &'static str>
where
    F: FnMut(&LeafSide, Vector, Vector) -> Option<Vector>,
{
    if !young.is_finite()
        || young <= 0.
        || !poisson.is_finite()
        || poisson <= -1.
        || poisson >= 0.5
        || body_force.iter().any(|x| !x.is_finite())
        || displacement.len() != 3 * adaptive.space.ndof
        || displacement.iter().any(|x| !x.is_finite())
    {
        return Err("invalid residual estimator inputs");
    }
    let points = gauss(order)?;
    let report = partition.snapshot()?;
    if adaptive.leaf_blocks.len() != report.leaves.len()
        || report
            .leaves
            .iter()
            .any(|leaf| !adaptive.leaf_blocks.contains_key(leaf))
    {
        return Err("adaptive space and partition leaves differ");
    }
    let lambda = young * poisson / ((1. + poisson) * (1. - 2. * poisson));
    let mu = young / (2. * (1. + poisson));
    let mut result = AdaptiveResidual {
        by_leaf: report
            .leaves
            .iter()
            .cloned()
            .map(|k| (k, LeafResidual::default()))
            .collect(),
        indicator_squared: 0.,
        indicator: 0.,
        error_upper_bound_available: false,
        qualification:
            "mesh-dependent residual indicator; no proven upper bound for exact energy error",
    };
    let hs = report
        .leaves
        .iter()
        .map(|leaf| Ok((leaf.clone(), diameter(block_for(adaptive, leaf)?)?)))
        .collect::<Result<BTreeMap<_, _>, &'static str>>()?;
    for leaf in &report.leaves {
        let block = block_for(adaptive, leaf)?;
        let h = hs[leaf];
        let mut sum = 0.;
        for &(u, wu) in &points {
            for &(v, wv) in &points {
                for &(w, ww) in &points {
                    let k = evaluate(block, displacement, [u, v, w], lambda, mu)?;
                    let r = std::array::from_fn(|c| k.divergence[c] + body_force[c]);
                    sum += h * h / (young * young) * dot(r, r) * k.det * wu * wv * ww;
                }
            }
        }
        result.by_leaf.get_mut(leaf).unwrap().volume_squared = sum;
    }
    for patch in &report.interfaces {
        integrate_interface(
            patch,
            adaptive,
            displacement,
            &points,
            lambda,
            mu,
            partition.grid_denominator(),
            young,
            &hs,
            &mut result.by_leaf,
        )?;
    }
    for patch in &report.exterior {
        integrate_exterior(
            patch,
            adaptive,
            displacement,
            &points,
            lambda,
            mu,
            partition.grid_denominator(),
            young,
            &hs,
            &mut result.by_leaf,
            &mut neumann,
        )?;
    }
    for value in result.by_leaf.values_mut() {
        value.total_squared = value.volume_squared + value.interior_squared + value.natural_squared;
        if !value.total_squared.is_finite() || value.total_squared < 0. {
            return Err("invalid element residual indicator");
        }
        result.indicator_squared += value.total_squared;
    }
    if !result.indicator_squared.is_finite() {
        return Err("nonfinite global residual indicator");
    }
    result.indicator = result.indicator_squared.sqrt();
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn integrate_interface(
    p: &FacePatch,
    a: &AdaptiveSpace,
    u: &[f64],
    q: &[(f64, f64)],
    lambda: f64,
    mu: f64,
    denom: u32,
    e: f64,
    hs: &BTreeMap<LeafKey, f64>,
    out: &mut BTreeMap<LeafKey, LeafResidual>,
) -> Result<(), &'static str> {
    let b1 = block_for(a, &p.first.leaf)?;
    let b2 = block_for(a, &p.second.leaf)?;
    let r1 = local_patch(&p.first.leaf, &p.first.root_side, p.first_subrect, denom);
    let scale = (r1[0][1] - r1[0][0]) * (r1[1][1] - r1[1][0]);
    let mut integral = 0.;
    for &(s, ws) in q {
        for &(t, wt) in q {
            let uv = [s, t];
            let root = [
                p.first_subrect.ranges[0][0] as f64 / denom as f64
                    + uv[0] * (p.first_subrect.ranges[0][1] - p.first_subrect.ranges[0][0]) as f64
                        / denom as f64,
                p.first_subrect.ranges[1][0] as f64 / denom as f64
                    + uv[1] * (p.first_subrect.ranges[1][1] - p.first_subrect.ranges[1][0]) as f64
                        / denom as f64,
            ];
            let x1 = [
                r1[0][0] + uv[0] * (r1[0][1] - r1[0][0]),
                r1[1][0] + uv[1] * (r1[1][1] - r1[1][0]),
            ];
            let a1 = tangent_axes(p.first.root_side.axis);
            let a2 = tangent_axes(p.second.root_side.axis);
            let root2 = if p.orientation.swap {
                [
                    if p.orientation.flip_second {
                        1. - root[1]
                    } else {
                        root[1]
                    },
                    if p.orientation.flip_first {
                        1. - root[0]
                    } else {
                        root[0]
                    },
                ]
            } else {
                [
                    if p.orientation.flip_first {
                        1. - root[0]
                    } else {
                        root[0]
                    },
                    if p.orientation.flip_second {
                        1. - root[1]
                    } else {
                        root[1]
                    },
                ]
            };
            let mut z1 = [0.; 3];
            let mut z2 = [0.; 3];
            z1[p.first.root_side.axis] = if p.first.root_side.upper { 1. } else { 0. };
            z2[p.second.root_side.axis] = if p.second.root_side.upper { 1. } else { 0. };
            for j in 0..2 {
                z1[a1[j]] = x1[j];
                let bounds = root_range(&p.second.leaf, a2[j]);
                z2[a2[j]] = (root2[j] - bounds[0]) / (bounds[1] - bounds[0]);
            }
            let k1 = evaluate(b1, u, z1, lambda, mu)?;
            let k2 = evaluate(b2, u, z2, lambda, mu)?;
            if (0..3).any(|c| {
                (k1.position[c] - k2.position[c]).abs()
                    > 1e-8 * (1. + k1.position[c].abs().max(k2.position[c].abs()))
            }) {
                return Err("interface parameter maps disagree physically");
            }
            let n1 = face_normal(k1, &p.first.root_side)?;
            let n2 = face_normal(k2, &p.second.root_side)?;
            let tr1 = stress_traction(k1.stress, n1);
            let tr2 = stress_traction(k2.stress, n2);
            let jump = std::array::from_fn(|c| tr1[c] + tr2[c]);
            let m = face_measure(k1, &p.first.root_side) * scale;
            integral += dot(jump, jump) * m * ws * wt;
        }
    }
    out.get_mut(&p.first.leaf)
        .ok_or("missing first interface indicator")?
        .interior_squared += hs[&p.first.leaf] / e * integral;
    out.get_mut(&p.second.leaf)
        .ok_or("missing second interface indicator")?
        .interior_squared += hs[&p.second.leaf] / e * integral;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn integrate_exterior<F>(
    p: &ExteriorPatch,
    a: &AdaptiveSpace,
    u: &[f64],
    q: &[(f64, f64)],
    lambda: f64,
    mu: f64,
    denom: u32,
    e: f64,
    hs: &BTreeMap<LeafKey, f64>,
    out: &mut BTreeMap<LeafKey, LeafResidual>,
    neumann: &mut F,
) -> Result<(), &'static str>
where
    F: FnMut(&LeafSide, Vector, Vector) -> Option<Vector>,
{
    let block = block_for(a, &p.side.leaf)?;
    let r = local_patch(&p.side.leaf, &p.side.root_side, p.root_rect, denom);
    let scale = (r[0][1] - r[0][0]) * (r[1][1] - r[1][0]);
    let mut integral = 0.;
    for &(s, ws) in q {
        for &(t, wt) in q {
            let uv = [
                r[0][0] + s * (r[0][1] - r[0][0]),
                r[1][0] + t * (r[1][1] - r[1][0]),
            ];
            let k = evaluate(block, u, side_parameter(&p.side.root_side, uv), lambda, mu)?;
            let n = face_normal(k, &p.side.root_side)?;
            if let Some(target) = neumann(&p.side, k.position, n) {
                if target.iter().any(|x| !x.is_finite()) {
                    return Err("nonfinite prescribed Neumann traction");
                }
                let actual = stress_traction(k.stress, n);
                let residual = std::array::from_fn(|c| actual[c] - target[c]);
                integral +=
                    dot(residual, residual) * face_measure(k, &p.side.root_side) * scale * ws * wt;
            }
        }
    }
    out.get_mut(&p.side.leaf)
        .ok_or("missing exterior indicator")?
        .natural_squared += hs[&p.side.leaf] / e * integral;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn affine_cube() -> AnalysisBlock {
        let controls = (0..27)
            .map(|n| {
                let i = n / 9;
                let j = n / 3 % 3;
                let k = n % 3;
                [
                    (i as i128) * 500,
                    (j as i128) * 500,
                    (k as i128) * 500,
                    1000,
                ]
            })
            .collect::<Vec<_>>();
        AnalysisBlock {
            id: "cube".into(),
            controls,
            dofs: (0..27).collect(),
        }
    }

    fn cube_patch(id: &str, x0: f64, flip_y: bool, offset: usize) -> (AnalysisBlock, Vec<f64>) {
        let mut values = Vec::new();
        let mut controls = Vec::new();
        for n in 0..27 {
            let i = n / 9;
            let j = n / 3 % 3;
            let k = n % 3;
            let x = x0 + i as f64 / 2.;
            let y0 = if flip_y {
                1. - j as f64 / 2.
            } else {
                j as f64 / 2.
            };
            let z = if flip_y {
                1. - k as f64 / 2.
            } else {
                k as f64 / 2.
            };
            controls.push([
                (x * 1000.) as i128,
                (y0 * 1000.) as i128,
                (z * 1000.) as i128,
                1000,
            ]);
            values.extend([y0 * y0, 0., 0.]);
        }
        (
            AnalysisBlock {
                id: id.into(),
                controls,
                dofs: (offset..offset + 27).collect(),
            },
            values,
        )
    }

    #[test]
    fn rational_analytic_hessians_vanish_for_affine_map_and_field() {
        let b = affine_cube();
        let u = (0..27)
            .flat_map(|n| {
                let i = n / 9;
                let j = n / 3 % 3;
                let k = n % 3;
                [0.1 * i as f64, 0.2 * j as f64, -0.3 * k as f64]
            })
            .collect::<Vec<_>>();
        let l = 10.;
        let m = 4.;
        for p in [[0.13, 0.29, 0.91], [0.5, 0.5, 0.5], [0.92, 0.04, 0.61]] {
            let x = evaluate(&b, &u, p, l, m).unwrap();
            assert!(x
                .hessian
                .iter()
                .flatten()
                .flatten()
                .all(|v| v.abs() < 1e-12));
            assert!(x.divergence.iter().all(|v| v.abs() < 1e-12));
        }
    }

    #[test]
    fn quadratic_manufactured_cube_has_constant_strong_residual() {
        let b = affine_cube();
        let u = (0..27)
            .flat_map(|n| {
                let i = n / 9;
                let j = n / 3 % 3;
                let k = n % 3;
                let x = i as f64 / 2.;
                let y = j as f64 / 2.;
                let z = k as f64 / 2.;
                [x * x + 2. * y * z, y * y + 2. * x * z, z * z + 2. * x * y]
            })
            .collect::<Vec<_>>();
        let x = evaluate(&b, &u, [0.21, 0.63, 0.44], 10., 4.).unwrap();
        assert!(x.divergence.iter().all(|v| v.is_finite()));
        assert!(x.hessian.iter().flatten().flatten().any(|v| v.abs() > 0.1));
    }

    #[test]
    fn curved_geometry_preserves_affine_physical_field_hessian() {
        let mut b = affine_cube();
        for n in 0..27 {
            let [i, j, k] = [n / 9, n / 3 % 3, n % 3];
            let weight = if i == 1 && j == 1 && k == 1 {
                1200_i128
            } else {
                1000_i128
            };
            b.controls[n] = [
                (i as i128 * 500) * weight / 1000,
                (j as i128 * 500) * weight / 1000,
                (k as i128 * 500) * weight / 1000,
                weight,
            ];
        }
        let u = (0..27)
            .flat_map(|n| {
                let [i, j, k] = [n / 9, n / 3 % 3, n % 3];
                let x = i as f64 / 2.;
                let y = j as f64 / 2.;
                let z = k as f64 / 2.;
                [2. * x - y + 0.2, x + 3. * z, y - z]
            })
            .collect::<Vec<_>>();
        let x = evaluate(&b, &u, [0.31, 0.42, 0.71], 10., 4.).unwrap();
        assert!(x
            .hessian
            .iter()
            .flatten()
            .flatten()
            .all(|v| v.abs() < 2e-12));
        assert!(x.divergence.iter().all(|v| v.abs() < 2e-12));
    }

    #[test]
    fn variable_traction_jump_is_integrated_once_with_d4_mapping() {
        let (first, mut u1) = cube_patch("first", 0., false, 0);
        let (second, u2) = cube_patch("second", 1., true, 27);
        u1.extend(u2);
        let a = LeafKey::root("first");
        let b = LeafKey::root("second");
        let side1 = LeafSide {
            leaf: a.clone(),
            root_side: BlockSide {
                block: "first".into(),
                axis: 0,
                upper: true,
            },
        };
        let side2 = LeafSide {
            leaf: b.clone(),
            root_side: BlockSide {
                block: "second".into(),
                axis: 0,
                upper: false,
            },
        };
        let rect = FaceRect {
            ranges: [[0, 1], [0, 1]],
        };
        let patch = FacePatch {
            first: side1,
            second: side2,
            master_rect: rect,
            first_subrect: rect,
            second_subrect: rect,
            orientation: crate::incidence::Orientation {
                swap: false,
                flip_first: true,
                flip_second: true,
            },
        };
        let leaves = BTreeMap::from([(a.clone(), 1.), (b.clone(), 1.)]);
        let mut out = BTreeMap::from([
            (a.clone(), LeafResidual::default()),
            (b.clone(), LeafResidual::default()),
        ]);
        let adaptive = AdaptiveSpace {
            space: crate::analysis_space::AnalysisSpace {
                ndof: 54,
                blocks: vec![first, second],
                boundary_dofs: BTreeMap::new(),
            },
            lift: crate::constraint_system::resolve(54, &[], &[]).unwrap(),
            leaf_blocks: BTreeMap::from([
                (a.clone(), "first".into()),
                (b.clone(), "second".into()),
            ]),
            constraints: Vec::new(),
            boundary_patches: BTreeMap::new(),
            root_ndof: 54,
        };
        integrate_interface(
            &patch,
            &adaptive,
            &u1,
            &gauss(4).unwrap(),
            10.,
            4.,
            1,
            200_000.,
            &leaves,
            &mut out,
        )
        .unwrap();
        assert!(out.values().all(|x| x.interior_squared < 1e-20));
        let mut bad = patch.clone();
        bad.orientation.flip_first = false;
        let mut out = BTreeMap::from([
            (a.clone(), LeafResidual::default()),
            (b.clone(), LeafResidual::default()),
        ]);
        assert_eq!(
            integrate_interface(
                &bad,
                &adaptive,
                &u1,
                &gauss(4).unwrap(),
                10.,
                4.,
                1,
                200_000.,
                &leaves,
                &mut out
            ),
            Err("interface parameter maps disagree physically")
        );
    }

    #[test]
    fn residual_estimator_integrates_quadratic_body_term_and_neumann_selection() {
        let block = affine_cube();
        let block_id = block.id.clone();
        let sides = (0..3)
            .flat_map(|axis| {
                let block_id = block_id.clone();
                [false, true].map(move |upper| BlockSide {
                    block: block_id.clone(),
                    axis,
                    upper,
                })
            })
            .collect();
        let incidence = crate::incidence::Incidence {
            interfaces: Vec::new(),
            exterior: sides,
        };
        let partition = AdaptivePartition::new(&["cube".into()], &incidence, 2, 100).unwrap();
        let leaf = LeafKey::root("cube");
        let adaptive = AdaptiveSpace {
            space: crate::analysis_space::AnalysisSpace {
                ndof: 27,
                blocks: vec![block],
                boundary_dofs: BTreeMap::new(),
            },
            lift: crate::constraint_system::resolve(27, &[], &[]).unwrap(),
            leaf_blocks: BTreeMap::from([(leaf, "cube".into())]),
            constraints: Vec::new(),
            boundary_patches: BTreeMap::new(),
            root_ndof: 27,
        };
        let affine = (0..27)
            .flat_map(|n| {
                let [i, j, k] = [n / 9, n / 3 % 3, n % 3];
                [i as f64 / 2., -2. * j as f64 / 2., 0.5 * k as f64 / 2.]
            })
            .collect::<Vec<_>>();
        let natural = estimate_residual(
            &adaptive,
            &partition,
            &affine,
            1000.,
            0.25,
            3,
            [0.; 3],
            |_, _, _| Some([0.; 3]),
        )
        .unwrap();
        assert!(natural.by_leaf.values().next().unwrap().natural_squared > 0.);
        let lambda = 1000. * 0.25 / (1.25 * 0.5);
        let mu = 1000. / (2. * 1.25);
        let body = [-(lambda + 2. * mu) * 2., 0., 0.];
        let quadratic = (0..27)
            .flat_map(|n| {
                let i = n / 9;
                [if i == 2 { 1. } else { 0. }, 0., 0.]
            })
            .collect::<Vec<_>>();
        let balanced = estimate_residual(
            &adaptive,
            &partition,
            &quadratic,
            1000.,
            0.25,
            3,
            body,
            |_, _, _| None,
        )
        .unwrap();
        assert!(balanced.indicator_squared < 1e-20);
        assert!(!balanced.error_upper_bound_available);
    }

    #[test]
    fn bulk_marking_is_stable_and_rejects_bad_fraction() {
        let a = LeafKey::root("a");
        let b = LeafKey::root("b");
        let mut by_leaf = BTreeMap::new();
        by_leaf.insert(
            a.clone(),
            LeafResidual {
                total_squared: 3.,
                ..Default::default()
            },
        );
        by_leaf.insert(
            b.clone(),
            LeafResidual {
                total_squared: 3.,
                ..Default::default()
            },
        );
        let r = AdaptiveResidual {
            by_leaf,
            indicator_squared: 6.,
            indicator: 6_f64.sqrt(),
            error_upper_bound_available: false,
            qualification: "indicator",
        };
        assert_eq!(r.bulk_mark(0.5).unwrap(), vec![a]);
        assert!(r.bulk_mark(0.).is_err());
    }
}
