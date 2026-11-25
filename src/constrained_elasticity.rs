//! Matrix-free solution of elasticity with exact affine scalar constraints.
//!
//! The three components use independent [`Lift`] maps. The reduced system is
//! applied as C^T K C and its right hand side is C^T(f-Kg), without forming a
//! dense constraint matrix or a dense reduced stiffness matrix.

use crate::constraint_system::Lift;
use crate::elasticity::Stiffness;

/// A displacement satisfying the supplied affine constraints.
#[derive(Debug, Clone)]
pub struct ConstrainedSolution {
    /// Reconstructed interleaved displacement, in mm.
    pub displacement_mm: Vec<f64>,
    /// Reduced unknowns, concatenated in x, y, z component order.
    pub reduced_values: Vec<f64>,
    /// Norm of the reduced residual divided by the reduced right-hand-side norm.
    pub relative_residual: f64,
    /// Number of conjugate-gradient iterations.
    pub iterations: usize,
    /// Strain energy in N mm, computed as 0.5 u^T K u.
    pub strain_energy_n_mm: f64,
}

#[derive(Clone, Copy)]
struct Term {
    index: usize,
    coefficient: f64,
}

/// Solve `K u = f` in the affine space defined independently for each
/// displacement component by `component_lifts`.
///
/// The CSR stiffness uses interleaved xyz vectors. Each lift describes the
/// exact affine map from its free roots to all scalar nodes. The reduced
/// system is solved by diagonally preconditioned conjugate gradients.
pub fn solve_constrained(
    stiffness: &Stiffness,
    component_lifts: [&Lift; 3],
    force_n: &[f64],
    tolerance: f64,
    max_iterations: usize,
) -> Result<ConstrainedSolution, &'static str> {
    if !tolerance.is_finite() || tolerance <= 0. || tolerance >= 1. {
        return Err("tolerance must be finite and in (0,1)");
    }
    if max_iterations == 0 {
        return Err("iteration limit must be positive");
    }
    let n = stiffness
        .offsets
        .len()
        .checked_sub(1)
        .ok_or("invalid CSR offsets")?;
    if n == 0 || stiffness.size() != force_n.len() || force_n.iter().any(|x| !x.is_finite()) {
        return Err("invalid stiffness or force dimensions");
    }
    stiffness.validate()?;

    let mut component_maps: Vec<Vec<Vec<Term>>> = Vec::with_capacity(3);
    let mut component_offsets: Vec<Vec<f64>> = Vec::with_capacity(3);
    let mut reduced_size = 0usize;
    for lift in component_lifts {
        if lift.lines.len() != n {
            return Err("constraint lift has wrong number of rows");
        }
        let base = reduced_size;
        reduced_size = reduced_size
            .checked_add(lift.free_roots.len())
            .ok_or("reduced size overflow")?;
        let mut maps = Vec::with_capacity(n);
        let mut offsets = Vec::with_capacity(n);
        for (dof, line) in lift.lines.iter().enumerate() {
            if line.dof != dof {
                return Err("constraint lift rows are not in DOF order");
            }
            let mut terms = Vec::with_capacity(line.terms.len());
            let mut previous = None;
            for &(index, ratio) in &line.terms {
                if index >= lift.free_roots.len() || previous.is_some_and(|p| index <= p) {
                    return Err("invalid or unsorted constraint lift terms");
                }
                previous = Some(index);
                let coefficient = ratio.numerator() as f64 / ratio.denominator() as f64;
                if !coefficient.is_finite() {
                    return Err("nonfinite constraint coefficient");
                }
                if coefficient != 0. {
                    terms.push(Term {
                        index: base + index,
                        coefficient,
                    });
                }
            }
            let offset = line.offset.numerator() as f64 / line.offset.denominator() as f64;
            if !offset.is_finite() {
                return Err("nonfinite constraint offset");
            }
            maps.push(terms);
            offsets.push(offset);
        }
        component_maps.push(maps);
        component_offsets.push(offsets);
    }

    let mut full_offsets = vec![0.; 3 * n];
    for c in 0..3 {
        for i in 0..n {
            full_offsets[3 * i + c] = component_offsets[c][i];
        }
    }
    let kg = stiffness.apply(&full_offsets)?;
    let shifted: Vec<f64> = force_n.iter().zip(kg).map(|(f, k)| f - k).collect();
    if shifted.iter().any(|x| !x.is_finite()) {
        return Err("nonfinite shifted force");
    }

    let mut rhs = vec![0.; reduced_size];
    for i in 0..n {
        for c in 0..3 {
            for term in &component_maps[c][i] {
                rhs[term.index] += term.coefficient * shifted[3 * i + c];
            }
        }
    }
    if rhs.iter().any(|x| !x.is_finite()) {
        return Err("nonfinite reduced force");
    }

    // Exact diagonal of the congruence, accumulated from sparse CSR entries.
    let mut diagonal = vec![0.; reduced_size];
    for i in 0..n {
        for entry in stiffness.offsets[i]..stiffness.offsets[i + 1] {
            let j = stiffness.columns[entry];
            let block = stiffness.values[entry];
            for c in 0..3 {
                for d in 0..3 {
                    let a = &component_maps[c][i];
                    let b = &component_maps[d][j];
                    for ta in a {
                        for tb in b {
                            if ta.index == tb.index {
                                diagonal[ta.index] += ta.coefficient * block[c][d] * tb.coefficient;
                            }
                        }
                    }
                }
            }
        }
    }
    if diagonal.iter().any(|d| !d.is_finite() || *d <= 0.) {
        return Err("reduced stiffness has nonpositive or nonfinite diagonal");
    }

    let mut reduced = vec![0.; reduced_size];
    let mut iterations = 0usize;
    if reduced_size > 0 {
        let rhs_norm = norm(&rhs);
        let target = tolerance * rhs_norm;
        let mut residual = rhs.clone();
        let mut z: Vec<f64> = residual.iter().zip(&diagonal).map(|(r, d)| r / d).collect();
        let mut direction = z.clone();
        let mut rz = dot(&residual, &z);
        if !rz.is_finite() || rz < 0. {
            return Err("invalid preconditioned residual");
        }
        while norm(&residual) > target && (rhs_norm == 0. || norm(&residual) != 0.) {
            if iterations >= max_iterations {
                return Err("conjugate-gradient iteration limit reached");
            }
            let applied = apply_reduced(stiffness, &component_maps, n, &direction)?;
            let curvature = dot(&direction, &applied);
            if !curvature.is_finite() || curvature <= 0. {
                return Err("reduced stiffness is singular or indefinite");
            }
            let alpha = rz / curvature;
            if !alpha.is_finite() {
                return Err("nonfinite conjugate-gradient step");
            }
            axpy(&mut reduced, alpha, &direction);
            axpy(&mut residual, -alpha, &applied);
            iterations += 1;
            if norm(&residual) <= target {
                break;
            }
            z = residual.iter().zip(&diagonal).map(|(r, d)| r / d).collect();
            let next_rz = dot(&residual, &z);
            if !next_rz.is_finite() || next_rz <= 0. {
                return Err("invalid conjugate-gradient residual");
            }
            let beta = next_rz / rz;
            if !beta.is_finite() {
                return Err("nonfinite conjugate-gradient recurrence");
            }
            for (p, &zi) in direction.iter_mut().zip(&z) {
                *p = zi + beta * *p;
            }
            rz = next_rz;
        }
    }
    if reduced.iter().any(|x| !x.is_finite()) {
        return Err("nonfinite reduced solution");
    }

    // Recompute the reduced residual to avoid relying on recurrence drift.
    let mut relative_residual = 0.;
    if reduced_size > 0 {
        let applied = apply_reduced(stiffness, &component_maps, n, &reduced)?;
        let res: Vec<f64> = rhs.iter().zip(applied).map(|(b, a)| b - a).collect();
        let denominator = norm(&rhs);
        relative_residual = if denominator == 0. {
            norm(&res)
        } else {
            norm(&res) / denominator
        };
        if !relative_residual.is_finite()
            || relative_residual > tolerance * (1. + 32. * f64::EPSILON)
        {
            return Err("reduced residual exceeds tolerance");
        }
    }

    let displacement = reconstruct(&component_maps, &component_offsets, &reduced, n);
    if displacement.iter().any(|x| !x.is_finite()) {
        return Err("nonfinite reconstructed displacement");
    }
    let kd = stiffness.apply(&displacement)?;
    let quadratic = dot(&displacement, &kd);
    if !quadratic.is_finite() || quadratic < -64. * f64::EPSILON * norm(&displacement) * norm(&kd) {
        return Err("invalid nonpositive strain energy");
    }
    Ok(ConstrainedSolution {
        displacement_mm: displacement,
        reduced_values: reduced,
        relative_residual,
        iterations,
        strain_energy_n_mm: 0.5 * quadratic.max(0.),
    })
}

fn reconstruct(
    maps: &[Vec<Vec<Term>>],
    offsets: &[Vec<f64>],
    reduced: &[f64],
    n: usize,
) -> Vec<f64> {
    let mut full = vec![0.; 3 * n];
    for i in 0..n {
        for c in 0..3 {
            full[3 * i + c] = offsets[c][i]
                + maps[c][i]
                    .iter()
                    .map(|t| t.coefficient * reduced[t.index])
                    .sum::<f64>();
        }
    }
    full
}

fn apply_reduced(
    stiffness: &Stiffness,
    maps: &[Vec<Vec<Term>>],
    n: usize,
    x: &[f64],
) -> Result<Vec<f64>, &'static str> {
    let mut full = vec![0.; 3 * n];
    for i in 0..n {
        for c in 0..3 {
            full[3 * i + c] = maps[c][i].iter().map(|t| t.coefficient * x[t.index]).sum();
        }
    }
    let product = stiffness.apply(&full)?;
    let mut result = vec![0.; x.len()];
    for i in 0..n {
        for c in 0..3 {
            for term in &maps[c][i] {
                result[term.index] += term.coefficient * product[3 * i + c];
            }
        }
    }
    if result.iter().any(|v| !v.is_finite()) {
        return Err("nonfinite reduced matrix product");
    }
    Ok(result)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}
fn axpy(y: &mut [f64], alpha: f64, x: &[f64]) {
    for (yi, xi) in y.iter_mut().zip(x) {
        *yi += alpha * xi;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraint_system::{LiftLine, Ratio};
    use crate::elasticity::Matrix;

    fn ratio(n: i128, d: i128) -> Ratio {
        Ratio::new(n, d).unwrap()
    }
    fn lift(rows: Vec<(Vec<(usize, Ratio)>, Ratio)>, roots: Vec<usize>) -> Lift {
        Lift {
            free_roots: roots,
            lines: rows
                .into_iter()
                .enumerate()
                .map(|(dof, (terms, offset))| LiftLine { dof, terms, offset })
                .collect(),
        }
    }
    fn identity_stiffness(n: usize) -> Stiffness {
        let eye: Matrix = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        Stiffness {
            offsets: (0..=n).collect(),
            columns: (0..n).collect(),
            values: vec![eye; n],
            integrated_volume: 1.,
        }
    }

    #[test]
    fn solves_componentwise_fractional_affine_lifts_and_offsets() {
        let half = ratio(1, 2);
        let three_quarters = ratio(3, 4);
        let x = lift(
            vec![
                (vec![(0, ratio(1, 1))], ratio(0, 1)),
                (vec![(0, half)], half),
                (vec![], ratio(1, 1)),
            ],
            vec![0],
        );
        let y = lift(
            vec![
                (vec![(0, three_quarters)], ratio(1, 2)),
                (vec![], ratio(2, 1)),
                (vec![(0, ratio(1, 1))], ratio(0, 1)),
            ],
            vec![2],
        );
        // The z component is entirely prescribed at nonzero offsets.
        let z = lift(
            vec![
                (vec![], ratio(1, 1)),
                (vec![], ratio(2, 1)),
                (vec![], ratio(3, 1)),
            ],
            vec![],
        );
        let f = vec![1., 1., 0., 2., 2., 0., 3., 4., 0.];
        let result =
            solve_constrained(&identity_stiffness(3), [&x, &y, &z], &f, 1e-12, 20).unwrap();
        let u = result.displacement_mm;
        assert!((u[0] - 1.4).abs() < 1e-12);
        assert!((u[3] - 1.2).abs() < 1e-12);
        assert!((u[6] - 1.).abs() < 1e-12);
        assert!((u[1] - 2.6).abs() < 1e-12);
        assert!((u[4] - 2.).abs() < 1e-12);
        assert!((u[7] - 2.8).abs() < 1e-12);
        assert_eq!([u[2], u[5], u[8]], [1., 2., 3.]);
        assert!(result.relative_residual <= 1e-12);
    }

    #[test]
    fn rejects_singular_excited_reduced_system() {
        let zero: Matrix = [[0.; 3]; 3];
        let k = Stiffness {
            offsets: vec![0, 1],
            columns: vec![0],
            values: vec![zero],
            integrated_volume: 1.,
        };
        let free = lift(vec![(vec![(0, ratio(1, 1))], ratio(0, 1))], vec![0]);
        let fixed = lift(vec![(vec![], ratio(0, 1))], vec![]);
        assert!(solve_constrained(&k, [&free, &fixed, &fixed], &[1., 0., 0.], 1e-10, 10).is_err());
    }

    #[test]
    fn refuses_invalid_force_and_accepts_fully_prescribed_state() {
        let k = identity_stiffness(1);
        let fixed = lift(vec![(vec![], ratio(2, 1))], vec![]);
        let solved = solve_constrained(&k, [&fixed, &fixed, &fixed], &[0.; 3], 1e-12, 1).unwrap();
        assert_eq!(solved.displacement_mm, vec![2.; 3]);
        assert_eq!(solved.iterations, 0);
        assert!(
            solve_constrained(&k, [&fixed, &fixed, &fixed], &[f64::NAN, 0., 0.], 1e-12, 1).is_err()
        );
    }
}
