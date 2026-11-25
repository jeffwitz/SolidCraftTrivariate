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
    solve_with_guess(
        stiffness,
        component_lifts,
        force_n,
        tolerance,
        max_iterations,
        None,
        None,
    )
}

/// Solve from a finite initial vector in the same reduced, component-major
/// ordering as `ConstrainedSolution::reduced_values`. Essential values are
/// applied through the supplied lifts independently of the initial guess.
pub fn solve_constrained_from(
    stiffness: &Stiffness,
    component_lifts: [&Lift; 3],
    force_n: &[f64],
    tolerance: f64,
    max_iterations: usize,
    initial_reduced: &[f64],
) -> Result<ConstrainedSolution, &'static str> {
    solve_with_guess(
        stiffness,
        component_lifts,
        force_n,
        tolerance,
        max_iterations,
        Some(initial_reduced),
        None,
    )
}

/// A fixed symmetric positive definite operator on the mounted reduced space.
pub trait ReducedPreconditioner {
    fn apply(&self, residual: &[f64]) -> Result<Vec<f64>, &'static str>;
}

/// Solve with a supplied reduced preconditioner in the same component-major
/// ordering as the returned reduced values. Convergence uses the full measured
/// force residual projected back through the exact affine constraints.
pub fn solve_preconditioned(
    stiffness: &Stiffness,
    component_lifts: [&Lift; 3],
    force_n: &[f64],
    tolerance: f64,
    max_iterations: usize,
    preconditioner: &dyn ReducedPreconditioner,
) -> Result<ConstrainedSolution, &'static str> {
    solve_with_guess(
        stiffness,
        component_lifts,
        force_n,
        tolerance,
        max_iterations,
        None,
        Some(preconditioner),
    )
}

fn solve_with_guess(
    stiffness: &Stiffness,
    component_lifts: [&Lift; 3],
    force_n: &[f64],
    tolerance: f64,
    max_iterations: usize,
    guess: Option<&[f64]>,
    preconditioner: Option<&dyn ReducedPreconditioner>,
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
    let shifted = stiffness.equilibrium_residual(&full_offsets, force_n)?;
    let rhs = project_equilibrium(&component_maps, &shifted, reduced_size)?;

    let diagonal = if preconditioner.is_none() {
        let mut diagonal = vec![0.; reduced_size];
        for i in 0..n {
            for entry in stiffness.offsets[i]..stiffness.offsets[i + 1] {
                let j = stiffness.columns[entry];
                let block = stiffness.values[entry];
                for (c, maps) in component_maps.iter().enumerate() {
                    for left in &maps[i] {
                        for right in &maps[j] {
                            if left.index == right.index {
                                diagonal[left.index] +=
                                    left.coefficient * block[c][c] * right.coefficient;
                            }
                        }
                    }
                }
            }
        }
        if diagonal.iter().any(|d| !d.is_finite() || *d <= 0.) {
            return Err("reduced stiffness has nonpositive or nonfinite diagonal");
        }
        diagonal
    } else {
        Vec::new()
    };
    let precondition = |r: &[f64]| -> Result<Vec<f64>, &'static str> {
        let z = match preconditioner {
            Some(prec) => prec.apply(r)?,
            None => r.iter().zip(&diagonal).map(|(x, d)| x / d).collect(),
        };
        if z.len() != r.len() || z.iter().any(|x| !x.is_finite()) {
            return Err("invalid reduced preconditioner output");
        }
        Ok(z)
    };

    if guess.is_some_and(|v| v.len() != reduced_size || v.iter().any(|x| !x.is_finite())) {
        return Err("invalid initial reduced vector");
    }
    let mut reduced = guess.map_or_else(|| vec![0.; reduced_size], |v| v.to_vec());
    let mut iterations = 0usize;
    let mut relative_residual = 0.;
    if reduced_size > 0 {
        let rhs_norm = norm(&rhs);
        if !rhs_norm.is_finite() {
            return Err("nonfinite reduced load norm");
        }
        let target = tolerance * rhs_norm;
        let equilibrium = |reduced: &[f64]| {
            let full = reconstruct(&component_maps, &component_offsets, reduced, n);
            let r = stiffness.equilibrium_residual(&full, force_n)?;
            project_equilibrium(&component_maps, &r, reduced_size)
        };
        let mut residual = if rhs_norm == 0. {
            reduced.fill(0.);
            vec![0.; reduced_size]
        } else if guess.is_some() {
            equilibrium(&reduced)?
        } else {
            rhs.clone()
        };
        if rhs_norm > 0. {
            relative_residual = norm(&residual) / rhs_norm;
        }
        let mut z: Vec<f64> = precondition(&residual)?;
        let mut direction = z.clone();
        let mut rz = dot(&residual, &z);
        if !rz.is_finite() || rz < 0. {
            return Err("invalid preconditioned residual");
        }
        while norm(&residual) > target {
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
                residual = equilibrium(&reduced)?;
                relative_residual = norm(&residual) / rhs_norm;
                if !relative_residual.is_finite() {
                    return Err("nonfinite measured reduced residual");
                }
                if relative_residual <= tolerance {
                    break;
                }
                z = precondition(&residual)?;
                direction.clone_from(&z);
                rz = dot(&residual, &z);
                if !rz.is_finite() || rz <= 0. {
                    return Err("nonpositive preconditioned equilibrium residual");
                }
                continue;
            }
            z = precondition(&residual)?;
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

/// Project a measured full force residual, retaining low-order terms where
/// neighboring hanging-constraint forces cancel on a common master.
fn project_equilibrium(
    maps: &[Vec<Vec<Term>>],
    full: &[f64],
    size: usize,
) -> Result<Vec<f64>, &'static str> {
    let mut sum = vec![0_f64; size];
    let mut correction = vec![0.; size];
    for (i, vector) in full.as_chunks::<3>().0.iter().enumerate() {
        for (component, &force) in vector.iter().enumerate() {
            for term in &maps[component][i] {
                let j = term.index;
                let product = term.coefficient * force;
                let next = sum[j] + product;
                let roundoff = if sum[j].abs() >= product.abs() {
                    (sum[j] - next) + product
                } else {
                    (product - next) + sum[j]
                };
                correction[j] += roundoff + term.coefficient.mul_add(force, -product);
                sum[j] = next;
            }
        }
    }
    for (s, c) in sum.iter_mut().zip(correction) {
        *s += c;
        if !s.is_finite() {
            return Err("nonfinite projected force residual");
        }
    }
    Ok(sum)
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
    fn initialized_reduced_solution_is_checked_against_the_full_affine_equilibrium() {
        let k = identity_stiffness(2);
        let x = lift(
            vec![
                (vec![(0, ratio(1, 1))], ratio(0, 1)),
                (vec![(0, ratio(1, 2))], ratio(3, 2)),
            ],
            vec![0],
        );
        let fixed = lift(vec![(vec![], ratio(0, 1)), (vec![], ratio(0, 1))], vec![]);
        let force = [2., 0., 0., 1., 0., 0.];
        let cold = solve_constrained(&k, [&x, &fixed, &fixed], &force, 1e-12, 10).unwrap();
        let warm = solve_constrained_from(
            &k,
            [&x, &fixed, &fixed],
            &force,
            1e-12,
            10,
            &cold.reduced_values,
        )
        .unwrap();
        assert_eq!(warm.iterations, 0);
        assert_eq!(cold.displacement_mm, warm.displacement_mm);
        assert!(warm.relative_residual <= 1e-12);
        assert!(solve_constrained_from(&k, [&x, &fixed, &fixed], &force, 1e-12, 10, &[]).is_err());
        assert!(
            solve_constrained_from(&k, [&x, &fixed, &fixed], &force, 1e-12, 10, &[f64::NAN])
                .is_err()
        );
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
