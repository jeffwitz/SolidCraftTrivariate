//! Cold and edited row mechanics with an affine manufactured displacement.
use solidcraft_trivariate::{
    constrained_elasticity::solve_constrained,
    constraint_system::{resolve, Lift, Ratio},
    elasticity::{basis_geometry, gauss},
    json::{object, Json},
    row_session::{RowAnalysisSession, RowRequest, RowState},
};
use std::collections::BTreeMap;
use std::time::Instant;

const A: [[i128; 3]; 3] = [[3, 1, -2], [-1, 2, 1], [2, 0, -1]];
const T: [i128; 3] = [7, -4, 3];
fn n(x: f64) -> Json {
    Json::Number(x.to_string())
}
fn stress(e: f64, nu: f64) -> ([[f64; 3]; 3], f64) {
    let mu = e / (2. * (1. + nu));
    let lam = e * nu / ((1. + nu) * (1. - 2. * nu));
    let eps: [[f64; 3]; 3] =
        std::array::from_fn(|i| std::array::from_fn(|j| (A[i][j] + A[j][i]) as f64 / 2000.));
    let tr = eps[0][0] + eps[1][1] + eps[2][2];
    let sigma: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| 2. * mu * eps[i][j] + if i == j { lam * tr } else { 0. })
    });
    (
        sigma,
        0.5 * sigma
            .iter()
            .flatten()
            .zip(eps.iter().flatten())
            .map(|(x, y)| x * y)
            .sum::<f64>(),
    )
}
fn physical_coefficients(state: &RowState) -> Result<Vec<[Ratio; 3]>, &'static str> {
    let mut field = vec![None; state.space.ndof];
    for b in &state.space.blocks {
        for (&i, h) in b.dofs.iter().zip(&b.controls) {
            let mut values = [Ratio::zero(); 3];
            for d in 0..3 {
                let t = T[d] * h[3] + (0..3).map(|j| A[d][j] * h[j]).sum::<i128>();
                values[d] = Ratio::new(t, 1000 * h[3])?;
            }
            if let Some(previous) = field[i] {
                if previous != values {
                    return Err("inconsistent common field control");
                }
            }
            field[i] = Some(values);
        }
    }
    field
        .into_iter()
        .map(|o| o.ok_or("unassigned field control"))
        .collect()
}
fn mechanics(state: &RowState) -> Result<Json, &'static str> {
    let t = Instant::now();
    let field = physical_coefficients(state)?;
    let r = &state.request;
    let (sigma, energy_density) = stress(r.young_mpa, r.poisson);
    let points = gauss(r.gauss_order)?;
    let blocks: BTreeMap<_, _> = state.space.blocks.iter().map(|b| (&b.id, b)).collect();
    let mut rhs = vec![0.; 3 * state.space.ndof];
    for (side, feature) in &state.certificate.design_faces {
        if feature == "left" {
            continue;
        }
        let b = blocks[&side.block];
        let tangent = match side.axis {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let orientation = if side.upper ^ (side.axis == 1) {
            1.
        } else {
            -1.
        };
        for &(u, wu) in &points {
            for &(v, wv) in &points {
                let mut uvw = [0.; 3];
                uvw[tangent[0]] = u;
                uvw[tangent[1]] = v;
                uvw[side.axis] = f64::from(u8::from(side.upper));
                let sample = basis_geometry(b, uvw)?;
                let [s, t] = tangent.map(|d| sample.jacobian_columns[d]);
                let area: [f64; 3] = std::array::from_fn(|i| {
                    orientation
                        * wu
                        * wv
                        * (s[(i + 1) % 3] * t[(i + 2) % 3] - s[(i + 2) % 3] * t[(i + 1) % 3])
                });
                let traction =
                    sigma.map(|row| row.iter().zip(area).map(|(x, y)| x * y).sum::<f64>());
                for (&id, value) in b.dofs.iter().zip(sample.values) {
                    for j in 0..3 {
                        rhs[3 * id + j] += traction[j] * value;
                    }
                }
            }
        }
    }
    let prescribed = &state.space.boundary_dofs["left"];
    let component = |c| {
        let bc = prescribed
            .iter()
            .map(|&d| (d, field[d][c]))
            .collect::<Vec<_>>();
        resolve(state.space.ndof, &[], &bc)
    };
    let lifts: [Lift; 3] = [component(0)?, component(1)?, component(2)?];
    let boundary_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let result = solve_constrained(
        &state.stiffness,
        [&lifts[0], &lifts[1], &lifts[2]],
        &rhs,
        1e-11,
        12000,
    )?;
    let solve_s = t.elapsed().as_secs_f64();
    let mut max: f64 = 0.;
    let mut err: f64 = 0.;
    for (analytic, actual) in field.iter().flatten().zip(&result.displacement_mm) {
        let expected = analytic.numerator() as f64 / analytic.denominator() as f64;
        max = max.max(expected.abs());
        err = err.max((actual - expected).abs());
    }
    let ideal_area = (r.half_widths[0] * r.half_widths[1]) as f64
        - std::f64::consts::PI
            * r.radii
                .iter()
                .map(|v| (f64::from(*v) / 2.).powi(2))
                .sum::<f64>();
    let energy_error = (result.strain_energy_n_mm / (8. * ideal_area * energy_density) - 1.).abs();
    if r.gauss_order == 8 && (err / max > 1e-7 || energy_error > 1e-7) {
        return Err("order-eight physical patch tolerance missed");
    }
    let update = state.update;
    Ok(object([
        ("integrated_blocks", n(update.integrated as f64)),
        ("reused_blocks", n(update.reused as f64)),
        ("build_s", n(update.build_s)),
        ("audit_numbering_s", n(update.audit_numbering_s)),
        ("element_s", n(update.element_s)),
        ("merge_csr_s", n(update.merge_s)),
        ("boundary_and_lift_s", n(boundary_s)),
        ("solve_s", n(solve_s)),
        ("solver_iterations", n(result.iterations as f64)),
        ("relative_residual", n(result.relative_residual)),
        ("max_coefficient_error_mm", n(err)),
        ("relative_displacement_error", n(err / max)),
        ("relative_energy_error", n(energy_error)),
        ("energy_n_mm", n(result.strain_energy_n_mm)),
        ("scalar_dofs", n(state.space.ndof as f64)),
    ]))
}
fn main() -> Result<(), &'static str> {
    let mut cases = Vec::new();
    for count in [1, 2, 4, 8, 16] {
        for order in [4, 6, 8] {
            let mut session = RowAnalysisSession::default();
            let mut r = RowRequest {
                half_widths: [10 * count, 20],
                centers: (0..count).map(|i| [20 * i + 10 - 10 * count, 0]).collect(),
                radii: vec![4; count as usize],
                young_mpa: 210000.,
                poisson: 0.3,
                gauss_order: order,
            };
            let initial = mechanics(session.update(r.clone())?)?;
            r.radii[count as usize / 2] += 1;
            let edit = mechanics(session.update(r.clone())?)?;
            let mut cold = RowAnalysisSession::default();
            let fresh = mechanics(cold.update(r)?)?;
            for key in [
                "energy_n_mm",
                "relative_residual",
                "relative_displacement_error",
            ] {
                if edit.field(key)? != fresh.field(key)? {
                    return Err("cold replay mismatch");
                }
            }
            cases.push(object([
                ("bores", n(count as f64)),
                ("order", n(order as f64)),
                ("initial", initial),
                ("edit", edit),
                ("cold_replay", fresh),
            ]));
        }
    }
    println!("{}", Json::Array(cases));
    Ok(())
}
