//! Integrated moving rational trace, geometry stiffness, and adjoint example.

use solidcraft_trivariate::analysis_space::{AnalysisBlock, AnalysisSpace};
use solidcraft_trivariate::constrained_elasticity::solve_constrained;
use solidcraft_trivariate::constraint_system::{resolve, Lift, Ratio};
use solidcraft_trivariate::elasticity::assemble;
use solidcraft_trivariate::elasticity_sensitivity::{assemble_sensitivity, MaterialJet};
use solidcraft_trivariate::face_constraints::{nested_constraints, FaceNet, Row, Span};
use solidcraft_trivariate::incidence::Orientation;
use solidcraft_trivariate::mechanical_sensitivity::{
    differentiate, LiftRate, MechanicalDerivative, SensitivityResult,
};
use solidcraft_trivariate::trace_sensitivity::nested_weight_jets;

const HEIGHT: i64 = 2;
const NU: f64 = 0.25;
const FORCE_DIRECTION: [f64; 3] = [0.0, 1.0, 0.0];
const G_DERIVATIVE: [f64; 3] = [1.0, 2.0, -1.0];
const G_COEFFICIENT: [i128; 3] = [1, 2, -1];

#[derive(Clone)]
struct TracePair {
    slave_dofs: Vec<usize>,
    master_dofs: Vec<usize>,
    master_net: FaceNet,
    master_dw: Vec<f64>,
    slave_net: FaceNet,
    slave_dw: Vec<f64>,
    span: Span,
}

struct MovingMesh {
    space: AnalysisSpace,
    jets: Vec<Vec<[f64; 4]>>,
    equations: Vec<(usize, Row)>,
    traces: Vec<TracePair>,
    outer_dofs: Vec<usize>,
    load_dof: usize,
    t: f64,
    denominator: i128,
}

type Curve = [[i128; 3]; 3];
type CurveJet = [[f64; 3]; 3];

fn split_curve(
    source: [[i128; 3]; 3],
    derivative: [[f64; 3]; 3],
) -> (Curve, Curve, CurveJet, CurveJet) {
    let mut lower = [[0; 3]; 3];
    let mut upper = [[0; 3]; 3];
    let mut lower_d = [[0.0; 3]; 3];
    let mut upper_d = [[0.0; 3]; 3];
    for c in 0..3 {
        let (a, b, d) = (source[0][c], source[1][c], source[2][c]);
        lower[0][c] = 4 * a;
        lower[1][c] = 2 * (a + b);
        lower[2][c] = a + 2 * b + d;
        upper[0][c] = lower[2][c];
        upper[1][c] = 2 * (b + d);
        upper[2][c] = 4 * d;
        let (a, b, d) = (derivative[0][c], derivative[1][c], derivative[2][c]);
        lower_d[0][c] = 4.0 * a;
        lower_d[1][c] = 2.0 * (a + b);
        lower_d[2][c] = a + 2.0 * b + d;
        upper_d[0][c] = lower_d[2][c];
        upper_d[1][c] = 2.0 * (b + d);
        upper_d[2][c] = 4.0 * d;
    }
    (lower, upper, lower_d, upper_d)
}

fn curve(t_num: i128, den: i128) -> ([[i128; 3]; 3], [[f64; 3]; 3]) {
    let den2 = den * den;
    let t = t_num as f64 / den as f64;
    let source = [
        [den2, 0, den2],
        [den2, t_num * den, den2],
        [den2 - t_num * t_num, 2 * t_num * den, den2 + t_num * t_num],
    ];
    let derivative = [
        [0.0, 0.0, 0.0],
        [0.0, den2 as f64, 0.0],
        [
            -2.0 * t * den2 as f64,
            2.0 * den2 as f64,
            2.0 * t * den2 as f64,
        ],
    ];
    (source, derivative)
}

fn build_block(
    id: &str,
    arc: [[i128; 3]; 3],
    arc_jet: [[f64; 3]; 3],
    radial_half: [i64; 2],
    start_dof: usize,
    shared_u0: Option<[usize; 9]>,
) -> (AnalysisBlock, Vec<[f64; 4]>, Vec<usize>, usize) {
    let mut linear_controls = Vec::with_capacity(12);
    let mut linear_jets = Vec::with_capacity(12);
    for i in 0..3 {
        for &radial_tick in &radial_half {
            for k in 0..2 {
                let [x, y, weight] = arc[i];
                let [dx, dy, dw] = arc_jet[i];
                let radius = i128::from(radial_tick);
                let z = if k == 0 { 0_i128 } else { i128::from(HEIGHT) };
                linear_controls.push([radius * x, radius * y, 2 * z * weight, 2 * weight]);
                linear_jets.push([
                    radius as f64 * dx,
                    radius as f64 * dy,
                    (2 * z) as f64 * dw,
                    2.0 * dw,
                ]);
            }
        }
    }
    let (controls, jets) = elevate_linear_axes(&linear_controls, &linear_jets);
    let mut dofs = Vec::with_capacity(27);
    let mut next = start_dof;
    let mut outer_face = Vec::with_capacity(9);
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let dof = if i == 0 {
                    shared_u0.map(|face| face[j * 3 + k]).unwrap_or_else(|| {
                        let value = next;
                        next += 1;
                        value
                    })
                } else {
                    let value = next;
                    next += 1;
                    value
                };
                dofs.push(dof);
                if j == 0 {
                    outer_face.push(dof);
                }
            }
        }
    }
    (
        AnalysisBlock {
            id: id.into(),
            controls,
            dofs,
        },
        jets,
        outer_face,
        next,
    )
}

fn face_net(
    block: &AnalysisBlock,
    jet: &[[f64; 4]],
    radial_index: usize,
) -> (FaceNet, Vec<f64>, Vec<usize>) {
    let mut controls = Vec::with_capacity(9);
    let mut dw = Vec::with_capacity(9);
    let mut dofs = Vec::with_capacity(9);
    for i in 0..3 {
        for k in 0..3 {
            let local = (i * 3 + radial_index) * 3 + k;
            controls.push(
                block.controls[local]
                    .map(|x| i64::try_from(x).expect("dyadic arc controls fit i64")),
            );
            dw.push(jet[local][3]);
            dofs.push(block.dofs[local]);
        }
    }
    (
        FaceNet {
            degrees: [2, 2],
            controls,
        },
        dw,
        dofs,
    )
}

fn elevation_weights(index: usize) -> &'static [(usize, i128)] {
    match index {
        0 => &[(0, 2)],
        1 => &[(0, 1), (1, 1)],
        _ => &[(1, 2)],
    }
}

fn elevate_linear_axes(
    controls: &[[i128; 4]],
    jets: &[[f64; 4]],
) -> (Vec<[i128; 4]>, Vec<[f64; 4]>) {
    let mut elevated = Vec::with_capacity(27);
    let mut elevated_jets = Vec::with_capacity(27);
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let mut point = [0_i128; 4];
                let mut jet = [0.0; 4];
                for &(old_j, v_factor) in elevation_weights(j) {
                    for &(old_k, w_factor) in elevation_weights(k) {
                        let source = (i * 2 + old_j) * 2 + old_k;
                        let factor = v_factor * w_factor;
                        for c in 0..4 {
                            point[c] += factor * controls[source][c];
                            jet[c] += factor as f64 * jets[source][c];
                        }
                    }
                }
                elevated.push(point);
                elevated_jets.push(jet);
            }
        }
    }
    (elevated, elevated_jets)
}

fn build_mesh(t_num: i128, den: i128) -> Result<MovingMesh, &'static str> {
    let (full, full_jet) = curve(t_num, den);
    let (low, high, low_jet, high_jet) = split_curve(full, full_jet);
    let (external, external_jet, outer_dofs, next) =
        build_block("external", full, full_jet, [4, 3], 0, None);
    let (inner_low, low_block_jet, _, next) =
        build_block("inner_low", low, low_jet, [3, 2], next, None);
    // Conforming internal face only: low.u1 and high.u0 share the same nine
    // radial/height Bernstein coefficients. The T1 master face stays separate.
    let shared = std::array::from_fn(|slot| inner_low.dofs[18 + slot]);
    let (inner_high, high_block_jet, _, next) =
        build_block("inner_high", high, high_jet, [3, 2], next, Some(shared));
    let load_dof = inner_low.dofs[19];
    let blocks = vec![external, inner_low, inner_high];
    let jets = vec![external_jet, low_block_jet, high_block_jet];
    let (master_net, master_dw, master_dofs) = face_net(&blocks[0], &jets[0], 2);
    let mut equations = Vec::new();
    let mut traces = Vec::new();
    for (block_index, span) in [(1usize, Span::Lower), (2usize, Span::Upper)] {
        let (slave_net, slave_dw, slave_dofs) =
            face_net(&blocks[block_index], &jets[block_index], 0);
        let rows = nested_constraints(
            &master_net,
            &slave_net,
            [span, Span::Full],
            Orientation {
                swap: false,
                flip_first: false,
                flip_second: false,
            },
        )?;
        for (row_index, row) in rows.into_iter().enumerate() {
            equations.push((
                slave_dofs[row_index],
                Row {
                    terms: row
                        .terms
                        .into_iter()
                        .map(|(m, weight)| (master_dofs[m], weight))
                        .collect(),
                    denominator: row.denominator,
                },
            ));
        }
        traces.push(TracePair {
            slave_dofs,
            master_dofs: master_dofs.clone(),
            master_net: master_net.clone(),
            master_dw: master_dw.clone(),
            slave_net,
            slave_dw,
            span,
        });
    }
    Ok(MovingMesh {
        space: AnalysisSpace {
            ndof: next,
            blocks,
            boundary_dofs: std::collections::BTreeMap::from([(
                "outer_r2".into(),
                outer_dofs.clone(),
            )]),
        },
        jets,
        equations,
        traces,
        outer_dofs,
        load_dof,
        t: t_num as f64 / den as f64,
        denominator: den,
    })
}

fn boundary_lifts(mesh: &MovingMesh, t_num: i128) -> Result<[Lift; 3], &'static str> {
    let prescribed: [Vec<_>; 3] = std::array::from_fn(|c| {
        mesh.outer_dofs
            .iter()
            .map(|&dof| {
                (
                    dof,
                    Ratio::new(G_COEFFICIENT[c] * t_num, mesh.denominator).unwrap(),
                )
            })
            .collect()
    });
    Ok([
        resolve(mesh.space.ndof, &mesh.equations, &prescribed[0])?,
        resolve(mesh.space.ndof, &mesh.equations, &prescribed[1])?,
        resolve(mesh.space.ndof, &mesh.equations, &prescribed[2])?,
    ])
}

fn lift_rates(mesh: &MovingMesh, lifts: [&Lift; 3]) -> Result<[Vec<LiftRate>; 3], &'static str> {
    let mut result: [Vec<LiftRate>; 3] =
        std::array::from_fn(|_| vec![LiftRate::default(); mesh.space.ndof]);
    for c in 0..3 {
        for &dof in &mesh.outer_dofs {
            result[c][dof].offset = G_DERIVATIVE[c];
        }
    }
    for trace in &mesh.traces {
        let rows = nested_weight_jets(
            &trace.master_net,
            &trace.master_dw,
            &trace.slave_net,
            &trace.slave_dw,
            [trace.span, Span::Full],
            Orientation {
                swap: false,
                flip_first: false,
                flip_second: false,
            },
        )?;
        for (row_index, row) in rows.iter().enumerate() {
            let dependent = trace.slave_dofs[row_index];
            let mut by_root = std::collections::BTreeMap::<usize, f64>::new();
            for coefficient in row {
                let master_dof = trace.master_dofs[coefficient.master];
                let master_line = &lifts[0].lines[master_dof];
                if !master_line.offset.eq(&Ratio::zero()) {
                    return Err("trace master unexpectedly carries a prescribed offset");
                }
                for &(root, ratio) in &master_line.terms {
                    *by_root.entry(root).or_default() += coefficient.derivative
                        * ratio.numerator() as f64
                        / ratio.denominator() as f64;
                }
            }
            for component in &mut result {
                let terms: Vec<_> = by_root
                    .iter()
                    .filter_map(|(&root, &rate)| (rate != 0.0).then_some((root, rate)))
                    .collect();
                let target = &mut component[dependent].terms;
                if target.is_empty() {
                    *target = terms;
                } else {
                    assert_eq!(target.len(), terms.len());
                    for ((left_root, left), (right_root, right)) in target.iter_mut().zip(terms) {
                        assert_eq!(*left_root, right_root);
                        assert!((*left - right).abs() <= 2e-13);
                        *left = 0.5 * (*left + right);
                    }
                }
            }
        }
    }
    Ok(result)
}

fn load(mesh: &MovingMesh) -> Vec<f64> {
    let mut force = vec![0.0; 3 * mesh.space.ndof];
    force[3 * mesh.load_dof + 1] = FORCE_DIRECTION[1];
    force
}

fn objective(mesh: &MovingMesh, order: usize) -> Result<(f64, Vec<f64>), &'static str> {
    let lifts = boundary_lifts(mesh, (mesh.t * mesh.denominator as f64).round() as i128)?;
    let force = load(mesh);
    let stiffness = assemble(&mesh.space, 10.0 + mesh.t, NU, order)?;
    let solution = solve_constrained(&stiffness, lifts.each_ref(), &force, 1e-13, 5000)?;
    Ok((
        solution.displacement_mm[3 * mesh.load_dof + 1],
        solution.displacement_mm,
    ))
}

fn bernstein2(x: f64) -> [f64; 3] {
    [(1.0 - x).powi(2), 2.0 * x * (1.0 - x), x.powi(2)]
}

fn evaluate_q2(block: &AnalysisBlock, [u, v, w]: [f64; 3]) -> [f64; 3] {
    let basis = [bernstein2(u), bernstein2(v), bernstein2(w)];
    let mut homogeneous = [0.0; 4];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let factor = basis[0][i] * basis[1][j] * basis[2][k];
                let control = block.controls[(i * 3 + j) * 3 + k];
                for c in 0..4 {
                    homogeneous[c] += factor * control[c] as f64;
                }
            }
        }
    }
    std::array::from_fn(|c| homogeneous[c] / homogeneous[3])
}

fn verify_physical_maps(mesh: &MovingMesh) {
    for (index, block) in mesh.space.blocks.iter().enumerate() {
        let (global_start, global_width) = match index {
            0 => (0.0, 1.0),
            1 => (0.0, 0.5),
            _ => (0.5, 0.5),
        };
        for (u, v, w) in [(0.173, 0.31, 0.61), (0.71, 0.83, 0.27)] {
            let actual = evaluate_q2(block, [u, v, w]);
            let global_u = global_start + global_width * u;
            let angle = 2.0 * (mesh.t * global_u).atan();
            let radius = if index == 0 {
                2.0 - 0.5 * v
            } else {
                1.5 - 0.5 * v
            };
            let expected = [radius * angle.cos(), radius * angle.sin(), 2.0 * w];
            for c in 0..3 {
                assert!(
                    (actual[c] - expected[c]).abs() <= 2e-13,
                    "block={} component={c} error={}",
                    block.id,
                    (actual[c] - expected[c]).abs()
                );
            }
        }
    }
}

fn derivative_result(
    mesh: &MovingMesh,
    order: usize,
    omit_geometry: bool,
    omit_trace: bool,
    omit_lifting: bool,
    omit_material: bool,
) -> Result<SensitivityResult, &'static str> {
    let lifts = boundary_lifts(mesh, 1)?;
    let rates = lift_rates(mesh, lifts.each_ref())?;
    let rates: [Vec<LiftRate>; 3] = std::array::from_fn(|c| {
        rates[c]
            .iter()
            .map(|row| LiftRate {
                terms: if omit_trace {
                    Vec::new()
                } else {
                    row.terms.clone()
                },
                offset: if omit_lifting { 0.0 } else { row.offset },
            })
            .collect()
    });
    let control_jets = if omit_geometry {
        vec![vec![[0.0; 4]; 27]; mesh.space.blocks.len()]
    } else {
        mesh.jets.clone()
    };
    let sensitivity = assemble_sensitivity(
        &mesh.space,
        &control_jets,
        MaterialJet {
            young_mpa: 11.0,
            poisson: NU,
            d_young_mpa: if omit_material { 0.0 } else { 1.0 },
            d_poisson: 0.0,
        },
        order,
    )?;
    let force = load(mesh);
    let force_rate = vec![0.0; force.len()];
    let mut objective_gradient = vec![0.0; force.len()];
    objective_gradient[3 * mesh.load_dof + 1] = 1.0;
    differentiate(
        MechanicalDerivative {
            stiffness: &sensitivity.stiffness,
            stiffness_rate: &sensitivity.derivative,
            force: &force,
            force_rate: &force_rate,
            lifts: lifts.each_ref(),
            lift_rates: rates.each_ref().map(Vec::as_slice),
            objective_gradient: &objective_gradient,
            objective_partial: 0.0,
            fixed_regime: true,
        },
        1e-13,
        5000,
    )
}

fn verify_stiffness_derivative(mesh: &MovingMesh, order: usize) -> Result<f64, &'static str> {
    let sensitivity = assemble_sensitivity(
        &mesh.space,
        &mesh.jets,
        MaterialJet {
            young_mpa: 11.0,
            poisson: NU,
            d_young_mpa: 1.0,
            d_poisson: 0.0,
        },
        order,
    )?;
    let h_num = 1_i128;
    let h_den = 8192_i128;
    let plus = build_mesh(h_den + h_num, h_den)?;
    let minus = build_mesh(h_den - h_num, h_den)?;
    let kp = assemble(&plus.space, 10.0 + plus.t, NU, order)?;
    let km = assemble(&minus.space, 10.0 + minus.t, NU, order)?;
    let h = 1.0 / h_den as f64;
    if kp.offsets != km.offsets
        || kp.offsets != sensitivity.derivative.offsets
        || kp.columns != km.columns
        || kp.columns != sensitivity.derivative.columns
    {
        return Err("stiffness finite-difference sparsity changed");
    }
    let mut error2 = 0.0;
    let mut reference2 = 0.0;
    for ((p, m), d) in kp
        .values
        .iter()
        .zip(&km.values)
        .zip(&sensitivity.derivative.values)
    {
        for row in 0..3 {
            for col in 0..3 {
                let fd = (p[row][col] - m[row][col]) / (2.0 * h);
                error2 += (fd - d[row][col]).powi(2);
                reference2 += d[row][col].powi(2);
            }
        }
    }
    Ok(error2.sqrt() / reference2.sqrt().max(1e-30))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mesh = build_mesh(1, 1)?;
    assert_eq!(mesh.space.blocks.len(), 3);
    assert_eq!(mesh.space.ndof, 72);
    assert_eq!(mesh.equations.len(), 18);
    verify_physical_maps(&mesh);
    let matrix_fd_error = verify_stiffness_derivative(&mesh, 8)?;
    println!("Kprime centered-FD relative Frobenius error={matrix_fd_error:.3e}");
    assert!(
        matrix_fd_error <= 1e-7,
        "Kprime FD error {matrix_fd_error:e}"
    );
    let mut records = Vec::new();
    let mut base_gradient = None;
    for order in [4usize, 6, 8] {
        let result = derivative_result(&mesh, order, false, false, false, false)?;
        let (base_objective, _) = objective(&mesh, order)?;
        let direct_adjoint_error = (result.objective_tangent - result.objective_adjoint).abs();
        assert!(
            direct_adjoint_error <= 1e-9,
            "q={order} direct/adjoint mismatch {direct_adjoint_error:e}"
        );
        assert!(
            result.tangent_relative_residual <= 1e-10 && result.adjoint_relative_residual <= 1e-10
        );
        println!("q={order} J={base_objective:.12e} Jprime_direct={:.12e} Jprime_adjoint={:.12e} agreement={direct_adjoint_error:.3e} tangent_res={:.3e} adjoint_res={:.3e}",
            result.objective_tangent, result.objective_adjoint, result.tangent_relative_residual, result.adjoint_relative_residual);
        let mut fd = Vec::new();
        let mut remainders = Vec::new();
        for den in [64_i128, 128, 256, 512, 1024] {
            let plus = build_mesh(den + 1, den)?;
            let minus = build_mesh(den - 1, den)?;
            let jp = objective(&plus, order)?.0;
            let jm = objective(&minus, order)?.0;
            let h = 1.0 / den as f64;
            let centered = (jp - jm) / (2.0 * h);
            let error = (centered - result.objective_adjoint).abs();
            let remainder = (jp - base_objective - h * result.objective_adjoint).abs();
            fd.push(error);
            remainders.push(remainder);
            println!("q={order} h={h:.9e} Jprime_fd={centered:.12e} derivative_error={error:.3e} taylor_remainder={remainder:.3e}");
        }
        let asymptotic_ratio = remainders[3] / remainders[4];
        assert!(
            asymptotic_ratio >= 3.5,
            "q={order} Taylor ratio {asymptotic_ratio}"
        );
        assert!(
            fd[4] <= 1e-6,
            "q={order} finest centered FD error {}",
            fd[4]
        );
        records.push((order, result.objective_adjoint, fd[4], asymptotic_ratio));
        if order == 6 {
            base_gradient = Some(result.objective_adjoint);
        }
    }
    let q6_gradient = base_gradient.unwrap();
    let q8_gradient = records.iter().find(|r| r.0 == 8).unwrap().1;
    println!(
        "quadrature q6_to_q8 objective derivative relative change={:.3e}",
        (q8_gradient - q6_gradient).abs() / q8_gradient.abs().max(1e-30)
    );
    for (label, result) in [
        (
            "omit_geometry_Kprime",
            derivative_result(&mesh, 6, true, false, false, false)?,
        ),
        (
            "omit_trace_Cprime",
            derivative_result(&mesh, 6, false, true, false, false)?,
        ),
        (
            "omit_lifting_gprime",
            derivative_result(&mesh, 6, false, false, true, false)?,
        ),
        (
            "omit_material_Eprime",
            derivative_result(&mesh, 6, false, false, false, true)?,
        ),
    ] {
        let delta = result.objective_adjoint - q6_gradient;
        println!(
            "diagnostic {label}: predicted_Jprime={:.12e} delta={delta:.6e}",
            result.objective_adjoint
        );
    }
    Ok(())
}
