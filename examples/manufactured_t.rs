//! Manufactured vector elasticity solution on a two-rate Q2 T-junction mesh.
//!
//! The mesh is built from affine rational Bernstein cells with one dyadic
//! homogeneous scale. The left half has n^3 cells and the right half n*2n*n.
//! The coarse interface is restricted exactly onto two fine slave faces.

use solidcraft_trivariate::analysis_space::{AnalysisBlock, AnalysisSpace};
use solidcraft_trivariate::constrained_elasticity::solve_constrained;
use solidcraft_trivariate::constraint_system::{resolve, Lift, Ratio};
use solidcraft_trivariate::elasticity::{assemble, basis_geometry, gauss};
use solidcraft_trivariate::face_constraints::{nested_constraints, FaceNet, Span};
use solidcraft_trivariate::incidence::Orientation;
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

type Point = [i64; 3];
type Equation = (usize, solidcraft_trivariate::face_constraints::Row);
const YOUNG: f64 = 10.0;
const POISSON: f64 = 0.25;
const ALPHA: f64 = std::f64::consts::FRAC_PI_2;
const AMPLITUDE: [f64; 3] = [1.0, 2.0, -1.0];
const AFFINE: [[i64; 3]; 3] = [[2, -1, 1], [1, 3, 2], [-1, 1, -2]];
const OFFSET: [i64; 3] = [1, -2, 3];

#[derive(Clone)]
struct Cell {
    half: usize,
    ix: usize,
    iy: usize,
    iz: usize,
    block: AnalysisBlock,
}

type Mesh = (AnalysisSpace, Vec<Point>, Vec<Cell>, Vec<Equation>);

fn force_at([x, y, z]: [f64; 3]) -> [f64; 3] {
    let (sx, cx) = (ALPHA * x).sin_cos();
    let (sy, cy) = (ALPHA * y).sin_cos();
    let (sz, cz) = (std::f64::consts::PI * z).sin_cos();
    let phi = sx * sy * sz;
    let hessian = [
        [
            -ALPHA * ALPHA * phi,
            ALPHA * ALPHA * cx * cy * sz,
            ALPHA * std::f64::consts::PI * cx * sy * cz,
        ],
        [
            ALPHA * ALPHA * cx * cy * sz,
            -ALPHA * ALPHA * phi,
            ALPHA * std::f64::consts::PI * sx * cy * cz,
        ],
        [
            ALPHA * std::f64::consts::PI * cx * sy * cz,
            ALPHA * std::f64::consts::PI * sx * cy * cz,
            -std::f64::consts::PI.powi(2) * phi,
        ],
    ];
    let mu = YOUNG / (2.0 * (1.0 + POISSON));
    let lambda = YOUNG * POISSON / ((1.0 + POISSON) * (1.0 - 2.0 * POISSON));
    let k2 = 2.0 * ALPHA * ALPHA + std::f64::consts::PI.powi(2);
    std::array::from_fn(|i| {
        mu * k2 * AMPLITUDE[i] * phi
            - (lambda + mu) * (0..3).map(|j| AMPLITUDE[j] * hessian[i][j]).sum::<f64>()
    })
}

fn exact([x, y, z]: [f64; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let (sx, cx) = (ALPHA * x).sin_cos();
    let (sy, cy) = (ALPHA * y).sin_cos();
    let (sz, cz) = (std::f64::consts::PI * z).sin_cos();
    let phi = sx * sy * sz;
    let dphi = [
        ALPHA * cx * sy * sz,
        ALPHA * sx * cy * sz,
        std::f64::consts::PI * sx * sy * cz,
    ];
    let u = std::array::from_fn(|c| {
        OFFSET[c] as f64
            + (0..3)
                .map(|a| AFFINE[c][a] as f64 * [x, y, z][a])
                .sum::<f64>()
            + AMPLITUDE[c] * phi
    });
    let grad = std::array::from_fn(|c| {
        std::array::from_fn(|a| AFFINE[c][a] as f64 + AMPLITUDE[c] * dphi[a])
    });
    (u, grad)
}

fn affine_value(point: Point, denominator: i128, component: usize) -> Ratio {
    let numerator = i128::from(OFFSET[component]) * denominator
        + (0..3)
            .map(|a| i128::from(AFFINE[component][a]) * i128::from(point[a]))
            .sum::<i128>();
    Ratio::new(numerator, denominator).unwrap()
}

fn lifts(
    space: &AnalysisSpace,
    points: &[Point],
    denominator: i128,
    equations: &[Equation],
) -> Result<[Lift; 3], &'static str> {
    let prescribed: Vec<Vec<(usize, Ratio)>> = (0..3)
        .map(|c| {
            points
                .iter()
                .enumerate()
                .filter_map(|(dof, p)| {
                    let boundary = p[0] == 0
                        || p[0] == 2 * denominator as i64
                        || p[1] == 0
                        || p[1] == 2 * denominator as i64
                        || p[2] == 0
                        || p[2] == denominator as i64;
                    boundary.then(|| (dof, affine_value(*p, denominator, c)))
                })
                .collect()
        })
        .collect();
    Ok([
        resolve(space.ndof, equations, &prescribed[0])?,
        resolve(space.ndof, equations, &prescribed[1])?,
        resolve(space.ndof, equations, &prescribed[2])?,
    ])
}

fn build_mesh(n: usize) -> Result<Mesh, &'static str> {
    let denominator = 2 * n as i64;
    let mut ids = BTreeMap::<(usize, Point), usize>::new();
    let mut points = Vec::<Point>::new();
    let mut cells = Vec::<Cell>::new();
    let mut half_nodes = [BTreeSet::new(), BTreeSet::new()];
    for (half, half_node_set) in half_nodes.iter_mut().enumerate() {
        let y_cells = if half == 0 { n } else { 2 * n };
        for ix in 0..n {
            for iy in 0..y_cells {
                for iz in 0..n {
                    let mut controls = Vec::with_capacity(27);
                    let mut dofs = Vec::with_capacity(27);
                    for i in 0..3 {
                        for j in 0..3 {
                            for k in 0..3 {
                                let x = (2 * (half * n + ix) + i) as i64;
                                let y = if half == 0 {
                                    (4 * iy + 2 * j) as i64
                                } else {
                                    (2 * iy + j) as i64
                                };
                                let z = (2 * iz + k) as i64;
                                let p = [x, y, z];
                                let dof = *ids.entry((half, p)).or_insert_with(|| {
                                    let next = points.len();
                                    points.push(p);
                                    next
                                });
                                half_node_set.insert(dof);
                                controls.push([
                                    i128::from(x),
                                    i128::from(y),
                                    i128::from(z),
                                    i128::from(denominator),
                                ]);
                                dofs.push(dof);
                            }
                        }
                    }
                    cells.push(Cell {
                        half,
                        ix,
                        iy,
                        iz,
                        block: AnalysisBlock {
                            id: format!("h{half}-x{ix}-y{iy}-z{iz}"),
                            controls,
                            dofs,
                        },
                    });
                }
            }
        }
    }
    let mut block_index = BTreeMap::new();
    for (index, cell) in cells.iter().enumerate() {
        block_index.insert((cell.half, cell.ix, cell.iy, cell.iz), index);
    }
    let mut equations = Vec::new();
    for iy in 0..n {
        for iz in 0..n {
            let master = &cells[block_index[&(0, n - 1, iy, iz)]].block;
            let master_net = FaceNet {
                degrees: [2, 2],
                controls: (0..3)
                    .flat_map(|j| {
                        (0..3).map(move |k| {
                            let p = master.controls[18 + 3 * j + k];
                            [p[0] as i64, p[1] as i64, p[2] as i64, p[3] as i64]
                        })
                    })
                    .collect(),
            };
            for child in 0..2 {
                let slave = &cells[block_index[&(1, 0, 2 * iy + child, iz)]].block;
                let slave_net = FaceNet {
                    degrees: [2, 2],
                    controls: (0..3)
                        .flat_map(|j| {
                            (0..3).map(move |k| {
                                let p = slave.controls[3 * j + k];
                                [p[0] as i64, p[1] as i64, p[2] as i64, p[3] as i64]
                            })
                        })
                        .collect(),
                };
                let span = if child == 0 { Span::Lower } else { Span::Upper };
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
                for (local, row) in rows.into_iter().enumerate() {
                    let dependent = slave.dofs[local];
                    let mapped = solidcraft_trivariate::face_constraints::Row {
                        terms: row
                            .terms
                            .into_iter()
                            .map(|(m, weight)| (master.dofs[18 + m], weight))
                            .collect(),
                        denominator: row.denominator,
                    };
                    equations.push((dependent, mapped));
                }
            }
        }
    }
    let boundary = points
        .iter()
        .enumerate()
        .filter_map(|(dof, p)| {
            (p[0] == 0
                || p[0] == 2 * denominator
                || p[1] == 0
                || p[1] == 2 * denominator
                || p[2] == 0
                || p[2] == denominator)
                .then_some(dof)
        })
        .collect();
    let space = AnalysisSpace {
        ndof: points.len(),
        blocks: cells.iter().map(|c| c.block.clone()).collect(),
        boundary_dofs: BTreeMap::from([("outer".to_string(), boundary)]),
    };
    let expected_left = (2 * n + 1).pow(3);
    let expected_right = (2 * n + 1) * (4 * n + 1) * (2 * n + 1);
    let expected_overlap = 0;
    assert_eq!(half_nodes[0].len(), expected_left);
    assert_eq!(half_nodes[1].len(), expected_right);
    assert_eq!(
        half_nodes[0].intersection(&half_nodes[1]).count(),
        expected_overlap
    );
    assert_eq!(space.ndof, expected_left + expected_right);
    Ok((space, points, cells, equations))
}

fn force(space: &AnalysisSpace, order: usize) -> Result<Vec<f64>, &'static str> {
    let mut rhs = vec![0.0; 3 * space.ndof];
    for block in &space.blocks {
        for (q, weight) in gauss(order)? {
            for (r, weight_r) in gauss(order)? {
                for (s, weight_s) in gauss(order)? {
                    let sample = basis_geometry(block, [q, r, s])?;
                    let f = force_at(sample.position);
                    let scale = sample.determinant * weight * weight_r * weight_s;
                    for (local, &dof) in block.dofs.iter().enumerate() {
                        for c in 0..3 {
                            rhs[3 * dof + c] += scale * sample.values[local] * f[c];
                        }
                    }
                }
            }
        }
    }
    Ok(rhs)
}

fn errors(
    space: &AnalysisSpace,
    solution: &[f64],
    order: usize,
) -> Result<(f64, f64), &'static str> {
    let points = gauss(order)?;
    let (mut l2, mut h1) = (0.0, 0.0);
    for block in &space.blocks {
        for &(u, wu) in &points {
            for &(v, wv) in &points {
                for &(w, ww) in &points {
                    let sample = basis_geometry(block, [u, v, w])?;
                    let (exact_u, exact_grad) = exact(sample.position);
                    let mut approx_u = [0.0; 3];
                    let mut approx_grad = [[0.0; 3]; 3];
                    for (local, &dof) in block.dofs.iter().enumerate() {
                        for c in 0..3 {
                            let value = solution[3 * dof + c];
                            approx_u[c] += sample.values[local] * value;
                            for (a, gradient) in approx_grad[c].iter_mut().enumerate() {
                                *gradient += sample.gradients[local][a] * value;
                            }
                        }
                    }
                    let measure = sample.determinant * wu * wv * ww;
                    l2 += measure
                        * (0..3)
                            .map(|c| (approx_u[c] - exact_u[c]).powi(2))
                            .sum::<f64>();
                    h1 += measure
                        * (0..3)
                            .flat_map(|c| {
                                (0..3).map(move |a| (approx_grad[c][a] - exact_grad[c][a]).powi(2))
                            })
                            .sum::<f64>();
                }
            }
        }
    }
    Ok((l2.sqrt(), h1.sqrt()))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut reference_n4 = None;
    for n in [1usize, 2, 4, 8] {
        let started = Instant::now();
        let (space, points, cells, equations) = build_mesh(n)?;
        let den = (2 * n) as i128;
        let lift = lifts(&space, &points, den, &equations)?;
        let lift_refs = [&lift[0], &lift[1], &lift[2]];
        let mut level = Vec::new();
        for q in [4usize, 6] {
            let assembled_at = Instant::now();
            let stiffness = assemble(&space, YOUNG, POISSON, q)?;
            let rhs = force(&space, q)?;
            let result = solve_constrained(&stiffness, lift_refs, &rhs, 1e-13, 200_000)?;
            let err8 = errors(&space, &result.displacement_mm, 8)?;
            let err10 = errors(&space, &result.displacement_mm, 10)?;
            assert!(
                result.relative_residual <= 1e-10,
                "n={n} q={q} residual {}",
                result.relative_residual
            );
            let error_variation = (err10.0 - err8.0).abs() / err10.0.max(1e-30);
            let h1_variation = (err10.1 - err8.1).abs() / err10.1.max(1e-30);
            assert!(error_variation < 0.01 && h1_variation < 0.01,
                    "error quadrature variation n={n} q={q}: L2={error_variation:e}, H1={h1_variation:e}");
            println!("n={n} q={q} ndof={} cells={} constraints={} free={}/{}/{} residual={:.3e} iters={} energy={:.9e} L2(q8)={:.9e} H1(q8)={:.9e} L2(q10)={:.9e} H1(q10)={:.9e} qseconds={:.3}",
                space.ndof, cells.len(), equations.len(), lift[0].free_roots.len(), lift[1].free_roots.len(), lift[2].free_roots.len(),
                result.relative_residual, result.iterations, result.strain_energy_n_mm,
                err8.0, err8.1, err10.0, err10.1, assembled_at.elapsed().as_secs_f64());
            level.push((q, err10, result));
        }
        let patch_rhs = vec![0.0; 3 * space.ndof];
        let patch_k = assemble(&space, YOUNG, POISSON, 4)?;
        let patch = solve_constrained(&patch_k, lift_refs, &patch_rhs, 1e-13, 200_000)?;
        let mut patch_error: f64 = 0.0;
        for (dof, p) in points.iter().enumerate() {
            for c in 0..3 {
                patch_error = patch_error.max(
                    (patch.displacement_mm[3 * dof + c]
                        - affine_value(*p, den, c).numerator() as f64
                            / affine_value(*p, den, c).denominator() as f64)
                        .abs(),
                );
            }
        }
        assert!(
            patch_error <= 1e-9,
            "affine patch error n={n}: {patch_error:e}"
        );
        println!(
            "n={n} affine_patch_max={patch_error:.3e} time_s={:.3}",
            started.elapsed().as_secs_f64()
        );
        if level.len() == 2 {
            let e = level[1].1;
            let previous = level[0].1;
            let l2_variation = (e.0 - previous.0).abs() / e.0.max(1e-30);
            let h1_variation = (e.1 - previous.1).abs() / e.1.max(1e-30);
            println!(
                "n={n} q4_to_q6 relative_error_variation L2={:.3e} H1={:.3e}",
                l2_variation, h1_variation
            );
            assert!(
                l2_variation < 0.01 && h1_variation < 0.01,
                "q4/q6 error variation exceeds 1% at n={n}"
            );
        }
        if n == 4 || n == 8 {
            let q6 = level[1].1;
            println!("n={n} final_error L2={:.9e} H1={:.9e}", q6.0, q6.1);
            if n == 4 {
                reference_n4 = Some(q6);
            } else if let Some(coarse) = reference_n4 {
                let l2_rate = (coarse.0 / q6.0).log2();
                let h1_rate = (coarse.1 / q6.1).log2();
                println!(
                    "n4_to_n8 rates L2={l2_rate:.6} (target >=2.7) H1={h1_rate:.6} (target >=1.8)"
                );
                assert!(
                    l2_rate >= 2.7,
                    "M1 L2 convergence rate {l2_rate} is below 2.7"
                );
                assert!(
                    h1_rate >= 1.8,
                    "M1 H1 convergence rate {h1_rate} is below 1.8"
                );
            }
        }
    }
    Ok(())
}
