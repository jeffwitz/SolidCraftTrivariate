//! Uniform cubic-spline convergence study on one affine cube parent.
use crate::{
    constraint_system::{Lift, LiftLine, Ratio},
    extracted_space::{AxisExtraction, ExtractedElement},
    json::{object, Json},
    parent_field::axis_rows,
    section_manufactured::{solve_sampled, IntegrationDomain},
};
use std::collections::BTreeSet;

type Result<T> = std::result::Result<T, &'static str>;
const DEGREE: usize = 3;

fn n(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}

fn open_knots(cells: usize, interior_multiplicity: usize) -> Result<Vec<Ratio>> {
    if cells == 0 || !(1..=DEGREE).contains(&interior_multiplicity) {
        return Err("cube grid needs positive cells and knot multiplicity in 1..=3");
    }
    let mut knots = vec![Ratio::zero(); DEGREE + 1];
    for knot in 1..cells {
        knots.extend(std::iter::repeat_n(
            Ratio::new(knot as i128, 1)?,
            interior_multiplicity,
        ));
    }
    knots.extend(std::iter::repeat_n(
        Ratio::new(cells as i128, 1)?,
        DEGREE + 1,
    ));
    Ok(knots)
}

fn axis_element(knots: &[Ratio], span: usize) -> Result<(AxisExtraction, Vec<usize>)> {
    let rows = axis_rows(
        knots,
        DEGREE,
        [
            Ratio::new(span as i128, 1)?,
            Ratio::new((span + 1) as i128, 1)?,
        ],
    )?;
    let active: Vec<usize> = rows
        .iter()
        .flat_map(|row| row.keys().copied())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut row_offsets = vec![0];
    let mut columns = Vec::new();
    let mut coefficients = Vec::new();
    for row in &rows {
        for (&global, &coefficient) in row {
            columns.push(
                active
                    .binary_search(&global)
                    .map_err(|_| "inactive spline coefficient")?,
            );
            coefficients.push(coefficient);
        }
        row_offsets.push(columns.len());
    }
    let extraction = AxisExtraction {
        degree: DEGREE,
        basis_count: active.len(),
        row_offsets,
        columns,
        coefficients,
    };
    extraction.validate()?;
    Ok((extraction, active))
}

fn greville(knots: &[Ratio], control: usize) -> Result<f64> {
    let sum = (1..=DEGREE).try_fold(Ratio::zero(), |acc, offset| {
        acc.checked_add(knots[control + offset])
    })?;
    Ok(sum.numerator() as f64 / sum.denominator() as f64 / DEGREE as f64)
}

fn identity_zero_lift(ndof: usize, shape: [usize; 3]) -> Lift {
    let constrained: BTreeSet<usize> = (0..shape[0])
        .flat_map(|i| (0..shape[1]).map(move |j| (i, j)))
        .map(|(i, j)| (i * shape[1] + j) * shape[2])
        .collect();
    let free_roots: Vec<usize> = (0..ndof).filter(|id| !constrained.contains(id)).collect();
    let lines = (0..ndof)
        .map(|dof| {
            let terms = free_roots
                .binary_search(&dof)
                .ok()
                .map(|free| vec![(free, Ratio::one())])
                .unwrap_or_default();
            LiftLine {
                dof,
                terms,
                offset: Ratio::zero(),
            }
        })
        .collect();
    Lift { free_roots, lines }
}

struct Mesh {
    elements: Vec<ExtractedElement>,
    scalar_dofs: usize,
    boundary: Vec<(usize, usize, bool)>,
    lift: Lift,
    min_det: f64,
    max_geometry_error: f64,
    max_jacobian_error: f64,
    max_affine_value_error: f64,
    max_affine_gradient_error: f64,
}

fn build_mesh(cells: usize, multiplicity: usize) -> Result<Mesh> {
    let knots = open_knots(cells, multiplicity)?;
    let basis_count = knots.len() - DEGREE - 1;
    let shape = [basis_count; 3];
    let scalar_dofs = basis_count
        .checked_pow(3)
        .ok_or("cube spline DOF overflow")?;
    let one_d_xyz: Vec<f64> = (0..basis_count)
        .map(|i| greville(&knots, i).map(|x| x / cells as f64))
        .collect::<Result<_>>()?;
    let mut dof_xyz = Vec::with_capacity(scalar_dofs);
    for &x in &one_d_xyz {
        for &y in &one_d_xyz {
            for &z in &one_d_xyz {
                dof_xyz.push([x, y, z]);
            }
        }
    }

    let mut axis_factors = Vec::with_capacity(cells);
    for span in 0..cells {
        axis_factors.push(axis_element(&knots, span)?);
    }
    let mut elements = Vec::with_capacity(cells.pow(3));
    let mut boundary = Vec::with_capacity(6 * cells * cells);
    let mut min_det = f64::INFINITY;
    let mut max_geometry_error = 0.0_f64;
    let mut max_jacobian_error = 0.0_f64;
    let mut max_affine_value_error = 0.0_f64;
    let mut max_affine_gradient_error = 0.0_f64;
    let gradient = [
        [0.013, -0.004, 0.002],
        [0.001, 0.009, -0.003],
        [-0.002, 0.005, 0.006],
    ];
    let offset = [0.7, -0.2, 0.3];
    let check_points = [0.0, 0.25, 0.5, 0.75, 1.0];
    for ix in 0..cells {
        for iy in 0..cells {
            for iz in 0..cells {
                let cell = elements.len();
                let factors = [&axis_factors[ix], &axis_factors[iy], &axis_factors[iz]];
                let axes = [
                    factors[0].0.clone(),
                    factors[1].0.clone(),
                    factors[2].0.clone(),
                ];
                let active = [&factors[0].1, &factors[1].1, &factors[2].1];
                let mut dofs = Vec::new();
                let mut geometry_controls = Vec::new();
                for &i in active[0] {
                    for &j in active[1] {
                        for &k in active[2] {
                            let id = (i * shape[1] + j) * shape[2] + k;
                            dofs.push(id);
                            let [x, y, z] = dof_xyz[id];
                            geometry_controls.push([x, y, z, 1.]);
                        }
                    }
                }
                let field_weights = vec![1.; dofs.len()];
                let element = ExtractedElement {
                    axes,
                    dofs,
                    field_weights,
                    geometry_controls,
                };
                element.validate()?;
                for &u in &check_points {
                    for &v in &check_points {
                        for &w in &check_points {
                            let (geom, basis) = element.sample([u, v, w])?;
                            min_det = min_det.min(geom.determinant);
                            let expected = [
                                (ix as f64 + u) / cells as f64,
                                (iy as f64 + v) / cells as f64,
                                (iz as f64 + w) / cells as f64,
                            ];
                            for (actual, target) in geom.position.iter().zip(expected) {
                                max_geometry_error =
                                    max_geometry_error.max((*actual - target).abs());
                            }
                            for (a, column) in geom.jacobian_columns.iter().enumerate() {
                                for (b, actual) in column.iter().enumerate() {
                                    let target = if a == b { 1.0 / cells as f64 } else { 0. };
                                    max_jacobian_error =
                                        max_jacobian_error.max((*actual - target).abs());
                                }
                            }
                            let mut value = [0.; 3];
                            let mut grad = [[0.; 3]; 3];
                            for local in 0..basis.dofs.len() {
                                let point = dof_xyz[basis.dofs[local]];
                                for a in 0..3 {
                                    let coefficient = offset[a]
                                        + (0..3).map(|b| gradient[a][b] * point[b]).sum::<f64>();
                                    value[a] += basis.values[local] * coefficient;
                                    for (b, component) in grad[a].iter_mut().enumerate() {
                                        *component +=
                                            basis.physical_gradients[local][b] * coefficient;
                                    }
                                }
                            }
                            for a in 0..3 {
                                let target = offset[a]
                                    + (0..3)
                                        .map(|b| gradient[a][b] * geom.position[b])
                                        .sum::<f64>();
                                max_affine_value_error =
                                    max_affine_value_error.max((value[a] - target).abs());
                                for b in 0..3 {
                                    max_affine_gradient_error = max_affine_gradient_error
                                        .max((grad[a][b] - gradient[a][b]).abs());
                                }
                            }
                        }
                    }
                }
                elements.push(element);
                if ix == 0 {
                    boundary.push((cell, 0, false));
                }
                if ix + 1 == cells {
                    boundary.push((cell, 0, true));
                }
                if iy == 0 {
                    boundary.push((cell, 1, false));
                }
                if iy + 1 == cells {
                    boundary.push((cell, 1, true));
                }
                if iz == 0 {
                    boundary.push((cell, 2, false));
                }
                if iz + 1 == cells {
                    boundary.push((cell, 2, true));
                }
            }
        }
    }
    let lift = identity_zero_lift(scalar_dofs, shape);
    Ok(Mesh {
        elements,
        scalar_dofs,
        boundary,
        lift,
        min_det,
        max_geometry_error,
        max_jacobian_error,
        max_affine_value_error,
        max_affine_gradient_error,
    })
}

/// Solve one maximum-continuity and one full-Bernstein uniform partition per level.
pub fn run(levels: &[usize], quadrature: usize) -> Result<Json> {
    if levels.is_empty() || levels.windows(2).any(|w| w[0] >= w[1]) || quadrature < 2 {
        return Err("cube study requires increasing levels and quadrature >= 2");
    }
    let mut rows = Vec::new();
    for &(continuity, multiplicity, label) in &[(2, 1, "C2_maximum"), (0, 3, "C0_bernstein")] {
        for &cells in levels {
            let mesh = build_mesh(cells, multiplicity)?;
            let domain = IntegrationDomain {
                scalar_dofs: mesh.scalar_dofs,
                cell_count: mesh.elements.len(),
                z_bounds: [0., 1.],
                boundary: &mesh.boundary,
                essential: &mesh.lift,
            };
            let result = solve_sampled(
                &domain,
                &object([
                    ("schema", Json::String("section_axial_sine_mms_v1".into())),
                    ("quadrature", n(quadrature)),
                ]),
                |cell, uvw| mesh.elements[cell].sample(uvw),
            )?;
            rows.push(object([
                ("continuity_label", Json::String(label.into())),
                ("interior_continuity", n(continuity)),
                ("cells_per_axis", n(cells)),
                ("cell_count", n(mesh.elements.len())),
                ("scalar_dofs", n(mesh.scalar_dofs)),
                ("boundary_face_observations", n(mesh.boundary.len())),
                ("min_sampled_determinant", n(mesh.min_det)),
                ("max_geometry_position_error_mm", n(mesh.max_geometry_error)),
                ("max_jacobian_error", n(mesh.max_jacobian_error)),
                ("max_affine_value_error_mm", n(mesh.max_affine_value_error)),
                (
                    "max_affine_gradient_error",
                    n(mesh.max_affine_gradient_error),
                ),
                ("result", result),
            ]));
        }
    }
    Ok(object([
        (
            "schema",
            Json::String("solidcraft_uniform_cubic_cube_convergence_v1".into()),
        ),
        (
            "geometry_parent",
            Json::String("single_affine_unit_cube".into()),
        ),
        ("spline_degree", n(DEGREE)),
        ("quadrature_order", n(quadrature)),
        (
            "levels_per_axis",
            Json::Array(levels.iter().map(|&x| n(x)).collect()),
        ),
        ("observations", Json::Array(rows)),
    ]))
}
