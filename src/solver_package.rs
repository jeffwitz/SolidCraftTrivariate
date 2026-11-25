//! Deterministic solver package for audited fixed and parameterized brackets.
//!
//! The source model is audited before any numerical arrays are exported. Its
//! exact homogeneous geometry is the solver input; elevated controls, topology,
//! boundary selections, quadrature data and the Dirichlet reduction are
//! provided as reproducible companion data.

use crate::{
    analysis_space::{self, AnalysisSpace},
    bracket::BracketBlock,
    bracket_json,
    bracket_solve::bracket_load,
    complex_skeleton::encode_skeleton,
    elasticity::{basis_geometry, gauss},
    interchange,
    json::{object, Json},
};
use std::collections::BTreeMap;

const YOUNG_MPA: f64 = 210_000.0;
const POISSON: f64 = 0.3;
const TOTAL_LOAD_N: f64 = 1_000.0;
const TARGET_AREA_MM2: f64 = 320.0;

fn text(s: &str) -> Json {
    Json::String(s.to_owned())
}
fn number(n: impl ToString) -> Json {
    Json::Number(n.to_string())
}
fn boolean(b: bool) -> Json {
    Json::Bool(b)
}
fn arr<T>(values: impl IntoIterator<Item = T>, mut f: impl FnMut(T) -> Json) -> Json {
    Json::Array(values.into_iter().map(&mut f).collect())
}
fn numeric_array<T: ToString>(values: impl IntoIterator<Item = T>) -> Json {
    arr(values, number)
}
fn f64_array(values: impl IntoIterator<Item = f64>) -> Json {
    arr(values, number)
}

pub(super) fn analysis_block_json(block: &analysis_space::AnalysisBlock) -> Json {
    object([
        ("id", text(&block.id)),
        ("degrees", numeric_array([2, 2, 2])),
        ("control_order", text("u-major,v-major,w-fast")),
        ("controls", arr(&block.controls, |p| numeric_array(*p))),
        ("dofs", numeric_array(&block.dofs)),
    ])
}

fn faces_for_feature(design_faces: &[Json], feature: &str) -> Vec<Json> {
    design_faces
        .iter()
        .filter(|face| face.field("feature").and_then(Json::text).ok() == Some(feature))
        .cloned()
        .collect()
}

fn gauss_json(order: usize) -> Result<Json, &'static str> {
    let rule = gauss(order)?;
    Ok(object([
        ("order", number(order)),
        ("domain", numeric_array([0, 1])),
        (
            "nodes_weights",
            arr(rule, |(node, weight)| f64_array([node, weight])),
        ),
        (
            "tensor_product_dimensions",
            object([("volume", number(3)), ("loaded_surface", number(2))]),
        ),
        (
            "accuracy_scope",
            text(
                "This is the requested Gauss rule; no order-independent accuracy guarantee is implied.",
            ),
        ),
    ]))
}

pub(super) fn integrate_volume(space: &AnalysisSpace, order: usize) -> Result<f64, &'static str> {
    let points = gauss(order)?;
    let mut volume = 0.0;
    for block in &space.blocks {
        for &(u, wu) in &points {
            for &(v, wv) in &points {
                for &(w, ww) in &points {
                    let sample = basis_geometry(block, [u, v, w])?;
                    let value = sample.determinant * wu * wv * ww;
                    if !value.is_finite() || value <= 0.0 {
                        return Err("invalid volume quadrature contribution");
                    }
                    volume += value;
                }
            }
        }
    }
    if !volume.is_finite() || volume <= 0.0 {
        return Err("invalid integrated volume");
    }
    Ok(volume)
}

fn material_json() -> Json {
    object([
        ("young_mpa", number(YOUNG_MPA)),
        ("poisson", number(POISSON)),
        (
            "units",
            object([
                ("young_modulus", text("MPa = N/mm^2")),
                ("length", text("mm")),
                ("force", text("N")),
            ]),
        ),
    ])
}

fn vectors_and_reduction(space: &AnalysisSpace) -> Result<(Json, Json), &'static str> {
    let fixed_scalar = space
        .boundary_dofs
        .get("y50")
        .ok_or("analysis space has no y50 selection")?;
    if fixed_scalar.is_empty() || fixed_scalar.iter().any(|&d| d >= space.ndof) {
        return Err("invalid y50 scalar DOFs");
    }
    let vector_size = space
        .ndof
        .checked_mul(3)
        .ok_or("vector DOF count overflow")?;
    let mut fixed_vector = Vec::with_capacity(fixed_scalar.len() * 3);
    for &scalar in fixed_scalar {
        for component in 0..3 {
            fixed_vector.push(3 * scalar + component);
        }
    }
    let mut fixed_mask = vec![false; vector_size];
    for &dof in &fixed_vector {
        fixed_mask[dof] = true;
    }
    let free_vector = fixed_mask
        .iter()
        .enumerate()
        .filter_map(|(i, fixed)| (!fixed).then_some(i))
        .collect::<Vec<_>>();
    let c_entries = free_vector
        .iter()
        .enumerate()
        .map(|(column, &row)| [row, column])
        .collect::<Vec<_>>();
    let clamping = object([
        ("feature", text("y50")),
        ("scalar_dofs", numeric_array(fixed_scalar)),
        ("vector_dofs", numeric_array(&fixed_vector)),
        (
            "component_order",
            Json::Array(["ux", "uy", "uz"].into_iter().map(text).collect()),
        ),
        ("value_mm", numeric_array([0, 0, 0])),
    ]);
    let reduction = object([
        ("full_vector_dofs", number(vector_size)),
        ("reduced_vector_dofs", number(free_vector.len())),
        ("ordering", text("scalar_node_major_xyz")),
        ("fixed_vector_dofs", numeric_array(&fixed_vector)),
        ("free_vector_dofs", numeric_array(&free_vector)),
        (
            "C",
            object([
                ("shape", numeric_array([vector_size, free_vector.len()])),
                ("representation", text("sparse_identity_injection")),
                ("unit_entries_row_column", arr(c_entries, numeric_array)),
                ("g_is_zero", boolean(true)),
            ]),
        ),
        ("clamped_nodes", clamping.clone()),
    ]);
    Ok((clamping, reduction))
}

fn derivative_json(model: &Json, input_blocks: &[BracketBlock]) -> Result<Json, &'static str> {
    let schema = model.field("schema")?.text()?;
    if schema == "filleted_bracket_v1" {
        return Ok(object([
            ("parameters", Json::Array(vec![])),
            (
                "scope",
                text("The fixed recipe exposes no radius derivatives."),
            ),
            ("blocks", Json::Array(vec![])),
            (
                "selection_derivative",
                text("zero: y50 and x60 selections are fixed"),
            ),
        ]));
    }
    if schema != "parameterized_bracket_v1" {
        return Err("unsupported derivative schema");
    }
    let radii = bracket_json::bracket_radii(model)?;
    let mut parameters = Vec::new();
    for (axis, name) in ["radius_x", "radius_y"].into_iter().enumerate() {
        let prefix = if axis == 0 {
            "hole_x_sector_"
        } else {
            "hole_y_sector_"
        };
        let (horizontal, center) = if axis == 0 { (0, 40) } else { (1, 32) };
        let mut jets = Vec::new();
        for block in input_blocks.iter().filter(|b| b.id.starts_with(prefix)) {
            // At fixed stored gauge, H_inner - W*center scales with radius.
            // A common denominator also represents arbitrary admissible gauges.
            let mut numerators = Vec::new();
            for (index, p) in block.controls.iter().enumerate() {
                let mut dh = [0i128; 4];
                if index / 2 % 2 == 1 {
                    let w = i128::from(p[3]);
                    dh[horizontal] = 16 * (i128::from(p[horizontal]) - center * w);
                    dh[2] = 16 * (i128::from(p[2]) - 20 * w);
                }
                numerators.push(dh);
            }
            let q = i128::from(radii[axis]);
            let divisor = numerators.iter().flatten().fold(q, |mut a, &x| {
                let mut b = x.abs();
                while b != 0 {
                    (a, b) = (b, a % b);
                }
                a
            });
            jets.push(object([
                ("block", text(&block.id)),
                ("degree", numeric_array(block.degrees)),
                ("control_order", text("u-major,v-major,w-fast")),
                ("derivative_denominator", number(q / divisor)),
                (
                    "d_homogeneous_controls_d_mm",
                    arr(numerators, |p| numeric_array(p.map(|x| x / divisor))),
                ),
            ]));
        }
        parameters.push(object([
            ("name",text(name)), ("input_ticks",number(radii[axis])), ("unit",text("mm")),
            ("control_derivatives",Json::Array(jets)),
            ("scope",text("Fixed topology and stored gauge. Divide each numerator by its block derivative_denominator.")),
        ]));
    }
    Ok(object([
        ("parameters", Json::Array(parameters)),
        (
            "scope",
            text(
                "Derivatives are with respect to physical radius in mm; each output jet uses the stored common gauge.",
            ),
        ),
        (
            "selection_derivative",
            text("zero: y50 and x60 selections and the identity reduction C are fixed"),
        ),
    ]))
}

/// Audit and package the exact fixed or R13 bracket for an independent solver.
pub fn export_analysis(model: &Json, order: usize) -> Result<Json, &'static str> {
    let schema = model.field("schema")?.text()?;
    if matches!(
        schema,
        "nurbs_multiloop_v1"
            | "g2_cross_plate_v1"
            | "primitive_sweep_v1"
            | "nurbs_edge_loops_v1"
            | "nurbs_sketch_strip_v1"
            | "nurbs_cross_bores_v1"
            | "nurbs_cross_bores_bounded_v1"
            | "nurbs_cross_oval_v1"
    ) {
        return crate::nurbs_space::analysis_package(model, order);
    }
    if matches!(
        schema,
        "sketch_layer_partition_v1" | "sketch_boolean_layers_v1" | "orthogonal_sketch_features_v1"
    ) {
        return export_sketch_analysis(model, order);
    }
    if schema != "filleted_bracket_v1" && schema != "parameterized_bracket_v1" {
        return Err("solver export supports only fixed and parameterized brackets");
    }
    let certificate = interchange::audit(model)?;
    if certificate.field("status")?.text()? != "valid" {
        return Err("solver input is not certified valid");
    }
    let mut blocks = bracket_json::decode_blocks(model)?;
    if blocks.len() != 36 {
        return Err("solver input must contain 36 bracket blocks");
    }
    blocks.sort_by(|a, b| a.id.cmp(&b.id));
    let space = if schema == "filleted_bracket_v1" {
        analysis_space::build_analysis_space(&blocks)?
    } else {
        let radii = bracket_json::bracket_radii(model)?;
        analysis_space::parameterized_bracket_space(radii[0], radii[1], &blocks)?
    };
    let skeleton = encode_skeleton(&space)?;
    let design_faces = model.field("design_faces")?.items()?.to_vec();
    let y50_faces = faces_for_feature(&design_faces, "y50");
    let x60_faces = faces_for_feature(&design_faces, "x60");
    if y50_faces.is_empty() || x60_faces.is_empty() {
        return Err("required y50/x60 design selections are missing");
    }
    let (y50_selection, reduction) = vectors_and_reduction(&space)?;
    let (rhs, area) = bracket_load(&space, order)?;
    if (area - TARGET_AREA_MM2).abs() > 1e-8 {
        return Err("loaded area does not match bracket selection");
    }
    let resultant = [
        rhs.iter().step_by(3).sum::<f64>(),
        rhs.iter().skip(1).step_by(3).sum::<f64>(),
        rhs.iter().skip(2).step_by(3).sum::<f64>(),
    ];
    if (resultant[1] + TOTAL_LOAD_N).abs() > 1e-8
        || resultant[0].abs() > 1e-10
        || resultant[2].abs() > 1e-10
    {
        return Err("integrated load resultant is invalid");
    }
    let gauss_rule = gauss_json(order)?;
    let volume = integrate_volume(&space, order)?;
    let (qx, qy) = if schema == "parameterized_bracket_v1" {
        let q = bracket_json::bracket_radii(model)?;
        (q[0], q[1])
    } else {
        (64, 64)
    };
    let rx = f64::from(qx) / 16.;
    let ry = f64::from(qy) / 16.;
    let ideal_volume =
        34080. - 360. * std::f64::consts::PI - 8. * std::f64::consts::PI * (rx * rx + ry * ry);
    let topology = model.field("topology")?.clone();
    let mut source_model = model.clone();
    if let Json::Object(fields) = &mut source_model {
        fields.insert("blocks".into(), bracket_json::encode_blocks(blocks.clone()));
    }
    let derivatives = derivative_json(model, &blocks)?;
    let scalar_ndof = space.ndof;
    let vector_ndof = scalar_ndof
        .checked_mul(3)
        .ok_or("vector DOF count overflow")?;
    let analysis = object([
        ("scalar_ndof", number(scalar_ndof)),
        ("vector_ndof", number(vector_ndof)),
        ("ordering", text("scalar_node_major_xyz")),
        ("block_control_order", text("u-major,v-major,w-fast")),
        (
            "basis",
            object([
                ("family", text("rational_tensor_bernstein_q2")),
                ("degrees", numeric_array([2, 2, 2])),
                (
                    "parametric_domain",
                    Json::Array(vec![numeric_array([0, 1]); 3]),
                ),
                (
                    "clamped_knot_vectors",
                    Json::Array(vec![numeric_array([0, 0, 0, 1, 1, 1]); 3]),
                ),
                (
                    "control_lattice_parameters",
                    Json::Array(vec![numeric_array([0., 0.5, 1.]); 3]),
                ),
                (
                    "extraction_identity",
                    object([
                        ("kind", text("identity")),
                        ("local_coefficients", number(27)),
                        ("global_mapping", text("analysis_space.blocks[].dofs")),
                        (
                            "rational_weights",
                            text("analysis_space.blocks[].controls[][3]"),
                        ),
                    ]),
                ),
            ]),
        ),
        ("blocks", arr(&space.blocks, analysis_block_json)),
        (
            "boundary_dofs",
            Json::Object(
                space
                    .boundary_dofs
                    .iter()
                    .map(|(name, dofs)| (name.clone(), numeric_array(dofs)))
                    .collect::<BTreeMap<_, _>>(),
            ),
        ),
    ]);
    let selections = object([
        (
            "y50",
            object([
                ("faces", Json::Array(y50_faces)),
                (
                    "scalar_dofs",
                    space
                        .boundary_dofs
                        .get("y50")
                        .map(numeric_array)
                        .unwrap_or(Json::Array(vec![])),
                ),
                ("value_mm", numeric_array([0, 0, 0])),
            ]),
        ),
        (
            "x60",
            object([
                ("faces", Json::Array(x60_faces)),
                ("area_mm2", number(area)),
                (
                    "traction_n_per_mm2",
                    numeric_array([0., -TOTAL_LOAD_N / TARGET_AREA_MM2, 0.]),
                ),
                ("resultant_n", f64_array(resultant)),
            ]),
        ),
    ]);
    let load = object([
        ("selection", text("x60")),
        (
            "traction_n_per_mm2",
            numeric_array([0., -TOTAL_LOAD_N / TARGET_AREA_MM2, 0.]),
        ),
        ("area_mm2", number(area)),
        ("resultant_n", f64_array(resultant)),
        ("nodal_rhs_n", f64_array(rhs)),
    ]);
    let measures = object([
        ("volume_quadrature_mm3", number(volume)),
        ("ideal_volume_target_mm3", number(ideal_volume)),
        ("x60_area_mm2", number(area)),
    ]);
    let geometry = object([
        ("length_unit", text("mm")),
        (
            "homogeneous_coordinate_order",
            Json::Array(["X", "Y", "Z", "W"].into_iter().map(text).collect()),
        ),
        ("control_order", text("u-major,v-major,w-fast")),
        ("source_blocks_json_pointer", text("/source_model/blocks")),
        (
            "source_coefficients",
            text("exact signed integers; source_model is authoritative"),
        ),
    ]);
    Ok(object([
        ("schema", text("solidcraft_solver_package_v1")),
        ("source_schema", text(schema)),
        ("quadrature_order", number(order)),
        ("source_model", source_model),
        ("certificate", certificate),
        ("topology", topology),
        ("design_faces", Json::Array(design_faces)),
        ("geometry", geometry),
        ("skeleton", skeleton),
        ("analysis_space", analysis),
        ("selections", selections),
        ("material", material_json()),
        ("load", load),
        ("quadrature", gauss_rule),
        ("reduction", reduction),
        ("clamping", y50_selection),
        ("measures", measures),
        ("derivatives", derivatives),
    ]))
}

fn side_from_json(value: &Json) -> Result<crate::incidence::BlockSide, &'static str> {
    let side = value.field("side")?.text()?;
    let (axis, upper) = match side {
        "u0" => (0, false),
        "u1" => (0, true),
        "v0" => (1, false),
        "v1" => (1, true),
        "w0" => (2, false),
        "w1" => (2, true),
        _ => return Err("unknown sketch face side"),
    };
    Ok(crate::incidence::BlockSide {
        block: value.field("block")?.text()?.to_owned(),
        axis,
        upper,
    })
}

fn sketch_design_faces(
    model: &Json,
) -> Result<Vec<(crate::incidence::BlockSide, String)>, &'static str> {
    model
        .field("design_faces")?
        .items()?
        .iter()
        .map(|face| {
            Ok((
                side_from_json(face)?,
                face.field("feature")?.text()?.to_owned(),
            ))
        })
        .collect()
}

fn control_positions(space: &AnalysisSpace) -> Result<Vec<[f64; 3]>, &'static str> {
    let mut positions: Vec<Option<[f64; 3]>> = vec![None; space.ndof];
    for block in &space.blocks {
        for (&dof, control) in block.dofs.iter().zip(&block.controls) {
            if dof >= space.ndof || control[3] <= 0 {
                return Err("invalid sketch analysis control mapping");
            }
            let point = std::array::from_fn(|axis| control[axis] as f64 / control[3] as f64);
            if let Some(previous) = positions[dof] {
                if (0..3).any(|axis| (previous[axis] - point[axis]).abs() > 1e-12) {
                    return Err("shared sketch Q2 control positions disagree");
                }
            } else {
                positions[dof] = Some(point);
            }
        }
    }
    positions
        .into_iter()
        .map(|point| point.ok_or("unreferenced sketch analysis control"))
        .collect()
}

fn vector_norm(vector: &[f64]) -> f64 {
    vector
        .iter()
        .fold(0.0_f64, |norm, value| norm.hypot(*value))
}

fn dot_vector(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn analysis_diagnostics(space: &AnalysisSpace) -> Result<Json, &'static str> {
    let positions = control_positions(space)?;
    let mut max_condition = 0.0_f64;
    let mut min_scaled_determinant = f64::INFINITY;
    let grid = [0.0, 0.5, 1.0];
    let mut samples = 0usize;
    for block in &space.blocks {
        for &u in &grid {
            for &v in &grid {
                for &w in &grid {
                    let sample = crate::elasticity::basis_geometry(block, [u, v, w])?;
                    let quality = crate::quality::at_jacobian(sample.jacobian_columns)?;
                    max_condition = max_condition.max(quality.condition_number);
                    min_scaled_determinant = min_scaled_determinant.min(quality.scaled_determinant);
                    samples += 1;
                }
            }
        }
    }

    // Use the public production assembler on the actual mounted Q2 space. The
    // small order is a patch/mode diagnostic, not a convergence study.
    let quadrature_order = 2;
    let young = 1.0;
    let poisson = 0.25;
    let stiffness = crate::elasticity::assemble(space, young, poisson, quadrature_order)?;
    let mut gradient = [[0.0; 3]; 3];
    gradient[0] = [0.013, -0.004, 0.002];
    gradient[1] = [0.001, 0.009, -0.003];
    gradient[2] = [-0.002, 0.005, 0.006];
    let offset = [0.7, -0.2, 0.3];
    let mut affine = vec![0.0; 3 * space.ndof];
    for (node, point) in positions.iter().enumerate() {
        for component in 0..3 {
            affine[3 * node + component] = offset[component]
                + (0..3)
                    .map(|axis| gradient[component][axis] * point[axis])
                    .sum::<f64>();
        }
    }
    let force = stiffness.apply(&affine)?;
    let energy = 0.5 * dot_vector(&affine, &force);
    let strain: [[f64; 3]; 3] =
        std::array::from_fn(|i| std::array::from_fn(|j| 0.5 * (gradient[i][j] + gradient[j][i])));
    let strain_norm2: f64 = strain.iter().flatten().map(|x| x * x).sum();
    let trace = strain[0][0] + strain[1][1] + strain[2][2];
    let mu = young / (2.0 * (1.0 + poisson));
    let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let energy_density = mu * strain_norm2 + 0.5 * lambda * trace * trace;
    let expected_energy = energy_density * stiffness.integrated_volume;
    let energy_error = (energy - expected_energy).abs() / expected_energy.abs().max(1e-30);

    let mut rigid_residuals = Vec::with_capacity(6);
    for mode in 0..6 {
        let mut displacement = vec![0.0; 3 * space.ndof];
        for (node, [x, y, z]) in positions.iter().copied().enumerate() {
            let value = match mode {
                0 => [1.0, 0.0, 0.0],
                1 => [0.0, 1.0, 0.0],
                2 => [0.0, 0.0, 1.0],
                3 => [0.0, -z, y],
                4 => [z, 0.0, -x],
                _ => [-y, x, 0.0],
            };
            displacement[3 * node..3 * node + 3].copy_from_slice(&value);
        }
        let residual = stiffness.apply(&displacement)?;
        let denominator =
            vector_norm(&displacement).max(1.0) * young * stiffness.integrated_volume.max(1.0);
        rigid_residuals.push(vector_norm(&residual) / denominator);
    }
    Ok(object([
        (
            "geometry_grid",
            object([
                ("axis_nodes", numeric_array(grid)),
                ("evaluations", number(samples)),
                ("max_condition_number", number(max_condition)),
                ("min_scaled_determinant", number(min_scaled_determinant)),
            ]),
        ),
        (
            "mounted_q2",
            object([
                ("assembler", text("elasticity::assemble")),
                ("quadrature_order", number(quadrature_order)),
                (
                    "material",
                    object([("young_mpa", number(young)), ("poisson", number(poisson))]),
                ),
                ("integrated_volume_mm3", number(stiffness.integrated_volume)),
                (
                    "affine_field_energy",
                    object([
                        (
                            "gradient",
                            Json::Array(gradient.into_iter().map(f64_array).collect()),
                        ),
                        ("offset_mm", f64_array(offset)),
                        ("assembled_n_mm", number(energy)),
                        ("constant_strain_expected_n_mm", number(expected_energy)),
                        ("relative_error", number(energy_error)),
                    ]),
                ),
                ("rigid_mode_relative_residuals", f64_array(rigid_residuals)),
            ]),
        ),
    ]))
}

/// Export PA3/PA5 geometry to the same production Q2 numbering and assembly
/// APIs, without attaching the legacy bracket-specific load and clamp.
/// Build a geometry-only Q2 package for PA3/PA5 sketches or bracket baselines.
/// The separate package schema leaves legacy solver packages unchanged.
pub fn export_sketch_analysis(model: &Json, order: usize) -> Result<Json, &'static str> {
    let schema = model.field("schema")?.text()?;
    if !matches!(
        schema,
        "sketch_layer_partition_v1"
            | "sketch_boolean_layers_v1"
            | "orthogonal_sketch_features_v1"
            | "filleted_bracket_v1"
            | "parameterized_bracket_v1"
    ) {
        return Err("unsupported model family for sketch analysis package");
    }
    if !(1..=16).contains(&order) {
        return Err("sketch analysis quadrature order must be in 1..=16");
    }
    let audit = interchange::audit(model)?;
    if audit.field("status")?.text()? != "valid" {
        return Err("sketch model is not certified valid");
    }
    let blocks = crate::bracket_json::decode_blocks(model)?;
    let incidence = crate::incidence::enumerate_faces(&blocks)?;
    let design_faces = sketch_design_faces(model)?;
    let space = analysis_space::from_certified_faces(&blocks, &incidence, &design_faces)?;
    let volume = integrate_volume(&space, order)?;
    let skeleton = encode_skeleton(&space)?;
    let audit_determinant = audit
        .field("determinant_lower_bound_mm3")
        .map_or(Json::Null, Clone::clone);
    let analysis = object([
        ("scalar_ndof", number(space.ndof)),
        (
            "vector_ndof",
            number(
                space
                    .ndof
                    .checked_mul(3)
                    .ok_or("vector DOF count overflow")?,
            ),
        ),
        ("ordering", text("scalar_node_major_xyz")),
        ("block_control_order", text("u-major,v-major,w-fast")),
        ("blocks", arr(&space.blocks, analysis_block_json)),
        (
            "boundary_dofs",
            Json::Object(
                space
                    .boundary_dofs
                    .iter()
                    .map(|(name, dofs)| (name.clone(), numeric_array(dofs)))
                    .collect(),
            ),
        ),
    ]);
    Ok(object([
        ("schema", text("solidcraft_sketch_analysis_package_v1")),
        ("source_schema", text(schema)),
        ("source_model", model.clone()),
        ("audit", audit),
        ("topology", model.field("topology")?.clone()),
        ("design_faces", model.field("design_faces")?.clone()),
        ("skeleton", skeleton),
        ("analysis_space", analysis),
        ("quadrature", gauss_json(order)?),
        (
            "measures",
            object([
                ("integrated_volume_mm3", number(volume)),
                ("audit_determinant_lower_bound_mm3", audit_determinant),
            ]),
        ),
        ("analysis_diagnostics", analysis_diagnostics(&space)?),
    ]))
}
