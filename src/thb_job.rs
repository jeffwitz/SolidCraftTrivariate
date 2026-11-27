//! Compact, replayable definition of the local field, with source-face lineage.
//! Detailed exact extraction is an optional export. A solve always rebuilds
//! the field from its audited source and the declared refinement request.
use crate::{
    constraint_system::Ratio,
    extracted_space::{assemble_sampled, compare_assembled},
    json::{object, Json},
    section_manufactured::{solve_sampled, IntegrationDomain},
    spline_support::recover_nurbs_strip_parent,
    thb_strip::{build_thb_strip_for_cells, ThbStrip},
};
type Result<T> = std::result::Result<T, &'static str>;
fn n(v: impl ToString) -> Json {
    Json::Number(v.to_string())
}
fn string(s: impl Into<String>) -> Json {
    Json::String(s.into())
}

fn mount(model: &Json, plan: &Json) -> Result<ThbStrip> {
    plan.keys(&["schema", "start", "end"])?;
    if plan.field("schema")?.text()? != "interior_thb_refinement_v1" {
        return Err("unsupported hierarchical refinement plan");
    }
    let index = |key: &str| -> Result<[usize; 3]> {
        let [x, y, z] = plan.field(key)?.items()? else {
            return Err("refinement bounds need three cell indices");
        };
        Ok([x.number()?, y.number()?, z.number()?])
    };
    build_thb_strip_for_cells(model, index("start")?, index("end")?)
}
fn summary(model: &Json, field: &ThbStrip) -> Result<Json> {
    let support = recover_nurbs_strip_parent(model)?;
    let mut sources = Vec::new();
    for leaf in &field.leaves {
        let mut found = None;
        for span in &support.source_spans {
            let mut inside = true;
            for (axis, [a, b]) in span.restriction.intervals.into_iter().enumerate() {
                inside &= !leaf.intervals[axis][0]
                    .cmp_value(&Ratio::new(i128::from(a), 1)?)
                    .is_lt()
                    && !leaf.intervals[axis][1]
                        .cmp_value(&Ratio::new(i128::from(b), 1)?)
                        .is_gt();
            }
            if inside && found.replace(span).is_some() {
                return Err("ambiguous integration leaf lineage");
            }
        }
        sources.push(found.ok_or("leaf crosses a construction span")?);
    }
    let mut boundary = Vec::new();
    let mut covered_faces = std::collections::BTreeSet::new();
    for &(i, axis, upper) in &field.boundary {
        let span = sources[i];
        let local = span
            .restriction
            .axes
            .iter()
            .position(|&a| a == axis)
            .ok_or("missing leaf boundary axis")?;
        let side = format!(
            "{}{}",
            ["u", "v", "w"][local],
            usize::from(upper ^ span.restriction.reversals[local])
        );
        let design = model
            .field("design_faces")?
            .items()?
            .iter()
            .find(|f| {
                f.field("block").and_then(Json::text) == Ok(span.source_block_id.as_str())
                    && f.field("side").and_then(Json::text) == Ok(side.as_str())
            })
            .ok_or("missing source design face")?;
        covered_faces.insert((span.source_block_id.clone(), side));
        boundary.push(object([
            ("leaf", string(field.leaves[i].id.clone())),
            ("leaf_axis", n(axis)),
            ("upper", Json::Bool(upper)),
            ("design_face", design.clone()),
        ]));
    }
    if covered_faces.len() != model.field("design_faces")?.items()?.len() {
        return Err("local refinement lost a source design face");
    }
    let fine = field.leaves.iter().filter(|c| c.level == 1).count();
    Ok(object([
        ("profile", string("interior_dyadic_thb_cubic_v1")),
        ("support_id", string(support.id)),
        (
            "parent_geometry_controls",
            n(support.homogeneous_controls.len()),
        ),
        ("initial_scalar_dofs", n(field.initial_ndof)),
        ("prepared_scalar_dofs", n(field.prepared_ndof)),
        ("local_scalar_dofs", n(field.ndof)),
        ("uniform_fine_scalar_dofs", n(field.uniform_fine_ndof)),
        ("prepared_transfer_rank", n(field.transfer_rank)),
        (
            "local_basis_rank_in_uniform_fine_space",
            n(field.uniform_fine_rank),
        ),
        ("fine_active", n(field.fine_active_count)),
        ("coarse_active", n(field.coarse_active_count)),
        ("coarse_leaves", n(field.leaves.len() - fine)),
        ("fine_leaves", n(fine)),
        (
            "common_partition_verified",
            Json::Bool(field.common_partition_verified),
        ),
        (
            "integration_sources",
            Json::Array(
                field
                    .leaves
                    .iter()
                    .zip(sources)
                    .map(|(l, s)| {
                        Json::Array(vec![
                            string(l.id.clone()),
                            string(s.source_block_id.clone()),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("boundary_sources", Json::Array(boundary)),
    ]))
}
/// A field refusal never changes the validity status of the source geometry.
pub fn refusal(model: &Json, reason: &str) -> Json {
    let validity = crate::interchange::audit(model)
        .ok()
        .and_then(|r| r.field("status").ok().cloned())
        .unwrap_or_else(|| string("unverified"));
    let outside = model.field("schema").and_then(Json::text) != Ok("nurbs_sketch_strip_v1");
    let status = if reason.contains("overflow")
        || reason.contains("budget")
        || reason.contains("uncertified")
    {
        "indeterminate"
    } else if outside
        || reason.starts_with("THB marked")
        || reason.contains("requires a cubic")
        || reason.contains("nonuniform")
        || reason.contains("repeated source")
    {
        "unsupported_profile"
    } else {
        "error"
    };
    object([
        ("status", string(status)),
        ("source_validity", validity),
        ("field_eligibility", string(status)),
        ("reason", string(reason)),
    ])
}

pub fn compile(model: &Json, plan: &Json) -> Result<Json> {
    let field = mount(model, plan)?;
    Ok(object([
        ("schema", string("local_thb_analysis_v1")),
        ("geometry", model.clone()),
        ("plan", plan.clone()),
        ("field", summary(model, &field)?),
    ]))
}
pub fn replay(package: &Json) -> Result<ThbStrip> {
    package.keys(&["schema", "geometry", "plan", "field"])?;
    if package.field("schema")?.text()? != "local_thb_analysis_v1" {
        return Err("invalid local THB package schema");
    }
    let geometry = package.field("geometry")?;
    let field = mount(geometry, package.field("plan")?)?;
    if package.field("field")? != &summary(geometry, &field)? {
        return Err("hierarchical field or source-face lineage differs from replay");
    }
    Ok(field)
}
pub fn solve(package: &Json, mode: &str, job: &Json) -> Result<Json> {
    let field = replay(package)?;
    let (count, lift) = match mode {
        "prepared" => (field.prepared_ndof, field.prepared_zero_bottom_lift()?),
        "local" => {
            let bottom: Vec<_> = field
                .boundary
                .iter()
                .filter(|(_, axis, upper)| *axis == 2 && !upper)
                .copied()
                .collect();
            (field.ndof, field.zero_face_lift(&bottom)?)
        }
        _ => return Err("expected prepared or local THB space"),
    };
    let domain = IntegrationDomain {
        scalar_dofs: count,
        cell_count: field.leaves.len(),
        z_bounds: field.z_bounds,
        boundary: &field.boundary,
        essential: &lift,
    };
    let mut result = solve_sampled(&domain, job, |i, p| {
        if mode == "local" {
            field.sample(i, p)
        } else {
            field.sample_prepared(i, p)
        }
    })?;
    if let Json::Object(fields) = &mut result {
        fields.insert("field_space".into(), string(mode));
    }
    Ok(result)
}
/// Restrict an assembled block matrix and force to unconstrained scalar roots
/// of a homogeneous, whole-feature displacement clamp.
fn principal_system(
    system: &crate::extracted_space::NativeAssembly,
    free: &[usize],
) -> crate::extracted_space::NativeAssembly {
    use crate::{elasticity::Stiffness, extracted_space::NativeAssembly};
    let k = &system.stiffness;
    let mut offsets = vec![0];
    let mut columns = Vec::new();
    let mut values = Vec::new();
    let mut force = Vec::new();
    for &i in free {
        for j in k.offsets[i]..k.offsets[i + 1] {
            if let Ok(col) = free.binary_search(&k.columns[j]) {
                columns.push(col);
                values.push(k.values[j]);
            }
        }
        offsets.push(columns.len());
        force.extend_from_slice(&system.force[3 * i..3 * i + 3]);
    }
    NativeAssembly {
        stiffness: Stiffness {
            offsets,
            columns,
            values,
            integrated_volume: k.integrated_volume,
        },
        force,
        quadrature_points: system.quadrature_points,
        quadrature_order: system.quadrature_order,
        cell_count: system.cell_count,
        quadrature_signature: system.quadrature_signature,
    }
}

pub fn galerkin(package: &Json, order: usize) -> Result<Json> {
    let field = replay(package)?;
    let load = |x: [f64; 3]| [1. + x[0], x[1] - x[2], x[0] * x[1] + 2.];
    let fine = assemble_sampled(
        field.ndof,
        field.leaves.len(),
        1.5,
        1.,
        order,
        |i, p| field.sample(i, p),
        load,
    )?;
    let coarse = assemble_sampled(
        field.prepared_ndof,
        field.leaves.len(),
        1.5,
        1.,
        order,
        |i, p| field.sample_prepared(i, p),
        load,
    )?;
    let cmp = compare_assembled(
        &fine,
        &coarse,
        field.ndof,
        field.prepared_ndof,
        &field.rational_prepared_to_thb,
    )?;
    let bottom: Vec<_> = field
        .boundary
        .iter()
        .copied()
        .filter(|(_, a, upper)| *a == 2 && !upper)
        .collect();
    let coarse_roots = field.prepared_zero_bottom_lift()?.free_roots;
    let fine_roots = field.zero_face_lift(&bottom)?.free_roots;
    let mut reduced_transfer = Vec::new();
    for (i, row) in field.rational_prepared_to_thb.iter().enumerate() {
        let reduced: std::collections::BTreeMap<_, _> = row
            .iter()
            .filter_map(|(&id, &x)| coarse_roots.binary_search(&id).ok().map(|j| (j, x)))
            .collect();
        if fine_roots.binary_search(&i).is_ok() {
            reduced_transfer.push(reduced);
        } else if !reduced.is_empty() {
            return Err("transfer violates homogeneous source-face clamp");
        }
    }
    let reduced = compare_assembled(
        &principal_system(&fine, &fine_roots),
        &principal_system(&coarse, &coarse_roots),
        fine_roots.len(),
        coarse_roots.len(),
        &reduced_transfer,
    )?;
    for check in [&cmp, &reduced] {
        if check.stiffness_relative_frobenius_defect > 1e-10
            || check.force_relative_l2_defect > 1e-10
        {
            return Err("THB Galerkin pullback exceeds tolerance");
        }
    }
    Ok(object([
        ("status", string("valid")),
        ("clamped_transfer_commutes_exactly", Json::Bool(true)),
        (
            "clamped_relative_stiffness_defect",
            n(reduced.stiffness_relative_frobenius_defect),
        ),
        (
            "clamped_relative_load_defect",
            n(reduced.force_relative_l2_defect),
        ),
        ("quadrature", n(order)),
        ("quadrature_points", n(fine.quadrature_points)),
        (
            "relative_stiffness_defect",
            n(cmp.stiffness_relative_frobenius_defect),
        ),
        ("relative_load_defect", n(cmp.force_relative_l2_defect)),
    ]))
}
