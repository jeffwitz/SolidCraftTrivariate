//! User-defined elasticity loads on compiled variable-size spline fields.
use crate::{
    constrained_elasticity::solve_constrained,
    constraint_system::Lift,
    elasticity::gauss,
    extracted_space::{assemble_sampled, ElementBasisSample, ElementGeometrySample},
    incidence::BlockSide,
    json::{object, Json},
    mechanical_job::{vector, Traction},
};
use std::collections::{BTreeMap, BTreeSet};
type Result<T> = std::result::Result<T, &'static str>;
type Face = (usize, usize, bool);
type Surfaces = BTreeMap<String, BTreeSet<Face>>;
fn n(v: impl ToString) -> Json {
    Json::Number(v.to_string())
}
fn string(v: &str) -> Json {
    Json::String(v.into())
}
fn select(surfaces: &Surfaces, name: &str) -> Result<Vec<Face>> {
    surfaces
        .get(name)
        .map(|v| v.iter().copied().collect())
        .ok_or("empty spline feature selection")
}

/// Both saved field types replay their exact coefficient maps before assembly.
pub fn solve(package: &Json, job: &Json) -> Result<Json> {
    solve_with_field(package, job, None)
}

/// Compile raw geometry once, then solve one feature-load job directly on the
/// in-memory automatic field. Persisted packages continue to use canonical replay.
pub fn solve_raw_inherited(model: &Json, job: &Json, degree: usize) -> Result<Json> {
    let compile_start = std::time::Instant::now();
    let plan = crate::inherited_space::automatic_plan(degree);
    let field = crate::inherited_space::space_from_plan(model, &plan)?;
    let field_compilation_ms = compile_start.elapsed().as_secs_f64() * 1000.;
    let package = object([
        ("schema", string("solidcraft_inherited_analysis_v1")),
        ("geometry", model.clone()),
        ("plan", plan),
        ("field", field.document()),
    ]);
    let mut result = solve_with_field(&package, job, Some(field))?;
    if let Json::Object(fields) = &mut result {
        fields.insert("field_compilation_ms".into(), n(field_compilation_ms));
    }
    Ok(result)
}

fn solve_with_field(
    package: &Json,
    job: &Json,
    mut prebuilt: Option<crate::inherited_space::InheritedSpace>,
) -> Result<Json> {
    if job.field("contact").is_ok() || job.field("contact_model").is_ok() {
        return Err("contact mechanics is unsupported by spline_feature_load_v1");
    }
    let clamp = job.field("clamp_feature")?.text()?;
    let mut surfaces = Surfaces::new();
    match package.field("schema")?.text()? {
        "local_thb_analysis_v1" => {
            let field = crate::thb_job::replay(package)?;
            let summary = package.field("field")?;
            for (&face, record) in field
                .boundary
                .iter()
                .zip(summary.field("boundary_sources")?.items()?)
            {
                let feature = record.field("design_face")?.field("feature")?.text()?;
                surfaces.entry(feature.into()).or_default().insert(face);
            }
            let lift = field.zero_face_lift(&select(&surfaces, clamp)?)?;
            solve_field(
                package,
                job,
                FieldDomain {
                    ndof: field.ndof,
                    cell_count: field.leaves.len(),
                    lift: &lift,
                    surfaces: &surfaces,
                    support_id: summary.field("support_id")?.clone(),
                },
                |i, p| field.sample(i, p),
            )
        }
        "solidcraft_inherited_analysis_v1" => {
            let field = match prebuilt.take() {
                Some(field) => field,
                None => crate::inherited_space::replay(package)?,
            };
            let model = package.field("geometry")?;
            let cells: BTreeMap<_, _> = field
                .cells
                .iter()
                .enumerate()
                .map(|(i, c)| (c.id.as_str(), i))
                .collect();
            let exterior: BTreeSet<_> =
                crate::incidence::enumerate_faces(&crate::bracket_json::decode_blocks(model)?)?
                    .exterior
                    .into_iter()
                    .collect();
            for record in model.field("design_faces")?.items()? {
                let block = record.field("block")?.text()?;
                let side = crate::mechanical_job::side(block, record.field("side")?.text()?)?;
                if !exterior.contains(&side) {
                    return Err("spline load selection is not an exterior face");
                }
                let idx = *cells.get(block).ok_or("unknown design-face block")?;
                surfaces
                    .entry(record.field("feature")?.text()?.into())
                    .or_default()
                    .insert((idx, side.axis, side.upper));
            }
            let fixed = select(&surfaces, clamp)?
                .into_iter()
                .map(|(i, axis, upper)| BlockSide {
                    block: field.cells[i].id.clone(),
                    axis,
                    upper,
                })
                .collect::<Vec<_>>();
            let lift = field.zero_face_lift(&fixed)?;
            solve_field(
                package,
                job,
                FieldDomain {
                    ndof: field.ndof,
                    cell_count: field.cells.len(),
                    lift: &lift,
                    surfaces: &surfaces,
                    support_id: Json::Null,
                },
                |i, p| field.cells[i].sample(p),
            )
        }
        _ => Err("unsupported compiled spline field"),
    }
}

struct FieldDomain<'a> {
    ndof: usize,
    cell_count: usize,
    lift: &'a Lift,
    surfaces: &'a Surfaces,
    support_id: Json,
}
fn solve_field<F>(package: &Json, job: &Json, domain: FieldDomain<'_>, sample: F) -> Result<Json>
where
    F: Fn(usize, [f64; 3]) -> Result<(ElementGeometrySample, ElementBasisSample)>,
{
    let FieldDomain {
        ndof,
        cell_count,
        lift,
        surfaces,
        support_id,
    } = domain;
    let select = |name| select(surfaces, name);
    if job.field("contact").is_ok() || job.field("contact_model").is_ok() {
        return Err("contact mechanics is unsupported by spline_feature_load_v1");
    }
    job.keys(&[
        "schema",
        "material",
        "clamp_feature",
        "tractions",
        "body_n_per_mm3",
        "quadrature",
        "relative_tolerance",
        "max_iterations",
    ])?;
    if job.field("schema")?.text()? != "spline_feature_load_v1" {
        return Err("unsupported local spline mechanical job");
    }
    let material = job.field("material")?;
    material.keys(&["young_mpa", "poisson"])?;
    let e = material.field("young_mpa")?.number::<f64>()?;
    let nu = material.field("poisson")?.number::<f64>()?;
    if !e.is_finite() || e <= 0. || !(-1. < nu && nu < 0.5) {
        return Err("invalid isotropic material");
    }
    let order = job.field("quadrature")?.number()?;
    let body = vector(job.field("body_n_per_mm3")?)?;
    let mu = e / (2. * (1. + nu));
    let lambda = 2. * mu * nu / (1. - 2. * nu);
    let mut assembly = assemble_sampled(ndof, cell_count, lambda, mu, order, &sample, |_| body)?;
    let q = gauss(order)?;
    for load in job.field("tractions")?.items()? {
        let selector = load.field("select")?;
        selector.keys(&["feature"])?;
        let traction = Traction::parse(load)?;
        for (i, axis, upper) in select(selector.field("feature")?.text()?)? {
            let (a, b) = ((axis + 1) % 3, (axis + 2) % 3);
            for &(u, wu) in &q {
                for &(v, wv) in &q {
                    let mut pt = [0.; 3];
                    pt[axis] = usize::from(upper) as f64;
                    pt[a] = u;
                    pt[b] = v;
                    let (g, basis) = sample(i, pt)?;
                    let (du, dv) = (g.jacobian_columns[a], g.jacobian_columns[b]);
                    let area: [f64; 3] = std::array::from_fn(|k| {
                        du[(k + 1) % 3] * dv[(k + 2) % 3] - du[(k + 2) % 3] * dv[(k + 1) % 3]
                    });
                    let norm = area.iter().map(|x| x * x).sum::<f64>().sqrt();
                    if !norm.is_finite() || norm <= 0. {
                        return Err("invalid spline load face");
                    }
                    let t = traction.value(area.map(|x| if upper { x / norm } else { -x / norm }));
                    for (&id, r) in basis.dofs.iter().zip(basis.values) {
                        for (c, &tc) in t.iter().enumerate() {
                            assembly.force[3 * id + c] += wu * wv * norm * r * tc;
                        }
                    }
                }
            }
        }
    }
    let tol = job.field("relative_tolerance")?.number()?;
    let limit = job.field("max_iterations")?.number()?;
    let solution = solve_constrained(&assembly.stiffness, [lift; 3], &assembly.force, tol, limit)?;
    let residual = assembly
        .stiffness
        .equilibrium_residual(&solution.displacement_mm, &assembly.force)?;
    let mut reactions = [0.; 3];
    for line in &lift.lines {
        let trace_weight = 1.
            - line
                .terms
                .iter()
                .map(|(_, x)| crate::inherited_space::val(*x))
                .sum::<f64>();
        for c in 0..3 {
            reactions[c] -= trace_weight * residual[3 * line.dof + c];
        }
    }
    let resultant = |v: &[f64]| -> Json {
        Json::Array(
            (0..3)
                .map(|c| n(v.iter().skip(c).step_by(3).sum::<f64>()))
                .collect(),
        )
    };
    let mut result = object([
        ("schema", string("spline_feature_result_v1")),
        ("support_id", support_id),
        ("plan", package.field("plan")?.clone()),
        ("job", job.clone()),
        ("scalar_dofs", n(ndof)),
        ("free_vector_dofs", n(3 * lift.free_roots.len())),
        ("volume_mm3", n(assembly.stiffness.integrated_volume)),
        ("force_n", resultant(&assembly.force)),
        (
            "reaction_n",
            Json::Array(reactions.into_iter().map(n).collect()),
        ),
        ("energy_n_mm", n(solution.strain_energy_n_mm)),
        ("relative_residual", n(solution.relative_residual)),
        ("iterations", n(solution.iterations)),
        (
            "coefficients_xyz_mm",
            Json::Array(
                solution
                    .displacement_mm
                    .as_chunks::<3>()
                    .0
                    .iter()
                    .map(|v| Json::Array(v.iter().map(n).collect()))
                    .collect(),
            ),
        ),
    ]);
    if let Ok(selection) = package
        .field("field")
        .and_then(|field| field.field("continuity_selection"))
    {
        if let Json::Object(fields) = &mut result {
            fields.insert("continuity_selection".into(), selection.clone());
            fields.insert(
                "interface_decisions".into(),
                package.field("field")?.field("interfaces")?.clone(),
            );
            fields.insert(
                "physics_scope".into(),
                Json::String("small_strain_linear_isotropic_elasticity; homogeneous young_mpa and poisson; constant body force and prescribed traction; no contact".into()),
            );
        }
    }
    Ok(result)
}
