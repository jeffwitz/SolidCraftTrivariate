//! Solve a versioned, linear mechanical job against an independently audited native model.
use crate::{
    adaptive_package::{lift_json, parse_refinement_plan},
    adaptive_partition::AdaptivePartition,
    adaptive_space::build_adaptive_space,
    analysis_space::{from_certified_faces, AnalysisBlock, AnalysisSpace},
    boundary_traction::{integrate_traction, LoadPatch},
    bracket_json::{decode_blocks, side_json},
    constrained_elasticity::solve_constrained,
    constraint_system::{resolve, Ratio},
    elasticity::{assemble, Vector},
    incidence::{enumerate_faces, BlockSide},
    interchange,
    json::{object, Json},
    prescribed_motion::Motion,
    projected_bulk::assemble_projected_bulk,
    solver_package::analysis_block_json,
};
use std::collections::{BTreeMap, BTreeSet};

fn n(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn text(x: &str) -> Json {
    Json::String(x.into())
}
fn numbers<T: ToString>(v: impl IntoIterator<Item = T>) -> Json {
    Json::Array(v.into_iter().map(n).collect())
}
fn finite(x: &Json) -> Result<f64, &'static str> {
    let v = x.number::<f64>()?;
    if v.is_finite() {
        Ok(v)
    } else {
        Err("nonfinite mechanical job scalar")
    }
}
pub(crate) fn vector(x: &Json) -> Result<Vector, &'static str> {
    let [a, b, c] = x.items()? else {
        return Err("expected mechanical vector of length three");
    };
    Ok([finite(a)?, finite(b)?, finite(c)?])
}
fn side(block: &str, s: &str) -> Result<BlockSide, &'static str> {
    let (axis, upper) = match s {
        "u0" => (0, false),
        "u1" => (0, true),
        "v0" => (1, false),
        "v1" => (1, true),
        "w0" => (2, false),
        "w1" => (2, true),
        _ => return Err("invalid mechanical face side"),
    };
    Ok(BlockSide {
        block: block.into(),
        axis,
        upper,
    })
}
fn face_indices(side: &BlockSide) -> impl Iterator<Item = usize> + '_ {
    (0..27).filter(move |&i| [i / 9, (i / 3) % 3, i % 3][side.axis] == 2 * usize::from(side.upper))
}

struct Surface<'a> {
    side: BlockSide,
    feature: &'a str,
    block: &'a AnalysisBlock,
}
fn select<'a>(
    input: &Json,
    faces: &'a [Surface<'a>],
) -> Result<Vec<&'a Surface<'a>>, &'static str> {
    let feature = input.field("feature")?.text()?;
    let plane = if let Ok(p) = input.field("plane") {
        input.keys(&["feature", "plane"])?;
        p.keys(&["axis", "coordinate_mm"])?;
        let axis = p.field("axis")?.number::<usize>()?;
        if axis >= 3 {
            return Err("invalid selector plane axis");
        }
        Some((axis, i128::from(p.field("coordinate_mm")?.number::<i64>()?)))
    } else {
        input.keys(&["feature"])?;
        None
    };
    let mut selection = Vec::new();
    for f in faces {
        if feature != "*" && f.feature != feature {
            continue;
        }
        if let Some((axis, offset)) = plane {
            if !face_indices(&f.side).all(|i| {
                let h = f.block.controls[i];
                h[3].checked_mul(offset) == Some(h[axis])
            }) {
                continue;
            }
        }
        selection.push(f);
    }
    if selection.is_empty() {
        return Err("empty mechanical face selection");
    }
    Ok(selection)
}

enum Traction {
    Vector(Vector),
    Pressure(f64),
    Stress([Vector; 3]),
}
impl Traction {
    fn parse(x: &Json) -> Result<Self, &'static str> {
        if let Ok(v) = x.field("vector_mpa") {
            x.keys(&["select", "vector_mpa"])?;
            return Ok(Self::Vector(vector(v)?));
        }
        if let Ok(v) = x.field("pressure_mpa") {
            x.keys(&["select", "pressure_mpa"])?;
            let p = finite(v)?;
            if p < 0. {
                return Err("negative compressive pressure");
            }
            return Ok(Self::Pressure(p));
        }
        x.keys(&["select", "stress_mpa"])?;
        let [a, b, c] = x.field("stress_mpa")?.items()? else {
            return Err("expected stress tensor");
        };
        let s = [vector(a)?, vector(b)?, vector(c)?];
        if (0..3).any(|i| (0..i).any(|j| s[i][j] != s[j][i])) {
            return Err("nonsymmetric prescribed stress tensor");
        }
        Ok(Self::Stress(s))
    }
    fn value(&self, normal: Vector) -> Vector {
        match self {
            Self::Vector(v) => *v,
            Self::Pressure(p) => normal.map(|v| -p * v),
            Self::Stress(s) => s.map(|row| dot(row, normal)),
        }
    }
    fn magnitude_bound(&self) -> f64 {
        match self {
            Self::Vector(x) => norm(*x),
            Self::Pressure(p) => *p,
            Self::Stress(s) => s.iter().map(|r| dot(*r, *r)).sum::<f64>().sqrt(),
        }
    }
}
fn dot(x: Vector, y: Vector) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}
fn norm(x: Vector) -> f64 {
    dot(x, x).sqrt()
}
fn cross(x: Vector, y: Vector) -> Vector {
    std::array::from_fn(|i| x[(i + 1) % 3] * y[(i + 2) % 3] - x[(i + 2) % 3] * y[(i + 1) % 3])
}
fn add(a: &mut Vector, b: Vector) {
    for (v, w) in a.iter_mut().zip(b) {
        *v += w;
    }
}
fn homogeneous_nodes(space: &AnalysisSpace) -> Vec<[i128; 4]> {
    let mut points = vec![[0; 4]; space.ndof];
    for b in &space.blocks {
        for (&id, h) in b.dofs.iter().zip(&b.controls) {
            points[id] = *h;
        }
    }
    points
}

/// Rank of the evaluation of physical rigid fields on prescribed coefficients.
/// Centering and uniform rescaling keep translations and rotations comparable.
fn mounting_rank(points: &[Vector], fixed: &[BTreeMap<usize, Ratio>; 3]) -> usize {
    let lo: Vector =
        std::array::from_fn(|j| points.iter().map(|x| x[j]).fold(f64::INFINITY, f64::min));
    let hi: Vector = std::array::from_fn(|j| {
        points
            .iter()
            .map(|x| x[j])
            .fold(f64::NEG_INFINITY, f64::max)
    });
    let size = (0..3).map(|j| hi[j] - lo[j]).fold(0., f64::max);
    if size <= 0. {
        return 0;
    }
    let mut pivots = [None::<[f64; 6]>; 6];
    for (c, ids) in fixed.iter().enumerate() {
        for &id in ids.keys() {
            let point = std::array::from_fn(|j| (points[id][j] - (lo[j] + hi[j]) * 0.5) / size);
            let unit: Vector = std::array::from_fn(|j| f64::from(c == j));
            let rotation = cross(point, unit);
            let mut row = [
                unit[0],
                unit[1],
                unit[2],
                rotation[0],
                rotation[1],
                rotation[2],
            ];
            for j in 0..6 {
                if let Some(b) = pivots[j] {
                    let factor = row[j];
                    for k in j..6 {
                        row[k] -= factor * b[k];
                    }
                } else if row[j].abs() > 1e-11 {
                    let scale = row[j];
                    for x in &mut row[j..] {
                        *x /= scale;
                    }
                    pivots[j] = Some(row);
                    break;
                }
            }
        }
    }
    pivots.iter().flatten().count()
}

/// Build the effective space, apply essential conditions, solve,
/// and return a reproducible result document. Publication is the caller's job.
pub fn solve_job(model: &Json, job: &Json) -> Result<Json, &'static str> {
    run_job(model, job, None)
}

/// Differentiate a fixed material selection and homogeneous mounting with respect
/// to a global affine deformation of the audited configuration.
pub fn differentiate_job(model: &Json, job: &Json, direction: &Json) -> Result<Json, &'static str> {
    let velocity = crate::affine_job_sensitivity::Velocity::parse(direction)?;
    if job.field("formulation")?.text()? != "displacement" {
        return Err("shape derivative requires displacement formulation");
    }
    for group in ["loads", "supports"] {
        for entry in job.field(group)?.items()? {
            if entry.field("select")?.field("plane").is_ok() {
                return Err("spatial plane selection requires event analysis");
            }
            if entry.field("affine_mm").is_ok() && velocity.exact.is_none() {
                return Err("prescribed motion derivative unavailable");
            }
        }
    }
    run_job(model, job, Some(&velocity))
}

fn run_job(
    model: &Json,
    job: &Json,
    velocity: Option<&crate::affine_job_sensitivity::Velocity>,
) -> Result<Json, &'static str> {
    let version = job.field("schema")?.text()?;
    let mut fields = vec![
        "schema",
        "material",
        "formulation",
        "order",
        "plan",
        "supports",
        "loads",
    ];
    let adaptive_surface = version == "solidcraft_mechanical_job_v5";
    let supports_body = version == "solidcraft_mechanical_job_v4"
        || (adaptive_surface && job.field("body_load").is_ok());
    let mut quadrature = if adaptive_surface {
        fields.push("surface_quadrature");
        let q = job.field("surface_quadrature")?;
        use crate::adaptive_traction::SurfaceSplit;
        let mut keys = vec!["relative_tolerance", "max_depth", "max_panels"];
        let freeze_primal = if let Ok(mode) = q.field("derivative") {
            keys.push("derivative");
            if mode.text()? != "freeze_primal" {
                return Err("unknown surface derivative plan");
            }
            true
        } else {
            false
        };
        if velocity.is_some() && !freeze_primal {
            return Err("shape derivative requires a fixed surface quadrature plan");
        }
        let strategy = if let Ok(split) = q.field("split") {
            keys.push("split");
            match split.text()? {
                "directional" => SurfaceSplit::Directional,
                "isotropic" => SurfaceSplit::Isotropic,
                _ => return Err("unsupported adaptive surface split"),
            }
        } else {
            SurfaceSplit::Isotropic
        };
        q.keys(&keys)?;
        let config = crate::adaptive_traction::SurfaceQuadrature {
            relative_tolerance: finite(q.field("relative_tolerance")?)?,
            max_depth: q.field("max_depth")?.number()?,
            max_panels: q.field("max_panels")?.number()?,
        };
        config.validate()?;
        Some((config, strategy, freeze_primal))
    } else {
        None
    };
    if supports_body {
        fields.push("body_load");
    }
    job.keys(&fields)?;
    let density = if supports_body {
        Some(crate::volume_load::Density::parse(job.field("body_load")?)?)
    } else {
        None
    };
    let eliminate_supports = matches!(
        version,
        "solidcraft_mechanical_job_v3"
            | "solidcraft_mechanical_job_v4"
            | "solidcraft_mechanical_job_v5"
    );
    let version_two = match version {
        "solidcraft_mechanical_job_v1" => false,
        "solidcraft_mechanical_job_v2"
        | "solidcraft_mechanical_job_v3"
        | "solidcraft_mechanical_job_v4"
        | "solidcraft_mechanical_job_v5" => true,
        _ => return Err("unknown mechanical job schema"),
    };
    let audit = interchange::audit(model)?;
    if audit.field("status")?.text()? != "valid"
        || !matches!(
            audit.field("method")?.text()?,
            "certified_fixed_extrusion_partition"
                | "certified_ruled_similarity_sweep"
                | "certified_coaxial_feature_program"
                | "certified_circular_bend"
                | "certified_orthogonal_feature_program"
                | "certified_polygon_prism"
                | "certified_perforated_polygon_prism"
                | crate::interval_bores_json::METHOD
                | crate::circular_polygon_json::METHOD
                | crate::affine_image::METHOD
                | "certified_parameterized_bracket"
        )
    {
        return Err("unsupported or unaudited mechanical model");
    }
    let material = job.field("material")?;
    material.keys(&["young_mpa", "poisson"])?;
    let young = finite(material.field("young_mpa")?)?;
    let poisson = finite(material.field("poisson")?)?;
    let order = job.field("order")?.number::<usize>()?;
    if !(2..=32).contains(&order) {
        return Err("invalid mechanical quadrature order");
    }
    let plan = parse_refinement_plan(job.field("plan")?)?;
    let blocks = decode_blocks(model)?;
    let incidence = enumerate_faces(&blocks)?;
    let mut neighbors: BTreeMap<&str, Vec<&str>> =
        blocks.iter().map(|b| (b.id.as_str(), Vec::new())).collect();
    for interface in &incidence.interfaces {
        let a = interface.first.block.as_str();
        let b = interface.second.block.as_str();
        neighbors
            .get_mut(a)
            .ok_or("unknown interface block")?
            .push(b);
        neighbors
            .get_mut(b)
            .ok_or("unknown interface block")?
            .push(a);
    }
    let mut pending = vec![*neighbors.keys().next().ok_or("empty mechanical model")?];
    let mut visited = BTreeSet::new();
    while let Some(id) = pending.pop() {
        if visited.insert(id) {
            pending.extend(neighbors[id].iter().copied());
        }
    }
    if visited.len() != blocks.len() {
        return Err("mechanical job requires one face-connected material domain");
    }
    let mut groups = BTreeMap::new();
    for face in model.field("design_faces")?.items()? {
        let s = side(face.field("block")?.text()?, face.field("side")?.text()?)?;
        if groups
            .insert(s, face.field("feature")?.text()?.to_owned())
            .is_some()
        {
            return Err("duplicate mechanical design face");
        }
    }
    if incidence.exterior.iter().collect::<BTreeSet<_>>() != groups.keys().collect() {
        return Err("incomplete mechanical exterior provenance");
    }
    let design = groups
        .iter()
        .map(|(s, t)| (s.clone(), t.clone()))
        .collect::<Vec<_>>();
    let parent = from_certified_faces(&blocks, &incidence, &design)?;
    let roots = parent
        .blocks
        .iter()
        .map(|b| b.id.clone())
        .collect::<Vec<_>>();
    let mut partition = AdaptivePartition::new(&roots, &incidence, plan.max_depth, plan.max_cells)?;
    for step in &plan.steps {
        partition.refine(step)?;
    }
    let adaptive = build_adaptive_space(&parent, &partition)?;
    let space = &adaptive.space;
    let lookup: BTreeMap<_, _> = space.blocks.iter().map(|b| (b.id.as_str(), b)).collect();
    let mut faces = Vec::new();
    for p in partition.exterior_patches() {
        let root = p.side.root_side;
        let id = &adaptive.leaf_blocks[&p.side.leaf];
        faces.push(Surface {
            side: BlockSide {
                block: id.clone(),
                axis: root.axis,
                upper: root.upper,
            },
            feature: &groups[&root],
            block: lookup[id.as_str()],
        });
    }
    let mut fixed: [BTreeMap<usize, Ratio>; 3] = Default::default();
    let exact_velocity = velocity.and_then(|v| v.exact);
    let mut fixed_rates: [BTreeMap<usize, Ratio>; 3] = Default::default();
    let mut frame_motion = Motion::zero();
    let mut frame_rows = [false; 3];
    let hom = homogeneous_nodes(space);
    let points = hom
        .iter()
        .map(|h| std::array::from_fn(|c| h[c] as f64 / h[3] as f64))
        .collect::<Vec<Vector>>();
    let mut supports = Vec::new();
    for support in job.field("supports")?.items()? {
        let explicit = support.field("affine_mm").ok();
        let motion = if version_two && explicit.is_some() {
            support.keys(&["select", "components", "affine_mm"])?;
            Motion::parse(support.field("affine_mm")?)?
        } else {
            support.keys(&["select", "components"])?;
            Motion::zero()
        };
        let components = support
            .field("components")?
            .items()?
            .iter()
            .map(Json::number::<usize>)
            .collect::<Result<Vec<_>, _>>()?;
        if components.is_empty()
            || components.iter().any(|&c| c >= 3)
            || components.iter().collect::<BTreeSet<_>>().len() != components.len()
        {
            return Err("invalid mechanical support components");
        }
        if explicit.is_some() {
            for &c in &components {
                if !frame_rows[c] {
                    frame_motion.set_component(c, motion);
                    frame_rows[c] = true;
                }
            }
        }
        let selected = select(support.field("select")?, &faces)?;
        let support_rate = exact_velocity.map(|v| motion.along(v)).transpose()?;
        for f in &selected {
            for &c in &components {
                for i in face_indices(&f.side) {
                    let id = f.block.dofs[i];
                    let value = motion.value(c, hom[id])?;
                    if let Some(rate) = support_rate {
                        let derivative = rate.value(c, hom[id])?;
                        if fixed_rates[c]
                            .insert(id, derivative)
                            .is_some_and(|old| old != derivative)
                        {
                            return Err("incompatible material derivatives on shared support");
                        }
                    }
                    if fixed[c].insert(id, value).is_some_and(|old| old != value) {
                        return Err("conflicting displacement on shared mechanical support");
                    }
                }
            }
        }
        supports.push(object([
            (
                "faces",
                Json::Array(selected.iter().map(|f| side_json(&f.side)).collect()),
            ),
            ("components", numbers(components)),
        ]));
    }
    let rank = mounting_rank(&points, &fixed);
    if rank != 6 {
        return Err("mechanical mounting leaves a rigid mode");
    }
    let make_lifts = |values_by_component: &[BTreeMap<usize, Ratio>; 3],
                      frame: &Motion|
     -> Result<_, &'static str> {
        let mut lifts = Vec::with_capacity(3);
        for (c, values) in values_by_component.iter().enumerate() {
            let prescribed = values
                .iter()
                .map(|(&id, &value)| {
                    let offset = frame.value(c, hom[id])?.checked_mul(Ratio::new(-1, 1)?)?;
                    Ok((id, value.checked_add(offset)?))
                })
                .collect::<Result<Vec<_>, &'static str>>()?;
            lifts.push(if eliminate_supports {
                crate::constraint_system::impose_on_lift(&adaptive.lift, &prescribed)?
            } else {
                resolve(space.ndof, &adaptive.constraints, &prescribed)?
            });
        }
        Ok(lifts)
    };
    let rigid_frame = frame_rows
        .iter()
        .any(|v| *v)
        .then(|| frame_motion.rigid_part())
        .transpose()?;
    let physical_lifts = make_lifts(&fixed, &Motion::zero())?;
    let relative_lifts = rigid_frame
        .as_ref()
        .map(|frame| make_lifts(&fixed, frame))
        .transpose()?;
    let lifts = relative_lifts.as_ref().unwrap_or(&physical_lifts);
    let exterior = faces
        .iter()
        .map(|f| f.side.clone())
        .collect::<BTreeSet<_>>();
    let mut seen = BTreeSet::new();
    let mut force = vec![0.; 3 * space.ndof];
    let mut applied = [0.; 3];
    let mut moment = [0.; 3];
    let mut scale = 0.;
    let mut loads = Vec::new();
    let mut force_rate = velocity.map(|_| vec![0.; 3 * space.ndof]);
    for specification in job.field("loads")?.items()? {
        let law = Traction::parse(specification)?;
        let selected = select(specification.field("select")?, &faces)?;
        let mut patches = Vec::new();
        for f in selected {
            if !seen.insert(f.side.clone()) {
                return Err("overlapping mechanical loads");
            }
            patches.push(LoadPatch::whole(f.side.clone()));
        }
        let rate_law = |velocity: &crate::affine_job_sensitivity::Velocity, normal| {
            let oriented = velocity.oriented_area_rate(normal);
            match &law {
                Traction::Vector(v) => v.map(|value| value * dot(normal, oriented)),
                Traction::Pressure(_) | Traction::Stress(_) => law.value(oriented),
            }
        };
        let (load, integration) = if let Some((config, strategy, freeze_primal)) = &mut quadrature {
            let out = crate::adaptive_traction::integrate_surface_traction(
                space,
                &exterior,
                &patches,
                order,
                (*config, *strategy),
                |_, normal| law.value(normal),
            )?;
            config.max_panels = config
                .max_panels
                .checked_sub(out.panels)
                .ok_or("adaptive surface quadrature panel budget exceeded")?;
            let mut info = object([
                ("panels", n(out.panels)),
                ("deepest", n(out.deepest)),
                (
                    "maximum_estimated_relative_error",
                    n(out.maximum_estimated_relative_error),
                ),
            ]);
            if *strategy == crate::adaptive_traction::SurfaceSplit::Directional || *freeze_primal {
                let plans = patches
                    .iter()
                    .zip(&out.patch_depths)
                    .map(|(patch, depths)| {
                        object([
                            ("block", Json::String(patch.side.block.clone())),
                            ("axis", n(patch.side.axis)),
                            ("upper", Json::Bool(patch.side.upper)),
                            ("depths", Json::Array(depths.iter().map(n).collect())),
                        ])
                    })
                    .collect();
                if let Json::Object(fields) = &mut info {
                    let split = match strategy {
                        crate::adaptive_traction::SurfaceSplit::Directional => "directional",
                        crate::adaptive_traction::SurfaceSplit::Isotropic => "isotropic",
                    };
                    fields.insert("split".into(), text(split));
                    fields.insert("patch_plans".into(), Json::Array(plans));
                }
            }
            if let (Some(velocity), Some(rhs)) = (velocity, &mut force_rate) {
                let (rate, count) = crate::adaptive_traction::replay_traction_rate(
                    space,
                    &patches,
                    &out.patch_depths,
                    order,
                    config.max_panels,
                    |_, normal| rate_law(velocity, normal),
                )?;
                config.max_panels -= count;
                for (a, b) in rhs.iter_mut().zip(rate) {
                    *a += b;
                }
                if let Json::Object(fields) = &mut info {
                    fields.insert("replay_panels".into(), n(count));
                }
            }
            (out.load, Some(info))
        } else {
            if let (Some(velocity), Some(rhs)) = (velocity, &mut force_rate) {
                let rate = integrate_traction(space, &exterior, &patches, order, |_, normal| {
                    rate_law(velocity, normal)
                })?;
                for (a, b) in rhs.iter_mut().zip(rate.nodal_force_n) {
                    *a += b;
                }
            }
            (
                integrate_traction(space, &exterior, &patches, order, |_, normal| {
                    law.value(normal)
                })?,
                None,
            )
        };
        add(&mut applied, load.resultant_n);
        add(&mut moment, load.moment_n_mm);
        scale += load.area_mm2 * law.magnitude_bound();
        for (f, v) in force.iter_mut().zip(load.nodal_force_n) {
            *f += v;
        }
        let mut record = object([
            (
                "faces",
                Json::Array(patches.iter().map(|p| side_json(&p.side)).collect()),
            ),
            ("area_mm2", n(load.area_mm2)),
            ("resultant_n", numbers(load.resultant_n)),
            ("moment_n_mm", numbers(load.moment_n_mm)),
        ]);
        if let (Some(info), Json::Object(fields)) = (integration, &mut record) {
            fields.insert("integration".into(), info);
        }
        loads.push(record);
    }
    let body = density
        .map(|d| d.integrate(space, order, velocity))
        .transpose()?;
    let surface_force = body.as_ref().map(|_| force.clone());
    if let Some(load) = &body {
        add(&mut applied, load.resultant);
        add(&mut moment, load.moment);
        scale += load.scale;
        for (total, f) in force.iter_mut().zip(&load.force) {
            *total += f;
        }
        if let (Some(total), Some(rate)) = (&mut force_rate, &load.force_rate) {
            for (f, df) in total.iter_mut().zip(rate) {
                *f += df;
            }
        }
    }
    let moving = fixed
        .iter()
        .any(|prescribed| prescribed.values().any(|&v| v != Ratio::zero()));
    if !scale.is_finite() || (scale <= 0. && !moving) {
        return Err("mechanical job needs a finite nonzero load scale");
    }
    if let (Some(velocity), Some(rate)) = (velocity, force_rate) {
        use crate::{
            affine_job_sensitivity::{self, SensitivityJob},
            elasticity_sensitivity::MaterialJet,
        };
        let support_rates = if velocity.exact.is_some() {
            let derivative_lifts = make_lifts(&fixed_rates, &Motion::zero())?;
            let mut rates: [Vec<crate::mechanical_sensitivity::LiftRate>; 3] = Default::default();
            for c in 0..3 {
                let derivative = &derivative_lifts[c];
                let primal = &physical_lifts[c];
                if derivative.free_roots != primal.free_roots {
                    return Err("material support derivative changes free roots");
                }
                for (a, b) in derivative.lines.iter().zip(&primal.lines) {
                    if a.dof != b.dof || a.terms != b.terms {
                        return Err("material support derivative changes trace weights");
                    }
                    rates[c].push(crate::mechanical_sensitivity::LiftRate {
                        terms: Vec::new(),
                        offset: a.offset.numerator() as f64 / a.offset.denominator() as f64,
                    });
                }
            }
            Some(rates)
        } else {
            None
        };
        let mut result = affine_job_sensitivity::solve(
            velocity,
            SensitivityJob {
                space,
                force: &force,
                force_rate: &rate,
                includes_body_load: supports_body,
                support_rates,
                lifts: [&physical_lifts[0], &physical_lifts[1], &physical_lifts[2]],
                material: MaterialJet {
                    young_mpa: young,
                    poisson,
                    d_young_mpa: 0.,
                    d_poisson: 0.,
                },
                order,
            },
        )?;
        if adaptive_surface {
            if let Json::Object(fields) = &mut result {
                fields.insert("surface_derivative_plan".into(), text("freeze_primal"));
                fields.insert("surface_loads".into(), Json::Array(loads));
            }
        }
        return Ok(result);
    }
    let formulation = job.field("formulation")?.text()?;
    let (full, projected) = match formulation {
        "displacement" => (Some(assemble(space, young, poisson, order)?), None),
        "projected_bulk_q1" => (
            None,
            Some(assemble_projected_bulk(space, young, poisson, order)?),
        ),
        _ => return Err("unknown mechanical formulation"),
    };
    let k = full
        .as_ref()
        .or_else(|| projected.as_ref().map(|p| &p.stiffness))
        .ok_or("missing mechanical stiffness")?;
    let mut solved = solve_constrained(k, [&lifts[0], &lifts[1], &lifts[2]], &force, 1e-10, 20000)?;
    let mut mean_stress = BTreeMap::new();
    if let Some(system) = &projected {
        for block in &space.blocks {
            let values = system.mean_stress_coefficients(&block.id, &solved.displacement_mm)?;
            mean_stress.insert(block.id.clone(), numbers(values));
        }
    }
    let residual = k.equilibrium_residual(&solved.displacement_mm, &force)?;
    let mut root_forces = vec![[0.; 3]; adaptive.lift.free_roots.len()];
    for line in &adaptive.lift.lines {
        for &(j, w) in &line.terms {
            let weight = w.numerator() as f64 / w.denominator() as f64;
            for c in 0..3 {
                root_forces[j][c] -= weight * residual[3 * line.dof + c];
            }
        }
    }
    // Rigid motions have identically zero strain. Adding this exact physical
    // field after the solve avoids subtracting nearly equal elastic energies.
    if let Some(frame) = rigid_frame {
        for (id, p) in hom.iter().enumerate() {
            for c in 0..3 {
                let r = frame.value(c, *p)?;
                solved.displacement_mm[3 * id + c] += r.numerator() as f64 / r.denominator() as f64;
            }
        }
    }
    let mut reactions = [0.; 3];
    let mut support_work = 0.;
    let mut support_scale = 0.;
    let mut reaction_moment = [0.; 3];
    // For q_e = C_e q_free + g_e, virtual work gives reaction
    // (lambda_e, -C_e^T lambda_e) in the unmounted trace-root coordinates.
    let mut constrained_forces = vec![[0.; 3]; adaptive.lift.free_roots.len()];
    for (j, &id) in adaptive.lift.free_roots.iter().enumerate() {
        for (c, mount) in physical_lifts.iter().enumerate() {
            if mount.free_roots.binary_search(&id).is_ok() {
                continue;
            }
            let reaction = root_forces[j][c];
            constrained_forces[j][c] += reaction;
            for &(index, weight) in &mount.lines[id].terms {
                let root = mount.free_roots[index];
                let master = adaptive
                    .lift
                    .free_roots
                    .binary_search(&root)
                    .map_err(|_| "invalid support reaction root")?;
                constrained_forces[master][c] -=
                    reaction * (weight.numerator() as f64 / weight.denominator() as f64);
            }
        }
    }
    for (&id, &r) in adaptive.lift.free_roots.iter().zip(&constrained_forces) {
        support_work += dot(
            r,
            std::array::from_fn(|c| solved.displacement_mm[3 * id + c]),
        );
        support_scale += norm(r);
        add(&mut reactions, r);
        add(&mut reaction_moment, cross(points[id], r));
    }
    if version_two {
        scale = scale.max(support_scale);
    }
    // Zero load and zero reactions are possible for a pure prescribed rigid motion.
    let normalize = |error: f64, scale: f64| {
        if scale == 0. && error == 0. {
            0.
        } else {
            error / scale
        }
    };
    let balance = normalize(
        norm(std::array::from_fn(|c| reactions[c] + applied[c])),
        scale,
    );
    let arm = points.iter().map(|&p| norm(p)).fold(1., f64::max);
    let moment_balance = normalize(
        norm(std::array::from_fn(|c| reaction_moment[c] + moment[c])),
        scale * arm,
    );
    if !balance.is_finite()
        || !moment_balance.is_finite()
        || balance > 1e-8
        || moment_balance > 1e-8
    {
        return Err("mechanical force or moment imbalance");
    }
    let applied_work = force
        .iter()
        .zip(&solved.displacement_mm)
        .map(|(a, b)| a * b)
        .sum::<f64>();
    let body_work = body.as_ref().map_or(0., |load| {
        load.force
            .iter()
            .zip(&solved.displacement_mm)
            .map(|(f, u)| f * u)
            .sum::<f64>()
    });
    let traction_work = surface_force.as_ref().map_or(applied_work, |surface| {
        surface
            .iter()
            .zip(&solved.displacement_mm)
            .map(|(f, u)| f * u)
            .sum::<f64>()
    });
    let work = applied_work + support_work;
    let twice_energy = 2. * solved.strain_energy_n_mm;
    let work_scale = work
        .abs()
        .max(twice_energy.abs())
        .max(traction_work.abs())
        .max(body_work.abs())
        .max(support_work.abs());
    let work_error = if work_scale == 0. {
        0.
    } else {
        (work - twice_energy).abs() / work_scale
    };
    if !work_error.is_finite() || work_error > 1e-8 {
        return Err("mechanical work-energy imbalance");
    }
    let mut result = object([
        ("relative_work_energy_error", n(work_error)),
        ("mean_stress_q1_mpa_by_block", Json::Object(mean_stress)),
        ("schema", text("solidcraft_mechanical_result_v1")),
        ("source_family", model.field("schema")?.clone()),
        ("source_audit", audit),
        ("job", job.clone()),
        (
            "analysis_blocks",
            Json::Array(space.blocks.iter().map(analysis_block_json).collect()),
        ),
        ("scalar_dofs", n(space.ndof)),
        (
            "component_lifts",
            Json::Array(physical_lifts.iter().map(lift_json).collect()),
        ),
        ("supports", Json::Array(supports)),
        ("loads", Json::Array(loads)),
        ("rigid_constraint_rank", n(rank)),
        ("force_n_xyz", numbers(force)),
        ("displacement_mm_xyz", numbers(solved.displacement_mm)),
        ("applied_resultant_n", numbers(applied)),
        ("applied_moment_n_mm", numbers(moment)),
        ("support_reaction_n", numbers(reactions)),
        ("support_moment_n_mm", numbers(reaction_moment)),
        ("relative_force_balance", n(balance)),
        ("relative_moment_balance", n(moment_balance)),
        ("volume_mm3", n(k.integrated_volume)),
        ("strain_energy_n_mm", n(solved.strain_energy_n_mm)),
        ("external_work_n_mm", n(work)),
        ("traction_work_n_mm", n(traction_work)),
        ("support_work_n_mm", n(support_work)),
        ("work_normalization_n_mm", n(work_scale)),
        ("equilibrium_normalization_n", n(scale)),
        ("iterations", n(solved.iterations)),
        ("relative_residual", n(solved.relative_residual)),
    ]);
    if let (Some(load), Json::Object(fields)) = (body, &mut result) {
        fields.insert(
            "body_load".into(),
            object([
                ("volume_mm3", n(load.volume)),
                ("resultant_n", numbers(load.resultant)),
                ("moment_n_mm", numbers(load.moment)),
                ("work_n_mm", n(body_work)),
            ]),
        );
    }
    Ok(result)
}
