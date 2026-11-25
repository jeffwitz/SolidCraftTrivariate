//! Compare two local fillet refinement ladders against the frozen J4b references.
use solidcraft_trivariate::{
    adaptive_partition::{AdaptivePartition, LeafKey},
    adaptive_space::{build_adaptive_space, AdaptiveSpace},
    analysis_space::parameterized_bracket_space,
    bracket_parameters::build_bracket_parameters,
    bracket_parameters_audit::audit_bracket_parameters,
    bracket_solve::{bracket_load, bracket_tip_displacement},
    constrained_elasticity::solve_constrained,
    constraint_system::{Lift, LiftLine, Ratio},
    elasticity::assemble,
    json::{object, Json},
    stress::fillet_limits,
};
use std::{collections::BTreeSet, path::Path, time::Instant};

#[path = "support/bracket_field.rs"]
mod bracket_field;

fn n(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn s(x: impl ToString) -> Json {
    Json::String(x.to_string())
}
fn vec_json(v: &[f64]) -> Json {
    Json::Array(v.iter().copied().map(n).collect())
}
fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}
fn scalar(
    block: &solidcraft_trivariate::analysis_space::AnalysisBlock,
    i: usize,
    j: usize,
    k: usize,
    c: usize,
) -> f64 {
    let p = block.controls[(i * 3 + j) * 3 + k];
    p[c] as f64 / p[3] as f64
}

fn target_parameters(
    parent: &solidcraft_trivariate::analysis_space::AnalysisSpace,
) -> Result<Vec<(String, [f64; 3])>, &'static str> {
    // Call the production probe to keep the requested physical location tied to
    // the exact IDs selected by stress::fillet_limits.
    let ids = fillet_limits(parent, &vec![0.; 3 * parent.ndof])?
        .into_iter()
        .map(|x| x.0)
        .collect::<BTreeSet<_>>();
    let q = 14. - 6. / 2_f64.sqrt();
    let mut points = Vec::new();
    for block in parent.blocks.iter().filter(|b| ids.contains(&b.id)) {
        for i in [0, 2] {
            for j in [0, 2] {
                let p0 = [
                    scalar(block, i, j, 0, 0),
                    scalar(block, i, j, 0, 1),
                    scalar(block, i, j, 0, 2),
                ];
                let p2 = [
                    scalar(block, i, j, 2, 0),
                    scalar(block, i, j, 2, 1),
                    scalar(block, i, j, 2, 2),
                ];
                if (p0[0] - q).abs() < 1e-10
                    && (p0[1] - q).abs() < 1e-10
                    && p0[2] <= 20.
                    && p2[2] >= 20.
                    && p2[2] > p0[2]
                {
                    points.push((
                        block.id.clone(),
                        [
                            i as f64 / 2.,
                            j as f64 / 2.,
                            (20. - p0[2]) / (p2[2] - p0[2]),
                        ],
                    ));
                }
            }
        }
    }
    points.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| a.1[0].total_cmp(&b.1[0]))
            .then_with(|| a.1[1].total_cmp(&b.1[1]))
    });
    points.dedup();
    if points.is_empty() {
        return Err("physical fillet point has no root-chart preimage");
    }
    Ok(points)
}

fn incident_leaves(
    partition: &AdaptivePartition,
    points: &[(String, [f64; 3])],
) -> BTreeSet<LeafKey> {
    let mut result = BTreeSet::new();
    for (root, xyz) in points {
        for leaf in partition.leaves().iter().filter(|leaf| &leaf.root == root) {
            let scale = 1_u32 << leaf.level;
            let inside = (0..3).all(|axis| {
                let lo = leaf.index[axis] as f64 / scale as f64;
                let hi = (leaf.index[axis] + 1) as f64 / scale as f64;
                xyz[axis] >= lo - 1e-14 && xyz[axis] <= hi + 1e-14
            });
            if inside {
                result.insert(leaf.clone());
            }
        }
    }
    result
}

fn refine_local(
    partition: &mut AdaptivePartition,
    points: &[(String, [f64; 3])],
) -> Result<(usize, usize, usize), &'static str> {
    let incident = incident_leaves(partition, points);
    if incident.is_empty() {
        return Err("fillet point is not incident to an active leaf");
    }
    let mut requested = incident.clone();
    for patch in partition.interfaces()? {
        if incident.contains(&patch.first.leaf) || incident.contains(&patch.second.leaf) {
            requested.insert(patch.first.leaf);
            requested.insert(patch.second.leaf);
        }
    }
    let report = partition.refine(&requested.iter().cloned().collect::<Vec<_>>())?;
    Ok((incident.len(), requested.len(), report.induced.len()))
}

fn json_fillet(values: &[(String, f64)]) -> Json {
    Json::Object(values.iter().map(|(id, v)| (id.clone(), n(*v))).collect())
}

fn mounted_zero_lift(build: &AdaptiveSpace, y50: &[usize]) -> Result<Lift, &'static str> {
    let mut fixed = BTreeSet::new();
    for &dof in y50 {
        let line = build
            .lift
            .lines
            .get(dof)
            .ok_or("y50 DOF outside adaptive lift")?;
        fixed.extend(line.terms.iter().map(|term| term.0));
    }
    let mut free_roots = Vec::new();
    let mut reduced_map = vec![None; build.lift.free_roots.len()];
    for (old, &root) in build.lift.free_roots.iter().enumerate() {
        if !fixed.contains(&old) {
            reduced_map[old] = Some(free_roots.len());
            free_roots.push(root);
        }
    }
    let lines = build
        .lift
        .lines
        .iter()
        .map(|line| LiftLine {
            dof: line.dof,
            terms: line
                .terms
                .iter()
                .filter_map(|(old, coefficient)| reduced_map[*old].map(|new| (new, *coefficient)))
                .collect(),
            offset: Ratio::zero(),
        })
        .collect::<Vec<_>>();
    if y50.iter().any(|&dof| !lines[dof].terms.is_empty()) {
        return Err("mounted lift leaves a prescribed y50 coefficient free");
    }
    Ok(Lift { free_roots, lines })
}

fn evaluate_state(
    build: &AdaptiveSpace,
    partition: &AdaptivePartition,
    q: usize,
    tolerance: f64,
    checkpoint: Option<(&Path, [u8; 2])>,
) -> Result<Json, Box<dyn std::error::Error>> {
    let space = &build.space;
    let t = Instant::now();
    let stiffness = assemble(space, 210000., 0.3, q)?;
    let assembly_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let (force, area) = bracket_load(space, q)?;
    let force_s = t.elapsed().as_secs_f64();
    let y50 = space.boundary_dofs.get("y50").ok_or("missing y50")?;
    let t = Instant::now();
    let mounted = mounted_zero_lift(build, y50)?;
    let free = &build.lift;
    let constraint_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let solved = solve_constrained(&stiffness, [&mounted; 3], &force, tolerance, 20000)?;
    let solve_s = t.elapsed().as_secs_f64();
    if let Some((path, radii)) = checkpoint {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        bracket_field::write(path, None, radii, q, space, &solved.displacement_mm)?;
    }
    let t = Instant::now();
    let ku = stiffness.apply(&solved.displacement_mm)?;
    let residual = ku
        .iter()
        .zip(&force)
        .map(|(a, b)| a - b)
        .collect::<Vec<_>>();
    let mut projected = vec![0.; 3 * free.free_roots.len()];
    let mut projected_force = vec![0.; 3 * free.free_roots.len()];
    for line in &free.lines {
        for c in 0..3 {
            for &(root, coefficient) in &line.terms {
                projected[3 * root + c] += residual[3 * line.dof + c]
                    * coefficient.numerator() as f64
                    / coefficient.denominator() as f64;
                projected_force[3 * root + c] += force[3 * line.dof + c]
                    * coefficient.numerator() as f64
                    / coefficient.denominator() as f64;
            }
        }
    }
    let y50_roots = space.boundary_dofs["y50"]
        .iter()
        .flat_map(|&d| free.lines[d].terms.iter().map(|x| x.0))
        .collect::<BTreeSet<_>>();
    let mut reaction = [0.; 3];
    for root in &y50_roots {
        for c in 0..3 {
            reaction[c] += projected[3 * root + c];
        }
    }
    let mut free_res = Vec::new();
    for root in 0..free.free_roots.len() {
        if !y50_roots.contains(&root) {
            free_res.extend_from_slice(&projected[3 * root..3 * root + 3]);
        }
    }
    let mut rhs = Vec::new();
    for root in 0..free.free_roots.len() {
        if !y50_roots.contains(&root) {
            rhs.extend_from_slice(&projected_force[3 * root..3 * root + 3]);
        }
    }
    let true_residual = norm(&free_res) / norm(&rhs).max(f64::MIN_POSITIVE);
    let tip = bracket_tip_displacement(space, &solved.displacement_mm)?;
    let fillet = fillet_limits(space, &solved.displacement_mm)?;
    let recovery_s = t.elapsed().as_secs_f64();
    let _ = partition;
    Ok(object([
        ("q", n(q)),
        ("blocks", n(space.blocks.len())),
        ("scalar_dofs", n(space.ndof)),
        ("masters", n(build.lift.free_roots.len())),
        ("volume_mm3", n(stiffness.integrated_volume)),
        ("loaded_area_mm2", n(area)),
        ("tip_mm", vec_json(&tip)),
        ("reaction_n", vec_json(&reaction)),
        ("fillet_limits_mpa", json_fillet(&fillet)),
        ("relative_residual", n(solved.relative_residual)),
        ("true_projected_residual", n(true_residual)),
        ("iterations", n(solved.iterations)),
        ("strain_energy_n_mm", n(solved.strain_energy_n_mm)),
        ("assembly_s", n(assembly_s)),
        ("force_s", n(force_s)),
        ("constraint_s", n(constraint_s)),
        ("solve_s", n(solve_s)),
        ("recovery_s", n(recovery_s)),
    ]))
}

fn save_report(
    out: &Path,
    rx: u8,
    ry: u8,
    geometry_audit_s: f64,
    targets: &[(String, [f64; 3])],
    paths: &[Json],
    total_wall_s: f64,
) -> std::io::Result<()> {
    if let Some(parent) = out.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let result = object([
        ("radii_ticks", Json::Array(vec![n(rx), n(ry)])),
        ("geometry_build_and_audit_s", n(geometry_audit_s)),
        (
            "fillet_parameter_targets",
            Json::Array(
                targets
                    .iter()
                    .map(|(id, uvw)| {
                        object([
                            ("root", s(id)),
                            ("uvw", Json::Array(uvw.iter().copied().map(n).collect())),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("paths", Json::Array(paths.to_vec())),
        ("total_wall_s", n(total_wall_s)),
    ]);
    let tmp = out.with_extension("json.tmp");
    std::fs::write(&tmp, format!("{result}\n"))?;
    std::fs::rename(tmp, out)
}

fn run_radius(rx: u8, ry: u8, out: &Path) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(parent) = out.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let total = Instant::now();
    let t = Instant::now();
    let model = build_bracket_parameters(rx, ry)?;
    let audit = audit_bracket_parameters(rx, ry, &model.blocks)?;
    let root_space = parameterized_bracket_space(rx, ry, &model.blocks)?;
    let geometry_audit_s = t.elapsed().as_secs_f64();
    let targets = target_parameters(&root_space)?;
    let roots = root_space
        .blocks
        .iter()
        .map(|b| b.id.clone())
        .collect::<Vec<_>>();
    let mut paths = Vec::new();
    for (uniform_depth, local_rounds, label) in [(1_u8, 4_usize, "from_h1"), (2, 3, "from_h2")] {
        let mut partition = AdaptivePartition::new(&roots, &audit.certificate.incidence, 7, 12000)?;
        let mut partition_s = 0.;
        for _ in 0..uniform_depth {
            let leaves = partition.leaves().iter().cloned().collect::<Vec<_>>();
            let t = Instant::now();
            partition.refine(&leaves)?;
            partition_s += t.elapsed().as_secs_f64();
        }
        let mut states = Vec::new();
        for step in 0..=local_rounds {
            let state_start = Instant::now();
            let mut local_counts = [0_usize; 3];
            if step > 0 {
                let t = Instant::now();
                match refine_local(&mut partition, &targets) {
                    Ok(counts) => {
                        partition_s += t.elapsed().as_secs_f64();
                        local_counts = [counts.0, counts.1, counts.2];
                    }
                    Err(e) => {
                        states.push(object([
                            ("step", n(step)),
                            ("status", s(format!("refused: {e}"))),
                            ("leaves", n(partition.leaves().len())),
                        ]));
                        let mut partial_paths = paths.clone();
                        partial_paths.push(object([
                            ("path", s(label)),
                            ("uniform_depth", n(uniform_depth)),
                            ("local_rounds", n(local_rounds)),
                            ("states", Json::Array(states.clone())),
                        ]));
                        save_report(
                            out,
                            rx,
                            ry,
                            geometry_audit_s,
                            &targets,
                            &partial_paths,
                            total.elapsed().as_secs_f64(),
                        )?;
                        break;
                    }
                }
            }
            let t = Instant::now();
            let build = build_adaptive_space(&root_space, &partition)?;
            let constraints_s = t.elapsed().as_secs_f64();
            let mut qresults = Vec::new();
            for q in [6_usize, 8] {
                qresults.push(evaluate_state(&build, &partition, q, 1e-11, None)?);
            }
            states.push(object([
                ("step", n(step)),
                ("status", s("completed")),
                ("leaves", n(partition.leaves().len())),
                (
                    "max_depth",
                    n(partition
                        .leaves()
                        .iter()
                        .map(|leaf| leaf.level)
                        .max()
                        .unwrap_or(0)),
                ),
                ("partition_cumulative_s", n(partition_s)),
                ("incident_point_leaves", n(local_counts[0])),
                ("requested_point_and_face_neighbours", n(local_counts[1])),
                ("balance_induced_splits", n(local_counts[2])),
                ("constraints_s", n(constraints_s)),
                ("wall_s", n(state_start.elapsed().as_secs_f64())),
                ("q6_q8", Json::Array(qresults)),
            ]));
            let mut partial_paths = paths.clone();
            partial_paths.push(object([
                ("path", s(label)),
                ("uniform_depth", n(uniform_depth)),
                ("local_rounds", n(local_rounds)),
                ("states", Json::Array(states.clone())),
            ]));
            save_report(
                out,
                rx,
                ry,
                geometry_audit_s,
                &targets,
                &partial_paths,
                total.elapsed().as_secs_f64(),
            )?;
            eprintln!(
                "radius {rx}/{ry} {label} step {step}: {} leaves",
                partition.leaves().len()
            );
        }
        paths.push(object([
            ("path", s(label)),
            ("uniform_depth", n(uniform_depth)),
            ("local_rounds", n(local_rounds)),
            ("states", Json::Array(states)),
        ]));
    }
    save_report(
        out,
        rx,
        ry,
        geometry_audit_s,
        &targets,
        &paths,
        total.elapsed().as_secs_f64(),
    )?;
    Ok(())
}

fn run_tier(
    rx: u8,
    ry: u8,
    path: &str,
    local_steps: usize,
    out: &Path,
    checkpoint: Option<&Path>,
    field_order: Option<usize>,
) -> Result<(), Box<dyn std::error::Error>> {
    if local_steps > 4 {
        return Err("invalid local tier".into());
    }
    let total = Instant::now();
    let t = Instant::now();
    let model = build_bracket_parameters(rx, ry)?;
    let audit = audit_bracket_parameters(rx, ry, &model.blocks)?;
    let parent = parameterized_bracket_space(rx, ry, &model.blocks)?;
    let geometry_audit_s = t.elapsed().as_secs_f64();
    let targets = target_parameters(&parent)?;
    let roots = parent
        .blocks
        .iter()
        .map(|b| b.id.clone())
        .collect::<Vec<_>>();
    let uniform_depth = match path {
        "from_h1" => 1_u8,
        "from_h2" => 2_u8,
        _ => return Err("invalid local tier path".into()),
    };
    let mut partition = AdaptivePartition::new(&roots, &audit.certificate.incidence, 7, 12000)?;
    let t = Instant::now();
    for _ in 0..uniform_depth {
        let leaves = partition.leaves().iter().cloned().collect::<Vec<_>>();
        partition.refine(&leaves)?;
    }
    let uniform_partition_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let mut local_counts = Vec::new();
    for _ in 0..local_steps {
        local_counts.push(refine_local(&mut partition, &targets)?);
    }
    let local_partition_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let build = build_adaptive_space(&parent, &partition)?;
    let constraints_s = t.elapsed().as_secs_f64();
    let q = field_order.unwrap_or(8);
    let tolerance = if field_order.is_some() { 5e-11 } else { 1e-11 };
    let result = evaluate_state(
        &build,
        &partition,
        q,
        tolerance,
        checkpoint.map(|p| (p, [rx, ry])),
    )?;
    let report = object([
        ("radii_ticks", Json::Array(vec![n(rx), n(ry)])),
        ("path", s(path)),
        ("uniform_depth", n(uniform_depth)),
        ("local_steps", n(local_steps)),
        ("leaves", n(partition.leaves().len())),
        (
            "max_depth",
            n(partition
                .leaves()
                .iter()
                .map(|leaf| leaf.level)
                .max()
                .unwrap_or(0)),
        ),
        ("geometry_build_and_audit_s", n(geometry_audit_s)),
        ("uniform_partition_s", n(uniform_partition_s)),
        ("local_partition_s", n(local_partition_s)),
        ("constraints_s", n(constraints_s)),
        (
            "local_refinement_counts",
            Json::Array(
                local_counts
                    .iter()
                    .map(|(incident, requested, induced)| {
                        object([
                            ("incident", n(*incident)),
                            ("requested", n(*requested)),
                            ("induced", n(*induced)),
                        ])
                    })
                    .collect(),
            ),
        ),
        (if q == 6 { "q6" } else { "q8" }, result),
        ("total_wall_s", n(total.elapsed().as_secs_f64())),
    ]);
    if let Some(parent) = out.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let tmp = out.with_extension("json.tmp");
    std::fs::write(&tmp, format!("{report}\n"))?;
    std::fs::rename(tmp, out)?;
    Ok(())
}

fn run_path(rx: u8, ry: u8, path: &str, out: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let (uniform_depth, local_rounds) = match path {
        "from_h1" => (1_u8, 4_usize),
        "from_h2" => (2_u8, 3_usize),
        _ => return Err("invalid local path".into()),
    };
    if let Some(parent) = out.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let total = Instant::now();
    let t = Instant::now();
    let model = build_bracket_parameters(rx, ry)?;
    let audit = audit_bracket_parameters(rx, ry, &model.blocks)?;
    let root_space = parameterized_bracket_space(rx, ry, &model.blocks)?;
    let geometry_audit_s = t.elapsed().as_secs_f64();
    let targets = target_parameters(&root_space)?;
    let roots = root_space
        .blocks
        .iter()
        .map(|b| b.id.clone())
        .collect::<Vec<_>>();
    let mut partition = AdaptivePartition::new(&roots, &audit.certificate.incidence, 7, 12000)?;
    let t = Instant::now();
    for _ in 0..uniform_depth {
        let leaves = partition.leaves().iter().cloned().collect::<Vec<_>>();
        partition.refine(&leaves)?;
    }
    let mut partition_s = t.elapsed().as_secs_f64();
    let mut states = Vec::new();
    for step in 0..=local_rounds {
        let state_start = Instant::now();
        let mut counts = [0_usize; 3];
        if step > 0 {
            let t = Instant::now();
            let c = refine_local(&mut partition, &targets)?;
            partition_s += t.elapsed().as_secs_f64();
            counts = [c.0, c.1, c.2];
        }
        let t = Instant::now();
        let build = build_adaptive_space(&root_space, &partition)?;
        let constraints_s = t.elapsed().as_secs_f64();
        let qresults = vec![
            evaluate_state(&build, &partition, 6, 1e-11, None)?,
            evaluate_state(&build, &partition, 8, 1e-11, None)?,
        ];
        states.push(object([
            ("step", n(step)),
            ("status", s("completed")),
            ("leaves", n(partition.leaves().len())),
            (
                "max_depth",
                n(partition
                    .leaves()
                    .iter()
                    .map(|leaf| leaf.level)
                    .max()
                    .unwrap_or(0)),
            ),
            ("partition_cumulative_s", n(partition_s)),
            ("incident_point_leaves", n(counts[0])),
            ("requested_point_and_face_neighbours", n(counts[1])),
            ("balance_induced_splits", n(counts[2])),
            ("constraints_s", n(constraints_s)),
            ("wall_s", n(state_start.elapsed().as_secs_f64())),
            ("q6_q8", Json::Array(qresults)),
        ]));
        let report = object([
            ("radii_ticks", Json::Array(vec![n(rx), n(ry)])),
            ("geometry_build_and_audit_s", n(geometry_audit_s)),
            (
                "fillet_parameter_targets",
                Json::Array(
                    targets
                        .iter()
                        .map(|(id, uvw)| {
                            object([
                                ("root", s(id)),
                                ("uvw", Json::Array(uvw.iter().copied().map(n).collect())),
                            ])
                        })
                        .collect(),
                ),
            ),
            (
                "paths",
                Json::Array(vec![object([
                    ("path", s(path)),
                    ("uniform_depth", n(uniform_depth)),
                    ("local_rounds", n(local_rounds)),
                    ("states", Json::Array(states.clone())),
                ])]),
            ),
            ("total_wall_s", n(total.elapsed().as_secs_f64())),
        ]);
        let tmp = out.with_extension("json.tmp");
        std::fs::write(&tmp, format!("{report}\n"))?;
        std::fs::rename(tmp, out)?;
        eprintln!(
            "radius {rx}/{ry} {path} step {step}: {} leaves",
            partition.leaves().len()
        );
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a = std::env::args().skip(1).collect::<Vec<_>>();
    if let [flag, rx, ry, path, out] = a.as_slice() {
        if flag == "--path" {
            return run_path(rx.parse()?, ry.parse()?, path, Path::new(out));
        }
    }
    if let [flag, rx, ry, path, step, out] = a.as_slice() {
        if flag == "--tier" {
            return run_tier(
                rx.parse()?,
                ry.parse()?,
                path,
                step.parse()?,
                Path::new(out),
                None,
                None,
            );
        }
    }
    if let [flag, rx, ry, path, step, out, checkpoint] = a.as_slice() {
        if flag == "--tier" {
            return run_tier(
                rx.parse()?,
                ry.parse()?,
                path,
                step.parse()?,
                Path::new(out),
                Some(Path::new(checkpoint)),
                None,
            );
        }
    }
    if let [flag, rx, ry, step, order, out, checkpoint] = a.as_slice() {
        if flag == "--field-tier" {
            let q = order.parse::<usize>()?;
            if ![6, 8].contains(&q) {
                return Err("field tier requires order 6 or 8".into());
            }
            return run_tier(
                rx.parse()?,
                ry.parse()?,
                "from_h2",
                step.parse()?,
                Path::new(out),
                (q == 8).then(|| Path::new(checkpoint)),
                Some(q),
            );
        }
    }
    let [rx, ry, out] = a.as_slice() else {
        return Err("usage: local_fillet_cost RADIUS_X_TICKS RADIUS_Y_TICKS OUTPUT_JSON".into());
    };
    run_radius(rx.parse()?, ry.parse()?, Path::new(out))
}
