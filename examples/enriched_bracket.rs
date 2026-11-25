//! Compare residual marking with marking from a uniformly enriched solution.
#[path = "support/bracket_field.rs"]
mod bracket_field;
#[path = "support/enrichment.rs"]
mod enrichment;

use enrichment::{bulk_mark, difference_energy_by_parent, galerkin_identity_relative_error};
use solidcraft_trivariate::{
    adaptive_estimator::estimate_residual,
    adaptive_partition::{AdaptivePartition, LeafKey},
    adaptive_space::{build_adaptive_space, AdaptiveSpace},
    adaptive_transfer::build_transfer,
    analysis_space::{parameterized_bracket_space, AnalysisSpace},
    bracket_parameters::build_bracket_parameters,
    bracket_parameters_audit::audit_bracket_parameters,
    bracket_solve::bracket_load,
    constrained_elasticity::solve_constrained,
    constraint_system::{resolve, Ratio},
    elasticity::{assemble, Stiffness},
    incidence::{BlockSide, Incidence},
    json::{object, Json},
};
use std::{collections::BTreeSet, path::Path, time::Instant};

type RootSpace = (
    AnalysisSpace,
    Incidence,
    std::collections::BTreeMap<BlockSide, String>,
);

const YOUNG_MPA: f64 = 210000.0;
const POISSON: f64 = 0.3;
const GAUSS_ORDER: usize = 8;
const MAX_ACCEPTED_CELLS: usize = 4000;
const MAX_TRIAL_CELLS: usize = 16000;
const MAX_DEPTH: u8 = 5;

fn n(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn leaf(key: &LeafKey) -> Json {
    object([
        ("root", Json::String(key.root.clone())),
        ("level", n(key.level)),
        ("ijk", Json::Array(key.index.map(n).into())),
    ])
}
fn save_json(path: &Path, value: &Json) -> std::io::Result<()> {
    let temp = path.with_extension("tmp");
    std::fs::write(&temp, format!("{value}\n"))?;
    std::fs::rename(temp, path)
}
fn rss_kib() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    status.lines().find_map(|line| {
        line.strip_prefix("VmRSS:")
            .and_then(|rest| rest.split_whitespace().next())
            .and_then(|value| value.parse().ok())
    })
}
fn vector_energy(stiffness: &Stiffness, field: &[f64]) -> Result<f64, &'static str> {
    let applied = stiffness.apply(field)?;
    Ok(field.iter().zip(applied).map(|(a, b)| a * b).sum())
}

struct SolvedSpace {
    build: AdaptiveSpace,
    stiffness: Stiffness,
    field: Vec<f64>,
    energy: f64,
    independent_scalar_dofs: usize,
    free_vector_dofs: usize,
    iterations: usize,
    relative_residual: f64,
    volume_mm3: f64,
    area_mm2: f64,
    prepare_s: f64,
    assembly_s: f64,
    load_s: f64,
    solve_s: f64,
    wall_s: f64,
    rss_kib: Option<u64>,
}

fn solve_space(
    parent: &AnalysisSpace,
    partition: &AdaptivePartition,
    radius_ticks: u8,
    round: usize,
    role: &str,
) -> Result<SolvedSpace, Box<dyn std::error::Error>> {
    let wall_start = Instant::now();
    eprintln!("A6_RESOLUTION_START|radius={radius_ticks}|round={round}|kind={role}");
    let start = Instant::now();
    let build = build_adaptive_space(parent, partition)?;
    let prepare_s = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let stiffness = assemble(&build.space, YOUNG_MPA, POISSON, GAUSS_ORDER)?;
    let assembly_s = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let (force, area_mm2) = bracket_load(&build.space, GAUSS_ORDER)?;
    let prescribed = build.space.boundary_dofs["y50"]
        .iter()
        .map(|&dof| (dof, Ratio::zero()))
        .collect::<Vec<_>>();
    let lift = resolve(build.space.ndof, &build.constraints, &prescribed)?;
    let load_s = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let solved = solve_constrained(&stiffness, [&lift; 3], &force, 1e-10, 20000)?;
    let solve_s = start.elapsed().as_secs_f64();
    let wall_s = wall_start.elapsed().as_secs_f64();
    eprintln!(
        "A6_RESOLUTION_END|radius={radius_ticks}|round={round}|kind={role}|elapsed_s={wall_s:.9}|rss_kib={}",
        rss_kib().map_or_else(|| "unknown".to_owned(), |value| value.to_string())
    );
    Ok(SolvedSpace {
        independent_scalar_dofs: build.lift.free_roots.len(),
        free_vector_dofs: lift.free_roots.len() * 3,
        energy: solved.strain_energy_n_mm,
        iterations: solved.iterations,
        relative_residual: solved.relative_residual,
        volume_mm3: stiffness.integrated_volume,
        area_mm2,
        field: solved.displacement_mm,
        build,
        stiffness,
        prepare_s,
        assembly_s,
        load_s,
        solve_s,
        wall_s,
        rss_kib: rss_kib(),
    })
}

fn ancestor(key: &LeafKey, level: u8) -> LeafKey {
    let shift = key.level - level;
    LeafKey {
        root: key.root.clone(),
        level,
        index: key.index.map(|coordinate| coordinate >> shift),
    }
}

fn same_partition_at_trial_depth(
    source: &AdaptivePartition,
    incidence: &Incidence,
) -> Result<AdaptivePartition, &'static str> {
    let mut target =
        AdaptivePartition::new(source.roots(), incidence, MAX_DEPTH + 1, MAX_TRIAL_CELLS)?;
    let desired = source.leaves().clone();
    for depth in 0..MAX_DEPTH {
        let requests = desired
            .iter()
            .filter(|leaf| leaf.level > depth)
            .map(|leaf| ancestor(leaf, depth))
            .filter(|candidate| target.leaves().contains(candidate))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        if !requests.is_empty() {
            target.refine(&requests)?;
        }
    }
    if target.leaves() != &desired {
        return Err("trial partition reconstruction changed the current topology");
    }
    Ok(target)
}

fn all_one_level_children(
    source: &AdaptivePartition,
    incidence: &Incidence,
) -> Result<AdaptivePartition, &'static str> {
    if source.leaves().iter().any(|leaf| leaf.level > MAX_DEPTH) {
        return Err("trial_depth_limit");
    }
    let required = source
        .leaves()
        .len()
        .checked_mul(8)
        .ok_or("trial leaf count overflow")?;
    if required > MAX_TRIAL_CELLS {
        return Err("trial_cell_limit");
    }
    let mut target = same_partition_at_trial_depth(source, incidence)?;
    let requests = target.leaves().iter().cloned().collect::<Vec<_>>();
    target.refine(&requests)?;
    if target.leaves().len() != required {
        return Err("uniform trial did not create exactly eight children per source leaf");
    }
    Ok(target)
}

fn root_radius_space(radius_ticks: u8) -> Result<RootSpace, Box<dyn std::error::Error>> {
    let model = build_bracket_parameters(radius_ticks, 64)?;
    let audit = audit_bracket_parameters(radius_ticks, 64, &model.blocks)?;
    let space = parameterized_bracket_space(radius_ticks, 64, &model.blocks)?;
    let tags = audit.certificate.design_faces.iter().cloned().collect();
    Ok((space, audit.certificate.incidence, tags))
}

fn write_field(
    directory: &Path,
    name: &str,
    radius_ticks: u8,
    solved: &SolvedSpace,
) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all(directory)?;
    let path = directory.join(format!("{name}.bin"));
    let temp = path.with_extension("part");
    bracket_field::write(
        &temp,
        None,
        [radius_ticks, 64],
        GAUSS_ORDER,
        &solved.build.space,
        &solved.field,
    )?;
    std::fs::rename(temp, path)?;
    Ok(())
}

fn run(radius_ticks: u8, output: &Path, rounds: usize) -> Result<(), Box<dyn std::error::Error>> {
    if !matches!(radius_ticks, 56 | 64) || !(1..=6).contains(&rounds) {
        return Err("expected radius ticks 56 or 64 and 1..6 rounds".into());
    }
    std::fs::create_dir_all(output)?;
    let fields = output.join("fields");
    std::fs::create_dir_all(&fields)?;
    let (parent, incidence, tags) = root_radius_space(radius_ticks)?;
    let roots = parent
        .blocks
        .iter()
        .map(|block| block.id.clone())
        .collect::<Vec<_>>();
    let mut current_partition =
        AdaptivePartition::new(&roots, &incidence, MAX_DEPTH, MAX_ACCEPTED_CELLS)?;
    let mut decisions = Vec::new();
    let mut stop_reason = String::from("round_limit");
    for round in 0..rounds {
        eprintln!(
            "A6 radius {radius_ticks}: round {round}, {} current leaves",
            current_partition.leaves().len()
        );
        let current = match solve_space(&parent, &current_partition, radius_ticks, round, "current")
        {
            Ok(value) => value,
            Err(error) => {
                stop_reason = format!("current_solve_failed:{error}");
                break;
            }
        };
        write_field(&fields, &format!("round-{round}"), radius_ticks, &current)?;
        let residual_start = Instant::now();
        let residual = estimate_residual(
            &current.build,
            &current_partition,
            &current.field,
            YOUNG_MPA,
            POISSON,
            4,
            [0.; 3],
            |side, _, _| match tags[&side.root_side].as_str() {
                "y50" => None,
                "x60" => Some([0., -1000. / 320., 0.]),
                _ => Some([0.; 3]),
            },
        )?;
        let residual_marked = residual.bulk_mark(0.5)?;
        let residual_s = residual_start.elapsed().as_secs_f64();
        let target_partition = match all_one_level_children(&current_partition, &incidence) {
            Ok(value) => value,
            Err(reason) => {
                stop_reason = reason.to_owned();
                decisions.push(object([
                    ("round", n(round)),
                    ("source_leaves", n(current_partition.leaves().len())),
                    ("trial_status", Json::String("not_started".into())),
                    ("stop_reason", Json::String(stop_reason.clone())),
                ]));
                break;
            }
        };
        let target = match solve_space(
            &parent,
            &target_partition,
            radius_ticks,
            round,
            "uniform_enrichment",
        ) {
            Ok(value) => value,
            Err(error) => {
                stop_reason = format!("trial_solve_failed:{error}");
                break;
            }
        };
        write_field(&fields, &format!("enriched-{round}"), radius_ticks, &target)?;
        let transfer_start = Instant::now();
        let transfer = build_transfer(
            &current.build,
            &current_partition,
            &target.build,
            &target_partition,
        )?;
        let prolonged = transfer.transfer_full_field(&current.field)?;
        let transfer_s = transfer_start.elapsed().as_secs_f64();
        let delta = target
            .field
            .iter()
            .zip(&prolonged)
            .map(|(fine, coarse)| fine - coarse)
            .collect::<Vec<_>>();
        let quadrature_start = Instant::now();
        let local = difference_energy_by_parent(
            &target.build,
            &target_partition,
            &current_partition,
            &delta,
            YOUNG_MPA,
            POISSON,
            GAUSS_ORDER,
        )?;
        let difference_energy = local.values().sum::<f64>();
        let global_difference_energy = vector_energy(&target.stiffness, &delta)?;
        let local_global_relative_error = (difference_energy - global_difference_energy).abs()
            / global_difference_energy.abs().max(1e-30);
        let identity_relative_error =
            galerkin_identity_relative_error(difference_energy, current.energy, target.energy)?;
        let quadrature_s = quadrature_start.elapsed().as_secs_f64();
        let marked = bulk_mark(&local, 0.5)?;
        let energy_mark_set = marked.iter().cloned().collect::<BTreeSet<_>>();
        let residual_mark_set = residual_marked.iter().cloned().collect::<BTreeSet<_>>();
        let overlap = energy_mark_set.intersection(&residual_mark_set).count();
        let union = energy_mark_set.union(&residual_mark_set).count();
        let residual_rows = residual
            .by_leaf
            .iter()
            .map(|(key, value)| {
                object([
                    ("leaf", leaf(key)),
                    ("volume_squared", n(value.volume_squared)),
                    ("interface_squared", n(value.interior_squared)),
                    ("natural_squared", n(value.natural_squared)),
                    ("total_squared", n(value.total_squared)),
                ])
            })
            .collect::<Vec<_>>();
        let contribution_rows = local
            .iter()
            .map(|(key, value)| {
                object([
                    ("leaf", leaf(key)),
                    ("difference_strain_energy_n_mm", n(*value)),
                ])
            })
            .collect::<Vec<_>>();
        let marked_energy = marked.iter().map(|key| local[key]).sum::<f64>();
        let record = object([
            ("round", n(round)),
            ("current_leaves", n(current_partition.leaves().len())),
            ("trial_leaves", n(target_partition.leaves().len())),
            ("current_scalar_dofs", n(current.build.space.ndof)),
            (
                "current_independent_scalar_dofs",
                n(current.independent_scalar_dofs),
            ),
            ("current_free_vector_dofs", n(current.free_vector_dofs)),
            ("trial_scalar_dofs", n(target.build.space.ndof)),
            (
                "trial_independent_scalar_dofs",
                n(target.independent_scalar_dofs),
            ),
            ("trial_free_vector_dofs", n(target.free_vector_dofs)),
            (
                "current_stiffness_block_entries",
                n(current.stiffness.columns.len()),
            ),
            (
                "trial_stiffness_block_entries",
                n(target.stiffness.columns.len()),
            ),
            ("current_volume_mm3", n(current.volume_mm3)),
            ("trial_volume_mm3", n(target.volume_mm3)),
            ("current_loaded_area_mm2", n(current.area_mm2)),
            ("trial_loaded_area_mm2", n(target.area_mm2)),
            ("current_energy_n_mm", n(current.energy)),
            ("trial_energy_n_mm", n(target.energy)),
            ("current_iterations", n(current.iterations)),
            ("trial_iterations", n(target.iterations)),
            (
                "current_true_residual_relative",
                n(current.relative_residual),
            ),
            ("trial_true_residual_relative", n(target.relative_residual)),
            ("difference_strain_energy_n_mm", n(difference_energy)),
            (
                "twice_energy_gain_n_mm",
                n(2.0 * (target.energy - current.energy)),
            ),
            (
                "galerkin_identity_relative_error",
                n(identity_relative_error),
            ),
            (
                "identity_qualified",
                Json::Bool(identity_relative_error <= 1e-5),
            ),
            (
                "local_vs_global_quadrature_relative_error",
                n(local_global_relative_error),
            ),
            ("marked", Json::Array(marked.iter().map(leaf).collect())),
            (
                "residual_marked",
                Json::Array(residual_marked.iter().map(leaf).collect()),
            ),
            ("residual_eta", n(residual.indicator)),
            (
                "residual_error_upper_bound_available",
                Json::Bool(residual.error_upper_bound_available),
            ),
            ("residual_marked_count", n(residual_marked.len())),
            ("energy_marked_count", n(marked.len())),
            ("marked_overlap_count", n(overlap)),
            ("marked_union_count", n(union)),
            (
                "marked_jaccard",
                n(if union == 0 {
                    1.0
                } else {
                    overlap as f64 / union as f64
                }),
            ),
            ("residual_contributions", Json::Array(residual_rows)),
            (
                "marked_fraction_energy",
                n(marked_energy / difference_energy.max(1e-30)),
            ),
            ("local_contributions", Json::Array(contribution_rows)),
            ("current_prepare_s", n(current.prepare_s)),
            ("current_assembly_s", n(current.assembly_s)),
            ("current_load_s", n(current.load_s)),
            ("current_solve_s", n(current.solve_s)),
            ("current_resolution_wall_s", n(current.wall_s)),
            (
                "current_rss_kib_after_solve",
                current.rss_kib.map_or(Json::Null, n),
            ),
            ("trial_prepare_s", n(target.prepare_s)),
            ("trial_assembly_s", n(target.assembly_s)),
            ("trial_load_s", n(target.load_s)),
            ("trial_solve_s", n(target.solve_s)),
            ("trial_resolution_wall_s", n(target.wall_s)),
            (
                "trial_rss_kib_after_solve",
                target.rss_kib.map_or(Json::Null, n),
            ),
            ("transfer_s", n(transfer_s)),
            ("difference_quadrature_s", n(quadrature_s)),
            ("residual_estimator_s", n(residual_s)),
        ]);
        decisions.push(record);
        let data = object([
            ("schema", Json::String("A6_enrichment_v1".into())),
            ("radius_x_ticks", n(radius_ticks)),
            ("radius_y_ticks", n(64)),
            ("max_accepted_leaves", n(MAX_ACCEPTED_CELLS)),
            ("max_trial_leaves", n(MAX_TRIAL_CELLS)),
            ("max_rounds", n(rounds)),
            ("decisions", Json::Array(decisions.clone())),
            ("stop_reason", Json::String(stop_reason.clone())),
        ]);
        save_json(&output.join("progress.json"), &data)?;
        if identity_relative_error > 1e-5 {
            eprintln!("A6 warning: Galerkin identity discrepancy {identity_relative_error:.3e}");
        }
        if round + 1 == rounds {
            break;
        }
        if marked.is_empty() {
            stop_reason = "zero_enrichment_difference".into();
            break;
        }
        match current_partition.refine(&marked) {
            Ok(update) => {
                if current_partition.leaves().len() > MAX_ACCEPTED_CELLS {
                    stop_reason = "accepted_cell_limit".into();
                    break;
                }
                if let Some(Json::Object(fields)) = decisions.last_mut() {
                    fields.insert(
                        "accepted_refinement_induced".into(),
                        Json::Array(update.induced.iter().map(leaf).collect()),
                    );
                    fields.insert(
                        "accepted_leaf_count_after_refine".into(),
                        n(current_partition.leaves().len()),
                    );
                }
            }
            Err(message) => {
                stop_reason = message.to_owned();
                break;
            }
        }
    }
    let complete = object([
        ("schema", Json::String("A6_enrichment_v1".into())),
        ("radius_x_ticks", n(radius_ticks)),
        ("radius_y_ticks", n(64)),
        ("young_mpa", n(YOUNG_MPA)),
        ("poisson", n(POISSON)),
        ("quadrature_order", n(GAUSS_ORDER)),
        ("max_accepted_leaves", n(MAX_ACCEPTED_CELLS)),
        ("max_trial_leaves", n(MAX_TRIAL_CELLS)),
        ("max_rounds", n(rounds)),
        ("decisions", Json::Array(decisions)),
        ("stop_reason", Json::String(stop_reason)),
    ]);
    save_json(&output.join("completed.json"), &complete)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let [radius, directory, rounds] = args.as_slice() else {
        return Err("usage: enriched_bracket RADIUS_X_TICKS OUTPUT_DIRECTORY ROUNDS".into());
    };
    run(radius.parse()?, Path::new(directory), rounds.parse()?)
}
