//! Refine a loaded bracket by its local residual indicators.
#[path = "support/bracket_field.rs"]
mod bracket_field;
use solidcraft_trivariate::{
    adaptive_estimator::estimate_residual,
    adaptive_partition::{AdaptivePartition, LeafKey},
    adaptive_space::{build_adaptive_space, AdaptiveSpace},
    adaptive_transfer::build_transfer,
    analysis_space::parameterized_bracket_space,
    bracket_parameters::build_bracket_parameters,
    bracket_parameters_audit::audit_bracket_parameters,
    bracket_solve::{bracket_load, bracket_tip_displacement},
    constrained_elasticity::{solve_constrained, solve_constrained_from},
    constraint_system::{resolve, Ratio},
    elasticity::{assemble, Stiffness},
    json::{object, Json},
};
use std::{collections::BTreeMap, path::Path, time::Instant};

fn n(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn leaf(k: &LeafKey) -> Json {
    object([
        ("root", Json::String(k.root.clone())),
        ("level", n(k.level)),
        ("ijk", Json::Array(k.index.map(n).into())),
    ])
}
fn save_json(path: &Path, data: Json) -> std::io::Result<()> {
    let tmp = path.with_extension("tmp");
    std::fs::write(&tmp, format!("{data}\n"))?;
    std::fs::rename(tmp, path)
}
struct Previous {
    build: AdaptiveSpace,
    partition: AdaptivePartition,
    displacement: Vec<f64>,
    energies: [f64; 2],
    works: [f64; 2],
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn energy(k: &Stiffness, u: &[f64]) -> Result<f64, &'static str> {
    Ok(0.5 * dot(u, &k.apply(u)?))
}
fn run(rx: u8, dir: &Path, rounds: usize, warm: bool) -> Result<(), Box<dyn std::error::Error>> {
    if !(1..=6).contains(&rounds) {
        return Err("expected 1..6 adaptive rounds".into());
    }
    std::fs::create_dir_all(dir)?;
    let blocks = build_bracket_parameters(rx, 64)?.blocks;
    let report = audit_bracket_parameters(rx, 64, &blocks)?;
    let tags: BTreeMap<_, _> = report.certificate.design_faces.iter().cloned().collect();
    let parent = parameterized_bracket_space(rx, 64, &blocks)?;
    let ids = parent
        .blocks
        .iter()
        .map(|b| b.id.clone())
        .collect::<Vec<_>>();
    let mut partition = AdaptivePartition::new(&ids, &report.certificate.incidence, 5, 4000)?;
    let mut history = Vec::new();
    let mut reason = String::from("round_limit");
    let mut previous: Option<Previous> = None;
    for round in 0..rounds {
        eprintln!(
            "radius {rx} adaptive round {round}: {} leaves",
            partition.leaves().len()
        );
        let start = Instant::now();
        let build = build_adaptive_space(&parent, &partition)?;
        let prepare_s = start.elapsed().as_secs_f64();
        let start = Instant::now();
        let stiffness = assemble(&build.space, 210000., 0.3, 8)?;
        let assembly_s = start.elapsed().as_secs_f64();
        let (force, area) = bracket_load(&build.space, 8)?;
        let prescribed = build.space.boundary_dofs["y50"]
            .iter()
            .map(|&d| (d, Ratio::zero()))
            .collect::<Vec<_>>();
        let lift = resolve(build.space.ndof, &build.constraints, &prescribed)?;
        let start = Instant::now();
        let mut warm_metrics = Json::Null;
        let checks = if warm {
            Some((
                assemble(&build.space, 210000., 0.3, 6)?,
                bracket_load(&build.space, 6)?.0,
            ))
        } else {
            None
        };
        let verification_assembly_s = start.elapsed().as_secs_f64();
        let start = Instant::now();
        let (solved, solve_s) = if let Some(old) = previous.take() {
            let transfer = build_transfer(&old.build, &old.partition, &build, &partition)?;
            let initial = transfer.transfer_full_field(&old.displacement)?;
            let mut initial_reduced = Vec::with_capacity(3 * lift.free_roots.len());
            for component in 0..3 {
                initial_reduced.extend(
                    lift.free_roots
                        .iter()
                        .map(|&id| initial[3 * id + component]),
                );
            }
            let transfer_s = start.elapsed().as_secs_f64();
            let start = Instant::now();
            let hot = solve_constrained_from(
                &stiffness,
                [&lift; 3],
                &force,
                1e-10,
                20000,
                &initial_reduced,
            )?;
            let hot_s = start.elapsed().as_secs_f64();
            let start = Instant::now();
            let cold = solve_constrained(&stiffness, [&lift; 3], &force, 1e-10, 20000)?;
            let cold_s = start.elapsed().as_secs_f64();
            let (k6, f6) = checks.as_ref().ok_or("missing verification operators")?;
            let inherited = [energy(k6, &initial)?, energy(&stiffness, &initial)?];
            let loads = [dot(f6, &initial), dot(&force, &initial)];
            let error = cold
                .displacement_mm
                .iter()
                .zip(&hot.displacement_mm)
                .map(|(a, b)| (a - b) * (a - b))
                .sum::<f64>()
                .sqrt()
                / dot(&cold.displacement_mm, &cold.displacement_mm).sqrt();
            warm_metrics = object([
                ("transfer_s", n(transfer_s)),
                ("warm_s", n(hot_s)),
                ("cold_s", n(cold_s)),
                ("warm_iterations", n(hot.iterations)),
                ("cold_iterations", n(cold.iterations)),
                ("cold_true_residual", n(cold.relative_residual)),
                ("cold_warm_displacement_relative", n(error)),
                (
                    "inherited_energy_q6_q8",
                    Json::Array(inherited.map(n).into()),
                ),
                (
                    "previous_energy_q6_q8",
                    Json::Array(old.energies.map(n).into()),
                ),
                ("inherited_work_q6_q8", Json::Array(loads.map(n).into())),
                ("previous_work_q6_q8", Json::Array(old.works.map(n).into())),
            ]);
            (hot, hot_s)
        } else {
            let cold = solve_constrained(&stiffness, [&lift; 3], &force, 1e-10, 20000)?;
            (cold, start.elapsed().as_secs_f64())
        };
        let tip = bracket_tip_displacement(&build.space, &solved.displacement_mm)?[1];
        let start = Instant::now();
        let residual = estimate_residual(
            &build,
            &partition,
            &solved.displacement_mm,
            210000.,
            0.3,
            4,
            [0.; 3],
            |side, _, _| match tags[&side.root_side].as_str() {
                "y50" => None,
                "x60" => Some([0., -1000. / 320., 0.]),
                _ => Some([0.; 3]),
            },
        )?;
        let marked = residual.bulk_mark(0.5)?;
        let indicator_s = start.elapsed().as_secs_f64();
        let totals = residual.by_leaf.values().fold([0.; 3], |[v, i, b], r| {
            [
                v + r.volume_squared,
                i + r.interior_squared,
                b + r.natural_squared,
            ]
        });
        let mass: f64 = marked
            .iter()
            .map(|key| residual.by_leaf[key].total_squared)
            .sum();
        let quantities = residual
            .by_leaf
            .iter()
            .map(|(key, v)| {
                object([
                    ("leaf", leaf(key)),
                    ("volume_squared", n(v.volume_squared)),
                    ("interface_squared", n(v.interior_squared)),
                    ("natural_squared", n(v.natural_squared)),
                    ("total_squared", n(v.total_squared)),
                ])
            })
            .collect();
        save_json(
            &dir.join(format!("indicators-{round}.json")),
            Json::Array(quantities),
        )?;
        let field_path = dir.join(format!("round-{round}.bin"));
        let temp = field_path.with_extension("part");
        bracket_field::write(
            &temp,
            None,
            [rx, 64],
            8,
            &build.space,
            &solved.displacement_mm,
        )?;
        std::fs::rename(temp, field_path)?;
        let record = object([
            ("round", n(round)),
            ("leaves", n(build.space.blocks.len())),
            ("full_scalar_dofs", n(build.space.ndof)),
            ("independent_scalar_dofs", n(build.lift.free_roots.len())),
            ("free_vector_dofs", n(lift.free_roots.len() * 3)),
            ("stiffness_scalar_entries", n(stiffness.columns.len() * 9)),
            ("volume_mm3", n(stiffness.integrated_volume)),
            ("loaded_area_mm2", n(area)),
            ("tip_y_mm", n(tip)),
            ("energy_n_mm", n(solved.strain_energy_n_mm)),
            ("true_residual_relative", n(solved.relative_residual)),
            ("cg_iterations", n(solved.iterations)),
            ("warm_comparison", warm_metrics),
            ("verification_assembly_s", n(verification_assembly_s)),
            ("eta", n(residual.indicator)),
            (
                "eta_upper_bound",
                Json::Bool(residual.error_upper_bound_available),
            ),
            ("volume_indicator_squared", n(totals[0])),
            ("interface_indicator_squared", n(totals[1])),
            ("natural_indicator_squared", n(totals[2])),
            ("marked", Json::Array(marked.iter().map(leaf).collect())),
            (
                "marked_fraction_squared",
                n(mass / residual.indicator_squared),
            ),
            ("prepare_s", n(prepare_s)),
            ("assembly_s", n(assembly_s)),
            ("solve_s", n(solve_s)),
            ("indicator_s", n(indicator_s)),
        ]);
        history.push(record);
        if let Some((k6, f6)) = checks {
            let u = solved.displacement_mm;
            previous = Some(Previous {
                energies: [energy(&k6, &u)?, energy(&stiffness, &u)?],
                works: [dot(&f6, &u), dot(&force, &u)],
                build,
                partition: partition.clone(),
                displacement: u,
            });
        }
        if round + 1 < rounds {
            match partition.refine(&marked) {
                Ok(update) => {
                    if let Some(Json::Object(fields)) = history.last_mut() {
                        fields.insert(
                            "induced_refinements".into(),
                            Json::Array(update.induced.iter().map(leaf).collect()),
                        );
                    }
                }
                Err(message) => {
                    reason = message.to_string();
                    break;
                }
            }
        }
        save_json(&dir.join("completed.json"), Json::Array(history.clone()))?;
    }
    save_json(&dir.join("completed.json"), Json::Array(history))?;
    save_json(
        &dir.join("stop.json"),
        object([("reason", Json::String(reason))]),
    )?;
    Ok(())
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let (args, warm) = match args.last().map(String::as_str) {
        Some("--warm") => (&args[..args.len() - 1], true),
        _ => (args.as_slice(), false),
    };
    let [ticks, dir, rounds] = args else {
        return Err("usage: adaptive_bracket RADIUS_X_TICKS DIRECTORY ROUNDS [--warm]".into());
    };
    run(ticks.parse()?, Path::new(dir), rounds.parse()?, warm)
}
