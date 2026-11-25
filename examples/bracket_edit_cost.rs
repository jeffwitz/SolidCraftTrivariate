//! Reproduce J4c transactional mechanical-edit measurements.
use solidcraft_trivariate::bracket_mechanics_session::{
    cache_preflight, BracketMechanicsConfig, BracketMechanicsSession,
};
use solidcraft_trivariate::{
    adaptive_package::parse_refinement_plan,
    json::{self, object, Json},
    stress::fillet_limits,
};
use std::{fmt::Write, fs, path::Path, sync::Arc, time::Instant};
#[path = "support/bracket_field.rs"]
mod bracket_field;
#[path = "support/bracket_geometry.rs"]
mod bracket_geometry;

fn metrics<W: Write>(
    out: &mut W,
    state: &solidcraft_trivariate::bracket_mechanics_session::BracketMechanicsState,
) -> std::fmt::Result {
    let t = state.timing;
    write!(out,
        "{{\"generation_s\":{:.9},\"audit_and_space_s\":{:.9},\"contribution_lookup_s\":{:.9},\"local_integration_s\":{:.9},\"merge_s\":{:.9},\"comparison_s\":{:.9},\"cold_assembly_s\":{:.9},\"load_s\":{:.9},\"cold_solve_s\":{:.9},\"warm_solve_s\":{:.9},\"reused_contributions\":{},\"integrated_contributions\":{},\"cache_entries\":{},\"cache_bytes_estimated\":{},\"candidate_cells\":{},\"ndof_scalar\":{},\"cold_solve_iterations\":{},\"warm_solve_iterations\":{},\"warm_start_used\":{},\"cold_relative_residual\":{:.12e},\"warm_relative_residual\":{},\"stiffness_action_relative_error\":{:.12e},\"warm_cold_displacement_relative_error\":{},\"warm_cold_energy_relative_error\":{},\"load_area_mm2\":{:.12},\"volume_mm3\":{:.12},\"energy_n_mm\":{:.12},\"tip_displacement_mm\":[{:.12},{:.12},{:.12}]}}",
        t.generation_s, t.audit_and_space_s, t.contribution_lookup_s, t.local_integration_s,
        t.merge_s, t.comparison_s, t.cold_assembly_s, t.load_s, t.cold_solve_s, t.warm_solve_s,
        t.reused_contributions, t.integrated_contributions, t.cache_entries,
        t.cache_bytes_estimated, t.candidate_cells, state.adaptive.space.ndof,
        t.cold_solve_iterations, t.warm_solve_iterations.map_or("null".to_owned(), |x| x.to_string()),
        t.warm_start_used,
        t.cold_relative_residual,
        t.warm_relative_residual.map_or("null".to_owned(), |x| format!("{x:.12e}")),
        t.stiffness_action_relative_error,
        t.warm_cold_displacement_relative_error.map_or("null".to_owned(), |x| format!("{x:.12e}")),
        t.warm_cold_energy_relative_error.map_or("null".to_owned(), |x| format!("{x:.12e}")),
        state.load_area_mm2, state.certificate.volume_mm3, state.strain_energy_n_mm,
        state.tip_displacement_mm[0], state.tip_displacement_mm[1], state.tip_displacement_mm[2])
}

fn run_level<W: Write>(out: &mut W, level: u8) -> Result<(), Box<dyn std::error::Error>> {
    let config = BracketMechanicsConfig::campaign(level);
    let mut session = BracketMechanicsSession::new(config)?;
    writeln!(out, "{{\"level\":{level},\"transactions\":[")?;
    let sequence = [[50_u8, 74_u8], [51, 74], [81, 74], [50, 74], [50, 73]];
    let mut first = true;
    for radii in sequence {
        if !first {
            writeln!(out, ",")?;
        }
        first = false;
        let started = Instant::now();
        if radii == [81, 74] {
            let prior = Arc::clone(session.state().expect("previous accepted state"));
            let bytes = session.cache_bytes_estimated();
            let entries = session.cache_entries();
            let error = session
                .update(radii)
                .expect_err("out-of-domain radius should be rejected");
            let preserved_state = Arc::ptr_eq(session.state().unwrap(), &prior);
            let preserved_cache =
                session.cache_bytes_estimated() == bytes && session.cache_entries() == entries;
            writeln!(out, "{{\"radii_sixteenths_mm\":[81,74],\"status\":\"rejected\",\"error\":\"{error}\",\"elapsed_s\":{:.9},\"prior_state_preserved\":{preserved_state},\"cache_preserved\":{preserved_cache}}}", started.elapsed().as_secs_f64())?;
        } else {
            let state = session.update(radii)?;
            write!(out, "{{\"radii_sixteenths_mm\":[{},{}],\"status\":\"accepted\",\"elapsed_s\":{:.9},\"metrics\":", radii[0], radii[1], started.elapsed().as_secs_f64())?;
            metrics(out, state)?;
            write!(out, "}}")?;
        }
    }
    writeln!(out, "]}}")?;
    Ok(())
}

fn run_matched(
    plans: &Path,
    folder: &Path,
    incremental: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all(folder)?;
    let mut session = BracketMechanicsSession::new(BracketMechanicsConfig {
        max_cells: 3000,
        ..BracketMechanicsConfig::campaign(2)
    })?;
    let sequence = [[50, 74], [51, 74], [50, 74], [81, 74], [40, 80], [80, 40]];
    for (index, radii) in sequence.into_iter().enumerate() {
        let [rx, ry] = radii;
        let start = Instant::now();
        let plan_id = if rx == 81 { [50, 74] } else { radii };
        let source = plans.join(format!("x{}-y{}.json", plan_id[0], plan_id[1]));
        let plan = parse_refinement_plan(&json::parse(&fs::read_to_string(source)?)?)?;
        let prior = session.state().cloned();
        let bytes = session.cache_bytes_estimated();
        let entries = session.cache_entries();
        let updated = if incremental {
            session.update_incremental_with_plan(radii, &plan)
        } else {
            session.update_with_plan(radii, &plan)
        };
        let mut row = match updated {
            Ok(state) => {
                let mut encoded = String::new();
                metrics(&mut encoded, state)?;
                let mut actual = json::parse(&encoded)?;
                if let Json::Object(fields) = &mut actual {
                    fields.insert(
                        "independent_cold_comparison".into(),
                        Json::Bool(state.timing.cold_comparison_performed),
                    );
                    fields.insert(
                        "accepted_solve_s".into(),
                        Json::Number(state.timing.accepted_solve_s.to_string()),
                    );
                    fields.insert(
                        "accepted_iterations".into(),
                        Json::Number(state.timing.accepted_iterations.to_string()),
                    );
                    fields.insert(
                        "accepted_relative_residual".into(),
                        Json::Number(state.timing.accepted_relative_residual.to_string()),
                    );
                    if incremental {
                        for k in [
                            "cold_assembly_s",
                            "cold_solve_s",
                            "cold_solve_iterations",
                            "cold_relative_residual",
                            "comparison_s",
                            "stiffness_action_relative_error",
                            "warm_cold_displacement_relative_error",
                            "warm_cold_energy_relative_error",
                        ] {
                            fields.remove(k);
                        }
                    }
                }
                let limits = fillet_limits(&state.adaptive.space, &state.displacement_mm)?;
                let field = folder.join(format!("event-{index}.bin"));
                bracket_field::write(
                    &field,
                    None,
                    radii,
                    8,
                    &state.adaptive.space,
                    &state.displacement_mm,
                )?;
                bracket_geometry::write(
                    &field.with_extension("geometry.bin"),
                    &state.adaptive.space,
                )?;
                object([
                    ("status", Json::String("accepted".into())),
                    ("metrics", actual),
                    (
                        "reaction_n",
                        Json::Array(
                            state
                                .support_reaction_n
                                .map(|v| Json::Number(v.to_string()))
                                .to_vec(),
                        ),
                    ),
                    (
                        "fillet_mpa",
                        Json::Object(
                            limits
                                .into_iter()
                                .map(|(k, v)| (k, Json::Number(v.to_string())))
                                .collect(),
                        ),
                    ),
                ])
            }
            Err(e) => {
                let preserved = match (&prior, session.state()) {
                    (Some(old), Some(new)) => Arc::ptr_eq(old, new),
                    (None, None) => true,
                    _ => false,
                };
                object([
                    ("status", Json::String("refused".into())),
                    ("diagnostic", Json::String(e.into())),
                    ("prior_state_preserved", Json::Bool(preserved)),
                    (
                        "cache_preserved",
                        Json::Bool(
                            session.cache_bytes_estimated() == bytes
                                && session.cache_entries() == entries,
                        ),
                    ),
                ])
            }
        };
        if let Json::Object(fields) = &mut row {
            fields.insert(
                "elapsed_with_checkpoint_s".into(),
                Json::Number(start.elapsed().as_secs_f64().to_string()),
            );
            fields.insert(
                "radii_ticks".into(),
                Json::Array(radii.map(|v| Json::Number(v.to_string())).to_vec()),
            );
        }
        fs::write(
            folder.join(format!("event-{index}.json")),
            format!("{row}\n"),
        )?;
        println!("event {index}, radii {rx}/{ry}");
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if let [mode, plans, folder] = args.as_slice() {
        if mode == "--verified-plan" || mode == "--incremental-plan" {
            return run_matched(
                Path::new(plans),
                Path::new(folder),
                mode == "--incremental-plan",
            );
        }
    }
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "docs/results/bracket-edit-cost-20261125.json".to_owned());
    let preflight = cache_preflight(3, 402_653_184)?;
    let mut levels = Vec::new();
    for level in 0..=2 {
        let mut result = String::new();
        run_level(&mut result, level)?;
        levels.push(result);
        let mut body = String::new();
        writeln!(body, "{{\"protocol\":\"bracket_mechanics_session_cost_v1\",\"base_revision\":\"fdc8653\",\"architecture\":\"{}-{}\",\"threads\":{},\"material\":{{\"young_mpa\":210000.0,\"poisson\":0.3}},\"quadrature_order\":8,\"cache_cap_bytes\":402653184,\"cache_accounting\":\"vector capacities plus key/control payload and 512 bytes per entry; excludes allocator metadata\",\"levels\":[{}],\"level_3_preflight\":{{\"cells\":{},\"bytes_per_entry_upper_bound\":{},\"total_bytes_upper_bound\":{},\"cap_bytes\":{},\"fits\":{},\"executed\":false}},\"oracle_thresholds\":{{\"stiffness_action_relative\":1e-11,\"warm_cold_solution_relative\":1e-8,\"true_relative_residual\":1e-10}}}}",
            std::env::consts::ARCH,
            std::env::consts::OS,
            std::thread::available_parallelism().map_or(1, usize::from),
            levels.join(","),
            preflight.cells,
            preflight.bytes_per_entry_upper_bound,
            preflight.total_bytes_upper_bound,
            preflight.cap_bytes,
            preflight.fits)?;
        fs::write(&path, body)?;
    }
    println!("wrote {path}");
    Ok(())
}
