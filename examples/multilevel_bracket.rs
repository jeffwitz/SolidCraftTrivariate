//! Measure additive geometric multilevel preconditioning on bracket spaces.
use solidcraft_trivariate::{
    adaptive_partition::{AdaptivePartition, LeafKey},
    adaptive_space::{build_adaptive_space, AdaptiveSpace},
    adaptive_transfer::build_transfer,
    analysis_space::{parameterized_bracket_space, AnalysisSpace},
    bracket_parameters::build_bracket_parameters,
    bracket_parameters_audit::audit_bracket_parameters,
    bracket_solve::{bracket_load, bracket_tip_displacement},
    constrained_elasticity::{solve_constrained, solve_preconditioned, ReducedPreconditioner},
    constraint_system::{resolve, Lift, Ratio},
    elasticity::{assemble, Stiffness},
    incidence::Incidence,
    json::{object, Json},
    multilevel::{MountedSpace, MultilevelPreconditioner},
};
use std::{
    cell::Cell,
    env, fs,
    mem::size_of,
    time::{Duration, Instant},
};

const YOUNG_MPA: f64 = 210_000.0;
const POISSON: f64 = 0.3;
const ORDER: usize = 8;
const TOLERANCE: f64 = 1e-10;
const MAX_ITERATIONS: usize = 20_000;
const MEMORY_CAP: usize = 6 * 1024 * 1024 * 1024;
const TIME_CAP_S: f64 = 1_800.0;
const COARSE_FACTOR_CAP: usize = MEMORY_CAP;

fn text(value: impl Into<String>) -> Json {
    Json::String(value.into())
}
fn number(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}
fn boolean(value: bool) -> Json {
    Json::Bool(value)
}
fn array(values: impl IntoIterator<Item = Json>) -> Json {
    Json::Array(values.into_iter().collect())
}
fn f64_array(values: impl IntoIterator<Item = f64>) -> Json {
    array(values.into_iter().map(number))
}

struct MeshLevel {
    name: String,
    partition: AdaptivePartition,
    adaptive: AdaptiveSpace,
    lift: Lift,
    stiffness: Stiffness,
    refinement_s: f64,
    space_s: f64,
    assembly_s: f64,
}

struct CaseTiming {
    root_setup_s: f64,
    cumulative_refinement_s: f64,
    cumulative_mesh_prep_s: f64,
    cumulative_assembly_s: f64,
    cumulative_transfer_s: f64,
    factor_s: f64,
    target_preflight_bytes: Option<usize>,
    transfer_adjoint: Option<(f64, f64)>,
}

struct TimedPreconditioner<'a> {
    inner: &'a MultilevelPreconditioner,
    count: Cell<usize>,
    elapsed: Cell<Duration>,
}

impl<'a> TimedPreconditioner<'a> {
    fn new(inner: &'a MultilevelPreconditioner) -> Self {
        Self {
            inner,
            count: Cell::new(0),
            elapsed: Cell::new(Duration::ZERO),
        }
    }
    fn reset(&self) {
        self.count.set(0);
        self.elapsed.set(Duration::ZERO);
    }
    fn measurements(&self) -> (usize, f64) {
        (self.count.get(), self.elapsed.get().as_secs_f64())
    }
}

impl ReducedPreconditioner for TimedPreconditioner<'_> {
    fn apply(&self, residual: &[f64]) -> Result<Vec<f64>, &'static str> {
        let started = Instant::now();
        let result = self.inner.apply(residual);
        self.count.set(self.count.get() + 1);
        self.elapsed.set(self.elapsed.get() + started.elapsed());
        result
    }
}

fn mounted_lift(adaptive: &AdaptiveSpace) -> Result<Lift, &'static str> {
    let prescribed = adaptive
        .space
        .boundary_dofs
        .get("y50")
        .ok_or("missing y50 clamp")?
        .iter()
        .map(|&dof| (dof, Ratio::zero()))
        .collect::<Vec<_>>();
    resolve(adaptive.space.ndof, &adaptive.constraints, &prescribed)
}

fn build_level(
    name: &str,
    parent: &AnalysisSpace,
    partition: AdaptivePartition,
    refinement_s: f64,
) -> Result<MeshLevel, &'static str> {
    let started = Instant::now();
    let adaptive = build_adaptive_space(parent, &partition)?;
    let space_s = started.elapsed().as_secs_f64();
    let lift = mounted_lift(&adaptive)?;
    let started = Instant::now();
    let stiffness = assemble(&adaptive.space, YOUNG_MPA, POISSON, ORDER)?;
    let assembly_s = started.elapsed().as_secs_f64();
    Ok(MeshLevel {
        name: name.to_owned(),
        partition,
        adaptive,
        lift,
        stiffness,
        refinement_s,
        space_s,
        assembly_s,
    })
}

fn refine_all(source: &AdaptivePartition) -> Result<(AdaptivePartition, f64), &'static str> {
    let started = Instant::now();
    let mut target = source.clone();
    let leaves = target.snapshot()?.leaves;
    target.refine(&leaves)?;
    Ok((target, started.elapsed().as_secs_f64()))
}

fn stiffness_bytes(stiffness: &Stiffness) -> usize {
    stiffness.offsets.capacity() * size_of::<usize>()
        + stiffness.columns.capacity() * size_of::<usize>()
        + stiffness.values.capacity() * size_of::<[[f64; 3]; 3]>()
}

fn space_bytes(adaptive: &AdaptiveSpace) -> usize {
    let blocks = adaptive
        .space
        .blocks
        .iter()
        .map(|block| {
            block.id.capacity()
                + block.controls.capacity() * size_of::<[i128; 4]>()
                + block.dofs.capacity() * size_of::<usize>()
        })
        .sum::<usize>();
    let lift = adaptive.lift.free_roots.capacity() * size_of::<usize>()
        + adaptive
            .lift
            .lines
            .iter()
            .map(|line| line.terms.capacity() * size_of::<(usize, Ratio)>() + size_of::<Ratio>())
            .sum::<usize>();
    blocks + lift
}

fn preconditioner_bytes(mg: &MultilevelPreconditioner, coarse: &Lift, fine: &Lift) -> usize {
    let dimension = 3 * coarse.free_roots.len();
    let factor = dimension * dimension * size_of::<f64>() + dimension * size_of::<usize>();
    let transfer = mg.transfer_nonzeros() * size_of::<(usize, f64)>()
        + 3 * fine.free_roots.len() * size_of::<f64>();
    factor + transfer
}

fn process_high_water_kib() -> Option<u64> {
    let status = fs::read_to_string("/proc/self/status").ok()?;
    status.lines().find_map(|line| {
        line.strip_prefix("VmHWM:")
            .and_then(|value| value.split_whitespace().next())
            .and_then(|value| value.parse().ok())
    })
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}
fn relative_difference(a: &[f64], b: &[f64]) -> f64 {
    let difference = a
        .iter()
        .zip(b)
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt();
    difference / norm(b).max(f64::MIN_POSITIVE)
}
fn deterministic(count: usize, phase: usize) -> Vec<f64> {
    (0..count)
        .map(|i| (((i * (37 + phase * 2) + 19) % 1009) as f64 - 504.) / 239.)
        .collect()
}

fn probe_preconditioner(
    mg: &MultilevelPreconditioner,
    free_roots: usize,
) -> Result<Json, &'static str> {
    let n = 3 * free_roots;
    let x = deterministic(n, 0);
    let y = deterministic(n, 1);
    let timed = TimedPreconditioner::new(mg);
    let bx = timed.apply(&x)?;
    let by = timed.apply(&y)?;
    let symmetry = (dot(&x, &by) - dot(&y, &bx)).abs()
        / (norm(&x) * norm(&by) + norm(&y) * norm(&bx)).max(f64::MIN_POSITIVE);
    let x_positive = dot(&x, &bx);
    let y_positive = dot(&y, &by);
    let (applications, apply_s) = timed.measurements();
    if symmetry > 1e-10
        || !x_positive.is_finite()
        || !y_positive.is_finite()
        || x_positive <= 0.0
        || y_positive <= 0.0
    {
        return Err("multilevel symmetry or positive quadratic check failed");
    }
    Ok(object([
        ("symmetry_relative", number(symmetry)),
        ("x_B_x", number(x_positive)),
        ("y_B_y", number(y_positive)),
        ("applications", number(applications)),
        ("application_s", number(apply_s)),
    ]))
}

fn mounted_transfer_adjoint_error(
    source: &MeshLevel,
    target: &MeshLevel,
) -> Result<f64, &'static str> {
    let transfer = build_transfer(
        &source.adaptive,
        &source.partition,
        &target.adaptive,
        &target.partition,
    )?;
    let source_map = source
        .adaptive
        .lift
        .free_roots
        .iter()
        .map(|&dof| source.lift.free_roots.binary_search(&dof).ok())
        .collect::<Vec<_>>();
    let target_map = target
        .adaptive
        .lift
        .free_roots
        .iter()
        .map(|&dof| target.lift.free_roots.binary_search(&dof).ok())
        .collect::<Vec<_>>();
    let mut rows = vec![Vec::<(usize, f64)>::new(); target.lift.free_roots.len()];
    for (row, target_index) in transfer.rows.iter().zip(target_map) {
        if let Some(target_index) = target_index {
            for &(source_index, weight) in row {
                let coefficient = weight.numerator() as f64 / weight.denominator() as f64;
                if let Some(source_index) = source_map[source_index] {
                    rows[target_index].push((source_index, coefficient));
                } else if coefficient != 0.0 {
                    let source_dof = source.adaptive.lift.free_roots[source_index];
                    let line = &source.lift.lines[source_dof];
                    if !line.terms.is_empty() || line.offset.numerator() != 0 {
                        return Err("mounted source root is not fixed to homogeneous zero");
                    }
                }
            }
        }
    }
    let ns = source.lift.free_roots.len();
    let nt = target.lift.free_roots.len();
    let x = deterministic(3 * nt, 3);
    let y = deterministic(3 * ns, 4);
    let mut py = vec![0.0; 3 * nt];
    let mut ptx = vec![0.0; 3 * ns];
    for component in 0..3 {
        for (i, row) in rows.iter().enumerate() {
            for &(j, weight) in row {
                py[component * nt + i] += weight * y[component * ns + j];
                ptx[component * ns + j] += weight * x[component * nt + i];
            }
        }
    }
    Ok((dot(&x, &py) - dot(&ptx, &y)).abs()
        / (norm(&x) * norm(&py) + norm(&y) * norm(&ptx)).max(f64::MIN_POSITIVE))
}

fn case_record(
    rx: u8,
    path: &str,
    level: &MeshLevel,
    mg: &MultilevelPreconditioner,
    coarse_lift: &Lift,
    timing: CaseTiming,
) -> Result<Json, &'static str> {
    let load_started = Instant::now();
    let (force_n, area) = bracket_load(&level.adaptive.space, ORDER)?;
    let load_s = load_started.elapsed().as_secs_f64();
    if (area - 320.0).abs() > 1e-8 {
        return Err("loaded x60 area differs from 320 mm^2");
    }
    let memory_estimate = preconditioner_bytes(mg, coarse_lift, &level.lift)
        + stiffness_bytes(&level.stiffness)
        + space_bytes(&level.adaptive)
        + 18 * 3 * level.adaptive.space.ndof * size_of::<f64>();
    if memory_estimate > MEMORY_CAP {
        return Err("estimated multilevel resolution memory exceeds 6 GiB");
    }
    let lifts = [&level.lift, &level.lift, &level.lift];
    let started = Instant::now();
    let jacobi = solve_constrained(&level.stiffness, lifts, &force_n, TOLERANCE, MAX_ITERATIONS);
    let jacobi_s = started.elapsed().as_secs_f64();

    let probe = probe_preconditioner(mg, level.lift.free_roots.len())?;
    let timed = TimedPreconditioner::new(mg);
    timed.reset();
    let started = Instant::now();
    let multilevel = solve_preconditioned(
        &level.stiffness,
        lifts,
        &force_n,
        TOLERANCE,
        MAX_ITERATIONS,
        &timed,
    );
    let multilevel_s = started.elapsed().as_secs_f64();
    let (applications, application_s) = timed.measurements();

    let jacobi_record = match &jacobi {
        Ok(result) => object([
            ("iterations", number(result.iterations)),
            ("relative_residual", number(result.relative_residual)),
            ("energy_n_mm", number(result.strain_energy_n_mm)),
            (
                "tip_displacement_mm",
                f64_array(bracket_tip_displacement(
                    &level.adaptive.space,
                    &result.displacement_mm,
                )?),
            ),
        ]),
        Err(error) => object([("error", text(*error))]),
    };
    let multilevel_record = match &multilevel {
        Ok(result) => object([
            ("iterations", number(result.iterations)),
            ("relative_residual", number(result.relative_residual)),
            ("energy_n_mm", number(result.strain_energy_n_mm)),
            (
                "tip_displacement_mm",
                f64_array(bracket_tip_displacement(
                    &level.adaptive.space,
                    &result.displacement_mm,
                )?),
            ),
        ]),
        Err(error) => object([("error", text(*error))]),
    };
    let (field_error, energy_error, iteration_ratio, residual_ok) = match (&jacobi, &multilevel) {
        (Ok(jacobi), Ok(multilevel)) => (
            Some(relative_difference(
                &multilevel.displacement_mm,
                &jacobi.displacement_mm,
            )),
            Some(
                (multilevel.strain_energy_n_mm - jacobi.strain_energy_n_mm).abs()
                    / jacobi.strain_energy_n_mm.abs().max(f64::MIN_POSITIVE),
            ),
            Some(jacobi.iterations as f64 / multilevel.iterations.max(1) as f64),
            jacobi.relative_residual <= 1e-10 && multilevel.relative_residual <= 1e-10,
        ),
        _ => (None, None, None, false),
    };
    let jacobi_total = timing.root_setup_s
        + timing.cumulative_refinement_s
        + level.space_s
        + level.assembly_s
        + load_s
        + jacobi_s;
    let multilevel_total = timing.root_setup_s
        + timing.cumulative_mesh_prep_s
        + timing.cumulative_assembly_s
        + timing.factor_s
        + timing.cumulative_transfer_s
        + load_s
        + multilevel_s;
    let hlevel = level
        .name
        .strip_prefix('h')
        .and_then(|x| x.parse::<u8>().ok());
    let iteration_target_applies = hlevel.is_some_and(|h| h >= 2);
    let elapsed_limit_ok = jacobi_s <= TIME_CAP_S && multilevel_s <= TIME_CAP_S;
    let field_and_energy_ok =
        field_error.is_some_and(|x| x <= 1e-8) && energy_error.is_some_and(|x| x <= 1e-8);
    let case_status = if residual_ok && field_and_energy_ok && elapsed_limit_ok {
        "qualified"
    } else if !elapsed_limit_ok {
        "over_time_budget"
    } else if jacobi.is_err() || multilevel.is_err() {
        "solver_failed"
    } else {
        "comparison_failed"
    };
    let mut entries = vec![
        (
            "radii_sixteenths_mm".to_owned(),
            array([number(rx), number(64)]),
        ),
        ("path".to_owned(), text(path)),
        ("level".to_owned(), text(&level.name)),
        ("status".to_owned(), text(case_status)),
        (
            "candidate_cells".to_owned(),
            number(level.adaptive.space.blocks.len()),
        ),
        ("scalar_dofs".to_owned(), number(level.adaptive.space.ndof)),
        (
            "mounted_free_scalar_masters".to_owned(),
            number(level.lift.free_roots.len()),
        ),
        ("load_area_mm2".to_owned(), number(area)),
        (
            "geometry_audit_space_s".to_owned(),
            number(timing.root_setup_s),
        ),
        ("refinement_plan_s".to_owned(), number(level.refinement_s)),
        ("analysis_space_s".to_owned(), number(level.space_s)),
        ("assembly_current_s".to_owned(), number(level.assembly_s)),
        (
            "assembly_hierarchy_cumulative_s".to_owned(),
            number(timing.cumulative_assembly_s),
        ),
        ("coarse_factor_s".to_owned(), number(timing.factor_s)),
        (
            "transfer_setup_cumulative_s".to_owned(),
            number(timing.cumulative_transfer_s),
        ),
        (
            "transfer_adjoint_relative_error".to_owned(),
            timing.transfer_adjoint.map_or(Json::Null, |x| number(x.0)),
        ),
        (
            "transfer_adjoint_oracle_s".to_owned(),
            timing.transfer_adjoint.map_or(Json::Null, |x| number(x.1)),
        ),
        ("operator_probe".to_owned(), probe),
        ("jacobi_solve_s".to_owned(), number(jacobi_s)),
        ("jacobi".to_owned(), jacobi_record),
        ("multilevel_solve_s".to_owned(), number(multilevel_s)),
        ("multilevel_apply_count".to_owned(), number(applications)),
        ("multilevel_apply_s".to_owned(), number(application_s)),
        ("multilevel".to_owned(), multilevel_record),
        (
            "warm_cold_displacement_relative_error".to_owned(),
            field_error.map_or(Json::Null, number),
        ),
        (
            "warm_cold_energy_relative_error".to_owned(),
            energy_error.map_or(Json::Null, number),
        ),
        ("true_residual_within_1e10".to_owned(), boolean(residual_ok)),
        (
            "field_and_energy_within_1e8".to_owned(),
            boolean(field_and_energy_ok),
        ),
        (
            "diagonal_to_multilevel_iteration_ratio".to_owned(),
            iteration_ratio.map_or(Json::Null, number),
        ),
        (
            "iteration_target_applicable".to_owned(),
            boolean(iteration_target_applies),
        ),
        (
            "iteration_target_met".to_owned(),
            if iteration_target_applies {
                iteration_ratio.map_or(Json::Null, |x| boolean(x >= 2.0))
            } else {
                Json::Null
            },
        ),
        ("jacobi_total_one_off_s".to_owned(), number(jacobi_total)),
        (
            "multilevel_total_one_off_s".to_owned(),
            number(multilevel_total),
        ),
        (
            "per_resolution_time_within_1800s".to_owned(),
            boolean(elapsed_limit_ok),
        ),
        ("memory_cap_bytes".to_owned(), number(MEMORY_CAP)),
        (
            "target_preflight_upper_bound_bytes".to_owned(),
            timing.target_preflight_bytes.map_or(Json::Null, number),
        ),
        (
            "hierarchy_memory_estimate_bytes".to_owned(),
            number(memory_estimate),
        ),
        (
            "process_high_water_kib".to_owned(),
            process_high_water_kib().map_or(Json::Null, number),
        ),
    ];
    let estimated = number(memory_estimate);
    entries.push(("preparation_included_in_total".to_owned(), boolean(true)));
    entries.push(("memory_basis".to_owned(), text("CSR and Q2 Vec capacities, lift rows, coarse factor, transfers and a conservative PCG-vector allowance; process VmHWM is also recorded")));
    let mut map = std::collections::BTreeMap::new();
    for (key, value) in entries {
        map.insert(key, value);
    }
    map.insert("estimated_bytes_check".to_owned(), estimated);
    Ok(Json::Object(map))
}

fn root_partition(
    ids: &[String],
    incidence: &Incidence,
) -> Result<AdaptivePartition, &'static str> {
    AdaptivePartition::new(ids, incidence, 3, 18_432)
}

fn preflight_target_memory(
    target_cells: usize,
    previous: &MeshLevel,
    mg: &MultilevelPreconditioner,
    coarse: &Lift,
) -> Result<usize, &'static str> {
    let target_bound = target_cells
        .checked_mul(128 * 1024)
        .ok_or("target memory estimate overflow")?;
    let retained_bound = preconditioner_bytes(mg, coarse, &previous.lift)
        + stiffness_bytes(&previous.stiffness)
        + 2 * space_bytes(&previous.adaptive)
        + 18 * 3 * previous.adaptive.space.ndof * size_of::<f64>();
    let total = target_bound
        .checked_add(retained_bound)
        .ok_or("target memory estimate overflow")?;
    if total > MEMORY_CAP {
        return Err("next level preflight exceeds the 6 GiB resolution cap");
    }
    Ok(total)
}

struct PathInput<'a> {
    rx: u8,
    path_name: &'a str,
    root_setup_s: f64,
    parent: &'a AnalysisSpace,
    incidence: &'a Incidence,
    output_path: &'a str,
}

fn run_path(
    input: PathInput<'_>,
    output: &mut Vec<Json>,
    errors: &mut Vec<Json>,
) -> Result<(), &'static str> {
    let PathInput {
        rx,
        path_name,
        root_setup_s,
        parent,
        incidence,
        output_path,
    } = input;
    let ids = parent
        .blocks
        .iter()
        .map(|block| block.id.clone())
        .collect::<Vec<_>>();
    let started = Instant::now();
    let p0 = root_partition(&ids, incidence)?;
    let partition_s = started.elapsed().as_secs_f64();
    let mut previous = build_level("h0", parent, p0, partition_s)?;
    let factor_start = Instant::now();
    let mut mg =
        MultilevelPreconditioner::new(&previous.stiffness, &previous.lift, COARSE_FACTOR_CAP)?;
    let factor_s = factor_start.elapsed().as_secs_f64();
    let coarse_lift = previous.lift.clone();
    let mut cumulative_refinement = previous.refinement_s;
    let mut cumulative_mesh = previous.refinement_s + previous.space_s;
    let mut cumulative_assembly = previous.assembly_s;
    let mut cumulative_transfer = 0.0;
    output.push(case_record(
        rx,
        path_name,
        &previous,
        &mg,
        &coarse_lift,
        CaseTiming {
            root_setup_s,
            cumulative_refinement_s: cumulative_refinement,
            cumulative_mesh_prep_s: cumulative_mesh,
            cumulative_assembly_s: cumulative_assembly,
            cumulative_transfer_s: cumulative_transfer,
            factor_s,
            target_preflight_bytes: None,
            transfer_adjoint: None,
        },
    )?);
    persist(output_path, output, errors).map_err(|_| "cannot persist campaign progress")?;

    let count = if path_name == "uniform" { 3 } else { 2 };
    for next in 1..=count {
        let (partition, refinement_s) = if path_name == "uniform" {
            refine_all(&previous.partition)?
        } else if next == 1 {
            let started = Instant::now();
            let mut partition = previous.partition.clone();
            partition.refine(&[LeafKey::root("hole_x_sector_east_z10_30")])?;
            (partition, started.elapsed().as_secs_f64())
        } else {
            let started = Instant::now();
            let mut partition = previous.partition.clone();
            partition.refine(&[LeafKey {
                root: "hole_x_sector_east_z10_30".to_owned(),
                level: 1,
                index: [1, 0, 0],
            }])?;
            (partition, started.elapsed().as_secs_f64())
        };
        let target_cells = partition.snapshot()?.leaves.len();
        let preflight_bytes =
            match preflight_target_memory(target_cells, &previous, &mg, &coarse_lift) {
                Ok(bytes) => bytes,
                Err(error) => {
                    errors.push(object([
                        ("radius_x_sixteenths_mm", number(rx)),
                        ("path", text(path_name)),
                        (
                            "level",
                            text(format!(
                                "{}{}",
                                if path_name == "uniform" { "h" } else { "local" },
                                next
                            )),
                        ),
                        ("candidate_cells", number(target_cells)),
                        ("status", text("preflight_refused")),
                        ("error", text(error)),
                    ]));
                    persist(output_path, output, errors)
                        .map_err(|_| "cannot persist campaign progress")?;
                    return Ok(());
                }
            };
        let name = if path_name == "uniform" {
            format!("h{next}")
        } else {
            format!("local{next}")
        };
        let target = build_level(&name, parent, partition, refinement_s)?;
        let transfer_started = Instant::now();
        let transfer_adjoint = mounted_transfer_adjoint_error(&previous, &target)?;
        let transfer_oracle_s = transfer_started.elapsed().as_secs_f64();
        if transfer_adjoint > 1e-11 {
            return Err("mounted transfer adjoint check exceeded tolerance");
        }
        let setup_started = Instant::now();
        mg.add_level(
            MountedSpace {
                adaptive: &previous.adaptive,
                partition: &previous.partition,
                lift: &previous.lift,
            },
            MountedSpace {
                adaptive: &target.adaptive,
                partition: &target.partition,
                lift: &target.lift,
            },
            &target.stiffness,
        )?;
        let transfer_setup = setup_started.elapsed().as_secs_f64();
        cumulative_transfer += transfer_setup;
        cumulative_refinement += target.refinement_s;
        cumulative_mesh += target.refinement_s + target.space_s;
        cumulative_assembly += target.assembly_s;
        output.push(case_record(
            rx,
            path_name,
            &target,
            &mg,
            &coarse_lift,
            CaseTiming {
                root_setup_s,
                cumulative_refinement_s: cumulative_refinement,
                cumulative_mesh_prep_s: cumulative_mesh,
                cumulative_assembly_s: cumulative_assembly,
                cumulative_transfer_s: cumulative_transfer,
                factor_s,
                target_preflight_bytes: Some(preflight_bytes),
                transfer_adjoint: Some((transfer_adjoint, transfer_oracle_s)),
            },
        )?);
        persist(output_path, output, errors).map_err(|_| "cannot persist campaign progress")?;
        previous = target;
    }
    Ok(())
}

fn persist(path: &str, records: &[Json], errors: &[Json]) -> std::io::Result<()> {
    fs::write(
        path,
        format!("{}\n", report(records.to_vec(), errors.to_vec())),
    )
}

fn report(records: Vec<Json>, errors: Vec<Json>) -> Json {
    object([
        ("protocol", text("bracket_multilevel_additive_cost_v1")),
        ("base_revision", text("2e0b22c")),
        (
            "architecture",
            text(format!("{}-{}", env::consts::ARCH, env::consts::OS)),
        ),
        (
            "threads",
            number(std::thread::available_parallelism().map_or(1, usize::from)),
        ),
        (
            "material",
            object([
                ("young_mpa", number(YOUNG_MPA)),
                ("poisson", number(POISSON)),
            ]),
        ),
        ("quadrature_order", number(ORDER)),
        ("tolerance", number(TOLERANCE)),
        ("per_resolution_memory_cap_bytes", number(MEMORY_CAP)),
        ("per_resolution_time_cap_s", number(TIME_CAP_S)),
        ("cases", array(records)),
        ("campaign_errors", array(errors)),
        (
            "process_high_water_kib",
            process_high_water_kib().map_or(Json::Null, number),
        ),
    ])
}

fn run_radius(
    rx: u8,
    records: &mut Vec<Json>,
    errors: &mut Vec<Json>,
    output_path: &str,
) -> Result<(), &'static str> {
    let started = Instant::now();
    let geometry = build_bracket_parameters(rx, 64)?;
    let certificate = audit_bracket_parameters(rx, 64, &geometry.blocks)?;
    let parent = parameterized_bracket_space(rx, 64, &geometry.blocks)?;
    let root_setup_s = started.elapsed().as_secs_f64();
    for path in ["uniform", "local"] {
        if let Err(error) = run_path(
            PathInput {
                rx,
                path_name: path,
                root_setup_s,
                parent: &parent,
                incidence: &certificate.certificate.incidence,
                output_path,
            },
            records,
            errors,
        ) {
            errors.push(object([
                ("radius_x_sixteenths_mm", number(rx)),
                ("path", text(path)),
                ("error", text(error)),
            ]));
            persist(output_path, records, errors).map_err(|_| "cannot persist campaign error")?;
        }
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output_path = env::args()
        .nth(1)
        .unwrap_or_else(|| "docs/results/multilevel-bracket-20261125.json".to_owned());
    let mut records = Vec::new();
    let mut errors = Vec::new();
    for rx in [64_u8, 56_u8] {
        if let Err(error) = run_radius(rx, &mut records, &mut errors, &output_path) {
            errors.push(object([
                ("radius_x_sixteenths_mm", number(rx)),
                ("error", text(error)),
            ]));
        }
        persist(&output_path, &records, &errors)?;
    }
    println!("wrote {output_path}");
    Ok(())
}
