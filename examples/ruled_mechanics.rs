//! Manufactured elasticity qualification for varying-section ruled sweeps.
use solidcraft_trivariate::{
    analysis_space::{ruled_sweep_analysis_space, AnalysisBlock, AnalysisSpace},
    boundary_traction::{integrate_traction, LoadPatch, SurfaceLoad},
    constrained_elasticity::solve_constrained,
    constraint_system::{resolve, Ratio},
    elasticity::{assemble, basis_geometry, gauss, BasisSample, Vector},
    h_refine::refine_uniform,
    incidence::BlockSide,
    json::{object, parse, Json},
    ruled_sweep::{build_ruled_sweep, SweepRecipe, SweepSection},
    ruled_sweep_audit::{audit_ruled_sweep, RuledSweepCertificate},
};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    time::Instant,
};

const YOUNG: f64 = 210_000.;
const POISSON: f64 = 0.3;
const LAMBDA: f64 = YOUNG * POISSON / ((1. + POISSON) * (1. - 2. * POISSON));
const MU: f64 = YOUNG / (2. * (1. + POISSON));
const C: f64 = LAMBDA + 2. * MU;
const SECTION_AREA: f64 = 400. - 16. * std::f64::consts::PI;

fn num(x: f64) -> Json {
    Json::Number(format!("{x:.12e}"))
}
fn integer(x: usize) -> Json {
    Json::Number(x.to_string())
}
fn dot(a: Vector, b: Vector) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: Vector, b: Vector) -> Vector {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: Vector) -> f64 {
    dot(a, a).sqrt()
}

#[derive(Clone)]
struct Program {
    recipe: SweepRecipe,
    base: AnalysisSpace,
    groups: BTreeMap<String, Vec<BlockSide>>,
    certificate: RuledSweepCertificate,
}

fn groups_for(
    recipe: SweepRecipe,
    certificate: &RuledSweepCertificate,
) -> Result<BTreeMap<String, Vec<BlockSide>>, &'static str> {
    let mut groups = BTreeMap::<String, Vec<BlockSide>>::new();
    for (side, tag) in &certificate.certificate.design_faces {
        let name = if side.axis == 2 {
            match tag.as_str() {
                "z0" => "bottom",
                "z8" => "top",
                _ => return Err("unexpected ruled cap selection"),
            }
        } else {
            tag.as_str()
        };
        groups
            .entry(name.to_owned())
            .or_default()
            .push(side.clone());
    }
    let expected = if recipe.section == SweepSection::PiercedSquare {
        7
    } else {
        6
    };
    if groups.len() != expected || !groups.contains_key("bottom") || !groups.contains_key("top") {
        return Err("ruled sweep selections are incomplete");
    }
    Ok(groups)
}

fn make_program(recipe: SweepRecipe) -> Result<Program, &'static str> {
    let blocks = build_ruled_sweep(recipe)?;
    let certificate = audit_ruled_sweep(recipe, &blocks)?;
    let base = ruled_sweep_analysis_space(recipe, &blocks)?;
    let groups = groups_for(recipe, &certificate)?;
    if base.blocks.len() != blocks.len()
        || base.boundary_dofs.keys().collect::<Vec<_>>() != groups.keys().collect::<Vec<_>>()
    {
        return Err("analysis space lost a certified feature");
    }
    Ok(Program {
        recipe,
        base,
        groups,
        certificate,
    })
}

fn refine(program: &Program, level: usize) -> Result<AnalysisSpace, &'static str> {
    let mut space = program.base.clone();
    for _ in 0..level {
        space = refine_uniform(&space)?;
    }
    Ok(space)
}

fn refined_groups(
    program: &Program,
    space: &AnalysisSpace,
) -> Result<BTreeMap<String, Vec<BlockSide>>, &'static str> {
    let mut result = BTreeMap::<String, Vec<BlockSide>>::new();
    for (group, sides) in &program.groups {
        for base in sides {
            for block in &space.blocks {
                let (root, path) = block
                    .id
                    .split_once('/')
                    .map_or((block.id.as_str(), ""), |(a, b)| (a, b));
                if root != base.block {
                    continue;
                }
                let bit = u8::from(base.upper);
                let on_boundary = path
                    .split('/')
                    .filter(|p| !p.is_empty())
                    .all(|p| p.as_bytes().get(base.axis).copied() == Some(b'0' + bit));
                if on_boundary {
                    result.entry(group.clone()).or_default().push(BlockSide {
                        block: block.id.clone(),
                        axis: base.axis,
                        upper: base.upper,
                    });
                }
            }
        }
    }
    for sides in result.values_mut() {
        sides.sort();
        sides.dedup();
    }
    if result.len() != program.groups.len() || result.values().any(Vec::is_empty) {
        return Err("failed to propagate ruled exterior selections");
    }
    Ok(result)
}

fn stress_from_axial_strain(e: f64) -> [[f64; 3]; 3] {
    [[LAMBDA * e, 0., 0.], [0., LAMBDA * e, 0.], [0., 0., C * e]]
}
fn traction_at(x: Vector, n: Vector, power: u32, coefficient: f64) -> Vector {
    let e = coefficient * power as f64 * x[2].powi(power as i32 - 1);
    let sigma = stress_from_axial_strain(e);
    [dot(sigma[0], n), dot(sigma[1], n), dot(sigma[2], n)]
}

fn integrate_surface_groups(
    space: &AnalysisSpace,
    groups: &BTreeMap<String, Vec<BlockSide>>,
    power: u32,
    coefficient: f64,
    order: usize,
) -> Result<BTreeMap<String, SurfaceLoad>, &'static str> {
    let exterior: BTreeSet<_> = groups.values().flatten().cloned().collect();
    let mut loads = BTreeMap::new();
    for (name, sides) in groups {
        let patches = sides
            .iter()
            .cloned()
            .map(LoadPatch::whole)
            .collect::<Vec<_>>();
        let load = integrate_traction(space, &exterior, &patches, order, |x, n| {
            traction_at(x, n, power, coefficient)
        })?;
        loads.insert(name.clone(), load);
    }
    Ok(loads)
}

fn body_load(
    space: &AnalysisSpace,
    power: u32,
    coefficient: f64,
    order: usize,
) -> Result<(Vec<f64>, Vector), &'static str> {
    let rule = gauss(order)?;
    let mut force = vec![0.; 3 * space.ndof];
    if power > 1 {
        for block in &space.blocks {
            for &(u, wu) in &rule {
                for &(v, wv) in &rule {
                    for &(w, ww) in &rule {
                        let s = basis_geometry(block, [u, v, w])?;
                        let density = -C
                            * coefficient
                            * power as f64
                            * (power - 1) as f64
                            * s.position[2].powi(power as i32 - 2);
                        let measure = s.determinant * wu * wv * ww;
                        for (local, &dof) in block.dofs.iter().enumerate() {
                            force[3 * dof + 2] += s.values[local] * density * measure;
                        }
                    }
                }
            }
        }
    }
    let result = [0., 0., force.iter().skip(2).step_by(3).sum()];
    Ok((force, result))
}

fn exact_moment(recipe: SweepRecipe, k: u32) -> f64 {
    let a = recipe.top_multiplier_eighths[0] as f64 / 8.;
    let b = recipe.top_multiplier_eighths[1] as f64 / 8.;
    let d = a - 1.;
    let q2 = d * d + b * b;
    let h = recipe.height_mm as f64;
    SECTION_AREA
        * h.powi(k as i32 + 1)
        * (1. / (k as f64 + 1.) + 2. * d / (k as f64 + 2.) + q2 / (k as f64 + 3.))
}

#[derive(Clone)]
struct Exact {
    volume: f64,
    l2_squared: f64,
    energy_seminorm_squared: f64,
    strain_energy: f64,
    external_work: f64,
    body_force: f64,
    bottom_force: f64,
    top_force: f64,
    lateral_force: f64,
    natural_force: f64,
    body_work: f64,
    natural_work: f64,
}
fn exact(recipe: SweepRecipe, power: u32, coefficient: f64) -> Exact {
    let volume = exact_moment(recipe, 0);
    let l2_squared = coefficient.powi(2) * exact_moment(recipe, 2 * power);
    let energy_seminorm_squared =
        (power as f64 * coefficient).powi(2) * C * exact_moment(recipe, 2 * power - 2);
    let strain_energy = energy_seminorm_squared / 2.;
    let external_work = energy_seminorm_squared;
    let body_force = if power >= 2 {
        -C * coefficient * power as f64 * (power - 1) as f64 * exact_moment(recipe, power - 2)
    } else {
        0.
    };
    let body_work = if power >= 2 {
        -C * coefficient.powi(2)
            * power as f64
            * (power - 1) as f64
            * exact_moment(recipe, 2 * power - 2)
    } else {
        0.
    };
    let bottom_force = if power == 1 {
        -C * coefficient * SECTION_AREA
    } else {
        0.
    };
    let [a, b] = recipe.top_multiplier_eighths.map(|x| x as f64 / 8.);
    let top_area = SECTION_AREA * (a * a + b * b);
    let top_force = C
        * coefficient
        * power as f64
        * (recipe.height_mm as f64).powi(power as i32 - 1)
        * top_area;
    let lateral_force = -body_force - bottom_force - top_force;
    let natural_force = -body_force - bottom_force;
    let natural_work = external_work - body_work;
    Exact {
        volume,
        l2_squared,
        energy_seminorm_squared,
        strain_energy,
        external_work,
        body_force,
        bottom_force,
        top_force,
        lateral_force,
        natural_force,
        body_work,
        natural_work,
    }
}

type EvaluatedField = ([f64; 3], [[f64; 3]; 3], BasisSample);

fn field_at(
    block: &AnalysisBlock,
    displacement: &[f64],
    uvw: Vector,
) -> Result<EvaluatedField, &'static str> {
    let sample = basis_geometry(block, uvw)?;
    let mut value = [0.; 3];
    let mut grad = [[0.; 3]; 3];
    for (i, &dof) in block.dofs.iter().enumerate() {
        for (c, component) in value.iter_mut().enumerate() {
            *component += sample.values[i] * displacement[3 * dof + c];
            for (axis, component_gradient) in grad[c].iter_mut().enumerate() {
                *component_gradient += sample.gradients[i][axis] * displacement[3 * dof + c];
            }
        }
    }
    Ok((value, grad, sample))
}
fn strain(g: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| 0.5 * (g[i][j] + g[j][i])))
}
fn energy_product(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> f64 {
    let ta = a[0][0] + a[1][1] + a[2][2];
    let tb = b[0][0] + b[1][1] + b[2][2];
    LAMBDA * ta * tb
        + 2. * MU
            * (0..3)
                .flat_map(|i| (0..3).map(move |j| a[i][j] * b[i][j]))
                .sum::<f64>()
}

#[derive(Clone, Debug)]
struct Metrics {
    l2: f64,
    energy: f64,
}
fn field_metrics(
    space: &AnalysisSpace,
    disp: &[f64],
    recipe: SweepRecipe,
    power: u32,
    coefficient: f64,
    order: usize,
) -> Result<Metrics, &'static str> {
    let rule = gauss(order)?;
    let mut s = [0.; 2];
    for block in &space.blocks {
        for &(u, wu) in &rule {
            for &(v, wv) in &rule {
                for &(w, ww) in &rule {
                    let (value, grad, g) = field_at(block, disp, [u, v, w])?;
                    let z = g.position[2];
                    let exact_u = coefficient * z.powi(power as i32);
                    let mut exact_g = [[0.; 3]; 3];
                    exact_g[2][2] = coefficient * power as f64 * z.powi(power as i32 - 1);
                    let d = [value[0], value[1], value[2] - exact_u];
                    let de = std::array::from_fn(|i| {
                        std::array::from_fn(|j| grad[i][j] - exact_g[i][j])
                    });
                    let weight = g.determinant * wu * wv * ww;
                    s[0] += dot(d, d) * weight;
                    s[1] += energy_product(strain(de), strain(de)) * weight;
                }
            }
        }
    }
    let ana = exact(recipe, power, coefficient);
    if s.iter().any(|x| !x.is_finite() || *x < 0.) {
        return Err("nonfinite error integral");
    }
    Ok(Metrics {
        l2: (s[0] / ana.l2_squared).sqrt(),
        energy: (s[1] / ana.energy_seminorm_squared).sqrt(),
    })
}

fn exact_natural_work(
    space: &AnalysisSpace,
    groups: &BTreeMap<String, Vec<BlockSide>>,
    power: u32,
    coefficient: f64,
    order: usize,
) -> Result<(f64, BTreeMap<String, f64>), &'static str> {
    let rule = gauss(order)?;
    let mut per_group = BTreeMap::new();
    for (group, sides) in groups {
        if group == "bottom" {
            continue;
        }
        let mut sum = 0.;
        for side in sides {
            let block = space
                .blocks
                .iter()
                .find(|x| x.id == side.block)
                .ok_or("missing exterior side block")?;
            let axes = match side.axis {
                0 => [1, 2],
                1 => [0, 2],
                _ => [0, 1],
            };
            let sign = if side.upper ^ (side.axis == 1) {
                1.
            } else {
                -1.
            };
            for &(u, wu) in &rule {
                for &(v, wv) in &rule {
                    let mut uvw = [0.; 3];
                    uvw[side.axis] = f64::from(u8::from(side.upper));
                    uvw[axes[0]] = u;
                    uvw[axes[1]] = v;
                    let g = basis_geometry(block, uvw)?;
                    let area_vec = cross(g.jacobian_columns[axes[0]], g.jacobian_columns[axes[1]]);
                    let area = norm(area_vec);
                    if !area.is_finite() || area <= 0. {
                        return Err("invalid natural face Jacobian");
                    }
                    let n = area_vec.map(|x| sign * x / area);
                    let f = traction_at(g.position, n, power, coefficient);
                    let uexact = coefficient * g.position[2].powi(power as i32);
                    sum += f[2] * uexact * area * wu * wv;
                }
            }
        }
        per_group.insert(group.clone(), sum);
    }
    Ok((per_group.values().sum(), per_group))
}

fn solve(
    program: &Program,
    level: usize,
    power: u32,
    coefficient: f64,
    order: usize,
) -> Result<Solved, &'static str> {
    let started = Instant::now();
    let space = refine(program, level)?;
    let groups = refined_groups(program, &space)?;
    let stiffness = assemble(&space, YOUNG, POISSON, order)?;
    let all_surfaces = integrate_surface_groups(&space, &groups, power, coefficient, order)?;
    let mut force = vec![0.; 3 * space.ndof];
    let mut natural_resultant = [0.; 3];
    for (group, load) in &all_surfaces {
        if group == "bottom" {
            continue;
        }
        for (f, l) in force.iter_mut().zip(&load.nodal_force_n) {
            *f += l;
        }
        for (d, resultant) in natural_resultant.iter_mut().enumerate() {
            *resultant += load.resultant_n[d];
        }
    }
    let (body, body_resultant) = body_load(&space, power, coefficient, order)?;
    for (f, b) in force.iter_mut().zip(body) {
        *f += b;
    }
    let prescribed = space
        .boundary_dofs
        .get("bottom")
        .ok_or("missing semantic bottom cap")?
        .iter()
        .map(|&d| (d, Ratio::zero()))
        .collect::<Vec<_>>();
    let lift = resolve(space.ndof, &[], &prescribed)?;
    let solved = solve_constrained(&stiffness, [&lift, &lift, &lift], &force, 1e-10, 30_000)?;
    let work = force
        .iter()
        .zip(&solved.displacement_mm)
        .map(|(f, u)| f * u)
        .sum::<f64>();
    let work_identity = (work / (2. * solved.strain_energy_n_mm) - 1.).abs();
    let metrics = field_metrics(
        &space,
        &solved.displacement_mm,
        program.recipe,
        power,
        coefficient,
        order,
    )?;
    let (natural_work, group_work) =
        exact_natural_work(&space, &groups, power, coefficient, order)?;
    let expected = exact(program.recipe, power, coefficient);
    let all_surface_z = all_surfaces.values().map(|x| x.resultant_n[2]).sum::<f64>();
    let net_z = all_surface_z + body_resultant[2];
    let flux_error = (net_z / (natural_resultant[2].abs() + body_resultant[2].abs() + 1.)).abs();
    let bottom_z = all_surfaces
        .get("bottom")
        .ok_or("missing bottom surface load")?
        .resultant_n[2];
    let top_z = all_surfaces
        .get("top")
        .ok_or("missing top surface load")?
        .resultant_n[2];
    let lateral_z = all_surfaces
        .iter()
        .filter(|(name, _)| name.as_str() != "bottom" && name.as_str() != "top")
        .map(|(_, load)| load.resultant_n[2])
        .sum::<f64>();
    let force_errors = [
        (body_resultant[2] - expected.body_force).abs(),
        (bottom_z - expected.bottom_force).abs(),
        (top_z - expected.top_force).abs(),
        (lateral_z - expected.lateral_force).abs(),
        (natural_resultant[2] - expected.natural_force).abs(),
    ];
    let force_scale = [
        expected.body_force,
        expected.bottom_force,
        expected.top_force,
        expected.lateral_force,
        expected.natural_force,
    ]
    .into_iter()
    .map(f64::abs)
    .fold(1.0, f64::max);
    let exact_body_work = body_work_integral(&space, power, coefficient, order)?;
    let exact_work_error = ((natural_work + exact_body_work) / expected.external_work - 1.).abs();
    if work_identity > 1e-7
        || flux_error > 1e-7
        || exact_work_error > 1e-7
        || force_errors.into_iter().any(|e| e > 1e-7 * force_scale)
    {
        return Err("swept manufactured balance or work check failed");
    }
    let mut group_forces = BTreeMap::new();
    let mut group_areas = BTreeMap::new();
    let mut group_moments = BTreeMap::new();
    for (name, load) in all_surfaces {
        group_forces.insert(
            name.clone(),
            Json::Array(load.resultant_n.map(num).to_vec()),
        );
        group_moments.insert(
            name.clone(),
            Json::Array(load.moment_n_mm.map(num).to_vec()),
        );
        group_areas.insert(name, num(load.area_mm2));
    }
    let volume_error = (stiffness.integrated_volume / expected.volume - 1.).abs();
    if order >= 8 && volume_error > 1e-8 {
        return Err("quadrature volume exceeds swept-sweep tolerance");
    }
    Ok(Solved {
        space,
        metrics,
        volume: stiffness.integrated_volume,
        natural_resultant,
        body_resultant,
        group_forces: Json::Object(group_forces),
        group_areas: Json::Object(group_areas),
        group_moments: Json::Object(group_moments),
        work_identity,
        exact_work_error,
        exact_load_work: natural_work,
        exact_body_work,
        exact_group_work: group_work,
        iterations: solved.iterations,
        residual: solved.relative_residual,
        elapsed: started.elapsed().as_secs_f64(),
        rss: peak_rss_kib(),
    })
}

fn body_work_integral(
    space: &AnalysisSpace,
    power: u32,
    coefficient: f64,
    order: usize,
) -> Result<f64, &'static str> {
    let rule = gauss(order)?;
    let mut work = 0.;
    if power == 1 {
        return Ok(0.);
    }
    for block in &space.blocks {
        for &(u, wu) in &rule {
            for &(v, wv) in &rule {
                for &(w, ww) in &rule {
                    let s = basis_geometry(block, [u, v, w])?;
                    let z = s.position[2];
                    work += (-C
                        * coefficient
                        * power as f64
                        * (power - 1) as f64
                        * z.powi(power as i32 - 2))
                        * (coefficient * z.powi(power as i32))
                        * s.determinant
                        * wu
                        * wv
                        * ww;
                }
            }
        }
    }
    Ok(work)
}

#[derive(Clone)]
struct Solved {
    space: AnalysisSpace,
    metrics: Metrics,
    volume: f64,
    natural_resultant: Vector,
    body_resultant: Vector,
    group_forces: Json,
    group_areas: Json,
    group_moments: Json,
    work_identity: f64,
    exact_work_error: f64,
    exact_load_work: f64,
    exact_body_work: f64,
    exact_group_work: BTreeMap<String, f64>,
    iterations: usize,
    residual: f64,
    elapsed: f64,
    rss: usize,
}
fn peak_rss_kib() -> usize {
    fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines().find_map(|l| {
                l.strip_prefix("VmHWM:")
                    .and_then(|v| v.split_whitespace().next()?.parse().ok())
            })
        })
        .unwrap_or(0)
}

fn affine_patch(program: &Program) -> Result<Json, &'static str> {
    let started = Instant::now();
    let power = 1;
    let coefficient = 0.001;
    let q = 8;
    let space = program.base.clone();
    let groups = program.groups.clone();
    let k = assemble(&space, YOUNG, POISSON, q)?;
    let loads = integrate_surface_groups(&space, &groups, power, coefficient, q)?;
    let mut force = vec![0.; 3 * space.ndof];
    let mut natural_resultant = [0.; 3];
    for (g, l) in &loads {
        if g != "bottom" {
            for (f, x) in force.iter_mut().zip(&l.nodal_force_n) {
                *f += x;
            }
            for (d, resultant) in natural_resultant.iter_mut().enumerate() {
                *resultant += l.resultant_n[d];
            }
        }
    }
    let prescribed = space.boundary_dofs["bottom"]
        .iter()
        .map(|&d| (d, Ratio::zero()))
        .collect::<Vec<_>>();
    let lift = resolve(space.ndof, &[], &prescribed)?;
    let sol = solve_constrained(&k, [&lift, &lift, &lift], &force, 1e-11, 10_000)?;
    let mut error: f64 = 0.;
    let mut scale: f64 = 0.;
    for block in &space.blocks {
        for (&dof, h) in block.dofs.iter().zip(&block.controls) {
            let z = h[2] as f64 / h[3] as f64;
            let e = coefficient * z;
            scale = scale.max(e.abs());
            error = error
                .max(sol.displacement_mm[3 * dof].abs())
                .max(sol.displacement_mm[3 * dof + 1].abs())
                .max((sol.displacement_mm[3 * dof + 2] - e).abs());
        }
    }
    let relative = error / scale.max(1e-30);
    let expected = exact(program.recipe, power, coefficient);
    let energy_expected = expected.strain_energy;
    let energy_error = (sol.strain_energy_n_mm / energy_expected - 1.).abs();
    let work = force
        .iter()
        .zip(&sol.displacement_mm)
        .map(|(f, u)| f * u)
        .sum::<f64>();
    let work_error = (work / expected.external_work - 1.).abs();
    if relative > 1e-7 || energy_error > 1e-7 || sol.relative_residual > 1e-10 {
        return Err("affine patch qualification failed");
    }
    Ok(object([
        ("relative_nodal_error", num(relative)),
        ("relative_energy_error", num(energy_error)),
        ("relative_external_work_error", num(work_error)),
        ("relative_residual", num(sol.relative_residual)),
        ("iterations", integer(sol.iterations)),
        ("natural_force_z_n", num(natural_resultant[2])),
        ("wall_s", num(started.elapsed().as_secs_f64())),
    ]))
}

fn record(level: usize, power: u32, c: f64, q: usize, s: &Solved) -> Json {
    object([
        ("level", integer(level)),
        ("power", integer(power as usize)),
        ("coefficient", num(c)),
        ("quadrature", integer(q)),
        ("blocks", integer(s.space.blocks.len())),
        ("scalar_dofs", integer(s.space.ndof)),
        ("l2_relative_error", num(s.metrics.l2)),
        ("energy_relative_error", num(s.metrics.energy)),
        ("volume_q_mm3", num(s.volume)),
        (
            "natural_resultant_n",
            Json::Array(s.natural_resultant.map(num).to_vec()),
        ),
        (
            "body_resultant_n",
            Json::Array(s.body_resultant.map(num).to_vec()),
        ),
        ("surface_group_forces_n", s.group_forces.clone()),
        ("surface_group_areas_mm2", s.group_areas.clone()),
        ("surface_group_moments_n_mm", s.group_moments.clone()),
        ("exact_natural_work_n_mm", num(s.exact_load_work)),
        ("exact_body_work_n_mm", num(s.exact_body_work)),
        (
            "exact_group_work_n_mm",
            Json::Object(
                s.exact_group_work
                    .iter()
                    .map(|(k, v)| (k.clone(), num(*v)))
                    .collect(),
            ),
        ),
        ("exact_load_work_error_relative", num(s.exact_work_error)),
        (
            "discrete_work_identity_error_relative",
            num(s.work_identity),
        ),
        ("residual_relative", num(s.residual)),
        ("iterations", integer(s.iterations)),
        ("wall_s", num(s.elapsed)),
        ("process_hwm_kib", integer(s.rss)),
    ])
}

fn reference_error(
    reference: &Json,
    recipe: SweepRecipe,
    power: u32,
    coefficient: f64,
    value: &Exact,
) -> Result<f64, &'static str> {
    for case in reference.field("cases")?.items()? {
        let rec = case.field("recipe")?.items()?;
        let vals = rec
            .iter()
            .map(|x| x.number::<i32>())
            .collect::<Result<Vec<_>, _>>()?;
        let expected_recipe = [
            recipe.top_multiplier_eighths[0],
            recipe.top_multiplier_eighths[1],
            recipe.height_mm,
            recipe.top_translation_eighths_mm[0],
            recipe.top_translation_eighths_mm[1],
        ];
        if vals.as_slice() == expected_recipe
            && case.field("power")?.number::<u32>()? == power
            && (case
                .field("coefficient")?
                .text()?
                .parse::<f64>()
                .map_err(|_| "bad oracle coefficient")?
                - coefficient)
                .abs()
                < 1e-14
        {
            let mut worst: f64 = 0.;
            for (key, val) in [
                ("volume_mm3", value.volume),
                ("l2_squared_mm5", value.l2_squared),
                ("h1_seminorm_squared_mm3", value.energy_seminorm_squared / C),
                ("strain_energy_n_mm", value.strain_energy),
                ("external_work_n_mm", value.external_work),
                ("body_force_z_n", value.body_force),
                ("bottom_force_z_n", value.bottom_force),
                ("top_force_z_n", value.top_force),
                ("lateral_force_z_n", value.lateral_force),
                ("natural_force_z_n", value.natural_force),
                ("body_work_n_mm", value.body_work),
                ("natural_work_n_mm", value.natural_work),
            ] {
                let refv = case
                    .field(key)?
                    .text()?
                    .parse::<f64>()
                    .map_err(|_| "bad Decimal reference")?;
                let err = if refv.abs() < 1e-10 {
                    (val - refv).abs()
                } else {
                    (val / refv - 1.).abs()
                };
                worst = worst.max(err);
            }
            return Ok(worst);
        }
    }
    Err("missing Decimal reference case")
}

fn main() -> Result<(), &'static str> {
    let started = Instant::now();
    let reference = parse(include_str!("../tests/references/ruled_mechanics.json"))?;
    let recipes = [
        SweepRecipe {
            section: SweepSection::PiercedSquare,
            top_multiplier_eighths: [8, 0],
            top_translation_eighths_mm: [0, 0],
            height_mm: 8,
        },
        SweepRecipe {
            section: SweepSection::PiercedSquare,
            top_multiplier_eighths: [0, 8],
            top_translation_eighths_mm: [4, -8],
            height_mm: 12,
        },
        SweepRecipe {
            section: SweepSection::PiercedSquare,
            top_multiplier_eighths: [-8, 1],
            top_translation_eighths_mm: [-16, 24],
            height_mm: 7,
        },
    ];
    let mut programs = Vec::new();
    let mut oracle_checks = Vec::new();
    let mut affine = Vec::new();
    let mut q8_records = Vec::new();
    let mut q6_records = Vec::new();
    let mut extra_records = Vec::new();
    for recipe in recipes {
        let p = make_program(recipe)?;
        let min_scale = p.certificate.minimum_squared_scale;
        let certificate = object([
            ("ideal_volume_mm3", num(p.certificate.ideal_volume_mm3)),
            (
                "determinant_lower_bound_mm3",
                num(p.certificate.certificate.determinant_lower_bound_mm3),
            ),
            (
                "minimum_squared_scale_ratio",
                Json::Array(min_scale.map(|v| Json::Number(v.to_string())).to_vec()),
            ),
            (
                "groups",
                Json::Array(p.groups.keys().cloned().map(Json::String).collect()),
            ),
        ]);
        affine.push(object([
            (
                "recipe",
                Json::Array(
                    [
                        recipe.top_multiplier_eighths[0],
                        recipe.top_multiplier_eighths[1],
                        recipe.height_mm,
                        recipe.top_translation_eighths_mm[0],
                        recipe.top_translation_eighths_mm[1],
                    ]
                    .map(|v| Json::Number(v.to_string()))
                    .to_vec(),
                ),
            ),
            ("result", affine_patch(&p)?),
        ]));
        programs.push((p, certificate));
    }
    for (program, certificate) in &programs {
        for (power, c) in [(1, 0.001), (3, 0.0001)] {
            let reference_exact = exact(program.recipe, power, c);
            let oracle_err =
                reference_error(&reference, program.recipe, power, c, &reference_exact)?;
            if oracle_err > 1e-12 {
                return Err("Decimal oracle mismatch");
            }
            oracle_checks.push(object([
                ("power", integer(power as usize)),
                ("coefficient", num(c)),
                ("maximum_relative_or_absolute_error", num(oracle_err)),
            ]));
            let mut per_level = Vec::new();
            let mut previous: Option<Metrics> = None;
            let mut q6_by_level = BTreeMap::new();
            let mut q6_case_records = Vec::new();
            for level in 0..=3 {
                eprintln!("S1 q8 recipe {:?} power {power} h{level}", program.recipe);
                let s8 = solve(program, level, power, c, 8)?;
                let r8 = record(level, power, c, 8, &s8);
                let rate_l2 = previous
                    .as_ref()
                    .map_or(Json::Null, |m| num((m.l2 / s8.metrics.l2).log2()));
                let rate_energy = previous
                    .as_ref()
                    .map_or(Json::Null, |m| num((m.energy / s8.metrics.energy).log2()));
                per_level.push(object([
                    ("solution", r8),
                    ("rate_l2", rate_l2),
                    ("rate_energy", rate_energy),
                ]));
                previous = Some(s8.metrics.clone());
                if power == 3 {
                    eprintln!("S1 q6 recipe {:?} power {power} h{level}", program.recipe);
                    let s6 = solve(program, level, power, c, 6)?;
                    q6_by_level.insert(level, (s6.metrics.clone(), s6.volume, s6.exact_work_error));
                    let row = record(level, power, c, 6, &s6);
                    q6_records.push(row.clone());
                    q6_case_records.push(row);
                }
            }
            let final_l2 = per_level[3].field("rate_l2")?.number::<f64>()?;
            let final_en = per_level[3].field("rate_energy")?.number::<f64>()?;
            let target_met = power == 3 && final_l2 >= 2.4 && final_en >= 1.6;
            let mut quadrature_comparisons = Vec::new();
            let mut max_quad_change: f64 = 0.;
            if power == 3 {
                for level in 0..=3 {
                    let q8rec = per_level[level].field("solution")?;
                    let l2_8 = q8rec.field("l2_relative_error")?.number::<f64>()?;
                    let en_8 = q8rec.field("energy_relative_error")?.number::<f64>()?;
                    let (m6, v6, w6) = &q6_by_level[&level];
                    let changes = [
                        (m6.l2 / l2_8 - 1.).abs(),
                        (m6.energy / en_8 - 1.).abs(),
                        (v6 / q8rec.field("volume_q_mm3")?.number::<f64>()? - 1.).abs(),
                        (*w6 - q8rec
                            .field("exact_load_work_error_relative")?
                            .number::<f64>()?)
                        .abs(),
                    ];
                    let worst = changes.into_iter().fold(0., f64::max);
                    max_quad_change = max_quad_change.max(worst);
                    quadrature_comparisons.push(object([
                        ("level", integer(level)),
                        ("relative_l2_error_change", num((m6.l2 / l2_8 - 1.).abs())),
                        (
                            "relative_energy_error_change",
                            num((m6.energy / en_8 - 1.).abs()),
                        ),
                        (
                            "relative_volume_change",
                            num((v6 / q8rec.field("volume_q_mm3")?.number::<f64>()? - 1.).abs()),
                        ),
                        (
                            "exact_work_error_change",
                            num((*w6
                                - q8rec
                                    .field("exact_load_work_error_relative")?
                                    .number::<f64>()?)
                            .abs()),
                        ),
                    ]));
                }
            }
            let additional = if power == 3 && max_quad_change > 1e-7 {
                let mut rows = Vec::new();
                for q in [12, 16] {
                    for level in 0..=3 {
                        eprintln!(
                            "S1 extra q{q} recipe {:?} power {power} h{level}",
                            program.recipe
                        );
                        let sx = solve(program, level, power, c, q)?;
                        rows.push(record(level, power, c, q, &sx));
                    }
                }
                rows
            } else {
                Vec::new()
            };
            extra_records.extend(additional.clone());
            let sweep_key = format!(
                "{:?}_{}_{}_{}_{}_{}",
                program.recipe.top_multiplier_eighths[0],
                program.recipe.top_multiplier_eighths[1],
                program.recipe.height_mm,
                program.recipe.top_translation_eighths_mm[0],
                program.recipe.top_translation_eighths_mm[1],
                power
            );
            let case = object([
                (
                    "recipe",
                    Json::Array(vec![
                        Json::Number(program.recipe.top_multiplier_eighths[0].to_string()),
                        Json::Number(program.recipe.top_multiplier_eighths[1].to_string()),
                        integer(program.recipe.height_mm as usize),
                        Json::Number(program.recipe.top_translation_eighths_mm[0].to_string()),
                        Json::Number(program.recipe.top_translation_eighths_mm[1].to_string()),
                    ]),
                ),
                ("power", integer(power as usize)),
                ("coefficient", num(c)),
                ("certificate", certificate.clone()),
                ("q8_h0_to_h3", Json::Array(per_level)),
                ("q6_h0_to_h3", Json::Array(q6_case_records)),
                ("quadrature_comparison", Json::Array(quadrature_comparisons)),
                ("max_q6_q8_relative_change", num(max_quad_change)),
                ("extra_q12_q16", Json::Array(additional)),
                (
                    "final_rate_threshold_met",
                    if power == 3 {
                        Json::Bool(target_met)
                    } else {
                        Json::Null
                    },
                ),
                (
                    "quasi_pinch_threshold_exceedance",
                    Json::Bool(
                        power == 3
                            && program.recipe.top_multiplier_eighths == [-8, 1]
                            && !target_met,
                    ),
                ),
                ("reference_key", Json::String(sweep_key)),
            ]);
            q8_records.push(case);
        }
    }
    let mut all_times = Vec::new();
    let mut max_rss = 0usize;
    for case in &q8_records {
        for level in case.field("q8_h0_to_h3")?.items()? {
            let solution = level.field("solution")?;
            all_times.push(solution.field("wall_s")?.number::<f64>()?);
            max_rss = max_rss.max(solution.field("process_hwm_kib")?.number::<usize>()?);
        }
    }
    for record in &q6_records {
        all_times.push(record.field("wall_s")?.number::<f64>()?);
        max_rss = max_rss.max(record.field("process_hwm_kib")?.number::<usize>()?);
    }
    let output=object([("schema",Json::String("ruled-manufactured-elasticity-s1-v1".into())),("material",object([("young_mpa",num(YOUNG)),("poisson",num(POISSON))])),("manufactured_fields",Json::Array(vec![Json::String("affine: u=(0,0,0.001 z) mm, zero body force; exact traction on all non-bottom faces".into()),Json::String("cubic: u=(0,0,0.0001 z^3) mm, f=-div(sigma), exact traction on all non-bottom faces".into())])),("affine_patch_q8",Json::Array(affine)),("cases",Json::Array(q8_records)),("q6_all_levels",Json::Array(q6_records)),("extra_q12_q16",Json::Array(extra_records)),("decimal_oracle_checks",Json::Array(oracle_checks)),("resource_summary",object([("q8_solve_wall_s_min",num(all_times.iter().copied().fold(f64::INFINITY,f64::min))), ("q8_solve_wall_s_max",num(all_times.iter().copied().fold(0.,f64::max))), ("process_hwm_kib_max",integer(max_rss)), ("campaign_wall_s",num(started.elapsed().as_secs_f64())),("rss_scope",Json::String("process high-water mark sampled after each solve; cumulative within this run".into()))])),("limitations",Json::String("The ruled-sweep audit establishes the fixed geometric family and injectivity; these results measure only this manufactured loading and the listed quadrature/mesh levels. The quasi-pinch case is reported even if rate or quadrature thresholds fail.".into()))]);
    fs::write(
        "docs/results/ruled-mechanics-20261125.json",
        format!("{output}\n"),
    )
    .map_err(|_| "cannot write ruled mechanics report")?;
    println!("{output}");
    Ok(())
}
