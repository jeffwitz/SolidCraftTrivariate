//! Three-dimensional pressure loading of a thick cylinder with restrained end faces.
use solidcraft_trivariate::projected_bulk::{assemble_projected_bulk, ProjectedElasticity};
use solidcraft_trivariate::{
    analysis_space::{section_analysis_space, AnalysisSpace},
    boundary_traction::{integrate_traction, LoadPatch},
    constrained_elasticity::solve_constrained,
    constraint_system::{resolve, Ratio},
    elasticity::{assemble, basis_geometry, gauss, BasisSample, Stiffness},
    h_refine::refine_cross_section,
    incidence::BlockSide,
    json::{object, Json},
    sections::section_blocks,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    time::Instant,
};

const YOUNG: f64 = 210000.;
fn n(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}
type Tensor = [[f64; 3]; 3];

#[derive(Clone, Copy)]
struct Pressure {
    inner: f64,
    outer: f64,
    nu: f64,
}
impl Pressure {
    fn a(self) -> f64 {
        (4. * self.inner - 16. * self.outer) / 12.
    }
    fn b(self) -> f64 {
        64. * (self.inner - self.outer) / 12.
    }
    fn radial(self, r: f64) -> f64 {
        (1. + self.nu) / YOUNG * ((1. - 2. * self.nu) * self.a() * r + self.b() / r)
    }
    fn energy(self) -> f64 {
        std::f64::consts::PI
            * 8.
            * (self.inner * 2. * self.radial(2.) - self.outer * 4. * self.radial(4.))
    }
    fn gradient(self, x: [f64; 3]) -> Tensor {
        let r = x[0].hypot(x[1]);
        let e = [x[0] / r, x[1] / r, 0.];
        let rr = (1. + self.nu) / YOUNG * ((1. - 2. * self.nu) * self.a() - self.b() / r.powi(2));
        let tt = self.radial(r) / r;
        std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                if i < 2 && j < 2 {
                    (if i == j { tt } else { 0. }) + (rr - tt) * e[i] * e[j]
                } else {
                    0.
                }
            })
        })
    }
    fn stress(self, g: Tensor) -> Tensor {
        let mu = YOUNG / (2. * (1. + self.nu));
        let lam = YOUNG * self.nu / ((1. + self.nu) * (1. - 2. * self.nu));
        let trace = g[0][0] + g[1][1] + g[2][2];
        std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                mu * (g[i][j] + g[j][i]) + if i == j { lam * trace } else { 0. }
            })
        })
    }
    fn norm_squared(self, g: Tensor) -> f64 {
        self.stress(g)
            .iter()
            .flatten()
            .zip(g.iter().flatten())
            .map(|(a, b)| a * b)
            .sum()
    }
}
fn field(ids: &[usize], s: &BasisSample, u: &[f64]) -> ([f64; 3], Tensor) {
    let mut value = [0.; 3];
    let mut gradient = [[0.; 3]; 3];
    for (i, &id) in ids.iter().enumerate() {
        for c in 0..3 {
            value[c] += s.values[i] * u[3 * id + c];
            for (out, derivative) in gradient[c].iter_mut().zip(s.gradients[i]) {
                *out += derivative * u[3 * id + c];
            }
        }
    }
    (value, gradient)
}

fn field_errors(
    space: &AnalysisSpace,
    p: Pressure,
    u: &[f64],
    order: usize,
) -> Result<(f64, f64, f64), &'static str> {
    let xy = gauss(order)?;
    let z = gauss(3)?;
    let mut sums = [0.; 5];
    for block in &space.blocks {
        for &(a, wa) in &xy {
            for &(b, wb) in &xy {
                for &(c, wc) in &z {
                    let s = basis_geometry(block, [a, b, c])?;
                    let x = s.position;
                    let r = x[0].hypot(x[1]);
                    let ur = p.radial(r);
                    let exact = [ur * x[0] / r, ur * x[1] / r, 0.];
                    let eg = p.gradient(x);
                    let (value, g) = field(&block.dofs, &s, u);
                    let de: [f64; 3] = std::array::from_fn(|i| value[i] - exact[i]);
                    let dg = std::array::from_fn(|i| std::array::from_fn(|j| g[i][j] - eg[i][j]));
                    let dv = s.determinant * wa * wb * wc;
                    sums[0] += dv * dot(&de, &de);
                    sums[1] += dv * dot(&exact, &exact);
                    sums[2] += dv * p.norm_squared(dg);
                    sums[3] += dv * p.norm_squared(eg);
                    sums[4] += 0.5 * dv * p.norm_squared(g);
                }
            }
        }
    }
    Ok((
        (sums[0] / sums[1]).sqrt(),
        (sums[2] / sums[3]).sqrt(),
        sums[4],
    ))
}

fn wall_patches(
    space: &AnalysisSpace,
    tag: &str,
    upper: bool,
) -> Result<Vec<LoadPatch>, &'static str> {
    let group: BTreeSet<_> = space
        .boundary_dofs
        .get(tag)
        .ok_or("missing wall")?
        .iter()
        .copied()
        .collect();
    let j = usize::from(upper) * 2;
    Ok(space
        .blocks
        .iter()
        .filter(|b| (0..3).all(|i| (0..3).all(|k| group.contains(&b.dofs[9 * i + 3 * j + k]))))
        .map(|b| {
            LoadPatch::whole(BlockSide {
                block: b.id.clone(),
                axis: 1,
                upper,
            })
        })
        .collect())
}
fn wall_metrics(
    space: &AnalysisSpace,
    patches: &[LoadPatch],
    p: Pressure,
    u: &[f64],
    radius: f64,
    projected: Option<&ProjectedElasticity>,
) -> Result<Json, &'static str> {
    let chosen: BTreeSet<_> = patches.iter().map(|p| p.side.block.as_str()).collect();
    let q = gauss(12)?;
    let axial = gauss(3)?;
    let mut area = 0.;
    let mut mean = 0.;
    let mut hoop_error: f64 = 0.;
    let mut primal_hoop_error: f64 = 0.;
    let exact = p.a() + p.b() / radius.powi(2);
    for b in space
        .blocks
        .iter()
        .filter(|b| chosen.contains(b.id.as_str()))
    {
        let pressure = projected
            .map(|m| m.mean_stress_coefficients(&b.id, u))
            .transpose()?;
        for &(t, wt) in &q {
            for &(w, ww) in &axial {
                let s = basis_geometry(b, [t, f64::from(radius == 2.), w])?;
                let [x, y, _] = s.position;
                let r = x.hypot(y);
                let hoop = [-y / r, x / r, 0.];
                let (value, g) = field(&b.dofs, &s, u);
                let sig = p.stress(g);
                let sh: f64 = (0..3).map(|i| hoop[i] * dot(&sig[i], &hoop)).sum();
                primal_hoop_error = primal_hoop_error.max((sh - exact).abs());
                let sh = if let Some(coefficients) = pressure {
                    let mean = ProjectedElasticity::mean_stress_at(
                        coefficients,
                        [t, f64::from(radius == 2.), w],
                    )?;
                    sh + mean - (sig[0][0] + sig[1][1] + sig[2][2]) / 3.
                } else {
                    sh
                };
                let du = s.jacobian_columns[0];
                let dw = s.jacobian_columns[2];
                let cross: [f64; 3] = std::array::from_fn(|i| {
                    du[(i + 1) % 3] * dw[(i + 2) % 3] - du[(i + 2) % 3] * dw[(i + 1) % 3]
                });
                let ds = dot(&cross, &cross).sqrt() * wt * ww;
                area += ds;
                mean += ds * (value[0] * x + value[1] * y) / r;
                hoop_error = hoop_error.max((sh - exact).abs());
            }
        }
    }
    Ok(object([
        ("area_mm2", n(area)),
        ("mean_radial_displacement_mm", n(mean / area)),
        ("exact_radial_displacement_mm", n(p.radial(radius))),
        ("exact_hoop_mpa", n(exact)),
        ("max_sampled_hoop_error_mpa", n(hoop_error)),
        ("primal_max_sampled_hoop_error_mpa", n(primal_hoop_error)),
        (
            "max_sampled_hoop_relative",
            n(hoop_error / exact.abs().max(1.)),
        ),
    ]))
}
fn run(
    level: usize,
    order: usize,
    p: Pressure,
    alternative: bool,
    projected: bool,
    capture_field: bool,
) -> Result<Json, &'static str> {
    if level > 4
        || ![6, 8, 12].contains(&order)
        || !(0.0..0.5).contains(&p.nu)
        || p.inner < 0.
        || p.outer < 0.
        || p.energy() <= 0.
    {
        return Err("invalid pressure cylinder case");
    }
    let start = Instant::now();
    let family = "full_annulus_v1";
    let mut space = section_analysis_space(family, &section_blocks(family)?)?;
    for _ in 0..level {
        space = refine_cross_section(&space)?;
    }
    let walls = [
        wall_patches(&space, "circle_r2", true)?,
        wall_patches(&space, "circle_r4", false)?,
    ];
    let exterior = walls.iter().flatten().map(|p| p.side.clone()).collect();
    let mut force = vec![0.; 3 * space.ndof];
    let mut area_error: f64 = 0.;
    for (patches, (pressure, radius)) in walls.iter().zip([(p.inner, 2.), (p.outer, 4.)]) {
        let surface = integrate_traction(&space, &exterior, patches, order, |_, normal| {
            normal.map(|v| -pressure * v)
        })?;
        area_error =
            area_error.max((surface.area_mm2 / (16. * std::f64::consts::PI * radius) - 1.).abs());
        for (f, v) in force.iter_mut().zip(surface.nodal_force_n) {
            *f += v;
        }
    }
    let mut axes: [Vec<(f64, usize)>; 2] = [vec![], vec![]];
    for b in &space.blocks {
        for (h, &id) in b.controls.iter().zip(&b.dofs) {
            if h[2] != 0 {
                continue;
            }
            for c in 0..2 {
                if h[c] == 0 {
                    axes[c].push((h[1 - c] as f64 / h[3] as f64, id));
                }
            }
        }
    }
    let mut conditions: [Vec<_>; 3] = [vec![], vec![], vec![]];
    for c in 0..2 {
        axes[c].sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let first = axes[c].first().ok_or("missing axial symmetry gauge")?.1;
        let last = axes[c].last().ok_or("missing axial symmetry gauge")?.1;
        conditions[c].push((first, Ratio::zero()));
        if c == usize::from(alternative) {
            conditions[c].push((last, Ratio::zero()));
        }
    }
    let zfixed: BTreeSet<_> = ["z0", "z8"]
        .into_iter()
        .flat_map(|tag| space.boundary_dofs[tag].iter().copied())
        .collect();
    conditions[2] = zfixed.into_iter().map(|id| (id, Ratio::zero())).collect();
    let [lx, ly, lz] = conditions.map(|values| resolve(space.ndof, &[], &values));
    let (lx, ly, lz) = (lx?, ly?, lz?);
    let prep = start.elapsed().as_secs_f64();
    let clock = Instant::now();
    let mixed = if projected {
        Some(assemble_projected_bulk(&space, YOUNG, p.nu, order)?)
    } else {
        None
    };
    let full;
    let k = if let Some(m) = &mixed {
        &m.stiffness
    } else {
        full = assemble(&space, YOUNG, p.nu, order)?;
        &full
    };
    let assembly = clock.elapsed().as_secs_f64();
    let clock = Instant::now();
    let solved = solve_constrained(k, [&lx, &ly, &lz], &force, 1e-10, 20000)?;
    let solve = clock.elapsed().as_secs_f64();
    let u = &solved.displacement_mm;
    let residual = k.equilibrium_residual(u, &force)?;
    let reactions = ["z0", "z8"].map(|name| {
        space.boundary_dofs[name]
            .iter()
            .map(|id| -residual[3 * id + 2])
            .sum::<f64>()
    });
    let exact_reaction = p.nu * 2. * p.a() * 12. * std::f64::consts::PI;
    let wall = [
        wall_metrics(&space, &walls[0], p, u, 2., mixed.as_ref())?,
        wall_metrics(&space, &walls[1], p, u, 4., mixed.as_ref())?,
    ];
    let errors = field_errors(&space, p, u, 12)?;
    let low = field_errors(&space, p, u, 8)?;
    let applied: [f64; 3] =
        std::array::from_fn(|a| force.as_chunks::<3>().0.iter().map(|x| x[a]).sum());
    let mut report = object([
        ("level", n(level)),
        ("order", n(order)),
        ("nu", n(p.nu)),
        ("inner_mpa", n(p.inner)),
        ("outer_mpa", n(p.outer)),
        ("alternative_gauge", Json::Bool(alternative)),
        ("projected_bulk", Json::Bool(projected)),
        ("physical_strain_energy_n_mm", n(errors.2)),
        ("vector_dofs", n(k.size())),
        ("blocks", n(space.blocks.len())),
        ("load_area_relative_error", n(area_error)),
        ("pressure_resultant_n", Json::Array(applied.map(n).into())),
        ("axial_reactions_n", Json::Array(reactions.map(n).into())),
        ("exact_top_reaction_n", n(exact_reaction)),
        ("relative_l2", n(errors.0)),
        ("relative_energy_norm", n(errors.1)),
        ("check_q8_relative_l2", n(low.0)),
        ("check_q8_relative_energy_norm", n(low.1)),
        ("pcg_iterations", n(solved.iterations)),
        ("relative_residual", n(solved.relative_residual)),
        ("work_n_mm", n(dot(&force, u))),
        ("strain_energy_n_mm", n(solved.strain_energy_n_mm)),
        ("exact_energy_n_mm", n(p.energy())),
        ("inner", wall[0].clone()),
        ("outer", wall[1].clone()),
        ("preparation_s", n(prep)),
        ("assembly_s", n(assembly)),
        ("solve_s", n(solve)),
        ("total_s", n(start.elapsed().as_secs_f64())),
    ]);
    if capture_field {
        let Json::Object(fields) = &mut report else {
            unreachable!()
        };
        fields.insert(
            "displacement_mm_xyz".into(),
            Json::Array(u.iter().copied().map(n).collect()),
        );
    }
    Ok(report)
}

struct PlaneSystem {
    space: AnalysisSpace,
    force: Vec<f64>,
    walls: [Vec<LoadPatch>; 2],
    fixed: [BTreeSet<usize>; 3],
    matrix: Stiffness,
    projected: Option<ProjectedElasticity>,
    preparation_s: f64,
    assembly_s: f64,
}

fn plane_system(
    level: usize,
    order: usize,
    p: Pressure,
    alternative: bool,
    projected: bool,
) -> Result<PlaneSystem, &'static str> {
    if level > 6
        || ![6, 8, 12].contains(&order)
        || !(0.0..0.5).contains(&p.nu)
        || p.inner < 0.
        || p.outer < 0.
        || p.energy() <= 0.
    {
        return Err("invalid pressure cylinder case");
    }
    let start = Instant::now();
    let family = "full_annulus_v1";
    let mut space = section_analysis_space(family, &section_blocks(family)?)?;
    for _ in 0..level {
        space = refine_cross_section(&space)?;
    }
    let walls = [
        wall_patches(&space, "circle_r2", true)?,
        wall_patches(&space, "circle_r4", false)?,
    ];
    let exterior = walls
        .iter()
        .flatten()
        .map(|patch| patch.side.clone())
        .collect();
    let mut force = vec![0.; 3 * space.ndof];
    for (patches, pressure) in walls.iter().zip([p.inner, p.outer]) {
        let load = integrate_traction(&space, &exterior, patches, order, |_, normal| {
            normal.map(|x| -pressure * x)
        })?;
        for (target, value) in force.iter_mut().zip(load.nodal_force_n) {
            *target += value;
        }
    }
    let mut axes: [Vec<(f64, usize)>; 2] = [vec![], vec![]];
    for block in &space.blocks {
        for (h, &id) in block.controls.iter().zip(&block.dofs) {
            if h[2] != 0 {
                continue;
            }
            for c in 0..2 {
                if h[c] == 0 {
                    axes[c].push((h[1 - c] as f64 / h[3] as f64, id));
                }
            }
        }
    }
    let mut prescribed: [Vec<(usize, Ratio)>; 3] = [vec![], vec![], vec![]];
    for c in 0..2 {
        axes[c].sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        prescribed[c].push((
            axes[c].first().ok_or("missing axial symmetry gauge")?.1,
            Ratio::zero(),
        ));
        if c == usize::from(alternative) {
            prescribed[c].push((
                axes[c].last().ok_or("missing axial symmetry gauge")?.1,
                Ratio::zero(),
            ));
        }
    }
    let axial_fixed: BTreeSet<_> = ["z0", "z8"]
        .into_iter()
        .flat_map(|tag| space.boundary_dofs[tag].iter().copied())
        .collect();
    prescribed[2] = axial_fixed
        .into_iter()
        .map(|id| (id, Ratio::zero()))
        .collect();
    let [lx, ly, lz] = prescribed.map(|values| resolve(space.ndof, &[], &values));
    let lifts = [lx?, ly?, lz?];
    let fixed = std::array::from_fn(|c| {
        lifts[c]
            .lines
            .iter()
            .filter(|line| line.terms.is_empty())
            .map(|line| line.dof)
            .collect()
    });
    let preparation_s = start.elapsed().as_secs_f64();
    let clock = Instant::now();
    let projected_matrix = if projected {
        Some(assemble_projected_bulk(&space, YOUNG, p.nu, order)?)
    } else {
        None
    };
    let matrix = if let Some(m) = &projected_matrix {
        m.stiffness.clone_for_example()?
    } else {
        assemble(&space, YOUNG, p.nu, order)?
    };
    let assembly_s = clock.elapsed().as_secs_f64();
    Ok(PlaneSystem {
        space,
        force,
        walls,
        fixed,
        matrix,
        projected: projected_matrix,
        preparation_s,
        assembly_s,
    })
}

trait CloneStiffnessForExample {
    fn clone_for_example(&self) -> Result<Stiffness, &'static str>;
}
impl CloneStiffnessForExample for Stiffness {
    fn clone_for_example(&self) -> Result<Stiffness, &'static str> {
        Ok(Stiffness {
            offsets: self.offsets.clone(),
            columns: self.columns.clone(),
            values: self.values.clone(),
            integrated_volume: self.integrated_volume,
        })
    }
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}
fn xy_key(h: [i128; 4]) -> Result<[i128; 3], &'static str> {
    if h[3] <= 0 {
        return Err("invalid plane-strain homogeneous weight");
    }
    let divisor = gcd_u128(
        gcd_u128(h[0].unsigned_abs(), h[1].unsigned_abs()),
        h[3] as u128,
    );
    let d = i128::try_from(divisor).map_err(|_| "plane-strain gauge overflow")?;
    Ok([h[0] / d, h[1] / d, h[3] / d])
}

fn plane_map(space: &AnalysisSpace) -> Result<(Vec<usize>, Vec<[i128; 3]>), &'static str> {
    let mut dof_keys = vec![None; space.ndof];
    let mut unique = BTreeSet::new();
    for block in &space.blocks {
        for (&id, &control) in block.dofs.iter().zip(&block.controls) {
            let key = xy_key(control)?;
            if dof_keys[id].is_some_and(|prior| prior != key) {
                return Err("one 3D scalar DOF has incompatible section coordinates");
            }
            dof_keys[id] = Some(key);
            unique.insert(key);
        }
    }
    if dof_keys.iter().any(Option::is_none) {
        return Err("incomplete 3D-to-plane DOF map");
    }
    let coordinates = unique.into_iter().collect::<Vec<_>>();
    let ids: BTreeMap<_, _> = coordinates
        .iter()
        .copied()
        .enumerate()
        .map(|(i, key)| (key, i))
        .collect();
    let mapping = dof_keys
        .into_iter()
        .map(|key| ids[&key.expect("checked map")])
        .collect();
    Ok((mapping, coordinates))
}

struct ReducedPlane {
    mapping: Vec<usize>,
    coordinates: Vec<[i128; 3]>,
    row_ptr: Vec<usize>,
    columns: Vec<usize>,
    values: Vec<f64>,
    rhs: Vec<f64>,
    fixed: Vec<usize>,
    free: Vec<usize>,
    reduction_s: f64,
}

fn reduce_to_plane(system: &PlaneSystem) -> Result<ReducedPlane, &'static str> {
    let started = Instant::now();
    let (mapping, coordinates) = plane_map(&system.space)?;
    let n = coordinates.len() * 2;
    let mut rows = vec![BTreeMap::<usize, f64>::new(); n];
    for i in 0..system.space.ndof {
        for pos in system.matrix.offsets[i]..system.matrix.offsets[i + 1] {
            let j = system.matrix.columns[pos];
            for a in 0..2 {
                for b in 0..2 {
                    *rows[2 * mapping[i] + a]
                        .entry(2 * mapping[j] + b)
                        .or_default() += system.matrix.values[pos][a][b];
                }
            }
        }
    }
    let mut row_ptr = vec![0];
    let mut columns = vec![];
    let mut values = vec![];
    for row in rows {
        for (col, value) in row {
            if value != 0. {
                columns.push(col);
                values.push(value);
            }
        }
        row_ptr.push(columns.len());
    }
    let mut rhs = vec![0.; n];
    for i in 0..system.space.ndof {
        for c in 0..2 {
            rhs[2 * mapping[i] + c] += system.force[3 * i + c];
        }
    }
    let mut fixed = BTreeSet::new();
    for c in 0..2 {
        for &dof in &system.fixed[c] {
            fixed.insert(2 * mapping[dof] + c);
        }
    }
    let free = (0..n).filter(|id| !fixed.contains(id)).collect::<Vec<_>>();
    let reduction_s = started.elapsed().as_secs_f64();
    Ok(ReducedPlane {
        mapping,
        coordinates,
        row_ptr,
        columns,
        values,
        rhs,
        fixed: fixed.into_iter().collect(),
        free,
        reduction_s,
    })
}

fn json_numbers<T: ToString>(v: impl IntoIterator<Item = T>) -> Json {
    Json::Array(v.into_iter().map(n).collect())
}
fn write_json(path: &str, value: &Json) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::write(path, format!("{value}\n"))?;
    Ok(())
}

fn plane_export(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    if args.len() != 8 {
        return Err(
            "usage: pressure_cylinder plane-export LEVEL ORDER NU PI PO full|projected OUTPUT.json"
                .into(),
        );
    }
    let level = args[1].parse()?;
    let order = args[2].parse()?;
    let p = Pressure {
        nu: args[3].parse()?,
        inner: args[4].parse()?,
        outer: args[5].parse()?,
    };
    let projected = match args[6].as_str() {
        "full" => false,
        "projected" => true,
        _ => return Err("expected full or projected formulation".into()),
    };
    let system = plane_system(level, order, p, false, projected)?;
    let ReducedPlane {
        mapping,
        coordinates,
        row_ptr,
        columns,
        values,
        rhs,
        fixed,
        free,
        reduction_s,
    } = reduce_to_plane(&system)?;
    let mut raw = Vec::new();
    for key in &coordinates {
        raw.push(Json::Array(key.iter().copied().map(n).collect()));
    }
    let matrix = object([
        ("dimension", n(rhs.len())),
        ("row_ptr", json_numbers(row_ptr)),
        ("column_indices", json_numbers(columns)),
        ("values", json_numbers(values)),
    ]);
    let output = object([
        (
            "schema",
            Json::String("plane_strain_sparse_system_v1".into()),
        ),
        ("level", n(level)),
        ("order", n(order)),
        ("nu", n(p.nu)),
        ("inner_mpa", n(p.inner)),
        ("outer_mpa", n(p.outer)),
        (
            "formulation",
            Json::String(
                if projected {
                    "projected_bulk_q1"
                } else {
                    "displacement"
                }
                .into(),
            ),
        ),
        ("space_dofs_3d", n(system.space.ndof)),
        ("scalar_dofs_2d", n(coordinates.len())),
        ("vector_dofs_2d", n(rhs.len())),
        ("plane_coordinates_homogeneous", Json::Array(raw)),
        ("map_3d_scalar_dof_to_2d", json_numbers(mapping)),
        ("matrix_csr", matrix),
        ("rhs_xy", json_numbers(rhs)),
        ("fixed_vector_dofs", json_numbers(fixed)),
        ("free_vector_dofs", json_numbers(free)),
        ("preparation_s", n(system.preparation_s)),
        ("assembly_s", n(system.assembly_s)),
        ("reduction_s", n(reduction_s)),
    ]);
    write_json(&args[7], &output)?;
    Ok(())
}

fn plane_read(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    if args.len() != 9 {
        return Err("usage: pressure_cylinder plane-read LEVEL ORDER NU PI PO full|projected SOLUTION.json OUTPUT.json".into());
    }
    let level = args[1].parse()?;
    let order = args[2].parse()?;
    let p = Pressure {
        nu: args[3].parse()?,
        inner: args[4].parse()?,
        outer: args[5].parse()?,
    };
    let projected = match args[6].as_str() {
        "full" => false,
        "projected" => true,
        _ => return Err("expected full or projected formulation".into()),
    };
    let solution = solidcraft_trivariate::json::parse(&std::fs::read_to_string(&args[7])?)?;
    if solution.field("schema")?.text()? != "plane_strain_solution_v1"
        || solution.field("level")?.number::<usize>()? != level
        || solution.field("order")?.number::<usize>()? != order
        || solution.field("nu")?.number::<f64>()? != p.nu
        || solution.field("inner_mpa")?.number::<f64>()? != p.inner
        || solution.field("outer_mpa")?.number::<f64>()? != p.outer
    {
        return Err("plane-strain solution metadata mismatch".into());
    }
    let system = plane_system(level, order, p, false, projected)?;
    let ReducedPlane {
        mapping,
        coordinates,
        free,
        reduction_s,
        ..
    } = reduce_to_plane(&system)?;
    let xy = solution
        .field("displacement_xy")?
        .items()?
        .iter()
        .map(Json::number::<f64>)
        .collect::<Result<Vec<_>, _>>()?;
    if xy.len() != 2 * coordinates.len() {
        return Err("plane-strain solution coefficient count mismatch".into());
    }
    let mut u = vec![0.; 3 * system.space.ndof];
    for (dof, &node) in mapping.iter().enumerate() {
        u[3 * dof] = xy[2 * node];
        u[3 * dof + 1] = xy[2 * node + 1];
    }
    let residual = system.matrix.equilibrium_residual(&u, &system.force)?;
    let residual_norm = (0..system.space.ndof)
        .flat_map(|node| (0..2).map(move |c| (node, c)))
        .filter(|(node, c)| !system.fixed[*c].contains(node))
        .map(|(node, c)| residual[3 * node + c].powi(2))
        .sum::<f64>()
        .sqrt();
    let force_norm = system
        .force
        .iter()
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt()
        .max(1.);
    let projected_residual = residual_norm / force_norm;
    let applied: [f64; 3] =
        std::array::from_fn(|c| system.force.as_chunks::<3>().0.iter().map(|x| x[c]).sum());
    let work = dot(&system.force, &u);
    let energy = 0.5 * dot(&u, &system.matrix.apply(&u)?);
    let wall = [
        wall_metrics(
            &system.space,
            &system.walls[0],
            p,
            &u,
            2.,
            system.projected.as_ref(),
        )?,
        wall_metrics(
            &system.space,
            &system.walls[1],
            p,
            &u,
            4.,
            system.projected.as_ref(),
        )?,
    ];
    let errors = field_errors(&system.space, p, &u, 12)?;
    let mut report = object([
        (
            "schema",
            Json::String("plane_strain_reconstruction_v1".into()),
        ),
        ("level", n(level)),
        ("order", n(order)),
        ("nu", n(p.nu)),
        ("inner_mpa", n(p.inner)),
        ("outer_mpa", n(p.outer)),
        ("projected_bulk", Json::Bool(projected)),
        ("relative_l2", n(errors.0)),
        ("relative_energy_norm", n(errors.1)),
        ("strain_energy_n_mm", n(errors.2)),
        ("physical_strain_energy_n_mm", n(errors.2)),
        ("work_n_mm", n(work)),
        ("strain_energy_matrix_n_mm", n(energy)),
        (
            "energy_work_relative",
            n((2. * energy - work).abs() / work.abs().max(1.)),
        ),
        ("pressure_resultant_n", json_numbers(applied)),
        ("relative_residual_3d_free", n(projected_residual)),
        ("inner", wall[0].clone()),
        ("outer", wall[1].clone()),
        ("space_dofs_3d", n(system.space.ndof)),
        ("scalar_dofs_2d", n(coordinates.len())),
        ("assembly_s", n(system.assembly_s)),
        ("reduction_s", n(reduction_s)),
        ("free_vector_dofs_2d", n(free.len())),
        ("volume_mm3", n(system.matrix.integrated_volume)),
    ]);
    if let Json::Object(fields) = &mut report {
        fields.insert("displacement_mm_xyz".into(), json_numbers(u));
    }
    write_json(&args[8], &report)?;
    Ok(())
}

fn full_export(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    if args.len() != 8 {
        return Err(
            "usage: pressure_cylinder full-export LEVEL ORDER NU PI PO full|projected OUTPUT.json"
                .into(),
        );
    }
    let level = args[1].parse()?;
    let order = args[2].parse()?;
    let p = Pressure {
        nu: args[3].parse()?,
        inner: args[4].parse()?,
        outer: args[5].parse()?,
    };
    let projected = match args[6].as_str() {
        "full" => false,
        "projected" => true,
        _ => return Err("expected full or projected formulation".into()),
    };
    let system = plane_system(level, order, p, false, projected)?;
    let started = Instant::now();
    let mut fixed = Vec::new();
    for c in 0..3 {
        for &dof in &system.fixed[c] {
            fixed.push(3 * dof + c);
        }
    }
    fixed.sort_unstable();
    let row_ptr = &system.matrix.offsets;
    let mut values = Vec::with_capacity(system.matrix.values.len() * 9);
    for block in &system.matrix.values {
        values.push(Json::Array(
            block
                .iter()
                .map(|row| json_numbers(row.iter().copied()))
                .collect(),
        ));
    }
    let matrix = object([
        ("block_size", n(3)),
        ("block_rows", n(system.space.ndof)),
        ("row_ptr", json_numbers(row_ptr.iter().copied())),
        (
            "column_indices",
            json_numbers(system.matrix.columns.iter().copied()),
        ),
        ("values_3x3", Json::Array(values)),
    ]);
    let export_s = started.elapsed().as_secs_f64();
    let output = object([
        (
            "schema",
            Json::String("full_cylinder_sparse_system_v1".into()),
        ),
        ("level", n(level)),
        ("order", n(order)),
        ("nu", n(p.nu)),
        ("inner_mpa", n(p.inner)),
        ("outer_mpa", n(p.outer)),
        (
            "formulation",
            Json::String(
                if projected {
                    "projected_bulk_q1"
                } else {
                    "displacement"
                }
                .into(),
            ),
        ),
        ("space_dofs_3d", n(system.space.ndof)),
        ("vector_dofs_3d", n(3 * system.space.ndof)),
        ("matrix_block_csr", matrix),
        ("rhs_xyz", json_numbers(system.force.iter().copied())),
        ("fixed_vector_dofs", json_numbers(fixed)),
        ("preparation_s", n(system.preparation_s)),
        ("assembly_s", n(system.assembly_s)),
        ("export_s", n(export_s)),
        ("integrated_volume_mm3", n(system.matrix.integrated_volume)),
    ]);
    write_json(&args[7], &output)?;
    Ok(())
}

fn full_read(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    if args.len() != 9 {
        return Err("usage: pressure_cylinder full-read LEVEL ORDER NU PI PO full|projected SOLUTION.json OUTPUT.json".into());
    }
    let level = args[1].parse()?;
    let order = args[2].parse()?;
    let p = Pressure {
        nu: args[3].parse()?,
        inner: args[4].parse()?,
        outer: args[5].parse()?,
    };
    let projected = match args[6].as_str() {
        "full" => false,
        "projected" => true,
        _ => return Err("expected full or projected formulation".into()),
    };
    let solution = solidcraft_trivariate::json::parse(&std::fs::read_to_string(&args[7])?)?;
    if solution.field("schema")?.text()? != "full_cylinder_solution_v1"
        || solution.field("level")?.number::<usize>()? != level
        || solution.field("order")?.number::<usize>()? != order
        || solution.field("nu")?.number::<f64>()? != p.nu
        || solution.field("inner_mpa")?.number::<f64>()? != p.inner
        || solution.field("outer_mpa")?.number::<f64>()? != p.outer
        || solution.field("formulation")?.text()?
            != if projected {
                "projected_bulk_q1"
            } else {
                "displacement"
            }
    {
        return Err("full-cylinder solution metadata mismatch".into());
    }
    let system = plane_system(level, order, p, false, projected)?;
    let u = solution
        .field("displacement_xyz")?
        .items()?
        .iter()
        .map(Json::number::<f64>)
        .collect::<Result<Vec<_>, _>>()?;
    if u.len() != 3 * system.space.ndof {
        return Err("full-cylinder solution coefficient count mismatch".into());
    }
    let residual = system.matrix.equilibrium_residual(&u, &system.force)?;
    let free_residual_squared = (0..system.space.ndof)
        .flat_map(|node| (0..3).map(move |c| (node, c)))
        .filter(|(node, c)| !system.fixed[*c].contains(node))
        .map(|(node, c)| residual[3 * node + c].powi(2))
        .sum::<f64>();
    let force_norm = system
        .force
        .iter()
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt()
        .max(1.);
    let free_residual = free_residual_squared.sqrt() / force_norm;
    let work = dot(&system.force, &u);
    let energy = 0.5 * dot(&u, &system.matrix.apply(&u)?);
    let applied: [f64; 3] =
        std::array::from_fn(|c| system.force.as_chunks::<3>().0.iter().map(|x| x[c]).sum());
    let wall = [
        wall_metrics(
            &system.space,
            &system.walls[0],
            p,
            &u,
            2.,
            system.projected.as_ref(),
        )?,
        wall_metrics(
            &system.space,
            &system.walls[1],
            p,
            &u,
            4.,
            system.projected.as_ref(),
        )?,
    ];
    let errors = field_errors(&system.space, p, &u, 12)?;
    let report = object([
        (
            "schema",
            Json::String("full_cylinder_reconstruction_v1".into()),
        ),
        ("level", n(level)),
        ("order", n(order)),
        ("nu", n(p.nu)),
        ("inner_mpa", n(p.inner)),
        ("outer_mpa", n(p.outer)),
        ("projected_bulk", Json::Bool(projected)),
        ("relative_l2", n(errors.0)),
        ("relative_energy_norm", n(errors.1)),
        ("physical_strain_energy_n_mm", n(errors.2)),
        ("work_n_mm", n(work)),
        ("strain_energy_matrix_n_mm", n(energy)),
        (
            "energy_work_relative",
            n((2. * energy - work).abs() / work.abs().max(1.)),
        ),
        ("pressure_resultant_n", json_numbers(applied)),
        ("relative_residual_3d_free", n(free_residual)),
        ("inner", wall[0].clone()),
        ("outer", wall[1].clone()),
        ("space_dofs_3d", n(system.space.ndof)),
        ("assembly_s", n(system.assembly_s)),
        ("volume_mm3", n(system.matrix.integrated_volume)),
    ]);
    write_json(&args[8], &report)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a = std::env::args().skip(1).collect::<Vec<_>>();
    if a.first().is_some_and(|arg| arg == "full-export") {
        return full_export(&a);
    }
    if a.first().is_some_and(|arg| arg == "full-read") {
        return full_read(&a);
    }
    if a.first().is_some_and(|arg| arg == "plane-export") {
        return plane_export(&a);
    }
    if a.first().is_some_and(|arg| arg == "plane-read") {
        return plane_read(&a);
    }
    if a.first().is_some_and(|arg| arg == "field") {
        if a.len() != 9 {
            return Err("usage: pressure_cylinder field LEVEL ORDER NU PI PO ALTERNATE full|projected OUTPUT.json".into());
        }
        let projected = match a[7].as_str() {
            "full" => false,
            "projected" => true,
            _ => return Err("expected full or projected formulation".into()),
        };
        let report = run(
            a[1].parse()?,
            a[2].parse()?,
            Pressure {
                nu: a[3].parse()?,
                inner: a[4].parse()?,
                outer: a[5].parse()?,
            },
            a[6].parse()?,
            projected,
            true,
        )?;
        write_json(&a[8], &report)?;
        println!("{{\"status\":\"field_exported\",\"path\":{:?}}}", a[8]);
        return Ok(());
    }
    if !(a.len() == 6 || (a.len() == 7 && a[6] == "--projected")) {
        return Err(
            "usage: pressure_cylinder level order nu inner_MPa outer_MPa alternate_gauge_boolean [--projected]"
                .into(),
        );
    }
    println!(
        "{}",
        run(
            a[0].parse()?,
            a[1].parse()?,
            Pressure {
                nu: a[2].parse()?,
                inner: a[3].parse()?,
                outer: a[4].parse()?
            },
            a[5].parse()?,
            a.len() == 7,
            false,
        )?
    );
    Ok(())
}
