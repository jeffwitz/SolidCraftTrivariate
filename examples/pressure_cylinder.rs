//! Three-dimensional pressure loading of a thick cylinder with restrained end faces.
use solidcraft_trivariate::projected_bulk::{assemble_projected_bulk, ProjectedElasticity};
use solidcraft_trivariate::{
    analysis_space::{section_analysis_space, AnalysisSpace},
    boundary_traction::{integrate_traction, LoadPatch},
    constrained_elasticity::solve_constrained,
    constraint_system::{resolve, Ratio},
    elasticity::{assemble, basis_geometry, gauss, BasisSample},
    h_refine::refine_cross_section,
    incidence::BlockSide,
    json::{object, Json},
    sections::section_blocks,
};
use std::{collections::BTreeSet, time::Instant};

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
    Ok(object([
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
    ]))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a = std::env::args().skip(1).collect::<Vec<_>>();
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
            a.len() == 7
        )?
    );
    Ok(())
}
