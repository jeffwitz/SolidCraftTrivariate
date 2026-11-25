//! Plane-strain restriction of the three-dimensional Q2 space on a square plate.
use solidcraft_trivariate::analysis_space::{section_analysis_space, AnalysisSpace};
use solidcraft_trivariate::constrained_elasticity::solve_constrained;
use solidcraft_trivariate::constraint_system::{resolve, Ratio};
use solidcraft_trivariate::elasticity::{assemble, basis_geometry, gauss};
use solidcraft_trivariate::face_constraints::Row;
use solidcraft_trivariate::h_refine::refine_cross_section;
use solidcraft_trivariate::sections::section_blocks;
use solidcraft_trivariate::stress::cauchy_stress;
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

const E: f64 = 210_000.;
const NU: f64 = 0.3;
const FAR: f64 = 100.;

/// Normal and shear components in the polar frame and their action on the
/// specified outer normal. Pressure on the circular inner boundary is zero.
fn traction(x: f64, y: f64, radius: f64, normal: [f64; 2]) -> [f64; 2] {
    let r = x.hypot(y);
    let (cx, sy) = (x / r, y / r);
    let d = (radius / r).powi(2);
    let cos2 = cx * cx - sy * sy;
    let sin2 = 2. * cx * sy;
    let rr = FAR / 2. * (1. - d + (1. - 4. * d + 3. * d * d) * cos2);
    let tt = FAR / 2. * (1. + d - (1. + 3. * d * d) * cos2);
    let rt = -FAR / 2. * (1. + 2. * d - 3. * d * d) * sin2;
    let nr = normal[0] * cx + normal[1] * sy;
    let nt = -normal[0] * sy + normal[1] * cx;
    let (fr, ft) = (rr * nr + rt * nt, rt * nr + tt * nt);
    [cx * fr - sy * ft, sy * fr + cx * ft]
}

fn run(level: usize, order: usize, radius: f64) -> Result<(), &'static str> {
    let start = Instant::now();
    if radius != 4. {
        return Err("fixed pierced-square template has radius 4 mm");
    }
    let family = "pierced_square_v1";
    let blocks = section_blocks(family)?;
    let mut space = section_analysis_space(family, &blocks)?;
    for _ in 0..level {
        space = refine_cross_section(&space)?;
    }
    let geometry = start.elapsed().as_secs_f64();
    run_space(&space, level, order, radius, geometry, start)
}

fn run_space(
    space: &AnalysisSpace,
    level: usize,
    order: usize,
    radius: f64,
    geometry: f64,
    start: Instant,
) -> Result<(), &'static str> {
    let mut columns = BTreeMap::new();
    let mut gauges: BTreeMap<(i128, i128), usize> = BTreeMap::new();
    let mut force = vec![0.; 3 * space.ndof];
    let mut force_scale = 0.;
    let mut area = 0.;
    let quad = gauss(order)?;
    let mut probes = vec![];
    for b in &space.blocks {
        for (controls, ids) in b
            .controls
            .as_chunks::<3>()
            .0
            .iter()
            .zip(b.dofs.as_chunks::<3>().0.iter())
        {
            for &id in &ids[1..] {
                if columns.insert(id, ids[0]).is_some_and(|old| old != ids[0]) {
                    return Err("inconsistent extrusion column");
                }
            }
            let p = controls[0];
            for xy in [(0, -10), (0, 10), (10, 0)] {
                if p[0] == xy.0 * p[3] && p[1] == xy.1 * p[3] {
                    gauges.insert(xy, ids[0]);
                }
            }
        }
        // The full outer square is at v=0 in the four original charts.
        let side: Vec<_> = (0..3)
            .flat_map(|i| (0..3).map(move |k| b.controls[9 * i + k]))
            .collect();
        for (axis, value) in [(0, -10), (0, 10), (1, -10), (1, 10)] {
            if !side.iter().all(|p| p[axis] == value * p[3]) {
                continue;
            }
            let mut normal = [0.; 2];
            normal[axis] = value as f64 / 10.;
            for &(u, wu) in &quad {
                for &(w, ww) in &quad {
                    let s = basis_geometry(b, [u, 0., w])?;
                    let [du, _, dw] = s.jacobian_columns;
                    let measure = du[0].hypot(du[1]) * dw[2] * wu * ww;
                    let load = traction(s.position[0], s.position[1], radius, normal);
                    area += measure;
                    force_scale += load[0].hypot(load[1]) * measure;
                    for (idx, &dof) in b.dofs.iter().enumerate() {
                        for c in 0..2 {
                            force[3 * dof + c] += load[c] * measure * s.values[idx];
                        }
                    }
                }
            }
        }
        for u in [0., 0.5, 1.] {
            let uvw = [u, 1., 0.5];
            let pos = basis_geometry(b, uvw)?.position;
            if pos[0].abs() < 1e-11 && (pos[1] - radius).abs() < 1e-11 {
                probes.push((b, uvw));
            }
        }
    }
    if (area / 640. - 1.).abs() > 1e-10 || gauges.len() != 3 || probes.is_empty() {
        return Err("invalid exterior, gauges or probe");
    }
    let eq: Vec<_> = columns
        .into_iter()
        .map(|(dependent, master)| {
            (
                dependent,
                Row {
                    terms: vec![(master, 1)],
                    denominator: 1,
                },
            )
        })
        .collect();
    let zero = Ratio::zero();
    let lx = resolve(
        space.ndof,
        &eq,
        &[(gauges[&(0, -10)], zero), (gauges[&(0, 10)], zero)],
    )?;
    let ly = resolve(space.ndof, &eq, &[(gauges[&(10, 0)], zero)])?;
    let roots: BTreeSet<_> = space
        .blocks
        .iter()
        .flat_map(|b| b.dofs.as_chunks::<3>().0.iter().map(|c| c[0]))
        .collect();
    let zfixed: Vec<_> = roots.iter().map(|&dof| (dof, zero)).collect();
    let lz = resolve(space.ndof, &eq, &zfixed)?;
    let boundary_done = start.elapsed().as_secs_f64();
    let k = assemble(space, E, NU, order)?;
    let assembly_done = start.elapsed().as_secs_f64();
    let solution = solve_constrained(&k, [&lx, &ly, &lz], &force, 1e-10, 20000)?;
    let solve_done = start.elapsed().as_secs_f64();
    let internal = k.apply(&solution.displacement_mm)?;
    // Evaluate gauge reactions after summing all members of each z column.
    let mut reactions = [0.; 2];
    let mut max_reaction: f64 = 0.;
    for (dir, idx) in [
        (0, gauges[&(0, -10)]),
        (0, gauges[&(0, 10)]),
        (1, gauges[&(10, 0)]),
    ] {
        let mut column = BTreeSet::new();
        for b in &space.blocks {
            for ids in b.dofs.as_chunks::<3>().0.iter() {
                if ids[0] == idx {
                    column.extend(ids.iter().copied());
                }
            }
        }
        let reaction: f64 = column
            .iter()
            .map(|i| internal[3 * i + dir] - force[3 * i + dir])
            .sum();
        max_reaction = max_reaction.max(reaction.abs());
        reactions[dir] += reaction;
    }
    let mut resultant = reactions;
    for v in force.as_chunks::<3>().0.iter() {
        resultant[0] += v[0];
        resultant[1] += v[1];
    }
    let balance = resultant[0].hypot(resultant[1]) / force_scale;
    if balance > 1e-6 || max_reaction / force_scale > 1e-6 {
        return Err("non-negligible artificial gauge reaction");
    }
    let stresses = probes
        .iter()
        .map(|(b, p)| cauchy_stress(b, &solution.displacement_mm, *p, E, NU).map(|s| s[0][0]))
        .collect::<Result<Vec<_>, _>>()?;
    let min = stresses.iter().copied().fold(f64::INFINITY, f64::min);
    let max = stresses.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    println!("{{\"radius_mm\":{radius},\"level\":{level},\"quadrature\":{order},\"blocks\":{},\"full_scalar_dofs\":{},\"reduced_dofs\":{},\"hoop_min_mpa\":{min},\"hoop_max_mpa\":{max},\"force_balance_relative\":{balance},\"residual\":{},\"iterations\":{},\"geometry_refinement_s\":{geometry},\"constraints_load_s\":{},\"assembly_s\":{},\"solve_s\":{},\"total_s\":{}}}",space.blocks.len(),space.ndof,solution.reduced_values.len(),solution.relative_residual,solution.iterations,boundary_done-geometry,assembly_done-boundary_done,solve_done-assembly_done,start.elapsed().as_secs_f64());
    Ok(())
}

fn main() -> Result<(), &'static str> {
    let mut args = std::env::args().skip(1);
    let level = args
        .next()
        .ok_or("missing level")?
        .parse()
        .map_err(|_| "invalid level")?;
    let order = args
        .next()
        .ok_or("missing order")?
        .parse()
        .map_err(|_| "invalid order")?;
    let radius = args
        .next()
        .ok_or("missing radius")?
        .parse()
        .map_err(|_| "invalid radius")?;
    run(level, order, radius)
}
