//! Plane-strain restriction of the three-dimensional Q2 space on a square plate.
use solidcraft_trivariate::analysis_space::{bore_analysis_space, AnalysisSpace};
use solidcraft_trivariate::constrained_elasticity::solve_constrained;
use solidcraft_trivariate::constraint_system::{resolve, Ratio};
use solidcraft_trivariate::elasticity::{assemble, basis_geometry, gauss};
use solidcraft_trivariate::face_constraints::Row;
use solidcraft_trivariate::h_refine::refine_cross_section;
use solidcraft_trivariate::movable_bore::build_movable_bore;
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

fn quadrant(space: AnalysisSpace) -> Result<AnalysisSpace, &'static str> {
    let original_count = space.blocks.len();
    let mut blocks = Vec::new();
    for block in space.blocks {
        if block.controls.iter().all(|p| p[0] >= 0 && p[1] >= 0) {
            blocks.push(block);
        } else if !(block.controls.iter().all(|p| p[0] <= 0)
            || block.controls.iter().all(|p| p[1] <= 0))
        {
            return Err("a cell straddles a symmetry plane");
        }
    }
    if 4 * blocks.len() != original_count {
        return Err("incorrect symmetry partition");
    }
    let used: BTreeSet<_> = blocks.iter().flat_map(|b| b.dofs.iter().copied()).collect();
    let compact: BTreeMap<_, _> = used
        .into_iter()
        .enumerate()
        .map(|(new, old)| (old, new))
        .collect();
    for b in &mut blocks {
        for id in &mut b.dofs {
            *id = compact[id];
        }
    }
    let boundary_dofs = space
        .boundary_dofs
        .into_iter()
        .filter_map(|(tag, nodes)| {
            let retained: Vec<_> = nodes
                .into_iter()
                .filter_map(|i| compact.get(&i).copied())
                .collect();
            (!retained.is_empty()).then_some((tag, retained))
        })
        .collect();
    Ok(AnalysisSpace {
        ndof: compact.len(),
        blocks,
        boundary_dofs,
    })
}

fn run(level: usize, order: usize, radius: f64, quarter: bool) -> Result<(), &'static str> {
    let start = Instant::now();
    let ticks = 2. * radius;
    if !(4.0..=10.0).contains(&ticks) || ticks.fract() != 0. {
        return Err("invalid half-millimeter bore radius");
    }
    let radius_ticks = ticks as i32;
    let candidate =
        build_movable_bore(0, 0, radius_ticks).map_err(|_| "bore construction failed")?;
    let blocks: Vec<_> = candidate.blocks.into_iter().map(|b| b.block).collect();
    let mut space = bore_analysis_space([0, 0], radius_ticks, &blocks)?;
    if quarter && level == 0 {
        return Err("quarter symmetry requires the first angular split");
    }
    for step in 0..level {
        space = refine_cross_section(&space)?;
        if quarter && step == 0 {
            space = quadrant(space)?;
        }
    }
    let geometry = start.elapsed().as_secs_f64();
    run_space(&space, level, order, radius, geometry, start, quarter)
}

fn run_space(
    space: &AnalysisSpace,
    level: usize,
    order: usize,
    radius: f64,
    geometry: f64,
    start: Instant,
    quarter: bool,
) -> Result<(), &'static str> {
    let mut columns = BTreeMap::new();
    let mut symmetry = [BTreeSet::new(), BTreeSet::new()];
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
            if quarter {
                for axis in 0..2 {
                    if p[axis] == 0 {
                        symmetry[axis].insert(ids[0]);
                    }
                }
            }
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
    let expected_area = if quarter { 160. } else { 640. };
    if (area / expected_area - 1.).abs() > 1e-10
        || (!quarter && gauges.len() != 3)
        || probes.is_empty()
    {
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
    if !quarter {
        symmetry[0].extend([gauges[&(0, -10)], gauges[&(0, 10)]]);
        symmetry[1].insert(gauges[&(10, 0)]);
    }
    let fixed = symmetry
        .each_ref()
        .map(|nodes| nodes.iter().map(|&i| (i, zero)).collect::<Vec<_>>());
    let lx = resolve(space.ndof, &eq, &fixed[0])?;
    let ly = resolve(space.ndof, &eq, &fixed[1])?;
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
    let mut root_of: Vec<_> = (0..space.ndof).collect();
    for &(id, ref row) in &eq {
        root_of[id] = row.terms[0].0;
    }
    for (dir, prescribed) in symmetry.iter().enumerate() {
        let mut by_root = BTreeMap::<usize, f64>::new();
        for (i, &root) in root_of.iter().enumerate() {
            if prescribed.contains(&root) {
                *by_root.entry(root).or_default() += internal[3 * i + dir] - force[3 * i + dir];
            }
        }
        for value in by_root.values() {
            max_reaction = max_reaction.max(value.abs());
            reactions[dir] += value;
        }
    }
    let mut resultant = reactions;
    for v in force.as_chunks::<3>().0.iter() {
        resultant[0] += v[0];
        resultant[1] += v[1];
    }
    let balance = resultant[0].hypot(resultant[1]) / force_scale;
    if balance > 1e-6 || (!quarter && max_reaction / force_scale > 1e-6) {
        return Err("non-negligible artificial gauge reaction");
    }
    let stresses = probes
        .iter()
        .map(|(b, p)| cauchy_stress(b, &solution.displacement_mm, *p, E, NU).map(|s| s[0][0]))
        .collect::<Result<Vec<_>, _>>()?;
    let min = stresses.iter().copied().fold(f64::INFINITY, f64::min);
    let max = stresses.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    println!("{{\"quarter_symmetry\":{quarter},\"radius_mm\":{radius},\"level\":{level},\"quadrature\":{order},\"blocks\":{},\"full_scalar_dofs\":{},\"reduced_dofs\":{},\"hoop_min_mpa\":{min},\"hoop_max_mpa\":{max},\"force_balance_relative\":{balance},\"residual\":{},\"iterations\":{},\"geometry_refinement_s\":{geometry},\"constraints_load_s\":{},\"assembly_s\":{},\"solve_s\":{},\"total_s\":{}}}",space.blocks.len(),space.ndof,solution.reduced_values.len(),solution.relative_residual,solution.iterations,boundary_done-geometry,assembly_done-boundary_done,solve_done-assembly_done,start.elapsed().as_secs_f64());
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
    let quarter = match args.next().as_deref() {
        Some("quarter") => true,
        Some("full") | None => false,
        _ => return Err("symmetry must be full or quarter"),
    };
    run(level, order, radius, quarter)
}
