//! Affine patch and loaded-bracket measurements on locally constrained Q2 spaces.
use solidcraft_trivariate::{
    adaptive_partition::{AdaptivePartition, LeafKey, LeafSide},
    adaptive_space::{build_adaptive_space, AdaptiveSpace},
    analysis_space::{parameterized_bracket_space, AnalysisBlock},
    boundary_traction::{integrate_traction, LoadPatch},
    bracket_parameters::build_bracket_parameters,
    bracket_parameters_audit::audit_bracket_parameters,
    bracket_solve::{bracket_load, bracket_tip_displacement, solve_bracket},
    constrained_elasticity::solve_constrained,
    constraint_system::{resolve, Ratio},
    elasticity::{assemble, basis_geometry},
    h_refine::refine_uniform,
    incidence::BlockSide,
    json::{object, Json},
};
use std::{
    collections::{BTreeMap, BTreeSet},
    time::Instant,
};

const A: [[i128; 3]; 3] = [[2, 1, -3], [-2, 4, 1], [1, 2, -1]];
const SHIFT: [i128; 3] = [300, -200, 70];
fn n(x: f64) -> Json {
    Json::Number(x.to_string())
}
fn value(r: Ratio) -> f64 {
    r.numerator() as f64 / r.denominator() as f64
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(u, v)| u * v).sum()
}

fn side_key(adaptive: &AdaptiveSpace, side: &LeafSide) -> BlockSide {
    BlockSide {
        block: adaptive.leaf_blocks[&side.leaf].clone(),
        axis: side.root_side.axis,
        upper: side.root_side.upper,
    }
}

type TraceProbe = ([f64; 3], [[f64; 3]; 2], f64);

fn trace_probe(
    b: &AnalysisBlock,
    side: &BlockSide,
    uv: [f64; 2],
    dofs: &[f64],
) -> Result<TraceProbe, &'static str> {
    let axes = match side.axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    };
    let mut p = [0.; 3];
    p[side.axis] = f64::from(u8::from(side.upper));
    p[axes[0]] = uv[0];
    p[axes[1]] = uv[1];
    let result = basis_geometry(b, p)?;
    let scalar = result
        .values
        .iter()
        .zip(&b.dofs)
        .map(|(v, &i)| v * dofs[i])
        .sum();
    Ok((
        result.position,
        axes.map(|i| result.jacobian_columns[i]),
        scalar,
    ))
}

fn verify_traces(
    a: &AdaptiveSpace,
    partition: &AdaptivePartition,
) -> Result<(usize, f64, f64), &'static str> {
    let mut state = 84291373_u64;
    let master = a
        .lift
        .free_roots
        .iter()
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((state >> 11) as f64) / (1_u64 << 53) as f64
        })
        .collect::<Vec<_>>();
    let dofs = a
        .lift
        .lines
        .iter()
        .map(|l| l.terms.iter().map(|&(i, r)| master[i] * value(r)).sum())
        .collect::<Vec<_>>();
    let by_id: BTreeMap<_, _> = a.space.blocks.iter().map(|b| (&b.id, b)).collect();
    let (mut count, mut position, mut jump) = (0, 0_f64, 0_f64);
    for patch in partition.interfaces()? {
        let (s, t) = if patch.first.leaf.level >= patch.second.leaf.level {
            (&patch.first, &patch.second)
        } else {
            (&patch.second, &patch.first)
        };
        let s = side_key(a, s);
        let t = side_key(a, t);
        for uv in [[0.1, 0.2], [0.37, 0.83], [0.9, 0.3], [0.7, 0.6]] {
            let (x, _, f) = trace_probe(by_id[&s.block], &s, uv, &dofs)?;
            let mut v = [0.5, 0.5];
            for _ in 0..32 {
                let (y, [p, q], _) = trace_probe(by_id[&t.block], &t, v, &dofs)?;
                let e = std::array::from_fn(|i| x[i] - y[i]);
                if dot(e, e) < 1e-26 {
                    break;
                }
                let (pp, qq, pq) = (dot(p, p), dot(q, q), dot(p, q));
                let det = pp * qq - pq * pq;
                if det <= 0. {
                    return Err("singular interface inversion in measurement");
                }
                let (ep, eq) = (dot(p, e), dot(q, e));
                v[0] = (v[0] + (qq * ep - pq * eq) / det).clamp(0., 1.);
                v[1] = (v[1] + (pp * eq - pq * ep) / det).clamp(0., 1.);
            }
            let (y, _, g) = trace_probe(by_id[&t.block], &t, v, &dofs)?;
            position = position.max(
                dot(
                    std::array::from_fn(|i| x[i] - y[i]),
                    std::array::from_fn(|i| x[i] - y[i]),
                )
                .sqrt(),
            );
            jump = jump.max((f - g).abs());
            count += 1;
        }
    }
    if position > 1e-10 || jump > 1e-10 {
        return Err("independent physical trace inversion or field continuity failed");
    }
    Ok((count, position, jump))
}

fn run(rx: u8, scenario: &str) -> Result<Json, &'static str> {
    let t = Instant::now();
    let geometry = build_bracket_parameters(rx, 64)?;
    let certified = audit_bracket_parameters(rx, 64, &geometry.blocks)?;
    let parent = parameterized_bracket_space(rx, 64, &geometry.blocks)?;
    let ids = parent
        .blocks
        .iter()
        .map(|b| b.id.clone())
        .collect::<Vec<_>>();
    let mut p = AdaptivePartition::new(&ids, &certified.certificate.incidence, 4, 1500)?;
    let focus = LeafKey::root("hole_x_sector_east_z10_30");
    let depth2 = LeafKey {
        root: focus.root.clone(),
        level: 1,
        index: [1, 0, 0],
    };
    match scenario {
        "h0" => (),
        "h1" => {
            p.refine(&ids.iter().map(LeafKey::root).collect::<Vec<_>>())?;
        }
        "local1" => {
            p.refine(&[focus])?;
        }
        "local2" => {
            p.refine(&[focus])?;
            p.refine(&[depth2])?;
        }
        _ => return Err("unknown refinement scenario"),
    }
    let adaptive = build_adaptive_space(&parent, &p)?;
    let space = &adaptive.space;
    for (col, &i) in adaptive.lift.free_roots.iter().enumerate() {
        if adaptive.lift.lines[i].terms != [(col, Ratio::one())] {
            return Err("missing identity minor in lift");
        }
    }
    let (probe_count, geom_err, trace_jump) = verify_traces(&adaptive, &p)?;
    let mut affine = vec![[Ratio::zero(); 3]; space.ndof];
    let mut coords = vec![[0.; 3]; space.ndof];
    for b in &space.blocks {
        for (&id, h) in b.dofs.iter().zip(&b.controls) {
            coords[id] = std::array::from_fn(|j| h[j] as f64 / h[3] as f64);
            for d in 0..3 {
                affine[id][d] = Ratio::new(
                    SHIFT[d] * h[3] + (0..3).map(|j| A[d][j] * h[j]).sum::<i128>(),
                    1000 * h[3],
                )?;
            }
        }
    }
    let mut affine_lift = 0_f64;
    for l in &adaptive.lift.lines {
        for (c, &coefficient) in affine[l.dof].iter().enumerate() {
            let expected = value(coefficient);
            let actual = l
                .terms
                .iter()
                .map(|&(i, r)| value(r) * value(affine[adaptive.lift.free_roots[i]][c]))
                .sum::<f64>();
            affine_lift = affine_lift.max((actual - expected).abs());
        }
    }
    if affine_lift > 1e-10 {
        return Err("affine field cannot be represented");
    }
    let startup = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let k = assemble(space, 210000., 0.3, 8)?;
    let assembly = t.elapsed().as_secs_f64();
    let diag = k
        .values
        .iter()
        .flatten()
        .flatten()
        .map(|x| x.abs())
        .fold(0_f64, f64::max);
    let mut rigid_max = 0_f64;
    for mode in 0..6 {
        let mut v = vec![0.; 3 * space.ndof];
        for (i, x) in coords.iter().enumerate() {
            if mode < 3 {
                v[3 * i + mode] = 1.;
            } else {
                let d = mode - 3;
                v[3 * i + (d + 1) % 3] = -x[(d + 2) % 3];
                v[3 * i + (d + 2) % 3] = x[(d + 1) % 3];
            }
        }
        let kv = k.apply(&v)?;
        let residual = kv.iter().map(|s| s * s).sum::<f64>().sqrt()
            / (diag * v.iter().map(|s| s * s).sum::<f64>().sqrt());
        rigid_max = rigid_max.max(residual);
    }
    if rigid_max > 1e-10 {
        return Err("discrete rigid motion has nonzero stiffness");
    }
    let exterior: BTreeSet<_> = p
        .exterior_patches()
        .iter()
        .map(|x| side_key(&adaptive, &x.side))
        .collect();
    let tags: BTreeMap<_, _> = certified.certificate.design_faces.iter().cloned().collect();
    let natural = p
        .exterior_patches()
        .into_iter()
        .filter(|x| tags[&x.side.root_side] != "y50")
        .map(|x| LoadPatch::whole(side_key(&adaptive, &x.side)))
        .collect::<Vec<_>>();
    let eps: [[f64; 3]; 3] =
        std::array::from_fn(|i| std::array::from_fn(|j| (A[i][j] + A[j][i]) as f64 / 2000.));
    let tr = (0..3).map(|d| eps[d][d]).sum::<f64>();
    let sigma: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            210000. / 1.3 * eps[i][j]
                + if i == j {
                    210000. * 0.3 / (1.3 * 0.4) * tr
                } else {
                    0.
                }
        })
    });
    let load = integrate_traction(space, &exterior, &natural, 8, |_, normal| {
        sigma.map(|r| dot(r, normal))
    })?;
    let bc = |d: usize| {
        resolve(
            space.ndof,
            &adaptive.constraints,
            &space.boundary_dofs["y50"]
                .iter()
                .map(|&i| (i, affine[i][d]))
                .collect::<Vec<_>>(),
        )
    };
    let (lx, ly, lz) = (bc(0)?, bc(1)?, bc(2)?);
    let t = Instant::now();
    let patch = solve_constrained(&k, [&lx, &ly, &lz], &load.nodal_force_n, 1e-11, 20000)?;
    let exact_norm = affine
        .iter()
        .flatten()
        .map(|r| value(*r).abs())
        .fold(0_f64, f64::max);
    let patch_err = affine
        .iter()
        .flatten()
        .zip(&patch.displacement_mm)
        .map(|(r, u)| (value(*r) - u).abs())
        .fold(0_f64, f64::max)
        / exact_norm;
    let ideal_energy =
        0.5 * certified.volume_mm3 * (0..3).map(|i| dot(eps[i], sigma[i])).sum::<f64>();
    let energy_err = (patch.strain_energy_n_mm / ideal_energy - 1.).abs();
    if energy_err > 1e-7 || patch_err > 1e-7 {
        return Err("physical affine patch failed");
    }
    let patch_s = t.elapsed().as_secs_f64();
    let fixed = space.boundary_dofs["y50"]
        .iter()
        .map(|&i| (i, Ratio::zero()))
        .collect::<Vec<_>>();
    let mounted = resolve(space.ndof, &adaptive.constraints, &fixed)?;
    let (f, _) = bracket_load(space, 8)?;
    let t = Instant::now();
    let loaded = solve_constrained(&k, [&mounted; 3], &f, 1e-10, 20000)?;
    let solve_s = t.elapsed().as_secs_f64();
    let tip = bracket_tip_displacement(space, &loaded.displacement_mm)?[1];
    let work: f64 = f
        .iter()
        .zip(&loaded.displacement_mm)
        .map(|(a, b)| a * b)
        .sum();
    let work_error = (work / (2. * loaded.strain_energy_n_mm) - 1.).abs();
    if work_error > 1e-8 {
        return Err("loaded bracket work discrepancy");
    }
    let reference = match scenario {
        "h0" => Some(solve_bracket(&parent, 8)?),
        "h1" => Some(solve_bracket(&refine_uniform(&parent)?, 8)?),
        _ => None,
    };
    let mut comparison = Json::Null;
    if let Some(r) = reference {
        let err = (tip / r.point_y_displacement_mm - 1.)
            .abs()
            .max((loaded.strain_energy_n_mm / r.strain_energy_n_mm - 1.).abs());
        if err > 1e-7 {
            return Err("uniform mechanical space differs from conforming reference");
        }
        comparison = n(err);
    }
    Ok(object([
        ("radius_x_mm", n(f64::from(rx) / 16.)),
        ("scenario", Json::String(scenario.into())),
        ("leaf_blocks", n(space.blocks.len() as f64)),
        ("full_scalar_dofs", n(space.ndof as f64)),
        (
            "independent_scalar_dofs",
            n(adaptive.lift.free_roots.len() as f64),
        ),
        ("trace_probes", n(probe_count as f64)),
        ("trace_position_error_mm", n(geom_err)),
        ("random_scalar_trace_jump", n(trace_jump)),
        ("affine_lift_error_mm", n(affine_lift)),
        ("normalized_rigid_residual_max", n(rigid_max)),
        ("affine_solution_error_relative", n(patch_err)),
        ("affine_energy_error_relative", n(energy_err)),
        ("tip_y_mm", n(tip)),
        ("loaded_strain_energy_n_mm", n(loaded.strain_energy_n_mm)),
        ("loaded_work_error_relative", n(work_error)),
        ("loaded_residual_relative", n(loaded.relative_residual)),
        ("uniform_reference_relative", comparison),
        ("build_and_trace_s", n(startup)),
        ("assembly_s", n(assembly)),
        ("patch_s", n(patch_s)),
        ("loaded_s", n(solve_s)),
    ]))
}

fn main() -> Result<(), &'static str> {
    let mut records = Vec::new();
    for rx in [64, 56] {
        for scenario in ["h0", "local1", "local2", "h1"] {
            eprintln!("adaptive {rx} {scenario}");
            records.push(run(rx, scenario)?);
        }
    }
    println!("{}", Json::Array(records));
    Ok(())
}
