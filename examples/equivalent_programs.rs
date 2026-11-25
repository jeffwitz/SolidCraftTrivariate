//! Compare equivalent extruded and axially-unioned G2 programs with a cubic manufactured field.
use solidcraft_trivariate::{
    analysis_space::{
        axial_boolean_analysis_space, section_analysis_space, AnalysisBlock, AnalysisSpace,
    },
    axial_boolean::{
        build_axial_boolean, AxialAxis, AxialBooleanOperation as Op, AxialOperand,
        AxialSectionFamily,
    },
    boundary_traction::{integrate_traction, LoadPatch},
    constrained_elasticity::solve_constrained,
    constraint_system::{resolve, Ratio},
    elasticity::{assemble, basis_geometry, gauss, BasisSample, Vector},
    h_refine::refine_uniform,
    incidence::BlockSide,
    json::{object, parse, Json},
    section_audit::audit_section,
    sections::section_blocks,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    time::Instant,
};

const YOUNG: f64 = 210_000.0;
const POISSON: f64 = 0.3;
const LAMBDA: f64 = YOUNG * POISSON / ((1.0 + POISSON) * (1.0 - 2.0 * POISSON));
const MU: f64 = YOUNG / (2.0 * (1.0 + POISSON));
const AMPLITUDE: f64 = 1.0e-4;
const SECTION_AREA: f64 = 400.0 - 16.0 * std::f64::consts::PI;

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
    name: &'static str,
    initial: AnalysisSpace,
    groups: BTreeMap<String, Vec<BlockSide>>,
    bottom_group: &'static str,
}

fn canonical_cap_group(source: &str) -> Result<String, &'static str> {
    match source {
        "z0" => Ok("bottom".into()),
        "z8" => Ok("top".into()),
        _ => Err("unknown axial cap provenance"),
    }
}

fn make_programs() -> Result<(Program, Program, f64), &'static str> {
    let family = "pierced_square_v1";
    let section = AxialSectionFamily::PiercedSquareV1;
    let source = section_blocks(family)?;
    let certificate = audit_section(family, &source)?;
    let a_space = section_analysis_space(family, &source)?;
    let mut a_groups = BTreeMap::<String, Vec<BlockSide>>::new();
    for (side, tag) in &certificate.design_faces {
        let key = if side.axis == 2 {
            canonical_cap_group(tag)?
        } else {
            tag.clone()
        };
        a_groups.entry(key).or_default().push(side.clone());
    }
    let a = Program {
        name: "extrusion_g2",
        initial: a_space,
        groups: a_groups,
        bottom_group: "z0",
    };

    let operands = [
        AxialOperand {
            section,
            axis: AxialAxis::Z,
            interval_mm: [0, 5],
        },
        AxialOperand {
            section,
            axis: AxialAxis::Z,
            interval_mm: [3, 8],
        },
    ];
    let model = build_axial_boolean(Op::Union, operands[0], operands[1])
        .map_err(|_| "axial union construction failed")?;
    let b_space = axial_boolean_analysis_space(Op::Union, operands[0], operands[1], &model)?;
    let mut b_groups = BTreeMap::<String, Vec<BlockSide>>::new();
    for face in &model.faces {
        let group = if face.side.axis == 2 {
            let names: BTreeSet<_> = face
                .provenance
                .iter()
                .map(|p| p.source_face.as_str())
                .collect();
            if names.len() != 1 {
                return Err("ambiguous cap source labels");
            }
            canonical_cap_group(names.into_iter().next().unwrap())?
        } else {
            let names: BTreeSet<_> = face
                .provenance
                .iter()
                .map(|p| p.source_face.as_str())
                .collect();
            if names.len() != 1 {
                return Err("ambiguous lateral source labels");
            }
            names.into_iter().next().unwrap().to_owned()
        };
        b_groups.entry(group).or_default().push(face.side.clone());
    }
    let b = Program {
        name: "axial_union_g6",
        initial: b_space,
        groups: b_groups,
        bottom_group: "bottom",
    };
    if a.groups.keys().collect::<Vec<_>>() != b.groups.keys().collect::<Vec<_>>() {
        return Err("physical design-face selections differ");
    }
    if a.initial.blocks.len() != 4
        || b.initial.blocks.len() != 12
        || model.selected_intervals.len() != 3
    {
        return Err("program block inventory differs from the M4 recipe");
    }
    Ok((a, b, model.volume_mm3))
}

fn refine(program: &Program, level: usize) -> Result<AnalysisSpace, &'static str> {
    let mut space = program.initial.clone();
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
    let blocks: BTreeMap<_, _> = space.blocks.iter().map(|b| (b.id.as_str(), b)).collect();
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
                let boundary_bit = u8::from(base.upper);
                let on_parent_boundary = path
                    .split('/')
                    .filter(|p| !p.is_empty())
                    .all(|p| p.as_bytes().get(base.axis).copied() == Some(b'0' + boundary_bit));
                if on_parent_boundary {
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
    if result.len() != program.groups.len()
        || result.values().any(Vec::is_empty)
        || blocks.len() != space.blocks.len()
    {
        return Err("failed to propagate exterior selections");
    }
    Ok(result)
}

fn surface_load(
    space: &AnalysisSpace,
    groups: &BTreeMap<String, Vec<BlockSide>>,
    order: usize,
) -> Result<(Vec<f64>, f64, Vector), &'static str> {
    let exterior: BTreeSet<_> = groups.values().flatten().cloned().collect();
    let patches: Vec<_> = groups
        .iter()
        .filter(|(group, _)| group.as_str() != "bottom")
        .flat_map(|(_, sides)| sides.iter().cloned().map(LoadPatch::whole))
        .collect();
    let load = integrate_traction(space, &exterior, &patches, order, |x, n| {
        let ezz = 3.0 * AMPLITUDE * x[2] * x[2];
        let stress = [LAMBDA * ezz, LAMBDA * ezz, (LAMBDA + 2.0 * MU) * ezz];
        [stress[0] * n[0], stress[1] * n[1], stress[2] * n[2]]
    })?;
    Ok((load.nodal_force_n, load.area_mm2, load.resultant_n))
}

fn body_force(space: &AnalysisSpace, order: usize) -> Result<(Vec<f64>, Vector), &'static str> {
    let rule = gauss(order)?;
    let mut force = vec![0.; 3 * space.ndof];
    for block in &space.blocks {
        for &(u, wu) in &rule {
            for &(v, wv) in &rule {
                for &(w, ww) in &rule {
                    let sample = basis_geometry(block, [u, v, w])?;
                    let density = -6.0 * AMPLITUDE * (LAMBDA + 2.0 * MU) * sample.position[2];
                    let measure = sample.determinant * wu * wv * ww;
                    for (local, &dof) in block.dofs.iter().enumerate() {
                        force[3 * dof + 2] += sample.values[local] * density * measure;
                    }
                }
            }
        }
    }
    let resultant = [0., 0., force.iter().skip(2).step_by(3).sum()];
    Ok((force, resultant))
}

fn analytical_load_work(
    space: &AnalysisSpace,
    groups: &BTreeMap<String, Vec<BlockSide>>,
    order: usize,
) -> Result<f64, &'static str> {
    let rule = gauss(order)?;
    let mut work = 0.0;
    for (group, sides) in groups {
        if group == "bottom" {
            continue;
        }
        for side in sides {
            let block = space
                .blocks
                .iter()
                .find(|b| b.id == side.block)
                .ok_or("missing exterior block")?;
            let axes: [usize; 2] = match side.axis {
                0 => [1, 2],
                1 => [0, 2],
                _ => [0, 1],
            };
            let normal_sign = if side.upper ^ (side.axis == 1) {
                1.0
            } else {
                -1.0
            };
            for &(u, wu) in &rule {
                for &(v, wv) in &rule {
                    let mut uvw = [0.0; 3];
                    uvw[side.axis] = f64::from(u8::from(side.upper));
                    uvw[axes[0]] = u;
                    uvw[axes[1]] = v;
                    let sample = basis_geometry(block, uvw)?;
                    let area_vector = cross(
                        sample.jacobian_columns[axes[0]],
                        sample.jacobian_columns[axes[1]],
                    );
                    let area = norm(area_vector);
                    if !area.is_finite() || area <= 0.0 {
                        return Err("invalid natural-boundary area element");
                    }
                    let normal = area_vector.map(|x| normal_sign * x / area);
                    let z = sample.position[2];
                    let exact_u = AMPLITUDE * z * z * z;
                    let ezz = 3.0 * AMPLITUDE * z * z;
                    let traction = [
                        LAMBDA * ezz * normal[0],
                        LAMBDA * ezz * normal[1],
                        (LAMBDA + 2.0 * MU) * ezz * normal[2],
                    ];
                    work += traction[2] * exact_u * area * wu * wv;
                }
            }
        }
    }
    for block in &space.blocks {
        for &(u, wu) in &rule {
            for &(v, wv) in &rule {
                for &(w, ww) in &rule {
                    let sample = basis_geometry(block, [u, v, w])?;
                    let z = sample.position[2];
                    let density = -6.0 * AMPLITUDE * (LAMBDA + 2.0 * MU) * z;
                    let exact_u = AMPLITUDE * z * z * z;
                    work += density * exact_u * sample.determinant * wu * wv * ww;
                }
            }
        }
    }
    if !work.is_finite() {
        return Err("nonfinite manufactured load work");
    }
    Ok(work)
}

fn strain(grad: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| 0.5 * (grad[i][j] + grad[j][i])))
}
fn elastic_product(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> f64 {
    let tra = a[0][0] + a[1][1] + a[2][2];
    let trb = b[0][0] + b[1][1] + b[2][2];
    LAMBDA * tra * trb
        + 2.0
            * MU
            * (0..3)
                .flat_map(|i| (0..3).map(move |j| 2.0 * a[i][j] * b[i][j]))
                .sum::<f64>()
}

#[derive(Clone)]
struct Metrics {
    l2_relative: f64,
    energy_relative: f64,
    exact_l2: f64,
    exact_energy: f64,
}
fn field_metrics(
    space: &AnalysisSpace,
    displacement: &[f64],
    order: usize,
) -> Result<Metrics, &'static str> {
    let rule = gauss(order)?;
    let mut sums = [0.; 4];
    for block in &space.blocks {
        for &(u, wu) in &rule {
            for &(v, wv) in &rule {
                for &(w, ww) in &rule {
                    let (value, grad, s) = eval_value_gradient(block, displacement, [u, v, w])?;
                    let z = s.position[2];
                    let exact = AMPLITUDE * z * z * z;
                    let mut exact_grad = [[0.; 3]; 3];
                    exact_grad[2][2] = 3.0 * AMPLITUDE * z * z;
                    let diff_grad = std::array::from_fn(|i| {
                        std::array::from_fn(|j| grad[i][j] - exact_grad[i][j])
                    });
                    let eps_diff = strain(diff_grad);
                    let eps_exact = strain(exact_grad);
                    let d = [value[0], value[1], value[2] - exact];
                    let measure = s.determinant * wu * wv * ww;
                    sums[0] += dot(d, d) * measure;
                    sums[1] += elastic_product(eps_diff, eps_diff) * measure;
                    sums[2] += exact * exact * measure;
                    sums[3] += elastic_product(eps_exact, eps_exact) * measure;
                }
            }
        }
    }
    if sums.iter().any(|x| !x.is_finite() || *x < 0.) {
        return Err("nonfinite field error metric");
    }
    Ok(Metrics {
        l2_relative: (sums[0] / sums[2]).sqrt(),
        energy_relative: (sums[1] / sums[3]).sqrt(),
        exact_l2: sums[2].sqrt(),
        exact_energy: sums[3].sqrt(),
    })
}

struct Solved {
    space: AnalysisSpace,
    displacement: Vec<f64>,
    metrics: Metrics,
    volume: f64,
    loaded_area: f64,
    surface_force: Vector,
    body_force: Vector,
    work_error: f64,
    iterations: usize,
    residual: f64,
    elapsed_s: f64,
    peak_rss_kib: usize,
    boundary_areas: BTreeMap<String, f64>,
    analytical_load_work: f64,
    analytical_load_work_error: f64,
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

fn solve(program: &Program, level: usize, order: usize) -> Result<Solved, &'static str> {
    let started = Instant::now();
    let space = refine(program, level)?;
    let groups = refined_groups(program, &space)?;
    let k = assemble(&space, YOUNG, POISSON, order)?;
    let (mut force, loaded_area, surface_resultant) = surface_load(&space, &groups, order)?;
    let (body, body_resultant) = body_force(&space, order)?;
    for (f, b) in force.iter_mut().zip(body) {
        *f += b;
    }
    let prescribed = space
        .boundary_dofs
        .get(program.bottom_group)
        .ok_or("missing fixed bottom group")?
        .iter()
        .map(|&d| (d, Ratio::zero()))
        .collect::<Vec<_>>();
    let lift = resolve(space.ndof, &[], &prescribed)?;
    let solution = solve_constrained(&k, [&lift, &lift, &lift], &force, 1e-10, 30_000)?;
    let work = force
        .iter()
        .zip(&solution.displacement_mm)
        .map(|(f, u)| f * u)
        .sum::<f64>();
    let work_error = (work / (2.0 * solution.strain_energy_n_mm) - 1.0).abs();
    let metrics = field_metrics(&space, &solution.displacement_mm, order)?;
    let analytical_load_work = analytical_load_work(&space, &groups, order)?;
    let expected_work =
        SECTION_AREA * 9.0 * AMPLITUDE * AMPLITUDE * (8.0_f64.powi(5) / 5.0) * (LAMBDA + 2.0 * MU);
    let analytical_load_work_error = (analytical_load_work / expected_work - 1.0).abs();
    let net = [
        surface_resultant[0] + body_resultant[0],
        surface_resultant[1] + body_resultant[1],
        surface_resultant[2] + body_resultant[2],
    ];
    let scale = surface_resultant
        .iter()
        .chain(body_resultant.iter())
        .map(|x| x.abs())
        .fold(0.0, f64::max);
    let exact_axial_force = 192.0 * AMPLITUDE * (LAMBDA + 2.0 * MU) * SECTION_AREA;
    if work_error > 1e-7
        || analytical_load_work_error > 1e-7
        || (surface_resultant[2] / exact_axial_force - 1.0).abs() > 1e-7
        || (body_resultant[2] / exact_axial_force + 1.0).abs() > 1e-7
        || norm(net) > 1e-7 * scale.max(1.0)
    {
        return Err("manufactured load or work balance failed");
    }
    let boundary_areas = group_areas(&groups, &space, order)?;
    Ok(Solved {
        space,
        displacement: solution.displacement_mm,
        metrics,
        volume: k.integrated_volume,
        loaded_area,
        surface_force: surface_resultant,
        body_force: body_resultant,
        work_error,
        iterations: solution.iterations,
        residual: solution.relative_residual,
        elapsed_s: started.elapsed().as_secs_f64(),
        peak_rss_kib: peak_rss_kib(),
        boundary_areas,
        analytical_load_work,
        analytical_load_work_error,
    })
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct PlaneKey {
    sector: String,
    xy_path: Vec<[u8; 2]>,
}
fn group_key(id: &str) -> Result<PlaneKey, &'static str> {
    let mut parts = id.split('/');
    let root = parts.next().ok_or("empty analysis block ID")?;
    let marker = root
        .find("pierced_sector_")
        .ok_or("block is not from pierced section")?
        + "pierced_sector_".len();
    let sector = root[marker..]
        .split("_z")
        .next()
        .ok_or("missing sector suffix")?
        .to_owned();
    let xy_path = parts
        .map(|p| {
            if p.len() != 3 || !p.as_bytes().iter().all(|b| matches!(b, b'0' | b'1')) {
                return Err("invalid refined block path");
            }
            Ok([p.as_bytes()[0] - b'0', p.as_bytes()[1] - b'0'])
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(PlaneKey { sector, xy_path })
}
fn z_bounds(block: &AnalysisBlock) -> Result<[f64; 2], &'static str> {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for p in &block.controls {
        if p[3] <= 0 {
            return Err("nonpositive analysis weight");
        }
        let z = p[2] as f64 / p[3] as f64;
        lo = lo.min(z);
        hi = hi.max(z);
    }
    if !(lo.is_finite() && hi.is_finite() && lo < hi) {
        return Err("invalid physical z span");
    }
    Ok([lo, hi])
}
type SampledField = ([f64; 3], [[f64; 3]; 3], BasisSample);

fn eval_value_gradient(
    block: &AnalysisBlock,
    disp: &[f64],
    uvw: [f64; 3],
) -> Result<SampledField, &'static str> {
    let sample = basis_geometry(block, uvw)?;
    let mut value = [0.; 3];
    let mut grad = [[0.; 3]; 3];
    for (i, &dof) in block.dofs.iter().enumerate() {
        for c in 0..3 {
            value[c] += sample.values[i] * disp[3 * dof + c];
            for (derivative, basis) in grad[c].iter_mut().zip(sample.gradients[i]) {
                *derivative += basis * disp[3 * dof + c];
            }
        }
    }
    Ok((value, grad, sample))
}

fn common_field_distance(a: &Solved, b: &Solved, order: usize) -> Result<(f64, f64), &'static str> {
    let mut ga = BTreeMap::<PlaneKey, Vec<(&AnalysisBlock, [f64; 2])>>::new();
    let mut gb = BTreeMap::<PlaneKey, Vec<(&AnalysisBlock, [f64; 2])>>::new();
    for block in &a.space.blocks {
        ga.entry(group_key(&block.id)?)
            .or_default()
            .push((block, z_bounds(block)?));
    }
    for block in &b.space.blocks {
        gb.entry(group_key(&block.id)?)
            .or_default()
            .push((block, z_bounds(block)?));
    }
    if ga.keys().collect::<Vec<_>>() != gb.keys().collect::<Vec<_>>() {
        return Err("programs have different common cross-section leaf partitions");
    }
    let rule = gauss(order)?;
    let mut difference = [0.; 2];
    for (key, aa) in &ga {
        let bb = &gb[key];
        let mut cuts: Vec<f64> = aa
            .iter()
            .chain(bb)
            .flat_map(|(_, z)| [z[0], z[1]])
            .collect();
        cuts.sort_by(f64::total_cmp);
        cuts.dedup_by(|x, y| x.to_bits() == y.to_bits());
        for span in cuts.windows(2) {
            let [z0, z1] = [span[0], span[1]];
            if z1 <= z0 {
                continue;
            }
            let zm = 0.5 * (z0 + z1);
            let (ba, za) = aa
                .iter()
                .find(|(_, z)| z[0] <= zm && zm < z[1])
                .ok_or("gap in extrusion program")?;
            let (bb, zb) = bb
                .iter()
                .find(|(_, z)| z[0] <= zm && zm < z[1])
                .ok_or("gap in axial union program")?;
            for &(u, wu) in &rule {
                for &(v, wv) in &rule {
                    for &(q, wq) in &rule {
                        let z = z0 + (z1 - z0) * q;
                        let wa = (z - za[0]) / (za[1] - za[0]);
                        let wb = (z - zb[0]) / (zb[1] - zb[0]);
                        let (ua, ga_, sa) = eval_value_gradient(ba, &a.displacement, [u, v, wa])?;
                        let (ub, gb_, _sb) = eval_value_gradient(bb, &b.displacement, [u, v, wb])?;
                        let area = norm(cross(sa.jacobian_columns[0], sa.jacobian_columns[1]));
                        let measure = area * (z1 - z0) * wu * wv * wq;
                        let dv = [ua[0] - ub[0], ua[1] - ub[1], ua[2] - ub[2]];
                        let gd =
                            std::array::from_fn(|i| std::array::from_fn(|j| ga_[i][j] - gb_[i][j]));
                        let ed = strain(gd);
                        difference[0] += dot(dv, dv) * measure;
                        difference[1] += elastic_product(ed, ed) * measure;
                    }
                }
            }
        }
    }
    let exact_l2_sq = a.metrics.exact_l2 * a.metrics.exact_l2;
    let exact_energy_sq = a.metrics.exact_energy * a.metrics.exact_energy;
    Ok((
        (difference[0] / exact_l2_sq).sqrt(),
        (difference[1] / exact_energy_sq).sqrt(),
    ))
}

fn group_areas(
    groups: &BTreeMap<String, Vec<BlockSide>>,
    space: &AnalysisSpace,
    order: usize,
) -> Result<BTreeMap<String, f64>, &'static str> {
    let exterior: BTreeSet<_> = groups.values().flatten().cloned().collect();
    let mut result = BTreeMap::new();
    for (group, sides) in groups {
        let patches = sides
            .iter()
            .cloned()
            .map(LoadPatch::whole)
            .collect::<Vec<_>>();
        let measure =
            integrate_traction(space, &exterior, &patches, order, |_, _| [0.; 3])?.area_mm2;
        result.insert(group.clone(), measure);
    }
    Ok(result)
}

fn json_solve(level: usize, q: usize, s: &Solved) -> Json {
    object([
        ("level", integer(level)),
        ("quadrature", integer(q)),
        ("blocks", integer(s.space.blocks.len())),
        ("scalar_dofs", integer(s.space.ndof)),
        ("relative_l2_error", num(s.metrics.l2_relative)),
        ("relative_energy_error", num(s.metrics.energy_relative)),
        ("volume_mm3", num(s.volume)),
        ("loaded_surface_area_mm2", num(s.loaded_area)),
        ("body_resultant_n_z", num(s.body_force[2])),
        ("surface_resultant_n_z", num(s.surface_force[2])),
        ("work_error_relative", num(s.work_error)),
        ("analytic_load_work_n_mm", num(s.analytical_load_work)),
        (
            "analytic_load_work_error_relative",
            num(s.analytical_load_work_error),
        ),
        ("linear_residual_relative", num(s.residual)),
        ("cg_iterations", integer(s.iterations)),
        ("wall_s", num(s.elapsed_s)),
        ("peak_rss_kib", integer(s.peak_rss_kib)),
        (
            "boundary_group_areas_mm2",
            Json::Object(
                s.boundary_areas
                    .iter()
                    .map(|(k, v)| (k.clone(), num(*v)))
                    .collect(),
            ),
        ),
    ])
}

fn main() -> Result<(), &'static str> {
    let campaign_started = Instant::now();
    let reference = parse(include_str!("../tests/references/equivalent_programs.json"))?;
    let ref_decimal = |key: &str| -> Result<f64, &'static str> {
        reference
            .field(key)?
            .text()?
            .parse()
            .map_err(|_| "invalid Decimal reference value")
    };
    let (a, b, analytic_volume) = make_programs()?;
    let exact_l2_squared = SECTION_AREA * AMPLITUDE.powi(2) * 8.0_f64.powi(7) / 7.0;
    let exact_energy_squared =
        SECTION_AREA * 9.0 * AMPLITUDE.powi(2) * (8.0_f64.powi(5) / 5.0) * (LAMBDA + 2.0 * MU);
    let expected_top_force = 192.0 * AMPLITUDE * (LAMBDA + 2.0 * MU) * SECTION_AREA;
    let expected_external_work = exact_energy_squared;
    let independent_errors = [
        (
            "section_area",
            (SECTION_AREA / ref_decimal("section_area_mm2")? - 1.0).abs(),
        ),
        (
            "volume",
            (analytic_volume / ref_decimal("volume_mm3")? - 1.0).abs(),
        ),
        (
            "L2_squared",
            (exact_l2_squared / ref_decimal("l2_squared_mm5")? - 1.0).abs(),
        ),
        (
            "energy",
            (0.5 * exact_energy_squared / ref_decimal("energy_n_mm")? - 1.0).abs(),
        ),
        (
            "external_work",
            (expected_external_work / ref_decimal("external_work_n_mm")? - 1.0).abs(),
        ),
        (
            "top_force",
            (expected_top_force / ref_decimal("top_force_z_n")? - 1.0).abs(),
        ),
        (
            "body_force",
            (-expected_top_force / ref_decimal("body_force_z_n")? - 1.0).abs(),
        ),
        (
            "tip",
            (AMPLITUDE * 8.0_f64.powi(3) / ref_decimal("tip_z_displacement_mm")? - 1.0).abs(),
        ),
    ];
    if independent_errors.iter().any(|(_, e)| *e > 2e-12) {
        return Err("independent Decimal manufactured oracle disagrees");
    }
    let mut q8 = Vec::new();
    let mut direct = Vec::new();
    let mut q6 = Vec::new();
    let mut prior_a = None;
    let mut prior_b = None;
    for level in 0..=3 {
        eprintln!("M4 q8 h{level}: {}", a.name);
        let sa = solve(&a, level, 8)?;
        eprintln!("M4 q8 h{level}: {}", b.name);
        let sb = solve(&b, level, 8)?;
        let area_a = &sa.boundary_areas;
        let area_b = &sb.boundary_areas;
        if area_a.keys().collect::<Vec<_>>() != area_b.keys().collect::<Vec<_>>()
            || area_a
                .iter()
                .any(|(k, x)| (x - area_b[k]).abs() > 1e-9 * x.abs().max(1.))
        {
            return Err("physical boundary area groups differ");
        }
        if (sa.volume - analytic_volume).abs() > 1e-8 * analytic_volume
            || (sb.volume - analytic_volume).abs() > 1e-8 * analytic_volume
            || (sa.volume - sb.volume).abs() > 1e-10 * analytic_volume
        {
            return Err("quadrature volume differs between equivalent programs");
        }
        let (l2, e) = common_field_distance(&sa, &sb, 8)?;
        direct.push(object([("level",integer(level)),("relative_l2_field_distance",num(l2)),("relative_energy_field_distance",num(e)),("integration",Json::String("shared physical domain, common xy leaves and union of both axial cell cuts, q8".into()))]));
        q8.push(object([
            ("program", Json::String(a.name.into())),
            ("solution", json_solve(level, 8, &sa)),
            (
                "convergence_rate_l2",
                prior_a.as_ref().map_or(Json::Null, |p: &Metrics| {
                    num((p.l2_relative / sa.metrics.l2_relative).log2())
                }),
            ),
            (
                "convergence_rate_energy",
                prior_a.as_ref().map_or(Json::Null, |p: &Metrics| {
                    num((p.energy_relative / sa.metrics.energy_relative).log2())
                }),
            ),
        ]));
        q8.push(object([
            ("program", Json::String(b.name.into())),
            ("solution", json_solve(level, 8, &sb)),
            (
                "convergence_rate_l2",
                prior_b.as_ref().map_or(Json::Null, |p: &Metrics| {
                    num((p.l2_relative / sb.metrics.l2_relative).log2())
                }),
            ),
            (
                "convergence_rate_energy",
                prior_b.as_ref().map_or(Json::Null, |p: &Metrics| {
                    num((p.energy_relative / sb.metrics.energy_relative).log2())
                }),
            ),
        ]));
        prior_a = Some(sa.metrics.clone());
        prior_b = Some(sb.metrics.clone());
        if level <= 1 {
            eprintln!("M4 q6 h{level}: {}", a.name);
            let qa = solve(&a, level, 6)?;
            eprintln!("M4 q6 h{level}: {}", b.name);
            let qb = solve(&b, level, 6)?;
            q6.push(object([
                ("program", Json::String(a.name.into())),
                ("level", integer(level)),
                ("relative_l2_error", num(qa.metrics.l2_relative)),
                ("relative_energy_error", num(qa.metrics.energy_relative)),
                ("volume_mm3", num(qa.volume)),
                ("loaded_surface_area_mm2", num(qa.loaded_area)),
                ("work_error_relative", num(qa.work_error)),
                ("scalar_dofs", integer(qa.space.ndof)),
                ("cg_iterations", integer(qa.iterations)),
                ("wall_s", num(qa.elapsed_s)),
                ("peak_rss_kib", integer(qa.peak_rss_kib)),
            ]));
            q6.push(object([
                ("program", Json::String(b.name.into())),
                ("level", integer(level)),
                ("relative_l2_error", num(qb.metrics.l2_relative)),
                ("relative_energy_error", num(qb.metrics.energy_relative)),
                ("volume_mm3", num(qb.volume)),
                ("loaded_surface_area_mm2", num(qb.loaded_area)),
                ("work_error_relative", num(qb.work_error)),
                ("scalar_dofs", integer(qb.space.ndof)),
                ("cg_iterations", integer(qb.iterations)),
                ("wall_s", num(qb.elapsed_s)),
                ("peak_rss_kib", integer(qb.peak_rss_kib)),
            ]));
        }
    }
    let final_rates = q8
        .iter()
        .filter_map(|r| {
            let p = r.field("program").ok()?.text().ok()?;
            let l = r
                .field("solution")
                .ok()?
                .field("level")
                .ok()?
                .number::<usize>()
                .ok()?;
            if l != 3 {
                return None;
            }
            Some(object([
                ("program", Json::String(p.to_owned())),
                ("l2", r.field("convergence_rate_l2").ok()?.clone()),
                ("energy", r.field("convergence_rate_energy").ok()?.clone()),
            ]))
        })
        .collect::<Vec<_>>();
    let scalar_dofs = q8
        .iter()
        .map(|r| {
            r.field("solution")
                .unwrap()
                .field("scalar_dofs")
                .unwrap()
                .number::<usize>()
                .unwrap()
        })
        .chain(
            q6.iter()
                .map(|r| r.field("scalar_dofs").unwrap().number::<usize>().unwrap()),
        )
        .collect::<Vec<_>>();
    let wall_times = q8
        .iter()
        .map(|r| {
            r.field("solution")
                .unwrap()
                .field("wall_s")
                .unwrap()
                .number::<f64>()
                .unwrap()
        })
        .chain(
            q6.iter()
                .map(|r| r.field("wall_s").unwrap().number::<f64>().unwrap()),
        )
        .collect::<Vec<_>>();
    let rss = q8
        .iter()
        .map(|r| {
            r.field("solution")
                .unwrap()
                .field("peak_rss_kib")
                .unwrap()
                .number::<usize>()
                .unwrap()
        })
        .chain(
            q6.iter()
                .map(|r| r.field("peak_rss_kib").unwrap().number::<usize>().unwrap()),
        )
        .collect::<Vec<_>>();
    let output=object([
        ("schema",Json::String("equivalent-programs-m4-v1".into())),("programs",Json::Array(vec![Json::String(a.name.into()),Json::String(b.name.into())])),
        ("geometry",object([("section",Json::String("pierced_square_v1".into())),("section_area_target_mm2",num(SECTION_AREA)),("height_mm",num(8.)),("radius_mm",num(4.)),("outer_width_mm",num(20.)),("analytic_volume_mm3",num(analytic_volume)),("analytic_l2_squared_mm5",num(exact_l2_squared)),("analytic_energy_seminorm_squared_n_mm",num(exact_energy_squared)),("analytic_external_work_n_mm",num(expected_external_work)),("analytic_top_force_z_n",num(expected_top_force)),("axial_union_layers_mm",Json::Array(vec![Json::Array(vec![integer(0),integer(3)]),Json::Array(vec![integer(3),integer(5)]),Json::Array(vec![integer(5),integer(8)])]))])),
        ("physics",object([("young_mpa",num(YOUNG)),("poisson",num(POISSON)),("manufactured_displacement",Json::String("u=(0,0,1e-4 z^3) mm".into())),("body_force",Json::String("f=-div(sigma)".into())),("traction",Json::String("sigma n on every exterior face except clamped z=0; internal caps excluded".into())),("quadrature",Json::String("Gauss-Legendre tensor product".into()))])),
        ("q8_h0_to_h3",Json::Array(q8)),("q6_h0_h1",Json::Array(q6)),("direct_field_distance_q8",Json::Array(direct)),("final_rates_h2_to_h3",Json::Array(final_rates)),
        ("independent_decimal_reference_max_relative_errors",Json::Object(independent_errors.iter().map(|(k,v)|(k.to_string(),num(*v))).collect())),
        ("resource_ranges",object([("scalar_dofs_minmax",Json::Array(vec![integer(*scalar_dofs.iter().min().unwrap()),integer(*scalar_dofs.iter().max().unwrap())])),("wall_s_minmax",Json::Array(vec![num(wall_times.iter().copied().fold(f64::INFINITY,f64::min)),num(wall_times.iter().copied().fold(0.,f64::max))])),("rss_kib_minmax",Json::Array(vec![integer(*rss.iter().min().unwrap()),integer(*rss.iter().max().unwrap())])),("resource_record_scope",Json::String("all q8 h0-h3 and q6 h0-h1 solves; RSS is process high-water mark observed after each solve, hence cumulative".into())),("campaign_wall_s",num(campaign_started.elapsed().as_secs_f64()))])),
        ("qualification",Json::String("Finite-mesh fields are reported separately; equality of the audited domain does not imply equality of the discrete solutions. Rates are observed on these four uniform levels only.".into()))
    ]);
    fs::write(
        "docs/results/equivalent-programs-20261125.json",
        format!("{output}\n"),
    )
    .map_err(|_| "cannot write M4 result file")?;
    println!("{output}");
    Ok(())
}
