//! Independent continuous audit for circular bores with exact axial intervals.
//!
//! The audit rebuilds the G18 guard triangulation, rational circle/guard and
//! circle/inner-diamond sectors, and the four quadrilaterals filling each
//! material disk. It certifies every planar Jacobian on the closed parameter
//! square and checks exact whole-face incidence on the full layered complex.
//! The proof is limited to this fixed guard-and-sector template.

use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    circular_polygon::CircularPolygonBore,
    incidence::{enumerate_faces, Incidence},
    interval_bores::{IntervalBore, IntervalBorePrism},
    perforated_prism::{LoopProfile, PerforatedPrism},
    perforated_prism_audit::audit_perforated_prism,
    planar::{certify_planar_refined, PlanarError},
    polygon_prism::Vertex,
    polygon_prism_audit::{expected_controls, next_down, parse_triangle_id},
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq)]
pub struct IntervalBoreCertificate {
    pub certificate: BracketCertificate,
    pub outer_area_twice_mm2: i64,
    pub removed_pi_volume_coefficient_mm3: i64,
    pub enclosed_voids: usize,
    pub through_bores: usize,
}

const QUADRANTS: [&str; 4] = ["east", "north", "west", "south"];

#[derive(Clone)]
enum Role {
    Outer {
        base_id: String,
        z0: i32,
        z1: i32,
    },
    Wall {
        bore: String,
        active: bool,
        z0: i32,
        z1: i32,
    },
    Fill {
        bore: String,
        z0: i32,
        z1: i32,
    },
    Core {
        bore: String,
        z0: i32,
        z1: i32,
    },
}

fn add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("interval_bore_audit_overflow")
}

fn mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("interval_bore_audit_overflow")
}

fn primitive(values: &[[i128; 4]]) -> Result<Vec<[i128; 4]>, &'static str> {
    let gcd = |mut a: u128, mut b: u128| {
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a
    };
    let common = values
        .iter()
        .flatten()
        .fold(0_u128, |g, x| gcd(g, x.unsigned_abs()));
    if common == 0 {
        return Err("interval_bore_zero_homogeneous_net");
    }
    let divisor = i128::try_from(common).map_err(|_| "interval_bore_audit_overflow")?;
    Ok(values.iter().map(|p| p.map(|x| x / divisor)).collect())
}

fn primitive_block(block: &BracketBlock) -> Result<Vec<[i128; 4]>, &'static str> {
    if block.controls.iter().any(|p| p[3] <= 0) {
        return Err("interval_bore_nonpositive_weight");
    }
    primitive(
        &block
            .controls
            .iter()
            .map(|p| p.map(i128::from))
            .collect::<Vec<_>>(),
    )
}

fn expected_block(
    id: String,
    degrees: [usize; 3],
    controls: Vec<[i128; 4]>,
) -> Result<BracketBlock, &'static str> {
    let controls = primitive(&controls)?;
    Ok(BracketBlock {
        id,
        degrees,
        controls: controls
            .into_iter()
            .map(|p| {
                Ok([
                    i64::try_from(p[0]).map_err(|_| "interval_bore_audit_overflow")?,
                    i64::try_from(p[1]).map_err(|_| "interval_bore_audit_overflow")?,
                    i64::try_from(p[2]).map_err(|_| "interval_bore_audit_overflow")?,
                    i64::try_from(p[3]).map_err(|_| "interval_bore_audit_overflow")?,
                ])
            })
            .collect::<Result<Vec<_>, _>>()?,
    })
}

fn signed_area2(vertices: &[Vertex]) -> i64 {
    (0..vertices.len())
        .map(|i| {
            let a = vertices[i].xy_mm;
            let b = vertices[(i + 1) % vertices.len()].xy_mm;
            i64::from(a[0]) * i64::from(b[1]) - i64::from(a[1]) * i64::from(b[0])
        })
        .sum()
}

fn guard_loop(bore: &CircularPolygonBore) -> Result<LoopProfile, &'static str> {
    let [cx, cy] = bore.center_mm.map(i64::from);
    let g = i64::from(bore.guard_mm);
    let points = [[cx + g, cy], [cx, cy + g], [cx - g, cy], [cx, cy - g]];
    if points.iter().flatten().any(|&x| !(-64..=64).contains(&x)) {
        return Err("interval_bore_guard_outside_domain");
    }
    Ok(LoopProfile {
        id: format!("guard_{}", bore.id),
        vertices: points
            .into_iter()
            .zip(["e", "n", "w", "s"])
            .map(|(xy_mm, d)| Vertex {
                id: format!("g_{}_{}", bore.id, d),
                xy_mm: xy_mm.map(|x| x as i32),
            })
            .collect(),
    })
}

fn normalized_loops(recipe: &IntervalBorePrism) -> Result<Vec<LoopProfile>, &'static str> {
    let mut loops = vec![recipe.outer.clone()];
    for bore in &recipe.bores {
        loops.push(guard_loop(&bore.bore)?);
    }
    for (i, loop_) in loops.iter_mut().enumerate() {
        let area = signed_area2(&loop_.vertices);
        if area == 0 {
            return Err("interval_bore_zero_area_profile");
        }
        if (area > 0) != (i == 0) {
            loop_.vertices.reverse();
        }
    }
    Ok(loops)
}

fn rotate_ring_profile(point: [i128; 3], quadrant: usize) -> [i128; 3] {
    let [x, y, w] = point;
    match quadrant {
        0 => [x, y, w],
        1 => [-x, y, w],
        2 => [-x, -y, w],
        _ => [x, -y, w],
    }
}

fn ring_source_index(u: usize, quadrant: usize) -> usize {
    if quadrant % 2 == 1 {
        2 - u
    } else {
        u
    }
}

#[derive(Clone, Copy)]
enum SectorKind {
    GuardToCircle,
    CircleToDiamond,
}

fn sector_profiles(bore: &CircularPolygonBore, kind: SectorKind) -> [[[i128; 3]; 2]; 3] {
    let r = i128::from(bore.radius_mm);
    let g = i128::from(bore.guard_mm);
    let guard = [[2 * g, 0, 2], [g, g, 2], [0, 2 * g, 2]];
    let circle = [[2 * r, 0, 2], [2 * r, 2 * r, 2], [0, 4 * r, 4]];
    // Degree-two elevation of the straight inner-diamond edge, with constant
    // weight so its restricted half-traces match the linear core patches.
    let diamond = [[2 * r, 0, 4], [r, r, 4], [0, 2 * r, 4]];
    std::array::from_fn(|u| match kind {
        SectorKind::GuardToCircle => [guard[u], circle[u]],
        SectorKind::CircleToDiamond => [circle[u], diamond[u]],
    })
}

fn power_quadratic(values: [i128; 3]) -> Result<[i128; 3], &'static str> {
    Ok([
        values[0],
        mul(
            2,
            values[1]
                .checked_sub(values[0])
                .ok_or("interval_bore_audit_overflow")?,
        )?,
        add(values[0], values[2])?
            .checked_sub(mul(2, values[1])?)
            .ok_or("interval_bore_audit_overflow")?,
    ])
}

fn polynomial_product(a: &[i128], b: &[i128]) -> Result<Vec<i128>, &'static str> {
    let mut product = vec![0_i128; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            product[i + j] = add(product[i + j], mul(x, y)?)?;
        }
    }
    Ok(product)
}

/// Prove the rational quadratic circle identity and strictly increasing arc
/// angle directly from its integer power coefficients.
fn prove_circle_profile(profile: &[[i128; 3]; 3], radius: i32) -> Result<(), &'static str> {
    let x = power_quadratic(std::array::from_fn(|i| profile[i][0]))?;
    let y = power_quadratic(std::array::from_fn(|i| profile[i][1]))?;
    let w = power_quadratic(std::array::from_fn(|i| profile[i][2]))?;
    let xx = polynomial_product(&x, &x)?;
    let yy = polynomial_product(&y, &y)?;
    let ww = polynomial_product(&w, &w)?;
    let r2 = mul(i128::from(radius), i128::from(radius))?;
    for k in 0..5 {
        if add(xx[k], yy[k])? != mul(r2, ww[k])? {
            return Err("interval_bore_not_exact_circle");
        }
    }
    let mut angular = [0_i128; 4];
    for i in 0..3 {
        for j in 1..3 {
            let k = i + j - 1;
            angular[k] = add(angular[k], mul(x[i], mul(j as i128, y[j])?)?)?;
            angular[k] = angular[k]
                .checked_sub(mul(y[i], mul(j as i128, x[j])?)?)
                .ok_or("interval_bore_audit_overflow")?;
        }
    }
    let bernstein_times_three = [
        mul(3, angular[0])?,
        add(mul(3, angular[0])?, angular[1])?,
        add(add(mul(3, angular[0])?, mul(2, angular[1])?)?, angular[2])?,
        mul(
            3,
            add(add(angular[0], angular[1])?, add(angular[2], angular[3])?)?,
        )?,
    ];
    if bernstein_times_three.iter().any(|&v| v <= 0) {
        return Err("interval_bore_nonmonotone_circle");
    }
    Ok(())
}

fn prove_bore_arcs(bore: &CircularPolygonBore) -> Result<(), &'static str> {
    for quadrant in 0..4 {
        for (kind, circle_v) in [
            (SectorKind::GuardToCircle, 1),
            (SectorKind::CircleToDiamond, 0),
        ] {
            let profiles = sector_profiles(bore, kind);
            let profile = std::array::from_fn(|u| {
                let source = ring_source_index(u, quadrant);
                rotate_ring_profile(profiles[source][circle_v], quadrant)
            });
            prove_circle_profile(&profile, bore.radius_mm)?;
        }
    }
    Ok(())
}

fn restricted(a: [i128; 4], b: [i128; 4], c: [i128; 4], half: usize, i: usize) -> [i128; 4] {
    match (half, i) {
        (0, 0) => a.map(|x| 4 * x),
        (0, 1) => std::array::from_fn(|k| 2 * (a[k] + b[k])),
        (0, 2) | (1, 0) => std::array::from_fn(|k| a[k] + 2 * b[k] + c[k]),
        (1, 1) => std::array::from_fn(|k| 2 * (b[k] + c[k])),
        (1, 2) => c.map(|x| 4 * x),
        _ => unreachable!("validated half-sector index"),
    }
}

fn expected_sector(
    bore: &CircularPolygonBore,
    recipe: &IntervalBorePrism,
    z0: i32,
    z1: i32,
    quadrant: usize,
    half: usize,
    kind: SectorKind,
) -> Result<BracketBlock, &'static str> {
    let profiles = sector_profiles(bore, kind);
    let [cx, cy] = bore.center_mm.map(i128::from);
    let h = i128::from(recipe.height_mm);
    let mut net = [[[[0_i128; 4]; 2]; 2]; 3];
    for (u, row) in net.iter_mut().enumerate() {
        let source_u = ring_source_index(u, quadrant);
        for (v, pair) in row.iter_mut().enumerate() {
            let [rx, ry, w] = rotate_ring_profile(profiles[source_u][v], quadrant);
            for (k, z) in [z0, z1].into_iter().enumerate() {
                let x = add(mul(cx, w)?, rx)?;
                let y = add(mul(cy, w)?, ry)?;
                pair[k] = [
                    add(
                        mul(h, x)?,
                        mul(mul(i128::from(z), i128::from(recipe.top_shift_mm[0]))?, w)?,
                    )?,
                    add(
                        mul(h, y)?,
                        mul(mul(i128::from(z), i128::from(recipe.top_shift_mm[1]))?, w)?,
                    )?,
                    mul(mul(h, i128::from(z))?, w)?,
                    mul(h, w)?,
                ];
            }
        }
    }
    let mut controls = Vec::with_capacity(12);
    for u in 0..3 {
        for ((row0, row1), row2) in net[0].iter().zip(&net[1]).zip(&net[2]) {
            for ((&a, &b), &c) in row0.iter().zip(row1).zip(row2) {
                controls.push(restricted(a, b, c, half, u));
            }
        }
    }
    let prefix = match kind {
        SectorKind::GuardToCircle => format!("bore_{}_sector", bore.id),
        SectorKind::CircleToDiamond => format!("bore_{}_fill_sector", bore.id),
    };
    expected_block(
        format!("{prefix}_{}_half_{half}_z{z0}_{z1}", QUADRANTS[quadrant]),
        [2, 1, 1],
        controls,
    )
}

fn expected_core(
    bore: &CircularPolygonBore,
    recipe: &IntervalBorePrism,
    z0: i32,
    z1: i32,
    quadrant: usize,
) -> Result<BracketBlock, &'static str> {
    let r = i128::from(bore.radius_mm);
    // Each cell is the convex quadrilateral around one cardinal vertex of the
    // inner diamond; neighboring cells meet on center-to-edge-midpoint lines.
    let e = [2 * r, 0, 4];
    let ne = [r, r, 4];
    let nw = [-r, r, 4];
    let w = [-2 * r, 0, 4];
    let sw = [-r, -r, 4];
    let se = [r, -r, 4];
    let n = [0, 2 * r, 4];
    let s = [0, -2 * r, 4];
    let origin = [0, 0, 4];
    let points = match quadrant {
        // u is the slow tensor index and v follows it, so these are
        // [P00, P01, P10, P11] in the flattened control order.
        0 => [e, se, ne, origin],
        1 => [n, ne, nw, origin],
        2 => [w, nw, sw, origin],
        _ => [s, sw, se, origin],
    };
    let [cx, cy] = bore.center_mm.map(i128::from);
    let h = i128::from(recipe.height_mm);
    let mut controls = Vec::with_capacity(8);
    for [x, y, weight] in points {
        for z in [z0, z1] {
            controls.push([
                add(
                    mul(h, add(mul(cx, weight)?, x)?)?,
                    mul(
                        mul(i128::from(z), i128::from(recipe.top_shift_mm[0]))?,
                        weight,
                    )?,
                )?,
                add(
                    mul(h, add(mul(cy, weight)?, y)?)?,
                    mul(
                        mul(i128::from(z), i128::from(recipe.top_shift_mm[1]))?,
                        weight,
                    )?,
                )?,
                mul(mul(h, i128::from(z))?, weight)?,
                mul(h, weight)?,
            ]);
        }
    }
    expected_block(
        format!(
            "bore_{}_core_sector_{}_z{z0}_{z1}",
            bore.id, QUADRANTS[quadrant]
        ),
        [1, 1, 1],
        controls,
    )
}

fn expected_outer_slab(
    base: &BracketBlock,
    recipe: &IntervalBorePrism,
    z0: i32,
    z1: i32,
) -> Result<BracketBlock, &'static str> {
    let h = i128::from(recipe.height_mm);
    let mut controls = Vec::with_capacity(8);
    for u in 0..2 {
        for v in 0..2 {
            let p = base.controls[(u * 2 + v) * 2].map(i128::from);
            for z in [z0, z1] {
                controls.push([
                    add(
                        mul(h, p[0])?,
                        mul(
                            mul(i128::from(z), i128::from(recipe.top_shift_mm[0]))?,
                            p[3],
                        )?,
                    )?,
                    add(
                        mul(h, p[1])?,
                        mul(
                            mul(i128::from(z), i128::from(recipe.top_shift_mm[1]))?,
                            p[3],
                        )?,
                    )?,
                    mul(mul(h, i128::from(z))?, p[3])?,
                    mul(h, p[3])?,
                ]);
            }
        }
    }
    expected_block(format!("{}_z{z0}_{z1}", base.id), [1, 1, 1], controls)
}

fn planar_bound(block: &BracketBlock) -> Result<f64, &'static str> {
    let [du, dv, dw] = block.degrees;
    if dw != 1 {
        return Err("interval_bore_not_extruded");
    }
    let mut controls = Vec::with_capacity((du + 1) * (dv + 1));
    for u in 0..=du {
        for v in 0..=dv {
            let p = block.controls[(u * (dv + 1) + v) * 2];
            controls.push([p[0], p[1], p[3]]);
        }
    }
    certify_planar_refined([du, dv], &controls).map_err(|error| match error {
        PlanarError::Invalid => "interval_bore_nonpositive_jacobian",
        PlanarError::Indeterminate => "interval_bore_indeterminate_jacobian",
    })
}

fn validate_recipe(recipe: &IntervalBorePrism) -> Result<Vec<i32>, &'static str> {
    if !(1..=4).contains(&recipe.bores.len())
        || !(1..=64).contains(&recipe.height_mm)
        || recipe.outer.vertices.len() + 4 * recipe.bores.len() > 32
        || recipe.top_shift_mm.iter().any(|x| !(-32..=32).contains(x))
    {
        return Err("interval_bore_recipe_range");
    }
    let mut ids = BTreeSet::<String>::new();
    if recipe.outer.id.is_empty() || !ids.insert(recipe.outer.id.clone()) {
        return Err("interval_bore_duplicate_profile_id");
    }
    for v in &recipe.outer.vertices {
        if v.id.is_empty() || !ids.insert(v.id.clone()) {
            return Err("interval_bore_duplicate_profile_id");
        }
    }
    let mut bore_ids = BTreeSet::new();
    let mut events = BTreeSet::from([0, recipe.height_mm]);
    for interval in &recipe.bores {
        let bore = &interval.bore;
        if bore.id.is_empty()
            || bore.id.len() > 18
            || !bore
                .id
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || c == b'_')
            || !bore_ids.insert(bore.id.as_str())
            || !ids.insert(bore.id.clone())
        {
            return Err("interval_bore_invalid_or_duplicate_id");
        }
        if bore.radius_mm < 1 || i64::from(bore.guard_mm) < 2 * i64::from(bore.radius_mm) {
            return Err("interval_bore_guard_radius");
        }
        let [lo, hi] = interval.z_interval_mm;
        if lo < 0 || lo >= hi || hi > recipe.height_mm {
            return Err("interval_bore_z_interval");
        }
        let guard = guard_loop(bore)?;
        if !ids.insert(guard.id.clone()) {
            return Err("interval_bore_duplicate_profile_id");
        }
        for vertex in &guard.vertices {
            if !ids.insert(vertex.id.clone()) {
                return Err("interval_bore_duplicate_profile_id");
            }
        }
        events.insert(lo);
        events.insert(hi);
    }
    Ok(events.into_iter().collect())
}

fn strip_layer_id(id: &str, z0: i32, z1: i32) -> Option<&str> {
    id.strip_suffix(&format!("_z{z0}_{z1}"))
}

/// Audit the exact layered circular-bore template, including cap provenance.
pub fn audit_interval_bore_prism(
    recipe: &IntervalBorePrism,
    blocks: &[BracketBlock],
) -> Result<IntervalBoreCertificate, &'static str> {
    let events = validate_recipe(recipe)?;
    let bores = {
        let mut v = recipe.bores.clone();
        v.sort_by(|a, b| a.bore.id.cmp(&b.bore.id));
        v
    };
    let guards = bores
        .iter()
        .map(|b| guard_loop(&b.bore))
        .collect::<Result<Vec<_>, _>>()?;
    let perforated = PerforatedPrism {
        outer: recipe.outer.clone(),
        holes: guards,
        height_mm: recipe.height_mm,
        top_shift_mm: recipe.top_shift_mm,
    };
    let mut loops = normalized_loops(recipe)?;
    let vertices = loops.drain(..).flat_map(|l| l.vertices).collect::<Vec<_>>();
    let idmap = vertices
        .iter()
        .enumerate()
        .map(|(i, v)| (v.id.as_str(), i))
        .collect::<BTreeMap<_, _>>();
    let mut actual = BTreeMap::<String, &BracketBlock>::new();
    for block in blocks {
        if actual.insert(block.id.clone(), block).is_some() {
            return Err("interval_bore_duplicate_block");
        }
        let _ = primitive_block(block)?;
    }

    let mut layers = Vec::new();
    for pair in events.windows(2) {
        let [z0, z1] = [pair[0], pair[1]];
        if z0 >= z1 {
            return Err("interval_bore_empty_layer");
        }
        layers.push((z0, z1));
    }

    // Recover the common G18 triangulation independently from its exact IDs,
    // then run the existing profile/coverage audit on synthesized full prisms.
    let mut base_outer_ids = None::<BTreeSet<String>>;
    let mut full_outer = BTreeMap::<String, BracketBlock>::new();
    for &(z0, z1) in &layers {
        let layer_ids = actual
            .keys()
            .filter_map(|id| {
                strip_layer_id(id, z0, z1)
                    .filter(|base| base.starts_with("tri/"))
                    .map(|base| (id, base.to_owned()))
            })
            .collect::<Vec<_>>();
        let current = layer_ids
            .iter()
            .map(|(_, base)| base.clone())
            .collect::<BTreeSet<_>>();
        if current.is_empty() || base_outer_ids.as_ref().is_some_and(|old| *old != current) {
            return Err("interval_bore_layer_outer_partition");
        }
        base_outer_ids.get_or_insert(current.clone());
        for (_, base_id) in layer_ids {
            let (tri, at) = parse_triangle_id(&base_id, &idmap)?;
            let controls =
                expected_controls(tri, at, &vertices, recipe.height_mm, recipe.top_shift_mm)
                    .ok_or("interval_bore_outer_expected_control")?;
            let base = BracketBlock {
                id: base_id.clone(),
                degrees: [1, 1, 1],
                controls,
            };
            if let Some(old) = full_outer.get(&base_id) {
                if old != &base {
                    return Err("interval_bore_duplicate_outer_patch");
                }
            } else {
                full_outer.insert(base_id.clone(), base.clone());
            }
            let expected = expected_outer_slab(&base, recipe, z0, z1)?;
            let actual_block = actual
                .get(&expected.id)
                .ok_or("interval_bore_missing_outer_patch")?;
            if actual_block.degrees != expected.degrees
                || primitive_block(actual_block)? != primitive_block(&expected)?
            {
                return Err("interval_bore_altered_outer_patch");
            }
        }
    }
    let outer_certificate = audit_perforated_prism(
        &perforated,
        &full_outer.values().cloned().collect::<Vec<_>>(),
    )?;
    let base_outer_ids = base_outer_ids.ok_or("interval_bore_empty_outer_partition")?;

    let mut expected = BTreeMap::<String, Role>::new();
    let mut min_det = f64::INFINITY;
    for &(z0, z1) in &layers {
        let dz = f64::from(z1 - z0);
        for base_id in &base_outer_ids {
            let base = full_outer
                .get(base_id)
                .ok_or("interval_bore_outer_expected_control")?;
            let block = expected_outer_slab(base, recipe, z0, z1)?;
            let det = outer_certificate.certificate.determinant_lower_bound_mm3 * dz
                / f64::from(recipe.height_mm);
            if !det.is_finite() || det <= 0.0 {
                return Err("interval_bore_indeterminate_jacobian");
            }
            min_det = min_det.min(crate::polygon_prism_audit::next_down(det));
            insert_expected(
                &mut expected,
                block,
                Role::Outer {
                    base_id: base_id.clone(),
                    z0,
                    z1,
                },
            )?;
        }
        for interval in &bores {
            let bore = &interval.bore;
            prove_bore_arcs(bore)?;
            let active = z0 >= interval.z_interval_mm[0] && z1 <= interval.z_interval_mm[1];
            for q in 0..4 {
                for half in 0..2 {
                    let block =
                        expected_sector(bore, recipe, z0, z1, q, half, SectorKind::GuardToCircle)?;
                    let det = next_down(planar_bound(&block)? * dz);
                    if !det.is_finite() || det <= 0.0 {
                        return Err("interval_bore_nonpositive_jacobian");
                    }
                    min_det = min_det.min(det);
                    insert_expected(
                        &mut expected,
                        block,
                        Role::Wall {
                            bore: bore.id.clone(),
                            active,
                            z0,
                            z1,
                        },
                    )?;
                }
            }
            if !active {
                for q in 0..4 {
                    for half in 0..2 {
                        let block = expected_sector(
                            bore,
                            recipe,
                            z0,
                            z1,
                            q,
                            half,
                            SectorKind::CircleToDiamond,
                        )?;
                        let det = next_down(planar_bound(&block)? * dz);
                        if !det.is_finite() || det <= 0.0 {
                            return Err("interval_bore_nonpositive_jacobian");
                        }
                        min_det = min_det.min(det);
                        insert_expected(
                            &mut expected,
                            block,
                            Role::Fill {
                                bore: bore.id.clone(),
                                z0,
                                z1,
                            },
                        )?;
                    }
                }
                for q in 0..4 {
                    let block = expected_core(bore, recipe, z0, z1, q)?;
                    let det = next_down(planar_bound(&block)? * dz);
                    if !det.is_finite() || det <= 0.0 {
                        return Err("interval_bore_nonpositive_jacobian");
                    }
                    min_det = min_det.min(det);
                    insert_expected(
                        &mut expected,
                        block,
                        Role::Core {
                            bore: bore.id.clone(),
                            z0,
                            z1,
                        },
                    )?;
                }
            }
        }
    }
    if expected.len() != blocks.len() || actual.keys().any(|id| !expected.contains_key(id)) {
        return Err("interval_bore_block_set_mismatch");
    }
    let expected_controls =
        expected_controls_by_id(&bores, recipe, &layers, &base_outer_ids, &full_outer)?;
    for (id, block) in &expected_controls {
        let actual_block = actual.get(id).ok_or("interval_bore_missing_block")?;
        if actual_block.degrees != block.degrees
            || primitive_block(actual_block)? != primitive_block(block)?
        {
            return Err("interval_bore_altered_block");
        }
    }

    let incidence: Incidence = enumerate_faces(blocks)?;
    let mut base_tags = BTreeMap::<(String, usize, bool), String>::new();
    for (side, tag) in &outer_certificate.certificate.design_faces {
        base_tags.insert((side.block.clone(), side.axis, side.upper), tag.clone());
    }
    let mut design_faces = Vec::new();
    for side in &incidence.exterior {
        let role = expected
            .get(&side.block)
            .ok_or("interval_bore_boundary_unknown_block")?;
        let tag = match role {
            Role::Outer { base_id, z0, z1 } => {
                if side.axis == 2 {
                    if side.upper && *z1 == recipe.height_mm {
                        "end".to_owned()
                    } else if !side.upper && *z0 == 0 {
                        "start".to_owned()
                    } else {
                        return Err("interval_bore_unexpected_outer_cap");
                    }
                } else {
                    let tag = base_tags
                        .get(&(base_id.clone(), side.axis, side.upper))
                        .ok_or("interval_bore_unexpected_outer_boundary")?;
                    if tag.starts_with("edge/guard_") {
                        return Err("interval_bore_open_guard_boundary");
                    }
                    tag.clone()
                }
            }
            Role::Wall {
                bore,
                active,
                z0,
                z1,
            } => {
                if side.axis == 1 && side.upper {
                    if !active {
                        return Err("interval_bore_exposed_filled_wall");
                    }
                    format!("bore/{bore}/wall")
                } else if side.axis == 2 {
                    if side.upper && *z1 == recipe.height_mm {
                        "end".to_owned()
                    } else if !side.upper && *z0 == 0 {
                        "start".to_owned()
                    } else {
                        return Err("interval_bore_unexpected_wall_cap");
                    }
                } else {
                    return Err("interval_bore_unexpected_wall_boundary");
                }
            }
            Role::Fill { bore, z0, z1 } | Role::Core { bore, z0, z1 } => {
                if side.axis != 2 {
                    return Err("interval_bore_unexpected_fill_boundary");
                }
                let interval = &bores
                    .iter()
                    .find(|b| b.bore.id == *bore)
                    .ok_or("interval_bore_unknown_bore")?
                    .z_interval_mm;
                if side.upper && *z1 == interval[0] && interval[0] > 0 {
                    format!("bore/{bore}/lower")
                } else if !side.upper && *z0 == interval[1] && interval[1] < recipe.height_mm {
                    format!("bore/{bore}/upper")
                } else if !side.upper && *z0 == 0 {
                    "start".to_owned()
                } else if side.upper && *z1 == recipe.height_mm {
                    "end".to_owned()
                } else {
                    return Err("interval_bore_unexpected_fill_cap");
                }
            }
        };
        design_faces.push((side.clone(), tag));
    }
    design_faces.sort_by(|a, b| a.0.cmp(&b.0));
    for bore in &bores {
        let wall_layers = layers
            .iter()
            .filter(|&&(z0, z1)| z0 >= bore.z_interval_mm[0] && z1 <= bore.z_interval_mm[1])
            .count();
        let count = |feature: &str| {
            design_faces
                .iter()
                .filter(|(_, tag)| tag == &format!("bore/{}/{feature}", bore.bore.id))
                .count()
        };
        if count("wall") != 8 * wall_layers
            || count("lower") != if bore.z_interval_mm[0] > 0 { 12 } else { 0 }
            || count("upper")
                != if bore.z_interval_mm[1] < recipe.height_mm {
                    12
                } else {
                    0
                }
        {
            return Err("interval_bore_face_provenance_mismatch");
        }
    }

    let outer_area_twice_mm2 = signed_area2(&recipe.outer.vertices)
        .checked_abs()
        .ok_or("interval_bore_area_overflow")?;
    let removed = bores.iter().try_fold(0_i64, |sum, b| {
        let r = i64::from(b.bore.radius_mm);
        let length = i64::from(b.z_interval_mm[1] - b.z_interval_mm[0]);
        sum.checked_add(
            r.checked_mul(r)
                .and_then(|x| x.checked_mul(length))
                .ok_or("interval_bore_volume_overflow")?,
        )
        .ok_or("interval_bore_volume_overflow")
    })?;
    let enclosed_voids = bores
        .iter()
        .filter(|b| b.z_interval_mm[0] > 0 && b.z_interval_mm[1] < recipe.height_mm)
        .count();
    let through_bores = bores
        .iter()
        .filter(|b| b.z_interval_mm == [0, recipe.height_mm])
        .count();
    let certificate = BracketCertificate {
        incidence,
        determinant_lower_bound_mm3: min_det,
        radial_error_upper_mm: 0.0,
        design_faces,
    };
    Ok(IntervalBoreCertificate {
        certificate,
        outer_area_twice_mm2,
        removed_pi_volume_coefficient_mm3: removed,
        enclosed_voids,
        through_bores,
    })
}

fn insert_expected(
    map: &mut BTreeMap<String, Role>,
    block: BracketBlock,
    role: Role,
) -> Result<(), &'static str> {
    if map.insert(block.id.clone(), role).is_some() {
        return Err("interval_bore_duplicate_expected_block");
    }
    Ok(())
}

fn expected_controls_by_id(
    bores: &[IntervalBore],
    recipe: &IntervalBorePrism,
    layers: &[(i32, i32)],
    base_outer_ids: &BTreeSet<String>,
    full_outer: &BTreeMap<String, BracketBlock>,
) -> Result<BTreeMap<String, BracketBlock>, &'static str> {
    let mut map = BTreeMap::new();
    for &(z0, z1) in layers {
        for bore in bores {
            let active = z0 >= bore.z_interval_mm[0] && z1 <= bore.z_interval_mm[1];
            for q in 0..4 {
                for half in 0..2 {
                    let wall = expected_sector(
                        &bore.bore,
                        recipe,
                        z0,
                        z1,
                        q,
                        half,
                        SectorKind::GuardToCircle,
                    )?;
                    map.insert(wall.id.clone(), wall);
                    if !active {
                        let fill = expected_sector(
                            &bore.bore,
                            recipe,
                            z0,
                            z1,
                            q,
                            half,
                            SectorKind::CircleToDiamond,
                        )?;
                        map.insert(fill.id.clone(), fill);
                    }
                }
            }
            if !active {
                for q in 0..4 {
                    let core = expected_core(&bore.bore, recipe, z0, z1, q)?;
                    map.insert(core.id.clone(), core);
                }
            }
        }
        for base_id in base_outer_ids {
            map.insert(
                format!("{}_z{z0}_{z1}", base_id),
                expected_outer_slab(
                    full_outer
                        .get(base_id)
                        .ok_or("interval_bore_outer_expected_control")?,
                    recipe,
                    z0,
                    z1,
                )?,
            );
        }
    }
    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::interval_bores::build_interval_bore_prism;

    fn recipe(interval: [i32; 2]) -> IntervalBorePrism {
        IntervalBorePrism {
            outer: LoopProfile {
                id: "outer".into(),
                vertices: [[-20, -20], [20, -20], [20, 20], [-20, 20]]
                    .into_iter()
                    .enumerate()
                    .map(|(i, xy_mm)| Vertex {
                        id: format!("v{i}"),
                        xy_mm,
                    })
                    .collect(),
            },
            bores: vec![IntervalBore {
                bore: CircularPolygonBore {
                    id: "pilot".into(),
                    center_mm: [0, 0],
                    radius_mm: 3,
                    guard_mm: 6,
                },
                z_interval_mm: interval,
            }],
            height_mm: 8,
            top_shift_mm: [2, -1],
        }
    }

    #[test]
    fn audit_through_blind_and_closed_voids() {
        for (interval, cavities, through) in [
            ([0, 8], 0, 1),
            ([0, 5], 0, 0),
            ([3, 8], 0, 0),
            ([2, 6], 1, 0),
        ] {
            let r = recipe(interval);
            let blocks = build_interval_bore_prism(&r).unwrap();
            let c = audit_interval_bore_prism(&r, &blocks).unwrap();
            assert_eq!(c.enclosed_voids, cavities);
            assert_eq!(c.through_bores, through);
            assert!(c.certificate.determinant_lower_bound_mm3 > 0.0);
            assert!(c.certificate.radial_error_upper_mm == 0.0);
            let labels = c
                .certificate
                .design_faces
                .iter()
                .map(|x| x.1.as_str())
                .collect::<Vec<_>>();
            assert_eq!(
                labels.iter().filter(|x| **x == "bore/pilot/lower").count(),
                if interval[0] > 0 { 12 } else { 0 }
            );
            assert_eq!(
                labels.iter().filter(|x| **x == "bore/pilot/upper").count(),
                if interval[1] < 8 { 12 } else { 0 }
            );
        }
    }

    #[test]
    fn mutations_refuse_and_common_projective_gauge_is_accepted() {
        let r = recipe([2, 6]);
        let mut blocks = build_interval_bore_prism(&r).unwrap();
        let cert = audit_interval_bore_prism(&r, &blocks).unwrap();
        assert!(cert.removed_pi_volume_coefficient_mm3 > 0);
        for p in &mut blocks[0].controls {
            for x in p {
                *x = x.checked_mul(2).unwrap();
            }
        }
        assert!(audit_interval_bore_prism(&r, &blocks).is_ok());
        blocks[0].controls[0][0] += 1;
        assert!(audit_interval_bore_prism(&r, &blocks).is_err());
        let mut missing = build_interval_bore_prism(&r).unwrap();
        missing.pop();
        assert!(audit_interval_bore_prism(&r, &missing).is_err());
    }
}
