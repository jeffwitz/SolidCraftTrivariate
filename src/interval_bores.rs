//! Finite-depth circular excavations in exact sheared polygonal prisms.
use crate::{
    bracket::BracketBlock,
    circular_polygon::CircularPolygonBore,
    perforated_prism::{build_perforated_prism, LoopProfile, PerforatedPrism},
    polygon_prism::Vertex,
};
use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntervalBore {
    pub bore: CircularPolygonBore,
    pub z_interval_mm: [i32; 2],
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntervalBorePrism {
    pub outer: LoopProfile,
    pub bores: Vec<IntervalBore>,
    pub height_mm: i32,
    pub top_shift_mm: [i32; 2],
}

fn add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("interval_bore_integer_overflow")
}

fn mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("interval_bore_integer_overflow")
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn block_from_homogeneous(
    id: String,
    degrees: [usize; 3],
    mut controls: Vec<[i128; 4]>,
) -> Result<BracketBlock, &'static str> {
    let expected = degrees.into_iter().try_fold(1_usize, |n, degree| {
        n.checked_mul(degree + 1)
            .ok_or("interval_bore_control_count_overflow")
    })?;
    if controls.len() != expected {
        return Err("interval_bore_control_count_mismatch");
    }
    if controls.iter().any(|p| p[3] <= 0) {
        return Err("interval_bore_nonpositive_weight");
    }
    let common = controls
        .iter()
        .flatten()
        .fold(0_u128, |g, x| gcd(g, x.unsigned_abs()));
    if common == 0 {
        return Err("interval_bore_zero_homogeneous_net");
    }
    for p in &mut controls {
        for x in p {
            *x /= i128::try_from(common).map_err(|_| "interval_bore_integer_overflow")?;
        }
    }
    let controls = controls
        .into_iter()
        .map(|p| {
            Ok([
                i64::try_from(p[0]).map_err(|_| "interval_bore_control_exceeds_i64")?,
                i64::try_from(p[1]).map_err(|_| "interval_bore_control_exceeds_i64")?,
                i64::try_from(p[2]).map_err(|_| "interval_bore_control_exceeds_i64")?,
                i64::try_from(p[3]).map_err(|_| "interval_bore_control_exceeds_i64")?,
            ])
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(BracketBlock {
        id,
        degrees,
        controls,
    })
}

fn layer_value(
    lower: [i64; 4],
    upper: [i64; 4],
    height: i32,
    z: i32,
) -> Result<[i128; 4], &'static str> {
    let h = i128::from(height);
    let z = i128::from(z);
    Ok(std::array::from_fn(|k| {
        (h - z) * i128::from(lower[k]) + z * i128::from(upper[k])
    }))
}

fn slice_prism_block(
    block: &BracketBlock,
    height: i32,
    z0: i32,
    z1: i32,
) -> Result<BracketBlock, &'static str> {
    let [du, dv, dw] = block.degrees;
    if dw != 1 {
        return Err("interval_bore_source_must_be_linearly_extruded");
    }
    let nv = dv + 1;
    let nw = 2;
    let mut controls = Vec::with_capacity(block.controls.len());
    for u in 0..=du {
        for v in 0..=dv {
            let lower = block.controls[(u * nv + v) * nw];
            let upper = block.controls[(u * nv + v) * nw + 1];
            controls.push(layer_value(lower, upper, height, z0)?);
            controls.push(layer_value(lower, upper, height, z1)?);
        }
    }
    block_from_homogeneous(
        format!("{}_z{}_{}", block.id, z0, z1),
        block.degrees,
        controls,
    )
}

fn guard_loop(bore: &CircularPolygonBore) -> Result<LoopProfile, &'static str> {
    let [cx, cy] = bore.center_mm.map(i64::from);
    let g = i64::from(bore.guard_mm);
    let points = [[cx + g, cy], [cx, cy + g], [cx - g, cy], [cx, cy - g]];
    if points.iter().flatten().any(|&v| !(-64..=64).contains(&v)) {
        return Err("interval_bore_guard_outside_domain");
    }
    Ok(LoopProfile {
        id: format!("guard_{}", bore.id),
        vertices: points
            .into_iter()
            .zip(["e", "n", "w", "s"])
            .map(|(xy, direction)| Vertex {
                id: format!("g_{}_{}", bore.id, direction),
                xy_mm: [xy[0] as i32, xy[1] as i32],
            })
            .collect(),
    })
}

fn rotate_profile(point: [i128; 3], quadrant: usize) -> [i128; 3] {
    let [x, y, w] = point;
    match quadrant {
        0 => [x, y, w],
        // The quarter-conic has endpoint weights 2,2,4. Odd sectors use a
        // paired parameter reversal and reflection so adjacent radial traces
        // retain the same projective parameterization while angle still grows
        // counter-clockwise.
        1 => [-x, y, w],
        2 => [-x, -y, w],
        _ => [x, -y, w],
    }
}

fn rotate_core_point(point: [i128; 3], quadrant: usize) -> [i128; 3] {
    let [x, y, w] = point;
    match quadrant {
        0 => [x, y, w],
        1 => [-y, x, w],
        2 => [-x, -y, w],
        _ => [y, -x, w],
    }
}

fn source_index(i: usize, quadrant: usize) -> usize {
    if quadrant % 2 == 1 {
        2 - i
    } else {
        i
    }
}

fn restricted_control(
    b0: [i128; 4],
    b1: [i128; 4],
    b2: [i128; 4],
    half: usize,
    i: usize,
) -> [i128; 4] {
    match (half, i) {
        (0, 0) => b0.map(|x| 4 * x),
        (0, 1) => std::array::from_fn(|k| 2 * (b0[k] + b1[k])),
        (0, 2) | (1, 0) => std::array::from_fn(|k| b0[k] + 2 * b1[k] + b2[k]),
        (1, 1) => std::array::from_fn(|k| 2 * (b1[k] + b2[k])),
        (1, 2) => b2.map(|x| 4 * x),
        _ => unreachable!("validated half sector and quadratic index"),
    }
}

#[derive(Clone, Copy)]
enum RingKind {
    GuardToCircle,
    CircleToInnerDiamond,
}

#[derive(Clone, Copy)]
struct Slab {
    height_mm: i32,
    shift_mm: [i32; 2],
    z0: i32,
    z1: i32,
}

fn ring_profiles(bore: &CircularPolygonBore, kind: RingKind) -> [[[i128; 3]; 2]; 3] {
    let r = i128::from(bore.radius_mm);
    let g = i128::from(bore.guard_mm);
    let guard = [[2 * g, 0, 2], [g, g, 2], [0, 2 * g, 2]];
    let circle = [[2 * r, 0, 2], [2 * r, 2 * r, 2], [0, 4 * r, 4]];
    // The inner diamond has half-diagonal r/2. Its middle control is the
    // exact degree-two elevation of a straight edge at common homogeneous scale.
    // Degree-elevate the straight inner-diamond edge at one common weight.
    // This keeps its trace parameterization identical to the Q1 core edge.
    let diamond = [[2 * r, 0, 4], [r, r, 4], [0, 2 * r, 4]];
    std::array::from_fn(|u| {
        let j = source_index(u, 0);
        match kind {
            RingKind::GuardToCircle => [guard[j], circle[j]],
            RingKind::CircleToInnerDiamond => [circle[j], diamond[j]],
        }
    })
}

fn ring_block(
    bore: &CircularPolygonBore,
    slab: Slab,
    quadrant: usize,
    half: usize,
    kind: RingKind,
) -> Result<BracketBlock, &'static str> {
    let profiles = ring_profiles(bore, kind);
    let center = bore.center_mm.map(i128::from);
    let h = i128::from(slab.height_mm);
    let mut layers = [[[[0_i128; 4]; 2]; 2]; 3]; // u, v, layer
    for (u, row) in layers.iter_mut().enumerate() {
        let source_u = source_index(u, quadrant);
        for (v, pair) in row.iter_mut().enumerate() {
            let [rx, ry, weight] = rotate_profile(profiles[source_u][v], quadrant);
            for (layer, z) in [slab.z0, slab.z1].into_iter().enumerate() {
                let zw = i128::from(z) * weight;
                pair[layer] = [
                    add(
                        mul(h, add(mul(center[0], weight)?, rx)?)?,
                        mul(zw, i128::from(slab.shift_mm[0]))?,
                    )?,
                    add(
                        mul(h, add(mul(center[1], weight)?, ry)?)?,
                        mul(zw, i128::from(slab.shift_mm[1]))?,
                    )?,
                    mul(mul(h, i128::from(z))?, weight)?,
                    mul(h, weight)?,
                ];
            }
        }
    }
    let mut controls = Vec::with_capacity(12);
    for u in 0..3 {
        for ((v0, v1), v2) in layers[0].iter().zip(&layers[1]).zip(&layers[2]) {
            for ((&b0, &b1), &b2) in v0.iter().zip(v1).zip(v2) {
                controls.push(restricted_control(b0, b1, b2, half, u));
            }
        }
    }
    let quadrant_name = ["east", "north", "west", "south"][quadrant];
    let prefix = match kind {
        RingKind::GuardToCircle => format!("bore_{}_sector", bore.id),
        RingKind::CircleToInnerDiamond => format!("bore_{}_fill_sector", bore.id),
    };
    block_from_homogeneous(
        format!(
            "{prefix}_{quadrant_name}_half_{half}_z{}_{}",
            slab.z0, slab.z1
        ),
        [2, 1, 1],
        controls,
    )
}

fn core_block(
    bore: &CircularPolygonBore,
    slab: Slab,
    quadrant: usize,
) -> Result<BracketBlock, &'static str> {
    let r = i128::from(bore.radius_mm);
    // Controls have weight four; the diamond's half-diagonal is r/2.
    let east = [2 * r, 0, 4];
    let north_east = [r, r, 4];
    let origin = [0, 0, 4];
    let south_east = [r, -r, 4];
    // Flattened tensor order is [P00, P01, P10, P11] with v varying first.
    // This gives positive orientation from the east diamond cell toward its
    // south-east and north-east neighbors; proper CCW rotations place the
    // remaining three cells without reflecting either parameter axis.
    let base = [east, south_east, north_east, origin];
    let center = bore.center_mm.map(i128::from);
    let h = i128::from(slab.height_mm);
    let mut controls = Vec::with_capacity(8);
    for point in base {
        let [x, y, weight] = rotate_core_point(point, quadrant);
        for z in [slab.z0, slab.z1] {
            controls.push([
                add(
                    mul(h, add(mul(center[0], weight)?, x)?)?,
                    mul(mul(i128::from(z), weight)?, i128::from(slab.shift_mm[0]))?,
                )?,
                add(
                    mul(h, add(mul(center[1], weight)?, y)?)?,
                    mul(mul(i128::from(z), weight)?, i128::from(slab.shift_mm[1]))?,
                )?,
                mul(mul(h, i128::from(z))?, weight)?,
                mul(h, weight)?,
            ]);
        }
    }
    let name = ["east", "north", "west", "south"][quadrant];
    block_from_homogeneous(
        format!(
            "bore_{}_core_sector_{name}_z{}_{}",
            bore.id, slab.z0, slab.z1
        ),
        [1, 1, 1],
        controls,
    )
}

pub(crate) fn validate(recipe: &IntervalBorePrism) -> Result<Vec<IntervalBore>, &'static str> {
    if !(1..=4).contains(&recipe.bores.len())
        || !(1..=64).contains(&recipe.height_mm)
        || recipe.outer.vertices.len() + 4 * recipe.bores.len() > 32
        || recipe.top_shift_mm.iter().any(|x| !(-32..=32).contains(x))
    {
        return Err("interval_bore_recipe_range");
    }
    let mut ids = BTreeSet::new();
    if recipe.outer.id.is_empty() || !ids.insert(recipe.outer.id.clone()) {
        return Err("interval_bore_invalid_outer_id");
    }
    for vertex in &recipe.outer.vertices {
        if vertex.id.is_empty() || !ids.insert(vertex.id.clone()) {
            return Err("interval_bore_duplicate_profile_id");
        }
    }
    let mut bores = recipe.bores.clone();
    bores.sort_by(|a, b| a.bore.id.cmp(&b.bore.id));
    for interval in &bores {
        let bore = &interval.bore;
        if bore.id.is_empty()
            || bore.id.len() > 18
            || !bore
                .id
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || c == b'_')
            || !ids.insert(bore.id.clone())
        {
            return Err("interval_bore_invalid_or_duplicate_id");
        }
        if bore.radius_mm < 1 || i64::from(bore.guard_mm) < 2 * i64::from(bore.radius_mm) {
            return Err("interval_bore_guard_radius");
        }
        if interval.z_interval_mm[0] < 0
            || interval.z_interval_mm[0] >= interval.z_interval_mm[1]
            || interval.z_interval_mm[1] > recipe.height_mm
        {
            return Err("interval_bore_z_interval");
        }
        let guard = guard_loop(bore)?;
        if !ids.insert(guard.id) {
            return Err("interval_bore_duplicate_profile_id");
        }
        for vertex in guard.vertices {
            if !ids.insert(vertex.id) {
                return Err("interval_bore_duplicate_profile_id");
            }
        }
    }
    Ok(bores)
}

/// Build a G21 partition with exact axial events, guard annuli and filled disk
/// sectors only on material slabs. This function constructs a candidate; a
/// caller must run the independent interval-bore audit before adoption.
pub fn build_interval_bore_prism(
    recipe: &IntervalBorePrism,
) -> Result<Vec<BracketBlock>, &'static str> {
    let bores = validate(recipe)?;
    let base = prepare_profile(recipe, &bores)?;
    let mut blocks = Vec::new();
    for [z0, z1] in layer_intervals(recipe.height_mm, &bores) {
        blocks.extend(build_layer(&base, recipe, &bores, z0, z1)?);
    }
    sort_partition(&mut blocks)?;
    Ok(blocks)
}

pub(crate) fn layer_intervals(height: i32, bores: &[IntervalBore]) -> Vec<[i32; 2]> {
    let mut events = BTreeSet::from([0, height]);
    for bore in bores {
        events.insert(bore.z_interval_mm[0]);
        events.insert(bore.z_interval_mm[1]);
    }
    let events = events.into_iter().collect::<Vec<_>>();
    events.windows(2).map(|z| [z[0], z[1]]).collect()
}

pub(crate) fn prepare_profile(
    recipe: &IntervalBorePrism,
    bores: &[IntervalBore],
) -> Result<Vec<BracketBlock>, &'static str> {
    let guards = bores
        .iter()
        .map(|bore| guard_loop(&bore.bore))
        .collect::<Result<Vec<_>, _>>()?;
    build_perforated_prism(&PerforatedPrism {
        outer: recipe.outer.clone(),
        holes: guards,
        height_mm: recipe.height_mm,
        top_shift_mm: recipe.top_shift_mm,
    })
}

pub(crate) fn build_layer(
    base: &[BracketBlock],
    recipe: &IntervalBorePrism,
    bores: &[IntervalBore],
    z0: i32,
    z1: i32,
) -> Result<Vec<BracketBlock>, &'static str> {
    let mut blocks = Vec::new();
    let slab = Slab {
        height_mm: recipe.height_mm,
        shift_mm: recipe.top_shift_mm,
        z0,
        z1,
    };
    for block in base {
        blocks.push(slice_prism_block(block, recipe.height_mm, z0, z1)?);
    }
    for interval in bores {
        let bore = &interval.bore;
        for quadrant in 0..4 {
            for half in 0..2 {
                blocks.push(ring_block(
                    bore,
                    slab,
                    quadrant,
                    half,
                    RingKind::GuardToCircle,
                )?);
            }
        }
        let [lo, hi] = interval.z_interval_mm;
        if z1 <= lo || z0 >= hi {
            for quadrant in 0..4 {
                for half in 0..2 {
                    blocks.push(ring_block(
                        bore,
                        slab,
                        quadrant,
                        half,
                        RingKind::CircleToInnerDiamond,
                    )?);
                }
            }
            for quadrant in 0..4 {
                blocks.push(core_block(bore, slab, quadrant)?);
            }
        }
    }
    Ok(blocks)
}

pub(crate) fn sort_partition(blocks: &mut [BracketBlock]) -> Result<(), &'static str> {
    blocks.sort_by(|a, b| a.id.cmp(&b.id));
    let mut ids = BTreeSet::new();
    if blocks.iter().any(|block| !ids.insert(block.id.clone())) {
        return Err("interval_bore_duplicate_block_id");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn finite_interval_adds_caps_only_on_material_sides() {
        let blocks = build_interval_bore_prism(&recipe([2, 6])).unwrap();
        let ids = blocks
            .iter()
            .map(|b| b.id.as_str())
            .collect::<BTreeSet<_>>();
        assert!(ids.contains("bore_pilot_sector_east_half_0_z2_6"));
        assert!(ids.contains("bore_pilot_fill_sector_east_half_0_z0_2"));
        assert!(ids.contains("bore_pilot_fill_sector_east_half_0_z6_8"));
        assert!(ids.contains("bore_pilot_core_sector_east_z0_2"));
        assert!(!ids.contains("bore_pilot_core_sector_east_z2_6"));
        assert!(blocks.iter().all(|b| b.controls.iter().all(|p| p[3] > 0)));
    }

    #[test]
    fn through_interval_has_no_disk_filler_and_invalid_intervals_refuse() {
        let blocks = build_interval_bore_prism(&recipe([0, 8])).unwrap();
        assert!(blocks
            .iter()
            .all(|b| !b.id.contains("fill_sector") && !b.id.contains("core_sector")));
        for interval in [[0, 0], [6, 2], [-1, 4], [4, 9]] {
            assert!(build_interval_bore_prism(&recipe(interval)).is_err());
        }
    }

    #[test]
    fn bore_order_does_not_change_the_partition() {
        let mut first = recipe([1, 7]);
        first.bores.push(IntervalBore {
            bore: CircularPolygonBore {
                id: "second".into(),
                center_mm: [10, 10],
                radius_mm: 2,
                guard_mm: 4,
            },
            z_interval_mm: [3, 6],
        });
        let mut second = first.clone();
        second.bores.reverse();
        assert_eq!(
            build_interval_bore_prism(&first).unwrap(),
            build_interval_bore_prism(&second).unwrap()
        );
    }
}
