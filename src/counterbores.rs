//! Exact polygon-prism counterbores built from the certified G21 guard complex.
use crate::{
    bracket::BracketBlock,
    circular_polygon::CircularPolygonBore,
    interval_bores::{build_interval_bore_prism, IntervalBore, IntervalBorePrism},
    interval_bores_audit::audit_interval_bore_prism,
    perforated_prism::LoopProfile,
};
use std::collections::BTreeSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CounterboreOpening {
    Start,
    End,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Counterbore {
    pub id: String,
    pub center_mm: [i32; 2],
    pub radius_mm: i32,
    pub guard_mm: i32,
    pub counterbore_radius_mm: i32,
    pub depth_mm: i32,
    pub opening: CounterboreOpening,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CounterborePrism {
    pub outer: LoopProfile,
    pub bores: Vec<Counterbore>,
    pub height_mm: i32,
    pub top_shift_mm: [i32; 2],
}

#[derive(Clone, Copy)]
struct Layer {
    z0: i32,
    z1: i32,
}

fn add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("counterbore_integer_overflow")
}

fn mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("counterbore_integer_overflow")
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn make_block(
    id: String,
    degrees: [usize; 3],
    mut controls: Vec<[i128; 4]>,
) -> Result<BracketBlock, &'static str> {
    let expected = degrees.into_iter().try_fold(1_usize, |count, degree| {
        count
            .checked_mul(degree + 1)
            .ok_or("counterbore_control_count_overflow")
    })?;
    if controls.len() != expected {
        return Err("counterbore_control_count_mismatch");
    }
    if controls.iter().any(|control| control[3] <= 0) {
        return Err("counterbore_nonpositive_weight");
    }
    let common = controls
        .iter()
        .flatten()
        .fold(0_u128, |current, value| gcd(current, value.unsigned_abs()));
    if common == 0 {
        return Err("counterbore_zero_homogeneous_net");
    }
    let divisor = i128::try_from(common).map_err(|_| "counterbore_integer_overflow")?;
    for control in &mut controls {
        for coordinate in control {
            *coordinate /= divisor;
        }
    }
    Ok(BracketBlock {
        id,
        degrees,
        controls: controls
            .into_iter()
            .map(|control| {
                Ok([
                    i64::try_from(control[0]).map_err(|_| "counterbore_control_exceeds_i64")?,
                    i64::try_from(control[1]).map_err(|_| "counterbore_control_exceeds_i64")?,
                    i64::try_from(control[2]).map_err(|_| "counterbore_control_exceeds_i64")?,
                    i64::try_from(control[3]).map_err(|_| "counterbore_control_exceeds_i64")?,
                ])
            })
            .collect::<Result<Vec<_>, _>>()?,
    })
}

fn split_block(
    block: &BracketBlock,
    height: i32,
    layer: Layer,
) -> Result<BracketBlock, &'static str> {
    let [du, dv, dw] = block.degrees;
    if dw != 1 {
        return Err("counterbore_source_block_must_be_linear_in_z");
    }
    let nv = dv + 1;
    let mut controls = Vec::with_capacity(block.controls.len());
    for u in 0..=du {
        for v in 0..=dv {
            let p0 = block.controls[(u * nv + v) * 2].map(i128::from);
            let p1 = block.controls[(u * nv + v) * 2 + 1].map(i128::from);
            for z in [layer.z0, layer.z1] {
                let mut point = [0_i128; 4];
                for k in 0..4 {
                    point[k] = add(
                        mul(i128::from(height - z), p0[k])?,
                        mul(i128::from(z), p1[k])?,
                    )?;
                }
                controls.push(point);
            }
        }
    }
    let suffix = format!("_z0_{height}");
    let prefix = block
        .id
        .strip_suffix(&suffix)
        .ok_or("counterbore_g21_block_has_unexpected_layer_id")?;
    make_block(
        format!("{prefix}_z{}_{}", layer.z0, layer.z1),
        block.degrees,
        controls,
    )
}

fn source_index(i: usize, quadrant: usize) -> usize {
    if quadrant % 2 == 1 {
        2 - i
    } else {
        i
    }
}

fn rotate_profile(point: [i128; 3], quadrant: usize) -> [i128; 3] {
    let [x, y, weight] = point;
    match quadrant {
        0 => [x, y, weight],
        1 => [-x, y, weight],
        2 => [-x, -y, weight],
        _ => [x, -y, weight],
    }
}

fn restrict_quadratic(
    p0: [i128; 4],
    p1: [i128; 4],
    p2: [i128; 4],
    half: usize,
    i: usize,
) -> [i128; 4] {
    match (half, i) {
        (0, 0) => p0.map(|x| 4 * x),
        (0, 1) => std::array::from_fn(|k| 2 * (p0[k] + p1[k])),
        (0, 2) | (1, 0) => std::array::from_fn(|k| p0[k] + 2 * p1[k] + p2[k]),
        (1, 1) => std::array::from_fn(|k| 2 * (p1[k] + p2[k])),
        (1, 2) => p2.map(|x| 4 * x),
        _ => unreachable!("validated quadratic restriction index"),
    }
}

fn concentric_sector(
    recipe: &CounterborePrism,
    bore: &Counterbore,
    layer: Layer,
    quadrant: usize,
    half: usize,
) -> Result<BracketBlock, &'static str> {
    let h = i128::from(recipe.height_mm);
    let center = bore.center_mm.map(i128::from);
    let shift = recipe.top_shift_mm.map(i128::from);
    let radii = [bore.counterbore_radius_mm, bore.radius_mm].map(i128::from);
    // The endpoint weights match the G21 quarter-circle chart exactly.
    let arc = |radius: i128| {
        [
            [2 * radius, 0, 2],
            [2 * radius, 2 * radius, 2],
            [0, 4 * radius, 4],
        ]
    };
    let outer = arc(radii[0]);
    let inner = arc(radii[1]);
    let mut net = [[[[0_i128; 4]; 2]; 2]; 3];
    for (u, row) in net.iter_mut().enumerate() {
        let source = source_index(u, quadrant);
        for (v, pair) in row.iter_mut().enumerate() {
            let point = if v == 0 { outer[source] } else { inner[source] };
            let [x, y, weight] = rotate_profile(point, quadrant);
            for (k, z) in [layer.z0, layer.z1].into_iter().enumerate() {
                let zw = mul(i128::from(z), weight)?;
                pair[k] = [
                    add(
                        mul(h, add(mul(center[0], weight)?, x)?)?,
                        mul(zw, shift[0])?,
                    )?,
                    add(
                        mul(h, add(mul(center[1], weight)?, y)?)?,
                        mul(zw, shift[1])?,
                    )?,
                    mul(mul(h, i128::from(z))?, weight)?,
                    mul(h, weight)?,
                ];
            }
        }
    }
    let mut controls = Vec::with_capacity(12);
    for u in 0..3 {
        for ((v0, v1), v2) in net[0].iter().zip(&net[1]).zip(&net[2]) {
            for ((&p0, &p1), &p2) in v0.iter().zip(v1).zip(v2) {
                controls.push(restrict_quadratic(p0, p1, p2, half, u));
            }
        }
    }
    let quadrant_name = ["east", "north", "west", "south"][quadrant];
    make_block(
        format!(
            "bore_{}_pilot_sector_{quadrant_name}_half_{half}_z{}_{}",
            bore.id, layer.z0, layer.z1
        ),
        [2, 1, 1],
        controls,
    )
}

fn shoulder_z(bore: &Counterbore, height: i32) -> i32 {
    match bore.opening {
        CounterboreOpening::Start => bore.depth_mm,
        CounterboreOpening::End => height - bore.depth_mm,
    }
}

fn pilot_only(bore: &Counterbore, height: i32, layer: Layer) -> bool {
    match bore.opening {
        CounterboreOpening::Start => layer.z0 >= bore.depth_mm,
        CounterboreOpening::End => layer.z1 <= height - bore.depth_mm,
    }
}

fn validate(recipe: &CounterborePrism) -> Result<Vec<Counterbore>, &'static str> {
    if !(1..=4).contains(&recipe.bores.len()) {
        return Err("counterbore_count");
    }
    let mut bores = recipe.bores.clone();
    bores.sort_by(|a, b| a.id.cmp(&b.id));
    let mut ids = BTreeSet::new();
    for bore in &bores {
        if bore.id.is_empty()
            || bore.id.len() > 18
            || !bore
                .id
                .bytes()
                .all(|byte| byte.is_ascii_alphanumeric() || byte == b'_')
            || !ids.insert(bore.id.clone())
        {
            return Err("counterbore_invalid_or_duplicate_id");
        }
        if bore.radius_mm < 1
            || bore.counterbore_radius_mm <= bore.radius_mm
            || i64::from(bore.guard_mm) < 2 * i64::from(bore.counterbore_radius_mm)
        {
            return Err("counterbore_radii_or_guard");
        }
        if !(1..recipe.height_mm).contains(&bore.depth_mm) {
            return Err("counterbore_depth");
        }
    }
    Ok(bores)
}

/// Build exact start- or end-opened counterbores using G21's certified guard maps.
pub fn build_counterbores(recipe: &CounterborePrism) -> Result<Vec<BracketBlock>, &'static str> {
    let bores = validate(recipe)?;
    let base_recipe = IntervalBorePrism {
        outer: recipe.outer.clone(),
        bores: bores
            .iter()
            .map(|bore| IntervalBore {
                bore: CircularPolygonBore {
                    id: bore.id.clone(),
                    center_mm: bore.center_mm,
                    radius_mm: bore.counterbore_radius_mm,
                    guard_mm: bore.guard_mm,
                },
                z_interval_mm: [0, recipe.height_mm],
            })
            .collect(),
        height_mm: recipe.height_mm,
        top_shift_mm: recipe.top_shift_mm,
    };
    let base_blocks = build_interval_bore_prism(&base_recipe)?;
    audit_interval_bore_prism(&base_recipe, &base_blocks)?;

    let mut events = BTreeSet::from([0, recipe.height_mm]);
    for bore in &bores {
        events.insert(shoulder_z(bore, recipe.height_mm));
    }
    let events = events.into_iter().collect::<Vec<_>>();
    let layers = events
        .windows(2)
        .map(|pair| Layer {
            z0: pair[0],
            z1: pair[1],
        })
        .collect::<Vec<_>>();

    let mut id_order = bores.clone();
    id_order.sort_by_key(|bore| std::cmp::Reverse(bore.id.len()));
    let id_prefixes = id_order
        .iter()
        .map(|bore| {
            (
                format!("bore_{}_sector_", bore.id),
                format!("bore_{}_counter_sector_", bore.id),
            )
        })
        .collect::<Vec<_>>();

    let mut blocks = Vec::new();
    for layer in &layers {
        for block in &base_blocks {
            let mut prefix = block
                .id
                .strip_suffix(&format!("_z0_{}", recipe.height_mm))
                .ok_or("counterbore_g21_block_has_unexpected_layer_id")?
                .to_owned();
            for (source, target) in &id_prefixes {
                if prefix.starts_with(source) {
                    prefix = prefix.replacen(source, target, 1);
                    break;
                }
            }
            let renamed = BracketBlock {
                id: format!("{prefix}_z0_{}", recipe.height_mm),
                degrees: block.degrees,
                controls: block.controls.clone(),
            };
            blocks.push(split_block(&renamed, recipe.height_mm, *layer)?);
        }
        for bore in &bores {
            if pilot_only(bore, recipe.height_mm, *layer) {
                for quadrant in 0..4 {
                    for half in 0..2 {
                        blocks.push(concentric_sector(recipe, bore, *layer, quadrant, half)?);
                    }
                }
            }
        }
    }
    blocks.sort_by(|a, b| a.id.cmp(&b.id));
    let mut ids = BTreeSet::new();
    if blocks.iter().any(|block| !ids.insert(block.id.clone())) {
        return Err("counterbore_duplicate_block_id");
    }
    Ok(blocks)
}
