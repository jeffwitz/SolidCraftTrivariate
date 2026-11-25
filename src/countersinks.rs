//! Exact polygon-prism countersinks built by varying certified G21 bore radii.
use crate::{
    bracket::BracketBlock,
    circular_polygon::CircularPolygonBore,
    counterbores::{CounterboreOpening, CounterborePrism},
    interval_bores::{build_interval_bore_prism, IntervalBore, IntervalBorePrism},
    interval_bores_audit::audit_interval_bore_prism,
};
use std::collections::BTreeSet;

#[derive(Clone, Copy)]
struct Layer {
    z0: i32,
    z1: i32,
}

fn add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("countersink_integer_overflow")
}

fn sub(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_sub(b).ok_or("countersink_integer_overflow")
}

fn mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("countersink_integer_overflow")
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
            .ok_or("countersink_control_count_overflow")
    })?;
    if controls.len() != expected || controls.iter().any(|p| p[3] <= 0) {
        return Err("countersink_invalid_homogeneous_net");
    }
    let common = controls
        .iter()
        .flatten()
        .fold(0_u128, |g, x| gcd(g, x.unsigned_abs()));
    if common == 0 {
        return Err("countersink_zero_homogeneous_net");
    }
    let divisor = i128::try_from(common).map_err(|_| "countersink_integer_overflow")?;
    for p in &mut controls {
        for x in p {
            *x /= divisor;
        }
    }
    Ok(BracketBlock {
        id,
        degrees,
        controls: controls
            .into_iter()
            .map(|p| {
                Ok([
                    i64::try_from(p[0]).map_err(|_| "countersink_control_exceeds_i64")?,
                    i64::try_from(p[1]).map_err(|_| "countersink_control_exceeds_i64")?,
                    i64::try_from(p[2]).map_err(|_| "countersink_control_exceeds_i64")?,
                    i64::try_from(p[3]).map_err(|_| "countersink_control_exceeds_i64")?,
                ])
            })
            .collect::<Result<Vec<_>, _>>()?,
    })
}

fn slice_block(
    block: &BracketBlock,
    height: i32,
    layer: Layer,
) -> Result<BracketBlock, &'static str> {
    let [du, dv, dw] = block.degrees;
    if dw != 1 {
        return Err("countersink_source_not_linear_in_z");
    }
    let mut controls = Vec::with_capacity(block.controls.len());
    for u in 0..=du {
        for v in 0..=dv {
            let lower = block.controls[(u * (dv + 1) + v) * 2].map(i128::from);
            let upper = block.controls[(u * (dv + 1) + v) * 2 + 1].map(i128::from);
            for z in [layer.z0, layer.z1] {
                controls.push(std::array::from_fn(|k| {
                    i128::from(height - z) * lower[k] + i128::from(z) * upper[k]
                }));
            }
        }
    }
    let prefix = block
        .id
        .strip_suffix(&format!("_z0_{height}"))
        .ok_or("countersink_unexpected_g21_layer_id")?;
    make_block(
        format!("{prefix}_z{}_{}", layer.z0, layer.z1),
        block.degrees,
        controls,
    )
}

fn radius_numerator(bore: &crate::counterbores::Counterbore, height: i32, z: i32) -> i128 {
    let r = i128::from(bore.radius_mm);
    let big = i128::from(bore.counterbore_radius_mm);
    let depth = i128::from(bore.depth_mm);
    let d = depth;
    let z = i128::from(z);
    match bore.opening {
        CounterboreOpening::Start if z <= d => big * d - (big - r) * z,
        CounterboreOpening::Start => r * d,
        CounterboreOpening::End if z <= i128::from(height) - d => r * d,
        CounterboreOpening::End => r * d + (big - r) * (z - (i128::from(height) - d)),
    }
}

fn reshape_ring(
    block: &mut BracketBlock,
    bore: &crate::counterbores::Counterbore,
    recipe: &CounterborePrism,
    layer: Layer,
) -> Result<(), &'static str> {
    if block.degrees != [2, 1, 1] {
        return Err("countersink_guard_ring_degree_mismatch");
    }
    let h = i128::from(recipe.height_mm);
    let big = i128::from(bore.counterbore_radius_mm);
    let depth = i128::from(bore.depth_mm);
    let scale = mul(mul(h, big)?, depth)?;
    let center = bore.center_mm.map(i128::from);
    let shift = recipe.top_shift_mm.map(i128::from);
    let mut controls = block
        .controls
        .iter()
        .map(|p| p.map(i128::from))
        .collect::<Vec<_>>();
    for u in 0..=2 {
        for v in 0..=1 {
            for w in 0..=1 {
                let index = (u * 2 + v) * 2 + w;
                let p = controls[index];
                let mut q = p
                    .map(|x| mul(scale, x))
                    .into_iter()
                    .collect::<Result<Vec<_>, _>>()?;
                if v == 1 {
                    let z = if w == 0 { layer.z0 } else { layer.z1 };
                    let rho = radius_numerator(bore, recipe.height_mm, z);
                    let complement = sub(big * depth, rho)?;
                    for axis in 0..2 {
                        q[axis] = add(
                            mul(mul(h, complement)?, mul(center[axis], p[3])?)?,
                            add(
                                mul(complement, mul(p[2], shift[axis])?)?,
                                mul(mul(h, rho)?, p[axis])?,
                            )?,
                        )?;
                    }
                }
                controls[index] = [q[0], q[1], q[2], q[3]];
            }
        }
    }
    // Rebuild once to reduce the shared gauge and verify the serialized range.
    *block = make_block(block.id.clone(), block.degrees, controls)?;
    Ok(())
}

/// Build exact linearly tapered bores from certified G21 guard annuli.
pub fn build_countersinks(recipe: &CounterborePrism) -> Result<Vec<BracketBlock>, &'static str> {
    if recipe
        .bores
        .iter()
        .any(|b| b.radius_mm < 1 || b.counterbore_radius_mm <= b.radius_mm)
    {
        return Err("countersink_radii");
    }
    let base_recipe = IntervalBorePrism {
        outer: recipe.outer.clone(),
        bores: recipe
            .bores
            .iter()
            .map(|b| IntervalBore {
                bore: CircularPolygonBore {
                    id: b.id.clone(),
                    center_mm: b.center_mm,
                    radius_mm: b.counterbore_radius_mm,
                    guard_mm: b.guard_mm,
                },
                z_interval_mm: [0, recipe.height_mm],
            })
            .collect(),
        height_mm: recipe.height_mm,
        top_shift_mm: recipe.top_shift_mm,
    };
    let base_blocks = build_interval_bore_prism(&base_recipe)?;
    audit_interval_bore_prism(&base_recipe, &base_blocks)?;

    let mut bores = recipe.bores.clone();
    bores.sort_by(|a, b| a.id.cmp(&b.id));
    let mut events = BTreeSet::from([0, recipe.height_mm]);
    for bore in &bores {
        if !(1..recipe.height_mm).contains(&bore.depth_mm) {
            return Err("countersink_depth");
        }
        events.insert(match bore.opening {
            CounterboreOpening::Start => bore.depth_mm,
            CounterboreOpening::End => recipe.height_mm - bore.depth_mm,
        });
    }
    let events = events.into_iter().collect::<Vec<_>>();
    let mut layers = Vec::new();
    for pair in events.windows(2) {
        layers.push(Layer {
            z0: pair[0],
            z1: pair[1],
        });
    }
    let mut result = Vec::new();
    for layer in layers {
        for base in &base_blocks {
            let mut block = slice_block(base, recipe.height_mm, layer)?;
            if let Some(bore) = bores
                .iter()
                .find(|b| block.id.starts_with(&format!("bore_{}_sector_", b.id)))
            {
                reshape_ring(&mut block, bore, recipe, layer)?;
            }
            result.push(block);
        }
    }
    result.sort_by(|a, b| a.id.cmp(&b.id));
    let mut ids = BTreeSet::new();
    if result.iter().any(|b| !ids.insert(b.id.clone())) {
        return Err("countersink_duplicate_block_id");
    }
    Ok(result)
}
