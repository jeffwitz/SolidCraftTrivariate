//! Independent exact-template audit for polygonal through bores with a blind
//! counterbore at one end.
//!
//! The proof is deliberately limited to the G21 guarded polygon partition and
//! the fixed concentric circular-sector template. G21 certifies the exterior
//! and guard-to-counterbore-radius partition; this module independently
//! checks every replacement annulus, its exact circular boundaries, positive
//! continuous planar Jacobian, incidence and face provenance.

use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    circular_polygon::CircularPolygonBore,
    counterbores::{Counterbore, CounterboreOpening, CounterborePrism},
    incidence::enumerate_faces,
    interval_bores::{build_interval_bore_prism, IntervalBore, IntervalBorePrism},
    interval_bores_audit::audit_interval_bore_prism,
    planar::{certify_planar_refined, PlanarError},
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq)]
pub struct CounterboreCertificate {
    pub certificate: BracketCertificate,
    pub outer_area_twice_mm2: i64,
    pub removed_pi_volume_coefficient_mm3: i64,
    pub enclosed_voids: usize,
    pub through_bores: usize,
}

const QUADRANTS: [&str; 4] = ["east", "north", "west", "south"];

#[derive(Clone)]
enum ExpectedRole {
    Outer,
    CounterRing {
        bore: String,
        active: bool,
        z0: i32,
        z1: i32,
    },
    PilotRing {
        bore: String,
        z0: i32,
        z1: i32,
    },
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn primitive(controls: &[[i128; 4]]) -> Result<Vec<[i128; 4]>, &'static str> {
    let common = controls
        .iter()
        .flatten()
        .fold(0_u128, |g, x| gcd(g, x.unsigned_abs()));
    if common == 0 {
        return Err("counterbore_zero_homogeneous_net");
    }
    let common = i128::try_from(common).map_err(|_| "counterbore_integer_overflow")?;
    Ok(controls.iter().map(|p| p.map(|x| x / common)).collect())
}

fn primitive_block(block: &BracketBlock) -> Result<Vec<[i128; 4]>, &'static str> {
    if block.controls.is_empty() || block.controls.iter().any(|p| p[3] <= 0) {
        return Err("counterbore_nonpositive_weight");
    }
    primitive(
        &block
            .controls
            .iter()
            .map(|p| p.map(i128::from))
            .collect::<Vec<_>>(),
    )
}

fn expected_block(id: String, controls: Vec<[i128; 4]>) -> Result<BracketBlock, &'static str> {
    let controls = primitive(&controls)?;
    Ok(BracketBlock {
        id,
        degrees: [2, 1, 1],
        controls: controls
            .into_iter()
            .map(|p| {
                Ok([
                    i64::try_from(p[0]).map_err(|_| "counterbore_control_overflow")?,
                    i64::try_from(p[1]).map_err(|_| "counterbore_control_overflow")?,
                    i64::try_from(p[2]).map_err(|_| "counterbore_control_overflow")?,
                    i64::try_from(p[3]).map_err(|_| "counterbore_control_overflow")?,
                ])
            })
            .collect::<Result<Vec<_>, _>>()?,
    })
}

fn rotate(point: [i128; 3], quadrant: usize) -> [i128; 3] {
    let [x, y, w] = point;
    match quadrant {
        0 => [x, y, w],
        1 => [-x, y, w],
        2 => [-x, -y, w],
        _ => [x, -y, w],
    }
}

fn source_index(i: usize, quadrant: usize) -> usize {
    if quadrant % 2 == 1 {
        2 - i
    } else {
        i
    }
}

fn circle_profile(radius: i32) -> [[i128; 3]; 3] {
    let r = i128::from(radius);
    [[2 * r, 0, 2], [2 * r, 2 * r, 2], [0, 4 * r, 4]]
}

fn restrict_half(a: [i128; 4], b: [i128; 4], c: [i128; 4], half: usize, i: usize) -> [i128; 4] {
    match (half, i) {
        (0, 0) => a.map(|x| 4 * x),
        (0, 1) => std::array::from_fn(|k| 2 * (a[k] + b[k])),
        (0, 2) | (1, 0) => std::array::from_fn(|k| a[k] + 2 * b[k] + c[k]),
        (1, 1) => std::array::from_fn(|k| 2 * (b[k] + c[k])),
        (1, 2) => c.map(|x| 4 * x),
        _ => unreachable!("validated half-sector index"),
    }
}

/// Reconstruct a pilot-only annulus R→r, including the exact axial shear.
fn expected_pilot_sector(
    bore: &Counterbore,
    recipe: &CounterborePrism,
    z0: i32,
    z1: i32,
    quadrant: usize,
    half: usize,
) -> Result<BracketBlock, &'static str> {
    let [cx, cy] = bore.center_mm.map(i128::from);
    let h = i128::from(recipe.height_mm);
    let profiles = [
        circle_profile(bore.counterbore_radius_mm),
        circle_profile(bore.radius_mm),
    ];
    let mut net = [[[0_i128; 4]; 2]; 3]; // [u][v][homogeneous 4-vector before axial restriction]
    let mut controls = Vec::with_capacity(12);
    for (u, row) in net.iter_mut().enumerate() {
        let source = source_index(u, quadrant);
        for (v, point) in row.iter_mut().enumerate() {
            let [rx, ry, w] = rotate(profiles[v][source], quadrant);
            let x = cx * w + rx;
            let y = cy * w + ry;
            *point = [x, y, w, w];
        }
    }
    // Assemble all homogeneous axial pairs first, then restrict each v/w line
    // in u. The base planar coordinates are unchanged by the axial shear.
    let mut full = [[[[0_i128; 4]; 2]; 2]; 3];
    for (source_row, row) in net.iter().zip(full.iter_mut()) {
        for (source, pair) in source_row.iter().zip(row.iter_mut()) {
            let [x, y, w, _] = *source;
            for (k, z) in [z0, z1].into_iter().enumerate() {
                pair[k] = [
                    h * x + i128::from(recipe.top_shift_mm[0]) * i128::from(z) * w,
                    h * y + i128::from(recipe.top_shift_mm[1]) * i128::from(z) * w,
                    h * i128::from(z) * w,
                    h * w,
                ];
            }
        }
    }
    let mut restricted = [[[[0_i128; 4]; 2]; 2]; 3];
    for (u, row) in restricted.iter_mut().enumerate() {
        for (v, pair) in row.iter_mut().enumerate() {
            for (w, out) in pair.iter_mut().enumerate() {
                *out = restrict_half(full[0][v][w], full[1][v][w], full[2][v][w], half, u);
            }
        }
    }
    for row in restricted {
        for pair in row {
            controls.extend(pair);
        }
    }
    expected_block(
        format!(
            "bore_{}_pilot_sector_{}_half_{half}_z{z0}_{z1}",
            bore.id, QUADRANTS[quadrant]
        ),
        controls,
    )
}

fn planar_lower_bound(block: &BracketBlock) -> Result<f64, &'static str> {
    let [du, dv, dw] = block.degrees;
    if [du, dv, dw] != [2, 1, 1] {
        return Err("counterbore_unexpected_ring_degrees");
    }
    let controls = (0..=du)
        .flat_map(|u| (0..=dv).map(move |v| block.controls[(u * 2 + v) * 2]))
        .map(|p| [p[0], p[1], p[3]])
        .collect::<Vec<_>>();
    certify_planar_refined([du, dv], &controls).map_err(|e| match e {
        PlanarError::Invalid => "counterbore_nonpositive_jacobian",
        PlanarError::Indeterminate => "counterbore_indeterminate_jacobian",
    })
}

fn layer_bounds(recipe: &CounterborePrism) -> BTreeSet<i32> {
    let mut events = BTreeSet::from([0, recipe.height_mm]);
    for bore in &recipe.bores {
        events.insert(match bore.opening {
            CounterboreOpening::Start => bore.depth_mm,
            CounterboreOpening::End => recipe.height_mm - bore.depth_mm,
        });
    }
    events
}

fn g21_comparator(recipe: &CounterborePrism, use_counter_radii: bool) -> IntervalBorePrism {
    let bores = recipe
        .bores
        .iter()
        .map(|bore| {
            let (radius, interval) = if use_counter_radii {
                let interval = match bore.opening {
                    CounterboreOpening::Start => [0, bore.depth_mm],
                    CounterboreOpening::End => [recipe.height_mm - bore.depth_mm, recipe.height_mm],
                };
                (bore.counterbore_radius_mm, interval)
            } else {
                (bore.radius_mm, [0, recipe.height_mm])
            };
            IntervalBore {
                bore: CircularPolygonBore {
                    id: bore.id.clone(),
                    center_mm: bore.center_mm,
                    radius_mm: radius,
                    guard_mm: bore.guard_mm,
                },
                z_interval_mm: interval,
            }
        })
        .collect();
    IntervalBorePrism {
        outer: recipe.outer.clone(),
        bores,
        height_mm: recipe.height_mm,
        top_shift_mm: recipe.top_shift_mm,
    }
}

fn mapped_counter_ring_id(reference_id: &str, bores: &[Counterbore]) -> Option<(String, String)> {
    // IDs may themselves contain `_sector_`; match the longest complete
    // `bore_{id}_sector_` prefix, as the constructor does.
    bores
        .iter()
        .filter_map(|bore| {
            let prefix = format!("bore_{}_sector_", bore.id);
            reference_id.strip_prefix(&prefix).map(|suffix| {
                (
                    bore.id.len(),
                    format!("bore_{}_counter_sector_{suffix}", bore.id),
                    bore.id.clone(),
                )
            })
        })
        .max_by_key(|(length, _, _)| *length)
        .map(|(_, mapped, id)| (mapped, id))
}

/// Independently audit a fixed-template polygonal counterbore model.
pub fn audit_counterbores(
    recipe: &CounterborePrism,
    blocks: &[BracketBlock],
) -> Result<CounterboreCertificate, &'static str> {
    if !(1..=4).contains(&recipe.bores.len()) || recipe.height_mm <= 1 {
        return Err("counterbore_recipe_range");
    }
    let ordered = {
        let mut v = recipe.bores.clone();
        v.sort_by(|a, b| a.id.cmp(&b.id));
        v
    };
    let events = layer_bounds(recipe).into_iter().collect::<Vec<_>>();
    if events.len() < 3 {
        return Err("counterbore_missing_shoulder_event");
    }

    // The G21 counter-radius comparator is itself independently audited. Its
    // whole-disk fills are replaced only on material slabs by the annulus
    // certified below; its outer partition and G→R ring are reused exactly.
    let outer_recipe = g21_comparator(recipe, true);
    let outer_blocks = build_interval_bore_prism(&outer_recipe)?;
    let outer_certificate = audit_interval_bore_prism(&outer_recipe, &outer_blocks)?;
    let pilot_recipe = g21_comparator(recipe, false);
    let pilot_blocks = build_interval_bore_prism(&pilot_recipe)?;
    audit_interval_bore_prism(&pilot_recipe, &pilot_blocks)?;

    let mut actual = BTreeMap::new();
    for block in blocks {
        if actual.insert(block.id.as_str(), block).is_some() {
            return Err("counterbore_duplicate_block");
        }
        let _ = primitive_block(block)?;
    }
    let mut expected_roles = BTreeMap::<String, ExpectedRole>::new();
    let mut min_det = f64::INFINITY;
    let mut expected_count = 0usize;
    for pair in events.windows(2) {
        let [z0, z1] = [pair[0], pair[1]];
        for reference in outer_blocks
            .iter()
            .filter(|block| block.id.ends_with(&format!("_z{z0}_{z1}")))
        {
            let ring = mapped_counter_ring_id(&reference.id, &ordered);
            let is_counter_ring = ring.is_some();
            let (id, role) = if let Some((mapped, bore_id)) = ring {
                let bore = ordered
                    .iter()
                    .find(|b| b.id == bore_id)
                    .ok_or("counterbore_unknown_counter_ring")?;
                let active = match bore.opening {
                    CounterboreOpening::Start => z0 < bore.depth_mm,
                    CounterboreOpening::End => z1 > recipe.height_mm - bore.depth_mm,
                };
                (
                    mapped,
                    ExpectedRole::CounterRing {
                        bore: bore_id,
                        active,
                        z0,
                        z1,
                    },
                )
            } else {
                (reference.id.clone(), ExpectedRole::Outer)
            };
            if is_counter_ring || reference.id.starts_with("tri/") {
                let a = actual
                    .get(id.as_str())
                    .ok_or("counterbore_missing_reused_block")?;
                if a.degrees != reference.degrees
                    || primitive_block(a)? != primitive_block(reference)?
                {
                    return Err("counterbore_altered_reused_g21_block");
                }
                expected_roles.insert(id, role);
                expected_count += 1;
                if is_counter_ring {
                    let det = planar_lower_bound(reference)? * f64::from(z1 - z0);
                    if !det.is_finite() || det <= 0.0 {
                        return Err("counterbore_nonpositive_jacobian");
                    }
                    min_det = min_det.min(det);
                }
            }
        }
        {
            for bore in &ordered {
                let is_pilot_layer = match bore.opening {
                    CounterboreOpening::Start => z0 >= bore.depth_mm,
                    CounterboreOpening::End => z1 <= recipe.height_mm - bore.depth_mm,
                };
                if !is_pilot_layer {
                    continue;
                }
                for q in 0..4 {
                    for half in 0..2 {
                        let expected = expected_pilot_sector(bore, recipe, z0, z1, q, half)?;
                        let actual_block = actual
                            .get(expected.id.as_str())
                            .ok_or("counterbore_missing_pilot_sector")?;
                        if actual_block.degrees != expected.degrees
                            || primitive_block(actual_block)? != primitive_block(&expected)?
                        {
                            return Err("counterbore_altered_pilot_sector");
                        }
                        let det = planar_lower_bound(&expected)? * f64::from(z1 - z0);
                        if !det.is_finite() || det <= 0.0 {
                            return Err("counterbore_nonpositive_jacobian");
                        }
                        min_det = min_det.min(det);
                        expected_roles.insert(
                            expected.id.clone(),
                            ExpectedRole::PilotRing {
                                bore: bore.id.clone(),
                                z0,
                                z1,
                            },
                        );
                        expected_count += 1;
                    }
                }
            }
        }
    }
    if expected_count != blocks.len() || actual.keys().any(|id| !expected_roles.contains_key(*id)) {
        return Err("counterbore_block_inventory_mismatch");
    }

    // The exact R and r circle profiles were certified by the two G21 audits;
    // the independently reconstructed pilot net matches the fixed conic
    // profiles coefficient-for-coefficient after a positive common gauge.
    let incidence = enumerate_faces(blocks)?;
    let mut g21_tags = BTreeMap::<(String, usize, bool), String>::new();
    for (side, label) in &outer_certificate.certificate.design_faces {
        g21_tags.insert((side.block.clone(), side.axis, side.upper), label.clone());
    }
    let mut design_faces = Vec::new();
    for side in &incidence.exterior {
        let role = expected_roles
            .get(&side.block)
            .ok_or("counterbore_unknown_boundary")?;
        let label = match role {
            ExpectedRole::Outer => g21_tags
                .get(&(side.block.clone(), side.axis, side.upper))
                .cloned()
                .ok_or("counterbore_unclassified_outer_face")?,
            ExpectedRole::CounterRing {
                bore,
                active,
                z0,
                z1,
            } => {
                if side.axis == 1 && side.upper {
                    if !active {
                        return Err("counterbore_exposed_internal_radius");
                    }
                    format!("bore/{bore}/counterbore_wall")
                } else if side.axis == 2 && !side.upper && *z0 == 0 {
                    "start".to_owned()
                } else if side.axis == 2 && side.upper && *z1 == recipe.height_mm {
                    "end".to_owned()
                } else {
                    return Err("counterbore_unexpected_counter_ring_boundary");
                }
            }
            ExpectedRole::PilotRing { bore, z0, z1 } => {
                if side.axis == 1 && side.upper {
                    format!("bore/{bore}/pilot_wall")
                } else if side.axis == 2 {
                    let shoulder = ordered
                        .iter()
                        .find(|b| b.id == *bore)
                        .ok_or("counterbore_unknown_bore")?;
                    let z = match shoulder.opening {
                        CounterboreOpening::Start => shoulder.depth_mm,
                        CounterboreOpening::End => recipe.height_mm - shoulder.depth_mm,
                    };
                    if side.upper && *z1 == z || !side.upper && *z0 == z {
                        format!("bore/{bore}/shoulder")
                    } else if !side.upper && *z0 == 0 {
                        "start".to_owned()
                    } else if side.upper && *z1 == recipe.height_mm {
                        "end".to_owned()
                    } else {
                        return Err("counterbore_unexpected_pilot_cap");
                    }
                } else {
                    return Err("counterbore_unexpected_pilot_boundary");
                }
            }
        };
        design_faces.push((side.clone(), label));
    }
    design_faces.sort_by(|a, b| a.0.cmp(&b.0));
    for bore in &ordered {
        let count = |feature: &str| {
            design_faces
                .iter()
                .filter(|(_, tag)| tag == &format!("bore/{}/{feature}", bore.id))
                .count()
        };
        let pilot_layers = events
            .windows(2)
            .filter(|pair| match bore.opening {
                CounterboreOpening::Start => pair[0] >= bore.depth_mm,
                CounterboreOpening::End => pair[1] <= recipe.height_mm - bore.depth_mm,
            })
            .count();
        let counter_layers = events
            .windows(2)
            .filter(|pair| match bore.opening {
                CounterboreOpening::Start => pair[1] <= bore.depth_mm,
                CounterboreOpening::End => pair[0] >= recipe.height_mm - bore.depth_mm,
            })
            .count();
        if count("pilot_wall") != 8 * pilot_layers
            || count("counterbore_wall") != 8 * counter_layers
            || count("shoulder") != 8
        {
            return Err("counterbore_face_provenance_inventory");
        }
    }
    let outer_area_twice_mm2 = outer_certificate.outer_area_twice_mm2;
    let removed = ordered.iter().try_fold(0_i64, |sum, bore| {
        let r = i64::from(bore.radius_mm);
        let big = i64::from(bore.counterbore_radius_mm);
        let h = i64::from(recipe.height_mm);
        let d = i64::from(bore.depth_mm);
        let r2 = r.checked_mul(r).ok_or("counterbore_volume_overflow")?;
        let big2 = big.checked_mul(big).ok_or("counterbore_volume_overflow")?;
        let value = r2
            .checked_mul(h)
            .and_then(|x| {
                big2.checked_sub(r2)
                    .and_then(|z| z.checked_mul(d))
                    .and_then(|y| x.checked_add(y))
            })
            .ok_or("counterbore_volume_overflow")?;
        sum.checked_add(value).ok_or("counterbore_volume_overflow")
    })?;
    min_det = min_det.min(outer_certificate.certificate.determinant_lower_bound_mm3);
    let certificate = BracketCertificate {
        incidence,
        determinant_lower_bound_mm3: min_det,
        radial_error_upper_mm: 0.0,
        design_faces,
    };
    Ok(CounterboreCertificate {
        certificate,
        outer_area_twice_mm2,
        removed_pi_volume_coefficient_mm3: removed,
        enclosed_voids: 0,
        through_bores: ordered.len(),
    })
}
