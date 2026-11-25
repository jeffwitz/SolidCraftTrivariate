//! Independent G2-transfer audit for variable-depth square pockets.
//!
//! Every accepted layer is mapped exactly back to the frozen G2 section
//! templates. Rational axial maps are evaluated in i128 and a single positive
//! homogeneous scale is removed per block before conversion back to i64.

use super::bracket::BracketBlock;
use super::bracket_audit::BracketCertificate;
use super::incidence::{enumerate_faces, BlockSide};
use super::section_audit::audit_section;
use std::collections::{BTreeMap, BTreeSet};
use std::convert::TryFrom;

const SECTORS: [&str; 4] = ["east", "north", "west", "south"];

#[derive(Clone, Copy, PartialEq, Eq)]
enum Regime {
    Solid,
    Blind,
    Through,
}

/// Certify one pocket model at an integer half-millimetre depth.
///
/// `depth_ticks` is in `0..=18`. The returned certificate covers local
/// regularity, exact interfaces, exterior faces, and the G2-qualified domain
/// partition for the selected solid, blind-pocket, or through-bore regime.
pub fn audit_depth(
    depth_ticks: i32,
    blocks: &[BracketBlock],
) -> Result<BracketCertificate, &'static str> {
    if !(0..=18).contains(&depth_ticks) {
        return Err("pocket_depth_out_of_range");
    }
    let regime = if depth_ticks == 0 {
        Regime::Solid
    } else if depth_ticks < 16 {
        Regime::Blind
    } else {
        Regime::Through
    };
    validate_ids(blocks, regime)?;
    let (canonical, axial_scale) = canonical_g2(blocks, depth_ticks, regime)?;
    let g2 = if regime == Regime::Through {
        audit_section("pierced_square_v1", &canonical)?
    } else {
        audit_section("layered_pocket_v1", &canonical)?
    };

    let incidence = enumerate_faces(blocks)?;
    let (expected_interfaces, expected_exterior) = match regime {
        Regime::Solid => (16, 22),
        Regime::Blind => (24, 30),
        Regime::Through => (4, 16),
    };
    if incidence.interfaces.len() != expected_interfaces
        || incidence.exterior.len() != expected_exterior
    {
        return Err("pocket_incidence_mismatch");
    }
    let design_faces = classify_faces(blocks, &incidence.exterior, regime)?;
    let determinant_lower_bound_mm3 = g2.determinant_lower_bound_mm3 * axial_scale;
    if !determinant_lower_bound_mm3.is_finite() || determinant_lower_bound_mm3 <= 0. {
        return Err("pocket_determinant_certificate");
    }
    Ok(BracketCertificate {
        incidence,
        determinant_lower_bound_mm3,
        radial_error_upper_mm: g2.radial_error_upper_mm,
        design_faces,
    })
}

fn validate_ids(blocks: &[BracketBlock], regime: Regime) -> Result<(), &'static str> {
    let expected: BTreeSet<String> = match regime {
        Regime::Solid => {
            let mut ids = BTreeSet::from(["pocket_core_lower".to_string()]);
            for sector in SECTORS {
                ids.insert(format!("pocket_disk_sector_lower_{sector}"));
                ids.insert(format!("pocket_outer_sector_lower_{sector}"));
            }
            ids
        }
        Regime::Blind => {
            let mut ids = BTreeSet::from(["pocket_core_lower".to_string()]);
            for sector in SECTORS {
                ids.insert(format!("pocket_disk_sector_lower_{sector}"));
                ids.insert(format!("pocket_outer_sector_lower_{sector}"));
                ids.insert(format!("pocket_outer_sector_upper_{sector}"));
            }
            ids
        }
        Regime::Through => SECTORS
            .iter()
            .map(|sector| format!("pocket_outer_sector_upper_{sector}"))
            .collect(),
    };
    if blocks.len() != expected.len() {
        return Err("pocket_block_inventory");
    }
    let mut actual = BTreeSet::new();
    for block in blocks {
        if !actual.insert(block.id.clone()) {
            return Err("duplicate_pocket_block_id");
        }
        if !expected.contains(&block.id) {
            return Err("unexpected_pocket_block_id");
        }
        let required_degree = if block.id == "pocket_core_lower" {
            [1, 1, 1]
        } else {
            [2, 1, 1]
        };
        if block.degrees != required_degree
            || block.controls.len() != required_degree.iter().map(|d| d + 1).product::<usize>()
        {
            return Err("pocket_block_dimensions");
        }
        if block.controls.iter().any(|p| p[3] <= 0) {
            return Err("nonpositive_pocket_weight");
        }
    }
    if actual != expected {
        return Err("pocket_block_inventory");
    }
    Ok(())
}

fn canonical_g2(
    blocks: &[BracketBlock],
    depth_ticks: i32,
    regime: Regime,
) -> Result<(Vec<BracketBlock>, f64), &'static str> {
    match regime {
        Regime::Solid => {
            let mut lower = Vec::with_capacity(9);
            for block in blocks {
                let id = lower_id(&block.id)?;
                lower.push(transform_lower(block, id, 16)?);
            }
            // Build the synthetic upper half from the same exact four square
            // sectors. Its addition lets the independent G2 layered-pocket
            // oracle qualify the full square cross-section at every height.
            let mut canonical = lower.clone();
            for sector in SECTORS {
                let source_id = format!("pocket_outer_sector_{sector}_z0_4");
                let source = lower
                    .iter()
                    .find(|block| block.id == source_id)
                    .ok_or("missing_solid_outer_sector")?;
                canonical.push(shift_z(
                    source,
                    format!("pocket_outer_sector_{sector}_z4_8"),
                    4,
                )?);
            }
            Ok((canonical, 2.0))
        }
        Regime::Blind => {
            let floor_ticks = 16 - depth_ticks;
            let mut canonical = Vec::with_capacity(13);
            for block in blocks {
                if block.id == "pocket_core_lower" || block.id.contains("_lower_") {
                    canonical.push(transform_lower(block, lower_id(&block.id)?, floor_ticks)?);
                } else if block.id.contains("_upper_") {
                    canonical.push(transform_upper(
                        block,
                        upper_id(&block.id)?,
                        floor_ticks,
                        depth_ticks,
                    )?);
                } else {
                    return Err("unknown_blind_pocket_layer");
                }
            }
            let scale = f64::from(floor_ticks.min(depth_ticks)) / 8.0;
            Ok((canonical, scale))
        }
        Regime::Through => {
            let mut canonical = Vec::with_capacity(4);
            for block in blocks {
                let sector = sector_from_id(&block.id)?;
                canonical.push(BracketBlock {
                    id: format!("pierced_sector_{sector}_z0_8"),
                    degrees: block.degrees,
                    controls: block.controls.clone(),
                });
            }
            Ok((canonical, 1.0))
        }
    }
}

fn lower_id(id: &str) -> Result<String, &'static str> {
    if id == "pocket_core_lower" {
        return Ok("pocket_core_z0_4".into());
    }
    if let Some(sector) = id.strip_prefix("pocket_disk_sector_lower_") {
        return Ok(format!("pocket_disk_sector_{sector}_z0_4"));
    }
    if let Some(sector) = id.strip_prefix("pocket_outer_sector_lower_") {
        return Ok(format!("pocket_outer_sector_{sector}_z0_4"));
    }
    Err("invalid_lower_id")
}

fn upper_id(id: &str) -> Result<String, &'static str> {
    let sector = id
        .strip_prefix("pocket_outer_sector_upper_")
        .ok_or("invalid_upper_id")?;
    Ok(format!("pocket_outer_sector_{sector}_z4_8"))
}

fn sector_from_id(id: &str) -> Result<&'static str, &'static str> {
    for sector in SECTORS {
        if id.contains(sector) {
            return Ok(sector);
        }
    }
    Err("unknown_pocket_sector")
}

fn checked_mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("pocket_transform_overflow")
}
fn checked_add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("pocket_transform_overflow")
}
fn checked_sub(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_sub(b).ok_or("pocket_transform_overflow")
}

fn normalize_controls(points: Vec<[i128; 4]>) -> Result<Vec<[i64; 4]>, &'static str> {
    let mut divisor = 0_i128;
    for point in &points {
        for &value in point {
            divisor = gcd_i128(divisor, value)?;
        }
    }
    if divisor <= 0 {
        return Err("zero_pocket_homogeneous_scale");
    }
    let mut result = Vec::with_capacity(points.len());
    for point in points {
        let mut converted = [0_i64; 4];
        for axis in 0..4 {
            converted[axis] =
                i64::try_from(point[axis] / divisor).map_err(|_| "pocket_transform_exceeds_i64")?;
        }
        result.push(converted);
    }
    Ok(result)
}

fn gcd_i128(mut a: i128, mut b: i128) -> Result<i128, &'static str> {
    a = a.checked_abs().ok_or("pocket_transform_overflow")?;
    b = b.checked_abs().ok_or("pocket_transform_overflow")?;
    while b != 0 {
        (a, b) = (b, a % b);
    }
    Ok(a)
}

fn transformed(
    block: &BracketBlock,
    id: String,
    transform: impl Fn([i128; 4]) -> Result<[i128; 4], &'static str>,
) -> Result<BracketBlock, &'static str> {
    let points = block
        .controls
        .iter()
        .map(|p| transform(p.map(i128::from)))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(BracketBlock {
        id,
        degrees: block.degrees,
        controls: normalize_controls(points)?,
    })
}

/// Map z ticks `[0,span]` to the G2 lower layer `[0,4]` mm.
fn transform_lower(
    block: &BracketBlock,
    id: String,
    span_ticks: i32,
) -> Result<BracketBlock, &'static str> {
    if span_ticks <= 0 {
        return Err("zero_height_pocket_layer");
    }
    transformed(block, id, |[x, y, z, w]| {
        Ok([
            checked_mul(i128::from(span_ticks), x)?,
            checked_mul(i128::from(span_ticks), y)?,
            checked_mul(8, z)?,
            checked_mul(i128::from(span_ticks), w)?,
        ])
    })
}

/// Map the blind upper layer `[floorTicks,16]` to G2's `[4,8]` mm.
fn transform_upper(
    block: &BracketBlock,
    id: String,
    floor_ticks: i32,
    span_ticks: i32,
) -> Result<BracketBlock, &'static str> {
    if span_ticks <= 0 {
        return Err("zero_height_pocket_layer");
    }
    transformed(block, id, |[x, y, z, w]| {
        let z = checked_add(
            checked_sub(
                checked_mul(8, z)?,
                checked_mul(4 * i128::from(floor_ticks), w)?,
            )?,
            checked_mul(4 * i128::from(span_ticks), w)?,
        )?;
        Ok([
            checked_mul(i128::from(span_ticks), x)?,
            checked_mul(i128::from(span_ticks), y)?,
            z,
            checked_mul(i128::from(span_ticks), w)?,
        ])
    })
}

fn shift_z(block: &BracketBlock, id: String, dz_mm: i128) -> Result<BracketBlock, &'static str> {
    transformed(block, id, |[x, y, z, w]| {
        Ok([x, y, checked_add(z, checked_mul(dz_mm, w)?)?, w])
    })
}

fn classify_faces(
    blocks: &[BracketBlock],
    exterior: &[BlockSide],
    regime: Regime,
) -> Result<Vec<(BlockSide, String)>, &'static str> {
    let by_id: BTreeMap<&str, &BracketBlock> = blocks
        .iter()
        .map(|block| (block.id.as_str(), block))
        .collect();
    let mut result = Vec::with_capacity(exterior.len());
    for side in exterior {
        let block = by_id
            .get(side.block.as_str())
            .ok_or("unknown_pocket_boundary_face")?;
        let id = block.id.as_str();
        let tag = if id.starts_with("pocket_outer_sector_") && side.axis == 1 {
            if side.upper {
                if regime == Regime::Solid {
                    return Err("unexpected_solid_bore_face");
                }
                "bore_wall".to_string()
            } else {
                let sector = sector_from_id(id)?;
                match sector {
                    "east" => "x10",
                    "north" => "y10",
                    "west" => "x-10",
                    "south" => "y-10",
                    _ => return Err("unknown_pocket_sector"),
                }
                .to_string()
            }
        } else if side.axis == 2 && !side.upper {
            "bottom".to_string()
        } else if side.axis == 2 && side.upper {
            if regime == Regime::Blind
                && (id == "pocket_core_lower" || id.starts_with("pocket_disk_sector_"))
            {
                "pocket_floor".to_string()
            } else {
                "top".to_string()
            }
        } else {
            return Err("unexpected_pocket_boundary_face");
        };
        result.push((side.clone(), tag));
    }
    result.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pocket_depth::{build_pocket_depth, PocketRegime, PocketSelectionEvent};

    fn audited(depth: i32) -> (Vec<BracketBlock>, BracketCertificate) {
        let model = build_pocket_depth(depth).unwrap();
        let blocks: Vec<_> = model.blocks.iter().map(|item| item.block.clone()).collect();
        let certificate = audit_depth(depth, &blocks)
            .unwrap_or_else(|error| panic!("audit rejected depth {depth}: {error}"));
        (blocks, certificate)
    }

    #[test]
    fn certifies_all_nineteen_depth_ticks_and_expected_volumes() {
        for depth in 0..=18 {
            let (blocks, cert) = audited(depth);
            assert!(cert.determinant_lower_bound_mm3.is_finite());
            assert!(cert.determinant_lower_bound_mm3 > 0.);
            assert!(cert.radial_error_upper_mm <= 1e-12);
            let effective = f64::from(depth.min(16)) * 0.5;
            let reference_volume = 3200. - 16. * std::f64::consts::PI * effective;
            let built = build_pocket_depth(depth).unwrap();
            assert_eq!(blocks.len(), built.blocks.len());
            assert!((built.volume_mm3 - reference_volume).abs() <= 4. * f64::EPSILON * 3200.);
            let certified_faces: BTreeSet<_> = cert
                .design_faces
                .iter()
                .map(|(side, label)| (side.clone(), label.as_str()))
                .collect();
            let constructed_faces: BTreeSet<_> = built
                .design_faces
                .iter()
                .map(|face| (face.side.clone(), face.label))
                .collect();
            assert_eq!(certified_faces, constructed_faces);
            for selected in &built.design_faces {
                assert!(cert
                    .incidence
                    .exterior
                    .iter()
                    .any(|side| side == &selected.side));
                assert!(!cert.incidence.interfaces.iter().any(|interface| {
                    interface.first == selected.side || interface.second == selected.side
                }));
            }
            match depth {
                0 => {
                    assert_eq!(built.regime, PocketRegime::Solid);
                    assert!(!cert
                        .design_faces
                        .iter()
                        .any(|(_, label)| label == "pocket_floor"));
                    assert!(!cert
                        .design_faces
                        .iter()
                        .any(|(_, label)| label == "bore_wall"));
                }
                1..=15 => {
                    assert_eq!(built.regime, PocketRegime::BlindPocket);
                    assert!(cert
                        .design_faces
                        .iter()
                        .any(|(_, label)| label == "pocket_floor"));
                    assert!(cert
                        .design_faces
                        .iter()
                        .any(|(_, label)| label == "bore_wall"));
                }
                _ => {
                    assert_eq!(built.regime, PocketRegime::ThroughBore);
                    assert!(!cert
                        .design_faces
                        .iter()
                        .any(|(_, label)| label == "pocket_floor"));
                    assert!(cert
                        .design_faces
                        .iter()
                        .any(|(_, label)| label == "bore_wall"));
                }
            }
        }
    }

    #[test]
    fn verifies_selection_events_and_saturated_domain_signature() {
        let zero = build_pocket_depth(0).unwrap();
        assert_eq!(zero.event, PocketSelectionEvent::ZeroDepth);
        assert!(!zero.selection_state.available.contains(&"pocket_floor"));
        assert!(!zero.selection_state.available.contains(&"bore_wall"));
        let below = build_pocket_depth(15).unwrap();
        let event = build_pocket_depth(16).unwrap();
        assert_eq!(event.event, PocketSelectionEvent::Breakthrough);
        assert!(!event.selection_state.available.contains(&"pocket_floor"));
        assert!(event.selection_state.available.contains(&"bore_wall"));
        let over = build_pocket_depth(18).unwrap();
        assert_eq!(over.event, PocketSelectionEvent::OverdrillSaturated);
        assert_eq!(over.partition_signature, event.partition_signature);
        assert_eq!(
            over.blocks.iter().map(|b| &b.block).collect::<Vec<_>>(),
            event.blocks.iter().map(|b| &b.block).collect::<Vec<_>>()
        );
        assert_ne!(below.partition_signature, event.partition_signature);
    }

    #[test]
    fn verifies_control_jets_within_the_blind_regime() {
        for depth in 2..=14 {
            let center = build_pocket_depth(depth).unwrap();
            let lower = build_pocket_depth(depth - 1).unwrap();
            let upper = build_pocket_depth(depth + 1).unwrap();
            for ((current, lo), hi) in center.blocks.iter().zip(&lower.blocks).zip(&upper.blocks) {
                assert_eq!(current.block.id, lo.block.id);
                assert_eq!(current.block.id, hi.block.id);
                let jet = current.dcontrols_dd.as_ref().unwrap();
                assert_eq!(jet.len(), current.block.controls.len());
                for ((derivative, before), after) in
                    jet.iter().zip(&lo.block.controls).zip(&hi.block.controls)
                {
                    // Two half-millimetre ticks span exactly one physical mm.
                    for axis in 0..4 {
                        assert_eq!(
                            i128::from(after[axis]) - i128::from(before[axis]),
                            i128::from(derivative[axis])
                        );
                    }
                }
            }
        }
        assert!(build_pocket_depth(0)
            .unwrap()
            .blocks
            .iter()
            .all(|b| b.dcontrols_dd.is_none()));
        assert!(build_pocket_depth(16)
            .unwrap()
            .blocks
            .iter()
            .all(|b| b.dcontrols_dd.is_none()));
    }

    #[test]
    fn rejects_wrong_depth_inventory_geometry_and_topology() {
        assert!(audit_depth(-1, &[]).is_err());
        assert!(audit_depth(19, &[]).is_err());
        let (mut blocks, _) = audited(7);
        blocks[0].controls[0][0] += 1;
        assert!(audit_depth(7, &blocks).is_err());
        let (mut blocks, _) = audited(8);
        blocks.pop();
        assert!(audit_depth(8, &blocks).is_err());
        let (mut blocks, _) = audited(16);
        blocks[1].id = blocks[0].id.clone();
        assert!(audit_depth(16, &blocks).is_err());
        let (mut blocks, _) = audited(7);
        blocks.reverse();
        assert!(audit_depth(7, &blocks).is_ok());
    }
}
