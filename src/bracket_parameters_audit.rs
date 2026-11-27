//! Independent audit for the fixed two-radius bracket parameter family.
//!
//! The audit inverts the exact inner-row homothety to the frozen E6 bracket
//! profile, then delegates only the restored reference net to its existing
//! independent fixed-family audit. The eight edited planar maps are certified
//! separately. Within qx,qy in `[40,80]` sixteenths of a millimetre, each inner
//! loop is a positive homothety (scale 5/8..5/4) of the simple E6 loop, remains
//! strictly inside its unchanged 10 mm square boundary, and joins that boundary
//! on the same four radial rays. Positive planar Jacobians and the simple
//! one-to-one boundary give degree one on each sector. The four sectors tile
//! the annular square once; extrusion over the independent 8 mm axis preserves
//! injectivity. This transfer is specific to this frozen template and range.
//! A positive local determinant alone would not establish injectivity.

use crate::bracket::BracketBlock;
use crate::bracket_audit::{audit_bracket, BracketCertificate};
use crate::incidence::{enumerate_faces, BlockSide, Incidence};
use crate::planar::certify_planar;
use std::collections::BTreeMap;
use std::convert::TryFrom;

const MIN_RADIUS_TICKS: u8 = 40;
const MAX_RADIUS_TICKS: u8 = 80;
const REFERENCE_RADIUS_TICKS: i128 = 64;
const SECTORS: [&str; 4] = ["east", "north", "south", "west"];
const LAYERS_AND_IDS: [(&str, &str); 2] = [("x", "hole_x_sector_"), ("y", "hole_y_sector_")];

/// Complete audit report for a parameterized bracket.
#[derive(Debug, Clone, PartialEq)]
pub struct ParameterizedBracketCertificate {
    /// Certificate for the restored fixed E6 reference and inherited regions.
    pub certificate: BracketCertificate,
    /// Exact family target volume, evaluated in mm³.
    pub volume_mm3: f64,
    pub radii_sixteenths_mm: [u8; 2],
    /// Canonical block IDs whose normalized coefficients change per radius.
    pub modified_blocks: [Vec<String>; 2],
}

/// Independently audit serialized bracket controls without invoking the
/// parameterized constructor. Input radii are sixteenths of a millimetre.
pub fn audit_bracket_parameters(
    radius_x_sixteenths_mm: u8,
    radius_y_sixteenths_mm: u8,
    blocks: &[BracketBlock],
) -> Result<ParameterizedBracketCertificate, &'static str> {
    let radii = [radius_x_sixteenths_mm, radius_y_sixteenths_mm];
    if radii
        .iter()
        .any(|q| !(MIN_RADIUS_TICKS..=MAX_RADIUS_TICKS).contains(q))
    {
        return Err("hole radius must be between 2.5 and 5 mm");
    }
    if blocks.len() != 36 {
        return Err("block_count");
    }

    let mut restored = Vec::with_capacity(blocks.len());
    let mut hole_counts = [0usize; 2];
    for block in blocks {
        let mut reference = block.clone();
        if let Some(axis) = hole_axis(&block.id) {
            hole_counts[axis] += 1;
            restore_inner_row(&mut reference, radii[axis])?;
        }
        restored.push(reference);
    }
    if hole_counts != [4, 4] {
        return Err("hole_sector_inventory");
    }

    // This certifies the inverse-transformed coefficient nets against E6's
    // independent profile constants, block inventory, topology and selections.
    let mut certificate = audit_bracket(&restored)?;
    let candidate_incidence = enumerate_faces(blocks).map_err(|_| "modified_partition")?;
    if candidate_incidence.interfaces.len() != 55 || candidate_incidence.exterior.len() != 106 {
        return Err("modified_partition_inventory");
    }
    if candidate_incidence != certificate.incidence {
        return Err("modified_partition_mismatch");
    }

    let candidate_selections = classify_design_faces(blocks, &candidate_incidence)?;
    if candidate_selections != certificate.design_faces {
        return Err("design_face_identity_mismatch");
    }

    // The E6 lower bound covers unchanged blocks. Certify every edited planar
    // patch independently and multiply by its exact 8 mm extrusion length.
    let mut edited_planar_lower = f64::INFINITY;
    for block in blocks.iter().filter(|b| hole_axis(&b.id).is_some()) {
        let profiles = planar_profiles(block)?;
        for profile in profiles {
            let lower =
                certify_planar([2, 1], &profile).map_err(|_| "modified_planar_certificate")?;
            if !lower.is_finite() || lower <= 0.0 {
                return Err("modified_planar_certificate");
            }
            edited_planar_lower = edited_planar_lower.min(lower);
        }
    }
    if !edited_planar_lower.is_finite() {
        return Err("modified_planar_certificate");
    }
    certificate.determinant_lower_bound_mm3 = certificate
        .determinant_lower_bound_mm3
        .min(edited_planar_lower * 8.0);
    if certificate.determinant_lower_bound_mm3 <= 0.0
        || !certificate.determinant_lower_bound_mm3.is_finite()
    {
        return Err("determinant_certificate");
    }
    let radius_factor = ((radii[0].max(radii[1]) as f64) / REFERENCE_RADIUS_TICKS as f64).max(1.0);
    certificate.radial_error_upper_mm *= radius_factor;

    let rx = radii[0] as f64 / 16.0;
    let ry = radii[1] as f64 / 16.0;
    let volume_mm3 =
        34_080.0 - 360.0 * std::f64::consts::PI - 8.0 * std::f64::consts::PI * (rx * rx + ry * ry);
    if !volume_mm3.is_finite() || volume_mm3 <= 0.0 {
        return Err("volume_range");
    }

    let modified_blocks = std::array::from_fn(|axis| {
        if radii[axis] == 64 {
            Vec::new()
        } else {
            hole_ids(axis)
        }
    });
    Ok(ParameterizedBracketCertificate {
        certificate,
        volume_mm3,
        radii_sixteenths_mm: radii,
        modified_blocks,
    })
}

/// Verify serialized analytic derivatives in millimetres per millimetre.
/// The sparse map format is the constructor's two-axis JSON representation.
pub fn audit_bracket_parameter_derivatives(
    radius_x_sixteenths_mm: u8,
    radius_y_sixteenths_mm: u8,
    blocks: &[BracketBlock],
    derivatives: &[BTreeMap<String, Vec<[i64; 4]>>; 2],
) -> Result<(), &'static str> {
    audit_bracket_parameters(radius_x_sixteenths_mm, radius_y_sixteenths_mm, blocks)?;
    let radii = [radius_x_sixteenths_mm, radius_y_sixteenths_mm];
    for axis in 0..2 {
        let ids = hole_ids(axis);
        if derivatives[axis].len() != ids.len()
            || derivatives[axis]
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>()
                != ids.iter().map(String::as_str).collect::<Vec<_>>()
        {
            return Err("parameter_derivative_inventory");
        }
        for id in &ids {
            let block = blocks
                .iter()
                .find(|block| block.id == *id)
                .ok_or("missing_hole_block")?;
            let jet = derivatives[axis]
                .get(id)
                .ok_or("missing_parameter_derivative")?;
            if jet.len() != 12 || block.degrees != [2, 1, 1] {
                return Err("parameter_derivative_shape");
            }
            let block_axis = hole_axis(&block.id).ok_or("unexpected_parameter_derivative_id")?;
            if block_axis != axis {
                return Err("parameter_derivative_axis_mismatch");
            }
            for (index, (control, actual)) in block.controls.iter().zip(jet).enumerate() {
                let radial_row = (index / 2) % 2;
                let k = index % 2;
                let weight = i128::from(control[3]);
                if weight <= 0 {
                    return Err("nonpositive_homogeneous_weight");
                }
                let center = center_times_weight(axis, k, weight)?;
                let expected = if radial_row == 0 {
                    [0i128; 4]
                } else {
                    let mut values = [0i128; 4];
                    for component in 0..3 {
                        let delta = i128::from(control[component]) - center[component];
                        let scaled = delta
                            .checked_mul(16)
                            .ok_or("parameter_derivative_overflow")?;
                        if scaled % i128::from(radii[axis]) != 0 {
                            return Err("parameter_derivative_not_integral");
                        }
                        values[component] = scaled / i128::from(radii[axis]);
                    }
                    values
                };
                for component in 0..4 {
                    let expected_i64 = i64::try_from(expected[component])
                        .map_err(|_| "parameter_derivative_range")?;
                    if actual[component] != expected_i64 {
                        return Err("parameter_derivative_mismatch");
                    }
                }
            }
        }
    }
    Ok(())
}

fn hole_axis(id: &str) -> Option<usize> {
    if id.starts_with("hole_x_sector_") {
        Some(0)
    } else if id.starts_with("hole_y_sector_") {
        Some(1)
    } else {
        None
    }
}

fn hole_ids(axis: usize) -> Vec<String> {
    let prefix = LAYERS_AND_IDS[axis].1;
    SECTORS
        .iter()
        .map(|sector| format!("{prefix}{sector}_z10_30"))
        .collect()
}

fn center_times_weight(axis: usize, k: usize, weight: i128) -> Result<[i128; 3], &'static str> {
    let axial = if k == 0 { 0i128 } else { 8i128 };
    let mul = |a: i128, b: i128| a.checked_mul(b).ok_or("homogeneous inverse overflow");
    Ok(if axis == 0 {
        [mul(40, weight)?, mul(8 - axial, weight)?, mul(20, weight)?]
    } else {
        [mul(axial, weight)?, mul(32, weight)?, mul(20, weight)?]
    })
}

fn restore_inner_row(block: &mut BracketBlock, radius_ticks: u8) -> Result<(), &'static str> {
    let axis = hole_axis(&block.id).ok_or("not a hole sector")?;
    if block.degrees != [2, 1, 1] || block.controls.len() != 12 {
        return Err("hole_block_shape");
    }
    for (index, control) in block.controls.iter_mut().enumerate() {
        let radial_row = (index / 2) % 2;
        if radial_row == 0 {
            continue;
        }
        let k = index % 2;
        let weight = i128::from(control[3]);
        if weight <= 0 {
            return Err("nonpositive_homogeneous_weight");
        }
        let center = center_times_weight(axis, k, weight)?;
        let radius = i128::from(radius_ticks);
        for component in 0..3 {
            let delta = i128::from(control[component])
                .checked_sub(center[component])
                .ok_or("homogeneous inverse overflow")?;
            let numerator = delta
                .checked_mul(64)
                .ok_or("homogeneous inverse overflow")?;
            if numerator % radius != 0 {
                return Err("inexact_inner_row_inverse");
            }
            let restored = center[component]
                .checked_add(numerator / radius)
                .ok_or("homogeneous inverse overflow")?;
            control[component] = i64::try_from(restored).map_err(|_| "restored_control_range")?;
        }
    }
    Ok(())
}

fn planar_profiles(block: &BracketBlock) -> Result<[Vec<[i64; 3]>; 2], &'static str> {
    let axis = hole_axis(&block.id).ok_or("not a hole sector")?;
    if block.degrees != [2, 1, 1] || block.controls.len() != 12 {
        return Err("hole_block_shape");
    }
    let mut profiles: [Vec<[i64; 3]>; 2] = std::array::from_fn(|_| Vec::with_capacity(6));
    for (k, profile) in profiles.iter_mut().enumerate() {
        for i in 0..3 {
            for j in 0..2 {
                let control = block.controls[(i * 2 + j) * 2 + k];
                let weight = i128::from(control[3]);
                if weight <= 0 {
                    return Err("nonpositive_homogeneous_weight");
                }
                let center = center_times_weight(axis, k, weight)?;
                let first_axis = if axis == 0 { 0 } else { 1 };
                let x = i64::try_from(i128::from(control[first_axis]) - center[first_axis])
                    .map_err(|_| "planar_control_range")?;
                let z = i64::try_from(i128::from(control[2]) - center[2])
                    .map_err(|_| "planar_control_range")?;
                profile.push([x, z, control[3]]);
            }
        }
    }
    if profiles[0] != profiles[1] {
        return Err("extruded_profile_mismatch");
    }
    Ok(profiles)
}

fn face_controls(block: &BracketBlock, side: &BlockSide) -> Vec<[i64; 4]> {
    let [du, dv, dw] = block.degrees;
    let dimensions = [du + 1, dv + 1, dw + 1];
    let fixed = if side.upper {
        dimensions[side.axis] - 1
    } else {
        0
    };
    let mut face = Vec::new();
    for i in 0..dimensions[0] {
        for j in 0..dimensions[1] {
            for k in 0..dimensions[2] {
                let index = [i, j, k];
                if index[side.axis] == fixed {
                    face.push(block.controls[(i * dimensions[1] + j) * dimensions[2] + k]);
                }
            }
        }
    }
    face
}

fn classify_design_faces(
    blocks: &[BracketBlock],
    incidence: &Incidence,
) -> Result<Vec<(BlockSide, String)>, &'static str> {
    let by_id: BTreeMap<&str, &BracketBlock> = blocks.iter().map(|b| (b.id.as_str(), b)).collect();
    let mut result = Vec::with_capacity(incidence.exterior.len());
    for side in &incidence.exterior {
        let block = by_id
            .get(side.block.as_str())
            .ok_or("unknown_exterior_block")?;
        let id = block.id.as_str();
        let feature = if side.axis == 1 && side.upper && id.starts_with("hole_x_sector_") {
            "hole_x".to_string()
        } else if side.axis == 1 && side.upper && id.starts_with("hole_y_sector_") {
            "hole_y".to_string()
        } else if side.axis == 1 && side.upper && id.starts_with("elbow_") {
            "inside_fillet".to_string()
        } else {
            let controls = face_controls(block, side);
            let planes: [(&str, usize, &[i64]); 3] = [
                ("x", 0, &[0, 8, 60]),
                ("y", 1, &[0, 8, 50]),
                ("z", 2, &[0, 40]),
            ];
            let mut matches = Vec::new();
            for (name, coordinate, levels) in planes {
                for level in levels {
                    if controls
                        .iter()
                        .all(|p| i128::from(p[coordinate]) == i128::from(*level) * i128::from(p[3]))
                    {
                        matches.push(format!("{name}{level}"));
                    }
                }
            }
            if matches.len() != 1 {
                return Err("unclassified_boundary");
            }
            matches[0].clone()
        };
        result.push((side.clone(), feature));
    }
    result.sort_by(|a, b| a.0.cmp(&b.0));
    if result.len() != 106 {
        return Err("design_face_inventory");
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bracket_parameters::build_bracket_parameters;

    #[test]
    fn accepts_nominal_and_independent_x_radius_edit() {
        for (rx, ry) in [(64, 64), (56, 64), (64, 56), (40, 80), (80, 40)] {
            let model = build_bracket_parameters(rx, ry).unwrap();
            let audit = audit_bracket_parameters(rx, ry, &model.blocks).unwrap();
            assert_eq!(audit.certificate.incidence.interfaces.len(), 55);
            assert_eq!(audit.certificate.incidence.exterior.len(), 106);
            assert_eq!(audit.radii_sixteenths_mm, [rx, ry]);
            assert_eq!(audit.modified_blocks[0].len(), if rx == 64 { 0 } else { 4 });
            assert_eq!(audit.modified_blocks[1].len(), if ry == 64 { 0 } else { 4 });
            let expected = 34_080.0
                - 360.0 * std::f64::consts::PI
                - 8.0
                    * std::f64::consts::PI
                    * ((rx as f64 / 16.0).powi(2) + (ry as f64 / 16.0).powi(2));
            assert!((audit.volume_mm3 - expected).abs() < 1e-10);
            assert!(audit.certificate.determinant_lower_bound_mm3 > 0.0);
        }
    }

    #[test]
    fn audits_jets_and_rejects_a_mutated_derivative() {
        let model = build_bracket_parameters(56, 64).unwrap();
        audit_bracket_parameter_derivatives(56, 64, &model.blocks, &model.derivatives).unwrap();
        let mut bad = model.derivatives.clone();
        let key = bad[0].keys().next().unwrap().clone();
        bad[0].get_mut(&key).unwrap()[3][0] += 1;
        assert_eq!(
            audit_bracket_parameter_derivatives(56, 64, &model.blocks, &bad),
            Err("parameter_derivative_mismatch")
        );
    }

    #[test]
    fn rejects_invalid_radius_and_mutated_blocks_or_interfaces() {
        let model = build_bracket_parameters(56, 64).unwrap();
        assert!(audit_bracket_parameters(39, 64, &model.blocks).is_err());
        let mut bad = model.blocks.clone();
        bad[0].controls[0][0] += 1;
        assert!(audit_bracket_parameters(56, 64, &bad).is_err());
        let mut bad = model.blocks.clone();
        let block = bad
            .iter_mut()
            .find(|b| b.id.starts_with("hole_x_sector_"))
            .unwrap();
        block.controls[2][0] += 7;
        assert!(audit_bracket_parameters(56, 64, &bad).is_err());
    }
}
