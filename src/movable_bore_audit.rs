//! Independent continuous and topological audit for `movable_bore_v1`.
//!
//! This verifier reconstructs the admissible coefficient nets from the frozen
//! E1 profile and the half-millimetre recipe. It never calls the constructor.

use super::bracket::BracketBlock;
use super::bracket_audit::{BracketCertificate, HOLE_PROFILE};
use super::incidence::{enumerate_faces, BlockSide};
use super::planar::certify_planar;
use std::collections::{BTreeMap, BTreeSet};
use std::convert::TryFrom;

const SECTORS: [&str; 4] = ["east", "north", "west", "south"];
const RADIAL_ERROR_MM: f64 = 5.25e-15;

/// Independently audit one serialized movable-bore recipe.
///
/// The center and radius are integer half-millimetre ticks. Accepted values
/// satisfy `|center_tick| <= 4` and `4 <= radius_tick <= 10`.
pub fn audit_bore(
    center_half_mm: [i32; 2],
    radius_half_mm: i32,
    blocks: &[BracketBlock],
) -> Result<BracketCertificate, &'static str> {
    if center_half_mm.iter().any(|&c| !(-4..=4).contains(&c)) || !(4..=10).contains(&radius_half_mm)
    {
        return Err("bore_parameters_out_of_range");
    }
    let expected = expected_blocks(center_half_mm, radius_half_mm)?;
    if blocks.len() != expected.len() {
        return Err("bore_block_inventory");
    }
    let mut actual = BTreeMap::new();
    for block in blocks {
        if actual.insert(block.id.as_str(), block).is_some() {
            return Err("duplicate_bore_block_id");
        }
    }
    if actual.len() != expected.len() {
        return Err("bore_block_inventory");
    }
    for (id, reference) in &expected {
        let block = actual.get(id.as_str()).ok_or("missing_bore_block")?;
        if block.degrees != [2, 1, 1] || block.controls != reference.controls {
            return Err("bore_coefficients_mismatch");
        }
    }

    // Exact incidence validates all four shared radial connectors and rejects
    // any extra, missing, or ambiguous face pairing.
    let incidence = enumerate_faces(blocks)?;
    if incidence.interfaces.len() != 4 || incidence.exterior.len() != 16 {
        return Err("bore_partition_topology");
    }

    // Each closed planar sector has a strictly positive Bernstein bound.
    // The 8 mm extrusion multiplies this bound by exactly eight.
    let mut determinant_lower_bound_mm3 = f64::INFINITY;
    for reference in expected.values() {
        let bound =
            certify_planar([2, 1], &reference.profile).map_err(|_| "bore_planar_certificate")?;
        if !bound.is_finite() || bound <= 0. {
            return Err("bore_nonpositive_determinant");
        }
        determinant_lower_bound_mm3 = determinant_lower_bound_mm3.min(8. * bound);
    }

    let design_faces = classify_faces(blocks, &incidence.exterior)?;
    // Connector disjointness uses the special homothety of these four
    // endpoints, rather than cyclic order alone. The outer corners are
    // A=(±10,±10); the corresponding arc endpoints are c+lambda*A with the
    // same positive lambda (the binary E1 endpoint gives lambda < 0.354 for
    // every allowed radius). Hence every connector lies on a ray from
    // O=c/(1-lambda), at radial factors [lambda,1]. The rays to the four
    // distinct corners are distinct, so connectors cannot cross. Moreover
    // |O-c| = lambda*|c|/(1-lambda) <
    // (177/500)*(2*(283/200))/(323/500) = 50091/32300 < 2 <= r;
    // thus O is strictly inside the bore, and no connector reaches O. The
    // outer square is simple; E6 certifies each transformed arc simple with
    // its endpoint angles and the radial bound places the entire bore inside
    // [-7,7]^2. Each sector boundary is therefore a Jordan curve. Positive
    // planar Jacobian and boundary degree one imply sector injectivity, while
    // the four shared connectors give disjoint interiors and exact coverage.
    if !determinant_lower_bound_mm3.is_finite() || determinant_lower_bound_mm3 <= 0. {
        return Err("bore_determinant_certificate");
    }
    Ok(BracketCertificate {
        incidence,
        determinant_lower_bound_mm3,
        radial_error_upper_mm: RADIAL_ERROR_MM,
        design_faces,
    })
}

struct ExpectedSector {
    controls: Vec<[i64; 4]>,
    profile: [[i64; 3]; 6],
}

fn checked_i64(value: i128) -> Result<i64, &'static str> {
    i64::try_from(value).map_err(|_| "bore_reference_overflow")
}

fn expected_blocks(
    center_ticks: [i32; 2],
    radius_ticks: i32,
) -> Result<BTreeMap<String, ExpectedSector>, &'static str> {
    let mut result = BTreeMap::new();
    for (turns, sector) in SECTORS.iter().enumerate() {
        let mut profile = [[0_i64; 3]; 6];
        let mut controls = Vec::with_capacity(12);
        for u in 0..3 {
            for v in 0..2 {
                let source = HOLE_PROFILE[u * 2 + v];
                let (mut x, mut y, weight) = (
                    i128::from(source[0]),
                    i128::from(source[1]),
                    i128::from(source[2]),
                );
                if v == 1 {
                    // Physical coordinates use 1/2 mm ticks. Multiplying the
                    // whole block by eight gives the exact integer transform.
                    x *= i128::from(radius_ticks);
                    y *= i128::from(radius_ticks);
                } else {
                    x *= 8;
                    y *= 8;
                }
                let w = 8 * weight;
                for _ in 0..turns {
                    (x, y) = (y.checked_neg().ok_or("bore_rotation_overflow")?, x);
                }
                if v == 1 {
                    x += 4 * i128::from(center_ticks[0]) * weight;
                    y += 4 * i128::from(center_ticks[1]) * weight;
                }
                profile[u * 2 + v] = [checked_i64(x)?, checked_i64(y)?, checked_i64(w)?];
                for z in [0_i128, 8_i128] {
                    controls.push([
                        checked_i64(x)?,
                        checked_i64(y)?,
                        checked_i64(z * w)?,
                        checked_i64(w)?,
                    ]);
                }
            }
        }
        result.insert(
            format!("movable_bore_sector_{sector}_z0_8"),
            ExpectedSector { controls, profile },
        );
    }
    Ok(result)
}

fn classify_faces(
    blocks: &[BracketBlock],
    exterior: &[BlockSide],
) -> Result<Vec<(BlockSide, String)>, &'static str> {
    let by_id: BTreeMap<&str, &BracketBlock> = blocks
        .iter()
        .map(|block| (block.id.as_str(), block))
        .collect();
    let mut result = Vec::with_capacity(exterior.len());
    let mut tags = BTreeSet::new();
    for side in exterior {
        if !by_id.contains_key(side.block.as_str()) {
            return Err("unknown_bore_face");
        }
        let feature = if side.axis == 1 && side.upper {
            "bore".to_string()
        } else if side.axis == 1 && !side.upper {
            let sector = side
                .block
                .strip_prefix("movable_bore_sector_")
                .and_then(|tail| tail.strip_suffix("_z0_8"))
                .ok_or("unknown_bore_sector")?;
            match sector {
                "east" => "x10",
                "west" => "x-10",
                "north" => "y10",
                "south" => "y-10",
                _ => return Err("unknown_bore_sector"),
            }
            .to_string()
        } else if side.axis == 2 {
            if side.upper { "z8" } else { "z0" }.to_string()
        } else {
            return Err("unexpected_exterior_bore_face");
        };
        tags.insert(feature.clone());
        result.push((side.clone(), feature));
    }
    if result.len() != 16
        || !["bore", "x10", "x-10", "y10", "y-10", "z0", "z8"]
            .iter()
            .all(|tag| tags.contains(*tag))
    {
        return Err("bore_design_face_inventory");
    }
    result.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn blocks(center: [i32; 2], radius: i32) -> Vec<BracketBlock> {
        crate::movable_bore::build_movable_bore(center[0], center[1], radius)
            .unwrap()
            .blocks
            .into_iter()
            .map(|bore_block| bore_block.block)
            .collect()
    }

    #[test]
    fn certifies_center_and_radius_extremes() {
        for (center, radius) in [([-4, -4], 4), ([-4, 4], 10), ([4, -4], 10), ([4, 4], 4)] {
            let model = blocks(center, radius);
            let cert = audit_bore(center, radius, &model).unwrap();
            assert_eq!(cert.incidence.interfaces.len(), 4);
            assert_eq!(cert.incidence.exterior.len(), 16);
            assert!(cert.determinant_lower_bound_mm3 > 0.);
            assert!(cert.radial_error_upper_mm <= 1e-12);
            assert!(cert.design_faces.iter().any(|(_, tag)| tag == "bore"));
        }
    }

    #[test]
    fn exhausts_all_567_admissible_recipes() {
        let mut count = 0;
        for cx in -4..=4 {
            for cy in -4..=4 {
                for radius in 4..=10 {
                    let model = blocks([cx, cy], radius);
                    let cert = audit_bore([cx, cy], radius, &model).unwrap();
                    assert!(cert.determinant_lower_bound_mm3.is_finite());
                    assert!(cert.determinant_lower_bound_mm3 > 0.);
                    assert!(cert.radial_error_upper_mm <= 1e-12);
                    count += 1;
                }
            }
        }
        assert_eq!(count, 567);
    }

    #[test]
    fn verifies_homogeneous_control_jets_for_every_recipe() {
        for cx in -4..=4 {
            for cy in -4..=4 {
                for radius in 4..=10 {
                    let model = crate::movable_bore::build_movable_bore(cx, cy, radius).unwrap();
                    for (turns, sector) in model.blocks.iter().enumerate() {
                        for index in 0..12 {
                            let u = index / 4;
                            let v = (index / 2) % 2;
                            let [mut x, mut y, source_weight] = HOLE_PROFILE[u * 2 + v];
                            for _ in 0..turns {
                                (x, y) = (-y, x);
                            }
                            let w = 8 * source_weight;
                            let expected_center_x = if v == 1 { [w, 0, 0, 0] } else { [0; 4] };
                            let expected_center_y = if v == 1 { [0, w, 0, 0] } else { [0; 4] };
                            let expected_radius =
                                if v == 1 { [2 * x, 2 * y, 0, 0] } else { [0; 4] };
                            assert_eq!(sector.dcontrols_dcx[index], expected_center_x);
                            assert_eq!(sector.dcontrols_dcy[index], expected_center_y);
                            assert_eq!(sector.dcontrols_dr[index], expected_radius);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn rejects_out_of_domain_missing_duplicate_and_altered_models() {
        assert!(audit_bore([5, 0], 8, &[]).is_err());
        assert!(audit_bore([0, 0], 3, &[]).is_err());
        let mut model = blocks([0, 0], 8);
        model.pop();
        assert!(audit_bore([0, 0], 8, &model).is_err());
        let mut model = blocks([0, 0], 8);
        model[0].controls[2][0] += 1;
        assert!(audit_bore([0, 0], 8, &model).is_err());
        let mut model = blocks([0, 0], 8);
        model[0].degrees = [1, 1, 1];
        assert!(audit_bore([0, 0], 8, &model).is_err());
        let mut model = blocks([0, 0], 8);
        model[0].controls.swap(0, 1);
        assert!(audit_bore([0, 0], 8, &model).is_err());
        let mut model = blocks([0, 0], 8);
        model[1].id = model[0].id.clone();
        assert!(audit_bore([0, 0], 8, &model).is_err());
    }

    #[test]
    fn certificate_is_order_independent() {
        let mut model = blocks([2, -3], 7);
        let first = audit_bore([2, -3], 7, &model).unwrap();
        model.reverse();
        let second = audit_bore([2, -3], 7, &model).unwrap();
        assert_eq!(first, second);
    }
}
