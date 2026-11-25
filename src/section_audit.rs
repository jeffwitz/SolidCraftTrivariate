//! Independent production qualification for the six fixed section templates.

use super::bracket::BracketBlock;
use super::incidence::{enumerate_faces, BlockSide, Incidence};
use super::planar::certify_planar;
use std::collections::{BTreeMap, BTreeSet};

/// Continuous and topological certificate for one fixed section family.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionCertificate {
    pub family: &'static str,
    pub incidence: Incidence,
    pub determinant_lower_bound_mm3: f64,
    pub radial_error_upper_mm: f64,
    /// Analytic target volume for the ideal circle/fillet geometry.
    pub volume_mm3: f64,
    pub design_faces: Vec<(BlockSide, String)>,
}

const SECTORS: [&str; 4] = ["east", "north", "west", "south"];
const ARC_ERROR_MM: f64 = 4.20e-15;

fn hole_profile() -> [[i128; 3]; 6] {
    super::bracket_audit::HOLE_PROFILE.map(|point| point.map(i128::from))
}

fn fillet_profile() -> [[i128; 3]; 6] {
    super::bracket_audit::FILLET_PROFILE.map(|point| point.map(i128::from))
}

#[derive(Clone)]
struct ExpectedBlock {
    block: BracketBlock,
    determinant_lower_bound_mm3: f64,
    radial_error_upper_mm: f64,
}

fn mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b)
        .ok_or("expected section coefficient overflow")
}

fn gcd(mut a: i128, mut b: i128) -> Result<i128, &'static str> {
    a = a.checked_abs().ok_or("coefficient norm overflow")?;
    b = b.checked_abs().ok_or("coefficient norm overflow")?;
    while b != 0 {
        (a, b) = (b, a % b);
    }
    Ok(a)
}

fn normalize(controls: &[[i128; 4]]) -> Result<Vec<[i128; 4]>, &'static str> {
    let mut divisor = 0_i128;
    for point in controls {
        for coefficient in point {
            divisor = gcd(divisor, *coefficient)?;
        }
    }
    if divisor <= 0 {
        return Err("zero homogeneous section net");
    }
    Ok(controls
        .iter()
        .map(|point| point.map(|value| value / divisor))
        .collect())
}

fn box_block(id: String, bounds: [(i64, i64); 3]) -> ExpectedBlock {
    let mut controls = Vec::with_capacity(8);
    for x in [bounds[0].0, bounds[0].1] {
        for y in [bounds[1].0, bounds[1].1] {
            for z in [bounds[2].0, bounds[2].1] {
                controls.push([x, y, z, 1]);
            }
        }
    }
    let volume = (bounds[0].1 - bounds[0].0) as f64
        * (bounds[1].1 - bounds[1].0) as f64
        * (bounds[2].1 - bounds[2].0) as f64;
    ExpectedBlock {
        block: BracketBlock {
            id,
            degrees: [1, 1, 1],
            controls,
        },
        determinant_lower_bound_mm3: volume,
        radial_error_upper_mm: 0.0,
    }
}

fn planar_bound(profile: &[[i128; 3]; 6]) -> Result<f64, &'static str> {
    let mut controls = [[0_i64; 3]; 6];
    for (target, source) in controls.iter_mut().zip(profile) {
        for component in 0..3 {
            target[component] =
                i64::try_from(source[component]).map_err(|_| "reference profile exceeds i64")?;
        }
    }
    certify_planar([2, 1], &controls).map_err(|_| "planar_profile_certificate")
}

fn extruded_profile(
    id: String,
    profile: &[[i128; 3]; 6],
    z0: i64,
    z1: i64,
) -> Result<ExpectedBlock, &'static str> {
    let planar = planar_bound(profile)?;
    let mut controls = Vec::with_capacity(12);
    for [x, y, weight] in profile {
        for z in [z0, z1] {
            controls.push([*x, *y, mul(i128::from(z), *weight)?, *weight]);
        }
    }
    let height = (z1 - z0) as f64;
    Ok(ExpectedBlock {
        block: BracketBlock {
            id,
            degrees: [2, 1, 1],
            controls: controls
                .iter()
                .map(|point| {
                    Ok([
                        i64::try_from(point[0]).map_err(|_| "reference control exceeds i64")?,
                        i64::try_from(point[1]).map_err(|_| "reference control exceeds i64")?,
                        i64::try_from(point[2]).map_err(|_| "reference control exceeds i64")?,
                        i64::try_from(point[3]).map_err(|_| "reference weight exceeds i64")?,
                    ])
                })
                .collect::<Result<Vec<_>, &'static str>>()?,
        },
        determinant_lower_bound_mm3: planar * height,
        radial_error_upper_mm: ARC_ERROR_MM,
    })
}

fn turned(profile: &[[i128; 3]; 6], turns: usize) -> Result<[[i128; 3]; 6], &'static str> {
    let mut result = *profile;
    for point in &mut result {
        for _ in 0..turns {
            (point[0], point[1]) = (
                point[1]
                    .checked_neg()
                    .ok_or("reference rotation overflow")?,
                point[0],
            );
        }
    }
    Ok(result)
}

fn append_sectors(
    expected: &mut Vec<ExpectedBlock>,
    prefix: &str,
    profile: &[[i128; 3]; 6],
    z0: i64,
    z1: i64,
) -> Result<(), &'static str> {
    for (turns, sector) in SECTORS.iter().enumerate() {
        expected.push(extruded_profile(
            format!("{prefix}_{sector}_z{z0}_{z1}"),
            &turned(profile, turns)?,
            z0,
            z1,
        )?);
    }
    Ok(())
}

fn core_disk_profile() -> Result<[[i128; 3]; 6], &'static str> {
    let hole = hole_profile();
    let mut profile = [[0_i128; 3]; 6];
    for angle in 0..3 {
        let circle = hole[angle * 2 + 1];
        let square = hole[angle * 2];
        if square[0] % 10 != 0 || square[1] % 10 != 0 {
            return Err("reference square cannot scale to unit core");
        }
        profile[angle * 2] = circle;
        profile[angle * 2 + 1] = [square[0] / 10, square[1] / 10, square[2]];
    }
    Ok(profile)
}

fn annulus_profile() -> Result<[[i128; 3]; 6], &'static str> {
    let hole = hole_profile();
    let mut profile = [[0_i128; 3]; 6];
    for angle in 0..3 {
        let outer = hole[angle * 2 + 1];
        if outer[0] % 2 != 0 || outer[1] % 2 != 0 {
            return Err("reference radius cannot be halved exactly");
        }
        profile[angle * 2] = outer;
        profile[angle * 2 + 1] = [outer[0] / 2, outer[1] / 2, outer[2]];
    }
    Ok(profile)
}

fn corner_profile(north: bool) -> [[i128; 3]; 6] {
    let fillet = fillet_profile();
    if !north {
        return fillet;
    }
    let mut result = [[0_i128; 3]; 6];
    for i in 0..3 {
        for j in 0..2 {
            let [x, y, weight] = fillet[(2 - i) * 2 + j];
            result[i * 2 + j] = [y, x, weight];
        }
    }
    result
}

fn expected_blocks(family: &str) -> Result<(Vec<ExpectedBlock>, f64), &'static str> {
    let mut expected = Vec::new();
    let hole = hole_profile();
    let volume = match family {
        "reference_rectangle_v1" => {
            expected.push(box_block(
                "rectangle_z0_8".into(),
                [(-10, 10), (-10, 10), (0, 8)],
            ));
            3200.0
        }
        "core_disk_v1" => {
            expected.push(box_block(
                "disk_core_z0_8".into(),
                [(-1, 1), (-1, 1), (0, 8)],
            ));
            append_sectors(&mut expected, "disk_sector", &core_disk_profile()?, 0, 8)?;
            128.0 * std::f64::consts::PI
        }
        "full_annulus_v1" => {
            append_sectors(&mut expected, "annulus_sector", &annulus_profile()?, 0, 8)?;
            96.0 * std::f64::consts::PI
        }
        "pierced_square_v1" => {
            append_sectors(&mut expected, "pierced_sector", &hole, 0, 8)?;
            3200.0 - 128.0 * std::f64::consts::PI
        }
        "rounded_corner_v1" => {
            for (side, north) in [("east", false), ("north", true)] {
                expected.push(extruded_profile(
                    format!("corner_{side}_z0_8"),
                    &corner_profile(north),
                    0,
                    8,
                )?);
            }
            1568.0 - 72.0 * std::f64::consts::PI
        }
        "layered_pocket_v1" => {
            expected.push(box_block(
                "pocket_core_z0_4".into(),
                [(-1, 1), (-1, 1), (0, 4)],
            ));
            append_sectors(
                &mut expected,
                "pocket_disk_sector",
                &core_disk_profile()?,
                0,
                4,
            )?;
            append_sectors(&mut expected, "pocket_outer_sector", &hole, 0, 4)?;
            append_sectors(&mut expected, "pocket_outer_sector", &hole, 4, 8)?;
            3200.0 - 64.0 * std::f64::consts::PI
        }
        _ => return Err("unsupported section family"),
    };
    Ok((expected, volume))
}

fn normalize_input(block: &BracketBlock) -> Result<Vec<[i128; 4]>, &'static str> {
    if block.controls.iter().any(|point| point[3] <= 0) {
        return Err("nonpositive homogeneous weight");
    }
    let controls: Vec<[i128; 4]> = block
        .controls
        .iter()
        .map(|point| point.map(i128::from))
        .collect();
    normalize(&controls)
}

fn validate_net(
    family: &str,
    blocks: &[BracketBlock],
) -> Result<(Vec<ExpectedBlock>, f64), &'static str> {
    let (expected, volume) = expected_blocks(family)?;
    if blocks.len() != expected.len() {
        return Err("section_block_count");
    }
    let expected_by_id: BTreeMap<&str, &ExpectedBlock> = expected
        .iter()
        .map(|item| (item.block.id.as_str(), item))
        .collect();
    let mut seen = BTreeSet::new();
    for block in blocks {
        if !seen.insert(block.id.as_str()) {
            return Err("duplicate_section_identity");
        }
        let reference = expected_by_id
            .get(block.id.as_str())
            .ok_or("unexpected_section_identity")?;
        if block.degrees != reference.block.degrees {
            return Err("section_degree_mismatch");
        }
        if block.controls.len() != reference.block.controls.len()
            || block.controls.len()
                != block
                    .degrees
                    .iter()
                    .map(|degree| degree + 1)
                    .product::<usize>()
        {
            return Err("section_control_count");
        }
        let actual = normalize_input(block)?;
        let expected_controls: Vec<[i128; 4]> = reference
            .block
            .controls
            .iter()
            .map(|point| point.map(i128::from))
            .collect();
        if actual != normalize(&expected_controls)? {
            return Err("section_coefficient_mismatch");
        }
    }
    Ok((expected, volume))
}

fn face_controls(block: &BracketBlock, side: &BlockSide) -> Vec<[i64; 4]> {
    let dimensions = block.degrees.map(|degree| degree + 1);
    let fixed = if side.upper {
        dimensions[side.axis] - 1
    } else {
        0
    };
    let mut controls = Vec::new();
    for i in 0..dimensions[0] {
        for j in 0..dimensions[1] {
            for k in 0..dimensions[2] {
                if [i, j, k][side.axis] == fixed {
                    controls.push(block.controls[(i * dimensions[1] + j) * dimensions[2] + k]);
                }
            }
        }
    }
    controls
}

fn tag_face(family: &str, block: &BracketBlock, side: &BlockSide) -> Result<String, &'static str> {
    let id = block.id.as_str();
    if family == "rounded_corner_v1" && side.axis == 1 && side.upper {
        return Ok("fillet_r6".into());
    }
    if family == "full_annulus_v1" && side.axis == 1 {
        return Ok(if side.upper { "circle_r2" } else { "circle_r4" }.into());
    }
    if family == "core_disk_v1" && id.starts_with("disk_sector_") && side.axis == 1 && !side.upper {
        return Ok("circle_r4".into());
    }
    if family == "pierced_square_v1" && side.axis == 1 && side.upper {
        return Ok("circle_r4".into());
    }
    if family == "layered_pocket_v1" {
        if id.starts_with("pocket_outer_sector_")
            && id.ends_with("_z4_8")
            && side.axis == 1
            && side.upper
        {
            return Ok("circle_r4".into());
        }
        if (id == "pocket_core_z0_4" || id.starts_with("pocket_disk_sector_"))
            && side.axis == 2
            && side.upper
        {
            return Ok("pocket_floor_z4".into());
        }
    }

    let plane_levels: [(&str, usize, &[i64]); 3] = if family == "rounded_corner_v1" {
        [("x", 0, &[0, 14]), ("y", 1, &[0, 14]), ("z", 2, &[0, 8])]
    } else if family == "layered_pocket_v1" {
        [
            ("x", 0, &[-10, 10]),
            ("y", 1, &[-10, 10]),
            ("z", 2, &[0, 4, 8]),
        ]
    } else {
        [
            ("x", 0, &[-10, 10]),
            ("y", 1, &[-10, 10]),
            ("z", 2, &[0, 8]),
        ]
    };
    let controls = face_controls(block, side);
    let mut tags = Vec::new();
    for (axis_name, coordinate, levels) in plane_levels {
        for level in levels {
            if controls.iter().all(|point| {
                i128::from(point[coordinate]) == i128::from(*level) * i128::from(point[3])
            }) {
                tags.push(format!("{axis_name}{level}"));
            }
        }
    }
    if tags.len() != 1 {
        return Err("unclassified_section_boundary");
    }
    Ok(tags.remove(0))
}

fn design_faces(
    family: &str,
    blocks: &[BracketBlock],
    incidence: &Incidence,
) -> Result<Vec<(BlockSide, String)>, &'static str> {
    let by_id: BTreeMap<&str, &BracketBlock> = blocks
        .iter()
        .map(|block| (block.id.as_str(), block))
        .collect();
    let mut result = Vec::with_capacity(incidence.exterior.len());
    for side in &incidence.exterior {
        let block = by_id
            .get(side.block.as_str())
            .ok_or("unknown section boundary block")?;
        result.push((side.clone(), tag_face(family, block, side)?));
    }
    result.sort_by(|left, right| left.0.cmp(&right.0));
    Ok(result)
}

/// Audit a fixed section family against independent profiles and topology.
///
/// The four curved templates are qualified only for these frozen profiles:
/// their planar contours are simple and nested without overlap, and their
/// positive planar Jacobian certificates preserve orientation. The radial
/// profile differs from its ideal circular arc by at most `ARC_ERROR_MM`;
/// the E6 arc's endpoint/tangent construction and this bound preserve the
/// contour ordering. Sector interiors therefore tile the disk/annulus once
/// (topological degree one), while the disjoint square/core regions fill the
/// remaining stated domains. Extrusion over disjoint z intervals preserves
/// that injective tiling. This is a template-specific transfer argument,
/// not a certificate for arbitrary blocks or arbitrary perturbed profiles.
pub fn audit_section(
    family: &str,
    blocks: &[BracketBlock],
) -> Result<SectionCertificate, &'static str> {
    let (expected, volume_mm3) = validate_net(family, blocks)?;
    let expected_nets: Vec<BracketBlock> = expected.iter().map(|item| item.block.clone()).collect();
    let incidence = enumerate_faces(blocks)?;
    let expected_incidence = enumerate_faces(&expected_nets)?;
    if incidence != expected_incidence {
        return Err("section_topology_mismatch");
    }

    let determinant_lower_bound_mm3 = expected
        .iter()
        .map(|item| item.determinant_lower_bound_mm3)
        .fold(f64::INFINITY, f64::min);
    let radial_error_upper_mm = expected
        .iter()
        .map(|item| item.radial_error_upper_mm)
        .fold(0.0, f64::max);
    if !determinant_lower_bound_mm3.is_finite() || determinant_lower_bound_mm3 <= 0.0 {
        return Err("section_determinant_certificate");
    }
    let design_faces = design_faces(family, blocks, &incidence)?;
    Ok(SectionCertificate {
        family: match family {
            "reference_rectangle_v1" => "reference_rectangle_v1",
            "core_disk_v1" => "core_disk_v1",
            "full_annulus_v1" => "full_annulus_v1",
            "pierced_square_v1" => "pierced_square_v1",
            "rounded_corner_v1" => "rounded_corner_v1",
            "layered_pocket_v1" => "layered_pocket_v1",
            _ => return Err("unsupported section family"),
        },
        incidence,
        determinant_lower_bound_mm3,
        radial_error_upper_mm,
        volume_mm3,
        design_faces,
    })
}
