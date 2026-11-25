//! Independent fixed-family audit for the filleted two-hole bracket.

use super::bracket::BracketBlock;
use super::incidence::{enumerate_faces, BlockSide, Incidence};
use super::planar::certify_planar;
use std::collections::{BTreeMap, BTreeSet};

/// Result of auditing the fixed bracket family.
#[derive(Debug, Clone, PartialEq)]
pub struct BracketCertificate {
    pub incidence: Incidence,
    pub determinant_lower_bound_mm3: f64,
    pub radial_error_upper_mm: f64,
    pub design_faces: Vec<(BlockSide, String)>,
}

const LAYERS: [(i64, i64); 3] = [(0, 10), (10, 30), (30, 40)];
const SECTORS: [&str; 4] = ["east", "north", "west", "south"];
const RADIAL_ERROR_BOUND_MM: f64 = 4.20e-15;

// E6-qualified local reference profiles, independently frozen from the
// analytical radius-6 fillet and radius-4 plate sector. Triples are (wx,wy,w).
// Their exact binary64 source values were scaled to integers by 2^52 and 2^53.
pub(crate) const FILLET_PROFILE: [[i64; 3]; 6] = [
    [0, 0, 4_503_599_627_370_496],
    [
        43_943_239_765_609_632,
        43_943_239_765_609_632,
        4_503_599_627_370_496,
    ],
    [31_525_197_391_593_472, 0, 4_503_599_627_370_496],
    [
        47_910_251_476_541_136,
        33_286_268_146_824_472,
        4_160_783_518_353_059,
    ],
    [63_050_394_783_186_944, 0, 4_503_599_627_370_496],
    [
        63_050_394_783_186_944,
        36_028_797_018_963_968,
        4_503_599_627_370_496,
    ],
];

pub(crate) const HOLE_PROFILE: [[i64; 3]; 6] = [
    [
        90_071_992_547_409_920,
        -90_071_992_547_409_920,
        9_007_199_254_740_992,
    ],
    [
        25_476_206_690_103_092,
        -25_476_206_690_103_092,
        9_007_199_254_740_992,
    ],
    [90_071_992_547_409_920, 0, 9_007_199_254_740_992],
    [36_028_797_018_963_976, 0, 6_369_051_672_525_773],
    [
        90_071_992_547_409_920,
        90_071_992_547_409_920,
        9_007_199_254_740_992,
    ],
    [
        25_476_206_690_103_092,
        25_476_206_690_103_092,
        9_007_199_254_740_992,
    ],
];

#[derive(Clone)]
struct ExpectedBlock {
    id: String,
    degrees: [usize; 3],
    controls: Vec<[i128; 4]>,
}

fn add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("expected coefficient overflow")
}

fn sub(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_sub(b).ok_or("expected coefficient overflow")
}

fn mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("expected coefficient overflow")
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
        return Err("zero homogeneous control net");
    }
    Ok(controls
        .iter()
        .map(|point| point.map(|value| value / divisor))
        .collect())
}

fn expected_box(id: String, x: (i64, i64), y: (i64, i64), z: (i64, i64)) -> ExpectedBlock {
    let mut controls = Vec::with_capacity(8);
    for xx in [x.0, x.1] {
        for yy in [y.0, y.1] {
            for zz in [z.0, z.1] {
                controls.push([i128::from(xx), i128::from(yy), i128::from(zz), 1_i128]);
            }
        }
    }
    ExpectedBlock {
        id,
        degrees: [1, 1, 1],
        controls,
    }
}

fn expected_elbow(side: &str, z0: i64, z1: i64) -> Result<ExpectedBlock, &'static str> {
    let mut controls = Vec::with_capacity(12);
    for i in 0..3 {
        let source_i = if side == "east" { i } else { 2 - i };
        for j in 0..2 {
            let [sx, sy, sw] = FILLET_PROFILE[source_i * 2 + j];
            let (x, y) = if side == "east" { (sx, sy) } else { (sy, sx) };
            let weight = i128::from(sw);
            for z in [z0, z1] {
                controls.push([
                    i128::from(x),
                    i128::from(y),
                    mul(i128::from(z), weight)?,
                    weight,
                ]);
            }
        }
    }
    Ok(ExpectedBlock {
        id: format!("elbow_{side}_z{z0}_{z1}"),
        degrees: [2, 1, 1],
        controls,
    })
}

fn rotate(x: i128, y: i128, turns: usize) -> Result<(i128, i128), &'static str> {
    let mut pair = (x, y);
    for _ in 0..turns {
        pair = (
            pair.1.checked_neg().ok_or("expected rotation overflow")?,
            pair.0,
        );
    }
    Ok(pair)
}

fn expected_hole(leg: &str, sector: &str, turns: usize) -> Result<ExpectedBlock, &'static str> {
    let mut controls = Vec::with_capacity(12);
    for [sx, sy, sw] in HOLE_PROFILE {
        let weight = i128::from(sw);
        let (a, b) = rotate(i128::from(sx), i128::from(sy), turns)?;
        for c in [0_i128, 8_i128] {
            let xyz = if leg == "x" {
                (
                    add(mul(40, weight)?, a)?,
                    sub(mul(8, weight)?, mul(c, weight)?)?,
                    add(mul(20, weight)?, b)?,
                )
            } else {
                (
                    mul(c, weight)?,
                    add(mul(32, weight)?, a)?,
                    add(mul(20, weight)?, b)?,
                )
            };
            controls.push([xyz.0, xyz.1, xyz.2, weight]);
        }
    }
    Ok(ExpectedBlock {
        id: format!("hole_{leg}_sector_{sector}_z10_30"),
        degrees: [2, 1, 1],
        controls,
    })
}

fn expected_blocks() -> Result<BTreeMap<String, ExpectedBlock>, &'static str> {
    let mut blocks = BTreeMap::new();
    let mut insert = |block: ExpectedBlock| -> Result<(), &'static str> {
        if blocks.insert(block.id.clone(), block).is_some() {
            return Err("duplicate expected block identity");
        }
        Ok(())
    };

    for (z0, z1) in LAYERS {
        insert(expected_elbow("east", z0, z1)?)?;
        insert(expected_elbow("north", z0, z1)?)?;

        let suffix = format!("z{z0}_{z1}");
        insert(expected_box(
            format!("connector_x_14_16_{suffix}"),
            (14, 16),
            (0, 8),
            (z0, z1),
        ))?;
        insert(expected_box(
            format!("connector_y_14_16_{suffix}"),
            (0, 8),
            (14, 16),
            (z0, z1),
        ))?;

        let x_intervals: &[(i64, i64)] = if (z0, z1) == (10, 30) {
            &[(16, 30), (50, 60)]
        } else {
            &[(16, 30), (30, 50), (50, 60)]
        };
        for &(x0, x1) in x_intervals {
            insert(expected_box(
                format!("arm_x_{x0}_{x1}_{suffix}"),
                (x0, x1),
                (0, 8),
                (z0, z1),
            ))?;
        }

        let y_intervals: &[(i64, i64)] = if (z0, z1) == (10, 30) {
            &[(16, 22), (42, 50)]
        } else {
            &[(16, 22), (22, 42), (42, 50)]
        };
        for &(y0, y1) in y_intervals {
            insert(expected_box(
                format!("arm_y_{y0}_{y1}_{suffix}"),
                (0, 8),
                (y0, y1),
                (z0, z1),
            ))?;
        }
    }

    for leg in ["x", "y"] {
        for (turns, sector) in SECTORS.iter().enumerate() {
            insert(expected_hole(leg, sector, turns)?)?;
        }
    }
    if blocks.len() != 36 {
        return Err("expected block inventory is inconsistent");
    }
    Ok(blocks)
}

fn normalized_input(block: &BracketBlock) -> Result<Vec<[i128; 4]>, &'static str> {
    let controls: Vec<[i128; 4]> = block
        .controls
        .iter()
        .map(|point| point.map(i128::from))
        .collect();
    if controls.iter().any(|point| point[3] <= 0) {
        return Err("nonpositive homogeneous weight");
    }
    normalize(&controls)
}

fn validate_blocks(input: &[BracketBlock]) -> Result<(), &'static str> {
    let expected = expected_blocks()?;
    if input.len() != expected.len() {
        return Err("block_count");
    }
    let mut seen = BTreeSet::new();
    for block in input {
        if !seen.insert(block.id.as_str()) {
            return Err("duplicate_block_identity");
        }
        let reference = expected.get(&block.id).ok_or("unexpected_block_identity")?;
        if block.degrees != reference.degrees {
            return Err("block_degree_mismatch");
        }
        if block.controls.len() != reference.controls.len()
            || block.controls.len()
                != block
                    .degrees
                    .iter()
                    .map(|degree| degree + 1)
                    .product::<usize>()
        {
            return Err("control_count");
        }
        let actual = normalized_input(block)?;
        let expected_controls = normalize(&reference.controls)?;
        if actual != expected_controls {
            return Err("coefficient_mismatch");
        }
    }
    Ok(())
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
    input: &[BracketBlock],
    incidence: &Incidence,
) -> Result<Vec<(BlockSide, String)>, &'static str> {
    let by_id: BTreeMap<&str, &BracketBlock> = input
        .iter()
        .map(|block| (block.id.as_str(), block))
        .collect();
    let mut result = Vec::with_capacity(incidence.exterior.len());
    for side in &incidence.exterior {
        let block = by_id
            .get(side.block.as_str())
            .ok_or("unknown exterior block")?;
        let name = block.id.as_str();
        let feature = if side.axis == 1 && side.upper && name.starts_with("hole_x_") {
            "hole_x".to_string()
        } else if side.axis == 1 && side.upper && name.starts_with("hole_y_") {
            "hole_y".to_string()
        } else if side.axis == 1 && side.upper && name.starts_with("elbow_") {
            "inside_fillet".to_string()
        } else {
            let controls = face_controls(block, side);
            let planes: [(&str, usize, &[i64]); 3] = [
                ("x", 0, &[0, 8, 60]),
                ("y", 1, &[0, 8, 50]),
                ("z", 2, &[0, 40]),
            ];
            let mut matches = Vec::new();
            for (axis_name, coordinate, levels) in planes {
                for level in levels {
                    if controls.iter().all(|point| {
                        i128::from(point[coordinate]) == i128::from(*level) * i128::from(point[3])
                    }) {
                        matches.push(format!("{axis_name}{level}"));
                    }
                }
            }
            if matches.len() != 1 {
                return Err("unclassified_boundary");
            }
            matches.remove(0)
        };
        result.push((side.clone(), feature));
    }
    result.sort_by(|left, right| left.0.cmp(&right.0));
    let feature_count: BTreeSet<&str> =
        result.iter().map(|(_, feature)| feature.as_str()).collect();
    if result.len() != 106 || feature_count.len() != 11 {
        return Err("design_face_inventory");
    }
    Ok(result)
}

fn profile_determinants() -> Result<(f64, f64), &'static str> {
    let fillet = certify_planar([2, 1], &FILLET_PROFILE).map_err(|_| "fillet_certificate")?;
    let hole = certify_planar([2, 1], &HOLE_PROFILE).map_err(|_| "hole_certificate")?;
    if !fillet.is_finite() || !hole.is_finite() || fillet <= 0.0 || hole <= 0.0 {
        return Err("nonpositive_profile_certificate");
    }
    Ok((fillet, hole))
}

/// Audit the fixed 36-block bracket without invoking its constructor.
pub fn audit_bracket(blocks: &[BracketBlock]) -> Result<BracketCertificate, &'static str> {
    validate_blocks(blocks)?;

    let incidence = enumerate_faces(blocks)?;
    if incidence.interfaces.len() != 55 || incidence.exterior.len() != 106 {
        return Err("partition_topology");
    }
    let design_faces = classify_design_faces(blocks, &incidence)?;

    let (fillet_det, hole_det) = profile_determinants()?;
    let mut determinant_lower_bound = f64::INFINITY;
    for (z0, z1) in LAYERS {
        determinant_lower_bound = determinant_lower_bound.min(fillet_det * (z1 - z0) as f64);
        let suffix = format!("_z{z0}_{z1}");
        for block in blocks.iter().filter(|block| block.id.ends_with(&suffix)) {
            if block.degrees == [1, 1, 1] {
                let mut lo = [i64::MAX; 3];
                let mut hi = [i64::MIN; 3];
                for point in &block.controls {
                    for axis in 0..3 {
                        let weight = i128::from(point[3]);
                        let coordinate = i128::from(point[axis]);
                        if coordinate % weight != 0 {
                            return Err("nonintegral_box_vertex");
                        }
                        let physical = i64::try_from(coordinate / weight)
                            .map_err(|_| "box_coordinate_range")?;
                        lo[axis] = lo[axis].min(physical);
                        hi[axis] = hi[axis].max(physical);
                    }
                }
                let volume =
                    (hi[0] - lo[0]) as f64 * (hi[1] - lo[1]) as f64 * (hi[2] - lo[2]) as f64;
                determinant_lower_bound = determinant_lower_bound.min(volume);
            }
        }
    }
    determinant_lower_bound = determinant_lower_bound.min(hole_det * 8.0);
    if !determinant_lower_bound.is_finite() || determinant_lower_bound <= 0.0 {
        return Err("determinant_certificate");
    }

    Ok(BracketCertificate {
        incidence,
        determinant_lower_bound_mm3: determinant_lower_bound,
        radial_error_upper_mm: RADIAL_ERROR_BOUND_MM,
        design_faces,
    })
}
