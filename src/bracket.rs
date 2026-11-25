//! Fixed 36-block filleted bracket, ported from the frozen E3/E5 profiles.

use std::convert::TryFrom;

/// A persisted homogeneous block with one positive integer scale per block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BracketBlock {
    pub id: String,
    pub degrees: [usize; 3],
    pub controls: Vec<[i64; 4]>,
}

const LAYERS: [(i64, i64); 3] = [(0, 10), (10, 30), (30, 40)];
const SECTORS: [&str; 4] = ["east", "north", "west", "south"];

pub(crate) type Profile = [[i128; 3]; 6];

// Frozen E3 profile derived from the radius-6 arc at (14,14) and the straight
// base in trivariate/elbow.py. Each triple is (wx, wy, w), scaled by 2^52.
const ELBOW_PROFILE_SCALE: i128 = 4_503_599_627_370_496;
pub(crate) const ELBOW_EAST_PROFILE: Profile = [
    [0, 0, ELBOW_PROFILE_SCALE],
    [
        43_943_239_765_609_632,
        43_943_239_765_609_632,
        ELBOW_PROFILE_SCALE,
    ],
    [31_525_197_391_593_472, 0, ELBOW_PROFILE_SCALE],
    [
        47_910_251_476_541_136,
        33_286_268_146_824_472,
        4_160_783_518_353_059,
    ],
    [63_050_394_783_186_944, 0, ELBOW_PROFILE_SCALE],
    [
        63_050_394_783_186_944,
        36_028_797_018_963_968,
        ELBOW_PROFILE_SCALE,
    ],
];

// Frozen E1 east-sector profile from trivariate/pierced.py. Triples are
// homogeneous angular/radial controls (wx, wy, w), scaled by 2^53.
const HOLE_PROFILE_SCALE: i128 = 9_007_199_254_740_992;
pub(crate) const HOLE_EAST_PROFILE: Profile = [
    [
        90_071_992_547_409_920,
        -90_071_992_547_409_920,
        HOLE_PROFILE_SCALE,
    ],
    [
        25_476_206_690_103_092,
        -25_476_206_690_103_092,
        HOLE_PROFILE_SCALE,
    ],
    [90_071_992_547_409_920, 0, HOLE_PROFILE_SCALE],
    [36_028_797_018_963_976, 0, 6_369_051_672_525_773],
    [
        90_071_992_547_409_920,
        90_071_992_547_409_920,
        HOLE_PROFILE_SCALE,
    ],
    [
        25_476_206_690_103_092,
        25_476_206_690_103_092,
        HOLE_PROFILE_SCALE,
    ],
];

fn checked_add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("transformed control exceeds i128")
}

fn checked_sub(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_sub(b).ok_or("transformed control exceeds i128")
}

fn checked_mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("transformed control exceeds i128")
}

fn push_control(
    controls: &mut Vec<[i64; 4]>,
    x: i128,
    y: i128,
    z: i128,
    weight: i128,
) -> Result<(), &'static str> {
    controls.push([
        i64::try_from(x).map_err(|_| "control exceeds i64")?,
        i64::try_from(y).map_err(|_| "control exceeds i64")?,
        i64::try_from(z).map_err(|_| "control exceeds i64")?,
        i64::try_from(weight).map_err(|_| "weight exceeds i64")?,
    ]);
    Ok(())
}

fn integer_box(id: String, x: (i64, i64), y: (i64, i64), z: (i64, i64)) -> BracketBlock {
    let mut controls = Vec::with_capacity(8);
    for xx in [x.0, x.1] {
        for yy in [y.0, y.1] {
            for zz in [z.0, z.1] {
                controls.push([xx, yy, zz, 1]);
            }
        }
    }
    BracketBlock {
        id,
        degrees: [1, 1, 1],
        controls,
    }
}

fn elbow_block(
    side: &str,
    z0: i64,
    z1: i64,
    source: &Profile,
) -> Result<BracketBlock, &'static str> {
    let mut controls = Vec::with_capacity(12);
    for i in 0..3 {
        let source_i = if side == "east" { i } else { 2 - i };
        for j in 0..2 {
            let [source_x, source_y, weight] = source[source_i * 2 + j];
            let (x, y) = if side == "east" {
                (source_x, source_y)
            } else {
                (source_y, source_x)
            };
            for z in [z0, z1] {
                push_control(
                    &mut controls,
                    x,
                    y,
                    checked_mul(i128::from(z), weight)?,
                    weight,
                )?;
            }
        }
    }
    Ok(BracketBlock {
        id: format!("elbow_{side}_z{z0}_{z1}"),
        degrees: [2, 1, 1],
        controls,
    })
}

fn rotate_quarter_turn(
    mut x: i128,
    mut y: i128,
    turns: usize,
) -> Result<(i128, i128), &'static str> {
    for _ in 0..turns {
        let next_x = y.checked_neg().ok_or("rotated control exceeds i128")?;
        y = x;
        x = next_x;
    }
    Ok((x, y))
}

fn hole_block(
    leg: &str,
    sector: &str,
    turns: usize,
    z0: &Profile,
) -> Result<BracketBlock, &'static str> {
    let mut controls = Vec::with_capacity(12);
    for point in z0 {
        let (local_x, local_y, weight) = (point[0], point[1], point[2]);
        let (a, b) = rotate_quarter_turn(local_x, local_y, turns)?;
        for c in [0_i128, 8_i128] {
            let (x, y, z) = if leg == "x" {
                (
                    checked_add(checked_mul(40, weight)?, a)?,
                    checked_sub(checked_mul(8, weight)?, checked_mul(c, weight)?)?,
                    checked_add(checked_mul(20, weight)?, b)?,
                )
            } else {
                (
                    checked_mul(c, weight)?,
                    checked_add(checked_mul(32, weight)?, a)?,
                    checked_add(checked_mul(20, weight)?, b)?,
                )
            };
            push_control(&mut controls, x, y, z, weight)?;
        }
    }
    Ok(BracketBlock {
        id: format!("hole_{leg}_sector_{sector}_z10_30"),
        degrees: [2, 1, 1],
        controls,
    })
}

fn push_layer_boxes(blocks: &mut Vec<BracketBlock>, z0: i64, z1: i64, middle: bool) {
    let suffix = format!("z{z0}_{z1}");
    blocks.push(integer_box(
        format!("connector_x_14_16_{suffix}"),
        (14, 16),
        (0, 8),
        (z0, z1),
    ));
    blocks.push(integer_box(
        format!("connector_y_14_16_{suffix}"),
        (0, 8),
        (14, 16),
        (z0, z1),
    ));

    let x_intervals: &[(i64, i64)] = if middle {
        &[(16, 30), (50, 60)]
    } else {
        &[(16, 30), (30, 50), (50, 60)]
    };
    for &(x0, x1) in x_intervals {
        blocks.push(integer_box(
            format!("arm_x_{x0}_{x1}_{suffix}"),
            (x0, x1),
            (0, 8),
            (z0, z1),
        ));
    }

    let y_intervals: &[(i64, i64)] = if middle {
        &[(16, 22), (42, 50)]
    } else {
        &[(16, 22), (22, 42), (42, 50)]
    };
    for &(y0, y1) in y_intervals {
        blocks.push(integer_box(
            format!("arm_y_{y0}_{y1}_{suffix}"),
            (0, 8),
            (y0, y1),
            (z0, z1),
        ));
    }
}

/// Build the fixed 60 x 50 x 40 mm bracket with two radius-4 through holes.
///
/// Each block uses one positive common dyadic scale, so every homogeneous
/// coefficient is exactly integral while preserving the prototype geometry.
pub fn bracket_blocks() -> Result<Vec<BracketBlock>, &'static str> {
    let mut blocks = Vec::with_capacity(36);

    for (z0, z1) in LAYERS {
        blocks.push(elbow_block("east", z0, z1, &ELBOW_EAST_PROFILE)?);
        blocks.push(elbow_block("north", z0, z1, &ELBOW_EAST_PROFILE)?);
        push_layer_boxes(&mut blocks, z0, z1, (z0, z1) == (10, 30));
    }

    for (leg, turns) in [
        ("x", 0),
        ("x", 1),
        ("x", 2),
        ("x", 3),
        ("y", 0),
        ("y", 1),
        ("y", 2),
        ("y", 3),
    ] {
        blocks.push(hole_block(leg, SECTORS[turns], turns, &HOLE_EAST_PROFILE)?);
    }

    if blocks.len() != 36 {
        return Err("bracket block inventory is inconsistent");
    }
    Ok(blocks)
}
