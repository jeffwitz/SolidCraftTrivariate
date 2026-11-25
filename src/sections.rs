//! Fixed planar section templates extruded as homogeneous trivariate blocks.

use crate::bracket::{BracketBlock, ELBOW_EAST_PROFILE, HOLE_EAST_PROFILE};
use std::convert::TryFrom;

const SECTORS: [&str; 4] = ["east", "north", "west", "south"];

type Profile = [[i128; 3]; 6];

fn checked_mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("section coefficient exceeds i128")
}

fn push_control(
    controls: &mut Vec<[i64; 4]>,
    x: i128,
    y: i128,
    z: i128,
    weight: i128,
) -> Result<(), &'static str> {
    controls.push([
        i64::try_from(x).map_err(|_| "section control exceeds i64")?,
        i64::try_from(y).map_err(|_| "section control exceeds i64")?,
        i64::try_from(z).map_err(|_| "section control exceeds i64")?,
        i64::try_from(weight).map_err(|_| "section weight exceeds i64")?,
    ]);
    Ok(())
}

fn box_block(id: &str, x: (i64, i64), y: (i64, i64), z: (i64, i64)) -> BracketBlock {
    let mut controls = Vec::with_capacity(8);
    for xx in [x.0, x.1] {
        for yy in [y.0, y.1] {
            for zz in [z.0, z.1] {
                controls.push([xx, yy, zz, 1]);
            }
        }
    }
    BracketBlock {
        id: id.to_string(),
        degrees: [1, 1, 1],
        controls,
    }
}

fn extruded_profile(
    id: String,
    profile: &Profile,
    z0: i64,
    z1: i64,
) -> Result<BracketBlock, &'static str> {
    let mut controls = Vec::with_capacity(12);
    for [x, y, weight] in profile {
        for z in [z0, z1] {
            push_control(
                &mut controls,
                *x,
                *y,
                checked_mul(i128::from(z), *weight)?,
                *weight,
            )?;
        }
    }
    Ok(BracketBlock {
        id,
        degrees: [2, 1, 1],
        controls,
    })
}

fn turned_profile(profile: &Profile, turns: usize) -> Result<Profile, &'static str> {
    let mut result = [[0_i128; 3]; 6];
    for (index, [mut x, mut y, weight]) in profile.iter().copied().enumerate() {
        for _ in 0..turns {
            (x, y) = (y.checked_neg().ok_or("section rotation exceeds i128")?, x);
        }
        result[index] = [x, y, weight];
    }
    Ok(result)
}

fn square_circle_profile() -> Profile {
    HOLE_EAST_PROFILE
}

fn circle_core_profile() -> Result<Profile, &'static str> {
    let mut profile = [[0_i128; 3]; 6];
    for angle in 0..3 {
        let outer_circle = HOLE_EAST_PROFILE[angle * 2 + 1];
        let square = HOLE_EAST_PROFILE[angle * 2];
        if square[0] % 10 != 0 || square[1] % 10 != 0 {
            return Err("square profile cannot be scaled to the unit core exactly");
        }
        profile[angle * 2] = outer_circle;
        profile[angle * 2 + 1] = [square[0] / 10, square[1] / 10, square[2]];
    }
    Ok(profile)
}

fn annulus_profile() -> Result<Profile, &'static str> {
    let mut profile = [[0_i128; 3]; 6];
    for angle in 0..3 {
        let radius_four = HOLE_EAST_PROFILE[angle * 2 + 1];
        if radius_four[0] % 2 != 0 || radius_four[1] % 2 != 0 {
            return Err("radius profile cannot be halved exactly");
        }
        profile[angle * 2] = radius_four;
        profile[angle * 2 + 1] = [radius_four[0] / 2, radius_four[1] / 2, radius_four[2]];
    }
    Ok(profile)
}

fn corner_profile(side: &str) -> Profile {
    let mut profile = ELBOW_EAST_PROFILE;
    if side == "north" {
        let source = profile;
        for i in 0..3 {
            for j in 0..2 {
                let [x, y, weight] = source[(2 - i) * 2 + j];
                profile[i * 2 + j] = [y, x, weight];
            }
        }
    }
    profile
}

fn append_sectors(
    blocks: &mut Vec<BracketBlock>,
    id_prefix: &str,
    profile: &Profile,
    z0: i64,
    z1: i64,
) -> Result<(), &'static str> {
    for (turns, sector) in SECTORS.iter().enumerate() {
        let rotated = turned_profile(profile, turns)?;
        blocks.push(extruded_profile(
            format!("{id_prefix}_{sector}_z{z0}_{z1}"),
            &rotated,
            z0,
            z1,
        )?);
    }
    Ok(())
}

/// Construct one of the six fixed section templates in the J1 geometry corpus.
///
/// Curved controls reuse the frozen E1/E3 binary64 profiles exposed as integer
/// dyadic nets by the bracket module. Every block is a local homogeneous net;
/// no model publication or geometric certificate is implied here.
pub fn section_blocks(family: &str) -> Result<Vec<BracketBlock>, &'static str> {
    let mut blocks = Vec::new();
    match family {
        "reference_rectangle_v1" => {
            blocks.push(box_block("rectangle_z0_8", (-10, 10), (-10, 10), (0, 8)));
        }
        "core_disk_v1" => {
            blocks.push(box_block("disk_core_z0_8", (-1, 1), (-1, 1), (0, 8)));
            append_sectors(&mut blocks, "disk_sector", &circle_core_profile()?, 0, 8)?;
        }
        "full_annulus_v1" => {
            append_sectors(&mut blocks, "annulus_sector", &annulus_profile()?, 0, 8)?;
        }
        "pierced_square_v1" => {
            append_sectors(
                &mut blocks,
                "pierced_sector",
                &square_circle_profile(),
                0,
                8,
            )?;
        }
        "rounded_corner_v1" => {
            for side in ["east", "north"] {
                blocks.push(extruded_profile(
                    format!("corner_{side}_z0_8"),
                    &corner_profile(side),
                    0,
                    8,
                )?);
            }
        }
        "layered_pocket_v1" => {
            blocks.push(box_block("pocket_core_z0_4", (-1, 1), (-1, 1), (0, 4)));
            append_sectors(
                &mut blocks,
                "pocket_disk_sector",
                &circle_core_profile()?,
                0,
                4,
            )?;
            append_sectors(
                &mut blocks,
                "pocket_outer_sector",
                &square_circle_profile(),
                0,
                4,
            )?;
            append_sectors(
                &mut blocks,
                "pocket_outer_sector",
                &square_circle_profile(),
                4,
                8,
            )?;
        }
        _ => return Err("unsupported section family"),
    }
    Ok(blocks)
}
