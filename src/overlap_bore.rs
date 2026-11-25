//! Fixed two-circle overlap bore in a square plate.
//!
//! The family is intentionally limited to two radius-2 mm circles centered
//! at (+/-1, 0) mm in [-10, 10]^2 x [0, 8] mm. The conic arcs use a frozen
//! rational approximation of sqrt(3); all emitted homogeneous coefficients
//! are exact integers and each radial generator shares its endpoint weights.

use crate::bracket::BracketBlock;
use std::convert::TryFrom;

const S_NUM: i128 = 1_732_050_807_568_877;
const S_DEN: i128 = 1_000_000_000_000_000;
const SCALE: i128 = S_DEN;
const HEIGHT_MM: i128 = 8;
const SECTORS: [&str; 2] = ["right", "left"];
// t intervals [-1,-1/3], [-1/3,0], [0,1/3], [1/3,1], as (p,r)/3.
const INTERVALS: [(i128, i128); 4] = [(-3, 2), (-1, 1), (0, 1), (1, 2)];

fn add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("overlap bore coefficient overflow")
}

fn sub(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_sub(b).ok_or("overlap bore coefficient overflow")
}

fn mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("overlap bore coefficient overflow")
}

fn output(value: i128) -> Result<i64, &'static str> {
    i64::try_from(value).map_err(|_| "overlap bore control exceeds i64")
}

/// Exact homogeneous controls for one restricted rational conic arc.
///
/// Each triple is `(wx, wy, w)`, and the local curve parameter runs from zero
/// to one over the corresponding fixed interval in the global conic parameter.
fn arc_controls(interval: usize) -> Result<[[i128; 3]; 3], &'static str> {
    let (p, r) = INTERVALS[interval];
    let l = SCALE;

    // Substitute t=(p+r*u)/3 into
    //   X = 1 + 2(1-3t^2)/(1+3t^2),
    //   Y = 4*s*t/(1+3t^2),
    //   W = 1.
    // Convert each power-basis quadratic to Bernstein controls. The middle
    // control is c0+c1/2; c1 is always even for these integer parameters.
    let p2 = mul(p, p)?;
    let pr2 = mul(2, mul(p, r)?)?;
    let r2 = mul(r, r)?;
    let x_coeff = [mul(l, sub(9, p2)?)?, mul(l, -pr2)?, mul(l, -r2)?];
    let y_coeff = [mul(4, mul(S_NUM, p)?)?, mul(4, mul(S_NUM, r)?)?, 0];
    let w_coeff = [mul(l, add(3, p2)?)?, mul(l, pr2)?, mul(l, r2)?];
    let coefficients = [x_coeff, y_coeff, w_coeff];
    let mut controls = [[0_i128; 3]; 3];
    for component in 0..3 {
        let [c0, c1, c2] = coefficients[component];
        if c1 % 2 != 0 {
            return Err("nonintegral overlap bore Bernstein control");
        }
        controls[0][component] = c0;
        controls[1][component] = add(c0, c1 / 2)?;
        controls[2][component] = add(add(c0, c1)?, c2)?;
    }
    Ok(controls)
}

fn outer_polyline() -> [[i128; 2]; 5] {
    [[0, -10], [10, -10], [10, 0], [10, 10], [0, 10]]
}

fn rotated_point(point: [i128; 2], half: usize) -> Result<[i128; 2], &'static str> {
    if half == 0 {
        Ok(point)
    } else {
        Ok([
            point[0]
                .checked_neg()
                .ok_or("overlap bore coefficient overflow")?,
            point[1]
                .checked_neg()
                .ok_or("overlap bore coefficient overflow")?,
        ])
    }
}

fn build_sector(half: usize, sector: usize) -> Result<BracketBlock, &'static str> {
    let arc = arc_controls(sector)?;
    let path = outer_polyline();
    let [path_start, path_end] = [path[sector], path[sector + 1]];
    let outer = [
        path_start,
        [
            add(path_start[0], path_end[0])? / 2,
            add(path_start[1], path_end[1])? / 2,
        ],
        path_end,
    ];
    let mut controls = Vec::with_capacity(12);
    for angle in 0..3 {
        let [x_inner, y_inner, weight] = arc[angle];
        let [outer_x, outer_y] = rotated_point(outer[angle], half)?;
        let (x_inner, y_inner) = if half == 0 {
            (x_inner, y_inner)
        } else {
            (
                x_inner
                    .checked_neg()
                    .ok_or("overlap bore coefficient overflow")?,
                y_inner
                    .checked_neg()
                    .ok_or("overlap bore coefficient overflow")?,
            )
        };
        let radial = [
            [mul(outer_x, weight)?, mul(outer_y, weight)?, weight],
            [x_inner, y_inner, weight],
        ];
        for [x, y, w] in radial {
            for z in [0_i128, HEIGHT_MM] {
                controls.push([output(x)?, output(y)?, output(mul(z, w)?)?, output(w)?]);
            }
        }
    }
    Ok(BracketBlock {
        id: format!("overlap_bore_{}_sector_{sector:02}_z0_8", SECTORS[half]),
        degrees: [2, 1, 1],
        controls,
    })
}

/// Build the fixed two-overlapping-circle bore as eight rational sectors.
///
/// The right half uses increasing conic parameter from the lower intersection
/// to the upper intersection. The left half is its direct pi rotation. For
/// each half, sector index 00..03 follows the outer polyline from the lower
/// square seam, across the bottom/right/top edges, to the upper square seam.
/// Controls are flattened in `(angle, radial, z)` order, with z fastest;
/// radial index 0 is the square-side curve and radial index 1 the bore arc.
///
/// The approximation `S_NUM/S_DEN` substitutes for sqrt(3); construction of
/// the represented rational model is exact, while the tiny boundary mismatch
/// to the ideal intersection is certified independently by the audit.
pub fn overlap_blocks() -> Result<Vec<BracketBlock>, &'static str> {
    if S_DEN != SCALE || S_NUM <= 0 || S_NUM >= 2 * SCALE {
        return Err("invalid frozen overlap bore constants");
    }
    let mut blocks = Vec::with_capacity(8);
    for half in 0..2 {
        for sector in 0..4 {
            blocks.push(build_sector(half, sector)?);
        }
    }
    if blocks.len() != 8
        || blocks.iter().any(|block| {
            block.degrees != [2, 1, 1]
                || block.controls.len() != 12
                || block.controls.iter().any(|control| control[3] <= 0)
        })
    {
        return Err("invalid overlap bore control net");
    }
    Ok(blocks)
}
