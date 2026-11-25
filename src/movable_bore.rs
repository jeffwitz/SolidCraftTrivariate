//! Movable circular bore in a fixed square plate, with analytic control jets.
//!
//! Inputs are integer half-millimeter ticks. Homogeneous controls use one
//! common scale (eight times the frozen E1 profile scale) for every block and
//! every parameter value. The returned geometry is a candidate net; no
//! continuous-domain certificate or rounding enclosure is implied.

use crate::bracket::{BracketBlock, HOLE_EAST_PROFILE};
use std::convert::TryFrom;

const SECTORS: [&str; 4] = ["east", "north", "west", "south"];
const HEIGHT_MM: i128 = 8;

/// Why a movable-bore request was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MovableBoreError {
    /// A center or radius lies outside the frozen parameter domain.
    ParameterOutOfDomain,
    /// An exact intermediate or persisted coefficient exceeds its integer type.
    NumericalRange,
}

/// A block and the physical-parameter derivatives of its homogeneous controls.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MovableBoreBlock {
    pub block: BracketBlock,
    /// Derivative of each `(wx, wy, wz, w)` control per physical mm of `cx`.
    pub dcontrols_dcx: Vec<[i64; 4]>,
    /// Derivative per physical mm of `cy`.
    pub dcontrols_dcy: Vec<[i64; 4]>,
    /// Derivative per physical mm of `r`.
    pub dcontrols_dr: Vec<[i64; 4]>,
}

/// Candidate geometry and exact analytic target-volume sensitivities.
#[derive(Debug, Clone, PartialEq)]
pub struct MovableBore {
    /// Center and radius in signed half-millimeter ticks.
    pub parameter_ticks: [i32; 3],
    pub blocks: Vec<MovableBoreBlock>,
    /// Ideal target volume in mm³. The materialized boundary is a rational
    /// approximation to the ideal circular bore.
    pub volume_mm3: f64,
    /// Ideal target derivatives with respect to physical `(cx, cy, r)` in mm.
    pub dvolume_dparameters: [f64; 3],
}

fn checked_mul(a: i128, b: i128) -> Result<i128, MovableBoreError> {
    a.checked_mul(b).ok_or(MovableBoreError::NumericalRange)
}

fn checked_add(a: i128, b: i128) -> Result<i128, MovableBoreError> {
    a.checked_add(b).ok_or(MovableBoreError::NumericalRange)
}

fn to_i64(value: i128) -> Result<i64, MovableBoreError> {
    i64::try_from(value).map_err(|_| MovableBoreError::NumericalRange)
}

fn rotate_profile(
    profile: &[[i128; 3]; 6],
    quarter_turns: usize,
) -> Result<[[i128; 3]; 6], MovableBoreError> {
    let mut result = [[0_i128; 3]; 6];
    for (index, [mut x, mut y, weight]) in profile.iter().copied().enumerate() {
        for _ in 0..quarter_turns {
            (x, y) = (y.checked_neg().ok_or(MovableBoreError::NumericalRange)?, x);
        }
        result[index] = [x, y, weight];
    }
    Ok(result)
}

fn build_sector(
    sector: &str,
    turns: usize,
    cx_ticks: i32,
    cy_ticks: i32,
    r_ticks: i32,
) -> Result<MovableBoreBlock, MovableBoreError> {
    let profile = rotate_profile(&HOLE_EAST_PROFILE, turns)?;
    let mut controls = Vec::with_capacity(12);
    let mut dcontrols_dcx = Vec::with_capacity(12);
    let mut dcontrols_dcy = Vec::with_capacity(12);
    let mut dcontrols_dr = Vec::with_capacity(12);

    for (index, [profile_x, profile_y, profile_weight]) in profile.into_iter().enumerate() {
        let radial = index % 2;
        for z in [0_i128, HEIGHT_MM] {
            let is_arc = radial == 1;
            // The square-side row is immutable. The inner circular row is
            // scaled about the bore center and translated in physical mm.
            let x = if is_arc {
                checked_add(
                    checked_mul(profile_x, i128::from(r_ticks))?,
                    checked_mul(checked_mul(i128::from(cx_ticks), 4)?, profile_weight)?,
                )?
            } else {
                checked_mul(profile_x, 8)?
            };
            let y = if is_arc {
                checked_add(
                    checked_mul(profile_y, i128::from(r_ticks))?,
                    checked_mul(checked_mul(i128::from(cy_ticks), 4)?, profile_weight)?,
                )?
            } else {
                checked_mul(profile_y, 8)?
            };
            let weight = checked_mul(profile_weight, 8)?;
            let zw = checked_mul(z, weight)?;
            controls.push([to_i64(x)?, to_i64(y)?, to_i64(zw)?, to_i64(weight)?]);

            // Ticks are 0.5 mm: d(radius scale)/dr = 2. Translation
            // derivatives equal the homogeneous weight at the common scale.
            dcontrols_dcx.push([to_i64(if is_arc { weight } else { 0 })?, 0, 0, 0]);
            dcontrols_dcy.push([0, to_i64(if is_arc { weight } else { 0 })?, 0, 0]);
            dcontrols_dr.push([
                to_i64(if is_arc {
                    checked_mul(profile_x, 2)?
                } else {
                    0
                })?,
                to_i64(if is_arc {
                    checked_mul(profile_y, 2)?
                } else {
                    0
                })?,
                0,
                0,
            ]);
        }
    }

    Ok(MovableBoreBlock {
        block: BracketBlock {
            id: format!("movable_bore_sector_{sector}_z0_8"),
            degrees: [2, 1, 1],
            controls,
        },
        dcontrols_dcx,
        dcontrols_dcy,
        dcontrols_dr,
    })
}

/// Construct four exact-profile sectors for a movable bore.
///
/// `cx_ticks`, `cy_ticks`, and `r_ticks` are half-millimeter units. Centers
/// satisfy `|cx|,|cy| <= 2 mm`; radii satisfy `2 <= r <= 5 mm`. The square
/// remains `[-10,10]²` and the plate occupies `z=[0,8]` mm. For `(0,0,8)`
/// (center zero, radius four mm), every persisted coefficient is exactly
/// eight times the corresponding `pierced_square_v1` coefficient.
pub fn build_movable_bore(
    cx_ticks: i32,
    cy_ticks: i32,
    r_ticks: i32,
) -> Result<MovableBore, MovableBoreError> {
    if !(-4..=4).contains(&cx_ticks)
        || !(-4..=4).contains(&cy_ticks)
        || !(4..=10).contains(&r_ticks)
    {
        return Err(MovableBoreError::ParameterOutOfDomain);
    }

    let mut blocks = Vec::with_capacity(4);
    for (turns, sector) in SECTORS.iter().enumerate() {
        blocks.push(build_sector(sector, turns, cx_ticks, cy_ticks, r_ticks)?);
    }

    let radius_mm = f64::from(r_ticks) * 0.5;
    let volume_mm3 = 8.0 * (400.0 - std::f64::consts::PI * radius_mm * radius_mm);
    let dvolume_dr = -16.0 * std::f64::consts::PI * radius_mm;
    if !volume_mm3.is_finite() || !dvolume_dr.is_finite() {
        return Err(MovableBoreError::NumericalRange);
    }

    Ok(MovableBore {
        parameter_ticks: [cx_ticks, cy_ticks, r_ticks],
        blocks,
        volume_mm3,
        dvolume_dparameters: [0.0, 0.0, dvolume_dr],
    })
}
