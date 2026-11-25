//! Exact-profile row of one circular bore per rectangular macrocell.
//!
//! Dimensions, centers, and radii use half-millimeter integer ticks. All
//! homogeneous coefficients are constructed with checked integer arithmetic;
//! the shared E1/E6 radius-four profile is uniformly scaled and translated.

use crate::bracket::{BracketBlock, HOLE_EAST_PROFILE};
use std::convert::TryFrom;

const MAX_CELLS: usize = 16;
const MAX_HALF_WIDTH_X_TICKS: i32 = 160;
const MAX_HALF_WIDTH_Y_TICKS: i32 = 20;
const SCALE: i128 = 8;
const SECTORS: [&str; 4] = ["east", "north", "west", "south"];

/// Specific failure classes for row-bore construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowBoresError {
    InvalidInput,
    NonRepresentableBand,
    DiskOutsideBand { cell: usize },
    DiskTangentToBoundary { cell: usize },
    DisksTangent { left: usize, right: usize },
    DisksIntersect { left: usize, right: usize },
    NumericalRange,
}

/// A constructed row of disjoint bores in equal-width rectangular cells.
#[derive(Debug, Clone, PartialEq)]
pub struct RowBores {
    pub half_width_ticks: [i32; 2],
    pub band_width_ticks: i32,
    pub centers_ticks: Vec<[i32; 2]>,
    pub radii_ticks: Vec<i32>,
    pub common_profile_scale: i64,
    /// Four-sector leaf blocks, ordered left-to-right then E/N/W/S.
    pub blocks: Vec<BracketBlock>,
    /// Ideal target volume for the disjoint circles, in mm³.
    pub ideal_volume_mm3: f64,
}

fn mul(a: i128, b: i128) -> Result<i128, RowBoresError> {
    a.checked_mul(b).ok_or(RowBoresError::NumericalRange)
}

fn add(a: i128, b: i128) -> Result<i128, RowBoresError> {
    a.checked_add(b).ok_or(RowBoresError::NumericalRange)
}

fn to_i64(value: i128) -> Result<i64, RowBoresError> {
    i64::try_from(value).map_err(|_| RowBoresError::NumericalRange)
}

fn rotated(point: [i128; 3], turns: usize) -> Result<[i128; 3], RowBoresError> {
    let [mut x, mut y, w] = point;
    for _ in 0..turns {
        (x, y) = (y.checked_neg().ok_or(RowBoresError::NumericalRange)?, x);
    }
    Ok([x, y, w])
}

fn outer_points_twice(
    sector: usize,
    left_ticks: i128,
    right_ticks: i128,
    half_height_ticks: i128,
) -> [[i128; 2]; 3] {
    let bottom = -half_height_ticks;
    let top = half_height_ticks;
    match sector {
        0 => [
            [2 * right_ticks, 2 * bottom],
            [2 * right_ticks, 0],
            [2 * right_ticks, 2 * top],
        ],
        1 => [
            [2 * right_ticks, 2 * top],
            [left_ticks + right_ticks, 2 * top],
            [2 * left_ticks, 2 * top],
        ],
        2 => [
            [2 * left_ticks, 2 * top],
            [2 * left_ticks, 0],
            [2 * left_ticks, 2 * bottom],
        ],
        _ => [
            [2 * left_ticks, 2 * bottom],
            [left_ticks + right_ticks, 2 * bottom],
            [2 * right_ticks, 2 * bottom],
        ],
    }
}

fn build_block(
    cell: usize,
    sector: usize,
    left_ticks: i128,
    right_ticks: i128,
    half_height_ticks: i128,
    center_ticks: [i128; 2],
    radius_ticks: i128,
) -> Result<BracketBlock, RowBoresError> {
    let outer = outer_points_twice(sector, left_ticks, right_ticks, half_height_ticks);
    let mut controls = Vec::with_capacity(12);
    for angle in 0..3 {
        let outer_weight = HOLE_EAST_PROFILE[2 * angle][2];
        let [inner_x, inner_y, inner_weight] = rotated(HOLE_EAST_PROFILE[2 * angle + 1], sector)?;
        for (radial, weight) in [(0, outer_weight), (1, inner_weight)] {
            let (x, y) = if radial == 0 {
                (
                    mul(mul(outer[angle][0], 2)?, outer_weight)?,
                    mul(mul(outer[angle][1], 2)?, outer_weight)?,
                )
            } else {
                (
                    add(
                        mul(inner_x, radius_ticks)?,
                        mul(mul(center_ticks[0], 4)?, inner_weight)?,
                    )?,
                    add(
                        mul(inner_y, radius_ticks)?,
                        mul(mul(center_ticks[1], 4)?, inner_weight)?,
                    )?,
                )
            };
            let output_weight = mul(weight, SCALE)?;
            for z_mm in [0_i128, 8_i128] {
                controls.push([
                    to_i64(x)?,
                    to_i64(y)?,
                    to_i64(mul(z_mm, output_weight)?)?,
                    to_i64(output_weight)?,
                ]);
            }
        }
    }
    Ok(BracketBlock {
        id: format!("row_bore_cell_{cell:02}_sector_{}_z0_8", SECTORS[sector]),
        degrees: [2, 1, 1],
        controls,
    })
}

fn validate_interactions(
    centers_ticks: &[[i32; 2]],
    radii_ticks: &[i32],
) -> Result<(), RowBoresError> {
    for left in 0..centers_ticks.len() {
        for right in left + 1..centers_ticks.len() {
            let dx = i128::from(centers_ticks[left][0]) - i128::from(centers_ticks[right][0]);
            let dy = i128::from(centers_ticks[left][1]) - i128::from(centers_ticks[right][1]);
            let distance_squared = add(mul(dx, dx)?, mul(dy, dy)?)?;
            let radius_sum = i128::from(radii_ticks[left]) + i128::from(radii_ticks[right]);
            let contact_squared = mul(radius_sum, radius_sum)?;
            if distance_squared == contact_squared {
                return Err(RowBoresError::DisksTangent { left, right });
            }
            if distance_squared < contact_squared {
                return Err(RowBoresError::DisksIntersect { left, right });
            }
        }
    }
    Ok(())
}

/// Construct a row of `n` equal-width rectangular macrocells.
///
/// Each dimension, center coordinate, and radius is specified in half-mm
/// ticks. `centers_ticks[i]` and `radii_ticks[i]` define the bore for band `i`;
/// order is therefore explicitly left-to-right. The rectangle is
/// `[-A,A] × [-B,B] × [0,8]` mm and must satisfy `A<=80`, `B<=10` mm.
pub fn build_row_bores(
    half_width_x_ticks: i32,
    half_width_y_ticks: i32,
    centers_ticks: &[[i32; 2]],
    radii_ticks: &[i32],
) -> Result<RowBores, RowBoresError> {
    let n = centers_ticks.len();
    if !(1..=MAX_CELLS).contains(&n)
        || radii_ticks.len() != n
        || !(1..=MAX_HALF_WIDTH_X_TICKS).contains(&half_width_x_ticks)
        || !(1..=MAX_HALF_WIDTH_Y_TICKS).contains(&half_width_y_ticks)
        || radii_ticks.iter().any(|radius| *radius <= 0)
    {
        return Err(RowBoresError::InvalidInput);
    }
    let total_width_ticks = i64::from(half_width_x_ticks) * 2;
    if total_width_ticks % n as i64 != 0 {
        return Err(RowBoresError::NonRepresentableBand);
    }
    let band_width_ticks = total_width_ticks / n as i64;
    // Give pairwise contact/intersection its own typed result before reporting
    // an individual disk's out-of-band condition.
    validate_interactions(centers_ticks, radii_ticks)?;

    let a = i128::from(half_width_x_ticks);
    let b = i128::from(half_width_y_ticks);
    let width = i128::from(band_width_ticks);
    for index in 0..n {
        let x = i128::from(centers_ticks[index][0]);
        let y = i128::from(centers_ticks[index][1]);
        let r = i128::from(radii_ticks[index]);
        let left = -a + index as i128 * width;
        let right = left + width;
        if x - r == left || x + r == right || y - r == -b || y + r == b {
            return Err(RowBoresError::DiskTangentToBoundary { cell: index });
        }
        if x - r < left || x + r > right || y - r < -b || y + r > b {
            return Err(RowBoresError::DiskOutsideBand { cell: index });
        }
    }

    let mut blocks = Vec::with_capacity(4 * n);
    for cell in 0..n {
        let left = -a + cell as i128 * width;
        let right = left + width;
        let center = centers_ticks[cell].map(i128::from);
        let radius = i128::from(radii_ticks[cell]);
        for sector in 0..4 {
            blocks.push(build_block(cell, sector, left, right, b, center, radius)?);
        }
    }

    let half_width_x_mm = f64::from(half_width_x_ticks) * 0.5;
    let half_width_y_mm = f64::from(half_width_y_ticks) * 0.5;
    let removed_area_mm2: f64 = radii_ticks
        .iter()
        .map(|radius| {
            let r = f64::from(*radius) * 0.5;
            std::f64::consts::PI * r * r
        })
        .sum();
    let ideal_volume_mm3 = 8.0 * (4.0 * half_width_x_mm * half_width_y_mm - removed_area_mm2);
    if !ideal_volume_mm3.is_finite() || ideal_volume_mm3 <= 0.0 {
        return Err(RowBoresError::NumericalRange);
    }

    Ok(RowBores {
        half_width_ticks: [half_width_x_ticks, half_width_y_ticks],
        band_width_ticks: i32::try_from(band_width_ticks)
            .map_err(|_| RowBoresError::NumericalRange)?,
        centers_ticks: centers_ticks.to_vec(),
        radii_ticks: radii_ticks.to_vec(),
        common_profile_scale: i64::try_from(SCALE).map_err(|_| RowBoresError::NumericalRange)?,
        blocks,
        ideal_volume_mm3,
    })
}

#[cfg(test)]
mod tests {
    use super::{build_row_bores, RowBoresError};
    use crate::row_bores_audit::audit_row;

    #[test]
    fn frozen_row_sizes_build_and_pass_the_independent_audit() {
        for n in [1_i32, 2, 4, 8, 16] {
            let half_width_x_ticks = 10 * n;
            let band_width_ticks = 20;
            let centers: Vec<[i32; 2]> = (0..n)
                .map(|i| {
                    [
                        -half_width_x_ticks + band_width_ticks / 2 + i * band_width_ticks,
                        0,
                    ]
                })
                .collect();
            let radii = vec![4; n as usize];
            let model = build_row_bores(half_width_x_ticks, 20, &centers, &radii).unwrap();
            assert_eq!(model.blocks.len(), 4 * n as usize);
            assert_eq!(model.band_width_ticks, band_width_ticks);
            let certificate = audit_row(
                model.half_width_ticks,
                &model.centers_ticks,
                &model.radii_ticks,
                &model.blocks,
            )
            .unwrap();
            assert_eq!(certificate.incidence.interfaces.len(), 5 * n as usize - 1);
            assert_eq!(certificate.incidence.exterior.len(), 14 * n as usize + 2);
            assert!(certificate.determinant_lower_bound_mm3 > 0.0);
            assert!(certificate.radial_error_upper_mm <= 1e-12);
        }
    }

    #[test]
    fn interaction_states_and_band_errors_are_distinguished() {
        for (separation_ticks, expected) in [
            (16, None),
            (12, None),
            (10, None),
            (8, Some(RowBoresError::DisksTangent { left: 0, right: 1 })),
            (4, Some(RowBoresError::DisksIntersect { left: 0, right: 1 })),
        ] {
            let result = build_row_bores(
                20,
                20,
                &[[-separation_ticks / 2, 0], [separation_ticks / 2, 0]],
                &[4, 4],
            );
            if let Some(error) = expected {
                assert_eq!(result.err(), Some(error));
            } else {
                assert!(result.is_ok());
            }
        }
        assert_eq!(
            build_row_bores(10, 20, &[[8, 0]], &[2]).err(),
            Some(RowBoresError::DiskTangentToBoundary { cell: 0 })
        );
        assert_eq!(
            build_row_bores(10, 20, &[[10, 0]], &[2]).err(),
            Some(RowBoresError::DiskOutsideBand { cell: 0 })
        );
        assert_eq!(
            build_row_bores(2, 20, &[[0, 0], [1, 0], [2, 0]], &[1, 1, 1]).err(),
            Some(RowBoresError::NonRepresentableBand)
        );
        assert_eq!(
            build_row_bores(20, 20, &[[0, 0]], &[0]).err(),
            Some(RowBoresError::InvalidInput)
        );
    }

    #[test]
    fn odd_tick_band_midpoints_remain_exact_with_the_common_profile_scale() {
        let model = build_row_bores(159, 20, &[[-80, 0], [80, 0]], &[4, 4]).unwrap();
        assert_eq!(model.band_width_ticks, 159);
        assert_eq!(model.blocks.len(), 8);
        audit_row(
            model.half_width_ticks,
            &model.centers_ticks,
            &model.radii_ticks,
            &model.blocks,
        )
        .unwrap();
    }
}
