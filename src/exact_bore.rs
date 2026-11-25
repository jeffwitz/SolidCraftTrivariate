//! Exact integer construction of a square plate with one circular bore.
//!
//! Inputs and spatial coordinates are integer micrometers. Each quadrant's
//! rational circle arc is split exactly at its Bezier parameter `u=1/2`; the
//! two outer charts follow the adjacent straight square edges. No floating
//! point operation participates in generating a control coefficient.

use crate::bracket::BracketBlock;
use std::convert::TryFrom;

const HALF_WIDTH_UM: i128 = 10_000;
const HEIGHT_MM: i128 = 8;
const COMMON_SCALE: i128 = 8;
const SECTORS: [&str; 4] = ["east", "north", "west", "south"];

/// Why an exact-bore request was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExactBoreError {
    /// The bore is empty, outside the plate, or tangent to the square boundary.
    GeometricDomain,
    /// An exact integer intermediate or persisted coefficient overflowed.
    NumericalRange,
}

/// An exact bore construction. Controls evaluate in millimeters; X/Y numerators
/// are in micrometers and weights include the exact factor of 1000.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExactBore {
    pub center_um: [i32; 2],
    pub radius_um: i32,
    pub common_scale: i64,
    /// Eight charts in deterministic E/N/W/S then half 0/1 order.
    pub blocks: Vec<BracketBlock>,
}

fn add(a: i128, b: i128) -> Result<i128, ExactBoreError> {
    a.checked_add(b).ok_or(ExactBoreError::NumericalRange)
}

fn sub(a: i128, b: i128) -> Result<i128, ExactBoreError> {
    a.checked_sub(b).ok_or(ExactBoreError::NumericalRange)
}

fn mul(a: i128, b: i128) -> Result<i128, ExactBoreError> {
    a.checked_mul(b).ok_or(ExactBoreError::NumericalRange)
}

fn as_i64(value: i128) -> Result<i64, ExactBoreError> {
    i64::try_from(value).map_err(|_| ExactBoreError::NumericalRange)
}

fn rotate_homogeneous(
    point: [i128; 3],
    center: [i128; 2],
    quarter_turns: usize,
) -> Result<[i128; 3], ExactBoreError> {
    let [x, y, weight] = point;
    let mut dx = sub(x, mul(center[0], weight)?)?;
    let mut dy = sub(y, mul(center[1], weight)?)?;
    for _ in 0..quarter_turns {
        (dx, dy) = (dy.checked_neg().ok_or(ExactBoreError::NumericalRange)?, dx);
    }
    Ok([
        add(mul(center[0], weight)?, dx)?,
        add(mul(center[1], weight)?, dy)?,
        weight,
    ])
}

fn split_at_half(points: [[i128; 3]; 3]) -> Result<[[[i128; 3]; 3]; 2], ExactBoreError> {
    let [p0, p1, p2] = points;
    let midpoint = [
        add(add(p0[0], mul(2, p1[0])?)?, p2[0])?,
        add(add(p0[1], mul(2, p1[1])?)?, p2[1])?,
        add(add(p0[2], mul(2, p1[2])?)?, p2[2])?,
    ];
    let left_middle = [
        mul(2, add(p0[0], p1[0])?)?,
        mul(2, add(p0[1], p1[1])?)?,
        mul(2, add(p0[2], p1[2])?)?,
    ];
    let right_middle = [
        mul(2, add(p1[0], p2[0])?)?,
        mul(2, add(p1[1], p2[1])?)?,
        mul(2, add(p1[2], p2[2])?)?,
    ];
    let left = [
        [mul(4, p0[0])?, mul(4, p0[1])?, mul(4, p0[2])?],
        left_middle,
        midpoint,
    ];
    let right = [
        midpoint,
        right_middle,
        [mul(4, p2[0])?, mul(4, p2[1])?, mul(4, p2[2])?],
    ];
    Ok([left, right])
}

fn circle_quadrant(
    center: [i128; 2],
    radius: i128,
    quarter_turns: usize,
) -> Result<[[i128; 3]; 3], ExactBoreError> {
    // Bernstein controls for ((1-u²), 2u, (1+u²)), u∈[0,1].
    let base = [
        [add(center[0], radius)?, center[1], 1],
        [add(center[0], radius)?, add(center[1], radius)?, 1],
        [
            mul(center[0], 2)?,
            add(mul(center[1], 2)?, mul(radius, 2)?)?,
            2,
        ],
    ];
    let mut result = [[0_i128; 3]; 3];
    for (index, point) in base.into_iter().enumerate() {
        result[index] = rotate_homogeneous(point, center, quarter_turns)?;
    }
    Ok(result)
}

fn outer_segments(center: [i128; 2], sector: usize) -> [[i128; 2]; 3] {
    let c = center;
    match sector {
        0 => [
            [HALF_WIDTH_UM, c[1]],
            [HALF_WIDTH_UM, HALF_WIDTH_UM],
            [c[0], HALF_WIDTH_UM],
        ],
        1 => [
            [c[0], HALF_WIDTH_UM],
            [-HALF_WIDTH_UM, HALF_WIDTH_UM],
            [-HALF_WIDTH_UM, c[1]],
        ],
        2 => [
            [-HALF_WIDTH_UM, c[1]],
            [-HALF_WIDTH_UM, -HALF_WIDTH_UM],
            [c[0], -HALF_WIDTH_UM],
        ],
        _ => [
            [c[0], -HALF_WIDTH_UM],
            [HALF_WIDTH_UM, -HALF_WIDTH_UM],
            [HALF_WIDTH_UM, c[1]],
        ],
    }
}

fn build_chart(
    center: [i128; 2],
    radius: i128,
    sector: usize,
    half: usize,
) -> Result<BracketBlock, ExactBoreError> {
    let arc = circle_quadrant(center, radius, sector)?;
    let arc_halves = split_at_half(arc)?;
    let arc_half = arc_halves[half];
    let outer = outer_segments(center, sector);
    let [outer_start, corner, outer_end] = outer;
    let (outer_a, outer_b) = if half == 0 {
        (outer_start, corner)
    } else {
        (corner, outer_end)
    };

    let mut controls = Vec::with_capacity(12);
    for (u, arc_control) in arc_half.into_iter().enumerate() {
        let weight = mul(arc_control[2], 2)?;
        if weight <= 0 || weight % 2 != 0 {
            return Err(ExactBoreError::NumericalRange);
        }
        let outer_xy = match u {
            0 => [mul(outer_a[0], weight)?, mul(outer_a[1], weight)?],
            1 => [
                mul(add(outer_a[0], outer_b[0])?, weight / 2)?,
                mul(add(outer_a[1], outer_b[1])?, weight / 2)?,
            ],
            _ => [mul(outer_b[0], weight)?, mul(outer_b[1], weight)?],
        };
        let inner_xy = [mul(arc_control[0], 2)?, mul(arc_control[1], 2)?];
        for [x, y] in [outer_xy, inner_xy] {
            let output_weight = mul(weight, 1_000)?;
            for z_mm in [0_i128, HEIGHT_MM] {
                controls.push([
                    as_i64(x)?,
                    as_i64(y)?,
                    as_i64(mul(z_mm, output_weight)?)?,
                    as_i64(output_weight)?,
                ]);
            }
        }
    }

    Ok(BracketBlock {
        id: format!("exact_bore_sector_{}_half_{half}_z0_8000", SECTORS[sector]),
        degrees: [2, 1, 1],
        controls,
    })
}

/// Construct an exact single-bore plate from integer micrometer coordinates.
///
/// Admissibility is strict: `r > 0`, `|cx| + r < 10000`, and
/// `|cy| + r < 10000`. Tangency to the square is rejected. Output uses 8 as
/// one common positive homogeneous scale: four for the exact half-parameter
/// restriction and two for the possible half-micrometer outer midpoint.
pub fn build_exact_bore(
    cx_um: i32,
    cy_um: i32,
    radius_um: i32,
) -> Result<ExactBore, ExactBoreError> {
    let cx = i128::from(cx_um);
    let cy = i128::from(cy_um);
    let radius = i128::from(radius_um);
    if radius <= 0 || cx.abs() + radius >= HALF_WIDTH_UM || cy.abs() + radius >= HALF_WIDTH_UM {
        return Err(ExactBoreError::GeometricDomain);
    }

    let center = [cx, cy];
    let mut blocks = Vec::with_capacity(8);
    for sector in 0..4 {
        for half in 0..2 {
            blocks.push(build_chart(center, radius, sector, half)?);
        }
    }
    if blocks.iter().any(|block| {
        block.degrees != [2, 1, 1]
            || block.controls.len() != 12
            || block.controls.iter().any(|point| point[3] <= 0)
    }) {
        return Err(ExactBoreError::NumericalRange);
    }

    Ok(ExactBore {
        center_um: [cx_um, cy_um],
        radius_um,
        common_scale: i64::try_from(COMMON_SCALE).map_err(|_| ExactBoreError::NumericalRange)?,
        blocks,
    })
}

#[cfg(test)]
mod tests {
    use super::{build_exact_bore, ExactBoreError};
    use crate::exact_bore_audit::audit_exact_bore;

    #[test]
    fn exact_bore_builds_eight_mm_charts_and_passes_independent_audit() {
        for (center, radius) in [
            ([0, 0], 4_000),
            ([5_999, 0], 4_000),
            ([-5_999, 0], 4_000),
            ([0, 5_999], 4_000),
            ([0, -5_999], 4_000),
        ] {
            let model = build_exact_bore(center[0], center[1], radius).unwrap();
            assert_eq!(model.common_scale, 8);
            assert_eq!(model.blocks.len(), 8);
            for block in &model.blocks {
                assert_eq!(block.degrees, [2, 1, 1]);
                assert_eq!(block.controls.len(), 12);
                for control in &block.controls {
                    assert!(control[3] > 0);
                    assert_eq!(control[3] % 1_000, 0);
                    assert!(control[0].abs() < 10_000 * control[3]);
                    assert!(control[1].abs() < 10_000 * control[3]);
                    assert!(control[2] == 0 || control[2] == 8 * control[3]);
                }
            }
            let certificate = audit_exact_bore(center, radius, &model.blocks).unwrap();
            assert_eq!(certificate.incidence.interfaces.len(), 8);
            assert_eq!(certificate.incidence.exterior.len(), 32);
            assert!(certificate.determinant_lower_bound_mm3.is_finite());
            assert!(certificate.determinant_lower_bound_mm3 > 0.0);
            assert_eq!(certificate.radial_error_upper_mm, 0.0);
        }
    }

    #[test]
    fn exact_bore_rejects_empty_invalid_and_tangent_inputs() {
        for args in [
            (0, 0, 0),
            (0, 0, -1),
            (6_000, 0, 4_000),
            (-6_000, 0, 4_000),
            (0, 6_000, 4_000),
            (0, -6_000, 4_000),
            (i32::MIN, 0, 1),
        ] {
            assert_eq!(
                build_exact_bore(args.0, args.1, args.2).err(),
                Some(ExactBoreError::GeometricDomain)
            );
        }
        // Strict inequality leaves a one-micrometer ligament and is admissible.
        assert!(build_exact_bore(0, 0, 9_999).is_ok());
    }
}
