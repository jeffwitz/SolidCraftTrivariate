//! Exact radius parameters for the frozen fixed-bracket recipe.
//!
//! Radius ticks are sixteenths of a millimetre. Only the four inner radial
//! control rows of each selected hole family move; the other coefficients and
//! every feature identity come from the canonical E7 builder.

use crate::bracket::{bracket_blocks, BracketBlock};
use std::collections::BTreeMap;
use std::convert::TryFrom;

const MIN_RADIUS_TICKS: u8 = 40;
const MAX_RADIUS_TICKS: u8 = 80;
const GAUGE: i128 = 16;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParameterizedBracket {
    pub radii_sixteenths_mm: [u8; 2],
    pub blocks: Vec<BracketBlock>,
    /// Sparse homogeneous derivatives with respect to physical millimetres.
    /// Each axis contains the four affected hole blocks in stable ID order.
    pub derivatives: [BTreeMap<String, Vec<[i64; 4]>>; 2],
}

fn mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b)
        .ok_or("bracket parameter arithmetic overflow")
}
fn add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b)
        .ok_or("bracket parameter arithmetic overflow")
}
fn sub(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_sub(b)
        .ok_or("bracket parameter arithmetic overflow")
}
fn out(value: i128) -> Result<i64, &'static str> {
    i64::try_from(value).map_err(|_| "bracket parameter control exceeds i64")
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

fn center_times_weight(axis: usize, k: usize, weight: i128) -> Result<[i128; 3], &'static str> {
    let axial_coordinate = if k == 0 { 0_i128 } else { 8_i128 };
    Ok(if axis == 0 {
        [
            mul(40, weight)?,
            mul(8 - axial_coordinate, weight)?,
            mul(20, weight)?,
        ]
    } else {
        [
            mul(axial_coordinate, weight)?,
            mul(32, weight)?,
            mul(20, weight)?,
        ]
    })
}

fn parameterize_hole(
    block: &mut BracketBlock,
    radius_ticks: u8,
) -> Result<Vec<[i64; 4]>, &'static str> {
    if block.degrees != [2, 1, 1] || block.controls.len() != 12 {
        return Err("invalid fixed-bracket hole block");
    }
    let axis = hole_axis(&block.id).ok_or("not a fixed-bracket hole block")?;
    let mut derivative = Vec::with_capacity(block.controls.len());
    for (index, control) in block.controls.iter_mut().enumerate() {
        let radial = (index / 2) % 2;
        let k = index % 2;
        let original = control.map(i128::from);
        let weight = original[3];
        if weight <= 0 {
            return Err("nonpositive fixed-bracket hole weight");
        }
        let center = center_times_weight(axis, k, weight)?;
        let delta = [
            sub(original[0], center[0])?,
            sub(original[1], center[1])?,
            sub(original[2], center[2])?,
        ];
        let mut d = [0_i64; 4];
        if radial == 1 {
            for component in 0..3 {
                if delta[component] % 4 != 0 {
                    return Err("E6 radial control is not divisible by four");
                }
                d[component] = out(mul(4, delta[component])?)?;
            }
        }
        derivative.push(d);

        let mut updated = [0_i64; 4];
        for component in 0..3 {
            let value = if radial == 0 {
                mul(GAUGE, original[component])?
            } else {
                let moving = mul(i128::from(radius_ticks), delta[component])?;
                if moving % 4 != 0 {
                    return Err("radius tick does not preserve integer controls");
                }
                add(mul(GAUGE, center[component])?, moving / 4)?
            };
            updated[component] = out(value)?;
        }
        updated[3] = out(mul(GAUGE, weight)?)?;
        control.copy_from_slice(&updated);
    }
    Ok(derivative)
}

/// Build the canonical fixed bracket at independent x/y hole radii.
///
/// This API intentionally takes no arbitrary model: it edits only the frozen
/// E7 recipe, never projects an unknown block complex onto the template.
pub fn build_bracket_parameters(
    radius_x_sixteenths_mm: u8,
    radius_y_sixteenths_mm: u8,
) -> Result<ParameterizedBracket, &'static str> {
    for radius in [radius_x_sixteenths_mm, radius_y_sixteenths_mm] {
        if !(MIN_RADIUS_TICKS..=MAX_RADIUS_TICKS).contains(&radius) {
            return Err("hole radius must be between 2.5 and 5 mm");
        }
    }
    let mut blocks = bracket_blocks()?;
    if blocks.len() != 36 {
        return Err("fixed-bracket block inventory mismatch");
    }
    let mut derivatives: [BTreeMap<String, Vec<[i64; 4]>>; 2] =
        std::array::from_fn(|_| BTreeMap::new());
    let mut hole_counts = [0_usize; 2];
    for block in &mut blocks {
        if let Some(axis) = hole_axis(&block.id) {
            let radius = if axis == 0 {
                radius_x_sixteenths_mm
            } else {
                radius_y_sixteenths_mm
            };
            let derivative = parameterize_hole(block, radius)?;
            if derivatives[axis]
                .insert(block.id.clone(), derivative)
                .is_some()
            {
                return Err("duplicate fixed-bracket hole ID");
            }
            hole_counts[axis] += 1;
        }
    }
    if hole_counts != [4, 4] {
        return Err("fixed-bracket hole-sector inventory mismatch");
    }
    Ok(ParameterizedBracket {
        radii_sixteenths_mm: [radius_x_sixteenths_mm, radius_y_sixteenths_mm],
        blocks,
        derivatives,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scaled(mut block: BracketBlock) -> BracketBlock {
        for p in &mut block.controls {
            for x in p {
                *x *= 16;
            }
        }
        block
    }

    #[test]
    fn base_radius_is_exactly_the_canonical_bracket_at_common_gauge() {
        let reference = bracket_blocks().unwrap();
        let model = build_bracket_parameters(64, 64).unwrap();
        assert_eq!(model.blocks.len(), 36);
        for (actual, source) in model.blocks.iter().zip(reference) {
            assert_eq!(actual.id, source.id);
            assert_eq!(actual.degrees, source.degrees);
            if hole_axis(&source.id).is_some() {
                assert_eq!(actual, &scaled(source));
            } else {
                assert_eq!(actual, &source);
            }
        }
        assert_eq!(model.derivatives.map(|d| d.len()), [4, 4]);
    }

    #[test]
    fn independent_radius_edit_changes_only_its_four_projective_blocks() {
        let base = build_bracket_parameters(64, 64).unwrap();
        let edited = build_bracket_parameters(56, 64).unwrap();
        let mut changed = Vec::new();
        for (a, b) in edited.blocks.iter().zip(&base.blocks) {
            if a != b {
                changed.push(a.id.as_str());
            }
        }
        assert_eq!(changed.len(), 4);
        assert!(changed.iter().all(|id| id.starts_with("hole_x_sector_")));
        assert_eq!(
            edited.derivatives[0].keys().cloned().collect::<Vec<_>>(),
            [
                "hole_x_sector_east_z10_30",
                "hole_x_sector_north_z10_30",
                "hole_x_sector_south_z10_30",
                "hole_x_sector_west_z10_30"
            ]
            .map(str::to_owned)
        );
        assert_eq!(
            edited.derivatives[1].keys().cloned().collect::<Vec<_>>(),
            [
                "hole_y_sector_east_z10_30",
                "hole_y_sector_north_z10_30",
                "hole_y_sector_south_z10_30",
                "hole_y_sector_west_z10_30"
            ]
            .map(str::to_owned)
        );
    }

    #[test]
    fn analytic_radius_jet_matches_exact_tick_differences() {
        let q = 55_u8;
        let minus = build_bracket_parameters(q - 1, 64).unwrap();
        let plus = build_bracket_parameters(q + 1, 64).unwrap();
        let center = build_bracket_parameters(q, 64).unwrap();
        for (id, derivative) in &center.derivatives[0] {
            let m = minus.blocks.iter().find(|b| &b.id == id).unwrap();
            let p = plus.blocks.iter().find(|b| &b.id == id).unwrap();
            for (i, row) in derivative.iter().enumerate() {
                for (c, value) in row.iter().enumerate() {
                    assert_eq!(
                        i128::from(p.controls[i][c] - m.controls[i][c]) * 8,
                        i128::from(*value)
                    );
                }
            }
        }
        for derivative in center.derivatives[0].values() {
            assert!(derivative.iter().all(|p| p[3] == 0));
        }
    }

    #[test]
    fn radius_bounds_and_i64_capacity_are_enforced() {
        assert!(build_bracket_parameters(39, 64).is_err());
        assert!(build_bracket_parameters(64, 81).is_err());
        for q in [40, 80] {
            let model = build_bracket_parameters(q, q).unwrap();
            assert!(model
                .blocks
                .iter()
                .flat_map(|b| &b.controls)
                .all(|p| p[3] > 0));
        }
    }
}
