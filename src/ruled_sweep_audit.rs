//! Independent certification for ruled section sweeps.

use crate::bracket::BracketBlock;
use crate::bracket_audit::BracketCertificate;
use crate::incidence::enumerate_faces;
use crate::ruled_sweep::SweepRecipe;
use crate::section_audit::audit_section;
use crate::sections::section_blocks;
use std::collections::BTreeMap;

/// Certificate for the exact ruled-sweep family.
#[derive(Debug, Clone, PartialEq)]
pub struct RuledSweepCertificate {
    pub certificate: BracketCertificate,
    pub ideal_volume_mm3: f64,
    /// Exact minimum of the squared planar scale along the sweep segment.
    pub minimum_squared_scale: [i64; 2],
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn canonical_net(controls: &[[i64; 4]]) -> Result<Vec<[i128; 4]>, &'static str> {
    let mut divisor = 0_i128;
    for point in controls {
        for coefficient in point {
            divisor = gcd(divisor, i128::from(*coefficient));
        }
    }
    if divisor <= 0 {
        return Err("ruled_sweep_zero_net");
    }
    if controls.iter().any(|point| point[3] <= 0) {
        return Err("ruled_sweep_nonpositive_weight");
    }
    Ok(controls
        .iter()
        .map(|point| point.map(|value| i128::from(value) / divisor))
        .collect())
}

fn canonical_wide_net(controls: &[[i128; 4]]) -> Result<Vec<[i128; 4]>, &'static str> {
    let mut divisor = 0_i128;
    for point in controls {
        for coefficient in point {
            divisor = gcd(divisor, *coefficient);
        }
    }
    if divisor <= 0 {
        return Err("ruled_sweep_zero_expected_net");
    }
    Ok(controls
        .iter()
        .map(|point| point.map(|value| value / divisor))
        .collect())
}

fn expected_net(block: &BracketBlock, recipe: SweepRecipe) -> Result<Vec<[i128; 4]>, &'static str> {
    if !matches!(block.degrees, [1, 1, 1] | [2, 1, 1])
        || block.controls.len() != (block.degrees[0] + 1) * 4
    {
        return Err("ruled_sweep_invalid_source_block");
    }
    let [a, b] = recipe.top_multiplier_eighths.map(i128::from);
    let [tx, ty] = recipe.top_translation_eighths_mm.map(i128::from);
    let mut expected = Vec::with_capacity(block.controls.len());
    for pair in block.controls.as_chunks::<2>().0 {
        let [x, y, _, weight] = pair[0].map(i128::from);
        if weight <= 0 {
            return Err("ruled_sweep_nonpositive_source_weight");
        }
        expected.push([8 * x, 8 * y, 0, 8 * weight]);
        expected.push([
            a * x - b * y + tx * weight,
            b * x + a * y + ty * weight,
            8 * i128::from(recipe.height_mm) * weight,
            8 * weight,
        ]);
    }
    canonical_wide_net(&expected)
}

fn minimum_squared_scale(recipe: SweepRecipe) -> [i64; 2] {
    let [a, b] = recipe.top_multiplier_eighths.map(i64::from);
    let dx = a - 8;
    let denominator = dx * dx + b * b;
    if denominator == 0 {
        return [1, 1];
    }
    let stationary_numerator = 8 * (8 - a);
    if stationary_numerator <= 0 {
        [1, 1]
    } else if stationary_numerator >= denominator {
        [a * a + b * b, 64]
    } else {
        [b * b, denominator]
    }
}

/// Audit exact control relations, incidence, selections, volume, and the
/// continuous planar scale bound for a ruled sweep of a certified section.
pub fn audit_ruled_sweep(
    recipe: SweepRecipe,
    blocks: &[BracketBlock],
) -> Result<RuledSweepCertificate, &'static str> {
    recipe.validate()?;
    let family = recipe.section.family();
    let source = section_blocks(family)?;
    let section = audit_section(family, &source)?;
    if blocks.len() != source.len() {
        return Err("ruled_sweep_block_count");
    }
    let mut candidate_by_id = BTreeMap::new();
    for block in blocks {
        if candidate_by_id.insert(block.id.as_str(), block).is_some() {
            return Err("ruled_sweep_duplicate_block_id");
        }
    }
    for base in &source {
        let candidate = candidate_by_id
            .get(base.id.as_str())
            .ok_or("ruled_sweep_missing_block")?;
        if candidate.degrees != base.degrees {
            return Err("ruled_sweep_degree_mismatch");
        }
        let expected = expected_net(base, recipe)?;
        let actual = canonical_net(&candidate.controls)?;
        if actual != expected {
            return Err("ruled_sweep_control_relation_mismatch");
        }
    }

    let incidence = enumerate_faces(blocks)?;
    if incidence != section.incidence {
        return Err("ruled_sweep_incidence_mismatch");
    }
    let minimum_squared_scale = minimum_squared_scale(recipe);
    if minimum_squared_scale[0] <= 0 || minimum_squared_scale[1] <= 0 {
        return Err("ruled_sweep_collapsed_section");
    }
    // Each binary64 operation in a reported bound is rounded outward.
    let scale_min = (minimum_squared_scale[0] as f64 / minimum_squared_scale[1] as f64).next_down();
    let [a, b] = recipe.top_multiplier_eighths.map(f64::from);
    let scale_max = (1.0_f64).max((a * a + b * b) / 64.0).sqrt();
    let height_ratio = f64::from(recipe.height_mm) / 8.0;
    let determinant_lower_bound_mm3 =
        ((section.determinant_lower_bound_mm3 * height_ratio).next_down() * scale_min).next_down();
    if !determinant_lower_bound_mm3.is_finite() || determinant_lower_bound_mm3 <= 0.0 {
        return Err("ruled_sweep_determinant_bound");
    }
    let ideal_volume_mm3 =
        section.volume_mm3 * height_ratio * (1.0 + a / 8.0 + (a * a + b * b) / 64.0) / 3.0;
    if !ideal_volume_mm3.is_finite() || ideal_volume_mm3 <= 0.0 {
        return Err("ruled_sweep_volume_overflow");
    }
    let radial_error_upper_mm = if section.radial_error_upper_mm == 0. {
        0.
    } else {
        (section.radial_error_upper_mm * scale_max.next_up()).next_up()
    };
    Ok(RuledSweepCertificate {
        certificate: BracketCertificate {
            incidence,
            determinant_lower_bound_mm3,
            radial_error_upper_mm,
            design_faces: section.design_faces,
        },
        ideal_volume_mm3,
        minimum_squared_scale,
    })
}
