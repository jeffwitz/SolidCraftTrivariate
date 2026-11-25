//! Independent audit for the rational, quadrant-aligned circular-bend family.
//!
//! Once the recovered section passes `section_audit`, the map is
//! `(ξ,η,t) -> ((R+ξ) cos θ,(R+ξ) sin θ,-η)` with
//! `θ=qπ/2+2 atan(t)`. The recovered section bound and the circular radial
//! error give a strict positive lower radius; `dθ/dt=2/(1+t²) >= 1`. The
//! Jacobian is therefore the section Jacobian times `(R+ξ)dθ/dt`, preserving
//! orientation. Polar inversion is unique on each open angular range, and the
//! four-quadrant endpoint is a single exact seam checked by whole-face
//! incidence. This proves injectivity and coverage for these audited section
//! templates; it does not apply to arbitrary profiles or unrestricted angles.

use crate::bracket::BracketBlock;
use crate::bracket_audit::BracketCertificate;
use crate::circular_bend::BendRecipe;
use crate::incidence::{enumerate_faces, BlockSide, Incidence};
use crate::section_audit::{audit_section, SectionCertificate};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq)]
pub struct BendCertificate {
    pub certificate: BracketCertificate,
    pub ideal_volume_mm3: f64,
    pub cavities: usize,
}

fn gcd_i128(mut a: i128, mut b: i128) -> Result<i128, &'static str> {
    a = a.checked_abs().ok_or("circular_bend_scale_overflow")?;
    b = b.checked_abs().ok_or("circular_bend_scale_overflow")?;
    while b != 0 {
        (a, b) = (b, a % b);
    }
    Ok(a)
}

fn normalized(points: &[[i128; 4]]) -> Result<Vec<[i128; 4]>, &'static str> {
    let mut g = 0;
    for p in points {
        for x in p {
            g = gcd_i128(g, *x)?;
        }
    }
    if g <= 0 {
        return Err("circular_bend_zero_net");
    }
    Ok(points.iter().map(|p| p.map(|x| x / g)).collect())
}

fn normalized_input(block: &BracketBlock) -> Result<Vec<[i128; 4]>, &'static str> {
    if block.controls.iter().any(|p| p[3] <= 0) {
        return Err("circular_bend_nonpositive_weight");
    }
    normalized(
        &block
            .controls
            .iter()
            .map(|p| p.map(i128::from))
            .collect::<Vec<_>>(),
    )
}

fn parse_id(id: &str, quadrants: u8) -> Result<(u8, String), &'static str> {
    let rest = id.strip_prefix("bend").ok_or("circular_bend_block_id")?;
    let (q, source) = rest.split_once('_').ok_or("circular_bend_block_id")?;
    let q = q.parse::<u8>().map_err(|_| "circular_bend_block_id")?;
    if q >= quadrants || source.is_empty() {
        return Err("circular_bend_block_id");
    }
    Ok((q, source.to_owned()))
}

fn checked_i64(value: i128) -> Result<i64, &'static str> {
    i64::try_from(value).map_err(|_| "circular_bend_control_range")
}

/// Invert the angle-zero chart to recover the represented certified section.
fn recover_section_block(
    recipe: BendRecipe,
    source_id: String,
    bent: &BracketBlock,
) -> Result<BracketBlock, &'static str> {
    if bent.degrees[2] != 2
        || bent.controls.len() != bent.degrees.iter().map(|d| d + 1).product::<usize>()
    {
        return Err("circular_bend_block_shape");
    }
    let points = normalized_input(bent)?;
    let mut source = Vec::with_capacity(points.len() / 3 * 2);
    for triplet in points.as_chunks::<3>().0 {
        let [x, y, z, w] = triplet[0];
        let [x1, y1, z1, w1] = triplet[1];
        let [x2, y2, z2, w2] = triplet[2];
        if w <= 0
            || y != 0
            || x2 != 0
            || y2 != 2 * x
            || w1 != w
            || w2 != 2 * w
            || x1 != x
            || y1 != x
            || z1 != z
            || z2 != 2 * z
        {
            return Err("circular_bend_angle_zero_chart");
        }
        let section_x = x
            .checked_sub(
                i128::from(recipe.major_radius_mm)
                    .checked_mul(w)
                    .ok_or("circular_bend_overflow")?,
            )
            .ok_or("circular_bend_overflow")?;
        let section_y = z.checked_neg().ok_or("circular_bend_overflow")?;
        source.push([section_x, section_y, 0, w]);
        source.push([section_x, section_y, 8 * w, w]);
    }
    let controls = source
        .into_iter()
        .map(|p| {
            Ok([
                checked_i64(p[0])?,
                checked_i64(p[1])?,
                checked_i64(p[2])?,
                checked_i64(p[3])?,
            ])
        })
        .collect::<Result<Vec<_>, &'static str>>()?;
    Ok(BracketBlock {
        id: source_id,
        degrees: [bent.degrees[0], bent.degrees[1], 1],
        controls,
    })
}

fn expected_bent_block(
    recipe: BendRecipe,
    quadrant: u8,
    section: &BracketBlock,
) -> Result<BracketBlock, &'static str> {
    let [p, q, r] = section.degrees;
    if r != 1 || section.controls.len() != (p + 1) * (q + 1) * 2 {
        return Err("circular_bend_source_shape");
    }
    let mut controls = Vec::with_capacity(section.controls.len() / 2 * 3);
    for pair in section.controls.as_chunks::<2>().0 {
        let [x, y, z, w] = pair[0].map(i128::from);
        let [x1, y1, z1, w1] = pair[1].map(i128::from);
        if w <= 0 || x1 != x || y1 != y || w1 != w || z != 0 || z1 != 8 * w {
            return Err("circular_bend_source_extrusion");
        }
        let radius = i128::from(recipe.major_radius_mm)
            .checked_mul(w)
            .and_then(|v| v.checked_add(x))
            .ok_or("circular_bend_overflow")?;
        for (c, s, d) in [(1_i128, 0_i128, 1_i128), (1, 1, 1), (0, 2, 2)] {
            let (mut a, mut b) = (radius * c, radius * s);
            for _ in 0..quadrant {
                (a, b) = (-b, a);
            }
            controls.push([
                checked_i64(a)?,
                checked_i64(b)?,
                checked_i64(
                    y.checked_neg()
                        .and_then(|v| v.checked_mul(d))
                        .ok_or("circular_bend_overflow")?,
                )?,
                checked_i64(w.checked_mul(d).ok_or("circular_bend_overflow")?)?,
            ]);
        }
    }
    Ok(BracketBlock {
        id: format!("bend{quadrant}_{}", section.id),
        degrees: [p, q, 2],
        controls,
    })
}

fn minimum_radius_lower(
    recipe: BendRecipe,
    section: &SectionCertificate,
) -> Result<f64, &'static str> {
    // The square walls are represented exactly. Circle-based section families
    // inherit the continuous radial-error bound from section_audit.
    let radial_error = match recipe.section {
        crate::circular_bend::BendSection::Disk | crate::circular_bend::BendSection::Annulus => {
            section.radial_error_upper_mm
        }
        crate::circular_bend::BendSection::Rectangle
        | crate::circular_bend::BendSection::PiercedSquare => 0.0,
    };
    let exact_margin = f64::from(recipe.major_radius_mm - recipe.section.extent_mm());
    let lower = next_down(exact_margin - radial_error);
    if !lower.is_finite() || lower <= 0.0 {
        return Err("circular_bend_axis_intersection");
    }
    Ok(lower)
}

fn next_down(value: f64) -> f64 {
    if value.is_nan() || value == f64::NEG_INFINITY {
        return value;
    }
    if value == 0.0 {
        return -f64::from_bits(1);
    }
    let bits = value.to_bits();
    f64::from_bits(if value > 0.0 { bits - 1 } else { bits + 1 })
}

fn transformed_design_faces(
    recipe: BendRecipe,
    incidence: &Incidence,
    section: &SectionCertificate,
) -> Result<Vec<(BlockSide, String)>, &'static str> {
    let source: BTreeMap<_, _> = section
        .design_faces
        .iter()
        .map(|(s, label)| ((s.block.clone(), s.axis, s.upper), label.as_str()))
        .collect();
    let mut out = Vec::with_capacity(incidence.exterior.len());
    for side in &incidence.exterior {
        let (quadrant, source_id) = parse_id(&side.block, recipe.quadrants)?;
        let label = if side.axis == 2 {
            if recipe.quadrants == 4 {
                return Err("circular_bend_unmatched_periodic_cap");
            }
            if quadrant == 0 && !side.upper {
                "start".to_owned()
            } else if quadrant + 1 == recipe.quadrants && side.upper {
                "end".to_owned()
            } else {
                return Err("circular_bend_unexpected_cap");
            }
        } else {
            let base = source
                .get(&(source_id, side.axis, side.upper))
                .ok_or("circular_bend_unclassified_side")?;
            format!("section_{base}")
        };
        out.push((side.clone(), label));
    }
    Ok(out)
}

/// Audit bend nets by recovering their source section and proving the
/// revolved chart equations independently of `build_circular_bend`.
pub fn audit_circular_bend(
    recipe: BendRecipe,
    blocks: &[BracketBlock],
) -> Result<BendCertificate, &'static str> {
    recipe.validate()?;
    let mut keyed = BTreeMap::<(u8, String), &BracketBlock>::new();
    for b in blocks {
        let key = parse_id(&b.id, recipe.quadrants)?;
        if keyed.insert(key, b).is_some() {
            return Err("circular_bend_duplicate_block");
        }
    }
    let source_ids: BTreeSet<_> = keyed
        .keys()
        .filter(|(q, _)| *q == 0)
        .map(|(_, id)| id.clone())
        .collect();
    if source_ids.is_empty() || keyed.len() != source_ids.len() * usize::from(recipe.quadrants) {
        return Err("circular_bend_block_count");
    }
    for q in 0..recipe.quadrants {
        for id in &source_ids {
            if !keyed.contains_key(&(q, id.clone())) {
                return Err("circular_bend_missing_block");
            }
        }
    }
    let mut sections = Vec::with_capacity(source_ids.len());
    for id in &source_ids {
        sections.push(recover_section_block(
            recipe,
            id.clone(),
            keyed[&(0, id.clone())],
        )?);
    }
    let section = audit_section(recipe.section.family(), &sections)?;
    let mut expected = Vec::with_capacity(blocks.len());
    for q in 0..recipe.quadrants {
        for s in &sections {
            expected.push(expected_bent_block(recipe, q, s)?);
        }
    }
    let expected_map: BTreeMap<_, _> = expected
        .iter()
        .map(|b| {
            Ok((
                b.id.clone(),
                normalized(
                    &b.controls
                        .iter()
                        .map(|p| p.map(i128::from))
                        .collect::<Vec<_>>(),
                )?,
            ))
        })
        .collect::<Result<_, &'static str>>()?;
    let actual_map: BTreeMap<_, _> = blocks
        .iter()
        .map(|b| Ok((b.id.clone(), normalized_input(b)?)))
        .collect::<Result<_, &'static str>>()?;
    if expected_map != actual_map {
        return Err("circular_bend_control_mismatch");
    }
    let incidence = enumerate_faces(blocks)?;
    let expected_incidence = enumerate_faces(&expected)?;
    if incidence != expected_incidence {
        return Err("circular_bend_incidence_mismatch");
    }
    if recipe.quadrants < 4 {
        if !incidence.exterior.iter().any(|s| s.axis == 2 && !s.upper)
            || !incidence.exterior.iter().any(|s| s.axis == 2 && s.upper)
        {
            return Err("circular_bend_caps_missing");
        }
    } else if incidence.exterior.iter().any(|s| s.axis == 2) {
        return Err("circular_bend_periodic_seam_open");
    }
    let radius_lower = minimum_radius_lower(recipe, &section)?;
    let determinant_lower = (section.determinant_lower_bound_mm3 / 8.0).next_down() * radius_lower;
    let determinant_lower = next_down(determinant_lower);
    if !determinant_lower.is_finite() || determinant_lower <= 0.0 {
        return Err("circular_bend_determinant");
    }
    let volume = (section.volume_mm3 / 8.0)
        * f64::from(recipe.major_radius_mm)
        * f64::from(recipe.quadrants)
        * std::f64::consts::FRAC_PI_2;
    let cavities = usize::from(
        recipe.quadrants == 4
            && matches!(
                recipe.section,
                crate::circular_bend::BendSection::Annulus
                    | crate::circular_bend::BendSection::PiercedSquare
            ),
    );
    let design_faces = transformed_design_faces(recipe, &incidence, &section)?;
    Ok(BendCertificate {
        certificate: BracketCertificate {
            incidence,
            determinant_lower_bound_mm3: determinant_lower,
            radial_error_upper_mm: section.radial_error_upper_mm,
            design_faces,
        },
        ideal_volume_mm3: volume,
        cavities,
    })
}
