//! Continuous audit of a guarded prism with linearly tapering circular bores.
//!
//! Two independently audited straight-through G21 complexes at the small and
//! large radii identify the template. Their canonical homogeneous weights,
//! outer guard geometry and angular parameter agree. Interpolating complete
//! nets in this common gauge therefore changes the circular radius linearly.
//! Every horizontal slice has the G20 guard-to-circle partition: the circle
//! stays strictly inside its diamond, the cardinal seams bound quadrants and
//! the two halves are exact restrictions of a single quadrant. The exterior
//! guard complex stays fixed modulo the prescribed axial shear. Positive
//! spatial Jacobians and increasing height give positive planar Jacobians;
//! the simple polygon/circle boundary cycles and paired internal seams give
//! degree one in the material and zero in the holes. Distinct heights cannot
//! overlap. These statements apply only after all coefficient identities pass.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    circular_polygon::CircularPolygonBore,
    constraint_system::Ratio,
    counterbores::{Counterbore, CounterboreOpening, CounterborePrism},
    incidence::{enumerate_faces, BlockSide},
    interval_bores::{build_interval_bore_prism, IntervalBore, IntervalBorePrism},
    interval_bores_audit::audit_interval_bore_prism,
    spatial_certificate::{certify_spatial, CertificateBudget, CertificateStatus},
};
use std::collections::{BTreeMap, BTreeSet};

pub struct CountersinkCertificate {
    pub certificate: BracketCertificate,
    pub outer_area_twice_mm2: i64,
    /// The volume removed is pi times this integer divided by three.
    pub removed_pi_numerator: i64,
}

fn comparator(recipe: &CounterborePrism, large: bool) -> IntervalBorePrism {
    IntervalBorePrism {
        outer: recipe.outer.clone(),
        height_mm: recipe.height_mm,
        top_shift_mm: recipe.top_shift_mm,
        bores: recipe
            .bores
            .iter()
            .map(|b| IntervalBore {
                bore: CircularPolygonBore {
                    id: b.id.clone(),
                    center_mm: b.center_mm,
                    radius_mm: if large {
                        b.counterbore_radius_mm
                    } else {
                        b.radius_mm
                    },
                    guard_mm: b.guard_mm,
                },
                z_interval_mm: [0, recipe.height_mm],
            })
            .collect(),
    }
}

fn distance_from_entry(bore: &Counterbore, h: i32, z: i32) -> i32 {
    match bore.opening {
        CounterboreOpening::Start => z,
        CounterboreOpening::End => h - z,
    }
}

fn blend(a: Ratio, b: Ratio, t: Ratio) -> Result<Ratio, &'static str> {
    a.checked_add(b.checked_sub(a)?.checked_mul(t)?)
}

fn coefficient_at_z(
    b: &BracketBlock,
    index: usize,
    axis: usize,
    z: Ratio,
) -> Result<Ratio, &'static str> {
    let gauge = i128::from(b.controls[0][3]);
    let at = |i: usize| Ratio::new(i128::from(b.controls[i][axis]), gauge);
    blend(at(index & !1)?, at(index | 1)?, z)
}

fn verify_layer(
    low: &BracketBlock,
    high: &BracketBlock,
    actual: &BracketBlock,
    bore: Option<&Counterbore>,
    layer: [i32; 2],
    h: i32,
) -> Result<f64, &'static str> {
    if low.degrees[2] != 1
        || low.degrees != high.degrees
        || actual.degrees != low.degrees
        || low.controls.len() != high.controls.len()
        || actual.controls.len() != low.controls.len()
        || actual.controls.iter().any(|p| p[3] <= 0)
    {
        return Err("countersink_layer_net_mismatch");
    }
    let gauge = i128::from(actual.controls[0][3]);
    for (i, p) in actual.controls.iter().enumerate() {
        let z = layer[i % 2];
        let t = Ratio::new(i128::from(z), i128::from(h))?;
        let opening = if let Some(b) = bore {
            Ratio::new(
                i128::from((b.depth_mm - distance_from_entry(b, h, z)).max(0)),
                i128::from(b.depth_mm),
            )?
        } else {
            Ratio::zero()
        };
        for (axis, &component) in p.iter().enumerate() {
            let lo = coefficient_at_z(low, i, axis, t)?;
            let hi = coefficient_at_z(high, i, axis, t)?;
            if (axis == 3 || bore.is_none()) && lo != hi {
                return Err("incompatible_countersink_reference_sections");
            }
            if blend(lo, hi, opening)? != Ratio::new(i128::from(component), gauge)? {
                return Err("countersink_stored_net_fails_reference_interpolation");
            }
        }
    }
    let c = certify_spatial(
        actual.degrees,
        &actual.controls,
        CertificateBudget::default(),
    );
    match c.status {
        CertificateStatus::Valid => c
            .determinant_lower_bound_mm3
            .ok_or("countersink_missing_positive_bound"),
        CertificateStatus::Invalid => Err("countersink_invalid_spatial_jacobian"),
        CertificateStatus::Indeterminate => Err("countersink_spatial_jacobian_indeterminate"),
    }
}

/// Verify the layer inventory, all stored nets, exterior selection and volume.
pub fn audit_countersinks(
    recipe: &CounterborePrism,
    blocks: &[BracketBlock],
) -> Result<CountersinkCertificate, &'static str> {
    let h = recipe.height_mm;
    if h <= 1 || !(1..=4).contains(&recipe.bores.len()) {
        return Err("countersink_recipe_range");
    }
    let mut cuts = BTreeSet::from([0, h]);
    let mut removed = 0_i128;
    for b in &recipe.bores {
        if !(1..h).contains(&b.depth_mm)
            || b.radius_mm <= 0
            || b.radius_mm >= b.counterbore_radius_mm
            || i64::from(b.guard_mm) < 2 * i64::from(b.counterbore_radius_mm)
        {
            return Err("countersink_invalid_profile");
        }
        cuts.insert(match b.opening {
            CounterboreOpening::Start => b.depth_mm,
            CounterboreOpening::End => h - b.depth_mm,
        });
        let (r, big, d) = (
            i128::from(b.radius_mm),
            i128::from(b.counterbore_radius_mm),
            i128::from(b.depth_mm),
        );
        removed += 3 * r * r * (i128::from(h) - d) + d * (big * big + big * r + r * r);
    }
    let cuts: Vec<_> = cuts.into_iter().collect();
    let low_recipe = comparator(recipe, false);
    let low = build_interval_bore_prism(&low_recipe)?;
    let reference = audit_interval_bore_prism(&low_recipe, &low)?;
    let high_recipe = comparator(recipe, true);
    let high = build_interval_bore_prism(&high_recipe)?;
    audit_interval_bore_prism(&high_recipe, &high)?;
    let index = |bs: &[BracketBlock]| {
        bs.iter()
            .enumerate()
            .map(|(i, b)| (b.id.clone(), i))
            .collect::<BTreeMap<_, _>>()
    };
    let actual = index(blocks);
    let upper = index(&high);
    if actual.len() != blocks.len()
        || actual.len() != low.len() * (cuts.len() - 1)
        || upper.len() != low.len()
    {
        return Err("countersink_block_inventory");
    }
    let suffix = format!("_z0_{h}");
    let base_tags: BTreeMap<_, _> = reference.certificate.design_faces.iter().cloned().collect();
    let mut expected_faces = BTreeMap::new();
    let mut lower_bound = f64::INFINITY;
    for pair in cuts.windows(2) {
        let layer = [pair[0], pair[1]];
        for small in &low {
            let stem = small
                .id
                .strip_suffix(&suffix)
                .ok_or("countersink_reference_suffix")?;
            let name = format!("{stem}_z{}_{}", layer[0], layer[1]);
            let big = &high[*upper
                .get(&small.id)
                .ok_or("countersink_missing_reference")?];
            let candidate = &blocks[*actual.get(&name).ok_or("countersink_missing_layer_block")?];
            let side = |axis, upper| BlockSide {
                block: small.id.clone(),
                axis,
                upper,
            };
            let bore = base_tags.get(&side(1, true)).and_then(|tag| {
                recipe
                    .bores
                    .iter()
                    .find(|b| tag == &format!("bore/{}/wall", b.id))
            });
            lower_bound = lower_bound.min(verify_layer(small, big, candidate, bore, layer, h)?);
            for axis in 0..3 {
                for top in [false, true] {
                    if axis == 2 && layer[usize::from(top)] != if top { h } else { 0 } {
                        continue;
                    }
                    if let Some(tag) = base_tags.get(&side(axis, top)) {
                        let label = if let Some(b) = bore.filter(|_| axis == 1 && top) {
                            let cone = layer
                                .iter()
                                .any(|&z| distance_from_entry(b, h, z) < b.depth_mm);
                            format!(
                                "bore/{}/{}",
                                b.id,
                                if cone { "cone_wall" } else { "pilot_wall" }
                            )
                        } else {
                            tag.clone()
                        };
                        expected_faces.insert(
                            BlockSide {
                                block: name.clone(),
                                axis,
                                upper: top,
                            },
                            label,
                        );
                    }
                }
            }
        }
    }
    let incidence = enumerate_faces(blocks)?;
    if incidence.exterior.len() != expected_faces.len()
        || incidence
            .exterior
            .iter()
            .any(|s| !expected_faces.contains_key(s))
    {
        return Err("countersink_boundary_inventory_mismatch");
    }
    Ok(CountersinkCertificate {
        outer_area_twice_mm2: reference.outer_area_twice_mm2,
        removed_pi_numerator: i64::try_from(removed).map_err(|_| "countersink_volume_range")?,
        certificate: BracketCertificate {
            incidence,
            design_faces: expected_faces.into_iter().collect(),
            determinant_lower_bound_mm3: lower_bound,
            radial_error_upper_mm: 0.,
        },
    })
}
