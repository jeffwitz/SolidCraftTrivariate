//! Analytic event and boundary-measure semantics for a translating circular bore.
//!
//! The fixed body is the 20×20×8 mm plate with a blind rectangular pocket
//! [-6,-2]×[-4,4] mm over z=[4,8] mm. A radius-2 mm vertical bore is centered
//! at (c,0). The integer input is c in micrometres. All event classification
//! uses integer comparisons; transcendental measures are reported as f64
//! evaluations of the exact circular-segment formulae below.
//!
//! This module describes the ideal physical set and its analytic measures. It
//! does not claim that a particular block decomposition represents every
//! event. At an exact tangency the physical volume remains differentiable, but
//! the representation signature changes and a representation jet is absent.

const RADIUS_UM: i64 = 2_000;
const MIN_CENTER_UM: i64 = -5_000;
const MAX_CENTER_UM: i64 = 14_000;
const CIRCLE_RADIUS_MM: f64 = 2.0;
const PLATE_X_MM: [f64; 2] = [-10.0, 10.0];
const POCKET_X_MM: [f64; 2] = [-6.0, -2.0];

/// Relation between the bore disk and the rectangular pocket footprint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BorePocketRelation {
    Disjoint,
    Tangent,
    Overlapping,
    /// The disk lies inside the pocket and touches both x walls.
    ContainedWithBoundaryContact,
}

/// Relation between the bore disk and the plate's right-hand side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BorePlateRelation {
    Interior,
    Tangent,
    OpenToExterior,
    Disjoint,
}

/// Exact status of a circular arc relative to a vertical feature line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LineRelation {
    Misses,
    Tangent,
    Crosses,
}

/// Physical-set classification, independent of how the set is partitioned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PhysicalSignature {
    pub bore_pocket: BorePocketRelation,
    pub bore_plate: BorePlateRelation,
}

/// Representation-relevant contacts against the fixed vertical arrangement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RepresentationSignature {
    pub plate_x_minus_10: LineRelation,
    pub pocket_x_minus_6: LineRelation,
    pub pocket_x_minus_2: LineRelation,
    pub plate_x_plus_10: LineRelation,
    /// At a line tangency the combinatorial cut changes and no representation
    /// derivative is supplied, even when the physical volume gradient exists.
    /// Alignment with the quadrant cut also changes the block inventory.
    /// This flag describes combinatorics, independently of integer realizability.
    pub locally_fixed_partition: bool,
    /// Position of the disk center relative to each vertical arrangement line.
    pub center_sides: [i8; 4],
    /// A transverse cut admits an exact integer-micrometre ordinate in G11b.
    pub integer_cut_ordinates: bool,
}

/// Status of the derivative of the ideal physical volume with respect to c.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhysicalGradientSignature {
    Exists,
}

/// Exact event at one of the vertical arrangement lines.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TangencyEvent {
    /// Feature name: `plate_x_minus_10`, `pocket_x_minus_6`, etc.
    pub feature: &'static str,
    /// Feature x-coordinate in micrometres.
    pub x_um: i64,
}

/// Surface measures grouped by their geometric provenance, in mm².
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundaryAreas {
    /// Lower cylindrical bore wall, z=[0,4].
    pub bore_wall_below_pocket_mm2: f64,
    /// Upper exposed cylindrical bore wall, z=[4,8], excluding the pocket.
    pub bore_wall_above_pocket_mm2: f64,
    /// Four vertical walls of the pocket, including bore truncation at x walls.
    pub pocket_vertical_walls_mm2: f64,
    /// Blind-pocket floor remaining outside the through bore.
    pub pocket_floor_mm2: f64,
    /// Plate top and bottom exposed faces after the cuts.
    pub plate_top_mm2: f64,
    pub plate_bottom_mm2: f64,
    /// Four outer vertical plate faces after the bore opens through x=+10.
    pub plate_outer_vertical_walls_mm2: f64,
}

/// Analytic metrics and three deliberately separate signatures.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TravellingBoreReport {
    pub center_um: i64,
    pub physical: PhysicalSignature,
    pub representation: RepresentationSignature,
    pub physical_gradient: PhysicalGradientSignature,
    pub events: [Option<TangencyEvent>; 4],
    /// Plate minus pocket minus in-plate bore plus bore/pocket overlap.
    pub volume_mm3: f64,
    /// Exact analytic derivative formula evaluated in mm² per mm.
    pub dvolume_dcenter_mm2: f64,
    pub areas: BoundaryAreas,
}

/// Compute the physical and representation signatures and analytic measures.
pub fn analyze_travelling_bore(center_um: i64) -> Result<TravellingBoreReport, &'static str> {
    if !(MIN_CENTER_UM..=MAX_CENTER_UM).contains(&center_um) {
        return Err("center_out_of_range");
    }
    let center_mm = center_um as f64 / 1_000.0;
    let event_lines = [
        ("plate_x_minus_10", -10_000_i64),
        ("pocket_x_minus_6", -6_000_i64),
        ("pocket_x_minus_2", -2_000_i64),
        ("plate_x_plus_10", 10_000_i64),
    ];
    let relations = event_lines.map(|(_, x)| line_relation(center_um, x));
    let events = std::array::from_fn(|index| {
        (relations[index] == LineRelation::Tangent).then_some(TangencyEvent {
            feature: event_lines[index].0,
            x_um: event_lines[index].1,
        })
    });
    let has_event = events.iter().any(Option::is_some);

    let center_sides = event_lines.map(|(_, x)| (center_um - x).signum() as i8);
    let integer_cut_ordinates = event_lines.iter().all(|&(_, x)| {
        let radicand = RADIUS_UM * RADIUS_UM - (center_um - x).pow(2);
        radicand <= 0 || radicand.isqrt().pow(2) == radicand
    });
    let pocket_relation = if center_um == -4000 {
        BorePocketRelation::ContainedWithBoundaryContact
    } else if center_um + RADIUS_UM < -6_000 || center_um - RADIUS_UM > -2_000 {
        BorePocketRelation::Disjoint
    } else if center_um + RADIUS_UM == -6_000 || center_um - RADIUS_UM == -2_000 {
        BorePocketRelation::Tangent
    } else {
        BorePocketRelation::Overlapping
    };
    let plate_relation = if center_um - RADIUS_UM >= 10_000 {
        if center_um - RADIUS_UM == 10_000 {
            BorePlateRelation::Tangent
        } else {
            BorePlateRelation::Disjoint
        }
    } else if center_um + RADIUS_UM == 10_000 {
        BorePlateRelation::Tangent
    } else if center_um + RADIUS_UM > 10_000 {
        BorePlateRelation::OpenToExterior
    } else {
        BorePlateRelation::Interior
    };

    let disk_in_plate = disk_strip_area(center_mm, PLATE_X_MM[0], PLATE_X_MM[1]);
    let overlap_area = disk_strip_area(center_mm, POCKET_X_MM[0], POCKET_X_MM[1]);
    let disk_derivative = disk_strip_area_derivative(center_mm, PLATE_X_MM[0], PLATE_X_MM[1]);
    let overlap_derivative = disk_strip_area_derivative(center_mm, POCKET_X_MM[0], POCKET_X_MM[1]);
    let volume_mm3 = 3_072.0 - 8.0 * disk_in_plate + 4.0 * overlap_area;
    let dvolume_dcenter_mm2 = -8.0 * disk_derivative + 4.0 * overlap_derivative;
    let bore_arc_plate = circle_arc_length_in_strip(center_mm, PLATE_X_MM[0], PLATE_X_MM[1]);
    let bore_arc_pocket = circle_arc_length_in_strip(center_mm, POCKET_X_MM[0], POCKET_X_MM[1]);
    let right_opening_chord = circle_chord_at_x(center_mm, PLATE_X_MM[1]);
    let pocket_left_cut = circle_chord_at_x(center_mm, POCKET_X_MM[0]);
    let pocket_right_cut = circle_chord_at_x(center_mm, POCKET_X_MM[1]);

    let areas = BoundaryAreas {
        bore_wall_below_pocket_mm2: 4.0 * bore_arc_plate,
        bore_wall_above_pocket_mm2: 4.0 * (bore_arc_plate - bore_arc_pocket).max(0.0),
        pocket_vertical_walls_mm2: 32.0
            + 4.0 * ((8.0 - pocket_left_cut).max(0.0) + (8.0 - pocket_right_cut).max(0.0)),
        pocket_floor_mm2: (32.0 - overlap_area).max(0.0),
        plate_top_mm2: (368.0 - disk_in_plate + overlap_area).max(0.0),
        plate_bottom_mm2: (400.0 - disk_in_plate).max(0.0),
        plate_outer_vertical_walls_mm2: 640.0 - 8.0 * right_opening_chord,
    };

    if !volume_mm3.is_finite()
        || !dvolume_dcenter_mm2.is_finite()
        || [
            areas.bore_wall_below_pocket_mm2,
            areas.bore_wall_above_pocket_mm2,
            areas.pocket_vertical_walls_mm2,
            areas.pocket_floor_mm2,
            areas.plate_top_mm2,
            areas.plate_bottom_mm2,
            areas.plate_outer_vertical_walls_mm2,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value < -1e-10)
    {
        return Err("analytic_measure_range");
    }

    Ok(TravellingBoreReport {
        center_um,
        physical: PhysicalSignature {
            bore_pocket: pocket_relation,
            bore_plate: plate_relation,
        },
        representation: RepresentationSignature {
            plate_x_minus_10: relations[0],
            pocket_x_minus_6: relations[1],
            pocket_x_minus_2: relations[2],
            plate_x_plus_10: relations[3],
            locally_fixed_partition: !has_event && !center_sides.contains(&0),
            center_sides,
            integer_cut_ordinates,
        },
        physical_gradient: PhysicalGradientSignature::Exists,
        events,
        volume_mm3,
        dvolume_dcenter_mm2,
        areas,
    })
}

fn line_relation(center_um: i64, line_um: i64) -> LineRelation {
    let distance = (center_um - line_um).abs();
    match distance.cmp(&RADIUS_UM) {
        std::cmp::Ordering::Less => LineRelation::Crosses,
        std::cmp::Ordering::Equal => LineRelation::Tangent,
        std::cmp::Ordering::Greater => LineRelation::Misses,
    }
}

/// Primitive of the disk's vertical chord width 2 sqrt(r²-x²), in mm².
fn disk_primitive(relative_x_mm: f64) -> f64 {
    let x = relative_x_mm.clamp(-CIRCLE_RADIUS_MM, CIRCLE_RADIUS_MM);
    let radicand = (CIRCLE_RADIUS_MM * CIRCLE_RADIUS_MM - x * x).max(0.0);
    x * radicand.sqrt() + CIRCLE_RADIUS_MM * CIRCLE_RADIUS_MM * (x / CIRCLE_RADIUS_MM).asin()
}

fn disk_strip_area(center_mm: f64, left_mm: f64, right_mm: f64) -> f64 {
    let left = (left_mm - center_mm).clamp(-CIRCLE_RADIUS_MM, CIRCLE_RADIUS_MM);
    let right = (right_mm - center_mm).clamp(-CIRCLE_RADIUS_MM, CIRCLE_RADIUS_MM);
    if right <= left {
        0.0
    } else {
        disk_primitive(right) - disk_primitive(left)
    }
}

fn chord_width(relative_x_mm: f64) -> f64 {
    if relative_x_mm.abs() >= CIRCLE_RADIUS_MM {
        0.0
    } else {
        2.0 * (CIRCLE_RADIUS_MM * CIRCLE_RADIUS_MM - relative_x_mm * relative_x_mm).sqrt()
    }
}

fn disk_strip_area_derivative(center_mm: f64, left_mm: f64, right_mm: f64) -> f64 {
    chord_width(left_mm - center_mm) - chord_width(right_mm - center_mm)
}

fn circle_arc_length_in_strip(center_mm: f64, left_mm: f64, right_mm: f64) -> f64 {
    let lower_cos = ((left_mm - center_mm) / CIRCLE_RADIUS_MM).clamp(-1.0, 1.0);
    let upper_cos = ((right_mm - center_mm) / CIRCLE_RADIUS_MM).clamp(-1.0, 1.0);
    if upper_cos <= lower_cos {
        0.0
    } else {
        CIRCLE_RADIUS_MM * 2.0 * (lower_cos.acos() - upper_cos.acos())
    }
}

fn circle_chord_at_x(center_mm: f64, x_mm: f64) -> f64 {
    chord_width(x_mm - center_mm)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn volume_formula(center_mm: f64) -> f64 {
        let plate = disk_strip_area(center_mm, -10.0, 10.0);
        let shared = disk_strip_area(center_mm, -6.0, -2.0);
        3_072.0 - 8.0 * plate + 4.0 * shared
    }

    #[test]
    fn exact_event_classification_uses_integer_supports() {
        let cases = [
            (-4_000, "pocket_x_minus_6"),
            (-4_000, "pocket_x_minus_2"),
            (0, "pocket_x_minus_2"),
            (8_000, "plate_x_plus_10"),
            (12_000, "plate_x_plus_10"),
        ];
        for (center, feature) in cases {
            let report = analyze_travelling_bore(center).unwrap();
            assert!(report
                .events
                .iter()
                .flatten()
                .any(|event| event.feature == feature));
            assert!(!report.representation.locally_fixed_partition);
            assert_eq!(report.physical_gradient, PhysicalGradientSignature::Exists);
        }
        assert_eq!(
            analyze_travelling_bore(-4_000)
                .unwrap()
                .physical
                .bore_pocket,
            BorePocketRelation::ContainedWithBoundaryContact
        );
        assert_eq!(
            analyze_travelling_bore(0).unwrap().physical.bore_pocket,
            BorePocketRelation::Tangent
        );
        assert_eq!(
            analyze_travelling_bore(8_000).unwrap().physical.bore_plate,
            BorePlateRelation::Tangent
        );
        assert_eq!(
            analyze_travelling_bore(12_000).unwrap().physical.bore_plate,
            BorePlateRelation::Tangent
        );
    }

    #[test]
    fn analytic_volume_and_derivative_match_formula_and_centered_difference() {
        for center_um in [
            -5_000, -4_001, -3_000, 1_000, 7_999, 8_001, 10_000, 11_999, 13_000,
        ] {
            let report = analyze_travelling_bore(center_um).unwrap();
            let center_mm = center_um as f64 / 1_000.0;
            assert!((report.volume_mm3 - volume_formula(center_mm)).abs() < 1e-12);
            let h = 1e-6;
            let fd = (volume_formula(center_mm + h) - volume_formula(center_mm - h)) / (2.0 * h);
            assert!(
                (report.dvolume_dcenter_mm2 - fd).abs() < 1e-5,
                "center={center_um}: analytic={}, fd={fd}",
                report.dvolume_dcenter_mm2
            );
        }
    }

    #[test]
    fn physical_gradient_survives_tangency_while_representation_jet_does_not() {
        for event in [-4_000, 0, 8_000, 12_000] {
            let exact = analyze_travelling_bore(event).unwrap();
            assert!(!exact.representation.locally_fixed_partition);
            assert!(exact.dvolume_dcenter_mm2.abs() < 1e-12);
            let left = analyze_travelling_bore(event - 1).unwrap();
            let right = analyze_travelling_bore(event + 1).unwrap();
            assert!(left.representation.locally_fixed_partition);
            assert!(right.representation.locally_fixed_partition);
            assert!(left.dvolume_dcenter_mm2.abs() < 1.1);
            assert!(right.dvolume_dcenter_mm2.abs() < 1.1);
            let event_mm = event as f64 / 1_000.0;
            let delta_mm = 0.001;
            let center_volume = volume_formula(event_mm);
            let left_secant = (center_volume - volume_formula(event_mm - delta_mm)) / delta_mm;
            let right_secant = (volume_formula(event_mm + delta_mm) - center_volume) / delta_mm;
            assert!(left_secant.abs() < 0.8, "event={event}, left={left_secant}");
            assert!(
                right_secant.abs() < 0.8,
                "event={event}, right={right_secant}"
            );
        }
    }

    #[test]
    fn boundary_areas_are_finite_nonnegative_and_report_opening_provenance() {
        let interior = analyze_travelling_bore(0).unwrap();
        assert!(
            (interior.areas.bore_wall_below_pocket_mm2 - 16.0 * std::f64::consts::PI).abs() < 1e-10
        );
        assert!(
            (interior.areas.bore_wall_above_pocket_mm2 - 16.0 * std::f64::consts::PI).abs() < 1e-10
        );
        assert!((interior.areas.pocket_vertical_walls_mm2 - 96.0).abs() < 1e-10);
        let contained = analyze_travelling_bore(-4_000).unwrap();
        assert_eq!(contained.areas.bore_wall_above_pocket_mm2, 0.0);
        let open = analyze_travelling_bore(10_000).unwrap();
        assert!(open.areas.plate_outer_vertical_walls_mm2 < 640.0);
        assert!(open.areas.bore_wall_below_pocket_mm2 < 16.0 * std::f64::consts::PI);
    }

    #[test]
    fn rejects_centers_outside_recipe_domain() {
        assert_eq!(analyze_travelling_bore(-5_001), Err("center_out_of_range"));
        assert_eq!(analyze_travelling_bore(14_001), Err("center_out_of_range"));
    }
}
