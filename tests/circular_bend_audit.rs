use solidcraft_trivariate::circular_bend::{build_circular_bend, BendRecipe, BendSection};
use solidcraft_trivariate::circular_bend_audit::audit_circular_bend;
use solidcraft_trivariate::sections::section_blocks;

fn sections() -> [BendSection; 4] {
    [
        BendSection::Rectangle,
        BendSection::Disk,
        BendSection::Annulus,
        BendSection::PiercedSquare,
    ]
}

fn area(section: BendSection) -> f64 {
    match section {
        BendSection::Rectangle => 400.0,
        BendSection::Disk => 16.0 * std::f64::consts::PI,
        BendSection::Annulus => 12.0 * std::f64::consts::PI,
        BendSection::PiercedSquare => 400.0 - 16.0 * std::f64::consts::PI,
    }
}

#[test]
fn frozen_section_and_quadrant_corpus_is_certified() {
    for section in sections() {
        for (quadrants, radius) in [(1, 16), (4, 16), (2, 24), (3, 24)] {
            let recipe = BendRecipe {
                section,
                major_radius_mm: radius,
                quadrants,
            };
            let blocks = build_circular_bend(recipe).unwrap();
            let report = audit_circular_bend(recipe, &blocks).unwrap();
            let expected = area(section)
                * f64::from(radius)
                * f64::from(quadrants)
                * std::f64::consts::FRAC_PI_2;
            assert!((report.ideal_volume_mm3 / expected - 1.0).abs() < 1e-12);
            assert!(report.certificate.determinant_lower_bound_mm3 > 0.0);
            let expected_cavities = usize::from(
                quadrants == 4
                    && matches!(section, BendSection::Annulus | BendSection::PiercedSquare),
            );
            assert_eq!(report.cavities, expected_cavities);
            if quadrants < 4 {
                assert!(report
                    .certificate
                    .design_faces
                    .iter()
                    .any(|(_, label)| label == "start"));
                assert!(report
                    .certificate
                    .design_faces
                    .iter()
                    .any(|(_, label)| label == "end"));
            } else {
                assert!(!report
                    .certificate
                    .incidence
                    .exterior
                    .iter()
                    .any(|s| s.axis == 2));
                let seam = report
                    .certificate
                    .incidence
                    .interfaces
                    .iter()
                    .filter(|i| {
                        let a = i.first.block.starts_with("bend0_")
                            && i.first.axis == 2
                            && !i.first.upper
                            && i.second.block.starts_with("bend3_")
                            && i.second.axis == 2
                            && i.second.upper;
                        let b = i.second.block.starts_with("bend0_")
                            && i.second.axis == 2
                            && !i.second.upper
                            && i.first.block.starts_with("bend3_")
                            && i.first.axis == 2
                            && i.first.upper;
                        a || b
                    })
                    .count();
                assert_eq!(seam, section_blocks(section.family()).unwrap().len());
            }
            assert!(report
                .certificate
                .design_faces
                .iter()
                .any(|(_, label)| label.starts_with("section_")));
        }
    }
}

#[test]
fn minimum_and_maximum_radii_and_turn_counts() {
    for (section, radius) in [
        (BendSection::Annulus, 5),
        (BendSection::PiercedSquare, 11),
        (BendSection::Annulus, 64),
        (BendSection::PiercedSquare, 64),
    ] {
        let recipe = BendRecipe {
            section,
            major_radius_mm: radius,
            quadrants: 4,
        };
        let blocks = build_circular_bend(recipe).unwrap();
        assert!(
            audit_circular_bend(recipe, &blocks)
                .unwrap()
                .certificate
                .determinant_lower_bound_mm3
                > 0.0
        );
    }
    for (section, radius, quadrants) in [
        (BendSection::Annulus, 4, 1),
        (BendSection::PiercedSquare, 10, 4),
        (BendSection::Rectangle, 65, 1),
        (BendSection::Disk, 16, 0),
        (BendSection::Annulus, 16, 5),
    ] {
        let recipe = BendRecipe {
            section,
            major_radius_mm: radius,
            quadrants,
        };
        assert!(build_circular_bend(recipe).is_err());
        assert!(audit_circular_bend(recipe, &[]).is_err());
    }
}

#[test]
fn control_corruption_and_fictitious_periodic_seam_are_rejected() {
    let recipe = BendRecipe {
        section: BendSection::Annulus,
        major_radius_mm: 24,
        quadrants: 4,
    };
    let mut blocks = build_circular_bend(recipe).unwrap();
    blocks[0].controls[3][0] += 1;
    assert!(audit_circular_bend(recipe, &blocks).is_err());

    let mut blocks = build_circular_bend(recipe).unwrap();
    blocks.retain(|b| !b.id.starts_with("bend3_"));
    assert!(audit_circular_bend(recipe, &blocks).is_err());

    let partial = BendRecipe {
        quadrants: 3,
        ..recipe
    };
    let blocks = build_circular_bend(partial).unwrap();
    assert!(audit_circular_bend(recipe, &blocks).is_err());
}

#[test]
fn each_block_accepts_a_global_projective_gauge_only() {
    let recipe = BendRecipe {
        section: BendSection::PiercedSquare,
        major_radius_mm: 16,
        quadrants: 2,
    };
    let mut blocks = build_circular_bend(recipe).unwrap();
    for block in &mut blocks {
        for point in &mut block.controls {
            for value in point {
                *value *= 2;
            }
        }
    }
    assert!(audit_circular_bend(recipe, &blocks).is_ok());
    blocks[0].controls[0][0] += 1;
    assert!(audit_circular_bend(recipe, &blocks).is_err());
}
