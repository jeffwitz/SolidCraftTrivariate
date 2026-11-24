use solidcraft_trivariate::{
    bracket::bracket_blocks,
    bracket_audit::audit_bracket,
};
use std::collections::HashSet;

fn translate(block: &mut solidcraft_trivariate::bracket::BracketBlock, axis: usize, amount: i64) {
    for point in &mut block.controls {
        point[axis] += amount * point[3];
    }
}

#[test]
fn fixed_bracket_audit_reports_continuous_margin_and_design_faces() {
    let blocks = bracket_blocks().unwrap();
    let certificate = audit_bracket(&blocks).expect("fixed bracket should certify");
    assert_eq!(certificate.incidence.interfaces.len(), 55);
    assert_eq!(certificate.incidence.exterior.len(), 106);
    assert_eq!(certificate.design_faces.len(), 106);
    assert_eq!(certificate.determinant_lower_bound_mm3, 160.0);
    assert!(certificate.radial_error_upper_mm <= 1e-12);
    let features: HashSet<_> = certificate
        .design_faces
        .iter()
        .map(|(_, feature)| feature.as_str())
        .collect();
    assert_eq!(features.len(), 11);
}

#[test]
fn audit_is_stable_under_input_block_reordering() {
    let blocks = bracket_blocks().unwrap();
    let expected = audit_bracket(&blocks).unwrap();
    let mut reversed = blocks;
    reversed.reverse();
    let actual = audit_bracket(&reversed).expect("reordering must not change the partition");
    assert_eq!(actual.incidence, expected.incidence);
    assert_eq!(actual.design_faces, expected.design_faces);
    assert_eq!(actual.determinant_lower_bound_mm3, expected.determinant_lower_bound_mm3);
    assert_eq!(actual.radial_error_upper_mm, expected.radial_error_upper_mm);
}

#[test]
fn audit_rejects_included_cell_translated_hole_twisted_extrusion_and_wrong_elbow_cut() {
    let original = bracket_blocks().unwrap();

    let mut included = original.clone();
    let cell = included.iter_mut().find(|block| block.id == "arm_x_50_60_z0_10").unwrap();
    translate(cell, 0, -20);
    assert!(audit_bracket(&included).is_err(), "contained box must be refused");

    let mut shifted_hole = original.clone();
    for block in shifted_hole.iter_mut().filter(|block| block.id.starts_with("hole_y_")) {
        translate(block, 1, 1);
    }
    assert!(audit_bracket(&shifted_hole).is_err(), "translated bore must be refused");

    let mut twisted = original.clone();
    let elbow = twisted.iter_mut().find(|block| block.id == "elbow_east_z0_10").unwrap();
    for (index, point) in elbow.controls.iter_mut().enumerate() {
        if index % 2 == 1 {
            point[0] += point[3];
        }
    }
    assert!(audit_bracket(&twisted).is_err(), "non-product extrusion must be refused");

    let mut wrong_cut = original;
    let elbow = wrong_cut.iter_mut().find(|block| block.id == "elbow_north_z0_10").unwrap();
    for j in 0..2 {
        for k in 0..2 {
            elbow.controls.swap(j * 2 + k, 8 + j * 2 + k);
        }
    }
    assert!(audit_bracket(&wrong_cut).is_err(), "reversed elbow cut must be refused");
}

#[test]
fn audit_rejects_missing_identity_and_handles_i64_extremes_without_panicking() {
    let original = bracket_blocks().unwrap();
    let mut missing_id = original.clone();
    missing_id[0].id.clear();
    assert!(audit_bracket(&missing_id).is_err());

    for extreme in [i64::MIN, i64::MAX] {
        let mut malformed = original.clone();
        malformed[0].controls[0][0] = extreme;
        let result = std::panic::catch_unwind(|| audit_bracket(&malformed));
        assert!(result.is_ok(), "auditing extreme coefficient must not panic");
        assert!(result.unwrap().is_err(), "unrepresentable profile must be refused");
    }
}
