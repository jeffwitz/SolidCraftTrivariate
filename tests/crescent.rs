use solidcraft_trivariate::{
    crescent_audit::audit_crescent_pocket, crescent_pocket::crescent_pocket_blocks,
};

#[test]
fn all_depths_have_the_expected_incidence_and_floor_selection() {
    for depth in 0..=16 {
        let blocks = crescent_pocket_blocks(depth).unwrap();
        let certificate = audit_crescent_pocket(depth, &blocks).unwrap();
        let expected = match depth {
            0 => (14, 20, 0),
            16 => (8, 8, 0),
            _ => (22, 36, 6),
        };
        assert_eq!(blocks.len(), expected.0);
        assert_eq!(certificate.incidence.interfaces.len(), expected.1);
        assert_eq!(
            certificate
                .design_faces
                .iter()
                .filter(|(_, f)| f == "pocket_floor")
                .count(),
            expected.2
        );
        assert!(certificate.determinant_lower_bound_mm3 > 0.);
    }
    for depth in [-1, 17] {
        assert!(crescent_pocket_blocks(depth).is_err());
    }
}

#[test]
fn the_audit_rejects_altered_arc_seam_height_and_inventory() {
    let blocks = crescent_pocket_blocks(8).unwrap();
    for (index, component) in [(0, 0), (6, 1), (1, 2), (2, 3)] {
        let mut changed = blocks.clone();
        let patch = changed
            .iter_mut()
            .find(|b| b.id == "g10_crescent_upper_a")
            .unwrap();
        patch.controls[index][component] += 1;
        if component != 2 {
            patch.controls[index ^ 1][component] += 1;
        }
        assert!(audit_crescent_pocket(8, &changed).is_err());
    }
    assert!(audit_crescent_pocket(7, &blocks).is_err());
    assert!(audit_crescent_pocket(8, &blocks[1..]).is_err());
}
