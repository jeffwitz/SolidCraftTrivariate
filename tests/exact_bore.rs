use solidcraft_trivariate::{exact_bore::build_exact_bore, exact_bore_audit::audit_exact_bore};

#[test]
fn audits_exact_circular_geometry_including_single_micrometre_clearance() {
    for (center, radius) in [
        ([0, 0], 4000),
        ([5999, 0], 4000),
        ([5999, -5999], 4000),
        ([-1, 9998], 1),
        ([0, 0], 9999),
    ] {
        let bore = build_exact_bore(center[0], center[1], radius).unwrap();
        let cert = audit_exact_bore(center, radius, &bore.blocks).unwrap();
        assert_eq!(
            (
                cert.incidence.interfaces.len(),
                cert.incidence.exterior.len()
            ),
            (8, 32)
        );
        assert_eq!(cert.radial_error_upper_mm, 0.);
        assert!(cert.determinant_lower_bound_mm3 > 0.);
    }
}

#[test]
fn endpoints_clearance_and_projective_tampering_have_separate_checks() {
    for (center, radius) in [([6000, 0], 4000), ([0, 0], 0), ([0, i32::MIN], 100)] {
        assert!(build_exact_bore(center[0], center[1], radius).is_err());
    }
    let mut bore = build_exact_bore(5999, 0, 4000).unwrap();
    for p in &mut bore.blocks[0].controls {
        for c in p {
            *c *= 3;
        }
    }
    assert!(audit_exact_bore([5999, 0], 4000, &bore.blocks).is_ok());
    bore.blocks[0].controls[6][0] += 1;
    bore.blocks[0].controls[7][0] += 1;
    assert!(audit_exact_bore([5999, 0], 4000, &bore.blocks).is_err());
}
