use solidcraft_trivariate::coaxial_audit::audit_coaxial_program;
use solidcraft_trivariate::coaxial_program::{
    build_coaxial_program, CoaxialProgram, CylinderFeature, CylinderOp,
};

fn program(features: &[(CylinderOp, i32, [i32; 2], &str)]) -> CoaxialProgram {
    CoaxialProgram {
        base_height_mm: 8,
        features: features
            .iter()
            .map(|(op, r, z, id)| CylinderFeature {
                id: (*id).into(),
                op: *op,
                radius_mm: *r,
                z_mm: *z,
            })
            .collect(),
    }
}

fn admitted(p: CoaxialProgram, volume: f64, cavities: usize) {
    let blocks = build_coaxial_program(&p).unwrap();
    let audit = audit_coaxial_program(&p, &blocks).unwrap();
    assert!(
        (audit.ideal_volume_mm3 - volume).abs() <= 1e-10 * volume.max(1.0),
        "volume {} differs from oracle {} for {:?}",
        audit.ideal_volume_mm3,
        volume,
        p
    );
    assert_eq!(audit.cavities, cavities);
    assert!(audit.certificate.determinant_lower_bound_mm3 > 0.0);
    if blocks.len() > 1 {
        assert!(!audit.certificate.incidence.interfaces.is_empty());
    }
    assert!(!audit.certificate.design_faces.is_empty());
}

#[test]
fn frozen_positive_programs_match_independent_volumes_and_cavities() {
    let pi = std::f64::consts::PI;
    admitted(program(&[]), 3200.0, 0);
    admitted(
        program(&[(CylinderOp::Cut, 4, [0, 8], "hole")]),
        3200.0 - 128.0 * pi,
        0,
    );
    admitted(
        program(&[(CylinderOp::Cut, 4, [4, 8], "pocket")]),
        3200.0 - 64.0 * pi,
        0,
    );
    admitted(
        program(&[
            (CylinderOp::Cut, 2, [0, 8], "hole"),
            (CylinderOp::Cut, 5, [5, 8], "seat"),
        ]),
        3200.0 - 95.0 * pi,
        0,
    );
    admitted(
        program(&[
            (CylinderOp::Add, 6, [8, 12], "boss"),
            (CylinderOp::Cut, 2, [0, 12], "hole"),
        ]),
        3200.0 + 96.0 * pi,
        0,
    );
    admitted(
        program(&[(CylinderOp::Cut, 3, [2, 6], "closed")]),
        3200.0 - 36.0 * pi,
        1,
    );
    admitted(
        program(&[
            (CylinderOp::Add, 6, [8, 12], "boss"),
            (CylinderOp::Cut, 4, [0, 12], "bore"),
        ]),
        3200.0 - 48.0 * pi,
        0,
    );
    admitted(
        program(&[
            (CylinderOp::Cut, 4, [0, 12], "bore"),
            (CylinderOp::Add, 6, [8, 12], "boss"),
        ]),
        3200.0 + 16.0 * pi,
        0,
    );
}

#[test]
fn rejects_corrupt_networks_and_invalid_topology() {
    let p = program(&[(CylinderOp::Cut, 4, [4, 8], "pocket")]);
    let mut blocks = build_coaxial_program(&p).unwrap();
    blocks[0].controls[0][0] += 1;
    assert!(audit_coaxial_program(&p, &blocks).is_err());
    let mut blocks = build_coaxial_program(&p).unwrap();
    blocks.pop();
    assert!(audit_coaxial_program(&p, &blocks).is_err());

    let floating = program(&[
        (CylinderOp::Cut, 5, [4, 8], "wide"),
        (CylinderOp::Add, 3, [6, 10], "plug"),
    ]);
    assert!(build_coaxial_program(&floating).is_err());
    let edge_contact = program(&[
        (CylinderOp::Cut, 5, [4, 8], "pocket"),
        (CylinderOp::Add, 5, [8, 10], "boss"),
    ]);
    assert!(build_coaxial_program(&edge_contact).is_err());

    let disconnected = program(&[(CylinderOp::Add, 2, [10, 12], "float")]);
    assert!(audit_coaxial_program(&disconnected, &[]).is_err());
    let duplicate = program(&[
        (CylinderOp::Cut, 2, [0, 8], "same"),
        (CylinderOp::Add, 4, [8, 10], "same"),
    ]);
    assert!(audit_coaxial_program(&duplicate, &[]).is_err());
    let zero = program(&[(CylinderOp::Cut, 0, [0, 8], "zero")]);
    assert!(audit_coaxial_program(&zero, &[]).is_err());
}

#[test]
fn floor_and_cylindrical_wall_keep_the_operation_that_created_them() {
    let p = program(&[(CylinderOp::Cut, 4, [4, 8], "pocket")]);
    let blocks = build_coaxial_program(&p).unwrap();
    let audit = audit_coaxial_program(&p, &blocks).unwrap();
    assert!(audit
        .certificate
        .design_faces
        .iter()
        .any(|(side, label)| { side.axis == 2 && side.upper && label == "pocket" }));
    assert!(audit
        .certificate
        .design_faces
        .iter()
        .any(|(side, label)| { side.axis == 1 && label == "pocket" }));
    assert!(audit
        .certificate
        .design_faces
        .iter()
        .any(|(_, label)| label == "base"));
}

#[test]
fn accepts_global_projective_gauge_and_rejects_local_perturbation() {
    let p = program(&[(CylinderOp::Cut, 4, [0, 8], "hole")]);
    let mut blocks = build_coaxial_program(&p).unwrap();
    for block in &mut blocks {
        for point in &mut block.controls {
            for x in point {
                *x *= 2;
            }
        }
    }
    assert!(audit_coaxial_program(&p, &blocks).is_ok());
    blocks[0].controls[0][0] += 1;
    assert!(audit_coaxial_program(&p, &blocks).is_err());
}
