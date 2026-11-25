use solidcraft_trivariate::coaxial_audit::audit_coaxial_program;
use solidcraft_trivariate::coaxial_program::{
    build_coaxial_program, CoaxialProgram, CylinderFeature, CylinderOp,
};
use solidcraft_trivariate::coaxial_session::{CoaxialSession, FaceSelection, SelectionStatus};

fn p(features: &[(CylinderOp, i32, [i32; 2], &str)]) -> CoaxialProgram {
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

fn selection(point_mm: [f64; 3], normal: [f64; 3], feature: &str) -> FaceSelection {
    FaceSelection {
        point_mm,
        normal,
        feature: feature.into(),
    }
}

fn cold_check(session: &CoaxialSession) {
    let state = session.state().unwrap();
    let cold = build_coaxial_program(&state.program).unwrap();
    assert_eq!(state.blocks, cold);
    let report = audit_coaxial_program(&state.program, &cold).unwrap();
    assert_eq!(state.certificate, report.certificate);
    assert_eq!(state.ideal_volume_mm3, report.ideal_volume_mm3);
    assert_eq!(state.cavities, report.cavities);
}

fn oracle_material(program: &CoaxialProgram, p: [f64; 3]) -> bool {
    let mut solid = p[0].abs() < 10.0
        && p[1].abs() < 10.0
        && p[2] > 0.0
        && p[2] < f64::from(program.base_height_mm);
    let radius = p[0].hypot(p[1]);
    for feature in &program.features {
        if radius < f64::from(feature.radius_mm)
            && p[2] > f64::from(feature.z_mm[0])
            && p[2] < f64::from(feature.z_mm[1])
        {
            solid = feature.op == CylinderOp::Add;
        }
    }
    solid
}

fn assert_probe_pair(
    program: &CoaxialProgram,
    selection: &FaceSelection,
    solid_on_negative_normal: bool,
) {
    let eps = 1.0e-6;
    let before = std::array::from_fn(|i| selection.point_mm[i] - eps * selection.normal[i]);
    let after = std::array::from_fn(|i| selection.point_mm[i] + eps * selection.normal[i]);
    assert_eq!(oracle_material(program, before), solid_on_negative_normal);
    assert_eq!(oracle_material(program, after), !solid_on_negative_normal);
}

#[test]
fn wall_selections_are_fixed_and_report_provenance_or_disappearance() {
    let through = p(&[(CylinderOp::Cut, 4, [0, 8], "hole")]);
    let mut session = CoaxialSession::new(through.clone()).unwrap();
    let wall = selection([4., 0., 4.], [-1., 0., 0.], "hole");
    session.select(wall.clone()).unwrap();
    assert_probe_pair(&through, &wall, true);

    let renamed = p(&[(CylinderOp::Cut, 4, [0, 8], "bore")]);
    let state = session.replace(renamed).unwrap();
    assert_eq!(
        state.selection_events[0].status,
        SelectionStatus::ProvenanceChanged {
            feature: "bore".into()
        }
    );
    cold_check(&session);

    let translated = p(&[(CylinderOp::Cut, 5, [0, 8], "bore")]);
    assert_eq!(
        session.replace(translated).unwrap().selection_events[0].status,
        SelectionStatus::Disappeared
    );
    cold_check(&session);

    let restored = p(&[(CylinderOp::Cut, 4, [0, 8], "hole")]);
    assert!(matches!(
        session.replace(restored.clone()).unwrap().selection_events[0].status,
        SelectionStatus::Attached { .. }
    ));
    session.replace(through).unwrap();
    assert!(matches!(
        session.state().unwrap().selection_events[0].status,
        SelectionStatus::Attached { .. }
    ));
    cold_check(&session);
}

#[test]
fn additions_cuts_order_and_noop_rebuilds_match_cold_models() {
    let mut session = CoaxialSession::new(p(&[])).unwrap();
    let boss = p(&[(CylinderOp::Add, 4, [8, 12], "boss")]);
    session.replace(boss.clone()).unwrap();
    let outer = selection([4., 0., 10.], [1., 0., 0.], "boss");
    session.select(outer.clone()).unwrap();
    assert_probe_pair(&boss, &outer, true);

    let noop = p(&[
        (CylinderOp::Add, 4, [8, 12], "boss"),
        (CylinderOp::Cut, 2, [20, 24], "noop"),
    ]);
    let state = session.replace(noop.clone()).unwrap();
    assert!(matches!(
        state.selection_events[0].status,
        SelectionStatus::Attached { .. }
    ));
    let prior_volume = build_coaxial_program(&boss).unwrap();
    let prior_certificate = audit_coaxial_program(&boss, &prior_volume).unwrap();
    assert_eq!(state.ideal_volume_mm3, prior_certificate.ideal_volume_mm3);
    cold_check(&session);

    let added_then_cut = p(&[
        (CylinderOp::Add, 6, [8, 12], "boss"),
        (CylinderOp::Cut, 2, [8, 12], "bore"),
    ]);
    let mut ordered = CoaxialSession::new(added_then_cut.clone()).unwrap();
    let inner = selection([2., 0., 10.], [-1., 0., 0.], "bore");
    ordered.select(inner.clone()).unwrap();
    assert_probe_pair(&added_then_cut, &inner, true);
    let reversed = p(&[
        (CylinderOp::Cut, 2, [8, 12], "bore"),
        (CylinderOp::Add, 6, [8, 12], "boss"),
    ]);
    assert_eq!(
        ordered.replace(reversed).unwrap().selection_events[0].status,
        SelectionStatus::Disappeared
    );
    cold_check(&ordered);
}

#[test]
fn cavity_creation_deletion_and_edge_ambiguity_are_explicit() {
    let mut session = CoaxialSession::new(p(&[])).unwrap();
    let cavity = p(&[(CylinderOp::Cut, 3, [2, 6], "closed")]);
    let state = session.replace(cavity.clone()).unwrap();
    assert_eq!(state.cavities, 1);
    let wall = selection([3., 0., 4.], [-1., 0., 0.], "closed");
    session.select(wall.clone()).unwrap();
    assert_probe_pair(&cavity, &wall, true);
    assert_eq!(
        session.replace(p(&[])).unwrap().selection_events[0].status,
        SelectionStatus::Disappeared
    );

    let through = p(&[(CylinderOp::Cut, 4, [0, 8], "hole")]);
    let mut edge = CoaxialSession::new(through).unwrap();
    let at_edge = selection([4., 0., 0.], [-1., 0., 0.], "hole");
    assert_eq!(
        edge.select(at_edge).unwrap_err(),
        "coaxial_selection_ambiguous"
    );
}

#[test]
fn rejected_edits_keep_geometry_and_selection_unchanged() {
    let through = p(&[(CylinderOp::Cut, 4, [0, 8], "hole")]);
    let mut session = CoaxialSession::new(through).unwrap();
    session
        .select(selection([4., 0., 4.], [-1., 0., 0.], "hole"))
        .unwrap();
    let before = session.state().unwrap().clone();
    assert!(session
        .replace(p(&[(CylinderOp::Cut, 0, [0, 8], "invalid")]))
        .is_err());
    assert_eq!(*session.state().unwrap(), before);
    cold_check(&session);
}
