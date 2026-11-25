use solidcraft_trivariate::{
    interchange,
    json::{parse, Json},
};

#[test]
fn unequal_leg_filleted_l_uses_certified_rectangular_local_chart() {
    let recipe = parse(include_str!("references/unequal_fillet_l.recipe.json")).unwrap();
    let model = interchange::build(&recipe).expect("unequal leg widths should be supported");
    let plan = model.field("partition_plan").unwrap();
    assert_eq!(
        plan.field("method").unwrap().text().unwrap(),
        "axis_planes_and_disjoint_local_sweeps"
    );
    assert_eq!(
        plan.field("fallback_reason").unwrap().text().unwrap(),
        "compact_bracket_unequal_leg_thickness"
    );
    let audit = interchange::audit(&model).expect("serialized partition must replay");
    assert_eq!(audit.field("status").unwrap().text().unwrap(), "valid");
    assert!(
        audit
            .field("determinant_lower_bound_mm3")
            .unwrap()
            .number::<f64>()
            .unwrap()
            > 0.0
    );
    let volume = audit.field("volume_exact").unwrap();
    assert_eq!(
        volume
            .field("twice_rational_part_mm3")
            .unwrap()
            .number::<i64>()
            .unwrap(),
        61_824
    );
    assert_eq!(
        volume.field("pi_coefficient_mm3").unwrap().items().unwrap()[0]
            .number::<i64>()
            .unwrap(),
        -192
    );
    let topology = audit
        .field("shared_faces")
        .unwrap()
        .number::<usize>()
        .unwrap();
    assert!(topology > 0);

    let mut legacy = model.clone();
    let Json::Object(model_fields) = &mut legacy else {
        unreachable!()
    };
    let Json::Object(plan_fields) = model_fields.get_mut("partition_plan").unwrap() else {
        unreachable!()
    };
    plan_fields.remove("fallback_reason");
    assert!(interchange::audit(&legacy).is_ok());

    let mut altered = model;
    let Json::Object(model_fields) = &mut altered else {
        unreachable!()
    };
    let Json::Object(plan_fields) = model_fields.get_mut("partition_plan").unwrap() else {
        unreachable!()
    };
    plan_fields.insert(
        "fallback_reason".into(),
        Json::String("different_reason".into()),
    );
    assert!(interchange::audit(&altered).is_err());
}

#[test]
fn unequal_arms_preserve_both_transverse_bores_on_grid_fallback() {
    let recipe = parse(include_str!(
        "references/unequal_fillet_cross_bores.recipe.json"
    ))
    .unwrap();
    let model = interchange::build(&recipe).unwrap();
    let report = interchange::audit(&model).unwrap();
    for (field, expected) in [("through_bores", 2), ("blocks", 48), ("shared_faces", 74)] {
        assert_eq!(
            report.field(field).unwrap().number::<usize>().unwrap(),
            expected
        );
    }
    let plan = model.field("partition_plan").unwrap();
    assert_eq!(
        plan.field("fallback_reason").unwrap().text().unwrap(),
        "compact_bracket_unequal_leg_thickness"
    );
}
