use solidcraft_trivariate::{
    inherited_space::{audit_package, build_space, compile},
    json::{parse, Json},
};
fn cylinder() -> Json {
    parse(include_str!(
        "../docs/site/box-cylinder/tall/cylinder.model.json"
    ))
    .unwrap()
}
fn pocket() -> Json {
    parse(include_str!(
        "../docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"
    ))
    .unwrap()
}

#[test]
fn layer_continuity_preserves_physical_affine_gradients() {
    for (model, baseline, restored, joins) in [(cylinder(), 119, 85, 10), (pocket(), 123, 99, 4)] {
        let space = build_space(&model, 2, 1, false, &[]).unwrap();
        assert_eq!(
            (space.c0_ndof, space.ndof, space.common_parent_faces),
            (baseline, restored, joins)
        );
        for cell in &space.cells {
            let (geo, basis) = cell.sample([0.19, 0.71, 0.43]).unwrap();
            let sum: f64 = basis.values.iter().sum();
            assert!((sum - 1.).abs() < 1e-12);
            for i in 0..3 {
                let position: f64 = basis
                    .dofs
                    .iter()
                    .zip(&basis.values)
                    .map(|(&j, r)| space.control_xyz[j][i] * r)
                    .sum();
                assert!((position - geo.position[i]).abs() < 1e-10);
                for k in 0..3 {
                    let derivative: f64 = basis
                        .dofs
                        .iter()
                        .zip(&basis.physical_gradients)
                        .map(|(&j, g)| space.control_xyz[j][i] * g[k])
                        .sum();
                    assert!((derivative - f64::from(i == k)).abs() < 1e-10);
                }
            }
        }
    }
}
#[test]
fn altered_field_is_rejected_and_unknown_joins_are_not_consumed() {
    let model = cylinder();
    let plan=parse(r#"{"schema":"inherited_field_plan_v1","axial_degree":1,"axial_continuity":1,"restore_strip":false,"transverse_interfaces":[]}"#).unwrap();
    let mut doc = compile(&model, &plan).unwrap();
    audit_package(&doc).unwrap();
    if let Json::Object(all) = &mut doc {
        if let Json::Object(f) = all.get_mut("field").unwrap() {
            f.insert("scalar_dofs".into(), Json::Number("1".into()));
        }
    }
    assert!(audit_package(&doc).is_err());
    assert!(build_space(&model, 1, 2, false, &[]).is_err());
    let unknown = ("unknown".to_owned(), "disk_lower_core".to_owned());
    assert!(build_space(&model, 1, 1, false, &[unknown]).is_err());
}
#[test]
fn failed_compile_preserves_a_previous_file() {
    let out =
        std::env::temp_dir().join(format!("sct-inherited-{}-atomic.json", std::process::id()));
    std::fs::write(&out, b"previous state").unwrap();
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
        .args([
            "export-inherited-space",
            "docs/site/box-cylinder/tall/cylinder.model.json",
            "tests/references/affine_pair.json",
        ])
        .arg(&out)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(&out).unwrap(), b"previous state");
    std::fs::remove_file(out).unwrap();
}

#[test]
fn one_affine_stock_is_a_complete_parent_and_a_tunnel_is_obstructed() {
    use solidcraft_trivariate::spline_support::{eligibility_report, recover_complete_parent};
    let box_model = parse(include_str!(
        "../docs/site/box-cylinder/tall/box.model.json"
    ))
    .unwrap();
    let parent = recover_complete_parent(&box_model).unwrap();
    assert_eq!(parent.control_shape, [2; 3]);
    assert_eq!(parent.source_spans.len(), 1);
    let tunnel = parse(include_str!(
        "../docs/site/box-cylinder/tall/a-minus-b.model.json"
    ))
    .unwrap();
    assert_eq!(
        eligibility_report(&tunnel)
            .unwrap()
            .field("common_parent_obstruction")
            .unwrap()
            .field("kind")
            .unwrap()
            .text()
            .unwrap(),
        "exact_through_tunnel_topology"
    );
}

#[test]
fn inherited_strip_uses_true_intervals_on_decreasing_nonuniform_sketch() {
    let mut recipe = parse(include_str!("references/nurbs_strip.recipe.json")).unwrap();
    let Json::Object(parts) = &mut recipe else {
        unreachable!()
    };
    for name in ["lower", "upper"] {
        let Some(Json::Object(curve)) = parts.get_mut(name) else {
            unreachable!()
        };
        let Some(Json::Array(controls)) = curve.get_mut("controls") else {
            unreachable!()
        };
        controls.reverse();
        curve.insert("knots".into(), parse("[-5,-5,-5,-5,0,4,4,4,4]").unwrap());
    }
    let model = solidcraft_trivariate::interchange::build(&recipe).unwrap();
    let field = build_space(&model, 1, 1, true, &[]).unwrap();
    assert_eq!(
        (field.c0_ndof, field.ndof, field.common_parent_faces),
        (28, 20, 1)
    );
}

#[test]
fn source_g2_sketch_only_restores_lateral_faces_with_complete_normal_jets() {
    let recipe = parse(include_str!("references/g2_cross_sketches.recipe.json")).unwrap();
    let model = solidcraft_trivariate::interchange::build(&recipe).unwrap();
    let plan = parse(include_str!("references/matched_chart_jets.plan.json")).unwrap();
    let field = solidcraft_trivariate::inherited_space::space_from_plan(&model, &plan).unwrap();
    assert_eq!(
        (
            field.common_parent_faces,
            field.transverse_faces,
            field.retained_faces
        ),
        (16, 4, 26)
    );
}
