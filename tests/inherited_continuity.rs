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
    let matched: Vec<_> = field
        .interface_evidence
        .iter()
        .filter(|x| {
            x.field("basis_evidence").and_then(Json::text).ok()
                == Some("exact_proportional_complete_homogeneous_normal_jet")
        })
        .collect();
    assert_eq!(matched.len(), 4);
    for join in matched {
        assert_eq!(
            join.field("normal_interval_ratio").unwrap(),
            &parse("[1,1]").unwrap()
        );
    }
}

#[test]
fn automatic_inherited_plan_reuses_only_the_qualified_radial_family() {
    let model = parse(include_str!(
        "../docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"
    ))
    .unwrap();
    let package = solidcraft_trivariate::inherited_space::compile_automatic(&model, 2).unwrap();
    let audit = solidcraft_trivariate::inherited_space::audit_package(&package).unwrap();
    assert_eq!(audit.field("status").unwrap().text().unwrap(), "valid");
    let field = package.field("field").unwrap();
    assert_eq!(
        field
            .field("continuity_selection")
            .unwrap()
            .field("selected_degree")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        2
    );
    assert_eq!(
        field
            .field("continuity_selection")
            .unwrap()
            .field("plane_degree")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        3
    );
    assert_eq!(
        field
            .field("continuity_selection")
            .unwrap()
            .field("qualified_radial_join_count")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        8
    );
    assert_eq!(
        field
            .field("transverse_c1_join_count")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        8
    );
    let interfaces = field.field("interfaces").unwrap().items().unwrap();
    assert_eq!(
        interfaces
            .iter()
            .filter(|x| x.field("role").unwrap().text().unwrap() == "distinct_support_join")
            .count(),
        8
    );
    assert!(interfaces.iter().all(|x| x.field("reason").is_ok()));
}

#[test]
fn automatic_inheritance_keeps_parent_c2_but_generic_lateral_joins_at_c0() {
    let recipe = parse(include_str!("references/nurbs_strip.recipe.json")).unwrap();
    let strip = solidcraft_trivariate::interchange::build(&recipe).unwrap();
    let strip_package =
        solidcraft_trivariate::inherited_space::compile_automatic(&strip, 2).unwrap();
    let strip_field = strip_package.field("field").unwrap();
    assert_eq!(strip_field.field("inherited_join_count").unwrap().number::<usize>().unwrap(), 1);
    let inherited = strip_field.field("interfaces").unwrap().items().unwrap().iter()
        .find(|v| v.field("role").unwrap().text().unwrap() == "common_parent_subdivision").unwrap();
    assert_eq!(inherited.field("field_order").unwrap().number::<usize>().unwrap(), 2);

    let model = cylinder();
    let generic = solidcraft_trivariate::inherited_space::compile_automatic(&model, 2).unwrap();
    let interfaces = generic.field("field").unwrap().field("interfaces").unwrap().items().unwrap();
    assert!(interfaces.iter().any(|v| {
        let first = v.field("first").unwrap().items().unwrap();
        first[1].text().unwrap().starts_with('u') || first[1].text().unwrap().starts_with('v')
    }));
    assert!(interfaces.iter().filter(|v| {
        let first = v.field("first").unwrap().items().unwrap();
        first[1].text().unwrap().starts_with('u') || first[1].text().unwrap().starts_with('v')
    }).all(|v| v.field("role").unwrap().text().unwrap() == "whole_face_join"));
}

#[test]
fn automatic_raw_load_rejects_declared_contact_and_unsupported_material_fields() {
    let model = pocket();
    let mut job = parse(include_str!("references/s3_inherited_stress.job.json")).unwrap();
    if let Json::Object(fields) = &mut job {
        fields.insert("clamp_feature".into(), Json::String("start".into()));
    }
    if let Json::Object(fields) = &mut job {
        fields.insert(
            "contact".into(),
            parse(r#"{"feature":"pocket/wall"}"#).unwrap(),
        );
    }
    assert_eq!(
        solidcraft_trivariate::mechanical_job::solve_job(&model, &job).unwrap_err(),
        "contact mechanics is unsupported by spline_feature_load_v1"
    );
    if let Json::Object(fields) = &mut job {
        fields.remove("contact");
        if let Json::Object(material) = fields.get_mut("material").unwrap() {
            material.insert("density_kg_per_mm3".into(), Json::Number("1".into()));
        }
    }
    assert_eq!(
        solidcraft_trivariate::mechanical_job::solve_job(&model, &job).unwrap_err(),
        "unexpected object fields"
    );
}

#[test]
fn manual_inherited_exports_keep_their_versioned_document_shape() {
    let recipe = parse(include_str!("references/g2_cross_sketches.recipe.json")).unwrap();
    let model = solidcraft_trivariate::interchange::build(&recipe).unwrap();
    let plan = parse(include_str!("references/matched_chart_jets.plan.json")).unwrap();
    let package = solidcraft_trivariate::inherited_space::compile(&model, &plan).unwrap();
    assert!(package
        .field("field")
        .unwrap()
        .field("continuity_selection")
        .is_err());
    assert!(package
        .field("field")
        .unwrap()
        .field("interfaces")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .all(|x| x.field("reason").is_err()));
    let interfaces = package
        .field("field")
        .unwrap()
        .field("interfaces")
        .unwrap()
        .items()
        .unwrap();
    let matched: Vec<_> = interfaces
        .iter()
        .filter(|x| {
            x.field("basis_evidence").and_then(Json::text).ok()
                == Some("exact_proportional_complete_homogeneous_normal_jet")
        })
        .collect();
    assert_eq!(matched.len(), 4);
    assert!(matched.iter().all(|x| {
        x.field("normal_interval_ratio").unwrap() == &parse("[1,1]").unwrap()
    }));
}
