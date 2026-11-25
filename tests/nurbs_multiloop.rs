use solidcraft_trivariate::{
    bracket::BracketBlock, incidence::enumerate_faces, json::parse, nurbs_multiloop::construct,
};

fn prospect() -> solidcraft_trivariate::json::Json {
    parse(include_str!("references/nurbs_multiloop.recipe.json")).unwrap()
}

#[test]
fn partitions_two_independent_profiles_across_global_z_events() {
    let recipe = prospect();
    let (blocks, faces) = construct(&recipe).unwrap();
    let incidence = enumerate_faces(&blocks).unwrap();
    assert_eq!(faces.len(), incidence.exterior.len());
    assert_eq!(
        faces
            .iter()
            .map(|(side, _)| side)
            .collect::<std::collections::BTreeSet<_>>(),
        incidence.exterior.iter().collect()
    );
    assert!(blocks.iter().any(|block| block.id.contains("_z0_")));
    assert!(blocks.iter().any(|block| block.id.contains("_z1_")));
    assert!(blocks.iter().all(|block| block.degrees[0] <= 3));
    for label in [
        "stock:xmin",
        "stock:xmax",
        "stock:ymin",
        "stock:ymax",
        "stock:zmin",
        "stock:zmax",
        "north-free-profile:wall_lower",
        "north-free-profile:wall_upper",
        "north-free-profile:end_lower",
        "south-curved-pocket:ceiling",
    ] {
        assert!(
            faces.iter().any(|(_, feature)| feature == label),
            "missing {label}"
        );
    }
}

#[test]
fn loop_order_and_names_do_not_select_geometry() {
    let original = prospect();
    let (expected, _) = construct(&original).unwrap();
    let mut loops = original.field("loops").unwrap().items().unwrap().to_vec();
    loops.reverse();
    let mut renamed = Vec::new();
    for (index, mut item) in loops.into_iter().enumerate() {
        if let solidcraft_trivariate::json::Json::Object(fields) = &mut item {
            fields.insert(
                "id".into(),
                solidcraft_trivariate::json::Json::String(format!("renamed_{index}")),
            );
        }
        renamed.push(item);
    }
    let mut fields = match original {
        solidcraft_trivariate::json::Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    fields.insert(
        "loops".into(),
        solidcraft_trivariate::json::Json::Array(renamed),
    );
    let changed = solidcraft_trivariate::json::Json::Object(fields);
    let (actual, _) = construct(&changed).unwrap();
    let normalize = |blocks: Vec<BracketBlock>| {
        blocks
            .into_iter()
            .map(|mut b| {
                b.id.clear();
                b
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(normalize(expected), normalize(actual));
}

#[test]
fn rejects_stationary_profile_without_weakening_monotonicity() {
    let stationary = parse(include_str!(
        "counterexamples/nurbs_multiloop_stationary.json"
    ))
    .unwrap();
    assert!(construct(&stationary).is_err());
}

#[test]
fn rejects_touching_and_overlapping_loop_projections() {
    let mut recipe = prospect();
    let mut fields = match recipe {
        solidcraft_trivariate::json::Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    let mut loops = fields.remove("loops").unwrap().items().unwrap().to_vec();
    let second = match loops[1].clone() {
        solidcraft_trivariate::json::Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    let mut second = second;
    second.insert("lower".into(), loops[0].field("lower").unwrap().clone());
    second.insert("upper".into(), loops[0].field("upper").unwrap().clone());
    loops[1] = solidcraft_trivariate::json::Json::Object(second);
    fields.insert(
        "loops".into(),
        solidcraft_trivariate::json::Json::Array(loops),
    );
    recipe = solidcraft_trivariate::json::Json::Object(fields);
    assert!(construct(&recipe).is_err());
}

#[test]
fn rejects_contours_touching_the_stock_and_degenerate_axial_intervals() {
    let original = prospect();
    let mut fields = match original.clone() {
        solidcraft_trivariate::json::Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    let mut loops = fields.remove("loops").unwrap().items().unwrap().to_vec();
    let mut lower = match loops[0].field("lower").unwrap().clone() {
        solidcraft_trivariate::json::Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    if let solidcraft_trivariate::json::Json::Array(controls) = lower.get_mut("controls").unwrap() {
        if let solidcraft_trivariate::json::Json::Array(first) = &mut controls[0] {
            first[0] = solidcraft_trivariate::json::Json::Number("0".into());
        }
    }
    let mut loop_fields = match loops[0].clone() {
        solidcraft_trivariate::json::Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    loop_fields.insert(
        "lower".into(),
        solidcraft_trivariate::json::Json::Object(lower),
    );
    loops[0] = solidcraft_trivariate::json::Json::Object(loop_fields);
    fields.insert(
        "loops".into(),
        solidcraft_trivariate::json::Json::Array(loops),
    );
    assert!(construct(&solidcraft_trivariate::json::Json::Object(fields)).is_err());

    let mut fields = match original {
        solidcraft_trivariate::json::Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    let mut loops = fields.remove("loops").unwrap().items().unwrap().to_vec();
    let mut loop_fields = match loops[0].clone() {
        solidcraft_trivariate::json::Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    loop_fields.insert(
        "z_interval_mm".into(),
        solidcraft_trivariate::json::Json::Array(vec![
            solidcraft_trivariate::json::Json::Number("4".into()),
            solidcraft_trivariate::json::Json::Number("4".into()),
        ]),
    );
    loops[0] = solidcraft_trivariate::json::Json::Object(loop_fields);
    fields.insert(
        "loops".into(),
        solidcraft_trivariate::json::Json::Array(loops),
    );
    assert!(construct(&solidcraft_trivariate::json::Json::Object(fields)).is_err());
}

#[test]
fn labels_both_caps_of_an_internal_pocket() {
    let original = prospect();
    let mut fields = match original {
        solidcraft_trivariate::json::Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    let mut loops = fields.remove("loops").unwrap().items().unwrap().to_vec();
    let mut loop_fields = match loops[0].clone() {
        solidcraft_trivariate::json::Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    loop_fields.insert(
        "z_interval_mm".into(),
        solidcraft_trivariate::json::Json::Array(vec![
            solidcraft_trivariate::json::Json::Number("4".into()),
            solidcraft_trivariate::json::Json::Number("12".into()),
        ]),
    );
    loops[0] = solidcraft_trivariate::json::Json::Object(loop_fields);
    fields.insert(
        "loops".into(),
        solidcraft_trivariate::json::Json::Array(loops),
    );
    let recipe = solidcraft_trivariate::json::Json::Object(fields);
    let (_, faces) = construct(&recipe).unwrap();
    assert!(faces
        .iter()
        .any(|(_, name)| name == "north-free-profile:floor"));
    assert!(faces
        .iter()
        .any(|(_, name)| name == "north-free-profile:ceiling"));
}
