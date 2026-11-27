use solidcraft_trivariate::{conforming_bending, g2_cross_plate, json};

fn model() -> json::Json {
    let recipe = json::parse(include_str!(
        "../tests/jobs/crossed_partial_rise3.recipe.json"
    ))
    .unwrap();
    g2_cross_plate::boolean::report(&recipe, "union").unwrap()
}

fn number(value: &json::Json, key: &str) -> f64 {
    value.field(key).unwrap().number().unwrap()
}

#[test]
fn reproduces_the_frozen_two_span_c0_control() {
    let result = conforming_bending::analyze(&model(), 2, 2, 8, "C0").unwrap();
    assert_eq!(number(&result, "free_vector_dofs"), 4128.);
    assert!((number(&result, "mean_ux_mm") - 0.29729400005065354).abs() < 2e-7);
    assert!((number(&result, "energy_N_mm") - 5.978912667685339).abs() < 2e-6);
    assert!(number(&result, "relative_residual") < 2e-8);
}

#[test]
fn reproduces_the_two_span_c1_control() {
    let result = conforming_bending::analyze(&model(), 2, 2, 8, "C1").unwrap();
    assert_eq!(number(&result, "free_vector_dofs"), 2520.);
    assert!((number(&result, "mean_ux_mm") - 0.2896992614239085).abs() < 3e-7);
    assert!((number(&result, "energy_N_mm") - 5.826174035303023).abs() < 3e-6);
    assert!(number(&result, "relative_residual") < 2e-8);
}

#[test]
fn reproduces_the_equal_dof_three_span_c1_control() {
    let result = conforming_bending::analyze(&model(), 3, 3, 8, "C1").unwrap();
    assert_eq!(number(&result, "free_vector_dofs"), 4128.);
    assert!((number(&result, "mean_ux_mm") - 0.29817755452912864).abs() < 3e-7);
    assert!((number(&result, "energy_N_mm") - 5.9966819299746685).abs() < 3e-6);
    assert!(number(&result, "relative_residual") < 2e-8);
}

#[test]
fn reproduces_the_matched_iga_six_by_two_c0_row() {
    let result = conforming_bending::analyze(&model(), 6, 2, 8, "C0").unwrap();
    assert_eq!(number(&result, "free_vector_dofs"), 11112.);
    assert!((number(&result, "mean_ux_mm") - 0.3055843317926726).abs() < 4e-7);
    assert!(number(&result, "relative_residual") < 2e-8);
}

#[test]
fn rejects_a_moved_source_chart_and_unsupported_parameters() {
    let mut moved = model();
    let json::Json::Object(root) = &mut moved else {
        panic!("expected model object")
    };
    let json::Json::Array(blocks) = root.get_mut("blocks").unwrap() else {
        panic!("expected chart list")
    };
    let json::Json::Object(block) = &mut blocks[0] else {
        panic!("expected chart object")
    };
    let json::Json::Array(controls) = block.get_mut("controls").unwrap() else {
        panic!("expected controls")
    };
    let json::Json::Array(point) = &mut controls[0] else {
        panic!("expected control point")
    };
    point[0] = json::Json::Number("86".into());
    assert!(conforming_bending::analyze(&moved, 2, 2, 8, "C0").is_err());

    let original = model();
    assert!(conforming_bending::analyze(&original, 2, 2, 8, "C2").is_err());
    assert!(conforming_bending::analyze(&original, 0, 2, 8, "C0").is_err());
}

#[test]
fn rejects_a_replayed_union_with_opposite_rise() {
    let mut source = model().field("source_recipe").unwrap().clone();
    let json::Json::Object(fields) = &mut source else {
        unreachable!()
    };
    fields.insert("transverse_rise_mm".into(), json::Json::Number("-3".into()));
    let opposite = g2_cross_plate::boolean::report(&source, "union").unwrap();
    assert!(conforming_bending::analyze(&opposite, 2, 2, 8, "C1").is_err());
}
