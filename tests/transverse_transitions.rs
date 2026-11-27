use solidcraft_trivariate::{json::parse, transverse_transitions::classify};

const MODEL_TEMPLATE: &str = r#"{
  "blocks":[
    {"id":"a","degrees":[1,1,1],"controls":[[0,0,0,1],[0,0,1,1],[0,1,0,1],[0,1,1,1],[1,0,0,1],[1,0,1,1],[1,1,0,1],[1,1,1,1]]},
    {"id":"b","degrees":[1,1,1],"controls":B_CONTROLS}
  ],
  "topology":{"interfaces":[{"first":["a","u1"],"second":["b","u0"],"orientation":{"swap":false,"flip_first":false,"flip_second":false}}]}
}"#;

fn as_g1_status(controls: &str) -> (String, String) {
    let model = parse(&MODEL_TEMPLATE.replace("B_CONTROLS", controls)).unwrap();
    let report = classify(&model).unwrap();
    let row = &report.items().unwrap()[0];
    (
        row.field("status").unwrap().text().unwrap().to_string(),
        row.field("as_g1")
            .unwrap()
            .field("status")
            .unwrap()
            .text()
            .unwrap()
            .to_string(),
    )
}

#[test]
fn affine_first_jet_transition_matches_rational_reference() {
    let model = parse(r#"{
      "blocks":[
        {"id":"a","degrees":[1,1,1],"controls":[[0,0,0,1],[0,0,1,1],[0,1,0,1],[0,1,1,1],[1,0,0,1],[1,0,1,1],[1,1,0,1],[1,1,1,1]]},
        {"id":"b","degrees":[1,1,1],"controls":[[1,0,0,1],[1,0,1,1],[1,1,0,1],[1,1,1,1],[2,0,0,1],[2,0,1,1],[2,1,0,1],[2,1,1,1]]}
      ],
      "topology":{"interfaces":[{"first":["a","u1"],"second":["b","u0"],"orientation":{"swap":false,"flip_first":false,"flip_second":false}}]}
    }"#).unwrap();
    let report = classify(&model).unwrap();
    let row = &report.items().unwrap()[0];
    assert_eq!(
        row.field("status").unwrap().text().unwrap(),
        "constant_transition"
    );
    assert_eq!(row.field("alpha").unwrap().to_string(), r#"["-1","0"]"#);
    assert_eq!(row.field("beta").unwrap().to_string(), r#"["0","0"]"#);
    assert_eq!(
        row.field("reversed_status").unwrap().text().unwrap(),
        "constant_transition"
    );
    assert_eq!(
        row.field("as_g1")
            .unwrap()
            .field("status")
            .unwrap()
            .text()
            .unwrap(),
        "as_g1"
    );
    assert_eq!(
        row.field("as_g1")
            .unwrap()
            .field("gluing_factors")
            .unwrap()
            .field("alpha_first")
            .unwrap()
            .to_string(),
        r#"["1"]"#
    );
    assert!(
        row.field("as_g1")
            .unwrap()
            .field("beta_split_verified")
            .unwrap()
            .bool()
            .unwrap()
    );
    assert!(
        row.field("as_g1")
            .unwrap()
            .field("syzygy_verified")
            .unwrap()
            .bool()
            .unwrap()
    );
}

#[test]
fn reversed_bilinear_orientation_is_not_certified_as_g1() {
    let model = parse(r#"{
      "blocks":[
        {"id":"a","degrees":[1,1,1],"controls":[[0,0,0,1],[0,0,1,1],[0,1,0,1],[0,1,1,1],[1,0,0,1],[1,0,1,1],[1,1,0,1],[1,1,1,1]]},
        {"id":"b","degrees":[1,1,1],"controls":[[1,0,0,1],[1,0,1,1],[1,1,0,1],[1,1,1,1],[0,0,0,1],[0,0,1,1],[0,1,0,1],[0,1,1,1]]}
      ],
      "topology":{"interfaces":[{"first":["a","u1"],"second":["b","u0"],"orientation":{"swap":false,"flip_first":false,"flip_second":false}}]}
    }"#).unwrap();
    let result = classify(&model).unwrap();
    let status = result.items().unwrap()[0]
        .field("as_g1")
        .unwrap()
        .field("status")
        .unwrap()
        .text()
        .unwrap();
    assert_eq!(status, "not_as_g1");
    assert!(
        !result.items().unwrap()[0]
            .field("as_g1")
            .unwrap()
            .field("orientable")
            .unwrap()
            .bool()
            .unwrap()
    );
}

#[test]
fn transition_classifier_rejects_unrepresentable_exact_difference() {
    let model = parse(
        r#"{
      "blocks":[
        {"id":"a","degrees":[1,1,1],"controls":[[0,0,-170141183460469231731687303715884105728,1],[0,0,170141183460469231731687303715884105727,1],[0,1,0,1],[0,1,1,1],[1,0,0,1],[1,0,1,1],[1,1,0,1],[1,1,1,1]]},
        {"id":"b","degrees":[1,1,1],"controls":[[1,0,0,1],[1,0,1,1],[1,1,0,1],[1,1,1,1],[2,0,0,1],[2,0,1,1],[2,1,0,1],[2,1,1,1]]}
      ],
      "topology":{"interfaces":[{"first":["a","u1"],"second":["b","u0"],"orientation":{"swap":false,"flip_first":false,"flip_second":false}}]}
    }"#,
    )
    .unwrap();
    assert_eq!(classify(&model).unwrap_err(), "rational overflow");
}

#[test]
fn equal_jets_do_not_certify_disjoint_faces() {
    let controls =
        r#"[[3,0,0,1],[3,0,1,1],[3,1,0,1],[3,1,1,1],[4,0,0,1],[4,0,1,1],[4,1,0,1],[4,1,1,1]]"#;
    let (transition, as_g1) = as_g1_status(controls);
    assert_eq!(transition, "constant_transition");
    assert_eq!(as_g1, "unsupported_trace");
}

#[test]
fn malformed_control_net_is_reported_without_panicking() {
    let controls = r#"[[2,0,0,1],[2,0,1,1],[2,1,0,1],[2,1,1,1],[3,0,0,1],[3,0,1,1],[3,1,0,1]]"#;
    let (transition, as_g1) = as_g1_status(controls);
    assert_eq!(transition, "unsupported_transition_frame");
    assert_eq!(as_g1, "unsupported_frame");
}

#[test]
fn nonconstant_projective_weights_are_not_certified() {
    let controls =
        r#"[[4,0,0,2],[4,0,2,2],[2,1,0,1],[2,1,1,1],[3,0,0,1],[3,0,1,1],[3,1,0,1],[3,1,1,1]]"#;
    let (_, as_g1) = as_g1_status(controls);
    assert_eq!(as_g1, "unsupported_projective_weights");
}

#[test]
fn nonplanar_section_is_not_certified_as_an_extruded_section() {
    let controls =
        r#"[[2,0,1,1],[2,0,2,1],[2,1,0,1],[2,1,1,1],[3,0,1,1],[3,0,2,1],[3,1,0,1],[3,1,1,1]]"#;
    let (_, as_g1) = as_g1_status(controls);
    assert_eq!(as_g1, "unsupported_extrusion_support");
}
