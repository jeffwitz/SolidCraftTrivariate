use solidcraft_trivariate::{json::parse, transverse_transitions::classify};

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
