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

fn curved_pair(degree: usize) -> String {
    let x = if degree == 2 { vec![0, 0, 1] } else { vec![0, 0, 0, 1] };
    let y: Vec<_> = (0..=degree).map(|i| i as i32).collect();
    let controls = |which: char| {
        let mut points = Vec::new();
        for i in 0..=degree {
            for v in 0..=1 {
                for z in 0..=1 {
                    let offset = match (which, v) {
                        ('a', 0) => -10,
                        ('a', 1) => 0,
                        ('b', 0) => 0,
                        ('b', 1) => 10,
                        _ => unreachable!(),
                    };
                    points.push(format!("[{}, {}, {z}, 1]", x[i] + offset, y[i]));
                }
            }
        }
        points.join(",")
    };
    format!(r#"{{"blocks":[
      {{"id":"a","degrees":[{degree},1,1],"controls":[{}]}},
      {{"id":"b","degrees":[{degree},1,1],"controls":[{}]}}
    ],"topology":{{"interfaces":[{{"first":["a","v1"],"second":["b","v0"],
      "orientation":{{"swap":false,"flip_first":false,"flip_second":false}}}}]}}}}"#,
      controls('a'), controls('b'))
}

fn curved_pair_on_first_side(degree: usize) -> String {
    let x = if degree == 2 { vec![0, 0, 1] } else { vec![0, 0, 0, 1] };
    let y: Vec<_> = (0..=degree).map(|i| (i * degree) as i32).collect();
    let controls = |which: char| {
        let mut points = Vec::new();
        for i in 0..=degree {
            for j in 0..=degree {
                for z in 0..=1 {
                    let offset = if which == 'a' { 0 } else { 10 * degree as i32 };
                    points.push(format!("[{}, {}, {z}, 1]", x[j] * degree as i32 + offset - 10 * i as i32, y[j]));
                }
            }
        }
        points.join(",")
    };
    format!(r#"{{"blocks":[
      {{"id":"a","degrees":[{degree},{degree},1],"controls":[{}]}},
      {{"id":"b","degrees":[{degree},{degree},1],"controls":[{}]}}
    ],"topology":{{"interfaces":[{{"first":["a","u0"],"second":["b","u1"],
      "orientation":{{"swap":false,"flip_first":false,"flip_second":false}}}}]}}}}"#,
      controls('a'), controls('b'))
}

#[test]
fn exact_curved_gluing_supports_quadratic_and_cubic_edges() {
    for degree in [2, 3] {
        for input in [curved_pair(degree), curved_pair_on_first_side(degree)] {
            let model = parse(&input).unwrap();
            let report = classify(&model).unwrap();
            let row = &report.items().unwrap()[0];
            let proof = row.field("as_g1").unwrap();
            assert_eq!(proof.field("status").unwrap().text().unwrap(), "as_g1");
            assert!(proof.field("face_association_verified").unwrap().bool().unwrap());
            assert!(proof.field("syzygy_verified").unwrap().bool().unwrap());
            assert!(proof.field("beta_split_verified").unwrap().bool().unwrap());
        }
    }
}

#[test]
fn cubic_coprime_quadratic_transverse_factors_are_disproven_as_g1() {
    // Both alpha factors have genuine quadratic degree and gcd 1. Their product
    // cannot be represented by the required linear beta-side factors.
    let model = parse(include_str!("references/as_g1_cubic_coprime_factors.json")).unwrap();
    let report = classify(&model).unwrap();
    let proof = report.items().unwrap()[0].field("as_g1").unwrap();
    assert_eq!(proof.field("status").unwrap().text().unwrap(), "not_as_g1");
    assert_eq!(proof.field("reason").unwrap().text().unwrap(), "no_polynomial_beta_split_within_degree");
}

#[test]
fn curved_trace_mismatch_with_matching_endpoints_is_not_a_glued_face() {
    let original = curved_pair(2);
    let (head, tail) = original.split_once("\"id\":\"b\"").unwrap();
    let altered = format!("{head}\"id\":\"b\"{}", tail.replace("[0, 1, 0, 1],[0, 1, 1, 1]", "[0, 2, 0, 1],[0, 2, 1, 1]"));
    let model = parse(&altered).unwrap();
    let report = classify(&model).unwrap();
    assert_eq!(report.items().unwrap()[0].field("as_g1").unwrap().field("status").unwrap().text().unwrap(), "unsupported_trace");
}

#[test]
fn lateral_face_association_checks_top_and_bottom_curved_traces() {
    let original = curved_pair(2);
    let (head, tail) = original.split_once("\"id\":\"b\"").unwrap();
    let unequal_height = format!("{head}\"id\":\"b\"{}", tail.replace(", 1, 1]", ", 2, 1]"));
    let model = parse(&unequal_height).unwrap();
    let report = classify(&model).unwrap();
    let proof = report.items().unwrap()[0].field("as_g1").unwrap();
    assert_eq!(proof.field("status").unwrap().text().unwrap(), "unsupported_frame");
}

#[test]
fn published_biquadratic_pocket_matches_independent_symbolic_gluing_corpus() {
    let model = parse(include_str!("../docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json")).unwrap();
    let expected = parse(include_str!("../docs/results/curved-pocket-fields/symbolic-gluing.json")).unwrap();
    let rows = classify(&model).unwrap();
    let lateral: Vec<_> = rows.items().unwrap().iter()
        .filter(|row| row.field("as_g1").unwrap().field("status").unwrap().text().unwrap() != "nonlateral_interface")
        .collect();
    assert_eq!(lateral.len(), expected.items().unwrap().len());
    let as_face = |v: &solidcraft_trivariate::json::Json| -> Vec<Vec<String>> {
        [v.field("first").unwrap(), v.field("second").unwrap()].into_iter()
            .map(|face| face.items().unwrap().iter().map(|s| s.text().unwrap().to_string()).collect()).collect()
    };
    let mut got = Vec::new();
    for row in lateral {
        let actual_faces: Vec<Vec<String>> = row.field("faces").unwrap().items().unwrap().iter()
            .map(|f| f.items().unwrap().iter().map(|x| x.text().unwrap().to_string()).collect()).collect();
        let reference = expected.items().unwrap().iter().find(|item| as_face(item) == actual_faces).expect("independent face association");
        let accepted = row.field("as_g1").unwrap().field("status").unwrap().text().unwrap() == "as_g1";
        assert_eq!(accepted, reference.field("as_g1").unwrap().bool().unwrap(), "face pair {actual_faces:?}");
        assert!(row.field("as_g1").unwrap().field("face_association_verified").unwrap().bool().unwrap());
        got.push((actual_faces, accepted));
    }
    let positives = got.iter().filter(|(_, ok)| *ok).count();
    let negatives = got.len() - positives;
    eprintln!("biquadratic pocket exact AS-G1 census: {positives} admitted, {negatives} disproven; 4 nonlateral links excluded");
    assert_eq!((positives, negatives), (8, 4));
}
