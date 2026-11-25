use solidcraft_trivariate::{
    interchange,
    json::{parse, Json},
    solidcraft_sketch_features::import_sketch_timeline,
};

const ORIGINAL: &str = r#"{
  "about":{"design":"original","format":"solidcraft/66","replay":"solidcraft-cli eval timeline.json","version":"0.48.1"},
  "commands":[
    {"command":"timeline.add","params":{"feature":{"id":1,"name":"Profile","plane":{"name":"XY","type":"origin"},"sketch":{"constraints":[{"id":"k1","p":0,"q":1,"type":"coincident"}],"counters":{"k":1,"l":6},"curves":[{"a":1,"b":2,"construction":false,"id":"l1","type":"line"},{"a":2,"b":3,"construction":false,"id":"l2","type":"line"},{"a":3,"b":4,"construction":false,"id":"l3","type":"line"},{"a":4,"b":5,"construction":false,"id":"l4","type":"line"},{"a":5,"b":6,"construction":false,"id":"l5","type":"line"},{"a":6,"b":1,"construction":false,"id":"l6","type":"line"}],"points":[{"fixed":true,"id":"origin","pos":[0.0,0.0]},{"fixed":false,"id":"l1.start","pos":[0.0,0.0]},{"fixed":false,"id":"l1.end","pos":[60.0,0.0]},{"fixed":false,"id":"l2.end","pos":[60.0,8.0]},{"fixed":false,"id":"l3.end","pos":[8.0,8.0]},{"fixed":false,"id":"l4.end","pos":[8.0,50.0]},{"fixed":false,"id":"l5.end","pos":[0.0,50.0]}]},"suppressed":false,"type":"sketch"}}},
    {"command":"timeline.add","params":{"feature":{"body_names":["Bracket"],"extent":{"convention":"total_distance","direction":"positive","distance":"40"},"id":2,"name":"Extrude1","operation":"new_body","param_names":["d1"],"profiles":{"type":"all"},"sketch":1,"suppressed":false,"type":"extrude"}}},
    {"command":"timeline.add","params":{"feature":{"edge_names":["F2:side:l3|F2:side:l4"],"edge_refs":[{"dir":[0.0,0.0,1.0],"max":[60.0,50.0,40.0],"min":[0.0,0.0,0.0]}],"edges":[[8.0,8.0,20.0]],"id":3,"name":"Fillet1","param_names":["d2"],"radius":"6","suppressed":false,"type":"fillet"}}},
    {"command":"timeline.add","params":{"feature":{"diameter":"8","direction":[-0.0,-1.0,0.0],"face_names":["F2:side:l3"],"hole":{"type":"simple"},"id":4,"name":"Hole1","param_names":["d3"],"position":[40.0,8.0,20.0],"suppressed":false,"type":"hole"}}},
    {"command":"timeline.add","params":{"feature":{"diameter":"8","direction":[-1.0,-0.0,-0.0],"face_names":["F2:side:l4"],"hole":{"type":"simple"},"id":5,"name":"Hole2","param_names":["d4"],"position":[8.0,32.0,20.0],"suppressed":false,"type":"hole"}}}
  ]
}"#;

fn original() -> Json {
    parse(ORIGINAL).unwrap()
}

fn object_mut(value: &mut Json) -> &mut std::collections::BTreeMap<String, Json> {
    let Json::Object(fields) = value else {
        panic!("expected object")
    };
    fields
}

fn array_mut(value: &mut Json) -> &mut Vec<Json> {
    let Json::Array(items) = value else {
        panic!("expected array")
    };
    items
}

fn set(value: &mut Json, key: &str, replacement: Json) {
    object_mut(value).insert(key.to_owned(), replacement);
}

fn feature_mut(script: &mut Json, command: usize) -> &mut Json {
    let commands = array_mut(object_mut(script).get_mut("commands").unwrap());
    let params = object_mut(commands.get_mut(command).unwrap())
        .get_mut("params")
        .unwrap();
    object_mut(params).get_mut("feature").unwrap()
}

#[test]
fn genuine_solidcraft_l_history_becomes_audited_48_block_sketch_model() {
    let model = import_sketch_timeline(&original()).unwrap();
    assert_eq!(
        model.field("schema").unwrap().text().unwrap(),
        "orthogonal_sketch_features_v1"
    );
    assert_eq!(model.field("blocks").unwrap().items().unwrap().len(), 48);
    assert_eq!(
        interchange::audit(&model)
            .unwrap()
            .field("status")
            .unwrap()
            .text()
            .unwrap(),
        "valid"
    );
    let report = interchange::audit(&model).unwrap();
    assert_eq!(
        report
            .field("through_bores")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        2
    );
    assert_eq!(
        report
            .field("shared_faces")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        74
    );
}

#[test]
fn moved_hole_in_timeline_changes_geometry_and_keeps_the_history_admissible() {
    let mut variant = original();
    let hole = feature_mut(&mut variant, 3);
    set(
        hole,
        "position",
        Json::Array(vec![
            Json::Number("44".into()),
            Json::Number("8".into()),
            Json::Number("20".into()),
        ]),
    );
    let original_model = import_sketch_timeline(&original()).unwrap();
    let moved_model = import_sketch_timeline(&variant).unwrap();
    assert_ne!(original_model, moved_model);
    assert_eq!(
        moved_model.field("blocks").unwrap().items().unwrap().len(),
        48
    );
    assert_eq!(
        interchange::audit(&moved_model)
            .unwrap()
            .field("status")
            .unwrap()
            .text()
            .unwrap(),
        "valid"
    );
}

#[test]
fn altered_or_unresolved_host_references_are_refused() {
    let mut bad = original();
    set(
        feature_mut(&mut bad, 2),
        "edge_names",
        Json::Array(vec![Json::String("F2:side:l1|F2:side:l2".into())]),
    );
    assert!(import_sketch_timeline(&bad).is_err());

    let mut bad = original();
    set(
        feature_mut(&mut bad, 3),
        "direction",
        Json::Array(vec![
            Json::Number("1".into()),
            Json::Number("0".into()),
            Json::Number("0".into()),
        ]),
    );
    assert!(import_sketch_timeline(&bad).is_err());

    let mut bad = original();
    set(
        feature_mut(&mut bad, 4),
        "diameter",
        Json::String("9.25".into()),
    );
    assert!(import_sketch_timeline(&bad).is_err());

    let mut bad = original();
    set(feature_mut(&mut bad, 1), "suppressed", Json::Bool(true));
    assert!(import_sketch_timeline(&bad).is_err());
}
