use solidcraft_trivariate::{
    interchange,
    json::{parse, Json},
    solidcraft_sketch_features::import_sketch_timeline,
};

fn original() -> Json {
    parse(include_str!("references/solidcraft_sketch_bracket.timeline.json")).unwrap()
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

#[test]
fn perpendicular_bores_share_a_guard_on_their_common_axial_layer() {
    // Radius pairs are in sixteenths of a millimeter: [x-axis bore, y-axis bore].
    // The 16x integer lattice preserves the source diameters exactly; the caller
    // can map the certified result back by an exact affine image.
    for [rx, ry] in [[80, 40], [74, 50], [74, 51], [64, 64], [40, 80]] {
        let recipe = interchange::build(
            &parse(&format!(
                r#"{{"family":"orthogonal_sketch_features_v1","outer":{{"id":"profile","vertices":[[0,0],[960,0],[960,128],[128,128],[128,800],[0,800]]}},"height_mm":640,"fillets":[{{"id":"inside_round","vertex":[128,128],"radius_mm":96}}],"bores":[{{"id":"north","axis":"x","center_mm":[64,512,320],"radius_mm":{rx}}},{{"id":"east","axis":"y","center_mm":[640,64,320],"radius_mm":{ry}}}]}}"#
            ))
            .unwrap(),
        )
        .unwrap();
        let report = interchange::audit(&recipe).unwrap();
        assert_eq!(report.field("status").unwrap().text().unwrap(), "valid");
        assert_eq!(report.field("through_bores").unwrap().number::<usize>().unwrap(), 2);
    }
}
