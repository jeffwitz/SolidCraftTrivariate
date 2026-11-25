use solidcraft_trivariate::{json, solidcraft_bridge::import_script};

fn interval(model: &json::Json) -> [i32; 2] {
    let bore = &model.field("bores").unwrap().items().unwrap()[0];
    bore.field("z_interval_mm")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|value| value.number().unwrap())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

fn center(model: &json::Json) -> [i32; 2] {
    let bore = &model.field("bores").unwrap().items().unwrap()[0];
    bore.field("center_mm")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|value| value.number().unwrap())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

fn features(model: &json::Json) -> Vec<(String, String)> {
    model
        .field("design_faces")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|face| {
            (
                face.field("block").unwrap().text().unwrap().to_owned(),
                face.field("feature").unwrap().text().unwrap().to_owned(),
            )
        })
        .collect()
}

#[test]
fn host_generated_blind_bore_histories_map_to_correct_g21_intervals_and_caps() {
    for (history, expected_center, expected_interval, cap_feature) in [
        (
            include_str!("references/solidcraft-blind-bore/centered_top.history.json"),
            [0, 0],
            [8, 12],
            "bore/blind_bore/lower",
        ),
        (
            include_str!("references/solidcraft-blind-bore/decentered_top.history.json"),
            [4, -3],
            [8, 12],
            "bore/blind_bore/lower",
        ),
        (
            include_str!("references/solidcraft-blind-bore/opposite_bottom.history.json"),
            [-5, 6],
            [0, 5],
            "bore/blind_bore/upper",
        ),
    ] {
        let script = json::parse(history).unwrap();
        let model = import_script(&script).unwrap();
        assert_eq!(center(&model), expected_center);
        assert_eq!(interval(&model), expected_interval);
        let selected = features(&model);
        assert!(selected.iter().any(|(_, feature)| feature == cap_feature));
        assert!(selected
            .iter()
            .any(|(_, feature)| feature == "bore/blind_bore/wall"));
        let wall_blocks = selected
            .iter()
            .filter(|(_, feature)| feature == "bore/blind_bore/wall")
            .map(|(block, _)| block.as_str())
            .collect::<Vec<_>>();
        assert!(!wall_blocks.is_empty());
        assert!(wall_blocks.iter().all(|id| id.ends_with(&format!(
            "_z{}_{}",
            expected_interval[0], expected_interval[1]
        ))));
    }
}

#[test]
fn direct_history_accepts_zero_decimal_spelling_from_admitted_host_recipe() {
    let direct = json::parse(include_str!(
        "references/solidcraft-blind-bore/centered_top.json"
    ))
    .unwrap();
    let from_direct = import_script(&direct).unwrap();
    let host_timeline = json::parse(include_str!(
        "references/solidcraft-blind-bore/centered_top.history.json"
    ))
    .unwrap();
    let from_timeline = import_script(&host_timeline).unwrap();
    assert_eq!(from_direct, from_timeline);
}

#[test]
fn unsupported_blind_hole_options_are_refused() {
    for (key, value) in [
        ("depth", json::Json::Number("12".into())),
        ("type", json::Json::String("drilled".into())),
        ("diameter", json::Json::Number("8.5".into())),
    ] {
        let mut script = json::parse(include_str!(
            "references/solidcraft-blind-bore/centered_top.json"
        ))
        .unwrap();
        let json::Json::Object(root) = &mut script else {
            unreachable!()
        };
        let json::Json::Array(commands) = root.get_mut("commands").unwrap() else {
            unreachable!()
        };
        let json::Json::Object(hole_command) = &mut commands[1] else {
            unreachable!()
        };
        let json::Json::Object(params) = hole_command.get_mut("params").unwrap() else {
            unreachable!()
        };
        params.insert(key.into(), value);
        assert!(import_script(&script).is_err(), "must refuse altered {key}");
    }

    let mut history = json::parse(include_str!(
        "references/solidcraft-blind-bore/centered_top.history.json"
    ))
    .unwrap();
    let json::Json::Object(root) = &mut history else {
        unreachable!()
    };
    let json::Json::Array(commands) = root.get_mut("commands").unwrap() else {
        unreachable!()
    };
    let json::Json::Object(hole_command) = &mut commands[1] else {
        unreachable!()
    };
    let json::Json::Object(params) = hole_command.get_mut("params").unwrap() else {
        unreachable!()
    };
    let json::Json::Object(feature) = params.get_mut("feature").unwrap() else {
        unreachable!()
    };
    feature.insert(
        "direction".into(),
        json::Json::Array(vec![
            json::Json::Number("0.0".into()),
            json::Json::Number("0.0".into()),
            json::Json::Number("1.0".into()),
        ]),
    );
    assert!(import_script(&history).is_err());
}
