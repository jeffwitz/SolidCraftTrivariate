use solidcraft_trivariate::{
    interchange,
    json::{object, parse, Json},
    solidcraft_box_timeline::import_box_timeline,
    solidcraft_boxes::import_boxes,
    solidcraft_bridge::import_script,
};
use std::{fs, path::Path};

fn n(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}

fn interval_box(value: &Json) -> [[i32; 2]; 3] {
    let axes = value.items().unwrap();
    std::array::from_fn(|axis| {
        let pair = axes[axis].items().unwrap();
        [pair[0].number().unwrap(), pair[1].number().unwrap()]
    })
}

fn half(twice: i32) -> Json {
    let sign = if twice < 0 { "-" } else { "" };
    let magnitude = twice.unsigned_abs();
    Json::Number(if magnitude.is_multiple_of(2) {
        format!("{sign}{}", magnitude / 2)
    } else {
        format!("{sign}{}.5", magnitude / 2)
    })
}

fn timeline_feature(
    index: usize,
    name: &str,
    operation: &str,
    bounds: [[i32; 2]; 3],
    centered: bool,
) -> Json {
    let length = [0, 1, 2].map(|axis| bounds[axis][1] - bounds[axis][0]);
    let point = if centered {
        std::array::from_fn(|axis| half(bounds[axis][0] + bounds[axis][1]))
    } else {
        bounds.map(|[lo, _]| n(lo))
    };
    let mut fields = vec![
        (
            "corner".to_owned(),
            Json::Array(point.into_iter().collect()),
        ),
        ("height".to_owned(), Json::String(length[2].to_string())),
        ("id".to_owned(), n(index)),
        ("length".to_owned(), Json::String(length[0].to_string())),
        ("name".to_owned(), Json::String(name.into())),
        ("operation".to_owned(), Json::String(operation.into())),
        (
            "param_names".to_owned(),
            Json::Array(
                (0..3)
                    .map(|axis| Json::String(format!("d{}", (index - 1) * 3 + axis + 1)))
                    .collect(),
            ),
        ),
        ("suppressed".to_owned(), Json::Bool(false)),
        ("type".to_owned(), Json::String("box".into())),
        ("width".to_owned(), Json::String(length[1].to_string())),
    ];
    if index == 1 {
        fields.push((
            "body_names".to_owned(),
            Json::Array(vec![Json::String("Body1".into())]),
        ));
    }
    if centered {
        fields.push(("centred".to_owned(), Json::Bool(true)));
    }
    object([
        ("command", Json::String("timeline.add".into())),
        (
            "params",
            object([("feature", Json::Object(fields.into_iter().collect()))]),
        ),
    ])
}

fn history(recipe: &Json) -> Json {
    let mut commands = vec![timeline_feature(
        1,
        "base",
        "new_body",
        interval_box(recipe.field("base_mm").unwrap()),
        false,
    )];
    for (index, feature) in recipe
        .field("features")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .enumerate()
    {
        commands.push(timeline_feature(
            index + 2,
            feature.field("id").unwrap().text().unwrap(),
            if feature.field("operation").unwrap().text().unwrap() == "add" {
                "join"
            } else {
                "cut"
            },
            interval_box(feature.field("bounds_mm").unwrap()),
            index % 2 == 0,
        ));
    }
    object([
        (
            "about",
            object([
                ("design", Json::String("editable".into())),
                ("format", Json::String("solidcraft/66".into())),
                ("replay", Json::String("canonical timeline".into())),
                ("version", Json::String("0.48.1".into())),
            ]),
        ),
        ("commands", Json::Array(commands)),
    ])
}

#[test]
fn canonical_timelines_match_all_frozen_native_programs() {
    let corpus = parse(include_str!("references/orthogonal_corpus.json")).unwrap();
    let mut accepted = 0;
    let mut rejected = 0;
    for case in corpus.items().unwrap() {
        let recipe = case.field("recipe").unwrap();
        let timeline = history(recipe);
        if case.field("admitted").unwrap() == &Json::Bool(true) {
            let model = import_box_timeline(&timeline).unwrap_or_else(|error| {
                panic!("{} should import, got {error}", case.field("id").unwrap())
            });
            assert_eq!(model, interchange::build(recipe).unwrap());
            accepted += 1;
        } else {
            assert!(
                import_box_timeline(&timeline).is_err(),
                "{}",
                case.field("id").unwrap()
            );
            rejected += 1;
        }
    }
    assert_eq!((accepted, rejected), (11, 6));
}

#[test]
fn resolved_center_and_corner_placements_preserve_the_same_program() {
    let corpus = parse(include_str!("references/orthogonal_corpus.json")).unwrap();
    let recipe = corpus
        .items()
        .unwrap()
        .iter()
        .find(|case| case.field("id").unwrap().text().unwrap() == "crossing_tunnels")
        .unwrap()
        .field("recipe")
        .unwrap();
    let model = import_box_timeline(&history(recipe)).unwrap();
    assert_eq!(model, interchange::build(recipe).unwrap());
}

#[test]
fn rejects_suppressed_deleted_unresolved_or_reordered_features() {
    let corpus = parse(include_str!("references/orthogonal_corpus.json")).unwrap();
    let recipe = corpus.items().unwrap()[3].field("recipe").unwrap();
    let base = history(recipe);

    let mut suppressed = base.clone();
    let Json::Object(root) = &mut suppressed else {
        unreachable!()
    };
    let Json::Array(commands) = root.get_mut("commands").unwrap() else {
        unreachable!()
    };
    let Json::Object(command) = &mut commands[1] else {
        unreachable!()
    };
    let Json::Object(params) = command.get_mut("params").unwrap() else {
        unreachable!()
    };
    let Json::Object(feature) = params.get_mut("feature").unwrap() else {
        unreachable!()
    };
    feature.insert("suppressed".into(), Json::Bool(true));
    assert!(import_box_timeline(&suppressed).is_err());

    let mut deleted = base.clone();
    let Json::Object(root) = &mut deleted else {
        unreachable!()
    };
    let Json::Array(commands) = root.get_mut("commands").unwrap() else {
        unreachable!()
    };
    commands.remove(1);
    assert!(import_box_timeline(&deleted).is_err());

    let mut expression = base.clone();
    let Json::Object(root) = &mut expression else {
        unreachable!()
    };
    let Json::Array(commands) = root.get_mut("commands").unwrap() else {
        unreachable!()
    };
    let Json::Object(command) = &mut commands[1] else {
        unreachable!()
    };
    let Json::Object(params) = command.get_mut("params").unwrap() else {
        unreachable!()
    };
    let Json::Object(feature) = params.get_mut("feature").unwrap() else {
        unreachable!()
    };
    feature.insert("length".into(), Json::String("length_expr + 1".into()));
    assert!(import_box_timeline(&expression).is_err());

    let mut reordered = base;
    let Json::Object(root) = &mut reordered else {
        unreachable!()
    };
    let Json::Array(commands) = root.get_mut("commands").unwrap() else {
        unreachable!()
    };
    let Json::Object(command) = &mut commands[1] else {
        unreachable!()
    };
    let Json::Object(params) = command.get_mut("params").unwrap() else {
        unreachable!()
    };
    let Json::Object(feature) = params.get_mut("feature").unwrap() else {
        unreachable!()
    };
    feature.insert(
        "param_names".into(),
        Json::Array(vec![
            Json::String("d5".into()),
            Json::String("d4".into()),
            Json::String("d6".into()),
        ]),
    );
    assert!(import_box_timeline(&reordered).is_err());
}

#[test]
fn production_archives_match_direct_scripts_when_host_artifacts_are_present() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let cases = root.join("build/solidcraft-boxes");
    if !cases.exists() {
        return;
    }
    let mut checked = 0;
    for entry in fs::read_dir(&cases).unwrap() {
        let dir = entry.unwrap().path();
        let archived = dir.join("canonical.json");
        let direct = dir.join("script.json");
        if !archived.exists() || !direct.exists() {
            continue;
        }
        let archived = parse(&fs::read_to_string(archived).unwrap()).unwrap();
        let direct = parse(&fs::read_to_string(direct).unwrap()).unwrap();
        let timeline_model = import_box_timeline(&archived)
            .unwrap_or_else(|error| panic!("{} archive rejected: {error}", dir.display()));
        assert_eq!(import_script(&archived).unwrap(), timeline_model);
        let direct_model = import_boxes(&direct)
            .unwrap_or_else(|error| panic!("{} direct script rejected: {error}", dir.display()));
        assert_eq!(timeline_model, direct_model, "{}", dir.display());
        checked += 1;
    }
    if checked != 0 {
        assert_eq!(checked, 11);
    }
}
