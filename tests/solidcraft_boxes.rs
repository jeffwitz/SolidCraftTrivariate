use solidcraft_trivariate::{
    interchange,
    json::{object, parse, Json},
    solidcraft_boxes::import_boxes,
};

fn number(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}

fn vector(values: impl IntoIterator<Item = i32>) -> Json {
    Json::Array(values.into_iter().map(number).collect())
}

fn bounds(value: &Json) -> [[i32; 2]; 3] {
    let axes = value.items().unwrap();
    std::array::from_fn(|axis| {
        let interval = axes[axis].items().unwrap();
        [interval[0].number().unwrap(), interval[1].number().unwrap()]
    })
}

fn decimal_half(twice: i32) -> Json {
    let magnitude = twice.unsigned_abs();
    let sign = if twice < 0 { "-" } else { "" };
    let token = if magnitude.is_multiple_of(2) {
        format!("{sign}{}", magnitude / 2)
    } else {
        format!("{sign}{}.5", magnitude / 2)
    };
    Json::Number(token)
}

fn command(name: &str, operation: &str, bounds_mm: [[i32; 2]; 3], center: bool) -> Json {
    let dimensions = bounds_mm.map(|[lo, hi]| hi - lo);
    let placement = if center {
        std::array::from_fn(|axis| decimal_half(bounds_mm[axis][0] + bounds_mm[axis][1]))
    } else {
        bounds_mm.map(|[lo, _]| number(lo))
    };
    let mut fields = vec![
        ("length", number(dimensions[0])),
        ("width", number(dimensions[1])),
        ("height", number(dimensions[2])),
        ("name", Json::String(name.into())),
        ("operation", Json::String(operation.into())),
    ];
    if center {
        fields.push(("center", Json::Array(placement.into())));
    } else {
        fields.push(("corner", Json::Array(placement.into())));
    }
    object([
        ("command", Json::String("solid.box".into())),
        (
            "params",
            Json::Object(
                fields
                    .into_iter()
                    .map(|(key, value)| (key.to_owned(), value))
                    .collect(),
            ),
        ),
    ])
}

fn script_for(recipe: &Json) -> Json {
    let base = bounds(recipe.field("base_mm").unwrap());
    let mut commands = vec![command("base", "new", base, false)];
    for (index, feature) in recipe
        .field("features")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .enumerate()
    {
        let operation = match feature.field("operation").unwrap().text().unwrap() {
            "add" => "join",
            "cut" => "cut",
            _ => unreachable!(),
        };
        commands.push(command(
            feature.field("id").unwrap().text().unwrap(),
            operation,
            bounds(feature.field("bounds_mm").unwrap()),
            index % 2 == 0,
        ));
    }
    object([("commands", Json::Array(commands))])
}

fn mutate_command(
    script: &mut Json,
    index: usize,
    change: impl FnOnce(&mut std::collections::BTreeMap<String, Json>),
) {
    let Json::Object(root) = script else {
        unreachable!()
    };
    let Json::Array(commands) = root.get_mut("commands").unwrap() else {
        unreachable!()
    };
    let Json::Object(command) = &mut commands[index] else {
        unreachable!()
    };
    let Json::Object(params) = command.get_mut("params").unwrap() else {
        unreachable!()
    };
    change(params);
}

#[test]
fn all_frozen_g16_command_histories_match_native_programs() {
    let corpus = parse(include_str!("references/orthogonal_corpus.json")).unwrap();
    let mut admitted = 0;
    let mut refused = 0;
    for case in corpus.items().unwrap() {
        let recipe = case.field("recipe").unwrap();
        let script = script_for(recipe);
        if case.field("admitted").unwrap() == &Json::Bool(true) {
            let imported = import_boxes(&script).unwrap_or_else(|error| {
                panic!(
                    "{} was expected to qualify, got {error}",
                    case.field("id").unwrap()
                )
            });
            let direct = interchange::build(recipe).unwrap();
            assert_eq!(imported, direct, "{}", case.field("id").unwrap());
            admitted += 1;
        } else {
            assert!(
                import_boxes(&script).is_err(),
                "{}",
                case.field("id").unwrap()
            );
            refused += 1;
        }
    }
    assert_eq!((admitted, refused), (11, 6));
}

#[test]
fn half_integer_centers_are_exact_when_their_bounds_are_integral() {
    let recipe = parse(include_str!("references/orthogonal_corpus.json")).unwrap();
    let case = recipe
        .items()
        .unwrap()
        .iter()
        .find(|case| case.field("id").unwrap().text().unwrap() == "crossing_tunnels")
        .unwrap();
    let script = script_for(case.field("recipe").unwrap());
    let Json::Object(root) = &script else {
        unreachable!()
    };
    let Json::Array(commands) = root.get("commands").unwrap() else {
        unreachable!()
    };
    let Json::Object(second) = &commands[1] else {
        unreachable!()
    };
    let Json::Object(params) = second.get("params").unwrap() else {
        unreachable!()
    };
    assert_eq!(
        params.get("center").unwrap(),
        &Json::Array(vec![
            number(6),
            Json::Number("5.5".into()),
            Json::Number("5.5".into())
        ])
    );
    assert!(import_boxes(&script).is_ok());
}

#[test]
fn center_object_form_is_equivalent_to_the_three_element_array() {
    let corpus = parse(include_str!("references/orthogonal_corpus.json")).unwrap();
    let case = corpus
        .items()
        .unwrap()
        .iter()
        .find(|case| case.field("id").unwrap().text().unwrap() == "crossing_tunnels")
        .unwrap();
    let mut script = script_for(case.field("recipe").unwrap());
    let Json::Object(root) = &mut script else {
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
    let Json::Array(center) = params.remove("center").unwrap() else {
        unreachable!()
    };
    params.insert(
        "center".into(),
        Json::Object(
            [
                ("x".into(), center[0].clone()),
                ("y".into(), center[1].clone()),
                ("z".into(), center[2].clone()),
            ]
            .into_iter()
            .collect(),
        ),
    );
    let array_form = script_for(case.field("recipe").unwrap());
    assert_eq!(
        import_boxes(&script).unwrap(),
        import_boxes(&array_form).unwrap()
    );
}

#[test]
fn expressions_options_and_invalid_body_transitions_are_rejected() {
    let corpus = parse(include_str!("references/orthogonal_corpus.json")).unwrap();
    let base_case = &corpus.items().unwrap()[0];
    let base = script_for(base_case.field("recipe").unwrap());

    let mut expression = base.clone();
    mutate_command(&mut expression, 0, |params| {
        params.insert("length".into(), Json::String("width_var + 1".into()));
    });
    assert!(import_boxes(&expression).is_err());

    let mut extra = base.clone();
    mutate_command(&mut extra, 0, |params| {
        params.insert(
            "targets".into(),
            Json::Array(vec![Json::String("Body1".into())]),
        );
    });
    assert!(import_boxes(&extra).is_err());

    let mut zero = base.clone();
    mutate_command(&mut zero, 0, |params| {
        params.insert("height".into(), number(0));
    });
    assert!(import_boxes(&zero).is_err());

    let mut dual_placement = base.clone();
    mutate_command(&mut dual_placement, 0, |params| {
        params.insert("center".into(), vector([6, 6, 6]));
    });
    assert!(import_boxes(&dual_placement).is_err());

    let mut next_body = base;
    let Json::Object(root) = &mut next_body else {
        unreachable!()
    };
    let Json::Array(commands) = root.get_mut("commands").unwrap() else {
        unreachable!()
    };
    commands.push(command("extra", "new", [[12, 14], [0, 2], [0, 2]], false));
    assert!(import_boxes(&next_body).is_err());
}
