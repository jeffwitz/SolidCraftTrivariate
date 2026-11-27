use solidcraft_trivariate::json::{object, parse, Json};
use solidcraft_trivariate::solidcraft_preview::descriptor;

const TIMELINE: &str = include_str!("references/solidcraft_sketch_bracket.timeline.json");

#[test]
fn sketch_only_preview_uses_the_planners_separator_before_extrusion() {
    let full = parse(TIMELINE).unwrap();
    let command = full.field("commands").unwrap().items().unwrap()[0].clone();
    let script = object([("commands", Json::Array(vec![command.clone()]))]);
    let snapshot = object([
        (
            "schema",
            Json::String("solidcraft_trivariate_timeline_v1".into()),
        ),
        ("omitted", Json::Array(vec![])),
        ("script", script.clone()),
    ]);
    let output = descriptor(&snapshot, 1).unwrap();
    assert_eq!(output.field("source_script").unwrap(), &script);
    assert_eq!(output.field("block_count").unwrap().number::<u32>(), Ok(2));
    assert_eq!(
        output.field("internal_edges").unwrap(),
        &parse("[[[[0,1],[0,1]],[[8,1],[8,1]]]]").unwrap()
    );
    let feature = output.field("feature").unwrap();
    assert_eq!(feature.field("id").unwrap().number::<u32>(), Ok(2));
    let curves = feature
        .field("sketch")
        .unwrap()
        .field("curves")
        .unwrap()
        .items()
        .unwrap();
    assert_eq!(curves.len(), 1);
    assert_eq!(curves[0].field("construction").unwrap().bool(), Ok(true));
}

#[test]
fn refuse_missing_source_feature_and_unsupported_plane() {
    assert!(descriptor(&parse(TIMELINE).unwrap(), 2).is_err());
    assert!(descriptor(&parse(TIMELINE).unwrap(), 6).is_err());
    assert!(descriptor(&parse(&TIMELINE.replace("\"XY\"", "\"XZ\"")).unwrap(), 1).is_err());
}

#[test]
fn extrusion_with_or_without_its_reference_sketch_builds_identical_cells() {
    let full = parse(TIMELINE).unwrap();
    let commands = full.field("commands").unwrap().items().unwrap();
    let timeline = |commands| {
        object([
            ("commands", Json::Array(commands)),
            ("about", full.field("about").unwrap().clone()),
        ])
    };
    let original = timeline(commands[..2].to_vec());
    let overlay = descriptor(&original, 1)
        .unwrap()
        .field("feature")
        .unwrap()
        .clone();
    let step = object([
        ("command", Json::String("timeline.add".into())),
        ("params", object([("feature", overlay)])),
    ]);
    let with_overlay = timeline(vec![commands[0].clone(), commands[1].clone(), step]);
    let model = solidcraft_trivariate::solidcraft_polygon::extrude_preview(&original).unwrap();
    assert_eq!(
        model,
        solidcraft_trivariate::solidcraft_polygon::extrude_preview(&with_overlay).unwrap()
    );
    let audit = solidcraft_trivariate::interchange::audit(&model).unwrap();
    assert_eq!(
        audit.field("volume_exact_mm3").unwrap(),
        &parse("[32640,1]").unwrap()
    );
    assert_eq!(audit.field("blocks").unwrap().number::<usize>(), Ok(2));
}
