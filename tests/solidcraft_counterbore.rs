use solidcraft_trivariate::{
    json::{parse, Json},
    solidcraft_bridge::import_script,
};

fn read(source: &str) -> Json {
    parse(source).expect("valid fixture JSON")
}

#[test]
fn saved_host_histories_and_direct_commands_reach_the_same_certified_g23_family() {
    for (direct_src, history_src, center, opening, shoulder_z) in [
        (
            include_str!("references/solidcraft-counterbore/centered_top.json"),
            include_str!("references/solidcraft-counterbore/centered_top.history.json"),
            [0, 0],
            "end",
            8,
        ),
        (
            include_str!("references/solidcraft-counterbore/decentered_top.json"),
            include_str!("references/solidcraft-counterbore/decentered_top.history.json"),
            [4, -3],
            "end",
            8,
        ),
        (
            include_str!("references/solidcraft-counterbore/opposite_bottom.json"),
            include_str!("references/solidcraft-counterbore/opposite_bottom.history.json"),
            [-5, 6],
            "start",
            4,
        ),
    ] {
        let direct = import_script(&read(direct_src)).unwrap();
        let timeline = import_script(&read(history_src)).unwrap();
        assert_eq!(direct, timeline);
        assert_eq!(direct.field("schema").unwrap(), &Json::String("polygon_counterbores_v1".into()));
        let bore = &direct.field("bores").unwrap().items().unwrap()[0];
        assert_eq!(bore.field("center_mm").unwrap(), &Json::Array(center.map(|x| Json::Number(x.to_string())).into()));
        assert_eq!(bore.field("opening").unwrap(), &Json::String(opening.into()));
        assert_eq!(bore.field("depth_mm").unwrap(), &Json::Number("4".into()));
        let audit = solidcraft_trivariate::interchange::audit(&direct).unwrap();
        assert_eq!(audit.field("status").unwrap(), &Json::String("valid".into()));
        assert_eq!(audit.field("through_bores").unwrap(), &Json::Number("1".into()));
        assert_eq!(shoulder_z, if opening == "end" { 8 } else { 4 });
        for feature in ["bore/hole_2/pilot_wall", "bore/hole_2/counterbore_wall", "bore/hole_2/shoulder"] {
            assert!(direct.field("design_faces").unwrap().items().unwrap().iter().any(|face| {
                face.field("feature").unwrap() == &Json::String(feature.into())
            }), "missing G23 selected boundary {feature}");
        }
    }
}

#[test]
fn host_accepted_odd_counterbore_diameter_is_rejected_by_the_integer_g23_contract() {
    let host_admitted = read(include_str!("references/solidcraft-counterbore/odd_counterbore.json"));
    assert!(import_script(&host_admitted).is_err());
}
