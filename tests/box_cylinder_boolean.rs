use solidcraft_trivariate::{
    interchange,
    json::{self, Json},
};
use std::collections::BTreeMap;

fn recipe() -> Json {
    json::parse(include_str!("references/box_minus_cylinder.recipe.json")).unwrap()
}

fn object_mut(value: &mut Json) -> &mut BTreeMap<String, Json> {
    match value {
        Json::Object(fields) => fields,
        _ => panic!("fixture object expected"),
    }
}

fn array_numbers(values: &[i32]) -> Json {
    Json::Array(values.iter().map(|v| Json::Number(v.to_string())).collect())
}

#[test]
fn centered_box_minus_cylinder_builds_and_replays() {
    let model = interchange::build(&recipe()).unwrap();
    assert_eq!(
        model.field("schema").unwrap().text().unwrap(),
        "box_cylinder_boolean_v1"
    );
    assert_eq!(
        model.field("operation").unwrap().text().unwrap(),
        "difference"
    );
    assert_eq!(model.field("blocks").unwrap().items().unwrap().len(), 4);
    assert_eq!(
        model
            .field("partition_plan")
            .unwrap()
            .field("method")
            .unwrap()
            .text()
            .unwrap(),
        "four_sector_local_sweeps"
    );
    assert_eq!(
        model
            .field("volume_mm3")
            .unwrap()
            .field("constant")
            .unwrap()
            .number::<i64>()
            .unwrap(),
        19_200
    );
    assert_eq!(
        model
            .field("volume_mm3")
            .unwrap()
            .field("pi_coefficient")
            .unwrap()
            .number::<i64>()
            .unwrap(),
        -300
    );
    let audit = interchange::audit(&model).unwrap();
    assert_eq!(audit.field("status").unwrap().text().unwrap(), "valid");
}

#[test]
fn rejects_non_centered_or_non_through_operands() {
    let mut input = recipe();
    let root = object_mut(&mut input);
    object_mut(root.get_mut("tool").unwrap()).insert("center_mm".into(), array_numbers(&[19, 20]));
    assert!(interchange::build(&input).is_err());

    let mut input = recipe();
    let root = object_mut(&mut input);
    object_mut(root.get_mut("tool").unwrap()).insert("interval_mm".into(), array_numbers(&[0, 14]));
    assert!(interchange::build(&input).is_err());
}

#[test]
fn exact_world_z_translation_and_model_tampering_are_checked() {
    let mut shifted = recipe();
    let root = object_mut(&mut shifted);
    let stock = object_mut(root.get_mut("stock").unwrap());
    stock.insert("min_mm".into(), array_numbers(&[0, 0, 10]));
    stock.insert("max_mm".into(), array_numbers(&[40, 40, 22]));
    object_mut(root.get_mut("tool").unwrap()).insert("interval_mm".into(), array_numbers(&[8, 24]));
    let base_model = interchange::build(&recipe()).unwrap();
    let shifted_model = interchange::build(&shifted).unwrap();
    let base_blocks = base_model.field("blocks").unwrap().items().unwrap();
    let shifted_blocks = shifted_model.field("blocks").unwrap().items().unwrap();
    assert_eq!(base_blocks.len(), shifted_blocks.len());
    for (base, moved) in base_blocks.iter().zip(shifted_blocks) {
        let base_controls = base.field("controls").unwrap().items().unwrap();
        let moved_controls = moved.field("controls").unwrap().items().unwrap();
        for (base_point, moved_point) in base_controls.iter().zip(moved_controls) {
            let p = base_point.items().unwrap();
            let q = moved_point.items().unwrap();
            let expected_z = p[2].number::<i64>().unwrap() + 10 * p[3].number::<i64>().unwrap();
            assert_eq!(q[2].number::<i64>().unwrap(), expected_z);
        }
    }

    let mut tampered = base_model.clone();
    object_mut(object_mut(&mut tampered).get_mut("tool").unwrap())
        .insert("radius_mm".into(), Json::Number("6".into()));
    assert!(interchange::audit(&tampered).is_err());
}
