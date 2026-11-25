use solidcraft_trivariate::{
    interchange::{audit, build},
    json::{parse, Json},
    solver_package::export_analysis,
};
use std::{fs, process::Command};

#[test]
fn reordered_input_reloads_and_derivatives_follow_the_serialized_gauge() {
    let recipe =
        parse(r#"{"family":"parameterized_bracket_v1","radii_sixteenths_mm":[56,64]}"#).unwrap();
    let mut model = build(&recipe).unwrap();
    let expected = export_analysis(&model, 4).unwrap();
    assert_eq!(
        audit(expected.field("source_model").unwrap()).unwrap(),
        audit(&model).unwrap()
    );
    let Json::Object(fields) = &mut model else {
        panic!()
    };
    let Some(Json::Array(blocks)) = fields.get_mut("blocks") else {
        panic!()
    };
    blocks.reverse();
    assert_eq!(expected, export_analysis(&model, 4).unwrap());
    // A common positive scale on a non-bore block leaves both its geometry and
    // the serialized radius jets on the other blocks unchanged.
    let Json::Object(fields) = &mut model else {
        panic!()
    };
    let Some(Json::Array(blocks)) = fields.get_mut("blocks") else {
        panic!()
    };
    let block = blocks
        .iter_mut()
        .find(|b| b.field("id").unwrap().text().unwrap().starts_with("arm_"))
        .unwrap();
    let Json::Object(fields) = block else {
        panic!()
    };
    let Some(Json::Array(net)) = fields.get_mut("controls") else {
        panic!()
    };
    for p in net {
        let Json::Array(h) = p else { panic!() };
        for x in h {
            *x = Json::Number((2 * x.number::<i64>().unwrap()).to_string());
        }
    }
    let scaled = export_analysis(&model, 4).unwrap();
    assert_eq!(
        scaled.field("derivatives").unwrap(),
        expected.field("derivatives").unwrap()
    );
}

#[test]
fn rejected_cli_export_preserves_the_last_solver_package() {
    let root = std::env::temp_dir().join(format!("solidcraft-x1-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let input = root.join("input.json");
    let output = root.join("package.json");
    fs::write(&input, r#"{"schema":"unknown_v1"}"#).unwrap();
    fs::write(&output, b"previous package").unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
        .arg("export-analysis")
        .arg(&input)
        .arg("8")
        .arg(&output)
        .output()
        .unwrap();
    assert!(!status.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"previous package");
    fs::remove_dir_all(root).unwrap();
}
