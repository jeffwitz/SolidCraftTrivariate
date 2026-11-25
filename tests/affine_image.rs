//! Exact affine transport, inverse audit and material shell identity.
use solidcraft_trivariate::{
    boundary_shell::boundary_shell,
    interchange::{audit, build},
    json::{object, parse, Json},
};

fn transformed(source: &str, matrix: &str) -> Json {
    object([
        ("family", Json::String("affine_image_v1".into())),
        ("source", parse(source).unwrap()),
        ("transform", parse(matrix).unwrap()),
    ])
}
const BOX: &str = r#"{"family":"orthogonal_feature_program_v1","base_mm":[[0,12],[0,10],[0,8]],"features":[{"id":"cavity","operation":"cut","bounds_mm":[[3,7],[2,9],[2,5]]}]}"#;
const MAP: &str =
    r#"{"matrix":[[0,-3,1],[4,1,2],[0,0,5]],"translation":[13,-7,4],"denominator":3}"#;

#[test]
fn material_shell_is_transportable_in_a_sheared_coordinate_system() {
    let model = build(&transformed(BOX, MAP)).unwrap();
    let certificate = audit(&model).unwrap();
    assert_eq!(
        certificate
            .field("determinant_numerator")
            .unwrap()
            .number::<u32>()
            .unwrap(),
        60
    );
    assert_eq!(
        certificate
            .field("hausdorff_error_upper_bound_mm")
            .unwrap()
            .number::<f64>()
            .unwrap(),
        0.0
    );
    let source = boundary_shell(model.field("source_model").unwrap()).unwrap();
    let image = boundary_shell(&model).unwrap();
    assert_eq!(source.components.len(), 2);
    assert_eq!(image.components, source.components);
}

fn change_path(j: &mut Json, path: &[&str]) {
    match (j, path) {
        (v, []) => *v = Json::Number((v.number::<i128>().unwrap_or(0) + 1).to_string()),
        (Json::Object(m), [k, rest @ ..]) => change_path(m.get_mut(*k).unwrap(), rest),
        (Json::Array(a), [i, rest @ ..]) => change_path(&mut a[i.parse::<usize>().unwrap()], rest),
        _ => panic!("unexpected mutation path"),
    }
}

#[test]
fn inverse_and_source_audits_reject_forged_data() {
    let model = build(&transformed(BOX, MAP)).unwrap();
    let paths: &[&[&str]] = &[
        &["source_model", "blocks", "0", "controls", "0", "1"],
        &["blocks", "0", "controls", "0", "0"],
        &["blocks", "0", "controls", "0", "3"],
        &["design_faces", "0", "feature"],
        &["topology"],
        &["transform", "translation", "2"],
    ];
    for path in paths {
        let mut forged = model.clone();
        change_path(&mut forged, path);
        assert!(audit(&forged).is_err(), "{path:?}");
    }
}

#[test]
fn refuse_negative_singular_out_of_range_and_recursive_images() {
    for matrix in [
        "[[1,0,0],[0,1,0],[0,0,-1]]",
        "[[1,2,3],[1,2,3],[0,0,1]]",
        "[[65,0,0],[0,1,0],[0,0,1]]",
    ] {
        let t = format!(r#"{{"matrix":{matrix},"translation":[0,0,0],"denominator":1}}"#);
        assert!(build(&transformed(BOX, &t)).is_err());
    }
    let inner = transformed(BOX, MAP);
    assert!(build(&transformed(&inner.to_string(), MAP)).is_err());
    assert!(build(&transformed(
        r#"{"family":"affine_pair_v1","widths_mm":[1,1]}"#,
        MAP
    ))
    .is_err());
}
