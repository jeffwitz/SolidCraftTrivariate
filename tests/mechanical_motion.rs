use solidcraft_trivariate::{
    interchange::build,
    json::{parse, Json},
    mechanical_job::solve_job,
};

fn fixture() -> (Json, Json) {
    let case = parse(include_str!(
        "references/affine_mechanical_jobs/rigid-rotation.json"
    ))
    .unwrap();
    (
        build(case.field("recipe").unwrap()).unwrap(),
        case.field("job").unwrap().clone(),
    )
}

fn set(target: &mut Json, key: &str, value: Json) {
    let Json::Object(fields) = target else {
        panic!("expected fixture object")
    };
    fields.insert(key.into(), value);
}

#[test]
fn independently_selected_rigid_components_have_identical_physical_lifts() {
    let (model, mut job) = fixture();
    let bundled = solve_job(&model, &job).unwrap();
    let base = job.field("supports").unwrap().items().unwrap()[0].clone();
    let split = (0..3)
        .rev()
        .map(|c| {
            let mut part = base.clone();
            set(
                &mut part,
                "components",
                Json::Array(vec![Json::Number(c.to_string())]),
            );
            part
        })
        .collect();
    set(&mut job, "supports", Json::Array(split));
    let result = solve_job(&model, &job).unwrap();
    for key in [
        "component_lifts",
        "displacement_mm_xyz",
        "strain_energy_n_mm",
        "support_work_n_mm",
        "support_reaction_n",
    ] {
        assert_eq!(result.field(key), bundled.field(key), "{key}");
    }
    assert_eq!(
        result
            .field("strain_energy_n_mm")
            .unwrap()
            .number::<f64>()
            .unwrap(),
        0.
    );
}

#[test]
fn older_job_schema_refuses_new_boundary_data_and_malformed_affines() {
    let (model, mut job) = fixture();
    set(
        &mut job,
        "schema",
        Json::String("solidcraft_mechanical_job_v1".into()),
    );
    assert!(solve_job(&model, &job).is_err());
    let (_, mut job) = fixture();
    let mut support = job.field("supports").unwrap().items().unwrap()[0].clone();
    set(
        &mut support,
        "affine_mm",
        parse(r#"{"denominator":1,"numerators":[[0],[0],[0]]}"#).unwrap(),
    );
    set(&mut job, "supports", Json::Array(vec![support]));
    assert!(solve_job(&model, &job).is_err());
}
