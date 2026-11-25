use solidcraft_trivariate::{
    interchange,
    json::{object, Json},
    mechanical_job::{differentiate_job, solve_job},
};

fn n(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}

fn recipe(t: [i32; 2]) -> Json {
    object([
        ("family", Json::String("rational_sector_prism_v1".into())),
        (
            "half_angle_tangent",
            Json::Array(t.into_iter().map(n).collect()),
        ),
        ("radii_mm", Json::Array(vec![n(2), n(6)])),
        ("height_mm", n(8)),
    ])
}

fn job(order: usize) -> Json {
    let refine = |level: u8, index: [u32; 3]| {
        object([
            ("root", Json::String("whole".into())),
            ("level", n(level)),
            ("index", Json::Array(index.into_iter().map(n).collect())),
        ])
    };
    object([
        (
            "schema",
            Json::String("solidcraft_mechanical_job_v1".into()),
        ),
        (
            "material",
            object([("young_mpa", n(210000)), ("poisson", n("0.3"))]),
        ),
        ("formulation", Json::String("displacement".into())),
        ("order", n(order)),
        (
            "plan",
            object([
                ("schema", Json::String("adaptive_refinement_plan_v1".into())),
                ("max_depth", n(2)),
                ("max_cells", n(128)),
                (
                    "steps",
                    Json::Array(vec![
                        object([("refine", Json::Array(vec![refine(0, [0, 0, 0])]))]),
                        object([("refine", Json::Array(vec![refine(1, [0, 0, 0])]))]),
                    ]),
                ),
            ]),
        ),
        (
            "supports",
            Json::Array(vec![object([
                (
                    "select",
                    object([("feature", Json::String("angle_start".into()))]),
                ),
                ("components", Json::Array(vec![n(0), n(1), n(2)])),
            ])]),
        ),
        (
            "loads",
            Json::Array(vec![object([
                (
                    "select",
                    object([("feature", Json::String("inner".into()))]),
                ),
                ("pressure_mpa", n("0.3")),
            ])]),
        ),
    ])
}

fn vector(result: &Json, key: &str) -> Vec<f64> {
    result
        .field(key)
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|v| v.number().unwrap())
        .collect()
}

fn layout(result: &Json) -> Vec<(String, Json)> {
    result
        .field("analysis_blocks")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|block| {
            (
                block.field("id").unwrap().text().unwrap().to_owned(),
                block.field("dofs").unwrap().clone(),
            )
        })
        .collect()
}

#[test]
fn local_dyadic_refinement_preserves_the_angular_shape_tangent() {
    let recipe_base = recipe([1, 2]);
    let model = interchange::build(&recipe_base).unwrap();
    let direction = object([(
        "schema",
        Json::String("solidcraft_sector_angle_velocity_v2".into()),
    )]);
    let plus = interchange::build(&recipe([128, 255])).unwrap();
    let minus = interchange::build(&recipe([127, 255])).unwrap();
    let step = 1. / 510.;
    let mut relative_errors = Vec::new();
    for order in [3, 4] {
        let job = job(order);
        let derivative = differentiate_job(&model, &job, &direction).unwrap();
        let base_solution = solve_job(&model, &job).unwrap();
        let plus_solution = solve_job(&plus, &job).unwrap();
        let minus_solution = solve_job(&minus, &job).unwrap();
        assert_eq!(layout(&base_solution), layout(&plus_solution));
        assert_eq!(layout(&base_solution), layout(&minus_solution));
        let up = vector(&plus_solution, "displacement_mm_xyz");
        let um = vector(&minus_solution, "displacement_mm_xyz");
        let analytic = vector(&derivative, "physical_displacement_rate_mm_xyz");
        let finite = up
            .iter()
            .zip(&um)
            .map(|(a, b)| (a - b) / (2. * step))
            .collect::<Vec<_>>();
        let error = analytic
            .iter()
            .zip(&finite)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        let scale = analytic.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-9);
        relative_errors.push(error / scale);
    }
    // The pre-rejection Q2 witness differed by 8.52%; q3/q4 converge.
    assert!(
        relative_errors[0] < 1e-4,
        "q3 tangent mismatch: {relative_errors:?}"
    );
    assert!(
        relative_errors[1] < 1e-4,
        "q4 tangent mismatch: {relative_errors:?}"
    );
}
