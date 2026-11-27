use solidcraft_trivariate::{interchange, json::parse, mechanical_job::solve_job, thb_job};

#[test]
fn declared_traction_follows_the_original_end_face_after_local_refinement() {
    let recipe = parse(include_str!("references/nurbs_strip.recipe.json")).unwrap();
    let plan = parse(include_str!("references/thb_interior.plan.json")).unwrap();
    let job = parse(include_str!("references/thb_end_traction.job.json")).unwrap();
    let model = interchange::build(&recipe).unwrap();
    let package = thb_job::compile(&model, &plan).unwrap();
    let result = solve_job(&package, &job).unwrap();
    let number = |key| result.field(key).unwrap().number::<f64>().unwrap();
    assert_eq!(number("scalar_dofs"), 544.);
    assert_eq!(number("free_vector_dofs"), 1401.);
    assert_eq!(
        result
            .field("coefficients_xyz_mm")
            .unwrap()
            .items()
            .unwrap()
            .len(),
        544
    );
    assert!((number("energy_n_mm") / (number("volume_mm3") / 20.) - 1.).abs() < 1e-9);
    assert!(number("relative_residual") < 1e-9);
    for (a, b) in result
        .field("force_n")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .zip(result.field("reaction_n").unwrap().items().unwrap())
    {
        assert!((a.number::<f64>().unwrap() + b.number::<f64>().unwrap()).abs() < 1e-8);
    }
}
