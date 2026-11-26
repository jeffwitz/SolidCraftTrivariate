use solidcraft_trivariate::{
    interchange::{audit, build, plan_sketch},
    json::{parse, Json},
};

#[test]
fn two_polygon_holes_need_no_construction_field() {
    let corpus = parse(include_str!(
        "../docs/results/sketch-planner-polygon-evaluation/cases.json"
    ))
    .unwrap();
    let mut recipe = corpus
        .field("cases")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .find(|case| case.field("id").unwrap().text().unwrap() == "seeded_convex_10")
        .unwrap()
        .field("recipe")
        .unwrap()
        .clone();
    if let Json::Object(fields) = &mut recipe {
        fields.remove("construction");
    }
    let plan = plan_sketch(&recipe).unwrap();
    assert_eq!(
        plan.field("method").unwrap().text().unwrap(),
        "separated_convex_quad_rings"
    );
    let model = build(&recipe).unwrap();
    assert_eq!(model.field("partition_plan").unwrap(), &plan);
    let report = audit(&model).unwrap();
    assert_eq!(
        report.field("volume_exact_mm3").unwrap(),
        &parse("[24954264,1]").unwrap()
    );
    // Both dodecagonal through holes require at least twelve ring cells each.
    assert_eq!(
        report.field("blocks").unwrap().number::<usize>().unwrap(),
        24
    );
}
