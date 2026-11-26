use solidcraft_trivariate::{interchange, json::parse};

fn recipe(bores: &str) -> solidcraft_trivariate::json::Json {
    parse(&format!(r#"{{"family":"orthogonal_sketch_features_v1","outer":{{"id":"profile","vertices":[[0,0],[120,0],[120,24],[24,24],[24,120],[0,120]]}},"height_mm":24,"fillets":[{{"id":"inside_round","vertex":[24,24],"radius_mm":6}}],"bores":[{bores}]}}"#)).unwrap()
}

fn admitted(bores: &str, expected: usize) {
    let model = interchange::build(&recipe(bores)).unwrap();
    let report = interchange::audit(&model).unwrap();
    assert_eq!(report.field("status").unwrap().text().unwrap(), "valid");
    assert_eq!(
        report
            .field("through_bores")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        2
    );
    assert_eq!(
        report.field("blocks").unwrap().number::<usize>().unwrap(),
        expected
    );
}

#[test]
fn stacked_z_bores_use_exact_host_strips_and_closed_partition() {
    admitted(
        r#"{"id":"z_low","axis":"z","center_mm":[70,12,6],"radius_mm":3},{"id":"z_high","axis":"z","center_mm":[95,12,18],"radius_mm":3}"#,
        11,
    );
}

#[test]
fn diagonally_displaced_parallel_bores_use_exact_midgap_partition() {
    admitted(
        r#"{"id":"diag_a","axis":"y","center_mm":[70,12,6],"radius_mm":2},{"id":"diag_b","axis":"y","center_mm":[95,12,18],"radius_mm":2}"#,
        11,
    );
}

#[test]
fn fieldcraft_stacks_close_the_cut_through_both_elbow_maps_and_opposite_arm() {
    for (path, opposite) in [
        (
            "references/fieldcraft_east_stack.recipe.json",
            "propagated_wing_1",
        ),
        (
            "references/fieldcraft_north_stack.recipe.json",
            "propagated_wing_0",
        ),
    ] {
        let source = match path {
            "references/fieldcraft_east_stack.recipe.json" => {
                include_str!("references/fieldcraft_east_stack.recipe.json")
            }
            _ => include_str!("references/fieldcraft_north_stack.recipe.json"),
        };
        let model = interchange::build(&parse(source).unwrap()).unwrap();
        let audit = interchange::audit(&model).unwrap();
        assert_eq!(audit.field("status").unwrap().text().unwrap(), "valid");
        assert_eq!(
            audit.field("blocks").unwrap().number::<usize>().unwrap(),
            17
        );
        assert_eq!(
            audit
                .field("through_bores")
                .unwrap()
                .number::<usize>()
                .unwrap(),
            3
        );
        let plan = model.field("partition_plan").unwrap();
        assert_eq!(
            plan.field("method").unwrap().text().unwrap(),
            "propagated_isoparametric_cuts"
        );
        let rounds = plan.field("cut_rounds").unwrap().items().unwrap();
        assert_eq!(rounds.len(), 1);
        assert_eq!(
            rounds[0]
                .field("requests")
                .unwrap()
                .number::<usize>()
                .unwrap(),
            2
        );
        assert_eq!(
            rounds[0]
                .field("closed_cuts")
                .unwrap()
                .number::<usize>()
                .unwrap(),
            3
        );
        assert_eq!(
            rounds[0]
                .field("blocks")
                .unwrap()
                .number::<usize>()
                .unwrap(),
            17
        );
        let sources = plan.field("cut_sources").unwrap().items().unwrap();
        assert_eq!(sources.len(), 3);
        assert!(sources
            .iter()
            .any(|s| s.field("block").unwrap().text().unwrap() == "elbow_east_fillet"));
        assert!(sources
            .iter()
            .any(|s| s.field("block").unwrap().text().unwrap() == "elbow_north_fillet"));
        assert!(sources.iter().any(|s| s
            .field("block")
            .unwrap()
            .text()
            .unwrap()
            .starts_with(opposite)));
    }
}

#[test]
fn fieldcraft_diagonal_three_bore_case_uses_disjoint_hosts_without_axial_cut() {
    let model = interchange::build(
        &parse(include_str!(
            "references/fieldcraft_diagonal_bores.recipe.json"
        ))
        .unwrap(),
    )
    .unwrap();
    let audit = interchange::audit(&model).unwrap();
    assert_eq!(audit.field("status").unwrap().text().unwrap(), "valid");
    assert_eq!(
        audit.field("blocks").unwrap().number::<usize>().unwrap(),
        14
    );
    assert_eq!(
        audit
            .field("through_bores")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        3
    );
    let plan = model.field("partition_plan").unwrap();
    assert_eq!(
        plan.field("method").unwrap().text().unwrap(),
        "local_bore_hosts"
    );
    assert!(plan.field("cut_rounds").is_err());
}
