use solidcraft_trivariate::{interchange, json::parse};

#[test]
fn slender_plate_uses_four_sector_host_and_two_box_infills() {
    let recipe=parse(r#"{"family":"orthogonal_sketch_features_v1","outer":{"id":"plate","vertices":[[0,0],[100,0],[100,12],[0,12]]},"height_mm":12,"fillets":[],"bores":[{"id":"hole","axis":"z","center_mm":[80,6,6],"radius_mm":3}]}"#).unwrap();
    let model = interchange::build(&recipe).unwrap();
    let report = interchange::audit(&model).unwrap();
    assert_eq!(
        report.field("blocks").unwrap().number::<usize>().unwrap(),
        6
    );
    assert_eq!(
        report
            .field("through_bores")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        1
    );
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
    assert!(model
        .field("partition_plan")
        .unwrap()
        .field("fallback_reason")
        .unwrap()
        .text()
        .unwrap()
        .contains("compact_local_host_selected"));
}

#[test]
fn host_refuses_a_bore_center_outside_its_extrusion_interval() {
    let recipe=parse(r#"{"family":"orthogonal_sketch_features_v1","outer":{"id":"plate","vertices":[[0,0],[100,0],[100,12],[0,12]]},"height_mm":12,"fillets":[],"bores":[{"id":"hole","axis":"z","center_mm":[80,6,12],"radius_mm":3}]}"#).unwrap();
    assert!(interchange::build(&recipe).is_err());
}

#[test]
fn two_bores_on_one_bracket_arm_use_disjoint_local_hosts_and_preserve_ten_block_case() {
    let single = r#"{"family":"orthogonal_sketch_features_v1","outer":{"id":"profile","vertices":[[0,0],[120,0],[120,24],[24,24],[24,120],[0,120]]},"height_mm":24,"fillets":[{"id":"inside_round","vertex":[24,24],"radius_mm":6}],"bores":[{"id":"a","axis":"y","center_mm":[70,12,12],"radius_mm":3},{"id":"b","axis":"y","center_mm":[95,12,12],"radius_mm":3}]}"#;
    let model = interchange::build(&parse(single).unwrap()).unwrap();
    let report = interchange::audit(&model).unwrap();
    assert_eq!(report.field("status").unwrap().text().unwrap(), "valid");
    assert_eq!(
        report.field("blocks").unwrap().number::<usize>().unwrap(),
        11
    );
    assert_eq!(
        report
            .field("through_bores")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        2
    );

    let legacy = r#"{"family":"orthogonal_sketch_features_v1","outer":{"id":"profile","vertices":[[0,0],[120,0],[120,24],[24,24],[24,120],[0,120]]},"height_mm":24,"fillets":[{"id":"inside_round","vertex":[24,24],"radius_mm":6}],"bores":[{"id":"north","axis":"x","center_mm":[12,70,12],"radius_mm":3},{"id":"east","axis":"y","center_mm":[70,12,12],"radius_mm":3}]}"#;
    let model = interchange::build(&parse(legacy).unwrap()).unwrap();
    assert_eq!(
        interchange::audit(&model)
            .unwrap()
            .field("blocks")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        10
    );
}

#[test]
fn overlapping_largest_hosts_use_mid_gap_strips_and_mixed_arm_groups_build() {
    let close = r#"{"family":"orthogonal_sketch_features_v1","outer":{"id":"profile","vertices":[[0,0],[120,0],[120,24],[24,24],[24,120],[0,120]]},"height_mm":24,"fillets":[{"id":"inside_round","vertex":[24,24],"radius_mm":6}],"bores":[{"id":"a","axis":"y","center_mm":[70,12,12],"radius_mm":2},{"id":"b","axis":"y","center_mm":[80,12,12],"radius_mm":2}]}"#;
    let model = interchange::build(&parse(close).unwrap()).unwrap();
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

    let mixed = r#"{"family":"orthogonal_sketch_features_v1","outer":{"id":"profile","vertices":[[0,0],[120,0],[120,24],[24,24],[24,120],[0,120]]},"height_mm":24,"fillets":[{"id":"inside_round","vertex":[24,24],"radius_mm":6}],"bores":[{"id":"north_a","axis":"x","center_mm":[12,70,12],"radius_mm":3},{"id":"north_b","axis":"x","center_mm":[12,95,12],"radius_mm":3},{"id":"east","axis":"y","center_mm":[70,12,12],"radius_mm":3}]}"#;
    let model = interchange::build(&parse(mixed).unwrap()).unwrap();
    let report = interchange::audit(&model).unwrap();
    assert_eq!(report.field("status").unwrap().text().unwrap(), "valid");
    assert_eq!(
        report
            .field("through_bores")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        3
    );
}

/// Regression for the 64-block grid fallback, which rejected this compact
/// 2+1-bore candidate on an unmatched wing/elbow face.
#[test]
fn fieldcraft_three_bore_arm_is_partitioned_with_full_section_strips() {
    let recipe = parse(include_str!(
        "regressions/fieldcraft_three_bores.recipe.json"
    ))
    .unwrap();
    let model = interchange::build(&recipe).unwrap();
    let report = interchange::audit(&model).unwrap();
    assert_eq!(report.field("status").unwrap().text().unwrap(), "valid");
    assert_eq!(
        report
            .field("through_bores")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        3
    );
    assert_eq!(
        report.field("blocks").unwrap().number::<usize>().unwrap(),
        14
    );
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
}
