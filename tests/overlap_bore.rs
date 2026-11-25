use solidcraft_trivariate::{
    interchange::{audit, build},
    json::{parse, Json},
    overlap_bore::overlap_blocks,
    overlap_bore_audit::audit_overlap,
};
#[test]
fn canonical_union_is_certified_and_all_four_limit_changes_are_detected() {
    let blocks = overlap_blocks().unwrap();
    let cert = audit_overlap(&blocks).unwrap();
    assert!(cert.determinant_lower_bound_mm3 > 0.);
    assert!(cert.radial_error_upper_mm < 1e-12);
    for i in 0..4 {
        let mut altered = blocks.clone();
        match i {
            0 => {
                altered[3].controls[10][1] += 1;
                altered[3].controls[11][1] += 1;
            }
            1 => {
                altered[1].controls[6][0] += 1;
                altered[1].controls[7][0] += 1;
            }
            2 => altered[6].controls[4][3] -= 1,
            _ => altered[4].id = altered[0].id.clone(),
        }
        assert!(audit_overlap(&altered).is_err());
    }
    let recipe = parse(r#"{"family":"overlap_bore_v1"}"#).unwrap();
    let mut model = build(&recipe).unwrap();
    assert_eq!(
        audit(&model)
            .unwrap()
            .field("status")
            .unwrap()
            .text()
            .unwrap(),
        "valid"
    );
    if let Json::Object(ref mut fields) = model {
        fields.insert("design_faces".into(), Json::Array(vec![]));
    }
    assert!(audit(&model).is_err());
    assert!(build(&parse(r#"{"family":"overlap_bore_v1","radius":3}"#).unwrap()).is_err());
}
