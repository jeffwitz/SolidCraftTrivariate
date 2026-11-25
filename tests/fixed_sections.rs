use solidcraft_trivariate::interchange::{audit, build, SECTION_FAMILIES};
use solidcraft_trivariate::json::{object, Json};
use solidcraft_trivariate::section_audit::audit_section;
use solidcraft_trivariate::sections::section_blocks;

#[test]
fn audit_transfers_certificates_only_for_the_complete_fixed_net() {
    for family in SECTION_FAMILIES {
        let blocks = section_blocks(family).unwrap();
        let certificate = audit_section(family, &blocks).unwrap();
        assert!(certificate.determinant_lower_bound_mm3 > 0.0);
        assert_eq!(
            certificate.design_faces.len(),
            certificate.incidence.exterior.len()
        );
        let mut reversed = blocks.clone();
        reversed.reverse();
        assert_eq!(
            certificate.incidence,
            audit_section(family, &reversed).unwrap().incidence
        );
        let mut altered = blocks.clone();
        altered[0].controls[0][0] += 1;
        assert!(audit_section(family, &altered).is_err());
        let mut contained = blocks.clone();
        let mut extra = blocks[0].clone();
        extra.id.push_str("_contained");
        contained.push(extra);
        assert!(audit_section(family, &contained).is_err());
        assert!(audit_section(family, &blocks[1..]).is_err());
    }
}

#[test]
fn saved_topology_and_boundary_labels_must_match_the_audited_geometry() {
    for family in SECTION_FAMILIES {
        let recipe = object([("family", Json::String(family.into()))]);
        let model = build(&recipe).unwrap();
        assert_eq!(
            audit(&model).unwrap().field("status").unwrap().text(),
            Ok("valid")
        );
        for field in ["topology", "design_faces"] {
            let mut changed = model.clone();
            let Json::Object(fields) = &mut changed else {
                panic!()
            };
            fields.insert(field.to_owned(), Json::Null);
            assert!(audit(&changed).is_err());
        }
    }
}
