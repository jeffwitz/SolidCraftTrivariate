use solidcraft_trivariate::{
    counterbores::{build_counterbores, Counterbore, CounterboreOpening, CounterborePrism},
    counterbores_audit::audit_counterbores,
    perforated_prism::LoopProfile,
    polygon_prism::Vertex,
};

fn recipe(id: &str, opening: CounterboreOpening) -> CounterborePrism {
    CounterborePrism {
        outer: LoopProfile {
            id: "outer".into(),
            vertices: [[-24, -16], [24, -16], [24, 16], [-24, 16]]
                .into_iter()
                .enumerate()
                .map(|(i, xy_mm)| Vertex {
                    id: format!("v{i}"),
                    xy_mm,
                })
                .collect(),
        },
        bores: vec![Counterbore {
            id: id.into(),
            center_mm: [0, 0],
            radius_mm: 3,
            guard_mm: 12,
            counterbore_radius_mm: 6,
            depth_mm: 4,
            opening,
        }],
        height_mm: 12,
        top_shift_mm: [3, -2],
    }
}

fn feature_count(
    faces: &[(solidcraft_trivariate::incidence::BlockSide, String)],
    feature: &str,
) -> usize {
    faces.iter().filter(|(_, label)| label == feature).count()
}

#[test]
fn audits_start_and_end_openings_with_exact_volume_and_faces() {
    for (id, opening) in [
        ("start", CounterboreOpening::Start),
        ("end", CounterboreOpening::End),
        ("a_sector_b", CounterboreOpening::Start),
    ] {
        let recipe = recipe(id, opening);
        let blocks = build_counterbores(&recipe).unwrap();
        let audit = audit_counterbores(&recipe, &blocks).unwrap();
        assert_eq!(audit.removed_pi_volume_coefficient_mm3, 216);
        assert_eq!(audit.outer_area_twice_mm2, 3072);
        assert_eq!(audit.through_bores, 1);
        assert_eq!(audit.enclosed_voids, 0);
        assert!(audit.certificate.determinant_lower_bound_mm3 > 0.0);
        assert_eq!(
            feature_count(
                &audit.certificate.design_faces,
                &format!("bore/{id}/pilot_wall")
            ),
            8
        );
        assert_eq!(
            feature_count(
                &audit.certificate.design_faces,
                &format!("bore/{id}/counterbore_wall")
            ),
            8
        );
        assert_eq!(
            feature_count(
                &audit.certificate.design_faces,
                &format!("bore/{id}/shoulder")
            ),
            8
        );

        let mut permuted = blocks.clone();
        permuted.reverse();
        assert_eq!(audit_counterbores(&recipe, &permuted).unwrap(), audit);
    }
}

#[test]
fn rejects_missing_and_altered_controls() {
    let recipe = recipe("bore", CounterboreOpening::Start);
    let blocks = build_counterbores(&recipe).unwrap();
    let mut missing = blocks.clone();
    missing.pop();
    assert!(audit_counterbores(&recipe, &missing).is_err());

    let mut altered = blocks;
    let block = altered
        .iter_mut()
        .find(|block| block.id.contains("pilot_sector"))
        .unwrap();
    block.controls[0][0] += 1;
    assert!(audit_counterbores(&recipe, &altered).is_err());
}

#[test]
fn audits_mixed_openings_and_coincident_layer_partition() {
    let mut recipe = recipe("left", CounterboreOpening::Start);
    recipe.bores = vec![
        Counterbore {
            id: "left".into(),
            center_mm: [-8, 0],
            radius_mm: 2,
            guard_mm: 6,
            counterbore_radius_mm: 3,
            depth_mm: 3,
            opening: CounterboreOpening::Start,
        },
        Counterbore {
            id: "right".into(),
            center_mm: [8, 0],
            radius_mm: 2,
            guard_mm: 6,
            counterbore_radius_mm: 3,
            depth_mm: 5,
            opening: CounterboreOpening::End,
        },
    ];
    let blocks = build_counterbores(&recipe).unwrap();
    let audit = audit_counterbores(&recipe, &blocks).unwrap();
    assert_eq!(audit.removed_pi_volume_coefficient_mm3, 136);
    assert_eq!(audit.through_bores, 2);
    assert_eq!(
        feature_count(&audit.certificate.design_faces, "bore/left/pilot_wall"),
        16
    );
    assert_eq!(
        feature_count(
            &audit.certificate.design_faces,
            "bore/right/counterbore_wall"
        ),
        8
    );
    assert_eq!(
        feature_count(&audit.certificate.design_faces, "bore/left/shoulder"),
        8
    );
    assert_eq!(
        feature_count(&audit.certificate.design_faces, "bore/right/shoulder"),
        8
    );
}
