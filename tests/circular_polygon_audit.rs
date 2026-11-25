use solidcraft_trivariate::{
    circular_polygon::{build_circular_polygon_prism, CircularPolygonBore, CircularPolygonPrism},
    circular_polygon_audit::audit_circular_polygon_prism,
    json::{parse, Json},
    perforated_prism::LoopProfile,
    polygon_prism::Vertex,
};

fn loop_profile(value: &Json) -> LoopProfile {
    let vertices = value
        .field("vertices")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|v| {
            let xy = v.field("xy_mm").unwrap().items().unwrap();
            Vertex {
                id: v.field("id").unwrap().text().unwrap().to_owned(),
                xy_mm: [xy[0].number().unwrap(), xy[1].number().unwrap()],
            }
        })
        .collect();
    LoopProfile {
        id: value.field("id").unwrap().text().unwrap().to_owned(),
        vertices,
    }
}

fn recipe(value: &Json) -> CircularPolygonPrism {
    let bores = value
        .field("bores")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|b| {
            let center = b.field("center_mm").unwrap().items().unwrap();
            CircularPolygonBore {
                id: b.field("id").unwrap().text().unwrap().to_owned(),
                center_mm: [center[0].number().unwrap(), center[1].number().unwrap()],
                radius_mm: b.field("radius_mm").unwrap().number().unwrap(),
                guard_mm: b.field("guard_mm").unwrap().number().unwrap(),
            }
        })
        .collect();
    let shift = value.field("top_shift_mm").unwrap().items().unwrap();
    CircularPolygonPrism {
        outer: loop_profile(value.field("outer").unwrap()),
        bores,
        height_mm: value.field("height_mm").unwrap().number().unwrap(),
        top_shift_mm: [shift[0].number().unwrap(), shift[1].number().unwrap()],
    }
}

#[test]
fn continuous_audit_accepts_frozen_circular_bore_corpus() {
    let corpus = parse(include_str!("references/circular_polygon_corpus.json")).unwrap();
    for case in corpus.items().unwrap() {
        let r = recipe(case.field("recipe").unwrap());
        if case.field("expected").unwrap() != &Json::Bool(true) {
            assert!(build_circular_polygon_prism(&r).is_err());
            continue;
        }
        let blocks = build_circular_polygon_prism(&r).unwrap();
        let cert = audit_circular_polygon_prism(&r, &blocks)
            .unwrap_or_else(|e| panic!("{} audit: {e}", case.field("id").unwrap().text().unwrap()));
        assert!(cert.certificate.determinant_lower_bound_mm3 > 0.0);
        assert_eq!(cert.certificate.radial_error_upper_mm, 0.0);
        let expected_area2 = r
            .outer
            .vertices
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let b = &r.outer.vertices[(i + 1) % r.outer.vertices.len()];
                i64::from(a.xy_mm[0]) * i64::from(b.xy_mm[1])
                    - i64::from(a.xy_mm[1]) * i64::from(b.xy_mm[0])
            })
            .sum::<i64>()
            .abs();
        assert_eq!(cert.outer_area_twice_mm2, expected_area2);
        assert_eq!(
            cert.sum_squared_radii_mm2,
            r.bores
                .iter()
                .map(|b| i64::from(b.radius_mm).pow(2))
                .sum::<i64>()
        );
        for bore in &r.bores {
            assert_eq!(
                cert.certificate
                    .design_faces
                    .iter()
                    .filter(|(_, tag)| tag == &format!("bore/{}", bore.id))
                    .count(),
                8
            );
        }
    }
}

#[test]
fn independent_audit_rejects_coefficient_topology_and_boundary_mutations() {
    let corpus = parse(include_str!("references/circular_polygon_corpus.json")).unwrap();
    let r = recipe(corpus.items().unwrap()[0].field("recipe").unwrap());
    let original = build_circular_polygon_prism(&r).unwrap();
    assert!(audit_circular_polygon_prism(&r, &original).is_ok());

    let mut changed = original.clone();
    let circle = changed
        .iter_mut()
        .find(|b| b.id.starts_with("bore_"))
        .unwrap();
    circle.controls[2][0] += 1;
    assert!(audit_circular_polygon_prism(&r, &changed).is_err());

    let mut changed = original.clone();
    changed.pop();
    assert!(audit_circular_polygon_prism(&r, &changed).is_err());

    let mut changed = original.clone();
    let outer = changed
        .iter_mut()
        .find(|b| b.id.starts_with("tri/"))
        .unwrap();
    outer.controls[0][0] += 1;
    assert!(audit_circular_polygon_prism(&r, &changed).is_err());

    let mut changed = original.clone();
    let ring = changed
        .iter_mut()
        .find(|b| b.id.starts_with("bore_"))
        .unwrap();
    ring.degrees[0] = 1;
    assert!(audit_circular_polygon_prism(&r, &changed).is_err());
}

#[test]
fn exact_report_is_invariant_to_bore_input_order_and_outer_winding() {
    let corpus = parse(include_str!("references/circular_polygon_corpus.json")).unwrap();
    let r = recipe(corpus.items().unwrap()[1].field("recipe").unwrap());
    let first = build_circular_polygon_prism(&r).unwrap();
    let first_cert = audit_circular_polygon_prism(&r, &first).unwrap();
    let mut reordered = r.clone();
    reordered.bores.reverse();
    reordered.outer.vertices.reverse();
    reordered.outer.vertices.rotate_left(1);
    let second = build_circular_polygon_prism(&reordered).unwrap();
    let second_cert = audit_circular_polygon_prism(&reordered, &second).unwrap();
    assert_eq!(first, second);
    assert_eq!(first_cert, second_cert);
}

#[test]
fn bore_ids_containing_sector_token_and_guard_prefixes_are_unambiguous() {
    let corpus = parse(include_str!("references/circular_polygon_corpus.json")).unwrap();
    let mut r = recipe(corpus.items().unwrap()[0].field("recipe").unwrap());
    r.bores[0].id = "a_sector_b".to_owned();
    let blocks = build_circular_polygon_prism(&r).unwrap();
    let certificate = audit_circular_polygon_prism(&r, &blocks).unwrap();
    assert_eq!(
        certificate
            .certificate
            .design_faces
            .iter()
            .filter(|(_, label)| label == "bore/a_sector_b")
            .count(),
        8
    );

    let mut r = recipe(corpus.items().unwrap()[1].field("recipe").unwrap());
    r.bores[0].id = "a".to_owned();
    r.bores[1].id = "ab".to_owned();
    let blocks = build_circular_polygon_prism(&r).unwrap();
    let certificate = audit_circular_polygon_prism(&r, &blocks).unwrap();
    for id in ["a", "ab"] {
        assert_eq!(
            certificate
                .certificate
                .design_faces
                .iter()
                .filter(|(_, label)| label == &format!("bore/{id}"))
                .count(),
            8
        );
    }
    assert!(certificate
        .certificate
        .design_faces
        .iter()
        .all(|(_, label)| !label.starts_with("edge/guard_")));

    let mut collision = r.clone();
    collision.outer.vertices[0].id = collision.outer.id.clone();
    assert!(audit_circular_polygon_prism(&collision, &blocks).is_err());
}
