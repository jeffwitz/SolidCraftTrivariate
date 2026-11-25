use solidcraft_trivariate::{
    circular_polygon::{build_circular_polygon_prism, CircularPolygonBore, CircularPolygonPrism},
    incidence::enumerate_faces,
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
fn frozen_corpus_admission_and_whole_face_incidence() {
    let corpus = parse(include_str!("references/circular_polygon_corpus.json")).unwrap();
    for case in corpus.items().unwrap() {
        let source = recipe(case.field("recipe").unwrap());
        let expected = case.field("expected").unwrap() == &Json::Bool(true);
        let built = build_circular_polygon_prism(&source);
        assert_eq!(
            built.is_ok(),
            expected,
            "{}: {:?}",
            case.field("id").unwrap().text().unwrap(),
            built.as_ref().err()
        );
        if let Ok(blocks) = built {
            let n = source.outer.vertices.len() + 4 * source.bores.len();
            let triangles = n + 2 * source.bores.len() - 2;
            assert_eq!(blocks.len(), 3 * triangles + 8 * source.bores.len());
            assert!(blocks.iter().all(|b| b.controls.iter().all(|p| p[3] > 0)));
            assert!(blocks
                .iter()
                .filter(|b| b.id.starts_with("bore_"))
                .all(|b| b.degrees == [2, 1, 1] && b.controls.len() == 12));
            let incidence = enumerate_faces(&blocks).unwrap_or_else(|e| {
                panic!(
                    "{} incidence: {e}",
                    case.field("id").unwrap().text().unwrap()
                )
            });
            assert_eq!(
                2 * incidence.interfaces.len() + incidence.exterior.len(),
                6 * blocks.len()
            );
            assert_eq!(
                incidence
                    .interfaces
                    .iter()
                    .filter(|i| i.first.block.starts_with("bore_")
                        || i.second.block.starts_with("bore_"))
                    .count(),
                16 * source.bores.len()
            );
        }
    }
}

#[test]
fn bore_order_and_outer_winding_preserve_the_partition() {
    let corpus = parse(include_str!("references/circular_polygon_corpus.json")).unwrap();
    let original = recipe(corpus.items().unwrap()[1].field("recipe").unwrap());
    let expected = build_circular_polygon_prism(&original).unwrap();
    let mut changed = original.clone();
    changed.outer.vertices.reverse();
    changed.outer.vertices.rotate_left(1);
    changed.bores.reverse();
    assert_eq!(build_circular_polygon_prism(&changed).unwrap(), expected);
}

#[test]
fn collisions_and_guard_limit_are_rejected() {
    let corpus = parse(include_str!("references/circular_polygon_corpus.json")).unwrap();
    let mut bad = recipe(corpus.items().unwrap()[0].field("recipe").unwrap());
    bad.bores[0].guard_mm = 2 * bad.bores[0].radius_mm - 1;
    assert!(build_circular_polygon_prism(&bad).is_err());
    let mut bad = recipe(corpus.items().unwrap()[0].field("recipe").unwrap());
    bad.bores[0].center_mm = [64, 64];
    assert!(build_circular_polygon_prism(&bad).is_err());
    let mut bad = recipe(corpus.items().unwrap()[0].field("recipe").unwrap());
    bad.bores[0].id = bad.outer.vertices[0].id.clone();
    assert!(build_circular_polygon_prism(&bad).is_err());
}
