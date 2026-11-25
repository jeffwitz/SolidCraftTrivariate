use solidcraft_trivariate::json::{parse, Json};
use solidcraft_trivariate::perforated_prism::{
    build_perforated_prism, LoopProfile, PerforatedPrism,
};
use solidcraft_trivariate::polygon_prism::Vertex;

fn lp(id: &str, points: &[(&str, i32, i32)]) -> LoopProfile {
    LoopProfile {
        id: id.into(),
        vertices: points
            .iter()
            .map(|(name, x, y)| Vertex {
                id: (*name).into(),
                xy_mm: [*x, *y],
            })
            .collect(),
    }
}
fn recipe(holes: Vec<LoopProfile>) -> PerforatedPrism {
    PerforatedPrism {
        outer: lp(
            "outer",
            &[
                ("a", -10, -10),
                ("b", 10, -10),
                ("c", 10, 10),
                ("d", -10, 10),
            ],
        ),
        holes,
        height_mm: 8,
        top_shift_mm: [2, -1],
    }
}

#[test]
fn square_with_one_hole_has_eight_triangles_and_exact_ids() {
    let model = recipe(vec![lp(
        "h",
        &[("p", -2, -2), ("q", -2, 2), ("r", 2, 2), ("s", 2, -2)],
    )]);
    let blocks = build_perforated_prism(&model).unwrap();
    assert_eq!(blocks.len(), 24);
    assert!(blocks
        .iter()
        .all(|b| b.degrees == [1, 1, 1] && b.controls.len() == 8));
    assert!(blocks.iter().all(|b| b.id.starts_with("tri/")));
}

#[test]
fn square_without_holes_matches_g17_count() {
    let blocks = build_perforated_prism(&recipe(vec![])).unwrap();
    assert_eq!(blocks.len(), 6);
}

#[test]
fn touching_hole_is_rejected() {
    let model = recipe(vec![lp(
        "h",
        &[("p", -2, -2), ("q", -2, 2), ("r", 2, 2), ("s", 2, -2)],
    )]);
    let mut bad = model;
    bad.holes[0].vertices[0].xy_mm = [-10, -10];
    assert!(build_perforated_prism(&bad).is_err());
}

#[test]
fn contour_cycling_and_reversal_preserve_block_ids_and_controls() {
    let original = recipe(vec![lp(
        "h",
        &[("p", -2, -2), ("q", -2, 2), ("r", 2, 2), ("s", 2, -2)],
    )]);
    let mut changed = original.clone();
    changed.outer.vertices.reverse();
    changed.outer.vertices.rotate_left(1);
    changed.holes[0].vertices.reverse();
    changed.holes[0].vertices.rotate_left(2);
    assert_eq!(
        build_perforated_prism(&original).unwrap(),
        build_perforated_prism(&changed).unwrap()
    );
}

#[test]
fn four_disjoint_triangular_holes_are_triangulated() {
    let holes = vec![
        lp("h0", &[("a0", -8, -7), ("b0", -5, -7), ("c0", -7, -4)]),
        lp("h1", &[("a1", 5, -7), ("b1", 8, -7), ("c1", 7, -4)]),
        lp("h2", &[("a2", -8, 7), ("b2", -7, 4), ("c2", -5, 7)]),
        lp("h3", &[("a3", 5, 7), ("b3", 7, 4), ("c3", 8, 7)]),
    ];
    let blocks = build_perforated_prism(&recipe(holes)).unwrap();
    // N=16, H=4 => 22 triangles, each represented by three quads.
    assert_eq!(blocks.len(), 66);
}

#[test]
fn hole_order_and_outer_winding_do_not_change_the_net() {
    let mut source = recipe(vec![
        lp(
            "left",
            &[("l0", -8, -2), ("l1", -5, -2), ("l2", -5, 2), ("l3", -8, 2)],
        ),
        lp(
            "right",
            &[("r0", 5, -2), ("r1", 8, -2), ("r2", 8, 2), ("r3", 5, 2)],
        ),
    ]);
    let expected = build_perforated_prism(&source).unwrap();
    source.holes.reverse();
    source.outer.vertices.reverse();
    source.holes[0].vertices.rotate_left(1);
    assert_eq!(build_perforated_prism(&source).unwrap(), expected);
}

#[test]
fn concave_outer_with_oblique_hole_builds() {
    let mut source = recipe(vec![lp(
        "oblique",
        &[("p", 4, 1), ("q", 6, 4), ("r", 8, 2)],
    )]);
    source.outer.vertices = [
        ("a", -10, -10),
        ("b", 10, -10),
        ("c", 10, 10),
        ("d", 2, 10),
        ("e", 2, 0),
        ("f", -2, 0),
        ("g", -2, 10),
        ("h", -10, 10),
    ]
    .into_iter()
    .map(|(id, x, y)| Vertex {
        id: id.into(),
        xy_mm: [x, y],
    })
    .collect();
    assert!(
        build_perforated_prism(&source).is_ok(),
        "{:?}",
        build_perforated_prism(&source)
    );
}

fn parse_loop(value: &Json) -> LoopProfile {
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

#[test]
fn frozen_corpus_admission_and_equivalence_cases() {
    let corpus = parse(include_str!("references/perforated_prism_corpus.json")).unwrap();
    let cases = corpus.items().unwrap();
    let mut outputs = std::collections::BTreeMap::new();
    for case in cases {
        let value = case.field("recipe").unwrap();
        let outer = parse_loop(value.field("outer").unwrap());
        let holes = value
            .field("holes")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .map(parse_loop)
            .collect();
        let shift = value.field("top_shift_mm").unwrap().items().unwrap();
        let model = PerforatedPrism {
            outer,
            holes,
            height_mm: value.field("height_mm").unwrap().number().unwrap(),
            top_shift_mm: [shift[0].number().unwrap(), shift[1].number().unwrap()],
        };
        let expected = case.field("expected").unwrap() == &Json::Bool(true);
        let result = build_perforated_prism(&model);
        assert_eq!(
            result.is_ok(),
            expected,
            "{}: {:?}",
            case.field("id").unwrap().text().unwrap(),
            result.as_ref().err()
        );
        if let Ok(blocks) = result {
            outputs.insert(case.field("id").unwrap().text().unwrap().to_owned(), blocks);
        }
    }
    for case in cases {
        if let Ok(name) = case.field("equivalent_to").and_then(Json::text) {
            let id = case.field("id").unwrap().text().unwrap();
            assert_eq!(
                outputs.get(id),
                outputs.get(name),
                "equivalence {id} -> {name}"
            );
        }
    }
}
