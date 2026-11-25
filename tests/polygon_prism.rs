use solidcraft_trivariate::polygon_prism::{build_polygon_prism, PolygonPrism, Vertex};

fn extrude(points: &[[i32; 2]]) -> PolygonPrism {
    PolygonPrism {
        vertices: points
            .iter()
            .enumerate()
            .map(|(i, &xy_mm)| Vertex {
                id: format!("v{i}"),
                xy_mm,
            })
            .collect(),
        height_mm: 4,
        top_shift_mm: [-1, 5],
    }
}

#[test]
fn cyclic_and_reversed_contours_emit_identical_quads() {
    let mut source = extrude(&[[0, 0], [11, 0], [11, 4], [3, 4], [3, 13], [0, 13]]);
    let initial = build_polygon_prism(&source).unwrap();
    assert_eq!(initial.len(), 12);
    for block in &initial {
        for pair in block.controls.as_chunks::<2>().0 {
            assert_eq!(pair[1], [pair[0][0] - 6, pair[0][1] + 30, 24, 6]);
        }
    }
    for _ in 0..source.vertices.len() {
        source.vertices.rotate_left(1);
        assert_eq!(build_polygon_prism(&source).unwrap(), initial);
    }
    source.vertices.reverse();
    assert_eq!(build_polygon_prism(&source).unwrap(), initial);
}

#[test]
fn rejects_crossing_and_touched_contours_before_tessellation() {
    for pts in [
        vec![[0, 0], [8, 8], [0, 8], [8, 0]],
        vec![[0, 0], [10, 0], [10, 8], [5, 0], [0, 8]],
        vec![[1, 2], [3, 2], [6, 2], [1, 8]],
    ] {
        assert!(build_polygon_prism(&extrude(&pts)).is_err());
    }
}
