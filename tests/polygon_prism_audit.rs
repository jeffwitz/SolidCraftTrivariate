use solidcraft_trivariate::{
    polygon_prism::{build_polygon_prism, PolygonPrism, Vertex},
    polygon_prism_audit::audit_polygon_prism,
};

fn recipe(points: &[[i32; 2]], height: i32, shift: [i32; 2]) -> PolygonPrism {
    PolygonPrism {
        vertices: points
            .iter()
            .enumerate()
            .map(|(i, &xy_mm)| Vertex {
                id: format!("v{i}"),
                xy_mm,
            })
            .collect(),
        height_mm: height,
        top_shift_mm: shift,
    }
}

#[test]
fn exact_audit_covers_convex_concave_and_thin_profiles() {
    let cases = [
        (vec![[0, 0], [8, 0], [0, 6]], 48_i64),
        (vec![[0, 0], [8, 0], [10, 5], [4, 9], [-1, 5]], 139),
        (vec![[0, 0], [12, 0], [12, 3], [4, 3], [4, 10], [0, 10]], 128),
        (vec![[0, 0], [12, 0], [12, 10], [9, 10], [9, 3], [3, 3], [3, 10], [0, 10]], 156),
        (vec![[0, 0], [20, 0], [20, 1], [1, 1], [1, 20], [0, 20]], 78),
    ];
    for (points, expected_area2) in cases {
        let source = recipe(&points, 7, [3, -2]);
        let blocks = build_polygon_prism(&source).unwrap();
        let result = audit_polygon_prism(&source, &blocks).unwrap();
        assert_eq!(result.projected_area_twice_mm2, expected_area2);
        assert_eq!(result.volume_sixths_mm3, expected_area2 * 3 * 7);
        assert_eq!(result.certificate.incidence.exterior.len(), 8 * points.len() - 12);
        assert_eq!(result.certificate.incidence.interfaces.len(), 5 * points.len() - 12);
        assert_eq!(result.certificate.design_faces.len(), 8 * points.len() - 12);
        assert_eq!(result.certificate.design_faces.iter().filter(|(_, feature)| *feature == "start").count(), 3 * (points.len() - 2));
        assert_eq!(result.certificate.design_faces.iter().filter(|(_, feature)| *feature == "end").count(), 3 * (points.len() - 2));
        assert_eq!(result.certificate.design_faces.iter().filter(|(_, feature)| feature.starts_with("edge/")).count(), 2 * points.len());
        assert!(result.certificate.determinant_lower_bound_mm3 > 0.0);
    }
}

#[test]
fn exact_audit_is_invariant_to_cycle_and_winding() {
    let mut source = recipe(&[[0, 0], [9, 0], [11, 4], [5, 8], [0, 5]], 5, [-2, 1]);
    let blocks = build_polygon_prism(&source).unwrap();
    let expected = audit_polygon_prism(&source, &blocks).unwrap();
    source.vertices.rotate_left(2);
    assert_eq!(audit_polygon_prism(&source, &blocks).unwrap(), expected);
    source.vertices.reverse();
    assert_eq!(audit_polygon_prism(&source, &blocks).unwrap(), expected);
}

#[test]
fn mutations_and_malformed_contours_are_rejected() {
    let source = recipe(&[[0, 0], [10, 0], [10, 8], [0, 8]], 6, [2, 3]);
    let original = build_polygon_prism(&source).unwrap();
    let mut altered = original.clone();
    altered[0].controls[0][0] += 1;
    assert!(audit_polygon_prism(&source, &altered).is_err());
    let mut altered = original.clone();
    altered[0].id.push_str("_bad");
    assert!(audit_polygon_prism(&source, &altered).is_err());
    let mut altered = original.clone();
    altered.pop();
    assert!(audit_polygon_prism(&source, &altered).is_err());

    for points in [
        vec![[0, 0], [8, 8], [0, 8], [8, 0]],
        vec![[0, 0], [10, 0], [10, 8], [5, 0], [0, 8]],
        vec![[0, 0], [4, 0], [8, 0], [0, 5]],
    ] {
        let bad = recipe(&points, 4, [0, 0]);
        assert!(audit_polygon_prism(&bad, &[]).is_err());
    }
}
