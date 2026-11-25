use solidcraft_trivariate::{
    perforated_prism::{build_perforated_prism, LoopProfile, PerforatedPrism},
    perforated_prism_audit::audit_perforated_prism,
    polygon_prism::Vertex,
};

fn lp(id: &str, points: &[[i32; 2]], prefix: &str) -> LoopProfile {
    LoopProfile {
        id: id.into(),
        vertices: points
            .iter()
            .enumerate()
            .map(|(i, &xy_mm)| Vertex {
                id: format!("{prefix}{i}"),
                xy_mm,
            })
            .collect(),
    }
}
fn source() -> PerforatedPrism {
    PerforatedPrism {
        outer: lp("outer", &[[0, 0], [20, 0], [20, 20], [0, 20]], "o"),
        holes: vec![
            lp("bore_a", &[[5, 5], [5, 10], [10, 10], [10, 5]], "a"),
            lp("bore_b", &[[13, 5], [13, 11], [17, 11], [17, 5]], "b"),
        ],
        height_mm: 8,
        top_shift_mm: [3, -2],
    }
}

#[test]
fn holes_are_open_through_tunnels_with_boundary_provenance() {
    let source = source();
    let blocks = build_perforated_prism(&source).unwrap();
    let result = audit_perforated_prism(&source, &blocks).unwrap();
    assert_eq!(blocks.len(), 3 * (12 + 2 * 2 - 2));
    // Outer area 400, holes 25 and 24 square mm.
    assert_eq!(result.projected_area_twice_mm2, 2 * (400 - 25 - 24));
    assert_eq!(result.volume_sixths_mm3, 3 * 2 * (400 - 25 - 24) * 8);
    let tags = result
        .certificate
        .design_faces
        .iter()
        .map(|(_, tag)| tag.as_str())
        .collect::<Vec<_>>();
    assert!(tags.iter().any(|tag| tag.starts_with("edge/bore_a/")));
    assert!(tags.iter().any(|tag| tag.starts_with("edge/bore_b/")));
    assert!(tags
        .iter()
        .all(|tag| *tag == "start" || *tag == "end" || tag.starts_with("edge/")));
    assert_eq!(
        result.certificate.incidence.exterior.len(),
        result.certificate.design_faces.len()
    );
    assert_eq!(
        result.certificate.incidence.interfaces.len(),
        3 * (12 + 2 * 2 - 2) + 2 * (12 + 3 * 2 - 3)
    );
}

#[test]
fn zero_hole_path_matches_simple_polygon_core() {
    let outer = vec![[0, 0], [12, 0], [12, 4], [5, 4], [5, 11], [0, 11]];
    let perforated = PerforatedPrism {
        outer: lp("outline", &outer, "p"),
        holes: vec![],
        height_mm: 6,
        top_shift_mm: [-2, 3],
    };
    let blocks = build_perforated_prism(&perforated).unwrap();
    let result = audit_perforated_prism(&perforated, &blocks).unwrap();
    assert_eq!(result.projected_area_twice_mm2, 2 * (12 * 4 + 5 * 7));
    assert_eq!(
        result.volume_sixths_mm3,
        3 * result.projected_area_twice_mm2 * 6
    );
    assert!(result
        .certificate
        .design_faces
        .iter()
        .all(|(_, tag)| tag == "start" || tag == "end" || tag.starts_with("edge/outline/")));
}

#[test]
fn rejects_mutations_and_invalid_loop_relationships() {
    let source = source();
    let blocks = build_perforated_prism(&source).unwrap();
    let mut changed = blocks.clone();
    changed[0].controls[0][0] += 1;
    assert!(audit_perforated_prism(&source, &changed).is_err());
    let mut missing = blocks;
    missing.pop();
    assert!(audit_perforated_prism(&source, &missing).is_err());

    let mut outside = source.clone();
    outside.holes[0] = lp("outside", &[[22, 4], [22, 8], [26, 8], [26, 4]], "x");
    assert!(audit_perforated_prism(&outside, &[]).is_err());
    let mut nested = source.clone();
    nested.holes[1] = lp("nested", &[[6, 6], [6, 8], [8, 8], [8, 6]], "n");
    assert!(audit_perforated_prism(&nested, &[]).is_err());
    let mut contact = source;
    contact.holes[0] = lp("contact", &[[0, 3], [0, 6], [3, 6], [3, 3]], "c");
    assert!(audit_perforated_prism(&contact, &[]).is_err());
}
