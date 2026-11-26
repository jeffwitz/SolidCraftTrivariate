use solidcraft_trivariate::sketch_arrangement::Segment;
use solidcraft_trivariate::sketch_layers::{partition, SketchLayers, SketchLoop};

fn contour(id: &str, vertices: &[[i32; 2]]) -> SketchLoop {
    SketchLoop {
        id: id.into(),
        vertices: vertices.to_vec(),
    }
}

fn job(outer: SketchLoop, cuts: Vec<(SketchLoop, [i32; 2])>, height: i32) -> SketchLayers {
    SketchLayers {
        outer,
        cuts,
        guides: Vec::<Segment>::new(),
        height,
        shift: [0, 0],
        strategy: solidcraft_trivariate::sketch_layers::PlannerStrategy::Automatic,
        force_triangulation: false,
        operations: None,
    }
}

#[test]
fn even_convex_polygon_hole_uses_conforming_ring_and_exact_partition() {
    let input = job(
        contour("stock", &[[0, 0], [100, 0], [100, 80], [0, 80]]),
        vec![(
            contour(
                "hex",
                &[[30, 30], [40, 25], [50, 30], [50, 40], [40, 45], [30, 40]],
            ),
            [0, 10],
        )],
        10,
    );
    let result = partition(&input).expect("a certified generic ring candidate");
    assert!(matches!(
        result.method,
        "convex_host_separators" | "rectangular_hosts_x" | "rectangular_hosts_y"
    ));
    assert_eq!(result.volume.numerator(), 77_000);
    assert_eq!(result.volume.denominator(), 1);
    assert!(result.certificate.determinant_lower_bound_mm3 > 0.0);
}

#[test]
fn planar_ring_lookahead_keeps_count_primary_and_certificate_gates() {
    use solidcraft_trivariate::sketch_layers::{planar_quality_score, PlannerStrategy};

    let input = job(
        contour("stock", &[[0, 0], [120, 0], [120, 80], [0, 80]]),
        vec![
            (
                contour("left", &[[22, 24], [32, 20], [38, 30], [30, 38]]),
                [0, 10],
            ),
            (
                contour("right", &[[82, 18], [94, 24], [90, 38], [78, 32]]),
                [0, 10],
            ),
        ],
        10,
    );
    let baseline = partition(&input).expect("automatic exact partition");
    let mut experimental = input.clone();
    experimental.strategy = PlannerStrategy::RingPlanarLookahead;
    let candidate = partition(&experimental).expect("lookahead exact partition");

    assert!(candidate.blocks.len() <= baseline.blocks.len());
    assert_eq!(candidate.volume, baseline.volume);
    assert!(candidate.certificate.determinant_lower_bound_mm3 > 0.0);
    if candidate.blocks.len() == baseline.blocks.len() {
        let old_quality = planar_quality_score(&input, &baseline).unwrap();
        let new_quality = planar_quality_score(&experimental, &candidate).unwrap();
        assert!(new_quality <= old_quality * (1.0 + 1e-12));
    }
}

#[test]
fn oblique_bisector_recursively_assigns_disjoint_convex_holes() {
    let input = job(
        contour("stock", &[[0, 0], [100, 0], [100, 100], [0, 100]]),
        vec![
            (
                contour("lower_left", &[[21, 36], [29, 36], [29, 44], [21, 44]]),
                [0, 10],
            ),
            (
                contour("upper_right", &[[71, 56], [79, 56], [79, 64], [71, 64]]),
                [0, 10],
            ),
        ],
        10,
    );
    let result = partition(&input).expect("oblique separator is exact and conforming");
    assert!(matches!(
        result.method,
        "convex_host_separators" | "rectangular_hosts_x" | "rectangular_hosts_y"
    ));
    assert_eq!(result.volume.numerator(), 98_720);
    assert_eq!(result.volume.denominator(), 1);
}

#[test]
fn convex_even_blind_pocket_gets_a_conforming_floor_cap() {
    let input = job(
        contour("stock", &[[0, 0], [100, 0], [100, 80], [0, 80]]),
        vec![(
            contour(
                "hex_pocket",
                &[[30, 30], [40, 25], [50, 30], [50, 40], [40, 45], [30, 40]],
            ),
            [3, 10],
        )],
        10,
    );
    let result = partition(&input).expect("pocket floor cap is conforming");
    assert!(matches!(
        result.method,
        "convex_host_separators" | "rectangular_hosts_x" | "rectangular_hosts_y"
    ));
    assert_eq!(result.blocks.len(), 14);
    assert_eq!(result.volume.numerator(), 77_900);
    assert_eq!(result.volume.denominator(), 1);
}

#[test]
fn odd_convex_through_cut_uses_bounded_exact_ring_subdivision() {
    let input = job(
        contour("stock", &[[0, 0], [100, 0], [100, 100], [0, 100]]),
        vec![(
            contour(
                "pentagon",
                &[[30, 30], [40, 25], [50, 30], [45, 45], [32, 43]],
            ),
            [0, 10],
        )],
        10,
    );
    let result = partition(&input).expect("exact scaled subdivision admits the odd ring");
    assert!(result.certificate.determinant_lower_bound_mm3 > 0.0);
    assert_eq!(result.volume.numerator(), 97_175);
    assert_eq!(result.volume.denominator(), 1);
}
