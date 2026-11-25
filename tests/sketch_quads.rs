use solidcraft_trivariate::{
    constraint_system::Ratio,
    sketch_quads::{quadrangulate, quadrangulate_with_provenance, BoundaryEdgeSource},
};

type Point = [Ratio; 2];

fn point(x: i128, y: i128) -> Point {
    [Ratio::new(x, 1).unwrap(), Ratio::new(y, 1).unwrap()]
}

fn ratio_add(a: Ratio, b: Ratio) -> Ratio {
    Ratio::new(
        a.numerator() * b.denominator() + b.numerator() * a.denominator(),
        a.denominator() * b.denominator(),
    )
    .unwrap()
}

fn ratio_sub(a: Ratio, b: Ratio) -> Ratio {
    Ratio::new(
        a.numerator() * b.denominator() - b.numerator() * a.denominator(),
        a.denominator() * b.denominator(),
    )
    .unwrap()
}

fn ratio_mul(a: Ratio, b: Ratio) -> Ratio {
    Ratio::new(
        a.numerator() * b.numerator(),
        a.denominator() * b.denominator(),
    )
    .unwrap()
}

fn orient(a: Point, b: Point, c: Point) -> Ratio {
    ratio_sub(
        ratio_mul(ratio_sub(b[0], a[0]), ratio_sub(c[1], a[1])),
        ratio_mul(ratio_sub(b[1], a[1]), ratio_sub(c[0], a[0])),
    )
}

fn area2(points: &[Point]) -> Ratio {
    let mut area = Ratio::zero();
    for i in 0..points.len() {
        let a = points[i];
        let b = points[(i + 1) % points.len()];
        area = ratio_add(
            area,
            ratio_sub(ratio_mul(a[0], b[1]), ratio_mul(a[1], b[0])),
        );
    }
    area
}

fn assert_positive_convex(quad: &[[Ratio; 2]; 4]) {
    for i in 0..4 {
        assert!(orient(quad[i], quad[(i + 1) % 4], quad[(i + 2) % 4]).numerator() > 0);
    }
}

fn assert_exact_area(face: &[Point]) {
    let quads = quadrangulate(face).unwrap();
    let total = quads
        .iter()
        .fold(Ratio::zero(), |sum, quad| ratio_add(sum, area2(quad)));
    assert_eq!(total, area2(face));
    for quad in &quads {
        assert_positive_convex(quad);
    }
}

#[test]
fn triangle_and_convex_quad_are_partitioned_into_positive_quads() {
    let triangle = [point(0, 0), point(4, 0), point(0, 3)];
    let tri_quads = quadrangulate(&triangle).unwrap();
    assert_eq!(tri_quads.len(), 3);
    assert_exact_area(&triangle);

    let square = [point(0, 0), point(4, 0), point(4, 4), point(0, 4)];
    let square_quads = quadrangulate(&square).unwrap();
    assert_eq!(square_quads.len(), 6);
    assert_exact_area(&square);
}

#[test]
fn concave_faces_and_collinear_boundary_vertices_preserve_area() {
    let l_shape = [
        point(0, 0),
        point(4, 0),
        point(4, 1),
        point(1, 1),
        point(1, 4),
        point(0, 4),
    ];
    let quads = quadrangulate(&l_shape).unwrap();
    assert_eq!(quads.len(), 12);
    assert_exact_area(&l_shape);

    let collinear = [
        point(0, 0),
        point(2, 0),
        point(4, 0),
        point(4, 2),
        point(0, 2),
    ];
    assert_eq!(quadrangulate(&collinear).unwrap().len(), 9);
    assert_exact_area(&collinear);
}

fn boundary_segments(
    quads: &[solidcraft_trivariate::sketch_quads::QuadWithProvenance],
    edge_index: usize,
) -> Vec<[Point; 2]> {
    let mut segments = Vec::new();
    for quad in quads {
        for edge in 0..4 {
            if quad.edge_sources[edge].is_some_and(|source| source.edge_index == edge_index) {
                let mut endpoints = [quad.vertices[edge], quad.vertices[(edge + 1) % 4]];
                if solidcraft_trivariate::sketch_quads::compare_points_exact(
                    endpoints[0],
                    endpoints[1],
                )
                .is_gt()
                {
                    endpoints.swap(0, 1);
                }
                segments.push(endpoints);
            }
        }
    }
    segments.sort_by(|a, b| {
        solidcraft_trivariate::sketch_quads::compare_points_exact(a[0], b[0])
            .then_with(|| solidcraft_trivariate::sketch_quads::compare_points_exact(a[1], b[1]))
    });
    segments
}

fn internal_segments(
    quads: &[solidcraft_trivariate::sketch_quads::QuadWithProvenance],
) -> Vec<[Point; 2]> {
    let mut segments = Vec::new();
    for quad in quads {
        for edge in 0..4 {
            if quad.edge_sources[edge].is_none() {
                let mut endpoints = [quad.vertices[edge], quad.vertices[(edge + 1) % 4]];
                if solidcraft_trivariate::sketch_quads::compare_points_exact(
                    endpoints[0],
                    endpoints[1],
                )
                .is_gt()
                {
                    endpoints.swap(0, 1);
                }
                segments.push(endpoints);
            }
        }
    }
    segments.sort_by(|a, b| {
        solidcraft_trivariate::sketch_quads::compare_points_exact(a[0], b[0])
            .then_with(|| solidcraft_trivariate::sketch_quads::compare_points_exact(a[1], b[1]))
    });
    segments
}

#[test]
fn internal_quadrilateral_edges_are_paired_exactly_twice() {
    let face = [
        point(0, 0),
        point(5, 0),
        point(5, 2),
        point(2, 2),
        point(2, 5),
        point(0, 5),
    ];
    let quads = quadrangulate_with_provenance(&face).unwrap();
    let segments = internal_segments(&quads);
    assert_eq!(segments.len() % 2, 0);
    let (pairs, remainder) = segments.as_chunks::<2>();
    assert!(remainder.is_empty());
    for pair in pairs {
        assert_eq!(pair[0], pair[1]);
    }
}

#[test]
fn shared_edges_split_identically_and_keep_edge_provenance() {
    let left = [point(0, 0), point(2, 0), point(2, 2), point(0, 2)];
    let right = [point(2, 0), point(4, 0), point(4, 2), point(2, 2)];
    let left_quads = quadrangulate_with_provenance(&left).unwrap();
    let right_quads = quadrangulate_with_provenance(&right).unwrap();
    assert_eq!(
        boundary_segments(&left_quads, 1),
        boundary_segments(&right_quads, 3)
    );
    assert_eq!(boundary_segments(&left_quads, 1).len(), 2);
    for (quads, face) in [(&left_quads, &left), (&right_quads, &right)] {
        for edge in 0..face.len() {
            let sources = quads
                .iter()
                .flat_map(|quad| quad.edge_sources)
                .flatten()
                .filter(|source| source.edge_index == edge)
                .collect::<Vec<BoundaryEdgeSource>>();
            assert_eq!(sources.len(), 2);
            assert!(sources.iter().any(|source| source.half == 0));
            assert!(sources.iter().any(|source| source.half == 1));
        }
    }
}

#[test]
fn rejects_clockwise_self_intersecting_and_degenerate_faces() {
    assert!(quadrangulate(&[point(0, 0), point(0, 2), point(2, 2), point(2, 0)]).is_err());
    assert!(quadrangulate(&[point(0, 0), point(3, 3), point(0, 3), point(3, 0)]).is_err());
    assert!(quadrangulate(&[point(0, 0), point(1, 0), point(2, 0)]).is_err());
    assert!(quadrangulate(&[point(0, 0), point(1, 0)]).is_err());
}
