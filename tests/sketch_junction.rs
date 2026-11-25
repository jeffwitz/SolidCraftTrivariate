use solidcraft_trivariate::sketch_junction::{classify, Continuity::*};

#[test]
fn one_sided_geometry_distinguishes_disconnect_corner_and_stationary_end() {
    let left = [[0, 0, 1], [1, 0, 1]];
    assert_eq!(classify(&left, &[[2, 0, 1], [3, 0, 1]]), Ok(Disconnected));
    assert_eq!(classify(&left, &[[1, 0, 1], [1, 1, 1]]), Ok(Corner));
    assert_eq!(classify(&left, &[[1, 0, 1], [0, 0, 1]]), Ok(Corner));
    assert_eq!(
        classify(&left, &[[1, 0, 1], [1, 0, 1], [2, 0, 1]]),
        Ok(Singular)
    );
    assert_eq!(classify(&left, &[[2, 0, 2], [30, 0, 3]]), Ok(Curvature));
    assert!(classify(&left, &[[1, 0, 0], [2, 0, 1]]).is_err());
}

#[test]
fn unequal_parameter_speeds_and_tangential_acceleration_preserve_curvature() {
    let left = [[-2, 2, 1], [-1, 0, 1], [0, 0, 1]];
    let right = [[0, 0, 1], [2, 0, 1], [5, 8, 1]];
    assert_eq!(classify(&left, &right), Ok(Curvature));
    // Positive geometric progressions in the homogeneous gauges apply
    // distinct non-affine Möbius reparameterizations to the same two curves.
    let weighted_left = [[-2, 2, 1], [-2, 0, 2], [0, 0, 4]];
    let weighted_right = [[0, 0, 3], [12, 0, 6], [60, 96, 12]];
    assert_eq!(classify(&weighted_left, &weighted_right), Ok(Curvature));
    assert_eq!(
        classify(&left, &[[0, 0, 1], [2, 0, 1], [5, 7, 1]]),
        Ok(Tangent)
    );
}

#[test]
fn multiplicity_three_can_hide_a_smooth_split_on_a_rectangle_side() {
    use solidcraft_trivariate::nurbs_curve::extract_bezier;
    let corners = [[0, 0], [3, 0], [6, 0], [6, 3], [0, 3], [0, 0]];
    let mut controls = vec![[0, 0, 1]];
    for segment in corners.windows(2) {
        for i in 1..=3 {
            controls.push([
                ((3 - i) * segment[0][0] + i * segment[1][0]) / 3,
                ((3 - i) * segment[0][1] + i * segment[1][1]) / 3,
                1,
            ]);
        }
    }
    let knots = [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5];
    let spans = extract_bezier(3, &knots, &controls).unwrap();
    let kinds: Vec<_> = (0..spans.len())
        .map(|i| classify(&spans[i].controls, &spans[(i + 1) % spans.len()].controls).unwrap())
        .collect();
    assert_eq!(kinds, vec![Curvature, Corner, Corner, Corner, Corner]);
}

#[test]
fn nonuniform_periodic_input_is_smooth_at_the_seam_and_all_three_knots() {
    use solidcraft_trivariate::sketch_junction::plan_loop;
    let poles = [
        [0, 0, 1],
        [8, 0, 2],
        [4, 4, 1],
        [0, 8, 2],
        [0, 0, 1],
        [8, 0, 2],
        [4, 4, 1],
    ];
    // Intervals [2,1,3,1] repeat one full period of length seven.
    let knots = [-5, -4, -1, 0, 2, 3, 6, 7, 9, 10, 13];
    let plan = plan_loop(3, &knots, &poles).unwrap();
    assert_eq!(plan.len(), 4);
    assert!(plan[0].closure);
    assert!(plan
        .iter()
        .all(|j| j.multiplicity == 1 && j.continuity == Curvature));
}

#[test]
fn active_span_blossoming_matches_clamped_knot_insertion_exactly() {
    use solidcraft_trivariate::nurbs_curve::{extract_active_bezier, extract_bezier};
    let knots = [0, 0, 0, 0, 1, 2, 2, 5, 5, 5, 5];
    let poles = [
        [0, 0, 1],
        [4, 2, 2],
        [2, 7, 1],
        [8, 6, 2],
        [5, 1, 1],
        [10, 0, 2],
        [6, 5, 1],
    ];
    assert_eq!(
        extract_active_bezier(3, &knots, &poles).unwrap(),
        extract_bezier(3, &knots, &poles).unwrap()
    );
}

#[test]
fn circle_and_ellipse_have_four_smooth_joins_despite_repeated_knots() {
    use solidcraft_trivariate::sketch_junction::plan_loop;
    // Exact rational quarter-circle (1-t², 2t, 1+t²). Rotate each quarter
    // and match its starting gauge to the preceding quarter's end gauge.
    let knots = [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4];
    for (rx, ry) in [(1, 1), (3, 2)] {
        let mut poles = vec![[rx, 0, 1]];
        for quarter in 0..4 {
            let scale = 1 << quarter;
            for (mut x, mut y, w) in [(1, 1, 1), (0, 2, 2)] {
                for _ in 0..quarter {
                    (x, y) = (-y, x);
                }
                poles.push([scale * rx * x, scale * ry * y, scale * w]);
            }
        }
        let plan = plan_loop(2, &knots, &poles).unwrap();
        assert_eq!(plan.len(), 4);
        assert!(plan.iter().all(|j| j.continuity == Curvature));
        assert!(plan[1..].iter().all(|j| j.multiplicity == 2));
    }
}
