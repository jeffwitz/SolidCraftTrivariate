use solidcraft_trivariate::{
    nurbs_curve::{extract_bezier, BezierSpan},
    nurbs_extrusion::{determinant_bound, extrude_strips},
};

fn paired_cubic() -> (Vec<BezierSpan>, Vec<BezierSpan>) {
    let weights = [1_i64, 2, 2, 1];
    let x = [0_i64, 1, 2, 3];
    let lower = (0..4).map(|i| [x[i] * weights[i], 0, weights[i]]).collect();
    let upper = (0..4)
        .map(|i| [x[i] * weights[i], 2 * weights[i], weights[i]])
        .collect();
    (
        vec![BezierSpan {
            interval: [0, 7],
            controls: lower,
        }],
        vec![BezierSpan {
            interval: [0, 7],
            controls: upper,
        }],
    )
}

#[test]
fn extrudes_a_cubic_rational_strip_with_positive_exact_bound() {
    let (lower, upper) = paired_cubic();
    let blocks = extrude_strips(&lower, &upper, 5).unwrap();
    assert_eq!(blocks.len(), 1);
    assert_eq!(blocks[0].degrees, [3, 1, 1]);
    assert_eq!(blocks[0].controls.len(), 16);
    assert!(determinant_bound(&blocks[0]).unwrap() > 0.0);
    for i in 0..4 {
        let bottom = blocks[0].controls[(i * 2) * 2];
        let top = blocks[0].controls[(i * 2) * 2 + 1];
        let upper_y = blocks[0].controls[(i * 2 + 1) * 2];
        assert_eq!(bottom[2], 0);
        assert_eq!(top[2], 5 * bottom[3]);
        assert_eq!(upper_y[1] - bottom[1], 2 * bottom[3]);
    }
}

#[test]
fn accepts_descending_parameterization_by_reversing_span_orientation() {
    let (mut lower, mut upper) = paired_cubic();
    lower[0].controls.reverse();
    upper[0].controls.reverse();
    let blocks = extrude_strips(&lower, &upper, 3).unwrap();
    assert!(determinant_bound(&blocks[0]).unwrap() > 0.0);
    assert_eq!(blocks[0].controls[0][0], 0);
}

#[test]
fn aligns_different_integer_homogeneous_gauges() {
    let (lower, mut upper) = paired_cubic();
    for point in &mut upper[0].controls {
        for coordinate in point {
            *coordinate *= 3;
        }
    }
    let blocks = extrude_strips(&lower, &upper, 4).unwrap();
    assert!(determinant_bound(&blocks[0]).unwrap() > 0.0);
}

#[test]
fn determinant_bound_rechecks_the_stored_extrusion_orientation() {
    let (lower, upper) = paired_cubic();
    let mut block = extrude_strips(&lower, &upper, 4).unwrap().remove(0);
    for u in 0..=3 {
        for w in 0..=1 {
            let a = (u * 2) * 2 + w;
            let b = (u * 2 + 1) * 2 + w;
            block.controls.swap(a, b);
        }
    }
    assert!(determinant_bound(&block).is_err());
}

#[test]
fn refuses_unsynchronized_or_nonseparated_curves() {
    let (lower, mut upper) = paired_cubic();
    upper[0].controls[1][0] += 1;
    assert!(extrude_strips(&lower, &upper, 5).is_err());

    let (lower, mut upper) = paired_cubic();
    upper[0].interval = [0, 8];
    assert!(extrude_strips(&lower, &upper, 5).is_err());

    let (lower, mut upper) = paired_cubic();
    upper[0].controls[1][1] = -10;
    assert!(extrude_strips(&lower, &upper, 5).is_err());
}

#[test]
fn admits_the_cubic_two_span_nurbs_recipe() {
    let knots = [-2, -2, -2, -2, 1, 4, 4, 4, 4];
    let lower_controls = [[0, 0, 1], [2, 0, 2], [6, 2, 3], [12, 3, 4], [20, 0, 5]];
    let upper_controls = [[0, 2, 1], [2, 6, 2], [6, 11, 3], [12, 11, 4], [20, 5, 5]];
    let lower = extract_bezier(3, &knots, &lower_controls).unwrap();
    let upper = extract_bezier(3, &knots, &upper_controls).unwrap();
    let blocks = extrude_strips(&lower, &upper, 5).unwrap();
    assert_eq!(blocks.len(), 2);
    assert!(blocks
        .iter()
        .all(|block| determinant_bound(block).unwrap() > 0.0));
}

#[test]
fn admits_all_three_frozen_positive_nurbs_strip_cases() {
    let cases = [
        (
            [0, 0, 0, 0, 1, 2, 2, 2, 2],
            [[0, 0, 1], [2, 0, 2], [2, 1, 1], [9, 3, 3], [4, 0, 1]],
            [[0, 2, 1], [2, 4, 2], [2, 3, 1], [9, 9, 3], [4, 2, 1]],
            8,
        ),
        (
            [-2, -2, -2, -2, 1, 4, 4, 4, 4],
            [[0, 0, 1], [2, 0, 2], [6, 2, 3], [12, 3, 4], [20, 0, 5]],
            [[0, 2, 1], [2, 6, 2], [6, 11, 3], [12, 11, 4], [20, 5, 5]],
            5,
        ),
        (
            [0, 0, 0, 0, 2, 5, 5, 5, 5],
            [
                [0, 0, 7],
                [21, 7, 14],
                [42, 0, 21],
                [63, 21, 28],
                [84, 0, 35],
            ],
            [
                [0, 7, 7],
                [21, 35, 14],
                [42, 42, 21],
                [63, 77, 28],
                [84, 35, 35],
            ],
            3,
        ),
    ];
    for (knots, lower_controls, upper_controls, height) in cases {
        let lower = extract_bezier(3, &knots, &lower_controls).unwrap();
        let upper = extract_bezier(3, &knots, &upper_controls).unwrap();
        let blocks = extrude_strips(&lower, &upper, height).unwrap();
        assert_eq!(blocks.len(), 2);
        assert!(blocks
            .iter()
            .all(|block| determinant_bound(block).unwrap() > 0.0));
    }
}
