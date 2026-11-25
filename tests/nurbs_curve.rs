use solidcraft_trivariate::nurbs_curve::{extract_bezier, BezierSpan};

fn de_boor(degree: usize, knots: &[i64], controls: &[[i64; 3]], u: f64) -> [f64; 2] {
    let n = controls.len() - 1;
    let span = if u == knots[n + 1] as f64 {
        n
    } else {
        (degree..=n)
            .find(|&i| knots[i] as f64 <= u && u < knots[i + 1] as f64)
            .unwrap()
    };
    let mut work = (0..=degree)
        .map(|j| controls[span - degree + j].map(|x| x as f64))
        .collect::<Vec<_>>();
    for level in 1..=degree {
        for j in (level..=degree).rev() {
            let index = span - degree + j;
            let denominator = (knots[index + degree + 1 - level] - knots[index]) as f64;
            let alpha = (u - knots[index] as f64) / denominator;
            work[j] = std::array::from_fn(|k| (1. - alpha) * work[j - 1][k] + alpha * work[j][k]);
        }
    }
    [
        work[degree][0] / work[degree][2],
        work[degree][1] / work[degree][2],
    ]
}

fn bezier(span: &BezierSpan, degree: usize, u: f64) -> [f64; 2] {
    let t = (u - span.interval[0] as f64) / (span.interval[1] - span.interval[0]) as f64;
    let mut h = [0.; 3];
    for (i, p) in span.controls.iter().enumerate() {
        let binomial = match (degree, i) {
            (1, _) => 1.,
            (2, 1) => 2.,
            (3, 1 | 2) => 3.,
            _ => 1.,
        };
        let basis = binomial * t.powi(i as i32) * (1. - t).powi((degree - i) as i32);
        for k in 0..3 {
            h[k] += basis * p[k] as f64;
        }
    }
    [h[0] / h[2], h[1] / h[2]]
}

#[test]
fn cubic_variable_weight_repeated_knot_matches_independent_de_boor() {
    let degree = 3;
    let knots = [0, 0, 0, 0, 1, 1, 2, 2, 2, 2];
    let controls = [
        [0, 0, 1],
        [2, 4, 2],
        [12, -3, 3],
        [20, 8, 1],
        [24, 0, 4],
        [28, 12, 2],
    ];
    let spans = extract_bezier(degree, &knots, &controls).unwrap();
    assert_eq!(
        spans.iter().map(|s| s.interval).collect::<Vec<_>>(),
        [[0, 1], [1, 2]]
    );
    assert!(spans
        .iter()
        .all(|s| s.controls.len() == degree + 1 && s.controls.iter().all(|p| p[2] > 0)));
    for span in &spans {
        for sample in 0..=20 {
            let t = sample as f64 / 20.;
            let u = span.interval[0] as f64 + t * (span.interval[1] - span.interval[0]) as f64;
            let expected = de_boor(degree, &knots, &controls, u);
            let actual = bezier(span, degree, u);
            assert!((expected[0] - actual[0]).abs() < 1e-12);
            assert!((expected[1] - actual[1]).abs() < 1e-12);
        }
    }
}

#[test]
fn linear_and_quadratic_clamped_curves_extract_without_extra_spans() {
    let linear = extract_bezier(1, &[0, 0, 5, 5], &[[0, 0, 2], [10, 4, 2]]).unwrap();
    assert_eq!(
        linear,
        vec![BezierSpan {
            interval: [0, 5],
            controls: vec![[0, 0, 1], [5, 2, 1]]
        }]
    );
    let quadratic =
        extract_bezier(2, &[0, 0, 0, 3, 3, 3], &[[0, 0, 1], [3, 6, 2], [12, 0, 3]]).unwrap();
    assert_eq!(quadratic.len(), 1);
    assert_eq!(quadratic[0].interval, [0, 3]);
    assert_eq!(quadratic[0].controls.len(), 3);
}

#[test]
fn rejects_invalid_knot_domains_weights_and_homogeneous_overflow() {
    assert!(extract_bezier(0, &[0, 0], &[[0, 0, 1]]).is_err());
    assert!(extract_bezier(4, &[0; 10], &[[0, 0, 1]; 6]).is_err());
    assert!(extract_bezier(2, &[0, 0, 1, 2, 2], &[[0, 0, 1]; 2]).is_err());
    assert!(extract_bezier(1, &[0, 1, 0, 1], &[[0, 0, 1], [1, 1, 1]]).is_err());
    assert!(extract_bezier(1, &[0, 0, 1, 1], &[[0, 0, 1], [1, 1, 0]]).is_err());
    assert!(extract_bezier(2, &[0, 0, 0, 2, 2, 2], &[[0, 0, 1], [1, 1, 1]]).is_err());

    let high = i64::MAX;
    let overflow = extract_bezier(
        3,
        &[0, 0, 0, 0, 1, 2, 2, 2, 2],
        &[[0, 0, 1], [high, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
    );
    assert_eq!(overflow, Err("NURBS coefficient exceeds i64"));
}
