use solidcraft_trivariate::{MapError, RationalBezierMap};
use std::f64::consts::FRAC_1_SQRT_2 as Q;

fn close(actual: [f64; 3], expected: [f64; 3]) {
    for (a, e) in actual.into_iter().zip(expected) {
        assert!(
            (a - e).abs() <= 1e-12,
            "{actual:?} differs from {expected:?}"
        );
    }
}

#[test]
fn annular_positions_and_derivatives_on_closed_cube() {
    let mut controls = Vec::new();
    for (x, y, weight) in [(1., 0., 1.), (1., 1., Q), (0., 1., 1.)] {
        for radius in [2., 1.] {
            for height in [0., 3.] {
                controls.push([
                    weight * radius * x,
                    weight * radius * y,
                    weight * height,
                    weight,
                ]);
            }
        }
    }
    let map = RationalBezierMap::new([2, 1, 1], controls).unwrap();
    let grid = [0., 0.25, 0.5, 0.75, 1.];
    for u in grid {
        // The direction of the rational circular arc is twice this half-angle.
        let (a, b) = (1. - u + Q * u, Q * u);
        let (sin, cos) = (2. * b.atan2(a)).sin_cos();
        let speed = 2. * Q / (a * a + b * b);
        for v in grid {
            for w in grid {
                let r = 2. - v;
                let value = map.evaluate([u, v, w]).unwrap();
                close(value.position, [r * cos, r * sin, 3. * w]);
                close(
                    value.jacobian_columns[0],
                    [-r * sin * speed, r * cos * speed, 0.],
                );
                close(value.jacobian_columns[1], [-cos, -sin, 0.]);
                close(value.jacobian_columns[2], [0., 0., 3.]);
            }
        }
    }
}

#[test]
fn affine_and_rejections() {
    let controls = (0..8)
        .map(|n| [(n / 4) as f64, ((n / 2) % 2) as f64, (n % 2) as f64, 1.])
        .collect();
    let map = RationalBezierMap::new([1; 3], controls).unwrap();
    let value = map.evaluate([0.2, 0.4, 0.7]).unwrap();
    close(value.position, [0.2, 0.4, 0.7]);
    for (axis, column) in value.jacobian_columns.into_iter().enumerate() {
        let mut expected = [0.; 3];
        expected[axis] = 1.;
        close(column, expected);
    }
    assert_eq!(
        map.evaluate([f64::NAN, 0., 0.]).unwrap_err(),
        MapError::OutsideDomain
    );
    assert_eq!(
        RationalBezierMap::new([1; 3], vec![]).unwrap_err(),
        MapError::WrongControlCount
    );
    assert_eq!(
        RationalBezierMap::new([1; 3], vec![[0.; 4]; 8]).unwrap_err(),
        MapError::InvalidControl
    );
    assert_eq!(
        RationalBezierMap::new([4, 1, 1], vec![[1.; 4]; 20]).unwrap_err(),
        MapError::UnsupportedDegrees
    );
}

#[test]
fn biquadratic_rational_chart_matches_separable_quotients_and_product_rule() {
    let mut controls = Vec::new();
    for (i, wu) in [1., 2., 1.].into_iter().enumerate() {
        for (j, wv) in [1., 3., 1.].into_iter().enumerate() {
            for k in 0..2 {
                let (x, y, weight) = (i as f64 / 2., j as f64 / 2., wu * wv);
                controls.push([
                    x * weight,
                    y * weight,
                    (k as f64 + 0.25 * x * y) * weight,
                    weight,
                ]);
            }
        }
    }
    let map = RationalBezierMap::new([2, 2, 1], controls).unwrap();
    fn coordinate(t: f64, weight: f64) -> (f64, f64) {
        let n = t * (weight + (1. - weight) * t);
        let d = 1. + 2. * (weight - 1.) * t * (1. - t);
        let dn = weight + 2. * (1. - weight) * t;
        let dd = 2. * (weight - 1.) * (1. - 2. * t);
        (n / d, (dn * d - n * dd) / (d * d))
    }
    for u in [0., 0.125, 0.5, 1.] {
        for v in [0., 0.375, 0.75, 1.] {
            let (x, dx) = coordinate(u, 2.);
            let (y, dy) = coordinate(v, 3.);
            let eval = map.evaluate([u, v, 0.625]).unwrap();
            close(eval.position, [x, y, 0.625 + 0.25 * x * y]);
            close(eval.jacobian_columns[0], [dx, 0., 0.25 * dx * y]);
            close(eval.jacobian_columns[1], [0., dy, 0.25 * x * dy]);
            close(eval.jacobian_columns[2], [0., 0., 1.]);
        }
    }
}

#[test]
fn every_tensor_degree_combination_has_the_correct_axis_strides() {
    for p in 1..=3 {
        for q in 1..=3 {
            for r in 1..=3 {
                let degrees = [p, q, r];
                let mut controls = Vec::new();
                for i in 0..=p {
                    for j in 0..=q {
                        for k in 0..=r {
                            let coordinates = [i, j, k];
                            let fun = |a: usize| f64::from(coordinates[a] == degrees[a]);
                            controls.push([fun(0), fun(1), fun(2), 1.]);
                        }
                    }
                }
                let map = RationalBezierMap::new(degrees, controls).unwrap();
                let uvw: [f64; 3] = [0.125, 0.375, 0.875];
                let e = map.evaluate(uvw).unwrap();
                for a in 0..3 {
                    assert!((e.position[a] - uvw[a].powi(degrees[a] as i32)).abs() < 1e-14);
                    for b in 0..3 {
                        let expected = if a == b {
                            degrees[a] as f64 * uvw[a].powi(degrees[a] as i32 - 1)
                        } else {
                            0.
                        };
                        assert!((e.jacobian_columns[b][a] - expected).abs() < 1e-14);
                    }
                }
            }
        }
    }
}
