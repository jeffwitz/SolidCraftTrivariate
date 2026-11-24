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
        RationalBezierMap::new([3, 1, 1], vec![[1.; 4]; 16]).unwrap_err(),
        MapError::UnsupportedDegrees
    );
}
