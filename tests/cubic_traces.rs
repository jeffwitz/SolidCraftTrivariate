use solidcraft_trivariate::{
    bracket::BracketBlock, incidence::enumerate_faces, MapError, RationalBezierMap,
};

#[test]
fn cubic_tensor_map_evaluates_position_and_all_derivatives() {
    let mut controls = Vec::new();
    for i in 0..=3 {
        for j in 0..=3 {
            for k in 0..=3 {
                let (u, v, w) = (i as f64 / 3., j as f64 / 3., k as f64 / 3.);
                let u_squared = [0., 0., 1. / 3., 1.][i];
                controls.push([u, v, w + u_squared * v, 1.]);
            }
        }
    }
    let map = RationalBezierMap::new([3, 3, 3], controls).unwrap();
    let [u, v, w] = [0.25_f64, 0.5, 0.75];
    let value = map.evaluate([u, v, w]).unwrap();
    assert!((value.position[0] - u).abs() < 1e-14);
    assert!((value.position[1] - v).abs() < 1e-14);
    assert!((value.position[2] - (w + u * u * v)).abs() < 1e-14);
    assert!((value.jacobian_columns[0][2] - 2. * u * v).abs() < 1e-14);
    assert!((value.jacobian_columns[1][2] - u * u).abs() < 1e-14);
    assert!((value.jacobian_columns[2][2] - 1.).abs() < 1e-14);
}

#[test]
fn cubic_and_quadratic_whole_face_traces_match_after_exact_elevation() {
    // Use identical linear-in-y traces, represented at degrees two and three.
    let linear_controls = |x0: i64, degree: usize| {
        let mut controls = Vec::new();
        for i in 0..=1 {
            for j in 0..=degree {
                for k in 0..=1 {
                    let (scale, y) = if degree == 2 {
                        (1, j as i64)
                    } else {
                        (3, 2 * j as i64)
                    };
                    controls.push([scale * (x0 + i as i64), y, scale * k as i64, scale]);
                }
            }
        }
        controls
    };
    let a = BracketBlock {
        id: "linear".into(),
        degrees: [1, 2, 1],
        controls: linear_controls(0, 2),
    };
    let b = BracketBlock {
        id: "cubic".into(),
        degrees: [1, 3, 1],
        controls: linear_controls(1, 3),
    };
    let result = enumerate_faces(&[a, b]).unwrap();
    let interface = result
        .interfaces
        .iter()
        .find(|f| f.first.block != f.second.block)
        .unwrap();
    assert_eq!(interface.first.block, "cubic");
    assert_eq!(interface.second.block, "linear");
    assert_eq!(result.interfaces.len(), 1);
}

#[test]
fn cubic_degree_is_supported_but_quartic_is_rejected() {
    assert!(RationalBezierMap::new([3, 1, 1], vec![[1., 1., 1., 1.]; 16]).is_ok());
    assert_eq!(
        RationalBezierMap::new([4, 1, 1], vec![[1., 1., 1., 1.]; 20]).unwrap_err(),
        MapError::UnsupportedDegrees
    );
}
