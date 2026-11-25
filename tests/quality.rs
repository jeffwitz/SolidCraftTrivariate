use solidcraft_trivariate::bracket::BracketBlock;
use solidcraft_trivariate::exact_bore::build_exact_bore;
use solidcraft_trivariate::exact_bore_audit::audit_exact_bore;
use solidcraft_trivariate::quality::{at_jacobian, bound_block};
use solidcraft_trivariate::RationalBezierMap;

#[test]
fn singular_values_distinguish_column_scaling_from_shear_and_orientation() {
    let anisotropic = [[1e6, 0., 0.], [0., 1., 0.], [0., 0., 1e-6]];
    let q = at_jacobian(anisotropic).unwrap();
    assert_eq!(q.singular_values, [1e6, 1., 1e-6]);
    assert_eq!(q.condition_number, 1e12);
    assert_eq!(q.scaled_determinant, 1.);
    let a = [[1., 0., 0.], [2., 1., 0.], [0., 0., 1.]];
    let q = at_jacobian(a).unwrap();
    let exact = [1. + 2f64.sqrt(), 1., 2f64.sqrt() - 1.];
    for (observed, expected) in q.singular_values.into_iter().zip(exact) {
        assert!((observed - expected).abs() < 1e-14);
    }
    assert!((q.scaled_determinant - 1. / 5f64.sqrt()).abs() < 1e-15);
    let rotated = a.map(|[x, y, z]| [-1000. * y, 1000. * x, 1000. * z]);
    let metric = at_jacobian(rotated).unwrap();
    assert!((q.condition_number / metric.condition_number - 1.).abs() < 1e-14);
    for (a, b) in q.singular_values.into_iter().zip(metric.singular_values) {
        assert!((b / (1000. * a) - 1.).abs() < 1e-14);
    }
    assert!(at_jacobian([[1., 0., 0.], [-1., 0., 0.], [0., 0., 1.]]).is_err());
}

#[test]
fn certified_bounds_contain_local_stretches_on_nearly_tangent_bore() {
    let bore = build_exact_bore(5999, 0, 4000).unwrap();
    let cert = audit_exact_bore([5999, 0], 4000, &bore.blocks).unwrap();
    let mut worst: f64 = 0.;
    for block in &bore.blocks {
        let bounds = bound_block(block, cert.determinant_lower_bound_mm3).unwrap();
        let map = RationalBezierMap::new(
            block.degrees,
            block.controls.iter().map(|p| p.map(|x| x as f64)).collect(),
        )
        .unwrap();
        for i in 0..=8 {
            for j in 0..=4 {
                let p = map.evaluate([i as f64 / 8., j as f64 / 4., 0.375]).unwrap();
                let q = at_jacobian(p.jacobian_columns).unwrap();
                assert!(q.singular_values[2] >= bounds.minimum_stretch_lower);
                assert!(q.singular_values[0] <= bounds.maximum_stretch_upper);
                assert!(q.condition_number <= bounds.condition_upper);
                assert!(q.scaled_determinant + 1e-14 >= bounds.scaled_determinant_lower);
                worst = worst.max(q.condition_number);
            }
        }
    }
    assert!(worst >= 7900.);
}

#[test]
fn certificate_is_translation_invariant_and_scales_in_physical_units() {
    let mut block = BracketBlock {
        id: "affine".into(),
        degrees: [1, 1, 1],
        controls: vec![],
    };
    for x in [0, 1] {
        for y in [0, 2] {
            for z in [0, 3] {
                block.controls.push([x, y, z, 1]);
            }
        }
    }
    let a = bound_block(&block, 6.).unwrap();
    for p in &mut block.controls {
        for c in &mut p[..3] {
            *c = 1000 * (*c + 42);
        }
    }
    let b = bound_block(&block, 6e9).unwrap();
    assert!((b.maximum_stretch_upper / a.maximum_stretch_upper / 1000. - 1.).abs() < 1e-14);
    assert!((b.minimum_stretch_lower / a.minimum_stretch_lower / 1000. - 1.).abs() < 1e-14);
    assert!((a.condition_upper / b.condition_upper - 1.).abs() < 1e-14);
    assert!((a.scaled_determinant_lower / b.scaled_determinant_lower - 1.).abs() < 1e-14);
}
