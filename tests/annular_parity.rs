use solidcraft_trivariate::{
    annular::build_annular, annular_audit::audit_annular, RationalBezierMap,
};
use std::f64::consts::FRAC_1_SQRT_2 as Q;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} differs from {expected:.17e} by more than {tolerance:.3e}"
    );
}

fn independent_quarter_annulus(u: f64, v: f64, z0: f64, height: f64, w: f64) -> [f64; 3] {
    // Closed-form rational arc, independent of the library's Bernstein sum.
    let a = (1.0 - u).powi(2) + 2.0 * Q * u * (1.0 - u);
    let b = 2.0 * Q * u * (1.0 - u) + u.powi(2);
    let denominator = (1.0 - u).powi(2) + 2.0 * Q * u * (1.0 - u) + u.powi(2);
    let radius = 2.0 - v;
    [
        radius * a / denominator,
        radius * b / denominator,
        z0 + height * w,
    ]
}

#[test]
fn annular_builder_matches_independent_closed_form_on_closed_cube() {
    let pair = build_annular([1, 2]).expect("valid fixed family heights");
    assert_eq!(pair.heights_mm, [1, 2]);
    assert_eq!(pair.interface, ["lower.w1", "upper.w0", "identity_uv"]);
    assert_eq!(pair.blocks.len(), 2);
    let grid = [0.0, 0.25, 0.5, 0.75, 1.0];

    for (index, block) in pair.blocks.iter().enumerate() {
        assert_eq!(block.degrees, [2, 1, 1]);
        assert_eq!(block.id, if index == 0 { "lower" } else { "upper" });
        let z0 = if index == 0 { 0.0 } else { 1.0 };
        let height = f64::from(pair.heights_mm[index]);
        let map = RationalBezierMap::new(block.degrees, block.controls.clone()).unwrap();
        for u in grid {
            for v in grid {
                for w in grid {
                    let actual = map.evaluate([u, v, w]).unwrap().position;
                    let expected = independent_quarter_annulus(u, v, z0, height, w);
                    for axis in 0..3 {
                        close(actual[axis], expected[axis], 1e-12);
                    }
                    close(actual[0].hypot(actual[1]), 2.0 - v, 1e-12);
                }
            }
        }
    }
}

#[test]
fn annular_certificate_reports_positive_weight_and_bounded_errors() {
    let pair = build_annular([1, 2]).unwrap();
    let certificate = audit_annular(&pair).expect("canonical coefficient family certifies");
    assert!(certificate.min_weight > 0.0);
    assert!(certificate.determinant_lower_bound_mm3 > 0.0);
    assert!(certificate.radial_error_upper_mm <= 1e-12);
    assert!(certificate.height_error_upper_mm <= 1e-12);
}

#[test]
fn annular_audit_rejects_one_ulp_coefficient_mutation() {
    let mut pair = build_annular([1, 2]).unwrap();
    let coefficient = &mut pair.blocks[0].controls[0][0];
    *coefficient = f64::from_bits(coefficient.to_bits() + 1);
    assert!(audit_annular(&pair).is_err());
}

#[test]
fn annular_builder_enforces_integer_height_limits() {
    assert!(build_annular([1, 1000]).is_ok());
    let maximum = build_annular([1000, 1000]).unwrap();
    assert_eq!(maximum.heights_mm, [1000, 1000]);
    assert!(audit_annular(&maximum).is_ok());
    assert!(build_annular([0, 1]).is_err());
    assert!(build_annular([1, 1001]).is_err());
    assert!(build_annular([u32::MAX, 1]).is_err());
}

#[test]
fn annular_binary64_coefficients_match_frozen_python_oracle() {
    let pair = build_annular([1, 2]).unwrap();
    let expected: Vec<[u64; 4]> = include_str!("references/annular_pair_bits.txt")
        .lines()
        .map(|line| {
            let words: Vec<u64> = line
                .split_whitespace()
                .map(|token| u64::from_str_radix(token, 16).expect("hex binary64 word"))
                .collect();
            words.try_into().expect("four coefficient words")
        })
        .collect();
    assert_eq!(expected.len(), 24);
    for (index, control) in pair.blocks.iter().flat_map(|block| &block.controls).enumerate() {
        let actual = control.map(f64::to_bits);
        assert_eq!(actual, expected[index], "control {index}");
    }
}

#[test]
fn annular_upper_block_supports_axial_level_1546_mm() {
    let pair = build_annular([1000, 546]).unwrap();
    let certificate = audit_annular(&pair).unwrap();
    assert!(certificate.height_error_upper_mm <= 1e-12);
    let upper = &pair.blocks[1];
    let end_control = upper.controls[11];
    assert_eq!(end_control[2], 1546.0);
    let map = RationalBezierMap::new(upper.degrees, upper.controls.clone()).unwrap();
    let point = map.evaluate([0.0, 0.0, 1.0]).unwrap().position;
    close(point[2], 1546.0, 1e-12);
}
