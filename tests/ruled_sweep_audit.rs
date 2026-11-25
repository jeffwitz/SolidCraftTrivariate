use solidcraft_trivariate::bracket::BracketBlock;
use solidcraft_trivariate::ruled_sweep::{build_ruled_sweep, SweepRecipe, SweepSection};
use solidcraft_trivariate::ruled_sweep_audit::audit_ruled_sweep;
use solidcraft_trivariate::section_audit::audit_section;
use solidcraft_trivariate::sections::section_blocks;
use solidcraft_trivariate::RationalBezierMap;

const CORPUS: [(i32, i32, i32, i32, i32); 6] = [
    (8, 0, 8, 0, 0),
    (0, 8, 12, 4, -8),
    (-8, 1, 7, -16, 24),
    (-8, -8, 16, 8, 4),
    (16, 0, 1, -32, 32),
    (4, -12, 32, 32, -32),
];

const GAUSS8_X: [f64; 8] = [
    0.0198550717512319,
    0.101666761293187,
    0.237233795041836,
    0.408282678752175,
    0.591717321247825,
    0.762766204958164,
    0.898333238706813,
    0.980144928248768,
];
const GAUSS8_W: [f64; 8] = [
    0.0506142681451881,
    0.111190517226687,
    0.156853322938944,
    0.181341891689181,
    0.181341891689181,
    0.156853322938944,
    0.111190517226687,
    0.0506142681451881,
];

fn recipe(section: SweepSection, case: (i32, i32, i32, i32, i32)) -> SweepRecipe {
    let (a, b, height, tx, ty) = case;
    SweepRecipe {
        section,
        top_multiplier_eighths: [a, b],
        top_translation_eighths_mm: [tx, ty],
        height_mm: height,
    }
}

fn as_map(block: &BracketBlock) -> RationalBezierMap {
    RationalBezierMap::new(
        block.degrees,
        block
            .controls
            .iter()
            .map(|point| point.map(|value| value as f64))
            .collect(),
    )
    .unwrap()
}

fn integrated_volume_q8(blocks: &[BracketBlock]) -> f64 {
    blocks
        .iter()
        .map(|block| {
            let map = as_map(block);
            let mut volume = 0.0;
            for (iu, u) in GAUSS8_X.iter().enumerate() {
                for (iv, v) in GAUSS8_X.iter().enumerate() {
                    for (iw, w) in GAUSS8_X.iter().enumerate() {
                        let jac = map.evaluate([*u, *v, *w]).unwrap().jacobian_columns;
                        let cross = [
                            jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1],
                            jac[1][2] * jac[2][0] - jac[1][0] * jac[2][2],
                            jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0],
                        ];
                        let determinant =
                            jac[0][0] * cross[0] + jac[0][1] * cross[1] + jac[0][2] * cross[2];
                        assert!(determinant > 0.0);
                        volume += GAUSS8_W[iu] * GAUSS8_W[iv] * GAUSS8_W[iw] * determinant;
                    }
                }
            }
            volume
        })
        .sum()
}

#[test]
fn corpus_matches_independent_similarity_and_volume_formulas() {
    for section in [SweepSection::Rectangle, SweepSection::PiercedSquare] {
        let family = section.family();
        let source = section_blocks(family).unwrap();
        let source_certificate = audit_section(family, &source).unwrap();
        let base_area = match section {
            SweepSection::Rectangle => 400.0,
            SweepSection::PiercedSquare => 400.0 - 16.0 * std::f64::consts::PI,
        };
        for case in CORPUS {
            let recipe = recipe(section, case);
            let blocks = build_ruled_sweep(recipe).unwrap();
            let audit = audit_ruled_sweep(recipe, &blocks).unwrap();
            assert_eq!(
                audit.certificate.design_faces,
                source_certificate.design_faces
            );
            assert_eq!(audit.certificate.incidence, source_certificate.incidence);
            let (a, b, height, _, _) = case;
            let expected_volume = base_area
                * f64::from(height)
                * (1.0 + f64::from(a) / 8.0 + f64::from(a * a + b * b) / 64.0)
                / 3.0;
            assert!((audit.ideal_volume_mm3 - expected_volume).abs() < 1e-10);
            let quadrature_volume = integrated_volume_q8(&blocks);
            assert!((quadrature_volume - expected_volume).abs() / expected_volume < 1e-8);

            for (source_block, sweep_block) in source.iter().zip(&blocks) {
                assert_eq!(source_block.id, sweep_block.id);
                let (src, dst) = (as_map(source_block), as_map(sweep_block));
                for sample in 0..100 {
                    let uvw = [
                        f64::from(sample % 11) / 10.0,
                        f64::from(sample % 7) / 6.0,
                        f64::from(sample % 17) / 16.0,
                    ];
                    let source_sample = src.evaluate(uvw).unwrap();
                    let target_sample = dst.evaluate(uvw).unwrap();
                    let point = source_sample.position;
                    let actual = target_sample.position;
                    let w = uvw[2];
                    let multiplier_re = 1.0 + w * (f64::from(a) / 8.0 - 1.0);
                    let multiplier_im = w * f64::from(b) / 8.0;
                    let expected = [
                        multiplier_re * point[0] - multiplier_im * point[1]
                            + w * f64::from(case.3) / 8.0,
                        multiplier_im * point[0]
                            + multiplier_re * point[1]
                            + w * f64::from(case.4) / 8.0,
                        w * f64::from(height),
                    ];
                    for axis in 0..3 {
                        assert!((actual[axis] - expected[axis]).abs() < 2e-12);
                    }
                    let complex = |x: f64, y: f64| {
                        [
                            multiplier_re * x - multiplier_im * y,
                            multiplier_im * x + multiplier_re * y,
                        ]
                    };
                    for axis in 0..2 {
                        let expected_column = complex(
                            source_sample.jacobian_columns[axis][0],
                            source_sample.jacobian_columns[axis][1],
                        );
                        assert!(
                            (target_sample.jacobian_columns[axis][0] - expected_column[0]).abs()
                                < 3e-11
                        );
                        assert!(
                            (target_sample.jacobian_columns[axis][1] - expected_column[1]).abs()
                                < 3e-11
                        );
                        assert!(target_sample.jacobian_columns[axis][2].abs() < 3e-12);
                    }
                    let derivative_re = f64::from(a) / 8.0 - 1.0;
                    let derivative_im = f64::from(b) / 8.0;
                    let expected_w = [
                        derivative_re * point[0] - derivative_im * point[1]
                            + f64::from(case.3) / 8.0,
                        derivative_im * point[0]
                            + derivative_re * point[1]
                            + f64::from(case.4) / 8.0,
                        f64::from(height),
                    ];
                    for (actual, expected) in
                        target_sample.jacobian_columns[2].iter().zip(expected_w)
                    {
                        assert!((actual - expected).abs() < 3e-11);
                    }
                    let mut inverse_point = actual;
                    if w == 1.0 {
                        inverse_point[2] = f64::from(height);
                    }
                    let inverse = recipe.inverse_coordinates(inverse_point).unwrap();
                    assert!((inverse[0] - point[0]).abs() < 3e-11);
                    assert!((inverse[1] - point[1]).abs() < 3e-11);
                    assert!((inverse[2] - w).abs() < 3e-14);
                }
            }
        }
    }
}

#[test]
fn audit_checks_exact_controls_incidence_selections_and_positive_block_scales() {
    let recipe = recipe(SweepSection::PiercedSquare, CORPUS[2]);
    let blocks = build_ruled_sweep(recipe).unwrap();
    let certificate = audit_ruled_sweep(recipe, &blocks).unwrap();
    assert_eq!(certificate.minimum_squared_scale, [1, 257]);
    assert!(certificate.certificate.determinant_lower_bound_mm3 > 0.0);
    assert_eq!(
        certificate.certificate.design_faces.len(),
        certificate.certificate.incidence.exterior.len()
    );

    let mut permuted = blocks.clone();
    permuted.reverse();
    assert!(audit_ruled_sweep(recipe, &permuted).is_ok());

    let mut commonly_scaled = blocks.clone();
    for point in &mut commonly_scaled[0].controls {
        for coefficient in point {
            *coefficient *= 7;
        }
    }
    assert!(audit_ruled_sweep(recipe, &commonly_scaled).is_ok());

    let mut coefficient_change = blocks.clone();
    coefficient_change[0].controls[0][0] += 1;
    assert!(audit_ruled_sweep(recipe, &coefficient_change).is_err());

    let mut weight_change = blocks.clone();
    weight_change[0].controls[0][3] += 1;
    assert!(audit_ruled_sweep(recipe, &weight_change).is_err());

    let mut degree_change = blocks.clone();
    degree_change[0].degrees[0] = 1;
    assert!(audit_ruled_sweep(recipe, &degree_change).is_err());

    let mut negative_scale = blocks.clone();
    for point in &mut negative_scale[0].controls {
        for coefficient in point {
            *coefficient = -*coefficient;
        }
    }
    assert!(audit_ruled_sweep(recipe, &negative_scale).is_err());
    assert!(audit_ruled_sweep(recipe, &blocks[1..]).is_err());
}

#[test]
fn collapsed_and_incompatible_sweeps_are_rejected() {
    let collapsed = SweepRecipe {
        section: SweepSection::Rectangle,
        top_multiplier_eighths: [-16, 0],
        top_translation_eighths_mm: [0, 0],
        height_mm: 8,
    };
    let blocks = build_ruled_sweep(collapsed);
    assert!(blocks.is_err());
    assert!(audit_ruled_sweep(collapsed, &[]).is_err());
    for [a, b] in [[-16, 0], [0, 0]] {
        let pinched = SweepRecipe {
            top_multiplier_eighths: [a, b],
            ..collapsed
        };
        assert!(build_ruled_sweep(pinched).is_err());
        assert!(audit_ruled_sweep(pinched, &[]).is_err());
    }

    for invalid in [
        SweepRecipe {
            height_mm: 0,
            ..collapsed
        },
        SweepRecipe {
            height_mm: 33,
            ..collapsed
        },
        SweepRecipe {
            top_multiplier_eighths: [17, 0],
            ..collapsed
        },
        SweepRecipe {
            top_translation_eighths_mm: [33, 0],
            ..collapsed
        },
    ] {
        assert!(audit_ruled_sweep(invalid, &[]).is_err());
    }

    let valid = recipe(SweepSection::Rectangle, CORPUS[0]);
    let blocks = build_ruled_sweep(valid).unwrap();
    let incompatible = recipe(SweepSection::Rectangle, CORPUS[1]);
    assert!(audit_ruled_sweep(incompatible, &blocks).is_err());
}
