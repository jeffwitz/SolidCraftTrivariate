use solidcraft_trivariate::{
    interchange,
    json::{object, Json},
    rational_sector::{build_rational_sector, RationalSector},
    rational_sector_audit::audit_rational_sector,
};

fn eval(block: &solidcraft_trivariate::bracket::BracketBlock, [u, v, w]: [f64; 3]) -> [f64; 3] {
    let bu = [(1.0 - u).powi(2), 2.0 * u * (1.0 - u), u * u];
    let bv = [1.0 - v, v];
    let bw = [1.0 - w, w];
    let mut h = [0.0; 4];
    for (i, &a) in bu.iter().enumerate() {
        for (j, &b) in bv.iter().enumerate() {
            for (k, &c) in bw.iter().enumerate() {
                let p = block.controls[(i * 2 + j) * 2 + k];
                for q in 0..4 {
                    h[q] += a * b * c * p[q] as f64;
                }
            }
        }
    }
    [h[0] / h[3], h[1] / h[3], h[2] / h[3]]
}

fn check_recipe(recipe: RationalSector) {
    let blocks = build_rational_sector(&recipe).unwrap();
    let certificate = audit_rational_sector(&recipe, &blocks).unwrap();
    assert!(certificate.determinant_lower_bound_mm3 > 0.0);
    let t = recipe.half_angle_tangent[0] as f64 / recipe.half_angle_tangent[1] as f64;
    let [ri, ro] = recipe.radii_mm.map(f64::from);
    let h = f64::from(recipe.height_mm);
    for block in &blocks {
        let [start, end] = match block.id.as_str() {
            "whole" => [0.0, 1.0],
            "first" => [0.0, 0.5],
            "second" => [0.5, 1.0],
            _ => panic!("unexpected chart ID"),
        };
        for local_u in [0.0, 0.17, 0.5, 0.83, 1.0] {
            for v in [0.0, 0.31, 1.0] {
                for w in [0.0, 0.63, 1.0] {
                    let global_u = start + (end - start) * local_u;
                    let theta = 2.0 * (t * global_u).atan();
                    let radius = ro + (ri - ro) * v;
                    let expected = [radius * theta.cos(), radius * theta.sin(), h * w];
                    let actual = eval(block, [local_u, v, w]);
                    for axis in 0..3 {
                        assert!((actual[axis] - expected[axis]).abs() < 2e-12);
                    }
                }
            }
        }
    }
}

#[test]
fn exact_charts_match_the_analytic_polar_domain() {
    for (t, radii, h) in [
        ([1, 4], [1, 64], 1),
        ([70, 29], [3, 6], 8),
        ([169, 70], [4, 5], 64),
        ([4, 1], [1, 2], 32),
    ] {
        check_recipe(RationalSector {
            half_angle_tangent: t,
            radii_mm: radii,
            height_mm: h,
        });
    }
}

#[test]
fn invalid_parameters_and_mutated_control_nets_fail_audit() {
    let recipe = RationalSector {
        half_angle_tangent: [1, 1],
        radii_mm: [2, 6],
        height_mm: 8,
    };
    assert!(build_rational_sector(&RationalSector {
        radii_mm: [6, 2],
        ..recipe
    })
    .is_err());
    let mut blocks = build_rational_sector(&recipe).unwrap();
    blocks[0].controls[5][0] += 1;
    assert!(audit_rational_sector(&recipe, &blocks).is_err());
}

#[test]
fn json_exchange_reports_exact_volume_and_analytic_volume_rate() {
    let recipe = object([
        ("family", Json::String("rational_sector_prism_v1".into())),
        (
            "half_angle_tangent",
            Json::Array(vec![Json::Number("169".into()), Json::Number("70".into())]),
        ),
        (
            "radii_mm",
            Json::Array(vec![Json::Number("3".into()), Json::Number("6".into())]),
        ),
        ("height_mm", Json::Number("8".into())),
    ]);
    let model = interchange::build(&recipe).unwrap();
    let audit = interchange::audit(&model).unwrap();
    assert_eq!(
        audit.field("status").unwrap(),
        &Json::String("valid".into())
    );
    assert_eq!(audit.field("blocks").unwrap(), &Json::Number("2".into()));
    assert_eq!(
        audit
            .field("volume_exact")
            .unwrap()
            .field("atan_coefficient_mm3")
            .unwrap(),
        &Json::Number("216".into())
    );
    let expected = 8.0 * (36.0 - 9.0) * 4900.0 / (169.0 * 169.0 + 70.0 * 70.0);
    let actual: f64 = audit.field("volume_rate_dt_mm3").unwrap().number().unwrap();
    assert!((actual - expected).abs() < 1e-12);
}
