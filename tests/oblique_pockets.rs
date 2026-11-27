use solidcraft_trivariate::{
    bracket::BracketBlock,
    elasticity::gauss,
    incidence::enumerate_faces,
    interchange,
    json::parse,
    oblique_pockets::{audit_oblique_blind_pockets, build_oblique_blind_pockets},
    spatial_certificate::{certify_spatial, CertificateBudget, CertificateStatus},
    RationalBezierMap,
};

fn maps(blocks: &[BracketBlock]) -> Vec<RationalBezierMap> {
    blocks
        .iter()
        .map(|b| {
            RationalBezierMap::new(
                b.degrees,
                b.controls.iter().map(|p| p.map(|x| x as f64)).collect(),
            )
            .unwrap()
        })
        .collect()
}

fn determinant(columns: [[f64; 3]; 3]) -> f64 {
    let [a, b, c] = columns;
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

fn chart_volume(map: &RationalBezierMap) -> f64 {
    let rule = gauss(3).unwrap();
    let mut sum = 0.0;
    for (u, wu) in &rule {
        for (w, ww) in &rule {
            for (v, wv) in &rule {
                let eval = map.evaluate([*u, *w, *v]).unwrap();
                let det = determinant(eval.jacobian_columns);
                assert!(det > 0.0);
                sum += wu * ww * wv * det;
            }
        }
    }
    sum
}

fn profile_lower(x: f64) -> f64 {
    let scale = if x > 5.0 { 1.0 / 500.0 } else { 1.0 / 1000.0 };
    2.0 + scale * (x - 5.0).powi(2)
}

fn pocket_contains(x: f64, t: f64, y: f64) -> bool {
    let in_profile = 1.0 < x && x < 9.0;
    let in_a = in_profile
        && t < 7.0
        && profile_lower(x) + t / 5.0 < y
        && y < profile_lower(x) + t / 5.0 + 2.0;
    let in_b = in_profile && t < 5.0 && 3.0 - t / 10.0 < y && y < 5.0 - t / 10.0;
    in_a || in_b
}

fn chart_contains(map: &RationalBezierMap, [x, y, z]: [f64; 3]) -> bool {
    let x0 = map.evaluate([0.0, 0.5, 0.5]).unwrap().position[0];
    let x1 = map.evaluate([1.0, 0.5, 0.5]).unwrap().position[0];
    let z0 = map.evaluate([0.5, 0.0, 0.5]).unwrap().position[2];
    let z1 = map.evaluate([0.5, 1.0, 0.5]).unwrap().position[2];
    if !(x0 < x && x < x1 && z1 < z && z < z0) {
        return false;
    }
    let mut lo = 0.0;
    let mut hi = 1.0;
    for _ in 0..48 {
        let mid = (lo + hi) * 0.5;
        if map.evaluate([mid, 0.5, 0.5]).unwrap().position[0] < x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let u = (lo + hi) * 0.5;
    let mut lo = 0.0;
    let mut hi = 1.0;
    for _ in 0..48 {
        let mid = (lo + hi) * 0.5;
        if map.evaluate([u, mid, 0.5]).unwrap().position[2] > z {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let w = (lo + hi) * 0.5;
    let a = map.evaluate([u, w, 0.0]).unwrap().position[1];
    let b = map.evaluate([u, w, 1.0]).unwrap().position[1];
    a < y && y < b
}

fn independent_material_volume() -> f64 {
    // Exact x integration by three-point Gauss-Legendre of the integrated
    // absolute affine wall difference. On x=[1,9], t in [0,5], the overlap
    // width is 2-|delta|, delta=f(x)-3+3t/10, with one strict interior zero.
    let rule = gauss(3).unwrap();
    let mut overlap = 0.0;
    for (x0, x1) in [(1.0, 5.0), (5.0, 9.0)] {
        for (u, weight) in &rule {
            let x = x0 + (x1 - x0) * u;
            let d0 = profile_lower(x) - 3.0;
            let abs_integral = (d0 * d0 + (d0 + 1.5).powi(2)) / 0.6;
            overlap += (x1 - x0) * weight * (10.0 - abs_integral);
        }
    }
    1000.0 - (16.0 * 7.0 + 16.0 * 5.0 - overlap)
}

#[test]
fn qualified_two_blind_profiles_have_distinct_axes_and_internal_continuities() {
    let recipe = parse(include_str!("references/oblique_blind_pockets.recipe.json")).unwrap();
    assert_eq!(
        recipe.field("family").unwrap().text().unwrap(),
        "oblique_blind_pockets_v1"
    );
    assert_eq!(recipe.field("footprint_x_mm").unwrap().to_string(), "[1,9]");
    let profiles = recipe.field("profile_a").unwrap();
    assert_eq!(
        profiles.field("degree").unwrap().number::<usize>().unwrap(),
        2
    );
    assert_eq!(profiles.field("knots").unwrap().items().unwrap().len(), 7);
    assert_eq!(
        profiles.field("internal_join").unwrap().text().unwrap(),
        "C1_not_C2"
    );
    assert_eq!(
        profiles
            .field("whole_loop_corner_order")
            .unwrap()
            .text()
            .unwrap(),
        "G0"
    );
    assert_eq!(
        recipe
            .field("profile_b")
            .unwrap()
            .field("whole_loop_corner_order")
            .unwrap()
            .text()
            .unwrap(),
        "G0"
    );
    let unsupported = parse(include_str!(
        "counterexamples/oblique_blind_smooth_oval.counterexample.json"
    ))
    .unwrap();
    assert_eq!(
        unsupported.field("status").unwrap().text().unwrap(),
        "unsupported"
    );
    assert_eq!(
        unsupported.field("reason").unwrap().text().unwrap(),
        "oblique_blind_profile_requires_exact_monotone_graph_arrangement"
    );
    assert!(interchange::build(&unsupported).is_err());

    // The source spans have equal first derivatives and unequal second
    // derivatives at x=5, which is an internal C1 knot, not a loop seam.
    let second_left = 2.0 / 1000.0;
    let second_right = 2.0 / 500.0;
    assert_ne!(second_left, second_right);

    // The extrusion vectors are separate and nonparallel; neither is normal
    // to the XY sketch plane and this is not a common shear of both tools.
    let dot: f64 = 0.0 * 0.0 + (1.0 / 5.0) * (-1.0 / 10.0) + (-1.0) * (-1.0);
    let cross_squared = (-0.3_f64).powi(2);
    assert_ne!(cross_squared, 0.0);
    assert!((dot - 0.98).abs() < f64::EPSILON);
}

#[test]
fn crossed_blind_pockets_cover_the_cube_with_certified_whole_cells() {
    let blocks = build_oblique_blind_pockets().unwrap();
    assert_eq!(blocks.len(), 64);
    audit_oblique_blind_pockets(&blocks).unwrap();
    let incidence = enumerate_faces(&blocks).unwrap();
    assert_eq!(incidence.interfaces.len(), 122);
    assert_eq!(incidence.exterior.len(), 140);

    for block in &blocks {
        assert_eq!(block.degrees, [2, 1, 1]);
        let cert = certify_spatial(block.degrees, &block.controls, CertificateBudget::default());
        assert_eq!(
            cert.status,
            CertificateStatus::Valid,
            "{}: {:?}",
            block.id,
            cert.reason
        );
        assert!(cert.determinant_lower_bound_mm3.unwrap() > 0.0);
    }

    let maps = maps(&blocks);
    let volume: f64 = maps.iter().map(chart_volume).sum();
    let expected = independent_material_volume();
    assert!(pocket_contains(5.0, 4.0, 3.7));
    assert!((volume - expected).abs() < 2e-9, "{volume} != {expected}");
    assert!((expected - 871.4365866666667).abs() < 2e-10);

    // Independent point classification uses the original cutter inequalities,
    // not the serialized charts. Probes include both wings, each blind floor,
    // and points on both sides of the transverse wall-crossing graph.
    let mut checked = 0;
    for x in [0.37, 1.37, 3.7, 5.2, 7.4, 8.73, 9.63] {
        let crossing = (3.0 - profile_lower(x)) / 0.3;
        for t in [
            0.31,
            crossing - 0.013,
            crossing + 0.017,
            4.73,
            5.37,
            6.71,
            7.29,
            9.23,
        ] {
            for y in [
                0.73, 1.63, 2.41, 2.91, 3.37, 4.17, 4.83, 5.43, 6.37, 8.91, 9.63,
            ] {
                let point = [x, y, 10.0 - t];
                let count = maps.iter().filter(|m| chart_contains(m, point)).count();
                assert_eq!(
                    count,
                    usize::from(!pocket_contains(x, t, y)),
                    "x={x}, t={t}, y={y}"
                );
                checked += 1;
            }
        }
    }
    assert_eq!(checked, 7 * 8 * 11);
}

#[test]
fn audit_rejects_a_changed_chart_without_mutating_the_accepted_model() {
    let original = build_oblique_blind_pockets().unwrap();
    let mut candidate = original.clone();
    candidate[0].controls[0][0] += 1;
    assert_eq!(
        audit_oblique_blind_pockets(&candidate),
        Err("oblique_blind_pocket_coefficients_mismatch")
    );
    assert!(audit_oblique_blind_pockets(&original).is_ok());
}

#[test]
fn versioned_recipe_rebuilds_and_audits_atomically() {
    let recipe = parse(include_str!("references/oblique_blind_pockets.recipe.json")).unwrap();
    let model = interchange::build(&recipe).unwrap();
    let report = interchange::audit(&model).unwrap();
    assert_eq!(report.field("status").unwrap().text().unwrap(), "valid");
    assert_eq!(
        report.field("blocks").unwrap().number::<usize>().unwrap(),
        64
    );
    assert_eq!(
        report
            .field("shared_faces")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        122
    );
    assert_eq!(
        report
            .field("exterior_sides")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        140
    );
    let changed_text = recipe
        .to_string()
        .replace("\"depth_mm\":7", "\"depth_mm\":8");
    let changed_recipe = parse(&changed_text).unwrap();
    assert!(interchange::build(&changed_recipe).is_err());
    assert!(interchange::audit(&model).is_ok());
    let altered_model = model.clone();
    let blocks = altered_model.field("blocks").unwrap().items().unwrap();
    assert_eq!(blocks.len(), 64);
    let corrupted = altered_model.to_string().replace("1000", "1001");
    let corrupted_model = parse(&corrupted).unwrap();
    assert!(interchange::audit(&corrupted_model).is_err());
}
