use solidcraft_trivariate::{
    elasticity::gauss,
    g2_cross_plate::boolean::{construct, Operation},
    json::{object, parse, Json},
    RationalBezierMap,
};

fn recipe(rise: i64) -> Json {
    let pair = parse(include_str!("references/g2_cross_sketches.recipe.json")).unwrap();
    object([
        ("family", Json::String("inclined_crossed_tools_v1".into())),
        ("sketch_pair", pair),
        ("transverse_rise_mm", Json::Number(rise.to_string())),
    ])
}

fn maps(input: &Json, op: Operation) -> Vec<RationalBezierMap> {
    construct(input, op)
        .unwrap()
        .into_iter()
        .map(|b| {
            RationalBezierMap::new(
                b.degrees,
                b.controls
                    .into_iter()
                    .map(|p| p.map(|v| v as f64))
                    .collect(),
            )
            .unwrap()
        })
        .collect()
}

// Physical-space inverse specialized to the monotone chart form. Probes are
// deliberately off knots, cut faces and the two collapsed end edges.
fn contains(m: &RationalBezierMap, [x, y, z]: [f64; 3]) -> bool {
    let pos = |u, v, w| m.evaluate([u, v, w]).unwrap().position;
    if y <= pos(0., 0., 0.)[1] || y >= pos(1., 0., 0.)[1] {
        return false;
    }
    let (mut lo, mut hi) = (0., 1.);
    for _ in 0..52 {
        let mid = (lo + hi) / 2.;
        if pos(mid, 0., 0.)[1] < y {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let u = (lo + hi) / 2.;
    let (right, left) = (pos(u, 0., 0.)[0], pos(u, 1., 0.)[0]);
    if x <= left || x >= right {
        return false;
    }
    let v = (right - x) / (right - left);
    z > pos(u, v, 0.)[2] && z < pos(u, v, 1.)[2]
}

fn volume(maps: &[RationalBezierMap]) -> f64 {
    let mut sum = 0.;
    for map in maps {
        // Each map is affine in v and w. Parallel swept caps give a height
        // affine in x, which is integrated by the midpoint rule in v.
        for (u, wu) in gauss(16).unwrap() {
            let [a, b, c] = map.evaluate([u, 0.5, 0.5]).unwrap().jacobian_columns;
            let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]);
            assert!(det > 0.);
            sum += wu * det;
        }
    }
    sum
}

#[test]
fn crossed_tools_match_independent_volumes_and_single_coverage() {
    // SolidCore reference measured before the native Boolean implementation.
    let common = 712.3493494134794;
    let a = 3489.731852993524;
    let b = 4022.2222222222226;
    for rise in [-3, 0, 3] {
        for (op, expected_volume) in [
            (Operation::Intersection, common),
            (Operation::ProfileMinusRectangle, a - common),
            (Operation::RectangleMinusProfile, b - common),
            (Operation::Union, a + b - common),
        ] {
            let maps = maps(&recipe(rise), op);
            assert!((volume(&maps) - expected_volume).abs() < 1e-6);
            // At y=24 the center line is well inside the free curved loop.
            for (x, in_a) in [(3., false), (75., true), (96., false)] {
                for z in [1., 3.9, 4.1, 6., 7.9, 8.1, 11.] {
                    let local_z = z - rise as f64 * x / 100.;
                    let in_b = (4.0..8.0).contains(&local_z);
                    let expected = match op {
                        Operation::Intersection => in_a && in_b,
                        Operation::ProfileMinusRectangle => in_a && !in_b,
                        Operation::RectangleMinusProfile => !in_a && in_b,
                        Operation::Union => in_a || in_b,
                    };
                    assert_eq!(
                        maps.iter().filter(|m| contains(m, [x, 24., z])).count(),
                        usize::from(expected)
                    );
                }
            }
        }
    }
}

#[test]
fn refuses_tilt_at_stock_cap_and_agrees_with_native_sketch() {
    for rise in [-4, 4, i64::MAX] {
        assert!(construct(&recipe(rise), Operation::Union).is_err());
    }
    let raw = parse(include_str!("references/g2_cross_native.recipe.json")).unwrap();
    assert_eq!(
        construct(&recipe(0), Operation::Intersection).unwrap(),
        construct(&raw, Operation::Intersection).unwrap()
    );
}

#[test]
fn cli_publishes_traces_and_refuses_atomically() {
    use std::{fs, process::Command};
    let directory = std::env::temp_dir().join(format!("crossed-tool-cli-{}", std::process::id()));
    fs::create_dir_all(&directory).unwrap();
    let input = directory.join("recipe.json");
    let output = directory.join("charts.json");
    let run = || {
        Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
            .arg("crossed-boolean")
            .arg(&input)
            .arg("union")
            .arg(&output)
            .output()
            .unwrap()
    };
    fs::write(&input, recipe(3).to_string()).unwrap();
    assert!(run().status.success());
    let bytes = fs::read_to_string(&output).unwrap();
    let payload = parse(&bytes).unwrap();
    assert_eq!(
        payload
            .field("profile_bands")
            .unwrap()
            .items()
            .unwrap()
            .len(),
        6
    );
    assert_eq!(
        payload
            .field("base_sketch_plan")
            .unwrap()
            .field("crossings")
            .unwrap()
            .items()
            .unwrap()
            .len(),
        4
    );
    fs::write(&input, recipe(4).to_string()).unwrap();
    assert!(!run().status.success());
    assert_eq!(bytes, fs::read_to_string(&output).unwrap());
    fs::remove_dir_all(directory).unwrap();
}
