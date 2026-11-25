use solidcraft_trivariate::{
    json::{object, Json},
    solidcraft_bridge,
};
use std::{fs, process::Command};

fn n(v: impl ToString) -> Json {
    Json::Number(v.to_string())
}
fn s(v: impl Into<String>) -> Json {
    Json::String(v.into())
}
fn run(program: &str, args: &[&str]) -> String {
    let out = Command::new(program).args(args).output().unwrap();
    assert!(
        out.status.success(),
        "{program} {args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

fn host_case(name: &str, outer: &[[i32; 2]], bores: &[([i32; 2], i32)], height: i32) {
    let host = std::env::var("SOLIDCRAFT_CLI").unwrap_or_else(|_| "solidcraft-cli".into());
    let freecad = std::env::var("FREECADCMD").unwrap_or_else(|_| "freecadcmd".into());
    if Command::new(&host).arg("--help").output().is_err()
        || Command::new(&freecad).arg("--help").output().is_err()
    {
        eprintln!("SC6 host/STEP differential skipped: host or FreeCAD is unavailable");
        return;
    }
    let work = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("build/solidcraft-sc6")
        .join(name);
    fs::create_dir_all(&work).unwrap();
    let outer_points = outer
        .iter()
        .map(|p| Json::Array(p.iter().map(|x| n(*x)).collect()))
        .collect();
    let mut commands = vec![
        object([
            ("command", s("sketch.create")),
            (
                "params",
                object([("plane", s("XY")), ("name", s("SC6 profile"))]),
            ),
        ]),
        object([
            ("command", s("sketch.line")),
            (
                "params",
                object([
                    ("points", Json::Array(outer_points)),
                    ("closed", Json::Bool(true)),
                ]),
            ),
        ]),
    ];
    for (center, radius) in bores {
        commands.push(object([
            ("command", s("sketch.circle.center")),
            (
                "params",
                object([
                    (
                        "center",
                        Json::Array(center.iter().map(|x| n(*x)).collect()),
                    ),
                    ("radius", n(*radius)),
                ]),
            ),
        ]));
    }
    commands.push(object([
        ("command", s("sketch.finish")),
        ("params", object([])),
    ]));
    commands.push(object([
        ("command", s("solid.extrude")),
        (
            "params",
            object([("distance", n(height)), ("body_name", s(name))]),
        ),
    ]));
    let commands_path = work.join("commands.json");
    let sc_path = work.join("editable.sc");
    let host_step = work.join("host.step");
    let history_path = work.join("history.json");
    let model_path = work.join("model.json");
    let native_step = work.join("native.step");
    fs::write(
        &commands_path,
        object([("commands", Json::Array(commands))]).to_string(),
    )
    .unwrap();
    run(
        &host,
        &[
            "run",
            commands_path.to_str().unwrap(),
            "--save",
            sc_path.to_str().unwrap(),
            "--out",
            host_step.to_str().unwrap(),
            "--quiet",
        ],
    );
    run(
        &host,
        &[
            "script",
            sc_path.to_str().unwrap(),
            "--out",
            history_path.to_str().unwrap(),
        ],
    );
    let history =
        solidcraft_trivariate::json::parse(&fs::read_to_string(&history_path).unwrap()).unwrap();
    let model = solidcraft_bridge::import_script(&history)
        .unwrap_or_else(|error| panic!("history import failed: {error}; {history}"));
    assert_eq!(
        model.field("schema").unwrap().text().unwrap(),
        "polygon_circular_bores_v1"
    );
    assert_eq!(
        model.field("bores").unwrap().items().unwrap().len(),
        bores.len()
    );
    run(
        env!("CARGO_BIN_EXE_solidcraft-trivariate"),
        &[
            "import-solidcraft",
            history_path.to_str().unwrap(),
            model_path.to_str().unwrap(),
        ],
    );
    run(
        env!("CARGO_BIN_EXE_solidcraft-trivariate"),
        &[
            "export-step",
            model_path.to_str().unwrap(),
            native_step.to_str().unwrap(),
        ],
    );
    let area2: i64 = (0..outer.len())
        .map(|i| {
            let a = outer[i];
            let b = outer[(i + 1) % outer.len()];
            i64::from(a[0]) * i64::from(b[1]) - i64::from(a[1]) * i64::from(b[0])
        })
        .sum();
    let sum_r2: i64 = bores
        .iter()
        .map(|(_, r)| i64::from(*r) * i64::from(*r))
        .sum();
    let expected = (area2 as f64 / 2.0 - std::f64::consts::PI * sum_r2 as f64) * height as f64;
    let z = f64::from(height) / 2.0;
    let mut probes = Vec::new();
    for (center, _) in bores {
        probes.push(vec![n(center[0]), n(center[1]), n(z), Json::Bool(false)]);
    }
    let material = if name == "concave" {
        [-40, 0]
    } else {
        [bores[0].0[0] + bores[0].1 + 1, bores[0].0[1]]
    };
    probes.push(vec![n(material[0]), n(material[1]), n(z), Json::Bool(true)]);
    probes.push(vec![n(70), n(70), n(z), Json::Bool(false)]);
    let macro_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/references/solidcraft_circular_polygon_oracle.FCMacro");
    let metrics_path = work.join("step-metrics.json");
    let out = Command::new(freecad)
        .arg(macro_path)
        .env("SC6_HOST_STEP", &host_step)
        .env("SC6_NATIVE_STEP", &native_step)
        .env("SC6_EXPECTED", expected.to_string())
        .env("SC6_METRICS", &metrics_path)
        .env(
            "SC6_PROBES",
            Json::Array(probes.into_iter().map(Json::Array).collect()).to_string(),
        )
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "FreeCAD STEP differential failed: {} {}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    let metrics_text = fs::read_to_string(&metrics_path).unwrap_or_else(|error| {
        panic!(
            "FreeCAD returned success but no metrics file ({error}); stdout={} stderr={}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        )
    });
    let metrics = solidcraft_trivariate::json::parse(&metrics_text).unwrap();
    assert_eq!(metrics.field("valid").unwrap(), &Json::Bool(true));
    println!("SC6 {name}: {metrics}");
    let protected = fs::read(&model_path).unwrap();
    let source = fs::read_to_string(&history_path).unwrap();
    let radius_token = format!("\"r\": {}.0", bores[0].1);
    let invalid = source.replacen(&radius_token, &format!("\"r\": {}.5", bores[0].1), 1);
    assert_ne!(source, invalid, "host radius mutation target not found");
    let invalid_path = work.join("invalid-history.json");
    fs::write(&invalid_path, invalid).unwrap();
    let rejected = Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
        .args([
            "import-solidcraft",
            invalid_path.to_str().unwrap(),
            model_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        !rejected.status.success(),
        "fractional-mm bore radius was accepted"
    );
    assert_eq!(
        fs::read(model_path).unwrap(),
        protected,
        "refusal replaced previous audited output"
    );
}

#[test]
fn real_solidcraft_circular_profiles_match_independent_step() {
    let outer = [[-48, -48], [48, -48], [48, 48], [-48, 48]];
    host_case("one", &outer, &[([0, 0], 6)], 8);
    host_case("two", &outer, &[([-18, 0], 4), ([18, 0], 4)], 8);
    host_case(
        "four",
        &outer,
        &[
            ([-20, -20], 3),
            ([20, -20], 3),
            ([-20, 20], 3),
            ([20, 20], 3),
        ],
        8,
    );
    host_case(
        "concave",
        &[
            [-48, -48],
            [48, -48],
            [48, -24],
            [0, -24],
            [0, 48],
            [-48, 48],
        ],
        &[([-24, 0], 4)],
        8,
    );
}

#[test]
fn frozen_real_host_fixture_imports_and_refusal_is_atomic() {
    let source = include_str!("references/solidcraft_circular_polygon_one_hole.history.json");
    let history = solidcraft_trivariate::json::parse(source).unwrap();
    let model = solidcraft_bridge::import_script(&history).unwrap();
    assert_eq!(
        model.field("schema").unwrap().text().unwrap(),
        "polygon_circular_bores_v1"
    );
    let work = std::env::temp_dir().join(format!("solidcraft-sc6-fixture-{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    fs::create_dir_all(&work).unwrap();
    let history_path = work.join("history.json");
    let model_path = work.join("model.json");
    fs::write(&history_path, source).unwrap();
    run(
        env!("CARGO_BIN_EXE_solidcraft-trivariate"),
        &[
            "import-solidcraft",
            history_path.to_str().unwrap(),
            model_path.to_str().unwrap(),
        ],
    );
    let protected = fs::read(&model_path).unwrap();
    let invalid = source.replace("\"r\": 6.0", "\"r\": 6.5");
    assert_ne!(
        source, invalid,
        "frozen host radius mutation target not found"
    );
    fs::write(&history_path, invalid).unwrap();
    let rejected = Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
        .args([
            "import-solidcraft",
            history_path.to_str().unwrap(),
            model_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        !rejected.status.success(),
        "fractional-mm radius was accepted"
    );
    assert_eq!(
        fs::read(&model_path).unwrap(),
        protected,
        "refusal replaced the previously audited output"
    );
    let _ = fs::remove_dir_all(work);
}
