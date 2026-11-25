use solidcraft_trivariate::{
    json::{self, Json},
    solidcraft_bridge,
};
use std::{fs, process::Command};

fn run(program: &str, args: &[&str]) -> String {
    let out = Command::new(program).args(args).output().unwrap();
    assert!(
        out.status.success(),
        "{program} {args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

struct Case {
    id: &'static str,
    plane: &'static str,
    hole: bool,
    bounds: [f64; 6],
}

fn reference(name: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/references/solidcraft-oriented-plane")
        .join(name)
}

fn world_probes(case: &Case) -> Json {
    let z_mid = 3.0;
    let (material, void, exterior) = match case.plane {
        "XZ" => ([11.0, 3.0, -11.0], [15.0, 3.0, -18.0], [30.0, 3.0, -30.0]),
        _ => ([3.0, 11.0, 11.0], [3.0, 15.0, 18.0], [3.0, 30.0, 30.0]),
    };
    let inside_at_void = !case.hole;
    Json::Array(vec![
        Json::Array(vec![
            Json::Number(material[0].to_string()),
            Json::Number(material[1].to_string()),
            Json::Number(material[2].to_string()),
            Json::Bool(true),
        ]),
        Json::Array(vec![
            Json::Number(void[0].to_string()),
            Json::Number(void[1].to_string()),
            Json::Number(void[2].to_string()),
            Json::Bool(inside_at_void),
        ]),
        Json::Array(vec![
            Json::Number(exterior[0].to_string()),
            Json::Number(exterior[1].to_string()),
            Json::Number(exterior[2].to_string()),
            Json::Bool(false),
        ]),
        Json::Array(vec![
            Json::Number(case.bounds[0].to_string()),
            Json::Number(case.bounds[2].to_string()),
            Json::Number(z_mid.to_string()),
            Json::Bool(false),
        ]),
    ])
}

fn host_step_case(case: &Case) -> Json {
    let host = std::env::var("SOLIDCRAFT_CLI").unwrap_or_else(|_| "solidcraft-cli".into());
    let freecad = std::env::var("FREECADCMD").unwrap_or_else(|_| "freecadcmd".into());
    if Command::new(&host).arg("--help").output().is_err()
        || Command::new(&freecad).arg("--help").output().is_err()
    {
        eprintln!("SC7 host/STEP oracle skipped: SolidCraft or FreeCAD is unavailable");
        return Json::Null;
    }
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let work = root.join("build/solidcraft-sc7").join(case.id);
    fs::create_dir_all(&work).unwrap();
    let design = reference(&format!("{}.sc", case.id));
    let frozen_history =
        json::parse(&fs::read_to_string(reference(&format!("{}.history.json", case.id))).unwrap())
            .unwrap();
    let host_step = work.join("host.step");
    let extracted_history = work.join("history.json");
    let native_model = work.join("native.json");
    let native_step = work.join("native.step");
    run(
        &host,
        &[
            "run",
            design.to_str().unwrap(),
            "--out",
            host_step.to_str().unwrap(),
            "--quiet",
        ],
    );
    run(
        &host,
        &[
            "script",
            design.to_str().unwrap(),
            "--out",
            extracted_history.to_str().unwrap(),
        ],
    );
    let extracted = json::parse(&fs::read_to_string(&extracted_history).unwrap()).unwrap();
    assert_eq!(
        extracted, frozen_history,
        "frozen history changed on host replay: {}",
        case.id
    );
    let model = solidcraft_bridge::import_script(&extracted).unwrap();
    assert_eq!(
        model.field("schema").unwrap().text().unwrap(),
        "affine_image_v1"
    );
    run(
        env!("CARGO_BIN_EXE_solidcraft-trivariate"),
        &[
            "import-solidcraft",
            extracted_history.to_str().unwrap(),
            native_model.to_str().unwrap(),
        ],
    );
    run(
        env!("CARGO_BIN_EXE_solidcraft-trivariate"),
        &[
            "export-step",
            native_model.to_str().unwrap(),
            native_step.to_str().unwrap(),
        ],
    );
    let expected_volume = if case.hole {
        1008.0 - 6.0 * std::f64::consts::PI
    } else {
        1008.0
    };
    let macro_path = root.join("tests/references/solidcraft_oriented_plane_oracle.FCMacro");
    let metrics_path = work.join("metrics.json");
    let output = Command::new(freecad)
        .arg(macro_path)
        .env("SC7_HOST_STEP", &host_step)
        .env("SC7_NATIVE_STEP", &native_step)
        .env("SC7_EXPECTED_VOLUME", expected_volume.to_string())
        .env(
            "SC7_EXPECTED_BOUNDS",
            Json::Array(
                case.bounds
                    .iter()
                    .map(|x| Json::Number(x.to_string()))
                    .collect(),
            )
            .to_string(),
        )
        .env("SC7_PROBES", world_probes(case).to_string())
        .env("SC7_METRICS", &metrics_path)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "SC7 STEP comparison failed for {}: {} {}",
        case.id,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let metrics = json::parse(&fs::read_to_string(&metrics_path).unwrap()).unwrap();
    assert_eq!(metrics.field("valid").unwrap(), &Json::Bool(true));
    metrics
}

#[test]
fn origin_xz_yz_histories_use_the_right_handed_g19_images() {
    let cases = [
        Case {
            id: "xz_polygon",
            plane: "XZ",
            hole: false,
            bounds: [10., 22., 0., 6., -28., -10.],
        },
        Case {
            id: "yz_polygon",
            plane: "YZ",
            hole: false,
            bounds: [0., 6., 10., 22., 10., 28.],
        },
        Case {
            id: "xz_hole",
            plane: "XZ",
            hole: true,
            bounds: [10., 22., 0., 6., -28., -10.],
        },
        Case {
            id: "yz_hole",
            plane: "YZ",
            hole: true,
            bounds: [0., 6., 10., 22., 10., 28.],
        },
    ];
    for case in cases {
        let history = json::parse(
            &fs::read_to_string(reference(&format!("{}.history.json", case.id))).unwrap(),
        )
        .unwrap();
        let model = solidcraft_bridge::import_script(&history)
            .unwrap_or_else(|error| panic!("{} import failed: {error}", case.id));
        assert_eq!(
            model.field("schema").unwrap().text().unwrap(),
            "affine_image_v1"
        );
        assert_eq!(
            model
                .field("source_model")
                .unwrap()
                .field("schema")
                .unwrap()
                .text()
                .unwrap(),
            if case.hole {
                "polygon_circular_bores_v1"
            } else {
                "polygon_prism_v1"
            }
        );
        let transform = model.field("transform").unwrap();
        let matrix = transform.field("matrix").unwrap();
        let expected = if case.plane == "XZ" {
            [[1, 0, 0], [0, 0, 1], [0, -1, 0]]
        } else {
            [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
        };
        for (row, values) in matrix.items().unwrap().iter().zip(expected) {
            assert_eq!(
                row,
                &Json::Array(
                    values
                        .into_iter()
                        .map(|x| Json::Number(x.to_string()))
                        .collect()
                )
            );
        }
        if !std::env::var("SC7_SKIP_HOST_ORACLE").is_ok_and(|v| v == "1") {
            let metrics = host_step_case(&case);
            if metrics != Json::Null {
                let directory =
                    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("build/solidcraft-sc7");
                fs::write(
                    directory.join(format!("{}-metrics.json", case.id)),
                    metrics.to_string(),
                )
                .unwrap();
            }
        }
    }
}

#[test]
fn real_host_negative_distance_and_custom_plane_refuse_atomically() {
    let work = std::env::temp_dir().join(format!("solidcraft-sc7-refusal-{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    fs::create_dir_all(&work).unwrap();
    let good = reference("xz_polygon.history.json");
    let output_model = work.join("existing.json");
    run(
        env!("CARGO_BIN_EXE_solidcraft-trivariate"),
        &[
            "import-solidcraft",
            good.to_str().unwrap(),
            output_model.to_str().unwrap(),
        ],
    );
    let protected = fs::read(&output_model).unwrap();
    for (id, reason) in [
        (
            "xz_negative",
            "SolidCraft extrusion height outside native range",
        ),
        ("xz_offset", "SC7 supports origin planes only"),
    ] {
        let history = reference(&format!("{id}.history.json"));
        let rejected = Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
            .args([
                "import-solidcraft",
                history.to_str().unwrap(),
                output_model.to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(
            !rejected.status.success(),
            "SC7 accepted refusal case {id} ({reason})"
        );
        let diagnostic = json::parse(&String::from_utf8_lossy(&rejected.stdout)).unwrap();
        assert_eq!(diagnostic.field("reason").unwrap().text().unwrap(), reason);
        assert_eq!(
            fs::read(&output_model).unwrap(),
            protected,
            "refusal {id} replaced the prior output"
        );
    }
    let _ = fs::remove_dir_all(work);
}
