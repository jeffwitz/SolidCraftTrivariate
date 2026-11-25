use solidcraft_trivariate::{
    json::{object, Json},
    solidcraft_bridge,
    solidcraft_polygon::import_polygon_timeline,
};
use std::{collections::BTreeMap, fs, process::Command};

fn num(n: impl ToString) -> Json {
    Json::Number(n.to_string())
}

fn run_ok(program: &str, args: &[&str]) -> String {
    let out = Command::new(program)
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("{program}: {e}"));
    assert!(
        out.status.success(),
        "{program} {args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}
fn host_step_case(
    name: &str,
    contours: &[Vec<[i32; 2]>],
    height: i32,
    expected_volume: f64,
    probes: &[[f64; 4]],
) {
    let host = std::env::var("SOLIDCRAFT_CLI").unwrap_or_else(|_| "solidcraft-cli".into());
    let freecad = std::env::var("FREECADCMD").unwrap_or_else(|_| "freecadcmd".into());
    if Command::new(&host).arg("--help").output().is_err()
        || Command::new(&freecad).arg("--help").output().is_err()
    {
        return;
    }
    let work = std::env::temp_dir().join(format!("solidcraft-sc5-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    fs::create_dir_all(&work).unwrap();
    let commands_json = contours
        .iter()
        .map(|contour| {
            Json::Array(
                contour
                    .iter()
                    .map(|p| Json::Array(p.iter().map(num).collect()))
                    .collect(),
            )
        })
        .collect::<Vec<_>>();
    let mut commands = vec![object([
        ("command", string("sketch.create")),
        (
            "params",
            object([("plane", string("XY")), ("name", string("Profile"))]),
        ),
    ])];
    for points in commands_json {
        commands.push(object([
            ("command", string("sketch.line")),
            (
                "params",
                object([("points", points), ("closed", Json::Bool(true))]),
            ),
        ]));
    }
    commands.push(object([
        ("command", string("sketch.finish")),
        ("params", object([])),
    ]));
    commands.push(object([
        ("command", string("solid.extrude")),
        (
            "params",
            object([("distance", num(height)), ("body_name", string(name))]),
        ),
    ]));
    let script = object([("commands", Json::Array(commands))]);
    let script_path = work.join("commands.json");
    let sc_path = work.join("editable.sc");
    let host_step = work.join("host.step");
    let timeline_path = work.join("history.json");
    let model_path = work.join("native.json");
    let native_step = work.join("native.step");
    fs::write(&script_path, script.to_string()).unwrap();
    run_ok(
        &host,
        &[
            "run",
            script_path.to_str().unwrap(),
            "--save",
            sc_path.to_str().unwrap(),
            "--out",
            host_step.to_str().unwrap(),
            "--quiet",
        ],
    );
    run_ok(
        &host,
        &[
            "script",
            sc_path.to_str().unwrap(),
            "--out",
            timeline_path.to_str().unwrap(),
        ],
    );
    let history =
        solidcraft_trivariate::json::parse(&fs::read_to_string(&timeline_path).unwrap()).unwrap();
    let model = solidcraft_bridge::import_script(&history).unwrap();
    let schema = model.field("schema").unwrap().text().unwrap();
    assert!(schema == "polygon_prism_v1" || schema == "polygonal_perforated_prism_v1");
    let native = env!("CARGO_BIN_EXE_solidcraft-trivariate");
    run_ok(
        native,
        &[
            "import-solidcraft",
            timeline_path.to_str().unwrap(),
            model_path.to_str().unwrap(),
        ],
    );
    run_ok(
        native,
        &[
            "export-step",
            model_path.to_str().unwrap(),
            native_step.to_str().unwrap(),
        ],
    );
    let protected = fs::read(&model_path).unwrap();
    let mut invalid = history.clone();
    feature_mut(&mut invalid, 1).insert("taper".into(), num(5));
    let invalid_path = work.join("invalid-history.json");
    fs::write(&invalid_path, invalid.to_string()).unwrap();
    let rejected = Command::new(native)
        .args([
            "import-solidcraft",
            invalid_path.to_str().unwrap(),
            model_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        !rejected.status.success(),
        "unsupported extrusion option was accepted"
    );
    assert_eq!(
        fs::read(&model_path).unwrap(),
        protected,
        "refusal changed the previous native model"
    );
    let macro_path = work.join("compare.FCMacro");
    let macro_source = r#"import FreeCAD as App, Part, os, json, math
a=Part.read(os.environ['SC5_HOST_STEP']); b=Part.read(os.environ['SC5_NATIVE_STEP'])
assert a.isValid() and b.isValid() and len(a.Solids)==1 and len(b.Solids)==1
assert abs(a.Volume-float(os.environ['SC5_EXPECTED_VOLUME']))<1e-8
assert abs(b.Volume-a.Volume)<1e-8 and abs(b.Area-a.Area)<1e-8
ca,cb=a.CenterOfMass,b.CenterOfMass
assert math.dist((ca.x,ca.y,ca.z),(cb.x,cb.y,cb.z))<1e-9
assert a.cut(b).Volume<1e-8 and b.cut(a).Volume<1e-8
probes=json.loads(os.environ['SC5_PROBES'])
for x,y,z,expected in probes:
    pa=Part.Vertex(App.Vector(x,y,z)); assert a.isInside(pa.Point,1e-8,False)==expected
    assert b.isInside(pa.Point,1e-8,False)==expected
print(json.dumps({'valid':True,'volume_mm3':a.Volume,'area_mm2':a.Area,'probes':len(probes)}))
"#;
    fs::write(&macro_path, macro_source).unwrap();
    let probe_json = Json::Array(
        probes
            .iter()
            .map(|p| {
                Json::Array(vec![
                    Json::Number(p[0].to_string()),
                    Json::Number(p[1].to_string()),
                    Json::Number(p[2].to_string()),
                    Json::Bool(p[3] != 0.0),
                ])
            })
            .collect(),
    )
    .to_string();
    let output = Command::new(&freecad)
        .arg(macro_path)
        .env("SC5_HOST_STEP", &host_step)
        .env("SC5_NATIVE_STEP", &native_step)
        .env("SC5_EXPECTED_VOLUME", expected_volume.to_string())
        .env("SC5_PROBES", probe_json)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "FreeCAD STEP comparison failed: {} {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let _ = fs::remove_dir_all(&work);
}
fn string(s: impl Into<String>) -> Json {
    Json::String(s.into())
}
fn make_history(contours: &[Vec<[i32; 2]>], height: i32) -> Json {
    let mut points = vec![object([
        ("fixed", Json::Bool(true)),
        ("id", string("origin")),
        ("pos", Json::Array(vec![num(0), num(0)])),
    ])];
    let mut curves = Vec::new();
    let mut curve_loops = Vec::new();
    let mut next_id = 1;
    for (loop_index, contour) in contours.iter().enumerate() {
        let start = next_id;
        for (i, xy) in contour.iter().enumerate() {
            let name = if i == 0 {
                format!("l{next_id}.start")
            } else {
                format!("l{}.end", next_id - 1)
            };
            points.push(object([
                ("fixed", Json::Bool(false)),
                ("id", string(name)),
                ("pos", Json::Array(xy.iter().map(num).collect())),
            ]));
            next_id += 1;
        }
        let ids = (0..contour.len())
            .map(|j| format!("l{}", start + j))
            .collect::<Vec<_>>();
        for (i, curve_id) in ids.iter().enumerate() {
            curves.push(object([
                ("a", num(start + i)),
                ("b", num(start + (i + 1) % contour.len())),
                ("construction", Json::Bool(false)),
                ("id", string(curve_id.clone())),
                ("type", string("line")),
            ]));
        }
        curve_loops.push(ids.into_iter().map(string).collect::<Vec<_>>());
        if loop_index == 0 && contour[0] == [0, 0] {
            // The host's origin snap is stored as a coincident constraint.
        }
    }
    let constraints = if !contours.is_empty() && contours[0].first() == Some(&[0, 0]) {
        vec![object([
            ("id", string("k1")),
            ("p", num(0)),
            ("q", num(1)),
            ("type", string("coincident")),
        ])]
    } else {
        vec![]
    };
    let sketch = object([
        ("constraints", Json::Array(constraints)),
        (
            "counters",
            object([
                (
                    "k",
                    num(
                        if contours.first().and_then(|c| c.first()) == Some(&[0, 0]) {
                            1
                        } else {
                            0
                        },
                    ),
                ),
                ("l", num(curves.len())),
            ]),
        ),
        ("curves", Json::Array(curves)),
        ("points", Json::Array(points)),
    ]);
    let sketch_feature = object([
        ("id", num(1)),
        ("name", string("Profile")),
        (
            "plane",
            object([("name", string("XY")), ("type", string("origin"))]),
        ),
        ("sketch", sketch),
        ("suppressed", Json::Bool(false)),
        ("type", string("sketch")),
    ]);
    let profiles = if contours.len() == 1 {
        object([("type", string("all"))])
    } else {
        object([
            (
                "loops",
                Json::Array(vec![Json::Array(curve_loops[0].clone())]),
            ),
            ("type", string("curves")),
        ])
    };
    let extrude = object([
        ("body_names", Json::Array(vec![string("Prism")])),
        (
            "extent",
            object([
                ("convention", string("total_distance")),
                ("direction", string("positive")),
                ("distance", string(height.to_string())),
            ]),
        ),
        ("id", num(2)),
        ("name", string("Extrude1")),
        ("operation", string("new_body")),
        ("param_names", Json::Array(vec![string("d1")])),
        ("profiles", profiles),
        ("sketch", num(1)),
        ("suppressed", Json::Bool(false)),
        ("type", string("extrude")),
    ]);
    object([
        (
            "about",
            object([
                ("design", string("editable")),
                ("format", string("solidcraft/66")),
                ("replay", string("host script")),
                ("version", string("0.48.1")),
            ]),
        ),
        (
            "commands",
            Json::Array(vec![
                object([
                    ("command", string("timeline.add")),
                    ("params", object([("feature", sketch_feature)])),
                ]),
                object([
                    ("command", string("timeline.add")),
                    ("params", object([("feature", extrude)])),
                ]),
            ]),
        ),
    ])
}

fn feature_mut(script: &mut Json, command_index: usize) -> &mut BTreeMap<String, Json> {
    let Json::Object(root) = script else {
        unreachable!()
    };
    let Json::Array(commands) = root.get_mut("commands").unwrap() else {
        unreachable!()
    };
    let Json::Object(command) = &mut commands[command_index] else {
        unreachable!()
    };
    let Json::Object(params) = command.get_mut("params").unwrap() else {
        unreachable!()
    };
    let Json::Object(feature) = params.get_mut("feature").unwrap() else {
        unreachable!()
    };
    feature
}

#[test]
fn imports_simple_and_perforated_histories_to_the_right_native_family() {
    let simple = make_history(&[vec![[0, 0], [16, 0], [16, 12], [0, 12]]], 8);
    let model = solidcraft_bridge::import_script(&simple).unwrap();
    assert_eq!(
        model.field("schema").unwrap().text().unwrap(),
        "polygon_prism_v1"
    );
    let perforated = make_history(
        &[
            vec![[0, 0], [20, 0], [20, 20], [0, 20]],
            vec![[5, 5], [5, 10], [10, 10], [10, 5]],
        ],
        8,
    );
    let model = import_polygon_timeline(&perforated).unwrap();
    assert_eq!(
        model.field("schema").unwrap().text().unwrap(),
        "polygonal_perforated_prism_v1"
    );
    assert_eq!(model.field("blocks").unwrap().items().unwrap().len(), 24);
}

#[test]
fn real_solidcraft_sources_export_matching_step_for_g17_and_g18() {
    host_step_case(
        "concave-g17",
        &[vec![[0, 0], [12, 0], [12, 3], [4, 3], [4, 10], [0, 10]]],
        8,
        512.0,
        &[
            [2.0, 8.0, 4.0, 1.0],
            [9.0, 8.0, 4.0, 0.0],
            [2.0, 2.0, 4.0, 1.0],
        ],
    );
    host_step_case(
        "square-hole-g18",
        &[
            vec![[0, 0], [20, 0], [20, 20], [0, 20]],
            vec![[5, 5], [5, 10], [10, 10], [10, 5]],
        ],
        8,
        3000.0,
        &[
            [2.0, 2.0, 4.0, 1.0],
            [7.5, 7.5, 4.0, 0.0],
            [18.0, 18.0, 4.0, 1.0],
            [22.0, 2.0, 4.0, 0.0],
        ],
    );
    host_step_case(
        "two-holes-g18",
        &[
            vec![[0, 0], [30, 0], [30, 20], [0, 20]],
            vec![[4, 4], [4, 8], [8, 8], [8, 4]],
            vec![[20, 4], [24, 4], [24, 8], [20, 8]],
        ],
        6,
        3408.0,
        &[
            [2.0, 2.0, 3.0, 1.0],
            [6.0, 6.0, 3.0, 0.0],
            [22.0, 6.0, 3.0, 0.0],
            [15.0, 10.0, 3.0, 1.0],
        ],
    );
}

#[test]
fn rejects_unclosed_fractional_and_unsupported_extrusions() {
    let good = make_history(&[vec![[0, 0], [16, 0], [16, 12], [0, 12]]], 8);
    let mut broken = good.clone();
    if let Json::Object(root) = &mut broken {
        if let Json::Array(commands) = root.get_mut("commands").unwrap() {
            if let Json::Object(command) = &mut commands[0] {
                if let Json::Object(params) = command.get_mut("params").unwrap() {
                    if let Json::Object(feature) = params.get_mut("feature").unwrap() {
                        if let Json::Object(sketch) = feature.get_mut("sketch").unwrap() {
                            if let Json::Array(curves) = sketch.get_mut("curves").unwrap() {
                                curves.pop();
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(import_polygon_timeline(&broken).is_err());
    let mut fractional = good.clone();
    if let Json::Object(root) = &mut fractional {
        if let Json::Array(commands) = root.get_mut("commands").unwrap() {
            if let Json::Object(command) = &mut commands[0] {
                if let Json::Object(params) = command.get_mut("params").unwrap() {
                    if let Json::Object(feature) = params.get_mut("feature").unwrap() {
                        if let Json::Object(sketch) = feature.get_mut("sketch").unwrap() {
                            if let Json::Array(points) = sketch.get_mut("points").unwrap() {
                                if let Json::Object(point) = &mut points[2] {
                                    point.insert(
                                        "pos".into(),
                                        Json::Array(vec![Json::Number("16.25".into()), num(0)]),
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(import_polygon_timeline(&fractional).is_err());
    let mut drafted = good;
    if let Json::Object(root) = &mut drafted {
        if let Json::Array(commands) = root.get_mut("commands").unwrap() {
            if let Json::Object(command) = &mut commands[1] {
                if let Json::Object(params) = command.get_mut("params").unwrap() {
                    if let Json::Object(feature) = params.get_mut("feature").unwrap() {
                        feature.insert("taper".into(), num(2));
                    }
                }
            }
        }
    }
    assert!(import_polygon_timeline(&drafted).is_err());
}

#[test]
fn rejects_expressions_constraints_arcs_and_inconsistent_selection() {
    let good = make_history(&[vec![[0, 0], [16, 0], [16, 12], [0, 12]]], 8);
    for mutation in [
        "arc",
        "construction",
        "constraint",
        "profile",
        "expression",
        "reverse",
        "wrong_sketch",
        "extra_step",
    ] {
        let mut bad = good.clone();
        match mutation {
            "arc" | "construction" => {
                let Json::Object(sketch) = feature_mut(&mut bad, 0).get_mut("sketch").unwrap()
                else {
                    unreachable!()
                };
                let Json::Array(curves) = sketch.get_mut("curves").unwrap() else {
                    unreachable!()
                };
                let Json::Object(curve) = &mut curves[0] else {
                    unreachable!()
                };
                curve.insert(
                    if mutation == "arc" {
                        "type"
                    } else {
                        "construction"
                    }
                    .into(),
                    if mutation == "arc" {
                        string("arc")
                    } else {
                        Json::Bool(true)
                    },
                );
            }
            "constraint" => {
                let Json::Object(sketch) = feature_mut(&mut bad, 0).get_mut("sketch").unwrap()
                else {
                    unreachable!()
                };
                let Json::Array(constraints) = sketch.get_mut("constraints").unwrap() else {
                    unreachable!()
                };
                constraints.push(object([
                    ("id", string("bad")),
                    ("p", num(1)),
                    ("q", num(2)),
                    ("type", string("distance")),
                ]));
            }
            "profile" => {
                feature_mut(&mut bad, 1).insert(
                    "profiles".into(),
                    object([
                        ("type", string("curves")),
                        (
                            "loops",
                            Json::Array(vec![Json::Array(vec![string("l1"), string("l3")])]),
                        ),
                    ]),
                );
            }
            "expression" => {
                let Json::Object(extent) = feature_mut(&mut bad, 1).get_mut("extent").unwrap()
                else {
                    unreachable!()
                };
                extent.insert("distance".into(), string("d1 / 2"));
            }
            "reverse" => {
                let Json::Object(extent) = feature_mut(&mut bad, 1).get_mut("extent").unwrap()
                else {
                    unreachable!()
                };
                extent.insert("direction".into(), string("negative"));
            }
            "wrong_sketch" => {
                feature_mut(&mut bad, 1).insert("sketch".into(), num(42));
            }
            "extra_step" => {
                let Json::Object(root) = &mut bad else {
                    unreachable!()
                };
                let Json::Array(commands) = root.get_mut("commands").unwrap() else {
                    unreachable!()
                };
                commands.push(commands[1].clone());
            }
            _ => unreachable!(),
        }
        assert!(
            import_polygon_timeline(&bad).is_err(),
            "unsupported history mutation accepted: {mutation}"
        );
    }
}

#[test]
fn non_two_feature_sketch_timeline_still_routes_to_sc2_adapter() {
    let mut source = make_history(&[vec![[0, 0], [16, 0], [16, 12], [0, 12]]], 8);
    let Json::Object(root) = &mut source else {
        unreachable!()
    };
    let Json::Array(commands) = root.get_mut("commands").unwrap() else {
        unreachable!()
    };
    commands.push(commands[1].clone());
    assert_eq!(
        solidcraft_bridge::import_script(&source).unwrap_err(),
        "SolidCraft bracket timeline must contain five active features"
    );
}
