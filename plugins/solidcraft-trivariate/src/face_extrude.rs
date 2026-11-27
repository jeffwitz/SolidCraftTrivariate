//! Resolve a real face attachment, construct the extrusion in its local frame,
//! and commit the result in the host's single plugin transaction.
use super::*;
use solidcraft_trivariate::face_profile_extrusion::{self, FaceProfileExtrusion};

type V = [f64; 3];
type Frame = [V; 4];

fn cmd(name: &str, args: Value) -> Result<Value, PluginError> {
    let reply = commands::run(name, &args.to_string())?;
    let result: Value = serde_json::from_str(&reply).map_err(|e| invalid(name, e))?;
    if let Some(err) = result.get("error").filter(|err| !err.is_null()) {
        return Err(invalid(name, err));
    }
    Ok(result)
}

fn xyz(v: &Value) -> Result<V, PluginError> {
    serde_json::from_value(v.clone())
        .map_err(|_| invalid("plane", "missing sketch-frame coordinates"))
}
fn dot(a: V, b: V) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn minus(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] - b[i])
}

fn frame(plane: &Value) -> Result<Frame, PluginError> {
    let (o, x, y) = if plane["type"] == "face" {
        let p = &plane["plane"];
        (xyz(&p["origin"])?, xyz(&p["x_dir"])?, xyz(&p["y_dir"])?)
    } else if plane["type"] == "origin" {
        let (x, y) = match plane["name"].as_str() {
            Some("XY") => ([1., 0., 0.], [0., 1., 0.]),
            Some("XZ") => ([1., 0., 0.], [0., 0., -1.]),
            Some("YZ") => ([0., 1., 0.], [0., 0., 1.]),
            _ => return Err(invalid("plane", "unknown origin plane")),
        };
        ([0.; 3], x, y)
    } else {
        return Err(invalid(
            "plane",
            "select an origin plane or an attached planar face",
        ));
    };
    if (dot(x, x) - 1.).abs() > 1e-9 || (dot(y, y) - 1.).abs() > 1e-9 || dot(x, y).abs() > 1e-9 {
        return Err(invalid("plane", "sketch frame is not orthonormal"));
    }
    Ok([o, x, y, cross(x, y)])
}

fn flip(frame: &mut Frame, spans: &mut [Vec<V>]) {
    for p in spans.iter_mut().flatten() {
        p[1] = -p[1];
    }
    frame[2] = frame[2].map(|v| -v);
    frame[3] = frame[3].map(|v| -v);
}

// Six outward planes, one at each bound, each with its full rectangular
// area; the selected support must be the upper face in the local frame.
fn box_section(faces: &[Value], f: &Frame) -> Option<([[f64; 2]; 2], f64)> {
    if faces.len() != 6 {
        return None;
    }
    let mut bounds = [[f64::NAN; 2]; 3];
    let mut sides = Vec::new();
    for face in faces {
        let normal = xyz(&face["plane_normal"]).ok()?;
        let axis = (0..3).find(|&i| (dot(normal, f[i + 1]).abs() - 1.).abs() < 1e-12)?;
        let end = usize::from(dot(normal, f[axis + 1]) > 0.);
        if bounds[axis][end].is_finite() {
            return None;
        }
        bounds[axis][end] = dot(minus(xyz(&face["centroid"]).ok()?, f[0]), f[axis + 1]);
        sides.push((axis, face["area"].as_f64()?));
    }
    let size = bounds.map(|[lo, hi]| hi - lo);
    if size.iter().any(|d| !d.is_finite() || *d <= 0.) || bounds[2][1].abs() > 1e-8 {
        return None;
    }
    for (i, area) in sides {
        let expected = size[(i + 1) % 3] * size[(i + 2) % 3];
        if (area - expected).abs() > expected * 1e-9 {
            return None;
        }
    }
    Some(([bounds[0], bounds[1]], size[2]))
}

fn model_body_for(id: u64) -> Result<String, PluginError> {
    let data = cmd("model.body_features", json!({}))?;
    let names: Vec<_> = data["bodies"]
        .as_array()
        .into_iter()
        .flatten()
        .filter(|row| {
            row["features"]
                .as_array()
                .is_some_and(|fs| fs.iter().any(|f| f["id"] == id))
        })
        .filter_map(|row| row["body"].as_str())
        .collect();
    match names.as_slice() {
        [name] => Ok(name.to_string()),
        _ => Err(invalid("body", "cannot resolve one result body")),
    }
}

fn volume(body: &str) -> Result<f64, PluginError> {
    let m = cmd("inspect.measure", json!({"bodies":[body]}))?;
    m["bodies"]
        .as_array()
        .filter(|b| b.len() == 1)
        .and_then(|b| {
            (b[0]["closed"] == true)
                .then(|| b[0]["volume_mm3"].as_f64())
                .flatten()
        })
        .ok_or_else(|| invalid("body", "expected a single closed material body"))
}

pub(super) fn run(params: &Value) -> Result<String, PluginError> {
    let timeline = export_timeline()?;
    let requested = params
        .get("source_feature_id")
        .and_then(Value::as_u64)
        .unwrap_or(0);
    let sketch = timeline
        .pointer("/script/commands")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .rev()
        .filter_map(|step| step.pointer("/params/feature"))
        .find(|f| {
            f["type"] == "sketch"
                && f["suppressed"] != true
                && !f["name"]
                    .as_str()
                    .is_some_and(|n| n.starts_with("SCT Partition "))
                && (requested == 0 || f["id"] == requested)
        })
        .ok_or_else(|| invalid("sketch", "no matching finished sketch"))?;
    let id = sketch["id"]
        .as_u64()
        .ok_or_else(|| invalid("sketch", "missing sketch id"))?;
    let plane = &sketch["plane"];
    let attached = plane["type"] == "face";
    let operation = match params
        .get("operation")
        .and_then(Value::as_str)
        .unwrap_or("auto")
    {
        "auto" if attached => "add",
        "auto" | "new" => "new",
        "add" if attached => "add",
        "cut" if attached => "cut",
        _ => {
            return Err(invalid(
                "operation",
                "add and cut require a sketch attached to a planar body face",
            ));
        }
    };
    let depth = params
        .get("height_mm")
        .and_then(Value::as_f64)
        .unwrap_or(12.);
    if !depth.is_finite() || depth <= 0. || depth > 4096. {
        return Err(invalid(
            "height_mm",
            "depth must be positive and at most 4096 mm",
        ));
    }
    // Retain the exact existing polygon path when it admits this input.
    if operation == "new" && plane["name"] == "XY" && depth.fract() == 0. && depth <= 64. {
        let input = solidcraft_trivariate::json::parse(&timeline.to_string())
            .map_err(|e| invalid("timeline", e))?;
        if let Ok(model) = solidcraft_trivariate::solidcraft_polygon::extrude_source_sketch(
            &input,
            u32::try_from(id).map_err(|_| invalid("sketch", "id overflow"))?,
            depth as i32,
        ) {
            return insert_audited_model(&model);
        }
    }
    let mut frame = frame(plane)?;
    let mut spans = face_spline::spans(&sketch["sketch"]).map_err(|e| invalid("contour", e))?;
    let mut section = None;
    let mut target = "";
    let mut initial_volume = 0.;
    if attached {
        target = plane
            .pointer("/anchor/target/body")
            .and_then(Value::as_str)
            .ok_or_else(|| invalid("face", "missing owner"))?;
        if plane
            .pointer("/anchor/path")
            .and_then(Value::as_array)
            .is_some_and(|p| !p.is_empty())
        {
            return Err(invalid(
                "face",
                "placed component sketches require a world transform",
            ));
        }
        let faces = cmd("model.faces", json!({"body":target}))?;
        let rows = faces["faces"]
            .as_array()
            .ok_or_else(|| invalid("face", "missing face table"))?;
        let name = plane
            .pointer("/anchor/target/name")
            .and_then(Value::as_str)
            .ok_or_else(|| invalid("face", "missing face identity"))?;
        let support = rows
            .iter()
            .find(|r| r["name"] == name)
            .or_else(|| {
                plane["name"]
                    .as_str()
                    .and_then(|alias| rows.iter().find(|r| r["name"] == alias))
            })
            .ok_or_else(|| invalid("face", "support face no longer exists"))?;
        let normal = xyz(&support["plane_normal"])?;
        let sign = dot(frame[3], normal);
        if (sign.abs() - 1.).abs() > 1e-8
            || dot(minus(xyz(&support["centroid"])?, frame[0]), normal).abs() > 1e-7
        {
            return Err(invalid(
                "face",
                "attachment frame differs from the current face",
            ));
        }
        if sign < 0. {
            flip(&mut frame, &mut spans);
        }
        initial_volume = volume(target)?;
        section = box_section(rows, &frame);
    }
    let mut stock_diagnostic = None;
    let mut stock_model = None;
    if let Some((bounds, thickness)) = section.filter(|_| operation != "new") {
        let recipe = FaceProfileExtrusion {
            spans: spans.clone(),
            bounds,
            thickness,
            depth,
            operation: operation.into(),
            frame,
            tolerance: 1e-6,
        };
        match face_profile_extrusion::build(&recipe) {
            Ok(model) => stock_model = Some(model),
            Err(message) => stock_diagnostic = Some(message),
        }
    }
    let full = stock_model.is_some();
    let model = if let Some(model) = stock_model {
        model
    } else {
        if operation == "cut" {
            flip(&mut frame, &mut spans);
        }
        let mut bounds = [[f64::INFINITY, f64::NEG_INFINITY]; 2];
        for p in spans.iter().flatten() {
            for j in 0..2 {
                let x = p[j] / p[2];
                bounds[j][0] = bounds[j][0].min(x - 1.);
                bounds[j][1] = bounds[j][1].max(x + 1.);
            }
        }
        let recipe = FaceProfileExtrusion {
            spans,
            bounds,
            thickness: depth,
            depth,
            operation: "tool".into(),
            frame,
            tolerance: 1e-6,
        };
        face_profile_extrusion::build(&recipe).map_err(|e| invalid("extrusion", e))?
    };
    // The STEP exporter replays the source recipe and verifies boundary incidence.
    let step = solidcraft_trivariate::step_boundary::export_step(&model)
        .map_err(|e| invalid("boundary", e))?;
    if full {
        cmd(
            "solid.remove",
            json!({"bodies":[target],"name":"Replace extrusion stock"}),
        )?;
    }
    let inserted: Value =
        serde_json::from_str(&insert_step_content(&step)?).map_err(|e| invalid("import", e))?;
    let imported = inserted["feature_id"]
        .as_u64()
        .ok_or_else(|| invalid("import", "missing result feature"))?;
    let body = model_body_for(imported)?;
    let mut final_body = body.clone();
    if !full && operation != "new" {
        let reply = cmd(
            "solid.combine",
            json!({"target":target,"tools":[body],"operation":if operation=="cut"{"cut"}else{"join"},
            "keep_tools":false,"name":format!("Trivariate face {operation}")}),
        )?;
        let fid = reply["feature"]
            .as_u64()
            .ok_or_else(|| invalid("boolean", "missing combined feature"))?;
        final_body = model_body_for(fid)?;
    }
    let final_volume = volume(&final_body)?;
    if operation != "new"
        && ((operation == "cut" && final_volume >= initial_volume)
            || (operation == "add" && final_volume <= initial_volume))
    {
        return Err(invalid(
            "boolean",
            "the extrusion did not change material in the requested direction",
        ));
    }
    let scope = if full || operation == "new" {
        "result"
    } else {
        "tool"
    };
    let file = format!("face-{id}-{operation}-{imported}.{scope}.model.json");
    files::write_file(GRANT, &file, model.to_string().as_bytes())?;
    Ok(
        json!({"sketch":id,"operation":operation,"trivariate_scope":scope,"stock_partition_diagnostic":stock_diagnostic,"model":file,
              "body":final_body,"volume_mm3":final_volume})
        .to_string(),
    )
}
