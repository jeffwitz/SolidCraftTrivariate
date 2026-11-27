//! Solid editing commands for the compact sketch ribbon.
use super::*;
type Model = solidcraft_trivariate::json::Json;

fn build(source: &Value) -> Result<Model, PluginError> {
    let recipe = solidcraft_trivariate::json::parse(&source.to_string())
        .map_err(|e| invalid("recipe", e))?;
    let model =
        solidcraft_trivariate::interchange::build(&recipe).map_err(|e| invalid("recipe", e))?;
    audit_model(&model)?;
    Ok(model)
}

fn boundary(model: &Model) -> Result<String, PluginError> {
    solidcraft_trivariate::step_boundary::export_step(model).map_err(|e| invalid("body", e))
}

fn feature_id(params: &Value, key: &str, default: u64) -> Result<u64, PluginError> {
    params
        .get(key)
        .unwrap_or(&json!(default))
        .as_u64()
        .filter(|n| *n > 0)
        .ok_or_else(|| invalid(key, "select a positive timeline feature id"))
}

// STEP exchange timestamps may differ between a desktop build and Wasm.
// Compare every DATA entity including incidences and coefficients.
fn verify_operand(timeline: &Value, id: u64, expected: &str) -> Result<(), PluginError> {
    let features = timeline
        .pointer("/script/commands")
        .and_then(Value::as_array)
        .ok_or_else(|| invalid("body", "missing design timeline"))?;
    let actual = features
        .iter()
        .filter_map(|c| c.pointer("/params/feature"))
        .find(|f| f["id"].as_u64() == Some(id) && f["type"] == "import" && f["suppressed"] == false)
        .and_then(|f| f["step"].as_str());
    let data = |s: &str| s.find("\nDATA;\n");
    if let Some((a, b)) = actual.and_then(|a| Some((&a[data(a)?..], &expected[data(expected)?..])))
    {
        if a == b {
            return Ok(());
        }
    }
    Err(invalid(
        "body",
        "the selected feature does not match the qualified trivariate operand",
    ))
}

fn replace(id: u64, step: String) -> Result<String, PluginError> {
    commands::run(
        "timeline.edit",
        &json!({"feature":id,"set":{"step":step}}).to_string(),
    )?;
    commands::run("view.visual_style", &json!({"style":"shaded"}).to_string())?;
    Ok(json!({"feature_id":id,"updated":true}).to_string())
}

fn suite() -> Value {
    serde_json::from_str(include_str!("../examples/box-cylinder-tall.recipe.json"))
        .expect("reference fixture")
}

pub(super) fn apply(entry: &str, params: &Value) -> Result<String, PluginError> {
    let id = feature_id(params, "body_id", 1)?;
    let timeline = export_timeline()?;
    if matches!(entry, "chamfer-bores" | "fillet-bores") {
        let recipe =
            serde_json::from_str(include_str!("../examples/fieldcraft_x64_y64.recipe.json"))
                .expect("reference fixture");
        let model = build(&recipe)?;
        verify_operand(&timeline, id, &boundary(&model)?)?;
        let state = if entry == "chamfer-bores" {
            "chamfer"
        } else {
            "fillet"
        };
        solidcraft_trivariate::fieldcraft_rims::construct(&model, state)
            .map_err(|e| invalid("rim", e))?;
        return replace(
            id,
            solidcraft_trivariate::fieldcraft_rims::export_step(&model, state)
                .map_err(|e| invalid("rim", e))?,
        );
    }

    let mut recipe = suite();
    recipe["operation"] = json!("box");
    let from_sketch = params
        .get("source_feature_id")
        .and_then(Value::as_u64)
        .unwrap_or(0);
    if entry == "drill" && from_sketch != 0 {
        let height = params["height_mm"].as_i64().unwrap_or(12);
        let input = solidcraft_trivariate::json::parse(&timeline.to_string())
            .map_err(|e| invalid("timeline", e))?;
        let model = solidcraft_trivariate::solidcraft_polygon::extrude_source_sketch(
            &input,
            u32::try_from(from_sketch).map_err(|_| invalid("sketch", "feature id too large"))?,
            i32::try_from(height).map_err(|_| invalid("height_mm", "height out of range"))?,
        )
        .map_err(|e| invalid("sketch", e))?;
        verify_operand(&timeline, id, &boundary(&model)?)?;
        let geometry: Value =
            serde_json::from_str(&model.to_string()).map_err(|e| invalid("body", e))?;
        let (lo, hi) = rectangular_bounds(&geometry).ok_or_else(|| {
            invalid(
                "body",
                "drilling requires one axis-aligned rectangular extrusion",
            )
        })?;
        recipe["stock"]["min_mm"] = json!(lo);
        recipe["stock"]["max_mm"] = json!(hi);
        recipe["tool"]["center_mm"] = json!([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2]);
        recipe["tool"]["interval_mm"] = json!([lo[2] - 1, hi[2] + 1]);
    } else {
        verify_operand(&timeline, id, &boundary(&build(&recipe)?)?)?;
    }
    let mut tool = None;
    let operation = if entry == "boolean" {
        let tool_id = feature_id(params, "tool_id", 2)?;
        recipe["operation"] = json!("cylinder");
        verify_operand(&timeline, tool_id, &boundary(&build(&recipe)?)?)?;
        tool = Some(tool_id);
        match params
            .get("operation")
            .and_then(Value::as_str)
            .unwrap_or("union")
        {
            op @ ("union" | "intersection" | "a_minus_b" | "b_minus_a") => op,
            _ => return Err(invalid("operation", "unknown Boolean operation")),
        }
    } else {
        recipe["tool"]["radius_mm"] = params.get("radius_mm").cloned().unwrap_or(json!(5));
        "a_minus_b"
    };
    recipe["operation"] = json!(operation);
    let step = boundary(&build(&recipe)?)?;
    let reply = replace(id, step)?;
    if let Some(tool_id) = tool {
        commands::run(
            "timeline.suppress",
            &json!({"feature":tool_id,"suppressed":true}).to_string(),
        )?;
    }
    Ok(reply)
}

fn rectangular_bounds(model: &Value) -> Option<([i64; 3], [i64; 3])> {
    let [block] = model["blocks"].as_array()?.as_slice() else {
        return None;
    };
    if block["degrees"] != json!([1, 1, 1]) {
        return None;
    }
    let mut corners = std::collections::BTreeSet::new();
    for control in block["controls"].as_array()? {
        let [x, y, z, w] = control.as_array()?.as_slice() else {
            return None;
        };
        let w = w.as_i64()?;
        let xyz = [x.as_i64()?, y.as_i64()?, z.as_i64()?];
        if w <= 0 || xyz.iter().any(|v| v % w != 0) {
            return None;
        }
        corners.insert(xyz.map(|v| v / w));
    }
    let lo = std::array::from_fn(|i| corners.iter().map(|p| p[i]).min().unwrap_or(0));
    let hi = std::array::from_fn(|i| corners.iter().map(|p| p[i]).max().unwrap_or(0));
    if corners.len() != 8
        || corners
            .iter()
            .any(|p| (0..3).any(|i| p[i] != lo[i] && p[i] != hi[i]))
    {
        return None;
    }
    Some((lo, hi))
}
