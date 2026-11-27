wit_bindgen::generate!({ path: "wit", world: "command-plugin" });

use craft::plugin::types::{FieldFailure, PluginError};
use craft::plugin::{commands, files};
use exports::craft::plugin::entry::Guest;
use serde_json::{Value, json};

const GRANT: &str = "companion";
const SNAPSHOT_NAME: &str = "timeline.json";
const DESCRIPTOR_SCHEMA: &str = "solidcraft_sketch_preview_v1";
const MAX_DESCRIPTOR_BYTES: usize = 1_048_576;

mod ribbon;

struct Plugin;

impl Guest for Plugin {
    fn execute(entry: String, params: String) -> Result<String, PluginError> {
        let params: Value = serde_json::from_str(&params).map_err(|e| invalid("params", e))?;
        match entry.as_str() {
            "drill" | "boolean" | "chamfer" | "fillet" => ribbon::apply(&entry, &params),
            "snapshot-timeline" => snapshot_timeline(),
            "show-sketch-stage" => set_sketch_stage(&params, true),
            "hide-sketch-stage" => set_sketch_stage(&params, false),
            "insert-result" => insert_result(),
            "extrude" | "extrude-active-sketch" => extrude_active_sketch(&params),
            "extrude-active-timeline" => extrude_active_timeline(),
            "build-recipe" => build_recipe(&params),
            _ => Err(invalid("entry", format!("unknown entry `{entry}`"))),
        }
    }
}

fn export_timeline() -> Result<Value, PluginError> {
    let raw = commands::run("timeline.export_script", "{}")?;
    let export: Value = serde_json::from_str(&raw).map_err(|e| invalid("timeline", e))?;
    if export
        .get("omitted")
        .and_then(Value::as_array)
        .is_some_and(|items| !items.is_empty())
    {
        return Err(invalid(
            "timeline",
            "the native companion requires a complete timeline snapshot; SolidCraft reports omitted design content",
        ));
    }
    if export.get("script").is_none() {
        return Err(invalid("timeline", "host export did not return a script"));
    }
    Ok(export)
}

fn snapshot_timeline() -> Result<String, PluginError> {
    let export = export_timeline()?;
    let script = export
        .get("script")
        .ok_or_else(|| invalid("timeline", "host export did not return a script"))?;
    let output = json!({
        "schema": "solidcraft_trivariate_timeline_v1",
        "source": "SolidCraft timeline.export_script",
        "steps": export.get("steps").cloned().unwrap_or(Value::Null),
        "omitted": export.get("omitted").cloned().unwrap_or_else(|| json!([])),
        "script": script,
    });
    let bytes = serde_json::to_vec_pretty(&output).map_err(|e| invalid("timeline", e))?;
    files::write_file(GRANT, SNAPSHOT_NAME, &bytes)?;
    Ok(json!({"path": SNAPSHOT_NAME, "bytes": bytes.len(), "steps": output["steps"], "omitted": output["omitted"]}).to_string())
}

fn set_sketch_stage(params: &Value, visible: bool) -> Result<String, PluginError> {
    let path = params
        .get("descriptor")
        .and_then(Value::as_str)
        .unwrap_or("partition.json");
    validate_leaf_name(path)?;
    let requested_source_id = params
        .get("source_feature_id")
        .and_then(Value::as_u64)
        .unwrap_or(1);
    let (mut descriptor, generated) = match files::read_file(GRANT, path) {
        Ok(bytes) => {
            if bytes.len() > MAX_DESCRIPTOR_BYTES {
                return Err(invalid("descriptor", "descriptor exceeds 1 MiB"));
            }
            (
                serde_json::from_slice(&bytes).map_err(|e| invalid("descriptor", e))?,
                false,
            )
        }
        Err(PluginError::InvalidInput(f)) if f.field == "path" && f.reason == "no such file" => {
            if !visible {
                return Ok(json!({"visible": false, "installed": false}).to_string());
            }
            (make_native_descriptor(requested_source_id)?, true)
        }
        Err(error) => return Err(error),
    };
    if descriptor.get("schema").and_then(Value::as_str) != Some(DESCRIPTOR_SCHEMA) {
        return Err(invalid(
            "descriptor",
            format!("expected schema `{DESCRIPTOR_SCHEMA}`"),
        ));
    }
    let source_id = descriptor
        .get("source_feature_id")
        .and_then(Value::as_u64)
        .ok_or_else(|| {
            invalid(
                "descriptor",
                "source_feature_id must be an unsigned feature id",
            )
        })?;
    let source_feature = descriptor
        .get("source_feature")
        .ok_or_else(|| invalid("descriptor", "source_feature is missing"))?;
    let source_script = descriptor
        .get("source_script")
        .ok_or_else(|| invalid("descriptor", "source_script is missing"))?;
    let preview_feature = descriptor
        .get("feature")
        .and_then(Value::as_object)
        .ok_or_else(|| {
            invalid(
                "descriptor",
                "feature must be a native sketch feature object",
            )
        })?;
    if preview_feature.get("type").and_then(Value::as_str) != Some("sketch") {
        return Err(invalid(
            "descriptor",
            "preview feature must have type `sketch`",
        ));
    }
    validate_construction_wire(preview_feature)?;
    if source_feature.get("id").and_then(Value::as_u64) != Some(source_id) {
        return Err(invalid(
            "descriptor",
            "source_feature id does not match source_feature_id",
        ));
    }
    if preview_feature.get("id").and_then(Value::as_u64).is_none()
        || preview_feature.get("name").and_then(Value::as_str)
            != Some(format!("SCT Partition {source_id}").as_str())
    {
        return Err(invalid(
            "descriptor",
            "preview sketch must have the native name `SCT Partition <source_feature_id>` and an assigned feature id",
        ));
    }

    let current_export = export_timeline()?;
    let mut current_script = current_export
        .get("script")
        .cloned()
        .ok_or_else(|| invalid("timeline", "host export did not return a script"))?;
    let installed_id = descriptor
        .get("installed_feature_id")
        .and_then(Value::as_u64);
    if let Some(id) = installed_id {
        if descriptor["feature"].get("id").and_then(Value::as_u64) != Some(id)
            || !script_has_feature(&current_script, id, &descriptor["feature"])
        {
            return Err(invalid(
                "descriptor",
                "the installed preview sketch no longer matches its descriptor",
            ));
        }
        strip_installed_feature(&mut current_script, id);
    }
    if &current_script != source_script {
        return Err(invalid(
            "descriptor",
            "descriptor source does not match the active SolidCraft timeline; refresh the snapshot and rebuild the preview",
        ));
    }
    if !script_has_feature(source_script, source_id, source_feature) {
        return Err(invalid(
            "descriptor",
            "source feature is not present in the recorded source timeline",
        ));
    }
    let descriptor_source = solidcraft_trivariate::json::parse(&source_script.to_string())
        .map_err(|e| invalid("descriptor", e))?;
    let expected_native = solidcraft_trivariate::solidcraft_preview::descriptor(
        &descriptor_source,
        u32::try_from(source_id)
            .map_err(|_| invalid("descriptor", "source feature id is too large"))?,
    )
    .map_err(|e| invalid("descriptor", e))?;
    let expected_native: Value =
        serde_json::from_str(&expected_native.to_string()).map_err(|e| invalid("descriptor", e))?;
    if expected_native.get("feature") != descriptor.get("feature")
        || expected_native.get("internal_edges") != descriptor.get("internal_edges")
    {
        return Err(invalid(
            "descriptor",
            "preview construction wires do not match the native reconstruction of the recorded source sketch",
        ));
    }

    let id = match installed_id {
        Some(id) => id,
        None if visible => {
            let reply = commands::run(
                "timeline.add",
                &json!({"feature": descriptor["feature"]}).to_string(),
            )?;
            let result: Value =
                serde_json::from_str(&reply).map_err(|e| invalid("timeline.add", e))?;
            let id = result
                .get("feature")
                .and_then(Value::as_u64)
                .ok_or_else(|| {
                    invalid(
                        "timeline.add",
                        "host did not return the new sketch feature id",
                    )
                })?;
            if descriptor["feature"].get("id").and_then(Value::as_u64) != Some(id) {
                return Err(invalid(
                    "timeline.add",
                    "host assigned a different feature id than the companion descriptor",
                ));
            }
            descriptor["installed_feature_id"] = json!(id);
            id
        }
        None => return Ok(json!({"visible": false, "installed": false}).to_string()),
    };

    commands::run(
        "browser.visibility",
        &json!({"visible": visible, "items": [id]}).to_string(),
    )?;
    if visible && (installed_id.is_none() || generated) {
        let updated =
            serde_json::to_vec_pretty(&descriptor).map_err(|e| invalid("descriptor", e))?;
        files::write_file(GRANT, path, &updated)?;
    }
    Ok(json!({"visible": visible, "installed": true, "sketch_feature_id": id, "source_feature_id": source_id}).to_string())
}

fn insert_result() -> Result<String, PluginError> {
    let bytes = files::read_file(GRANT, "result.step")?;
    let step = std::str::from_utf8(&bytes).map_err(|e| invalid("result.step", e))?;
    insert_step_content(step)
}

fn insert_step_content(step: &str) -> Result<String, PluginError> {
    let export = export_timeline()?;
    let command_list = export
        .pointer("/script/commands")
        .and_then(Value::as_array)
        .ok_or_else(|| invalid("timeline", "host export did not return script commands"))?;
    let next_id = command_list
        .iter()
        .filter_map(|command| {
            command
                .pointer("/params/feature/id")
                .and_then(Value::as_u64)
        })
        .max()
        .unwrap_or(0)
        .checked_add(1)
        .ok_or_else(|| invalid("timeline", "feature id space is exhausted"))?;
    let feature = json!({
        "type": "import",
        "name": "Trivariate result",
        "file": "result.step",
        "step": step,
        "id": next_id,
        "suppressed": false
    });
    let reply = commands::run("timeline.add", &json!({"feature": feature}).to_string())?;
    let result: Value = serde_json::from_str(&reply).map_err(|e| invalid("timeline.add", e))?;
    if result.get("feature").and_then(Value::as_u64) != Some(next_id) {
        return Err(invalid(
            "timeline.add",
            "host did not confirm the imported result feature id",
        ));
    }
    commands::run("view.visual_style", &json!({"style": "shaded"}).to_string())?;
    Ok(json!({"feature_id": next_id, "name": "Trivariate result"}).to_string())
}

fn make_native_descriptor(source_id: u64) -> Result<Value, PluginError> {
    let mut export = export_timeline()?;
    let source_id_u32 = u32::try_from(source_id)
        .map_err(|_| invalid("source_feature_id", "feature id is too large"))?;
    let mut existing = None;
    if let Some(steps) = export
        .pointer_mut("/script/commands")
        .and_then(Value::as_array_mut)
    {
        let matches = steps
            .iter()
            .enumerate()
            .filter(|(_, step)| {
                step.get("command").and_then(Value::as_str) == Some("timeline.add")
                    && step.pointer("/params/feature/name").and_then(Value::as_str)
                        == Some(format!("SCT Partition {source_id}").as_str())
            })
            .map(|(index, step)| (index, step.pointer("/params/feature").cloned()))
            .collect::<Vec<_>>();
        if matches.len() > 1 {
            return Err(invalid(
                "descriptor",
                "multiple generated partition sketches are present",
            ));
        }
        if let Some((index, Some(feature))) = matches.into_iter().next() {
            let id = feature.get("id").and_then(Value::as_u64).ok_or_else(|| {
                invalid("descriptor", "generated partition sketch has no feature id")
            })?;
            steps.remove(index);
            existing = Some((id, feature));
        }
    }
    let input = solidcraft_trivariate::json::parse(&export.to_string())
        .map_err(|e| invalid("timeline", e))?;
    let native = solidcraft_trivariate::solidcraft_preview::descriptor(&input, source_id_u32)
        .map_err(|e| invalid("source_feature_id", e))?;
    let mut descriptor: Value =
        serde_json::from_str(&native.to_string()).map_err(|e| invalid("descriptor", e))?;
    if let Some((id, feature)) = existing {
        if descriptor.get("feature") != Some(&feature) {
            return Err(invalid(
                "descriptor",
                "existing partition sketch does not match the native preview reconstruction",
            ));
        }
        descriptor["installed_feature_id"] = json!(id);
    }
    Ok(descriptor)
}

fn extrude_active_sketch(params: &Value) -> Result<String, PluginError> {
    let source_id = params
        .get("source_feature_id")
        .and_then(Value::as_u64)
        .unwrap_or(1);
    let source_id = u32::try_from(source_id)
        .map_err(|_| invalid("source_feature_id", "feature id is too large"))?;
    let height = params
        .get("height_mm")
        .and_then(Value::as_i64)
        .unwrap_or(12);
    let height = i32::try_from(height)
        .map_err(|_| invalid("height_mm", "height must be an integer from 1 to 64 mm"))?;
    if !(1..=64).contains(&height) {
        return Err(invalid(
            "height_mm",
            "height must be an integer from 1 to 64 mm",
        ));
    }
    let export = export_timeline()?;
    let input = solidcraft_trivariate::json::parse(&export.to_string())
        .map_err(|e| invalid("timeline", e))?;
    let model =
        solidcraft_trivariate::solidcraft_polygon::extrude_source_sketch(&input, source_id, height)
            .map_err(|e| invalid("timeline", e))?;
    insert_audited_model(&model)
}

fn extrude_active_timeline() -> Result<String, PluginError> {
    let export = export_timeline()?;
    let input = solidcraft_trivariate::json::parse(&export.to_string())
        .map_err(|e| invalid("timeline", e))?;
    let model = solidcraft_trivariate::solidcraft_polygon::extrude_preview(&input)
        .map_err(|e| invalid("timeline", e))?;
    insert_audited_model(&model)
}

fn build_recipe(params: &Value) -> Result<String, PluginError> {
    let path = params
        .get("recipe")
        .and_then(Value::as_str)
        .unwrap_or("recipe.json");
    validate_leaf_name(path)?;
    let bytes = files::read_file(GRANT, path)?;
    if bytes.len() > MAX_DESCRIPTOR_BYTES {
        return Err(invalid("recipe", "recipe exceeds 1 MiB"));
    }
    let mut source: Value = serde_json::from_slice(&bytes).map_err(|e| invalid("recipe", e))?;
    let operation = params
        .get("operation")
        .and_then(Value::as_str)
        .unwrap_or("from_recipe");
    match operation {
        "from_recipe" => {}
        "union" | "intersection" | "a_minus_b" | "b_minus_a" => {
            if source.get("family").and_then(Value::as_str) != Some("box_cylinder_boolean_suite_v1")
            {
                return Err(invalid(
                    "operation",
                    "operation overrides require family `box_cylinder_boolean_suite_v1`; no feature was added",
                ));
            }
            source["operation"] = json!(operation);
        }
        _ => {
            return Err(invalid(
                "operation",
                "choose from_recipe, union, intersection, a_minus_b or b_minus_a",
            ));
        }
    }
    let recipe = solidcraft_trivariate::json::parse(&source.to_string())
        .map_err(|e| invalid("recipe", e))?;
    let model =
        solidcraft_trivariate::interchange::build(&recipe).map_err(|e| invalid("recipe", e))?;
    let rim = params.get("rim").and_then(Value::as_str).unwrap_or("sharp");
    let step = match rim {
        "sharp" => {
            audit_model(&model)?;
            solidcraft_trivariate::step_boundary::export_step(&model)
        }
        "chamfer" | "fillet" => {
            solidcraft_trivariate::fieldcraft_rims::construct(&model, rim)
                .map_err(|e| invalid("rim", e))?;
            solidcraft_trivariate::fieldcraft_rims::export_step(&model, rim)
        }
        _ => return Err(invalid("rim", "choose sharp, chamfer or fillet")),
    }
    .map_err(|e| invalid("model", e))?;
    files::write_file(GRANT, "result.step", step.as_bytes())?;
    let inserted = insert_step_content(&step)?;
    let mut result: Value =
        serde_json::from_str(&inserted).map_err(|e| invalid("timeline.add", e))?;
    result["rim"] = json!(rim);
    Ok(result.to_string())
}

fn insert_audited_model(model: &solidcraft_trivariate::json::Json) -> Result<String, PluginError> {
    audit_model(model)?;
    let step = solidcraft_trivariate::step_boundary::export_step(model)
        .map_err(|e| invalid("model", e))?;
    files::write_file(GRANT, "result.step", step.as_bytes())?;
    insert_step_content(&step)
}

fn audit_model(model: &solidcraft_trivariate::json::Json) -> Result<(), PluginError> {
    let audit =
        solidcraft_trivariate::interchange::audit(model).map_err(|e| invalid("model", e))?;
    if audit
        .field("status")
        .and_then(solidcraft_trivariate::json::Json::text)
        != Ok("valid")
    {
        return Err(invalid(
            "model",
            "native trivariate audit did not return valid",
        ));
    }
    Ok(())
}

fn validate_construction_wire(feature: &serde_json::Map<String, Value>) -> Result<(), PluginError> {
    if feature.get("plane").is_none() {
        return Err(invalid(
            "descriptor",
            "preview sketch must declare its plane",
        ));
    }
    let curves = feature
        .get("sketch")
        .and_then(|sketch| sketch.get("curves"))
        .and_then(Value::as_array)
        .ok_or_else(|| invalid("descriptor", "preview sketch must contain sketch.curves"))?;
    if curves.is_empty()
        || curves.iter().any(|curve| {
            curve.get("type").and_then(Value::as_str) != Some("line")
                || curve.get("construction").and_then(Value::as_bool) != Some(true)
        })
    {
        return Err(invalid(
            "descriptor",
            "preview geometry must contain only construction lines",
        ));
    }
    Ok(())
}

fn strip_installed_feature(script: &mut Value, id: u64) {
    let Some(commands) = script.get_mut("commands").and_then(Value::as_array_mut) else {
        return;
    };
    commands.retain(|step| {
        let is_preview = step.get("command").and_then(Value::as_str) == Some("timeline.add")
            && step.pointer("/params/feature/id").and_then(Value::as_u64) == Some(id);
        !is_preview
    });
}

fn script_has_feature(script: &Value, id: u64, expected: &Value) -> bool {
    script
        .get("commands")
        .and_then(Value::as_array)
        .is_some_and(|steps| {
            steps.iter().any(|step| {
                step.get("command").and_then(Value::as_str) == Some("timeline.add")
                    && step.pointer("/params/feature/id").and_then(Value::as_u64) == Some(id)
                    && step.pointer("/params/feature") == Some(expected)
            })
        })
}

fn validate_leaf_name(path: &str) -> Result<(), PluginError> {
    let valid = !path.is_empty()
        && path.len() <= 128
        && path != "."
        && path != ".."
        && !path.starts_with('.')
        && path
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-'));
    if valid {
        Ok(())
    } else {
        Err(invalid(
            "descriptor",
            "use a simple filename without folders or path traversal",
        ))
    }
}

fn invalid(field: &str, reason: impl ToString) -> PluginError {
    PluginError::InvalidInput(FieldFailure {
        field: field.into(),
        reason: reason.to_string(),
    })
}

export!(Plugin);
