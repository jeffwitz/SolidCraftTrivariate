//! Queue a trusted desktop-side IGA bridge request without granting Wasm process execution.
use super::*;

const REQUEST: &str = "fieldcraft-analysis-request.json";
const STATUS: &str = "fieldcraft-analysis-status.json";
const MAX_INPUT_BYTES: usize = 1_048_576;

fn read_small_json(path: &str, field: &str) -> Result<Value, PluginError> {
    validate_leaf_name(path)?;
    let bytes = files::read_file(GRANT, path)?;
    if bytes.len() > MAX_INPUT_BYTES {
        return Err(invalid(field, "input exceeds 1 MiB"));
    }
    serde_json::from_slice(&bytes).map_err(|e| invalid(field, e))
}

pub(super) fn queue(params: &Value) -> Result<String, PluginError> {
    let recipe = params
        .get("recipe")
        .and_then(Value::as_str)
        .unwrap_or("fieldcraft_x64_y64.recipe.json");
    let job = params
        .get("job")
        .and_then(Value::as_str)
        .unwrap_or("fieldcraft_bracket_hex_pa5.json");
    let case_name = params
        .get("case_name")
        .and_then(Value::as_str)
        .unwrap_or("solidcraft-analysis");
    validate_leaf_name(recipe)?;
    validate_leaf_name(job)?;
    if case_name.is_empty()
        || case_name.len() > 64
        || !case_name
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'-' || b == b'_')
    {
        return Err(invalid(
            "case_name",
            "use 1-64 ASCII letters, digits, '-' or '_'",
        ));
    }
    let status = match files::read_file(GRANT, STATUS) {
        Ok(bytes) if bytes.len() <= 65_536 => {
            serde_json::from_slice::<Value>(&bytes).map_err(|e| invalid("status", e))?
        }
        Ok(_) => return Err(invalid("status", "status record exceeds 64 KiB")),
        Err(PluginError::InvalidInput(f)) if f.field == "path" && f.reason == "no such file" => {
            return Err(invalid("FieldCraft", "desktop bridge is not active; launch the compact SolidCraft desktop with FieldCraft installed"));
        }
        Err(error) => return Err(error),
    };
    if status.get("available").and_then(Value::as_bool) != Some(true)
        || !matches!(
            status.get("state").and_then(Value::as_str),
            Some("ready" | "complete" | "failed")
        )
    {
        return Err(invalid(
            "FieldCraft",
            status.get("message").and_then(Value::as_str)
                .unwrap_or("FieldCraft solve and display integration is unavailable; CAD editing remains available"),
        ));
    }
    let recipe_value = read_small_json(recipe, "recipe")?;
    let job_value = read_small_json(job, "job")?;
    if recipe_value.get("family").and_then(Value::as_str).is_none() {
        return Err(invalid(
            "recipe",
            "expected a versioned recipe with a family identifier",
        ));
    }
    if job_value.as_object().is_none() {
        return Err(invalid("job", "expected a JSON object"));
    }
    if files::read_file(GRANT, REQUEST).is_ok() {
        return Err(invalid(
            "request",
            "another FieldCraft analysis request is pending",
        ));
    }
    let queued = json!({
        "schema": "solidcraft_fieldcraft_analysis_request_v1",
        "recipe_source": recipe,
        "job_source": job,
        "recipe": recipe_value,
        "job": job_value,
        "case_name": case_name,
    });
    let bytes = serde_json::to_vec(&queued).map_err(|e| invalid("request", e))?;
    if bytes.len() > MAX_INPUT_BYTES {
        return Err(invalid(
            "request",
            "combined recipe and job snapshots exceed 1 MiB",
        ));
    }
    files::write_file(GRANT, REQUEST, &bytes)?;
    Ok(json!({
        "queued": true,
        "case_name": case_name,
        "state": "queued",
        "provenance": "SolidCraft native IGA result; FieldCraft archive is a visualization carrier, not a FieldCraft solve",
    }).to_string())
}

pub(super) fn status() -> Result<String, PluginError> {
    match files::read_file(GRANT, STATUS) {
        Ok(bytes) if bytes.len() <= 65_536 => {
            let value: Value = serde_json::from_slice(&bytes)
                .map_err(|e| invalid("status", e))?;
            Ok(value.to_string())
        }
        Ok(_) => Err(invalid("status", "status record exceeds 64 KiB")),
        Err(PluginError::InvalidInput(f)) if f.field == "path" && f.reason == "no such file" => {
            Ok(json!({
                "available": false,
                "state": "unavailable",
                "message": "Run the SolidCraft compact desktop launcher with FieldCraft installed to enable solve and display. CAD editing remains available without FieldCraft."
            }).to_string())
        }
        Err(error) => Err(error),
    }
}
