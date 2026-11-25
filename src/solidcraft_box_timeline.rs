//! Translate the canonical solidcraft/66 box feature timeline.
//!
//! Only literal axis-aligned boxes with a single initial body followed by
//! ordered joins and cuts are accepted. Feature identities, parameter names,
//! suppression state and center/corner placement are checked before the
//! history is converted to the direct `solid.box` form audited by SC3.

use crate::{
    json::{object, Json},
    solidcraft_boxes,
};

fn number(value: i128) -> Json {
    Json::Number(value.to_string())
}

fn exact_integer(value: &Json) -> Result<i128, &'static str> {
    let Json::String(value) = value else {
        return Err("box timeline dimensions must be resolved string literals");
    };
    let (numerator, denominator) = solidcraft_boxes::decimal_rational(value)?;
    if numerator % denominator != 0 {
        return Err("box timeline dimensions must be integer millimetres");
    }
    Ok(numerator / denominator)
}

fn point(value: &Json) -> Result<(), &'static str> {
    let values = value.items()?;
    if values.len() != 3 {
        return Err("box timeline placement needs three coordinates");
    }
    for item in values {
        let Json::Number(token) = item else {
            return Err("box timeline placement must contain numeric coordinates");
        };
        let (numerator, denominator) = solidcraft_boxes::decimal_rational(token)?;
        if numerator
            .checked_mul(2)
            .ok_or("box timeline coordinate range")?
            % denominator
            != 0
        {
            return Err("box timeline placement is outside the exact half-mm grid");
        }
    }
    Ok(())
}

fn feature(command: &Json, expected_id: usize, base: bool) -> Result<&Json, &'static str> {
    command.keys(&["command", "params"])?;
    if command.field("command")?.text()? != "timeline.add" {
        return Err("box timeline contains an unsupported operation");
    }
    let params = command.field("params")?;
    params.keys(&["feature"])?;
    let feature = params.field("feature")?;
    let mut keys = vec![
        "corner",
        "height",
        "id",
        "length",
        "name",
        "operation",
        "param_names",
        "suppressed",
        "type",
        "width",
    ];
    if base {
        keys.push("body_names");
    }
    let has_centred = feature.field("centred").is_ok();
    if has_centred {
        keys.push("centred");
    }
    feature.keys(&keys)?;
    if feature.field("type")?.text()? != "box" {
        return Err("box timeline feature type mismatch");
    }
    if feature.field("id")?.number::<usize>()? != expected_id {
        return Err("box timeline feature ids must be consecutive");
    }
    if feature.field("suppressed")? != &Json::Bool(false) {
        return Err("suppressed box timeline features are not supported");
    }
    let name = feature.field("name")?.text()?;
    if name.is_empty()
        || name.len() > 48
        || !name
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || byte == b'_')
    {
        return Err("box timeline feature name is not a stable native id");
    }
    if base && name != "base" {
        return Err("the initial box timeline feature must be named base");
    }
    let operation = feature.field("operation")?.text()?;
    if (base && operation != "new_body") || (!base && !matches!(operation, "join" | "cut")) {
        return Err("box timeline operation is outside the base/join/cut subset");
    }
    let names = feature.field("param_names")?.items()?;
    if names.len() != 3 {
        return Err("box timeline feature needs three resolved dimension parameters");
    }
    for (axis, name) in names.iter().enumerate() {
        let expected = format!("d{}", (expected_id - 1) * 3 + axis + 1);
        if name.text()? != expected {
            return Err("box timeline parameter references are missing, duplicated or reordered");
        }
    }
    for key in ["length", "width", "height"] {
        if exact_integer(feature.field(key)?)? <= 0 {
            return Err("box timeline dimensions must be positive");
        }
    }
    point(feature.field("corner")?)?;
    if has_centred
        && feature.field("centred")? != &Json::Bool(true)
        && feature.field("centred")? != &Json::Bool(false)
    {
        return Err("box timeline centered flag must be boolean");
    }
    if base {
        let body_names = feature.field("body_names")?.items()?;
        if body_names.len() != 1 || body_names[0].text()?.is_empty() {
            return Err("base box must create exactly one named body");
        }
    }
    Ok(feature)
}

/// Convert a canonical box timeline into the native audited orthogonal model.
pub fn import_box_timeline(script: &Json) -> Result<Json, &'static str> {
    script.keys(&["about", "commands"])?;
    let about = script.field("about")?;
    about.keys(&["design", "format", "replay", "version"])?;
    if !matches!(about.field("design")?.text()?, "editable" | "original")
        || about.field("format")?.text()? != "solidcraft/66"
        || about.field("version")?.text()? != "0.48.1"
        || about.field("replay")?.text()?.is_empty()
    {
        return Err("unsupported SolidCraft box timeline format or version");
    }
    let commands = script.field("commands")?.items()?;
    if commands.is_empty() || commands.len() > 9 {
        return Err("box timeline requires one base and at most eight operations");
    }
    let mut names = std::collections::BTreeSet::new();
    let mut direct_commands = Vec::with_capacity(commands.len());
    for (index, command) in commands.iter().enumerate() {
        let base = index == 0;
        let feature = feature(command, index + 1, base)?;
        let name = feature.field("name")?.text()?.to_owned();
        if !names.insert(name.clone()) {
            return Err("box timeline feature names must be unique");
        }
        let mut params = vec![
            (
                "length".to_owned(),
                number(exact_integer(feature.field("length")?)?),
            ),
            (
                "width".to_owned(),
                number(exact_integer(feature.field("width")?)?),
            ),
            (
                "height".to_owned(),
                number(exact_integer(feature.field("height")?)?),
            ),
            ("name".to_owned(), Json::String(name)),
            (
                "operation".to_owned(),
                Json::String(
                    match feature.field("operation")?.text()? {
                        "new_body" => "new",
                        "join" => "join",
                        "cut" => "cut",
                        _ => return Err("unsupported box timeline operation"),
                    }
                    .into(),
                ),
            ),
        ];
        let centered = feature
            .field("centred")
            .is_ok_and(|value| value == &Json::Bool(true));
        params.push((
            if centered { "center" } else { "corner" }.into(),
            feature.field("corner")?.clone(),
        ));
        direct_commands.push(object([
            ("command", Json::String("solid.box".into())),
            ("params", Json::Object(params.into_iter().collect())),
        ]));
    }
    solidcraft_boxes::import_boxes(&object([("commands", Json::Array(direct_commands))]))
}
