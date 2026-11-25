//! Import the qualified SolidCraft box plus simple blind-hole subset through G21.
use crate::{
    interval_bores_json,
    json::{object, Json},
    solidcraft_boxes::decimal_rational,
};

#[derive(Clone, Copy)]
struct BoxGeometry {
    lo: [i32; 3],
    size: [i32; 3],
}

#[derive(Clone, Copy)]
struct HoleGeometry {
    center: [i32; 2],
    diameter: i32,
    depth: i32,
    face_z: i32,
    direction_z: i32,
}

fn text(value: &Json) -> Result<&str, &'static str> {
    value.text()
}

fn numeric_literal(value: &Json, allow_string: bool) -> Result<(i128, i128), &'static str> {
    match value {
        Json::Number(token) => decimal_rational(token),
        Json::String(token) if allow_string => decimal_rational(token),
        _ => Err("SC8 values must be literal numbers"),
    }
}

fn integer(value: &Json, allow_string: bool) -> Result<i32, &'static str> {
    let (numerator, denominator) = numeric_literal(value, allow_string)?;
    if numerator % denominator != 0 {
        return Err("SC8 geometry must use integer millimetres");
    }
    i32::try_from(numerator / denominator).map_err(|_| "SC8 integer range")
}

fn vector3(value: &Json, allow_string: bool) -> Result<[i32; 3], &'static str> {
    let values = value.items()?;
    if values.len() != 3 {
        return Err("SC8 point must have three coordinates");
    }
    Ok([
        integer(&values[0], allow_string)?,
        integer(&values[1], allow_string)?,
        integer(&values[2], allow_string)?,
    ])
}

fn dimensions(length: i32, width: i32, height: i32) -> Result<[i32; 3], &'static str> {
    if !(1..=64).contains(&length) || !(1..=64).contains(&width) || !(1..=64).contains(&height) {
        return Err("SC8 box dimensions outside native range");
    }
    Ok([length, width, height])
}

fn direct_box(feature: &Json) -> Result<BoxGeometry, &'static str> {
    feature.keys(&["command", "params"])?;
    if text(feature.field("command")?)? != "solid.box" {
        return Err("SC8 requires solid.box as the first direct command");
    }
    let params = feature.field("params")?;
    params.keys(&["corner", "height", "length", "name", "operation", "width"])?;
    if text(params.field("operation")?)? != "new" {
        return Err("SC8 requires a new base box");
    }
    let _ = text(params.field("name")?)?;
    let lo = vector3(params.field("corner")?, false)?;
    let size = dimensions(
        integer(params.field("length")?, false)?,
        integer(params.field("width")?, false)?,
        integer(params.field("height")?, false)?,
    )?;
    if lo[2] != 0
        || !(-32..=32).contains(&lo[0])
        || !(-32..=32).contains(&lo[1])
        || lo[0].checked_add(size[0]).is_none_or(|x| x > 64)
        || lo[1].checked_add(size[1]).is_none_or(|y| y > 64)
    {
        return Err("SC8 base box must be an origin-based bounded XY prism");
    }
    Ok(BoxGeometry { lo, size })
}

fn direct_hole(feature: &Json) -> Result<HoleGeometry, &'static str> {
    feature.keys(&["command", "params"])?;
    if text(feature.field("command")?)? != "solid.hole" {
        return Err("SC8 requires one solid.hole after the base box");
    }
    let params = feature.field("params")?;
    params.keys(&["depth", "diameter", "name", "position", "type"])?;
    if text(params.field("type")?)? != "simple" {
        return Err("SC8 supports simple flat-bottom holes only");
    }
    let _ = text(params.field("name")?)?;
    let position = vector3(params.field("position")?, false)?;
    Ok(HoleGeometry {
        center: [position[0], position[1]],
        diameter: integer(params.field("diameter")?, false)?,
        depth: integer(params.field("depth")?, false)?,
        face_z: position[2],
        direction_z: 0,
    })
}

fn timeline_feature<'a>(command: &'a Json, expected: &str) -> Result<&'a Json, &'static str> {
    command.keys(&["command", "params"])?;
    if text(command.field("command")?)? != "timeline.add" {
        return Err("SC8 accepts timeline.add features only");
    }
    let params = command.field("params")?;
    params.keys(&["feature"])?;
    let feature = params.field("feature")?;
    if text(feature.field("type")?)? != expected {
        return Err("SC8 feature sequence type mismatch");
    }
    Ok(feature)
}

fn timeline_box(feature: &Json) -> Result<BoxGeometry, &'static str> {
    feature.keys(&[
        "body_names",
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
    ])?;
    if integer(feature.field("id")?, false)? != 1
        || feature.field("suppressed")? != &Json::Bool(false)
        || text(feature.field("operation")?)? != "new_body"
    {
        return Err("SC8 requires one unsuppressed initial box feature");
    }
    if feature
        .field("param_names")?
        .items()?
        .iter()
        .map(text)
        .collect::<Result<Vec<_>, _>>()?
        != ["d1", "d2", "d3"]
    {
        return Err("SC8 box dimensions are not resolved canonical parameters");
    }
    let bodies = feature.field("body_names")?.items()?;
    if bodies.len() != 1 || text(&bodies[0])?.is_empty() {
        return Err("SC8 base box must create exactly one body");
    }
    let _ = text(feature.field("name")?)?;
    let lo = vector3(feature.field("corner")?, false)?;
    let size = dimensions(
        integer(feature.field("length")?, true)?,
        integer(feature.field("width")?, true)?,
        integer(feature.field("height")?, true)?,
    )?;
    if lo[2] != 0 || !(-32..=32).contains(&lo[0]) || !(-32..=32).contains(&lo[1]) {
        return Err("SC8 base box must start at z=0 within the bounded XY domain");
    }
    if lo[0].checked_add(size[0]).is_none_or(|x| x > 64)
        || lo[1].checked_add(size[1]).is_none_or(|y| y > 64)
    {
        return Err("SC8 base box exceeds the native domain");
    }
    Ok(BoxGeometry { lo, size })
}

fn timeline_hole(feature: &Json) -> Result<HoleGeometry, &'static str> {
    feature.keys(&[
        "depth",
        "diameter",
        "direction",
        "face_names",
        "hole",
        "id",
        "name",
        "param_names",
        "position",
        "suppressed",
        "type",
    ])?;
    if integer(feature.field("id")?, false)? != 2
        || feature.field("suppressed")? != &Json::Bool(false)
    {
        return Err("SC8 requires one unsuppressed hole feature after the box");
    }
    let params = feature.field("param_names")?.items()?;
    if params.len() != 2 || text(&params[0])? != "d4" || text(&params[1])? != "d5" {
        return Err("SC8 hole depth and diameter must be resolved parameters");
    }
    let hole = feature.field("hole")?;
    hole.keys(&["type"])?;
    if text(hole.field("type")?)? != "simple" {
        return Err("SC8 supports simple flat-bottom holes only");
    }
    let faces = feature.field("face_names")?.items()?;
    if faces.len() != 1 {
        return Err("SC8 hole must select one end face of its box");
    }
    let face = text(&faces[0])?;
    let direction = vector3(feature.field("direction")?, false)?;
    let position = vector3(feature.field("position")?, false)?;
    let (face_z, direction_z) = match face {
        "F1:+z" if direction == [0, 0, -1] => (position[2], -1),
        "F1:-z" if direction == [0, 0, 1] => (position[2], 1),
        _ => return Err("SC8 hole face and direction are unsupported or inconsistent"),
    };
    Ok(HoleGeometry {
        center: [position[0], position[1]],
        diameter: integer(feature.field("diameter")?, true)?,
        depth: integer(feature.field("depth")?, true)?,
        face_z,
        direction_z,
    })
}

fn recipe(box_: BoxGeometry, hole: HoleGeometry) -> Result<Json, &'static str> {
    let [x0, y0, z0] = box_.lo;
    let [length, width, height] = box_.size;
    let x1 = x0
        .checked_add(length)
        .ok_or("SC8 box coordinate overflow")?;
    let y1 = y0.checked_add(width).ok_or("SC8 box coordinate overflow")?;
    if z0 != 0 {
        return Err("SC8 box must start at z=0");
    }
    if hole.diameter <= 0 || hole.diameter % 2 != 0 || hole.depth <= 0 || hole.depth >= height {
        return Err("SC8 requires a positive even diameter and strictly blind depth");
    }
    let radius = hole.diameter / 2;
    let guard = radius.checked_mul(2).ok_or("SC8 guard overflow")?;
    if i64::from(hole.center[0]) - i64::from(guard) <= i64::from(x0)
        || i64::from(hole.center[0]) + i64::from(guard) >= i64::from(x1)
        || i64::from(hole.center[1]) - i64::from(guard) <= i64::from(y0)
        || i64::from(hole.center[1]) + i64::from(guard) >= i64::from(y1)
    {
        return Err("SC8 bore guard must lie strictly inside the box profile");
    }
    let interval = if hole.direction_z == -1 {
        if hole.face_z != height {
            return Err("SC8 top bore position is not on the top face");
        }
        [height - hole.depth, height]
    } else if hole.direction_z == 1 {
        if hole.face_z != z0 {
            return Err("SC8 bottom bore position is not on the bottom face");
        }
        [z0, z0 + hole.depth]
    } else if hole.face_z == height {
        [height - hole.depth, height]
    } else if hole.face_z == z0 {
        [z0, z0 + hole.depth]
    } else {
        return Err("SC8 direct hole must start on a box end face");
    };
    let [cx, cy] = hole.center;
    let points = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
    let vertices = points
        .into_iter()
        .enumerate()
        .map(|(i, xy)| {
            object([
                ("id", Json::String(format!("sc8_v{i}"))),
                (
                    "xy_mm",
                    Json::Array(
                        xy.into_iter()
                            .map(|v| Json::Number(v.to_string()))
                            .collect(),
                    ),
                ),
            ])
        })
        .collect();
    Ok(object([
        ("family", Json::String("polygon_interval_bores_v1".into())),
        (
            "outer",
            object([
                ("id", Json::String("sc8_outer".into())),
                ("vertices", Json::Array(vertices)),
            ]),
        ),
        (
            "bores",
            Json::Array(vec![object([
                ("id", Json::String("blind_bore".into())),
                (
                    "center_mm",
                    Json::Array([cx, cy].map(|v| Json::Number(v.to_string())).into()),
                ),
                ("radius_mm", Json::Number(radius.to_string())),
                ("guard_mm", Json::Number(guard.to_string())),
                (
                    "z_interval_mm",
                    Json::Array(interval.map(|v| Json::Number(v.to_string())).into()),
                ),
            ])]),
        ),
        ("height_mm", Json::Number(height.to_string())),
        (
            "top_shift_mm",
            Json::Array(vec![Json::Number("0".into()), Json::Number("0".into())]),
        ),
    ]))
}

fn build(box_: BoxGeometry, hole: HoleGeometry) -> Result<Json, &'static str> {
    interval_bores_json::build(&recipe(box_, hole)?)
}

/// Import a SolidCraft direct command script containing exactly one box and one simple blind hole.
pub fn import_blind_bore_direct(script: &Json) -> Result<Json, &'static str> {
    script.keys(&["commands"])?;
    let [box_command, hole_command] = script.field("commands")?.items()? else {
        return Err("SC8 direct script requires exactly one box and one hole");
    };
    build(direct_box(box_command)?, direct_hole(hole_command)?)
}

/// Import a saved SolidCraft 0.48.1 box-plus-blind-hole feature history.
pub fn import_blind_bore_timeline(script: &Json) -> Result<Json, &'static str> {
    script.keys(&["about", "commands"])?;
    let about = script.field("about")?;
    about.keys(&["design", "format", "replay", "version"])?;
    if text(about.field("design")?)?.is_empty()
        || text(about.field("format")?)? != "solidcraft/66"
        || text(about.field("version")?)? != "0.48.1"
        || text(about.field("replay")?)?.is_empty()
    {
        return Err("unsupported SC8 SolidCraft history metadata");
    }
    let commands = script.field("commands")?.items()?;
    if commands.len() != 2 {
        return Err("SC8 timeline requires exactly one box and one hole");
    }
    let box_feature = timeline_feature(&commands[0], "box")?;
    let hole_feature = timeline_feature(&commands[1], "hole")?;
    build(timeline_box(box_feature)?, timeline_hole(hole_feature)?)
}

/// Route a direct script or saved feature history through the SC8 subset.
pub fn import_blind_bore(script: &Json) -> Result<Json, &'static str> {
    if script.field("about").is_ok() {
        import_blind_bore_timeline(script)
    } else {
        import_blind_bore_direct(script)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    fn interval_and_center(model: &Json) -> Result<([i32; 2], [i32; 2]), &'static str> {
        let bore = &model.field("bores")?.items()?[0];
        Ok((
            [
                integer(&bore.field("z_interval_mm")?.items()?[0], false)?,
                integer(&bore.field("z_interval_mm")?.items()?[1], false)?,
            ],
            [
                integer(&bore.field("center_mm")?.items()?[0], false)?,
                integer(&bore.field("center_mm")?.items()?[1], false)?,
            ],
        ))
    }

    fn features(model: &Json) -> Result<BTreeSet<String>, &'static str> {
        model
            .field("design_faces")?
            .items()?
            .iter()
            .map(|face| face.field("feature")?.text().map(str::to_owned))
            .collect()
    }

    #[test]
    fn direct_zero_decimal_and_host_top_histories_map_exactly_to_g21() {
        let direct = crate::json::parse(include_str!(
            "../tests/references/solidcraft-blind-bore/centered_top.json"
        ))
        .unwrap();
        let direct_model = import_blind_bore(&direct).unwrap();
        assert_eq!(interval_and_center(&direct_model), Ok(([8, 12], [0, 0])));

        let centered = crate::json::parse(include_str!(
            "../tests/references/solidcraft-blind-bore/centered_top.history.json"
        ))
        .unwrap();
        assert_eq!(import_blind_bore(&centered).unwrap(), direct_model);
        let centered_features = features(&direct_model).unwrap();
        assert!(centered_features.contains("bore/blind_bore/lower"));
        assert!(centered_features.contains("bore/blind_bore/wall"));

        let decentered = crate::json::parse(include_str!(
            "../tests/references/solidcraft-blind-bore/decentered_top.history.json"
        ))
        .unwrap();
        let decentered = import_blind_bore(&decentered).unwrap();
        assert_eq!(interval_and_center(&decentered), Ok(([8, 12], [4, -3])));
        assert!(features(&decentered)
            .unwrap()
            .contains("bore/blind_bore/lower"));
    }

    #[test]
    fn bottom_entry_places_the_cap_at_its_interval_upper_end() {
        let history = crate::json::parse(include_str!(
            "../tests/references/solidcraft-blind-bore/opposite_bottom.history.json"
        ))
        .unwrap();
        let model = import_blind_bore(&history).unwrap();
        assert_eq!(interval_and_center(&model), Ok(([0, 5], [-5, 6])));
        let tags = features(&model).unwrap();
        assert!(tags.contains("bore/blind_bore/upper"));
        assert!(tags.contains("bore/blind_bore/wall"));
        assert!(!tags.contains("bore/blind_bore/lower"));
    }

    #[test]
    fn direct_import_rejects_nonblind_or_unsupported_holes() {
        let mut script = crate::json::parse(include_str!(
            "../tests/references/solidcraft-blind-bore/centered_top.json"
        ))
        .unwrap();
        let Json::Object(root) = &mut script else {
            unreachable!()
        };
        let Json::Array(commands) = root.get_mut("commands").unwrap() else {
            unreachable!()
        };
        let Json::Object(hole) = &mut commands[1] else {
            unreachable!()
        };
        let Json::Object(params) = hole.get_mut("params").unwrap() else {
            unreachable!()
        };
        params.insert("depth".into(), Json::Number("12".into()));
        assert!(import_blind_bore(&script).is_err());
    }
}
