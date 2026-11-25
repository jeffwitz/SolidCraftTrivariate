//! Import the qualified SolidCraft native counterbore subset through G23.
use crate::{
    interchange,
    json::{object, Json},
    solidcraft_boxes::decimal_rational,
};

#[derive(Clone, Copy)]
struct BoxGeometry {
    lo: [i32; 3],
    size: [i32; 3],
}

#[derive(Clone, Copy)]
struct CounterboreGeometry {
    center: [i32; 2],
    pilot_diameter: i32,
    cb_diameter: i32,
    cb_depth: i32,
    opening_z: i32,
    from_top: bool,
}

fn text(value: &Json) -> Result<&str, &'static str> {
    value.text()
}

fn integer(value: &Json, allow_string: bool) -> Result<i32, &'static str> {
    let token = match value {
        Json::Number(value) => value.as_str(),
        Json::String(value) if allow_string => value.as_str(),
        _ => return Err("SC9 values must be literal numbers"),
    };
    let (numerator, denominator) = decimal_rational(token)?;
    if numerator % denominator != 0 {
        return Err("SC9 geometry must use integer millimetres");
    }
    i32::try_from(numerator / denominator).map_err(|_| "SC9 integer range")
}

fn vector3(value: &Json, allow_string: bool) -> Result<[i32; 3], &'static str> {
    let [x, y, z] = value.items()? else {
        return Err("SC9 point must have three coordinates");
    };
    Ok([
        integer(x, allow_string)?,
        integer(y, allow_string)?,
        integer(z, allow_string)?,
    ])
}

fn box_dimensions(length: i32, width: i32, height: i32) -> Result<[i32; 3], &'static str> {
    if !(1..=64).contains(&length) || !(1..=64).contains(&width) || !(1..=64).contains(&height) {
        return Err("SC9 box dimensions outside native range");
    }
    Ok([length, width, height])
}

fn direct_box(command: &Json) -> Result<BoxGeometry, &'static str> {
    command.keys(&["command", "params"])?;
    if text(command.field("command")?)? != "solid.box" {
        return Err("SC9 requires solid.box first");
    }
    let params = command.field("params")?;
    params.keys(&["corner", "height", "length", "name", "operation", "width"])?;
    if text(params.field("operation")?)? != "new" {
        return Err("SC9 requires a new base box");
    }
    let _ = text(params.field("name")?)?;
    let lo = vector3(params.field("corner")?, false)?;
    let size = box_dimensions(
        integer(params.field("length")?, false)?,
        integer(params.field("width")?, false)?,
        integer(params.field("height")?, false)?,
    )?;
    validate_box(BoxGeometry { lo, size })
}

fn validate_box(box_: BoxGeometry) -> Result<BoxGeometry, &'static str> {
    let [x, y, z] = box_.lo;
    let [length, width, _] = box_.size;
    if z != 0
        || !(-32..=32).contains(&x)
        || !(-32..=32).contains(&y)
        || x.checked_add(length).is_none_or(|v| v > 64)
        || y.checked_add(width).is_none_or(|v| v > 64)
    {
        return Err("SC9 base box must be a bounded origin-based XY prism");
    }
    Ok(box_)
}

fn direct_counterbore(command: &Json) -> Result<CounterboreGeometry, &'static str> {
    command.keys(&["command", "params"])?;
    if text(command.field("command")?)? != "solid.hole" {
        return Err("SC9 requires one solid.hole after the box");
    }
    let params = command.field("params")?;
    params.keys(&[
        "cb_depth",
        "cb_diameter",
        "diameter",
        "name",
        "position",
        "type",
    ])?;
    if text(params.field("type")?)? != "counterbore" {
        return Err("SC9 supports native counterbore holes only");
    }
    let _ = text(params.field("name")?)?;
    let [x, y, z] = vector3(params.field("position")?, false)?;
    Ok(CounterboreGeometry {
        center: [x, y],
        pilot_diameter: integer(params.field("diameter")?, false)?,
        cb_diameter: integer(params.field("cb_diameter")?, false)?,
        cb_depth: integer(params.field("cb_depth")?, false)?,
        opening_z: z,
        from_top: z != 0,
    })
}

fn timeline_feature<'a>(command: &'a Json, expected: &str) -> Result<&'a Json, &'static str> {
    command.keys(&["command", "params"])?;
    if text(command.field("command")?)? != "timeline.add" {
        return Err("SC9 accepts timeline.add commands only");
    }
    let params = command.field("params")?;
    params.keys(&["feature"])?;
    let feature = params.field("feature")?;
    if text(feature.field("type")?)? != expected {
        return Err("SC9 feature sequence type mismatch");
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
        || feature
            .field("param_names")?
            .items()?
            .iter()
            .map(text)
            .collect::<Result<Vec<_>, _>>()?
            != ["d1", "d2", "d3"]
    {
        return Err("SC9 requires one canonical unsuppressed initial box");
    }
    let bodies = feature.field("body_names")?.items()?;
    if bodies.len() != 1 || text(&bodies[0])?.is_empty() {
        return Err("SC9 base feature must create exactly one body");
    }
    let _ = text(feature.field("name")?)?;
    validate_box(BoxGeometry {
        lo: vector3(feature.field("corner")?, false)?,
        size: box_dimensions(
            integer(feature.field("length")?, true)?,
            integer(feature.field("width")?, true)?,
            integer(feature.field("height")?, true)?,
        )?,
    })
}

fn timeline_counterbore(feature: &Json) -> Result<CounterboreGeometry, &'static str> {
    feature.keys(&[
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
        || feature
            .field("param_names")?
            .items()?
            .iter()
            .map(text)
            .collect::<Result<Vec<_>, _>>()?
            != ["d4", "d5", "d6"]
    {
        return Err("SC9 requires one unsuppressed canonical counterbore feature");
    }
    let _ = text(feature.field("name")?)?;
    let hole = feature.field("hole")?;
    hole.keys(&["cb_depth", "cb_diameter", "type"])?;
    if text(hole.field("type")?)? != "counterbore" {
        return Err("SC9 supports native counterbore holes only");
    }
    let faces = feature.field("face_names")?.items()?;
    if faces.len() != 1 {
        return Err("SC9 counterbore must select one box end face");
    }
    let [x, y, z] = vector3(feature.field("position")?, false)?;
    let direction = vector3(feature.field("direction")?, false)?;
    let from_top = match (text(&faces[0])?, direction) {
        ("F1:+z", [0, 0, -1]) => true,
        ("F1:-z", [0, 0, 1]) => false,
        _ => return Err("SC9 counterbore face and direction are inconsistent"),
    };
    Ok(CounterboreGeometry {
        center: [x, y],
        pilot_diameter: integer(feature.field("diameter")?, true)?,
        cb_diameter: integer(hole.field("cb_diameter")?, true)?,
        cb_depth: integer(hole.field("cb_depth")?, true)?,
        opening_z: z,
        from_top,
    })
}

fn make_recipe(box_: BoxGeometry, bore: CounterboreGeometry) -> Result<Json, &'static str> {
    let [x0, y0, z0] = box_.lo;
    let [length, width, height] = box_.size;
    let x1 = x0
        .checked_add(length)
        .ok_or("SC9 box coordinate overflow")?;
    let y1 = y0.checked_add(width).ok_or("SC9 box coordinate overflow")?;
    if z0 != 0
        || bore.pilot_diameter <= 0
        || bore.pilot_diameter % 2 != 0
        || bore.cb_diameter <= 0
        || bore.cb_diameter % 2 != 0
    {
        return Err("SC9 diameters must be positive even integer millimetres");
    }
    let r = bore.pilot_diameter / 2;
    let big_r = bore.cb_diameter / 2;
    if big_r <= r || bore.cb_depth <= 0 || bore.cb_depth >= height {
        return Err("SC9 counterbore radius and depth are outside the native domain");
    }
    let guard = big_r.checked_mul(2).ok_or("SC9 guard overflow")?;
    let [cx, cy] = bore.center;
    if i64::from(cx) - i64::from(guard) <= i64::from(x0)
        || i64::from(cx) + i64::from(guard) >= i64::from(x1)
        || i64::from(cy) - i64::from(guard) <= i64::from(y0)
        || i64::from(cy) + i64::from(guard) >= i64::from(y1)
    {
        return Err("SC9 counterbore guard must lie strictly inside the box profile");
    }
    let opening = if bore.from_top {
        if bore.opening_z != height {
            return Err("SC9 top counterbore is not on the top face");
        }
        "end"
    } else {
        if bore.opening_z != 0 {
            return Err("SC9 bottom counterbore is not on the bottom face");
        }
        "start"
    };
    let points = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
    let vertices = points
        .into_iter()
        .enumerate()
        .map(|(i, xy)| {
            object([
                ("id", Json::String(format!("sc9_v{i}"))),
                (
                    "xy_mm",
                    Json::Array(xy.map(|n| Json::Number(n.to_string())).into()),
                ),
            ])
        })
        .collect();
    let recipe = object([
        ("family", Json::String("polygon_counterbores_v1".into())),
        (
            "outer",
            object([
                ("id", Json::String("sc9_outer".into())),
                ("vertices", Json::Array(vertices)),
            ]),
        ),
        (
            "bores",
            Json::Array(vec![object([
                ("id", Json::String("hole_2".into())),
                (
                    "center_mm",
                    Json::Array([cx, cy].map(|n| Json::Number(n.to_string())).into()),
                ),
                ("radius_mm", Json::Number(r.to_string())),
                ("guard_mm", Json::Number(guard.to_string())),
                ("counterbore_radius_mm", Json::Number(big_r.to_string())),
                ("depth_mm", Json::Number(bore.cb_depth.to_string())),
                ("opening", Json::String(opening.into())),
            ])]),
        ),
        ("height_mm", Json::Number(height.to_string())),
        (
            "top_shift_mm",
            Json::Array(vec![Json::Number("0".into()), Json::Number("0".into())]),
        ),
    ]);
    interchange::build(&recipe)
}

/// Import a literal script containing exactly one base box and one counterbore.
pub fn import_counterbore_direct(script: &Json) -> Result<Json, &'static str> {
    script.keys(&["commands"])?;
    let [box_command, hole_command] = script.field("commands")?.items()? else {
        return Err("SC9 direct script requires exactly one box and one counterbore");
    };
    make_recipe(direct_box(box_command)?, direct_counterbore(hole_command)?)
}

/// Import a saved SolidCraft 0.48.1 box plus native counterbore feature history.
pub fn import_counterbore_timeline(script: &Json) -> Result<Json, &'static str> {
    script.keys(&["about", "commands"])?;
    let about = script.field("about")?;
    about.keys(&["design", "format", "replay", "version"])?;
    if text(about.field("design")?)?.is_empty()
        || text(about.field("format")?)? != "solidcraft/66"
        || text(about.field("version")?)? != "0.48.1"
        || text(about.field("replay")?)?.is_empty()
    {
        return Err("unsupported SC9 SolidCraft history metadata");
    }
    let commands = script.field("commands")?.items()?;
    if commands.len() != 2 {
        return Err("SC9 timeline requires one box and one counterbore");
    }
    make_recipe(
        timeline_box(timeline_feature(&commands[0], "box")?)?,
        timeline_counterbore(timeline_feature(&commands[1], "hole")?)?,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn real_host_top_and_bottom_histories_map_to_g23() {
        for (name, expected_center, opening, shoulder, feature) in [
            ("centered_top", [0, 0], "end", 8, "bore/hole_2/shoulder"),
            ("decentered_top", [4, -3], "end", 8, "bore/hole_2/shoulder"),
            (
                "opposite_bottom",
                [-5, 6],
                "start",
                4,
                "bore/hole_2/shoulder",
            ),
        ] {
            let source = match name {
                "centered_top" => crate::json::parse(include_str!(
                    "../tests/references/solidcraft-counterbore/centered_top.history.json"
                )),
                "decentered_top" => crate::json::parse(include_str!(
                    "../tests/references/solidcraft-counterbore/decentered_top.history.json"
                )),
                _ => crate::json::parse(include_str!(
                    "../tests/references/solidcraft-counterbore/opposite_bottom.history.json"
                )),
            }
            .unwrap();
            let model = import_counterbore_timeline(&source).unwrap();
            let bore = &model.field("bores").unwrap().items().unwrap()[0];
            assert_eq!(
                bore.field("center_mm").unwrap(),
                &Json::Array(expected_center.map(|v| Json::Number(v.to_string())).into())
            );
            assert_eq!(
                bore.field("opening").unwrap(),
                &Json::String(opening.into())
            );
            assert_eq!(bore.field("depth_mm").unwrap(), &Json::Number("4".into()));
            assert!(model
                .field("design_faces")
                .unwrap()
                .items()
                .unwrap()
                .iter()
                .any(|face| face.field("feature").unwrap() == &Json::String(feature.into())));
            assert_eq!(
                model.field("schema").unwrap(),
                &Json::String("polygon_counterbores_v1".into())
            );
            assert_eq!(shoulder, if opening == "end" { 8 } else { 4 });
        }
    }

    #[test]
    fn direct_and_timeline_are_identical_and_out_of_contract_inputs_fail() {
        let direct = crate::json::parse(include_str!(
            "../tests/references/solidcraft-counterbore/centered_top.json"
        ))
        .unwrap();
        let timeline = crate::json::parse(include_str!(
            "../tests/references/solidcraft-counterbore/centered_top.history.json"
        ))
        .unwrap();
        assert_eq!(
            import_counterbore_direct(&direct).unwrap(),
            import_counterbore_timeline(&timeline).unwrap()
        );
        assert!(import_counterbore_direct(&direct).is_ok());
        let Json::Object(mut script) = direct else {
            panic!("fixture object")
        };
        let Json::Array(commands) = script.get_mut("commands").unwrap() else {
            panic!("commands array")
        };
        let Json::Object(params) = commands[1].field("params").unwrap() else {
            panic!("params object")
        };
        let mut params = params.clone();
        params.insert("cb_diameter".into(), Json::Number("13".into()));
        commands[1] = object([
            ("command", Json::String("solid.hole".into())),
            ("params", Json::Object(params)),
        ]);
        assert!(import_counterbore_direct(&Json::Object(script)).is_err());
    }
}
