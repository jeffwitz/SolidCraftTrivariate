//! Translate the verified solidcraft/66 bracket timeline subset.
//!
//! This module consumes the canonical history emitted by `solidcraft-cli script`.
//! It never runs a user script in the host kernel. Before constructing anything,
//! it resolves every feature, curve and generated face reference against the
//! sketch's directed six-edge contour and verifies the geometric support.

use crate::{
    interchange,
    json::{object, Json},
    solidcraft_bridge,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Segment {
    Bottom,
    Right,
    Top,
    InnerVertical,
    UpperLeft,
    Left,
}

impl Segment {
    fn directed_endpoints(self) -> ([i64; 2], [i64; 2]) {
        match self {
            Self::Bottom => ([0, 0], [60, 0]),
            Self::Right => ([60, 0], [60, 8]),
            Self::Top => ([60, 8], [8, 8]),
            Self::InnerVertical => ([8, 8], [8, 50]),
            Self::UpperLeft => ([8, 50], [0, 50]),
            Self::Left => ([0, 50], [0, 0]),
        }
    }
}

fn id(value: &Json, expected: i64) -> Result<(), &'static str> {
    if value.number::<i64>()? == expected {
        Ok(())
    } else {
        Err("unexpected timeline feature identity")
    }
}

fn text(value: &Json, expected: &str) -> Result<(), &'static str> {
    if value.text()? == expected {
        Ok(())
    } else {
        Err("unexpected timeline semantic value")
    }
}

fn bool_value(value: &Json, expected: bool) -> Result<(), &'static str> {
    if value == &Json::Bool(expected) {
        Ok(())
    } else {
        Err("suppressed or altered timeline feature")
    }
}

/// Parse a finite decimal token and multiply it by an integer scale exactly.
/// A non-integral scaled result is rejected rather than rounded.
fn scaled_decimal(value: &Json, scale: i128, allow_string: bool) -> Result<i128, &'static str> {
    let token = match value {
        Json::Number(s) => s.as_str(),
        Json::String(s) if allow_string => s.as_str(),
        _ => return Err("timeline value must be a literal number"),
    };
    let (mantissa, exponent) = match token.split_once(['e', 'E']) {
        Some((m, e)) => (m, e.parse::<i32>().map_err(|_| "timeline numeric range")?),
        None => (token, 0),
    };
    let (negative, mantissa) = if let Some(rest) = mantissa.strip_prefix('-') {
        (true, rest)
    } else {
        (false, mantissa)
    };
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    if whole.is_empty()
        || !whole.bytes().all(|b| b.is_ascii_digit())
        || !fraction.bytes().all(|b| b.is_ascii_digit())
    {
        return Err("invalid timeline decimal literal");
    }
    let mut digits = format!("{whole}{fraction}");
    let mut power = exponent
        .checked_sub(i32::try_from(fraction.len()).map_err(|_| "timeline numeric range")?)
        .ok_or("timeline numeric range")?;
    while digits.ends_with('0') {
        digits.pop();
        power = power.checked_add(1).ok_or("timeline numeric range")?;
    }
    if digits.is_empty() {
        return Ok(0);
    }
    let mut result = digits
        .parse::<i128>()
        .map_err(|_| "timeline numeric range")?
        .checked_mul(scale)
        .ok_or("timeline numeric range")?;
    if power >= 0 {
        let factor = 10_i128
            .checked_pow(power as u32)
            .ok_or("timeline numeric range")?;
        result = result.checked_mul(factor).ok_or("timeline numeric range")?;
    } else {
        let divisor = 10_i128
            .checked_pow(power.unsigned_abs())
            .ok_or("timeline numeric range")?;
        if result % divisor != 0 {
            return Err("timeline value is outside the exact parameter lattice");
        }
        result /= divisor;
    }
    if negative {
        result.checked_neg().ok_or("timeline numeric range")
    } else {
        Ok(result)
    }
}

fn scalar(value: &Json, expected: i64) -> Result<(), &'static str> {
    if scaled_decimal(value, 1, false)? == i128::from(expected) {
        Ok(())
    } else {
        Err("unexpected sketch coordinate or direction")
    }
}

fn vector<const N: usize>(value: &Json, expected: [i64; N]) -> Result<(), &'static str> {
    let values = value.items()?;
    if values.len() != N {
        return Err("unexpected timeline vector length");
    }
    for (value, expected) in values.iter().zip(expected) {
        scalar(value, expected)?;
    }
    Ok(())
}

fn feature<'a>(command: &'a Json, expected_type: &str) -> Result<&'a Json, &'static str> {
    command.keys(&["command", "params"])?;
    text(command.field("command")?, "timeline.add")?;
    let params = command.field("params")?;
    params.keys(&["feature"])?;
    let feature = params.field("feature")?;
    text(feature.field("type")?, expected_type)?;
    Ok(feature)
}

fn validate_sketch(sketch_feature: &Json) -> Result<[String; 6], &'static str> {
    sketch_feature.keys(&["id", "name", "plane", "sketch", "suppressed", "type"])?;
    id(sketch_feature.field("id")?, 1)?;
    text(sketch_feature.field("name")?, "Profile")?;
    text(sketch_feature.field("type")?, "sketch")?;
    bool_value(sketch_feature.field("suppressed")?, false)?;

    let plane = sketch_feature.field("plane")?;
    plane.keys(&["name", "type"])?;
    text(plane.field("name")?, "XY")?;
    text(plane.field("type")?, "origin")?;

    let sketch = sketch_feature.field("sketch")?;
    sketch.keys(&["constraints", "counters", "curves", "points"])?;
    let counters = sketch.field("counters")?;
    counters.keys(&["k", "l"])?;
    id(counters.field("k")?, 1)?;
    id(counters.field("l")?, 6)?;

    let points = sketch.field("points")?.items()?;
    if points.len() != 7 {
        return Err("bracket sketch must contain the origin and six vertices");
    }
    let origin = &points[0];
    origin.keys(&["fixed", "id", "pos"])?;
    text(origin.field("id")?, "origin")?;
    bool_value(origin.field("fixed")?, true)?;
    vector(origin.field("pos")?, [0, 0])?;

    let mut point_positions = [[0_i64; 2]; 7];
    let mut point_ids = Vec::with_capacity(7);
    for (index, point) in points.iter().enumerate() {
        point.keys(&["fixed", "id", "pos"])?;
        let point_id = point.field("id")?.text()?.to_owned();
        if point_id.is_empty() || point_ids.contains(&point_id) {
            return Err("sketch point ids must be nonempty and unique");
        }
        point_ids.push(point_id);
        let is_origin = index == 0;
        bool_value(point.field("fixed")?, is_origin)?;
        let pos = point.field("pos")?.items()?;
        if pos.len() != 2 {
            return Err("sketch point must have two coordinates");
        }
        for axis in 0..2 {
            point_positions[index][axis] = i64::try_from(scaled_decimal(&pos[axis], 1, false)?)
                .map_err(|_| "sketch coordinate range")?;
        }
    }

    let constraints = sketch.field("constraints")?.items()?;
    if constraints.len() != 1 {
        return Err("unsupported or missing sketch constraint");
    }
    let constraint = &constraints[0];
    constraint.keys(&["id", "p", "q", "type"])?;
    text(constraint.field("id")?, "k1")?;
    text(constraint.field("type")?, "coincident")?;
    id(constraint.field("p")?, 0)?;
    let origin_vertex = (1..=6)
        .find(|&index| point_positions[index] == [0, 0])
        .ok_or("sketch has no origin vertex")?;
    if (1..=6)
        .filter(|&index| point_positions[index] == [0, 0])
        .count()
        != 1
    {
        return Err("sketch origin vertex is ambiguous");
    }
    id(
        constraint.field("q")?,
        i64::try_from(origin_vertex).map_err(|_| "sketch point index range")?,
    )?;

    let curves = sketch.field("curves")?.items()?;
    if curves.len() != 6 {
        return Err("bracket sketch must contain six line segments");
    }
    let expected_segments = [
        Segment::Bottom,
        Segment::Right,
        Segment::Top,
        Segment::InnerVertical,
        Segment::UpperLeft,
        Segment::Left,
    ];
    let mut semantic_ids: [Option<String>; 6] = std::array::from_fn(|_| None);
    let mut seen_ids = Vec::<String>::new();
    let mut outgoing = [0_u8; 7];
    let mut incoming = [0_u8; 7];
    for curve in curves {
        curve.keys(&["a", "b", "construction", "id", "type"])?;
        text(curve.field("type")?, "line")?;
        bool_value(curve.field("construction")?, false)?;
        let curve_id = curve.field("id")?.text()?.to_owned();
        if curve_id.is_empty() || seen_ids.contains(&curve_id) {
            return Err("sketch curve ids must be nonempty and unique");
        }
        seen_ids.push(curve_id.clone());
        let a = usize::try_from(curve.field("a")?.number::<i64>()?)
            .map_err(|_| "invalid sketch endpoint index")?;
        let b = usize::try_from(curve.field("b")?.number::<i64>()?)
            .map_err(|_| "invalid sketch endpoint index")?;
        if !(1..=6).contains(&a) || !(1..=6).contains(&b) || a == b {
            return Err("sketch line endpoint is outside the vertex graph");
        }
        outgoing[a] += 1;
        incoming[b] += 1;
        let pa = point_positions[a];
        let pb = point_positions[b];
        let Some((semantic, _)) = expected_segments
            .iter()
            .enumerate()
            .find(|(_, segment)| segment.directed_endpoints() == (pa, pb))
        else {
            return Err("sketch lines do not form the supported directed L contour");
        };
        if semantic_ids[semantic].replace(curve_id).is_some() {
            return Err("duplicate geometric segment in sketch contour");
        }
    }
    for index in 1..=6 {
        if outgoing[index] != 1 || incoming[index] != 1 {
            return Err("sketch line graph is not one directed six-edge loop");
        }
    }
    if semantic_ids.iter().any(Option::is_none) {
        return Err("sketch line graph omits a contour segment");
    }
    Ok(semantic_ids.map(Option::unwrap))
}

fn validate_extrude(feature: &Json) -> Result<(), &'static str> {
    feature.keys(&[
        "body_names",
        "extent",
        "id",
        "name",
        "operation",
        "param_names",
        "profiles",
        "sketch",
        "suppressed",
        "type",
    ])?;
    id(feature.field("id")?, 2)?;
    text(feature.field("name")?, "Extrude1")?;
    text(feature.field("type")?, "extrude")?;
    bool_value(feature.field("suppressed")?, false)?;
    text(feature.field("operation")?, "new_body")?;
    id(feature.field("sketch")?, 1)?;
    let body_names = feature.field("body_names")?.items()?;
    if body_names.len() != 1 {
        return Err("extrusion must create one bracket body");
    }
    text(&body_names[0], "Bracket")?;
    let params = feature.field("param_names")?.items()?;
    if params.len() != 1 {
        return Err("extrusion parameter reference is unresolved");
    }
    text(&params[0], "d1")?;
    let profiles = feature.field("profiles")?;
    profiles.keys(&["type"])?;
    text(profiles.field("type")?, "all")?;
    let extent = feature.field("extent")?;
    extent.keys(&["convention", "direction", "distance"])?;
    text(extent.field("convention")?, "total_distance")?;
    text(extent.field("direction")?, "positive")?;
    if scaled_decimal(extent.field("distance")?, 1, true)? != 40 {
        return Err("unsupported bracket extrusion distance");
    }
    Ok(())
}

fn side_name(value: &Json) -> Result<&str, &'static str> {
    let name = value.text()?;
    name.strip_prefix("F2:side:")
        .ok_or("generated face or edge does not reference the extrusion")
}

fn curve_index(name: &str, curve_ids: &[String; 6]) -> Result<usize, &'static str> {
    curve_ids
        .iter()
        .position(|id| id == name)
        .ok_or("timeline references an unknown sketch curve")
}

fn validate_fillet(feature: &Json, curve_ids: &[String; 6]) -> Result<(), &'static str> {
    feature.keys(&[
        "edge_names",
        "edge_refs",
        "edges",
        "id",
        "name",
        "param_names",
        "radius",
        "suppressed",
        "type",
    ])?;
    id(feature.field("id")?, 3)?;
    text(feature.field("name")?, "Fillet1")?;
    text(feature.field("type")?, "fillet")?;
    bool_value(feature.field("suppressed")?, false)?;
    let params = feature.field("param_names")?.items()?;
    if params.len() != 1 {
        return Err("fillet parameter reference is unresolved");
    }
    text(&params[0], "d2")?;
    if scaled_decimal(feature.field("radius")?, 8, true)? != 48 {
        return Err("unsupported inside-corner fillet radius");
    }
    let edges = feature.field("edges")?.items()?;
    if edges.len() != 1 {
        return Err("fillet must select the single inside vertical edge");
    }
    vector(&edges[0], [8, 8, 20])?;
    let edge_names = feature.field("edge_names")?.items()?;
    if edge_names.len() != 1 {
        return Err("fillet edge reference count mismatch");
    }
    let names = edge_names[0].text()?.split('|').collect::<Vec<_>>();
    if names.len() != 2 {
        return Err("fillet must reference the two contour sides at the inside corner");
    }
    let mut semantic = Vec::new();
    for name in names {
        let index = curve_index(side_name(&Json::String(name.to_owned()))?, curve_ids)?;
        semantic.push(index);
    }
    semantic.sort_unstable();
    if semantic != [2, 3] {
        return Err("fillet reference is not the inside contour corner");
    }
    let refs = feature.field("edge_refs")?.items()?;
    if refs.len() != 1 {
        return Err("fillet analytic edge reference count mismatch");
    }
    let edge = &refs[0];
    edge.keys(&["dir", "max", "min"])?;
    vector(edge.field("dir")?, [0, 0, 1])?;
    vector(edge.field("min")?, [0, 0, 0])?;
    vector(edge.field("max")?, [60, 50, 40])?;
    Ok(())
}

#[derive(Clone, Copy)]
struct Hole {
    support: Segment,
    diameter_eighths: i128,
}

fn validate_hole(
    feature: &Json,
    expected_id: i64,
    expected_param: &str,
    curve_ids: &[String; 6],
) -> Result<Hole, &'static str> {
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
    id(feature.field("id")?, expected_id)?;
    let expected_name = if expected_id == 4 { "Hole1" } else { "Hole2" };
    text(feature.field("name")?, expected_name)?;
    text(feature.field("type")?, "hole")?;
    bool_value(feature.field("suppressed")?, false)?;
    let params = feature.field("param_names")?.items()?;
    if params.len() != 1 {
        return Err("hole parameter reference is unresolved");
    }
    text(&params[0], expected_param)?;
    let hole_type = feature.field("hole")?;
    hole_type.keys(&["type"])?;
    text(hole_type.field("type")?, "simple")?;

    let faces = feature.field("face_names")?.items()?;
    if faces.len() != 1 {
        return Err("hole must select one generated side face");
    }
    let face = side_name(&faces[0])?;
    let line_index = curve_index(face, curve_ids)?;
    let diameter_eighths = scaled_decimal(feature.field("diameter")?, 8, true)?;
    if !(40..=80).contains(&diameter_eighths) {
        return Err("hole diameter is outside the SC1 qualified lattice");
    }

    let support = if line_index == 2 {
        vector(feature.field("position")?, [40, 8, 20])?;
        vector(feature.field("direction")?, [0, -1, 0])?;
        Segment::Top
    } else if line_index == 3 {
        vector(feature.field("position")?, [8, 32, 20])?;
        vector(feature.field("direction")?, [-1, 0, 0])?;
        Segment::InnerVertical
    } else {
        return Err("hole support is not one of the two qualified bracket faces");
    };
    Ok(Hole {
        support,
        diameter_eighths,
    })
}

fn decimal_eighths(ticks: i128) -> String {
    let whole = ticks.div_euclid(8);
    match ticks.rem_euclid(8) {
        0 => whole.to_string(),
        1 => format!("{whole}.125"),
        2 => format!("{whole}.25"),
        3 => format!("{whole}.375"),
        4 => format!("{whole}.5"),
        5 => format!("{whole}.625"),
        6 => format!("{whole}.75"),
        _ => format!("{whole}.875"),
    }
}

fn point(values: &[i64]) -> Json {
    Json::Array(
        values
            .iter()
            .map(|value| Json::Number(value.to_string()))
            .collect(),
    )
}

fn direct_script(curve_ids: &[String; 6], holes: [Hole; 2]) -> Result<Json, &'static str> {
    let _ = curve_ids;
    let mut ordered = holes;
    ordered.sort_by_key(|hole| match hole.support {
        Segment::Top => 0,
        Segment::InnerVertical => 1,
        _ => 2,
    });
    if ordered[0].support != Segment::Top || ordered[1].support != Segment::InnerVertical {
        return Err("timeline must contain one hole on each qualified face");
    }
    let command = |name: &str, params: Json| {
        object([
            ("command", Json::String(name.to_owned())),
            ("params", params),
        ])
    };
    Ok(object([(
        "commands",
        Json::Array(vec![
            command(
                "sketch.create",
                object([
                    ("plane", Json::String("XY".into())),
                    ("name", Json::String("Profile".into())),
                ]),
            ),
            command(
                "sketch.line",
                object([
                    (
                        "points",
                        Json::Array(
                            [[0, 0], [60, 0], [60, 8], [8, 8], [8, 50], [0, 50]]
                                .iter()
                                .map(|values| point(values))
                                .collect(),
                        ),
                    ),
                    ("closed", Json::Bool(true)),
                ]),
            ),
            command("sketch.finish", object([])),
            command(
                "solid.extrude",
                object([
                    ("distance", Json::Number("40".into())),
                    ("body_name", Json::String("Bracket".into())),
                ]),
            ),
            command(
                "solid.fillet",
                object([
                    ("edges", Json::Array(vec![point(&[8, 8, 20])])),
                    ("radius", Json::Number("6".into())),
                ]),
            ),
            command(
                "solid.hole",
                object([
                    ("position", point(&[40, 8, 20])),
                    (
                        "diameter",
                        Json::Number(decimal_eighths(ordered[0].diameter_eighths)),
                    ),
                ]),
            ),
            command(
                "solid.hole",
                object([
                    ("position", point(&[8, 32, 20])),
                    (
                        "diameter",
                        Json::Number(decimal_eighths(ordered[1].diameter_eighths)),
                    ),
                ]),
            ),
        ]),
    )]))
}

/// Translate a validated canonical solidcraft/66 bracket timeline into the
/// native audited parameterized bracket JSON model.
pub fn translate_timeline(script: &Json) -> Result<Json, &'static str> {
    script.keys(&["about", "commands"])?;
    let about = script.field("about")?;
    about.keys(&["design", "format", "replay", "version"])?;
    let design = about.field("design")?.text()?;
    if design != "editable" && design != "original" {
        return Err("timeline is not an editable or original SolidCraft design");
    }
    text(about.field("format")?, "solidcraft/66")?;
    if about.field("replay")?.text()?.is_empty() || about.field("version")?.text()?.is_empty() {
        return Err("incomplete SolidCraft timeline metadata");
    }
    let commands = script.field("commands")?.items()?;
    if commands.len() != 5 {
        return Err("SolidCraft bracket timeline must contain five active features");
    }
    let sketch = feature(&commands[0], "sketch")?;
    let curve_ids = validate_sketch(sketch)?;
    let extrude = feature(&commands[1], "extrude")?;
    validate_extrude(extrude)?;
    let fillet = feature(&commands[2], "fillet")?;
    validate_fillet(fillet, &curve_ids)?;
    let first = feature(&commands[3], "hole")?;
    let first = validate_hole(first, 4, "d3", &curve_ids)?;
    let second = feature(&commands[4], "hole")?;
    let second = validate_hole(second, 5, "d4", &curve_ids)?;
    if first.support == second.support {
        return Err("timeline duplicates a bracket hole support");
    }
    let direct = direct_script(&curve_ids, [first, second])?;
    let recipe = solidcraft_bridge::bracket_recipe(&direct)?;
    let model = interchange::build(&recipe)?;
    let audit = interchange::audit(&model)?;
    if audit.field("status")?.text()? != "valid" {
        return Err("translated timeline failed the native model audit");
    }
    Ok(model)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    const HOST_SCRIPT: &str = r#"{
      "about":{"design":"editable","format":"solidcraft/66","replay":"canonical","version":"0.48.1"},
      "commands":[
        {"command":"timeline.add","params":{"feature":{"id":1,"name":"Profile","plane":{"name":"XY","type":"origin"},"sketch":{"constraints":[{"id":"k1","p":0,"q":1,"type":"coincident"}],"counters":{"k":1,"l":6},"curves":[{"a":1,"b":2,"construction":false,"id":"l1","type":"line"},{"a":2,"b":3,"construction":false,"id":"l2","type":"line"},{"a":3,"b":4,"construction":false,"id":"l3","type":"line"},{"a":4,"b":5,"construction":false,"id":"l4","type":"line"},{"a":5,"b":6,"construction":false,"id":"l5","type":"line"},{"a":6,"b":1,"construction":false,"id":"l6","type":"line"}],"points":[{"fixed":true,"id":"origin","pos":[0,0]},{"fixed":false,"id":"l1.start","pos":[0,0]},{"fixed":false,"id":"l1.end","pos":[60,0]},{"fixed":false,"id":"l2.end","pos":[60,8]},{"fixed":false,"id":"l3.end","pos":[8,8]},{"fixed":false,"id":"l4.end","pos":[8,50]},{"fixed":false,"id":"l5.end","pos":[0,50]}]},"suppressed":false,"type":"sketch"}}},
        {"command":"timeline.add","params":{"feature":{"body_names":["Bracket"],"extent":{"convention":"total_distance","direction":"positive","distance":"40"},"id":2,"name":"Extrude1","operation":"new_body","param_names":["d1"],"profiles":{"type":"all"},"sketch":1,"suppressed":false,"type":"extrude"}}},
        {"command":"timeline.add","params":{"feature":{"edge_names":["F2:side:l3|F2:side:l4"],"edge_refs":[{"dir":[0,0,1],"max":[60,50,40],"min":[0,0,0]}],"edges":[[8,8,20]],"id":3,"name":"Fillet1","param_names":["d2"],"radius":"6","suppressed":false,"type":"fillet"}}},
        {"command":"timeline.add","params":{"feature":{"diameter":"8","direction":[0,-1,0],"face_names":["F2:side:l3"],"hole":{"type":"simple"},"id":4,"name":"Hole1","param_names":["d3"],"position":[40,8,20],"suppressed":false,"type":"hole"}}},
        {"command":"timeline.add","params":{"feature":{"diameter":"8","direction":[-1,0,0],"face_names":["F2:side:l4"],"hole":{"type":"simple"},"id":5,"name":"Hole2","param_names":["d4"],"position":[8,32,20],"suppressed":false,"type":"hole"}}}
      ]
    }"#;

    fn host() -> Json {
        crate::json::parse(HOST_SCRIPT).unwrap()
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
    fn translates_semantic_timeline_to_native_model() {
        let script = host();
        let model = translate_timeline(&script).unwrap();
        let direct = direct_script(
            &["l1", "l2", "l3", "l4", "l5", "l6"].map(str::to_owned),
            [
                Hole {
                    support: Segment::Top,
                    diameter_eighths: 64,
                },
                Hole {
                    support: Segment::InnerVertical,
                    diameter_eighths: 64,
                },
            ],
        )
        .unwrap();
        let direct_model =
            interchange::build(&solidcraft_bridge::bracket_recipe(&direct).unwrap()).unwrap();
        assert_eq!(model, direct_model);
        assert_eq!(
            interchange::audit(&model)
                .unwrap()
                .field("status")
                .unwrap()
                .text()
                .unwrap(),
            "valid"
        );
    }

    #[test]
    fn rejects_suppression_unresolved_constraints_and_wrong_references() {
        let mut script = host();
        feature_mut(&mut script, 0).insert("suppressed".into(), Json::Bool(true));
        assert!(translate_timeline(&script).is_err());

        let mut script = host();
        let feature = feature_mut(&mut script, 0);
        let Json::Object(sketch) = feature.get_mut("sketch").unwrap() else {
            unreachable!()
        };
        let Json::Array(constraints) = sketch.get_mut("constraints").unwrap() else {
            unreachable!()
        };
        constraints.push(object([
            ("id", Json::String("k2".into())),
            ("p", Json::Number("1".into())),
            ("q", Json::Number("2".into())),
            ("type", Json::String("distance".into())),
        ]));
        assert!(translate_timeline(&script).is_err());

        let mut script = host();
        let feature = feature_mut(&mut script, 2);
        let Json::Array(names) = feature.get_mut("edge_names").unwrap() else {
            unreachable!()
        };
        names[0] = Json::String("F2:side:l2|F2:side:l4".into());
        assert!(translate_timeline(&script).is_err());

        let mut script = host();
        let feature = feature_mut(&mut script, 3);
        feature.insert(
            "direction".into(),
            Json::Array(vec![
                Json::Number("1".into()),
                Json::Number("0".into()),
                Json::Number("0".into()),
            ]),
        );
        assert!(translate_timeline(&script).is_err());
    }

    #[test]
    fn rejects_fractional_parameters_outside_direct_lattice() {
        let mut script = host();
        let feature = feature_mut(&mut script, 3);
        feature.insert(
            "diameter".into(),
            Json::String("8.00000000000000001".into()),
        );
        assert!(translate_timeline(&script).is_err());
        let mut script = host();
        let feature = feature_mut(&mut script, 1);
        let Json::Object(extent) = feature.get_mut("extent").unwrap() else {
            unreachable!()
        };
        extent.insert("distance".into(), Json::String("d1".into()));
        assert!(translate_timeline(&script).is_err());
    }
}
