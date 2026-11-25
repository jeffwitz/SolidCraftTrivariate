//! Translate a narrow, explicit SolidCraft sketch/extrude/fillet/hole history.
//!
//! The adapter reads the archived feature graph and its geometric references;
//! it does not execute host commands and does not recognize a part by name or
//! dimensions. Every selected edge and face is checked against the source sketch.

use crate::{
    interchange,
    json::{object, Json},
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone)]
struct Line {
    id: String,
    a: usize,
    b: usize,
    pa: [i32; 2],
    pb: [i32; 2],
}

fn literal_ratio(value: &Json) -> Result<(i128, i128), &'static str> {
    let token = match value {
        Json::Number(s) | Json::String(s) => s.as_str(),
        _ => return Err("solidcraft sketch requires numeric literals"),
    };
    let (mantissa, exponent) = match token.split_once(['e', 'E']) {
        Some((m, e)) => (
            m,
            e.parse::<i32>()
                .map_err(|_| "solidcraft numeric exponent")?,
        ),
        None => (token, 0),
    };
    let (negative, mantissa) = match mantissa.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, mantissa),
    };
    if mantissa.starts_with('+') || mantissa.is_empty() {
        return Err("invalid solidcraft numeric literal");
    }
    let (whole, fractional) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    if whole.is_empty()
        || !whole.bytes().all(|b| b.is_ascii_digit())
        || !fractional.bytes().all(|b| b.is_ascii_digit())
    {
        return Err("invalid solidcraft numeric literal");
    }
    let mut digits = format!("{whole}{fractional}");
    let first_nonzero = digits.bytes().position(|b| b != b'0');
    let Some(first_nonzero) = first_nonzero else {
        return Ok((0, 1));
    };
    digits.drain(..first_nonzero);
    let mut power = exponent
        .checked_sub(i32::try_from(fractional.len()).map_err(|_| "solidcraft numeric range")?)
        .ok_or("solidcraft numeric range")?;
    while digits.ends_with('0') {
        digits.pop();
        power = power.checked_add(1).ok_or("solidcraft numeric range")?;
    }
    let magnitude = digits
        .parse::<i128>()
        .map_err(|_| "solidcraft numeric range")?;
    let (numerator, denominator) = if power >= 0 {
        let scale = 10_i128
            .checked_pow(power as u32)
            .ok_or("solidcraft numeric range")?;
        (
            magnitude
                .checked_mul(scale)
                .ok_or("solidcraft numeric range")?,
            1,
        )
    } else {
        let denominator = 10_i128
            .checked_pow(power.unsigned_abs())
            .ok_or("solidcraft numeric range")?;
        (magnitude, denominator)
    };
    Ok((if negative { -numerator } else { numerator }, denominator))
}

fn exact_i32(value: &Json) -> Result<i32, &'static str> {
    let (n, d) = literal_ratio(value)?;
    if n % d != 0 {
        return Err("solidcraft coordinate is not an integer millimeter value");
    }
    i32::try_from(n / d).map_err(|_| "solidcraft coordinate range")
}

fn exact_vector<const N: usize>(value: &Json) -> Result<[i32; N], &'static str> {
    value
        .items()?
        .iter()
        .map(exact_i32)
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| "solidcraft vector length")
}

fn literal_length(value: &Json) -> Result<i32, &'static str> {
    exact_i32(value)
        .map_err(|_| "solidcraft feature dimension must be an integer millimeter literal")
}

fn strict_false(value: &Json) -> Result<(), &'static str> {
    if value == &Json::Bool(false) {
        Ok(())
    } else {
        Err("suppressed SolidCraft feature")
    }
}

fn feature<'a>(step: &'a Json, required_type: &str) -> Result<&'a Json, &'static str> {
    step.keys(&["command", "params"])?;
    if step.field("command")?.text()? != "timeline.add" {
        return Err("unsupported SolidCraft timeline command");
    }
    let params = step.field("params")?;
    params.keys(&["feature"])?;
    let f = params.field("feature")?;
    if f.field("type")?.text()? != required_type {
        return Err("unsupported SolidCraft feature ordering");
    }
    Ok(f)
}

fn feature_id(f: &Json) -> Result<i64, &'static str> {
    f.field("id")?
        .number::<i64>()
        .map_err(|_| "invalid SolidCraft feature id")
}

fn feature_name(prefix: &str, id: i64) -> String {
    format!("{prefix}_{id}")
}

type SketchOutline = (String, Vec<[i32; 2]>, Vec<Line>);

fn sketch_loop(sketch_feature: &Json) -> Result<SketchOutline, &'static str> {
    sketch_feature.keys(&["id", "name", "plane", "sketch", "suppressed", "type"])?;
    strict_false(sketch_feature.field("suppressed")?)?;
    let plane = sketch_feature.field("plane")?;
    plane.keys(&["name", "type"])?;
    if plane.field("type")?.text()? != "origin" || plane.field("name")?.text()? != "XY" {
        return Err("SolidCraft sketch must lie on the origin XY plane");
    }
    let sk = sketch_feature.field("sketch")?;
    sk.keys(&["constraints", "counters", "curves", "points"])?;
    let points = sk
        .field("points")?
        .items()?
        .iter()
        .map(|p| {
            p.keys(&["fixed", "id", "pos"])?;
            let id = p.field("id")?.text()?;
            if id == "origin" && p.field("fixed")? != &Json::Bool(true) {
                return Err("SolidCraft sketch origin must be fixed");
            }
            Ok((id.to_owned(), exact_vector::<2>(p.field("pos")?)?))
        })
        .collect::<Result<Vec<_>, &'static str>>()?;
    if points
        .iter()
        .any(|(_, p)| p.iter().any(|v| !(-4096..=4096).contains(v)))
    {
        return Err("SolidCraft outline exceeds sketch coordinate range");
    }
    if points.is_empty() {
        return Err("empty SolidCraft sketch points");
    }
    let point_ids = points.iter().map(|p| p.0.clone()).collect::<BTreeSet<_>>();
    if point_ids.len() != points.len() {
        return Err("duplicate SolidCraft sketch point id");
    }
    let count = sk.field("counters")?;
    count.keys(&["k", "l"])?;
    let curves = sk.field("curves")?.items()?;
    if count.field("l")?.number::<usize>()? != curves.len() {
        return Err("SolidCraft sketch curve counter mismatch");
    }
    let mut lines = Vec::new();
    let mut line_names = BTreeSet::new();
    for c in curves {
        c.keys(&["a", "b", "construction", "id", "type"])?;
        if c.field("type")?.text()? != "line" || c.field("construction")? != &Json::Bool(false) {
            return Err("only active straight SolidCraft sketch lines are supported");
        }
        let id = c.field("id")?.text()?.to_owned();
        if !line_names.insert(id.clone()) {
            return Err("duplicate SolidCraft sketch curve id");
        }
        let a = c.field("a")?.number::<usize>()?;
        let b = c.field("b")?.number::<usize>()?;
        if a == b || a >= points.len() || b >= points.len() {
            return Err("SolidCraft sketch curve point reference out of range");
        }
        lines.push(Line {
            id,
            a,
            b,
            pa: points[a].1,
            pb: points[b].1,
        });
    }
    let mut constraint_ids = BTreeSet::new();
    for c in sk.field("constraints")?.items()? {
        c.keys(&["id", "p", "q", "type"])?;
        if !constraint_ids.insert(c.field("id")?.text()?.to_owned())
            || c.field("type")?.text()? != "coincident"
        {
            return Err("unsupported or duplicate SolidCraft sketch constraint");
        }
        let p = c.field("p")?.number::<usize>()?;
        let q = c.field("q")?.number::<usize>()?;
        if p >= points.len() || q >= points.len() || points[p].1 != points[q].1 {
            return Err("SolidCraft coincident constraint is not satisfied by stored points");
        }
    }
    if lines.len() < 3 {
        return Err("SolidCraft profile needs three or more line segments");
    }
    let mut incidence: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (i, line) in lines.iter().enumerate() {
        incidence.entry(line.a).or_default().push(i);
        incidence.entry(line.b).or_default().push(i);
    }
    if incidence.values().any(|edges| edges.len() != 2) {
        return Err("SolidCraft profile is not a single closed contour");
    }
    let start_edge = 0usize;
    let start_point = lines[start_edge].a;
    let mut current_point = start_point;
    let mut current_edge = start_edge;
    let mut used_edges = BTreeSet::new();
    let mut vertices = Vec::with_capacity(lines.len());
    loop {
        if !used_edges.insert(current_edge) {
            return Err("SolidCraft profile repeats a curve before closure");
        }
        vertices.push(points[current_point].1);
        let line = &lines[current_edge];
        let next_point = if line.a == current_point {
            line.b
        } else {
            line.a
        };
        if next_point == start_point {
            break;
        }
        let next_edge = incidence[&next_point]
            .iter()
            .copied()
            .find(|edge| *edge != current_edge)
            .ok_or("open SolidCraft profile")?;
        current_point = next_point;
        current_edge = next_edge;
        if used_edges.len() > lines.len() {
            return Err("SolidCraft profile cycle budget");
        }
    }
    if used_edges.len() != lines.len() {
        return Err("SolidCraft sketch contains multiple contours");
    }
    let constraint_count = sk.field("constraints")?.items()?.len();
    if count.field("k")?.number::<usize>()? != constraint_count {
        return Err("SolidCraft sketch constraint counter mismatch");
    }
    Ok((
        feature_name("sketch", feature_id(sketch_feature)?),
        vertices,
        lines,
    ))
}

fn reference_line<'a>(
    name: &str,
    prefix: &str,
    lines: &'a [Line],
) -> Result<&'a Line, &'static str> {
    let id = name
        .strip_prefix(prefix)
        .ok_or("unsupported SolidCraft face reference")?;
    lines
        .iter()
        .find(|line| line.id == id)
        .ok_or("SolidCraft face references an unknown sketch edge")
}

fn signed_area2(points: &[[i32; 2]]) -> i64 {
    points
        .iter()
        .copied()
        .zip(points.iter().copied().cycle().skip(1))
        .map(|(a, b)| i64::from(a[0]) * i64::from(b[1]) - i64::from(a[1]) * i64::from(b[0]))
        .sum()
}

fn validate_extrusion(f: &Json, sketch_id: i64) -> Result<i32, &'static str> {
    f.keys(&[
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
    strict_false(f.field("suppressed")?)?;
    if f.field("sketch")?.number::<i64>()? != sketch_id
        || f.field("operation")?.text()? != "new_body"
        || f.field("profiles")?.field("type")?.text()? != "all"
        || f.field("body_names")?.items()?.len() != 1
    {
        return Err("SolidCraft extrusion references unsupported sketch or body operation");
    }
    let extent = f.field("extent")?;
    extent.keys(&["convention", "direction", "distance"])?;
    if extent.field("convention")?.text()? != "total_distance"
        || extent.field("direction")?.text()? != "positive"
    {
        return Err("SolidCraft extrusion must use positive total distance");
    }
    let height = literal_length(extent.field("distance")?)?;
    if height <= 0 {
        return Err("SolidCraft extrusion height must be positive");
    }
    Ok(height)
}

fn bounds(vertices: &[[i32; 2]], height: i32) -> [[i32; 3]; 2] {
    let mut low = [i32::MAX, i32::MAX, 0];
    let mut high = [i32::MIN, i32::MIN, height];
    for p in vertices {
        for axis in 0..2 {
            low[axis] = low[axis].min(p[axis]);
            high[axis] = high[axis].max(p[axis]);
        }
    }
    [low, high]
}

fn parse_fillet(
    f: &Json,
    vertices: &[[i32; 2]],
    lines: &[Line],
    prefix: &str,
    height: i32,
) -> Result<(String, [i32; 2], i32), &'static str> {
    f.keys(&[
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
    strict_false(f.field("suppressed")?)?;
    let id = feature_id(f)?;
    let names = f.field("edge_names")?.items()?;
    let edges = f.field("edges")?.items()?;
    let refs = f.field("edge_refs")?.items()?;
    if names.len() != 1 || edges.len() != 1 || refs.len() != 1 {
        return Err("exactly one vertical SolidCraft fillet edge is supported");
    }
    let selected = names[0].text()?;
    let (first, second) = selected
        .split_once('|')
        .ok_or("SolidCraft fillet edge pair required")?;
    let first_line = reference_line(first, prefix, lines)?;
    let second_line = reference_line(second, prefix, lines)?;
    let p = exact_vector::<3>(&edges[0])?;
    if i64::from(p[2]) * 2 != i64::from(height) {
        return Err("SolidCraft fillet edge point must lie at extrusion mid-height");
    }
    let vertex = [p[0], p[1]];
    if !vertices.contains(&vertex) {
        return Err("SolidCraft fillet edge is not at a sketch vertex");
    }
    let incident = lines
        .iter()
        .filter(|line| line.pa == vertex || line.pb == vertex)
        .map(|line| line.id.as_str())
        .collect::<BTreeSet<_>>();
    if incident.len() != 2
        || !incident.contains(first_line.id.as_str())
        || !incident.contains(second_line.id.as_str())
    {
        return Err("SolidCraft fillet edge references do not meet at selected vertex");
    }
    let edge_ref = &refs[0];
    edge_ref.keys(&["dir", "max", "min"])?;
    if exact_vector::<3>(edge_ref.field("dir")?)? != [0, 0, 1]
        || exact_vector::<3>(edge_ref.field("min")?)? != bounds(vertices, height)[0]
        || exact_vector::<3>(edge_ref.field("max")?)? != bounds(vertices, height)[1]
    {
        return Err("SolidCraft fillet edge reference bounds or direction mismatch");
    }
    let radius = literal_length(f.field("radius")?)?;
    if radius <= 0 {
        return Err("SolidCraft fillet radius must be positive");
    }
    Ok((feature_name("fillet", id), vertex, radius))
}

fn parse_hole(
    f: &Json,
    vertices: &[[i32; 2]],
    lines: &[Line],
    prefix: &str,
    height: i32,
) -> Result<(String, String, [i32; 3], i32), &'static str> {
    f.keys(&[
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
    strict_false(f.field("suppressed")?)?;
    if f.field("hole")?.field("type")?.text()? != "simple" {
        return Err("only simple through holes are supported");
    }
    let names = f.field("face_names")?.items()?;
    if names.len() != 1 {
        return Err("SolidCraft hole must reference exactly one side face");
    }
    let line = reference_line(names[0].text()?, prefix, lines)?;
    let position = exact_vector::<3>(f.field("position")?)?;
    let direction = exact_vector::<3>(f.field("direction")?)?;
    let (a, b) = (line.pa, line.pb);
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let orientation = signed_area2(vertices).signum();
    let inward = [
        -dy * i32::try_from(orientation).unwrap(),
        dx * i32::try_from(orientation).unwrap(),
    ];
    let normal_axis = if inward[0] != 0 {
        0
    } else if inward[1] != 0 {
        1
    } else {
        return Err("degenerate hole support edge");
    };
    let mut expected_direction = [0, 0, 0];
    expected_direction[normal_axis] = inward[normal_axis].signum();
    if direction != expected_direction || position[2] <= 0 || position[2] >= height {
        return Err("SolidCraft hole direction or axial position disagrees with selected face");
    }
    if normal_axis == 0 {
        if position[0] != a[0] || position[1] <= a[1].min(b[1]) || position[1] >= a[1].max(b[1]) {
            return Err("SolidCraft hole center is outside its selected vertical face");
        }
    } else if position[1] != a[1] || position[0] <= a[0].min(b[0]) || position[0] >= a[0].max(b[0])
    {
        return Err("SolidCraft hole center is outside its selected horizontal face");
    }
    let axis_coords = vertices
        .iter()
        .map(|p| p[normal_axis])
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let mut axial_interval = None;
    for pair in axis_coords.windows(2) {
        let (lo, hi) = (pair[0], pair[1]);
        let touches_face = if direction[normal_axis] < 0 {
            hi == position[normal_axis]
        } else {
            lo == position[normal_axis]
        };
        if !touches_face {
            continue;
        }
        let mid = crate::constraint_system::Ratio::new(i128::from(lo) + i128::from(hi), 2)?;
        let fixed = crate::constraint_system::Ratio::new(i128::from(position[1 - normal_axis]), 1)?;
        let sample = if normal_axis == 0 {
            [mid, fixed]
        } else {
            [fixed, mid]
        };
        if crate::sketch_layers::inside(
            &crate::sketch_layers::SketchLoop {
                id: "profile".into(),
                vertices: vertices.to_vec(),
            },
            sample,
        )? && axial_interval.replace([lo, hi]).is_some()
        {
            return Err("ambiguous SolidCraft hole through-run");
        }
    }
    let [axial_lo, axial_hi] =
        axial_interval.ok_or("SolidCraft hole direction does not enter the sketch material")?;
    let axial_center = axial_lo
        .checked_add(axial_hi)
        .ok_or("SolidCraft hole center overflow")?;
    if axial_center % 2 != 0 {
        return Err("SolidCraft hole through-run has a half-millimeter axis center");
    }
    let mut center = position;
    center[normal_axis] = axial_center / 2;
    let diameter = literal_length(f.field("diameter")?)?;
    if diameter <= 0 || diameter % 2 != 0 {
        return Err("SolidCraft hole diameter must define an integer-mm radius");
    }
    let axis_name = ["x", "y", "z"][normal_axis];
    Ok((
        feature_name("hole", feature_id(f)?),
        axis_name.into(),
        center,
        diameter / 2,
    ))
}

/// Translate an archived SolidCraft 0.48.1 L-sketch history into the generic
/// orthogonal sketch feature family, then construct and audit the native model.
pub fn import_sketch_timeline(script: &Json) -> Result<Json, &'static str> {
    script.keys(&["about", "commands"])?;
    let about = script.field("about")?;
    about.keys(&["design", "format", "replay", "version"])?;
    if about.field("format")?.text()? != "solidcraft/66"
        || about.field("version")?.text()? != "0.48.1"
    {
        return Err("unsupported SolidCraft timeline format or version");
    }
    let commands = script.field("commands")?.items()?;
    if commands.len() != 5 {
        return Err("expected sketch, extrusion, fillet and two hole features");
    }
    let sketch_f = feature(&commands[0], "sketch")?;
    let (outer_id, vertices, lines) = sketch_loop(sketch_f)?;
    let sketch_id = feature_id(sketch_f)?;
    let extrude = feature(&commands[1], "extrude")?;
    let height = validate_extrusion(extrude, sketch_id)?;
    if feature_id(extrude)? <= sketch_id {
        return Err("SolidCraft feature ids are not topologically ordered");
    }
    let prefix = format!("F{}:side:", feature_id(extrude)?);
    let fillet_f = feature(&commands[2], "fillet")?;
    let (fillet_id, fillet_vertex, radius) =
        parse_fillet(fillet_f, &vertices, &lines, &prefix, height)?;
    let mut bores = Vec::new();
    let mut used_axes = BTreeSet::new();
    for step in &commands[3..] {
        let hole_f = feature(step, "hole")?;
        let (id, axis, center, radius) = parse_hole(hole_f, &vertices, &lines, &prefix, height)?;
        if !used_axes.insert(axis.clone()) {
            return Err("duplicate SolidCraft hole axis");
        }
        bores.push(object([
            ("id", Json::String(id)),
            ("axis", Json::String(axis)),
            (
                "center_mm",
                Json::Array(
                    center
                        .into_iter()
                        .map(|v| Json::Number(v.to_string()))
                        .collect(),
                ),
            ),
            ("radius_mm", Json::Number(radius.to_string())),
        ]));
    }
    let recipe = object([
        (
            "family",
            Json::String("orthogonal_sketch_features_v1".into()),
        ),
        (
            "outer",
            object([
                ("id", Json::String(outer_id)),
                (
                    "vertices",
                    Json::Array(
                        vertices
                            .iter()
                            .map(|p| {
                                Json::Array(p.iter().map(|v| Json::Number(v.to_string())).collect())
                            })
                            .collect(),
                    ),
                ),
            ]),
        ),
        ("height_mm", Json::Number(height.to_string())),
        (
            "fillets",
            Json::Array(vec![object([
                ("id", Json::String(fillet_id)),
                (
                    "vertex",
                    Json::Array(
                        fillet_vertex
                            .into_iter()
                            .map(|v| Json::Number(v.to_string()))
                            .collect(),
                    ),
                ),
                ("radius_mm", Json::Number(radius.to_string())),
            ])]),
        ),
        ("bores", Json::Array(bores)),
    ]);
    let candidate = interchange::build(&recipe)?;
    if interchange::audit(&candidate)?.field("status")?.text()? != "valid" {
        return Err("SolidCraft sketch candidate failed native audit");
    }
    Ok(candidate)
}
