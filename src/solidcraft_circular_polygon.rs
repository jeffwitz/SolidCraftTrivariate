//! Translate a SolidCraft sketch with circular voids and one extrusion to G20.
use crate::{
    circular_polygon_json,
    json::{object, Json},
    polygon_prism::Vertex,
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone)]
struct LineCurve {
    id: String,
    a: usize,
    b: usize,
}

#[derive(Clone)]
struct CircleCurve {
    id: String,
    center: usize,
    radius_mm: i32,
}

type BoreDescription = (String, [i32; 2], i32, i32);
type ParsedProfile = (Vec<Vertex>, BTreeSet<String>, Vec<BoreDescription>, usize);

struct DisjointSet(Vec<usize>);
impl DisjointSet {
    fn new(n: usize) -> Self {
        Self((0..n).collect())
    }
    fn find(&mut self, x: usize) -> usize {
        if self.0[x] != x {
            self.0[x] = self.find(self.0[x]);
        }
        self.0[x]
    }
    fn union(&mut self, a: usize, b: usize) {
        let (a, b) = (self.find(a), self.find(b));
        if a != b {
            self.0[b] = a;
        }
    }
}

fn text(value: &Json) -> Result<&str, &'static str> {
    value.text()
}
fn exact_i32(value: &Json, allow_string: bool) -> Result<i32, &'static str> {
    let token = match value {
        Json::Number(s) => s.as_str(),
        Json::String(s) if allow_string => s.as_str(),
        _ => return Err("SolidCraft circular profile values must be literal numbers"),
    };
    // SolidCraft serializes integral dimensions as decimals (including 0.0).
    // The shared box parser intentionally rejects some zero spellings, so
    // recognize a valid all-zero JSON decimal before delegating other values.
    let mantissa = token.split_once(['e', 'E']).map_or(token, |(m, _)| m);
    if mantissa
        .bytes()
        .all(|byte| matches!(byte, b'0' | b'.' | b'-'))
    {
        return Ok(0);
    }
    let (numerator, denominator) = crate::solidcraft_boxes::decimal_rational(token)?;
    if denominator != 1 {
        return Err("SolidCraft circular profile values must be integral millimetres");
    }
    i32::try_from(numerator).map_err(|_| "SolidCraft circular profile value outside integer range")
}
fn index(value: &Json, limit: usize) -> Result<usize, &'static str> {
    let index: usize = value.number()?;
    if index >= limit {
        return Err("SolidCraft sketch point index out of range");
    }
    Ok(index)
}
fn safe_vertex_id(source: &str) -> Result<String, &'static str> {
    if source.is_empty()
        || source.len() > 22
        || !source
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'.')
    {
        return Err("unsupported SolidCraft polygon point id");
    }
    Ok(format!("v_{}", source.replace('.', "_")))
}
fn timeline_feature<'a>(command: &'a Json, kind: &str) -> Result<&'a Json, &'static str> {
    command.keys(&["command", "params"])?;
    if text(command.field("command")?)? != "timeline.add" {
        return Err("unsupported SolidCraft timeline command");
    }
    let params = command.field("params")?;
    params.keys(&["feature"])?;
    let feature = params.field("feature")?;
    if text(feature.field("type")?)? != kind {
        return Err("unsupported SolidCraft feature order");
    }
    Ok(feature)
}
fn feature_id(feature: &Json) -> Result<usize, &'static str> {
    let id: usize = feature.field("id")?.number()?;
    if id == 0 {
        return Err("invalid SolidCraft feature id");
    }
    Ok(id)
}
fn signed_area2(vertices: &[Vertex]) -> i64 {
    (0..vertices.len())
        .map(|i| {
            let a = vertices[i].xy_mm;
            let b = vertices[(i + 1) % vertices.len()].xy_mm;
            i64::from(a[0]) * i64::from(b[1]) - i64::from(a[1]) * i64::from(b[0])
        })
        .sum()
}

fn parse_profile(sketch_feature: &Json) -> Result<ParsedProfile, &'static str> {
    sketch_feature.keys(&["id", "name", "plane", "sketch", "suppressed", "type"])?;
    let sketch_id = feature_id(sketch_feature)?;
    if sketch_feature.field("name")?.text()?.is_empty()
        || sketch_feature.field("suppressed")? != &Json::Bool(false)
    {
        return Err("invalid or suppressed SolidCraft sketch");
    }
    let plane = sketch_feature.field("plane")?;
    plane.keys(&["name", "type"])?;
    if text(plane.field("type")?)? != "origin" || text(plane.field("name")?)? != "XY" {
        return Err("SolidCraft circular profile must lie on origin XY");
    }
    let sketch = sketch_feature.field("sketch")?;
    sketch.keys(&["constraints", "counters", "curves", "points"])?;
    let points = sketch.field("points")?.items()?;
    let curves = sketch.field("curves")?.items()?;
    let constraints = sketch.field("constraints")?.items()?;
    if !(4..=32).contains(&points.len()) || curves.len() < 4 || curves.len() > 36 {
        return Err("SolidCraft circular sketch size out of range");
    }
    let counters = sketch.field("counters")?;
    let nk: usize = if counters.field("k").is_ok() {
        counters.keys(&["c", "k", "l"])?;
        counters.field("k")?.number()?
    } else {
        counters.keys(&["c", "l"])?;
        if !constraints.is_empty() {
            return Err("missing SolidCraft sketch constraint counter");
        }
        0
    };
    let (nc, nl): (usize, usize) = (
        counters.field("c")?.number()?,
        counters.field("l")?.number()?,
    );
    if nk < constraints.len() {
        return Err("invalid SolidCraft sketch constraint counter");
    }
    let mut coordinates = Vec::with_capacity(points.len());
    let mut point_ids = BTreeSet::new();
    for point in points {
        point.keys(&["fixed", "id", "pos"])?;
        if !matches!(point.field("fixed")?, Json::Bool(_)) {
            return Err("invalid SolidCraft point fixed flag");
        }
        let id = text(point.field("id")?)?;
        if id.is_empty() || !point_ids.insert(id.to_owned()) {
            return Err("duplicate SolidCraft point id");
        }
        let pos = point.field("pos")?.items()?;
        if pos.len() != 2 {
            return Err("SolidCraft sketch point must have two coordinates");
        }
        coordinates.push([exact_i32(&pos[0], false)?, exact_i32(&pos[1], false)?]);
    }
    let mut dsu = DisjointSet::new(points.len());
    let mut constraint_ids = BTreeSet::new();
    for constraint in constraints {
        constraint.keys(&["id", "p", "q", "type"])?;
        let id = text(constraint.field("id")?)?;
        if id.is_empty() || !constraint_ids.insert(id.to_owned()) {
            return Err("duplicate SolidCraft constraint id");
        }
        if text(constraint.field("type")?)? != "coincident" {
            return Err("unsupported SolidCraft sketch constraint");
        }
        let p = index(constraint.field("p")?, points.len())?;
        let q = index(constraint.field("q")?, points.len())?;
        if coordinates[p] != coordinates[q] {
            return Err("inconsistent SolidCraft coincident constraint");
        }
        dsu.union(p, q);
    }
    let mut point_names = BTreeMap::<usize, String>::new();
    let mut point_coordinates = BTreeMap::<usize, [i32; 2]>::new();
    for (i, point) in points.iter().enumerate() {
        let root = dsu.find(i);
        let source_name = text(point.field("id")?)?;
        if let Some(existing) = point_coordinates.insert(root, coordinates[i]) {
            if existing != coordinates[i] {
                return Err("coincident point coordinates differ");
            }
        }
        point_names
            .entry(root)
            .and_modify(|name| {
                if source_name < name {
                    *name = source_name.to_owned();
                }
            })
            .or_insert_with(|| source_name.to_owned());
    }
    let mut line_curves = Vec::new();
    let mut circles = Vec::new();
    let mut curve_ids = BTreeSet::new();
    let mut used_roots = BTreeSet::new();
    for curve in curves {
        let kind = text(curve.field("type")?)?;
        let id = text(curve.field("id")?)?;
        if id.is_empty() || !curve_ids.insert(id.to_owned()) {
            return Err("duplicate SolidCraft curve id");
        }
        if curve.field("construction")? != &Json::Bool(false) {
            return Err("construction geometry is outside SC6");
        }
        match kind {
            "line" => {
                curve.keys(&["a", "b", "construction", "id", "type"])?;
                if !id
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'.')
                {
                    return Err("invalid SolidCraft line id");
                }
                let a = dsu.find(index(curve.field("a")?, points.len())?);
                let b = dsu.find(index(curve.field("b")?, points.len())?);
                if a == b {
                    return Err("degenerate SolidCraft profile edge");
                }
                used_roots.insert(a);
                used_roots.insert(b);
                line_curves.push(LineCurve {
                    id: id.to_owned(),
                    a,
                    b,
                });
            }
            "circle" => {
                curve.keys(&["c", "construction", "id", "r", "type"])?;
                let center = dsu.find(index(curve.field("c")?, points.len())?);
                used_roots.insert(center);
                let radius_mm = exact_i32(curve.field("r")?, false)?;
                if radius_mm < 1 {
                    return Err("circle radius must be a positive integral millimetre");
                }
                circles.push(CircleCurve {
                    id: id.to_owned(),
                    center,
                    radius_mm,
                });
            }
            _ => return Err("SC6 supports only straight lines and circle sketch curves"),
        }
    }
    if circles.is_empty() || circles.len() > 4 || nc < circles.len() {
        return Err("SC6 requires one to four stored circle curves");
    }
    if nl < line_curves.len() {
        return Err("invalid SolidCraft sketch curve counter");
    }
    for (i, point) in points.iter().enumerate() {
        let root = dsu.find(i);
        if !used_roots.contains(&root) && text(point.field("id")?)? != "origin" {
            return Err("unused SolidCraft sketch point");
        }
    }
    let mut adjacency = BTreeMap::<usize, Vec<(usize, usize)>>::new();
    let mut unique_edges = BTreeSet::new();
    for (ci, line) in line_curves.iter().enumerate() {
        if !unique_edges.insert((line.a.min(line.b), line.a.max(line.b))) {
            return Err("duplicate SolidCraft profile edge");
        }
        adjacency.entry(line.a).or_default().push((line.b, ci));
        adjacency.entry(line.b).or_default().push((line.a, ci));
    }
    if adjacency.values().any(|neighbors| neighbors.len() != 2) {
        return Err("outer SolidCraft profile must be one closed line cycle");
    }
    let mut consumed = BTreeSet::new();
    let mut cycle_roots = Vec::new();
    let mut selected_lines = BTreeSet::new();
    let start = *adjacency
        .keys()
        .min_by_key(|root| point_names[root].clone())
        .ok_or("missing polygon profile")?;
    let mut current = start;
    let mut previous = None;
    loop {
        cycle_roots.push(current);
        let next = adjacency[&current]
            .iter()
            .filter(|(_, edge)| Some(*edge) != previous)
            .min_by_key(|(root, edge)| (point_names[root].clone(), line_curves[*edge].id.clone()))
            .copied()
            .ok_or("open polygon profile")?;
        if !consumed.insert(next.1) {
            if next.0 == start {
                break;
            }
            return Err("invalid polygon curve graph");
        }
        selected_lines.insert(line_curves[next.1].id.clone());
        previous = Some(next.1);
        current = next.0;
        if current == start {
            break;
        }
        if cycle_roots.len() > line_curves.len() {
            return Err("invalid polygon profile cycle");
        }
    }
    if consumed.len() != line_curves.len() || cycle_roots.len() < 3 {
        return Err("SC6 requires exactly one simple outer line contour");
    }
    let vertices = cycle_roots
        .iter()
        .map(|root| {
            Ok(Vertex {
                id: safe_vertex_id(&point_names[root])?,
                xy_mm: point_coordinates[root],
            })
        })
        .collect::<Result<Vec<_>, &'static str>>()?;
    if signed_area2(&vertices) == 0 {
        return Err("zero-area SolidCraft outer profile");
    }
    let bores = circles
        .into_iter()
        .map(|circle| {
            let center = point_coordinates[&circle.center];
            let guard_mm = circle
                .radius_mm
                .checked_mul(2)
                .ok_or("circular bore guard overflow")?;
            Ok((circle.id, center, circle.radius_mm, guard_mm))
        })
        .collect::<Result<Vec<_>, &'static str>>()?;
    let _ = nl;
    Ok((vertices, selected_lines, bores, sketch_id))
}

fn validate_extrude(
    feature: &Json,
    sketch_id: usize,
    line_ids: &BTreeSet<String>,
) -> Result<i32, &'static str> {
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
    feature_id(feature)?;
    if text(feature.field("type")?)? != "extrude"
        || feature.field("suppressed")? != &Json::Bool(false)
        || feature.field("name")?.text()?.is_empty()
    {
        return Err("invalid SolidCraft extrusion feature");
    }
    if feature.field("sketch")?.number::<usize>()? != sketch_id {
        return Err("SolidCraft extrusion references another sketch");
    }
    if text(feature.field("operation")?)? != "new_body" {
        return Err("SolidCraft circular extrusion must create a body");
    }
    let bodies = feature.field("body_names")?.items()?;
    if bodies.len() != 1 || bodies[0].text()?.is_empty() {
        return Err("SolidCraft circular extrusion must target one new body");
    }
    let params = feature.field("param_names")?.items()?;
    let mut names = BTreeSet::new();
    if params.is_empty()
        || params.iter().any(|p| {
            p.text()
                .map_or(true, |s| s.is_empty() || !names.insert(s.to_owned()))
        })
    {
        return Err("invalid SolidCraft extrusion parameter names");
    }
    let extent = feature.field("extent")?;
    extent.keys(&["convention", "direction", "distance"])?;
    if text(extent.field("convention")?)? != "total_distance"
        || text(extent.field("direction")?)? != "positive"
    {
        return Err("unsupported SolidCraft extrusion direction");
    }
    let height = exact_i32(extent.field("distance")?, true)?;
    if !(1..=64).contains(&height) {
        return Err("SolidCraft extrusion height outside native range");
    }
    let profiles = feature.field("profiles")?;
    profiles.keys(&["loops", "type"])?;
    if text(profiles.field("type")?)? != "curves" {
        return Err("circular profile selection must be explicit");
    }
    let loops = profiles.field("loops")?.items()?;
    if loops.len() != 1 {
        return Err("SolidCraft extrusion must select one outer contour");
    }
    let selected = loops[0]
        .items()?
        .iter()
        .map(|v| v.text().map(str::to_owned))
        .collect::<Result<BTreeSet<_>, _>>()?;
    if selected.len() != loops[0].items()?.len() || selected != *line_ids {
        return Err("SolidCraft extrusion selection differs from outer contour");
    }
    Ok(height)
}

/// Import a canonical two-feature SolidCraft circle-sketch extrusion as G20.
pub fn import_circular_polygon_timeline(script: &Json) -> Result<Json, &'static str> {
    script.keys(&["about", "commands"])?;
    let about = script.field("about")?;
    about.keys(&["design", "format", "replay", "version"])?;
    if text(about.field("design")?)?.is_empty()
        || text(about.field("format")?)? != "solidcraft/66"
        || text(about.field("version")?)? != "0.48.1"
        || text(about.field("replay")?)?.is_empty()
    {
        return Err("unsupported SolidCraft history metadata");
    }
    let commands = script.field("commands")?.items()?;
    if commands.len() != 2 {
        return Err("SC6 accepts exactly one sketch and one extrusion");
    }
    let sketch = timeline_feature(&commands[0], "sketch")?;
    let (vertices, line_ids, circles, sketch_id) = parse_profile(sketch)?;
    let extrude = timeline_feature(&commands[1], "extrude")?;
    let height = validate_extrude(extrude, sketch_id, &line_ids)?;
    let bores = circles
        .into_iter()
        .map(|(id, center_mm, radius_mm, guard_mm)| {
            object([
                ("id", Json::String(id)),
                (
                    "center_mm",
                    Json::Array(
                        center_mm
                            .into_iter()
                            .map(|v| Json::Number(v.to_string()))
                            .collect(),
                    ),
                ),
                ("radius_mm", Json::Number(radius_mm.to_string())),
                ("guard_mm", Json::Number(guard_mm.to_string())),
            ])
        })
        .collect();
    let recipe = object([
        ("family", Json::String("polygon_circular_bores_v1".into())),
        (
            "outer",
            object([
                ("id", Json::String("outer".into())),
                (
                    "vertices",
                    Json::Array(
                        vertices
                            .into_iter()
                            .map(|v| {
                                object([
                                    ("id", Json::String(v.id)),
                                    (
                                        "xy_mm",
                                        Json::Array(
                                            v.xy_mm
                                                .into_iter()
                                                .map(|n| Json::Number(n.to_string()))
                                                .collect(),
                                        ),
                                    ),
                                ])
                            })
                            .collect(),
                    ),
                ),
            ]),
        ),
        ("bores", Json::Array(bores)),
        ("height_mm", Json::Number(height.to_string())),
        (
            "top_shift_mm",
            Json::Array(vec![Json::Number("0".into()), Json::Number("0".into())]),
        ),
    ]);
    circular_polygon_json::build(&recipe)
}
