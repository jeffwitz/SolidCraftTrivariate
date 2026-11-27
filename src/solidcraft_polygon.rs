//! Import a narrowly qualified SolidCraft sketch plus extrusion history into G17/G18.
use crate::json::{object, Json};
use crate::polygon_json::{self, FAMILY, PERFORATED};
use crate::polygon_prism::Vertex;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone)]
struct Curve {
    id: String,
    a: usize,
    b: usize,
}

#[derive(Clone)]
struct ProfileLoop {
    vertices: Vec<Vertex>,
    curve_ids: BTreeSet<String>,
}

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
fn boolean(value: &Json, expected: bool) -> Result<(), &'static str> {
    if value == &Json::Bool(expected) {
        Ok(())
    } else {
        Err("unsupported SolidCraft boolean option")
    }
}
fn exact_i32(value: &Json, allow_string: bool) -> Result<i32, &'static str> {
    let token = match value {
        Json::Number(value) => value.as_str(),
        Json::String(value) if allow_string => value.as_str(),
        _ => return Err("SolidCraft polygon values must be literal numbers"),
    };
    let (n, d) = crate::solidcraft_boxes::decimal_rational(token)?;
    if d != 1 {
        return Err("SolidCraft polygon values must be integral millimetres");
    }
    i32::try_from(n).map_err(|_| "SolidCraft polygon value outside native integer range")
}
fn index(value: &Json, limit: usize) -> Result<usize, &'static str> {
    let n: usize = value.number()?;
    if n >= limit {
        return Err("SolidCraft sketch point index out of range");
    }
    Ok(n)
}
fn valid_name(s: &str) -> bool {
    !s.is_empty()
        && s.len() <= 22
        && s.bytes()
            .all(|c| c.is_ascii_alphanumeric() || c == b'_' || c == b'.')
}
fn native_vertex_id(source: &str) -> Result<String, &'static str> {
    if !valid_name(source) {
        return Err("unsupported SolidCraft sketch point id");
    }
    let id = format!("v_{}", source.replace('.', "_"));
    if id.len() > 24 {
        return Err("SolidCraft sketch point id too long");
    }
    Ok(id)
}
fn feature<'a>(command: &'a Json, kind: &str) -> Result<&'a Json, &'static str> {
    command.keys(&["command", "params"])?;
    if text(command.field("command")?)? != "timeline.add" {
        return Err("unsupported SolidCraft timeline command");
    }
    let params = command.field("params")?;
    params.keys(&["feature"])?;
    let f = params.field("feature")?;
    if text(f.field("type")?)? != kind {
        return Err("unsupported SolidCraft timeline feature order");
    }
    Ok(f)
}
fn feature_id(feature: &Json) -> Result<usize, &'static str> {
    let id: usize = feature.field("id")?.number()?;
    if id == 0 {
        return Err("invalid SolidCraft feature id");
    }
    Ok(id)
}
fn point_inside(point: [i32; 2], vertices: &[Vertex]) -> bool {
    let [px, py] = point.map(|v| 2 * i64::from(v));
    let mut inside = false;
    for i in 0..vertices.len() {
        let a = vertices[i].xy_mm.map(i64::from);
        let b = vertices[(i + 1) % vertices.len()].xy_mm.map(i64::from);
        let cross = (b[0] - a[0]) * (py - 2 * a[1]) - (b[1] - a[1]) * (px - 2 * a[0]);
        if cross == 0
            && 2 * a[0].min(b[0]) <= px
            && px <= 2 * a[0].max(b[0])
            && 2 * a[1].min(b[1]) <= py
            && py <= 2 * a[1].max(b[1])
        {
            return false;
        }
        if (2 * a[1] > py) != (2 * b[1] > py) {
            let lhs = (px - 2 * a[0]) * (b[1] - a[1]);
            let rhs = (py - 2 * a[1]) * (b[0] - a[0]);
            if (b[1] > a[1]) == (lhs < rhs) {
                inside = !inside;
            }
        }
    }
    inside
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

fn parse_loops(sketch_feature: &Json) -> Result<Vec<ProfileLoop>, &'static str> {
    sketch_feature.keys(&["id", "name", "plane", "sketch", "suppressed", "type"])?;
    let sketch_id = feature_id(sketch_feature)?;
    if sketch_feature.field("name")?.text()?.is_empty() {
        return Err("empty SolidCraft sketch name");
    }
    boolean(sketch_feature.field("suppressed")?, false)?;
    let plane = sketch_feature.field("plane")?;
    plane.keys(&["name", "type"])?;
    if text(plane.field("type")?)? != "origin" || text(plane.field("name")?)? != "XY" {
        return Err("SolidCraft polygon sketch must lie on the origin XY plane");
    }
    let sketch = sketch_feature.field("sketch")?;
    sketch.keys(&["constraints", "counters", "curves", "points"])?;
    let points = sketch.field("points")?.items()?;
    let curves_json = sketch.field("curves")?.items()?;
    let constraints = sketch.field("constraints")?.items()?;
    if points.len() < 4 || curves_json.len() < 3 || curves_json.len() > 32 {
        return Err("SolidCraft polygon sketch size out of range");
    }
    let counters = sketch.field("counters")?;
    counters.keys(&["k", "l"])?;
    let kc: usize = counters.field("k")?.number()?;
    let lc: usize = counters.field("l")?.number()?;
    if kc < constraints.len() || lc < curves_json.len() {
        return Err("invalid SolidCraft sketch counters");
    }

    let mut coords = Vec::with_capacity(points.len());
    let mut point_ids = BTreeSet::new();
    for p in points {
        p.keys(&["fixed", "id", "pos"])?;
        if !matches!(p.field("fixed")?, Json::Bool(_)) {
            return Err("invalid SolidCraft point fixed flag");
        }
        let id = text(p.field("id")?)?;
        if id.is_empty() || !point_ids.insert(id.to_owned()) {
            return Err("duplicate SolidCraft sketch point id");
        }
        let pos = p.field("pos")?.items()?;
        if pos.len() != 2 {
            return Err("SolidCraft sketch point must have two coordinates");
        }
        coords.push([exact_i32(&pos[0], false)?, exact_i32(&pos[1], false)?]);
    }
    let mut dsu = DisjointSet::new(points.len());
    let mut constraint_ids = BTreeSet::new();
    for c in constraints {
        c.keys(&["id", "p", "q", "type"])?;
        let id = text(c.field("id")?)?;
        if id.is_empty() || !constraint_ids.insert(id.to_owned()) {
            return Err("duplicate SolidCraft sketch constraint id");
        }
        if text(c.field("type")?)? != "coincident" {
            return Err("unsupported SolidCraft sketch constraint");
        }
        let p = index(c.field("p")?, points.len())?;
        let q = index(c.field("q")?, points.len())?;
        if coords[p] != coords[q] {
            return Err("inconsistent SolidCraft coincident constraint");
        }
        dsu.union(p, q);
    }
    let mut curves = Vec::with_capacity(curves_json.len());
    let mut curve_ids = BTreeSet::new();
    let mut used = BTreeSet::new();
    for c in curves_json {
        c.keys(&["a", "b", "construction", "id", "type"])?;
        if text(c.field("type")?)? != "line" {
            return Err("SolidCraft polygon sketch only supports straight lines");
        }
        boolean(c.field("construction")?, false)?;
        let id = text(c.field("id")?)?;
        if !valid_name(id) || !curve_ids.insert(id.to_owned()) {
            return Err("invalid or duplicate SolidCraft curve id");
        }
        let a = index(c.field("a")?, points.len())?;
        let b = index(c.field("b")?, points.len())?;
        used.insert(a);
        used.insert(b);
        curves.push(Curve {
            id: id.to_owned(),
            a,
            b,
        });
    }
    let used_roots = used.iter().map(|&i| dsu.find(i)).collect::<BTreeSet<_>>();
    for (i, point) in points.iter().enumerate() {
        let root = dsu.find(i);
        if !used_roots.contains(&root) && text(point.field("id")?)? != "origin" {
            return Err("unused SolidCraft sketch point");
        }
    }
    let mut root_coords = BTreeMap::<usize, [i32; 2]>::new();
    let mut root_names = BTreeMap::<usize, String>::new();
    for (i, point) in points.iter().enumerate() {
        let root = dsu.find(i);
        let source_name = text(point.field("id")?)?;
        if let Some(existing) = root_coords.insert(root, coords[i]) {
            if existing != coords[i] {
                return Err("coincident SolidCraft point coordinates differ");
            }
        }
        root_names
            .entry(root)
            .and_modify(|name| {
                if source_name < name {
                    *name = source_name.to_owned();
                }
            })
            .or_insert_with(|| source_name.to_owned());
    }
    let mut graph = BTreeMap::<usize, Vec<(usize, usize)>>::new();
    let mut unique_edges = BTreeSet::new();
    for (ci, c) in curves.iter().enumerate() {
        let a = dsu.find(c.a);
        let b = dsu.find(c.b);
        if a == b || !unique_edges.insert((a.min(b), a.max(b))) {
            return Err("degenerate or duplicate SolidCraft sketch edge");
        }
        graph.entry(a).or_default().push((b, ci));
        graph.entry(b).or_default().push((a, ci));
    }
    if graph.values().any(|neighbors| neighbors.len() != 2) {
        return Err("SolidCraft sketch curves must form closed loops");
    }
    let mut used_edges = BTreeSet::<usize>::new();
    let mut loops = Vec::new();
    while used_edges.len() < curves.len() {
        let start = *graph
            .iter()
            .filter(|(_, edges)| edges.iter().any(|(_, e)| !used_edges.contains(e)))
            .map(|(v, _)| v)
            .min_by_key(|&&root| root_names[&root].clone())
            .ok_or("SolidCraft sketch cycle extraction failed")?;
        let mut cycle = Vec::new();
        let mut edge_set = BTreeSet::new();
        let mut previous = None;
        let mut current = start;
        loop {
            cycle.push(current);
            let next = graph[&current]
                .iter()
                .filter(|(_, e)| Some(*e) != previous && !used_edges.contains(e))
                .min_by_key(|(v, e)| (root_names[v].clone(), curves[*e].id.clone()))
                .copied()
                .ok_or("open SolidCraft sketch cycle")?;
            used_edges.insert(next.1);
            edge_set.insert(curves[next.1].id.clone());
            previous = Some(next.1);
            current = next.0;
            if current == start {
                break;
            }
            if cycle.len() > curves.len() {
                return Err("invalid SolidCraft sketch cycle");
            }
        }
        let vertices = cycle
            .iter()
            .map(|root| {
                Ok(Vertex {
                    id: native_vertex_id(&root_names[root])?,
                    xy_mm: root_coords[root],
                })
            })
            .collect::<Result<Vec<_>, &'static str>>()?;
        if vertices.len() < 3 {
            return Err("SolidCraft sketch loop has fewer than three vertices");
        }
        let area2 = signed_area2(&vertices);
        if area2 == 0 {
            return Err("zero-area SolidCraft sketch loop");
        }
        loops.push(ProfileLoop {
            vertices,
            curve_ids: edge_set,
        });
    }
    if loops.len() > 5 {
        return Err("SolidCraft sketch has too many closed loops");
    }
    let outer_candidates = (0..loops.len())
        .filter(|&i| {
            (0..loops.len()).all(|j| {
                i == j
                    || loops[j]
                        .vertices
                        .iter()
                        .all(|v| point_inside(v.xy_mm, &loops[i].vertices))
            })
        })
        .collect::<Vec<_>>();
    if outer_candidates.len() != 1 {
        return Err("SolidCraft sketch must contain one outer loop and strict holes");
    }
    let outer_i = outer_candidates[0];
    if loops.iter().enumerate().any(|(i, l)| {
        i != outer_i
            && l.vertices
                .iter()
                .any(|v| !point_inside(v.xy_mm, &loops[outer_i].vertices))
    }) {
        return Err("SolidCraft hole lies outside outer loop or touches it");
    }
    let mut ordered = Vec::with_capacity(loops.len());
    ordered.push(loops.remove(outer_i));
    loops.sort_by_key(|l| {
        l.vertices
            .iter()
            .map(|v| v.id.clone())
            .min()
            .unwrap_or_default()
    });
    ordered.extend(loops);
    let _ = sketch_id;
    Ok(ordered)
}

fn validate_extrude(
    extrude: &Json,
    sketch_id: usize,
    loops: &[ProfileLoop],
) -> Result<i32, &'static str> {
    extrude.keys(&[
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
    let _ = feature_id(extrude)?;
    if extrude.field("name")?.text()?.is_empty() || text(extrude.field("type")?)? != "extrude" {
        return Err("invalid SolidCraft extrusion feature");
    }
    boolean(extrude.field("suppressed")?, false)?;
    let reference: usize = extrude.field("sketch")?.number()?;
    if reference != sketch_id {
        return Err("SolidCraft extrusion references another sketch");
    }
    if text(extrude.field("operation")?)? != "new_body" {
        return Err("SolidCraft extrusion must create a new body");
    }
    let bodies = extrude.field("body_names")?.items()?;
    if bodies.len() != 1 || bodies[0].text()?.is_empty() {
        return Err("SolidCraft extrusion must target one new body");
    }
    let params = extrude.field("param_names")?.items()?;
    let mut param_names = BTreeSet::new();
    if params.is_empty()
        || params.iter().any(|p| {
            p.text().map_or(true, |name| {
                name.is_empty() || !param_names.insert(name.to_owned())
            })
        })
    {
        return Err("invalid SolidCraft extrusion parameter names");
    }
    let extent = extrude.field("extent")?;
    extent.keys(&["convention", "direction", "distance"])?;
    if text(extent.field("convention")?)? != "total_distance"
        || text(extent.field("direction")?)? != "positive"
    {
        return Err("unsupported SolidCraft extrusion extent");
    }
    let height = exact_i32(extent.field("distance")?, true)?;
    if !(1..=64).contains(&height) {
        return Err("SolidCraft extrusion height outside native range");
    }
    let profile = extrude.field("profiles")?;
    if loops.len() == 1 && profile.keys(&["type"]).is_ok() && text(profile.field("type")?)? == "all"
    {
        return Ok(height);
    }
    profile.keys(&["loops", "type"])?;
    if text(profile.field("type")?)? != "curves" {
        return Err("ambiguous SolidCraft extrusion profile selection");
    }
    let selected = profile.field("loops")?.items()?;
    if selected.len() != 1 {
        return Err("SolidCraft extrusion must select one material profile");
    }
    let ids = selected[0]
        .items()?
        .iter()
        .map(|id| id.text().map(str::to_owned))
        .collect::<Result<BTreeSet<_>, _>>()?;
    if ids.len() != selected[0].items()?.len() || ids != loops[0].curve_ids {
        return Err("SolidCraft extrusion selection is not the outer profile");
    }
    Ok(height)
}

/// Import an extracted `solidcraft/66` sketch/extrude history as an audited native prism.
pub fn import_polygon_timeline(script: &Json) -> Result<Json, &'static str> {
    let (mut loops, height) = extrusion_profiles(script, false)?;
    let outer = loops.remove(0);
    let recipe = if loops.is_empty() {
        if outer.vertices.len() > 24 {
            return Err("SolidCraft simple polygon exceeds native G17 limit");
        }
        object([
            ("family", Json::String(FAMILY.into())),
            (
                "vertices",
                Json::Array(
                    outer
                        .vertices
                        .into_iter()
                        .map(|v| {
                            object([
                                ("id", Json::String(v.id)),
                                (
                                    "xy_mm",
                                    Json::Array(
                                        v.xy_mm
                                            .into_iter()
                                            .map(|x| Json::Number(x.to_string()))
                                            .collect(),
                                    ),
                                ),
                            ])
                        })
                        .collect(),
                ),
            ),
            ("height_mm", Json::Number(height.to_string())),
            (
                "top_shift_mm",
                Json::Array(vec![Json::Number("0".into()), Json::Number("0".into())]),
            ),
        ])
    } else {
        let holes = loops
            .into_iter()
            .enumerate()
            .map(|(i, l)| {
                object([
                    ("id", Json::String(format!("hole{i}"))),
                    (
                        "vertices",
                        Json::Array(
                            l.vertices
                                .into_iter()
                                .map(|v| {
                                    object([
                                        ("id", Json::String(v.id)),
                                        (
                                            "xy_mm",
                                            Json::Array(
                                                v.xy_mm
                                                    .into_iter()
                                                    .map(|x| Json::Number(x.to_string()))
                                                    .collect(),
                                            ),
                                        ),
                                    ])
                                })
                                .collect(),
                        ),
                    ),
                ])
            })
            .collect();
        object([
            ("family", Json::String(PERFORATED.into())),
            (
                "outer",
                object([
                    ("id", Json::String("outer".into())),
                    (
                        "vertices",
                        Json::Array(
                            outer
                                .vertices
                                .into_iter()
                                .map(|v| {
                                    object([
                                        ("id", Json::String(v.id)),
                                        (
                                            "xy_mm",
                                            Json::Array(
                                                v.xy_mm
                                                    .into_iter()
                                                    .map(|x| Json::Number(x.to_string()))
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
            ("holes", Json::Array(holes)),
            ("height_mm", Json::Number(height.to_string())),
            (
                "top_shift_mm",
                Json::Array(vec![Json::Number("0".into()), Json::Number("0".into())]),
            ),
        ])
    };
    polygon_json::build(&recipe)
}

fn check_preview_coordinates(feature: &Json) -> Result<(), &'static str> {
    for point in feature.field("sketch")?.field("points")?.items()? {
        for value in point.field("pos")?.items()? {
            if !(-4096..=4096).contains(&exact_i32(value, false)?) {
                return Err("sketch preview exceeds planner coordinate range");
            }
        }
    }
    Ok(())
}

pub(crate) fn preview_contours(feature: &Json) -> Result<Vec<Vec<[i32; 2]>>, &'static str> {
    check_preview_coordinates(feature)?;
    Ok(parse_loops(feature)?
        .into_iter()
        .map(|profile| profile.vertices.into_iter().map(|v| v.xy_mm).collect())
        .collect())
}

fn extrusion_profiles(
    script: &Json,
    preview: bool,
) -> Result<(Vec<ProfileLoop>, i32), &'static str> {
    script.keys(&["about", "commands"])?;
    let about = script.field("about")?;
    about.keys(&["design", "format", "replay", "version"])?;
    if text(about.field("format")?)? != "solidcraft/66"
        || text(about.field("version")?)? != "0.48.1"
        || text(about.field("design")?)?.is_empty()
        || text(about.field("replay")?)?.is_empty()
    {
        return Err("unsupported SolidCraft timeline metadata");
    }
    let commands = script.field("commands")?.items()?;
    let first = commands.first().ok_or("empty SolidCraft timeline")?;
    let sketch = feature(first, "sketch")?;
    let sketch_id = feature_id(sketch)?;
    if preview {
        check_preview_coordinates(sketch)?;
    }
    let loops = parse_loops(sketch)?;
    let extrude = match commands {
        [_, extrusion] => feature(extrusion, "extrude")?,
        [_, a, b] if preview => {
            let (overlay, extrusion, prefix) = if feature(a, "sketch").is_ok() {
                (a, b, vec![first.clone()])
            } else {
                (b, a, vec![first.clone(), a.clone()])
            };
            let prefix = object([("commands", Json::Array(prefix))]);
            let expected = crate::solidcraft_preview::descriptor(
                &prefix,
                sketch_id.try_into().map_err(|_| "invalid sketch id")?,
            )?;
            if feature(overlay, "sketch")? != expected.field("feature")? {
                return Err("additional sketch differs from the generated construction preview");
            }
            feature(extrusion, "extrude")?
        }
        _ => return Err("SolidCraft polygon import requires a sketch and extrusion"),
    };
    let height = validate_extrude(extrude, sketch_id, &loops)?;
    Ok((loops, height))
}

/// Extrude the same planar layout shown by the native sketch preview.
pub fn extrude_preview(input: &Json) -> Result<Json, &'static str> {
    let script = crate::solidcraft_preview::snapshot_script(input)?;
    let (profiles, height) = extrusion_profiles(script, true)?;
    let n = |i: i32| Json::Number(i.to_string());
    let mut contours = profiles.into_iter().enumerate().map(|(i, profile)| {
        object([
            ("id", Json::String(format!("loop{i}"))),
            (
                "vertices",
                Json::Array(
                    profile
                        .vertices
                        .iter()
                        .map(|v| Json::Array(v.xy_mm.map(n).to_vec()))
                        .collect(),
                ),
            ),
        ])
    });
    let outer = contours.next().ok_or("empty extrusion profile")?;
    let mut cuts = Vec::new();
    for contour in contours {
        let Json::Object(mut values) = contour else {
            unreachable!()
        };
        values.insert("z_interval_mm".into(), Json::Array(vec![n(0), n(height)]));
        cuts.push(Json::Object(values));
    }
    crate::interchange::build(&object([
        ("family", Json::String("sketch_layer_partition_v1".into())),
        ("outer", outer),
        ("cuts", Json::Array(cuts)),
        ("height_mm", n(height)),
    ]))
}
