//! Adapt SolidCraft edge and face selections to the stock-corner kernel or host blend.
//! The host owns the B-rep topology and executes the edit as a timeline transaction.
//! Tangent propagation is configurable for host blends across split rims.

use super::*;
use std::collections::{BTreeSet, HashMap};

#[derive(Clone)]
struct EdgeRow {
    index: u64,
    name: Option<String>,
    length: f64,
    start: Option<[f64; 3]>,
    end: Option<[f64; 3]>,
    mid: Option<[f64; 3]>,
}

#[derive(Clone)]
struct FaceRow {
    index: u64,
    name: Option<String>,
    point: Option<[f64; 3]>,
    normal: Option<[f64; 3]>,
}

fn bad(field: &str, reason: impl ToString) -> PluginError {
    invalid(field, reason)
}

fn vec3(v: &Value, field: &str) -> Result<[f64; 3], PluginError> {
    let a = v
        .as_array()
        .filter(|a| a.len() == 3)
        .ok_or_else(|| bad(field, "point must be a 3D coordinate"))?;
    let mut p = [0.0; 3];
    for (i, x) in a.iter().enumerate() {
        p[i] = x
            .as_f64()
            .filter(|v| v.is_finite())
            .ok_or_else(|| bad(field, "point coordinates must be finite numbers"))?;
    }
    Ok(p)
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let den = ab.iter().map(|x| x * x).sum::<f64>();
    let t = if den > 0.0 {
        (ap.iter().zip(ab).map(|(x, y)| x * y).sum::<f64>() / den).clamp(0.0, 1.0)
    } else {
        0.0
    };
    distance(p, [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]])
}

fn one_body(selection: &[Value]) -> Result<String, PluginError> {
    if selection.is_empty() {
        return Err(bad("selection", "select at least one edge or face"));
    }
    let mut owner: Option<&str> = None;
    for item in selection {
        match item.get("type").and_then(Value::as_str) {
            Some("edge" | "face") => {}
            Some("body") => {
                return Err(bad("selection", "select edges or faces, not a whole body"));
            }
            _ => return Err(bad("selection", "selection items must be edges or faces")),
        }
        let body = item
            .get("body")
            .and_then(Value::as_str)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| {
                bad(
                    "selection",
                    "each selected edge or face must name its owning body",
                )
            })?;
        if owner.is_some_and(|current| current != body) {
            return Err(bad(
                "selection",
                "select edges or faces from one body at a time",
            ));
        }
        owner = Some(body);
    }
    Ok(owner.expect("nonempty selection").to_string())
}

fn body_exists(body: &str) -> Result<(), PluginError> {
    let raw = commands::run("model.body_features", "{}")?;
    let data: Value = serde_json::from_str(&raw).map_err(|e| bad("selection", e))?;
    let found = data
        .get("bodies")
        .and_then(Value::as_array)
        .is_some_and(|rows| {
            rows.iter()
                .any(|row| row.get("body").and_then(Value::as_str) == Some(body))
        });
    if !found {
        return Err(bad(
            "selection",
            format!("body `{body}` is not present in this design"),
        ));
    }
    Ok(())
}

fn edges(body: &str) -> Result<Vec<EdgeRow>, PluginError> {
    let raw = commands::run("model.edges", &json!({"body":body}).to_string())?;
    let data: Value = serde_json::from_str(&raw).map_err(|e| bad("selection", e))?;
    let rows = data
        .get("edges")
        .and_then(Value::as_array)
        .ok_or_else(|| bad("selection", "SolidCraft did not return the body's edges"))?;
    rows.iter()
        .map(|e| {
            Ok(EdgeRow {
                index: e.get("index").and_then(Value::as_u64).ok_or_else(|| {
                    bad("selection", "SolidCraft returned an edge without an index")
                })?,
                name: e.get("name").and_then(Value::as_str).map(str::to_string),
                length: e
                    .get("length")
                    .and_then(Value::as_f64)
                    .filter(|x| x.is_finite() && *x >= 0.0)
                    .ok_or_else(|| {
                        bad(
                            "selection",
                            "SolidCraft returned an edge with invalid length",
                        )
                    })?,
                start: e
                    .get("start")
                    .filter(|v| !v.is_null())
                    .map(|v| vec3(v, "selection"))
                    .transpose()?,
                end: e
                    .get("end")
                    .filter(|v| !v.is_null())
                    .map(|v| vec3(v, "selection"))
                    .transpose()?,
                mid: e
                    .get("mid")
                    .filter(|v| !v.is_null())
                    .map(|v| vec3(v, "selection"))
                    .transpose()?,
            })
        })
        .collect()
}

fn faces(body: &str) -> Result<Vec<FaceRow>, PluginError> {
    let raw = commands::run("model.faces", &json!({"body":body}).to_string())?;
    let data: Value = serde_json::from_str(&raw).map_err(|e| bad("selection", e))?;
    let rows = data
        .get("faces")
        .and_then(Value::as_array)
        .ok_or_else(|| bad("selection", "SolidCraft did not return the body's faces"))?;
    rows.iter()
        .map(|f| {
            Ok(FaceRow {
                index: f.get("index").and_then(Value::as_u64).ok_or_else(|| {
                    bad("selection", "SolidCraft returned a face without an index")
                })?,
                name: f.get("name").and_then(Value::as_str).map(str::to_string),
                point: f
                    .get("point")
                    .filter(|v| !v.is_null())
                    .map(|v| vec3(v, "selection"))
                    .transpose()?,
                normal: f
                    .get("plane_normal")
                    .filter(|v| !v.is_null())
                    .map(|v| vec3(v, "selection"))
                    .transpose()?,
            })
        })
        .collect()
}

fn point_witness(item: &Value) -> Result<[f64; 3], PluginError> {
    vec3(
        item.get("point").ok_or_else(|| {
            bad(
                "selection",
                "each selected edge or face needs a point witness",
            )
        })?,
        "selection",
    )
}

fn verify_face_witness(body: &str, face_name: &str, point: [f64; 3]) -> Result<(), PluginError> {
    let raw = commands::run(
        "inspect.classify_points",
        &json!({"body":body,"face":face_name,"points":[point]}).to_string(),
    )?;
    let data: Value = serde_json::from_str(&raw).map_err(|e| bad("selection", e))?;
    let result = data.pointer("/results/0").ok_or_else(|| {
        bad(
            "selection",
            "SolidCraft did not classify the selected face witness",
        )
    })?;
    if result.get("class").and_then(Value::as_str) != Some("on")
        || result.get("face").and_then(Value::as_str) != Some(face_name)
    {
        return Err(bad(
            "selection",
            format!("point witness does not lie on selected face `{face_name}`"),
        ));
    }
    Ok(())
}

fn resolve_edge(item: &Value, edges: &[EdgeRow], body: &str) -> Result<u64, PluginError> {
    let index = item.get("index").and_then(Value::as_u64).ok_or_else(|| {
        bad(
            "selection",
            "each selected edge needs its model.edges index",
        )
    })?;
    let edge = edges.iter().find(|e| e.index == index).ok_or_else(|| {
        bad(
            "selection",
            format!("edge index {index} is not present on `{body}`"),
        )
    })?;
    if let (Some(point), Some(a), Some(b)) = (
        item.get("point")
            .map(|v| vec3(v, "selection"))
            .transpose()?,
        edge.start,
        edge.end,
    ) {
        let chord = distance(a, b);
        let tolerance = (edge.length.max(chord).max(1.0) * 1e-8).max(1e-8);
        // The midpoint is not a projection oracle for a curved edge. Check a witness only
        // when the host's measured length certifies that this edge is straight.
        if (edge.length - chord).abs() <= tolerance && segment_distance(point, a, b) > tolerance {
            return Err(bad(
                "selection",
                format!("point witness does not lie on straight edge {index} of `{body}`"),
            ));
        }
    }
    Ok(index)
}

fn coplanar(a: &FaceRow, b: &FaceRow, tolerance: f64) -> bool {
    let (Some(na), Some(nb), Some(pa), Some(pb)) = (a.normal, b.normal, a.point, b.point) else {
        return false;
    };
    let dot = na.iter().zip(nb).map(|(x, y)| x * y).sum::<f64>();
    if dot < 1.0 - 1e-10 {
        return false;
    }
    let delta = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
    delta.iter().zip(na).map(|(x, n)| x * n).sum::<f64>().abs() <= tolerance
}

fn edge_has_face(name: &str, face: &str) -> bool {
    // Edge names append `~N` only when a face pair has multiple pieces.
    name.split('|').any(|part| {
        let canonical = part
            .rsplit_once('~')
            .filter(|(_, suffix)| !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_digit()))
            .map_or(part, |(base, _)| base);
        canonical == face
    })
}

fn edge_face_sides(name: &str) -> Vec<&str> {
    name.split('|')
        .map(|part| {
            part.rsplit_once('~')
                .filter(|(_, suffix)| {
                    !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_digit())
                })
                .map_or(part, |(base, _)| base)
        })
        .collect()
}

type BodyMetrics = solidcraft_trivariate::stock_corner::BodyMetrics;

fn find_named<'a>(value: &'a Value, name: &str) -> Option<&'a Value> {
    if value.get("name").and_then(Value::as_str) == Some(name) {
        return Some(value);
    }
    match value {
        Value::Array(rows) => rows.iter().find_map(|row| find_named(row, name)),
        Value::Object(rows) => rows.values().find_map(|row| find_named(row, name)),
        _ => None,
    }
}

fn body_metrics(body: &str) -> Result<BodyMetrics, PluginError> {
    let raw = commands::run("inspect.measure", &json!({"bodies":[body]}).to_string())?;
    let data: Value = serde_json::from_str(&raw).map_err(|e| bad("selection", e))?;
    let row = find_named(&data, body)
        .ok_or_else(|| bad("selection", "SolidCraft did not measure the selected body"))?;
    // Current host releases wrap measurements in `outcome`; accepting both shapes keeps the
    // plugin compatible with the established solidcraft-cli and the embedded host adapter.
    let measured = row.get("outcome").unwrap_or(row);
    let bounds = measured.get("bbox").ok_or_else(|| {
        bad(
            "selection",
            "SolidCraft did not return selected-body bounds",
        )
    })?;
    let min = vec3(
        bounds
            .get("min")
            .ok_or_else(|| bad("selection", "missing bbox minimum"))?,
        "selection",
    )?;
    let max = vec3(
        bounds
            .get("max")
            .ok_or_else(|| bad("selection", "missing bbox maximum"))?,
        "selection",
    )?;
    let volume = measured
        .get("volume_mm3")
        .and_then(Value::as_f64)
        .filter(|v| v.is_finite() && *v > 0.0)
        .ok_or_else(|| bad("selection", "selected body is not a measured closed solid"))?;
    let closed = measured
        .get("closed")
        .and_then(Value::as_bool)
        .unwrap_or(false);
    let shells = measured
        .pointer("/kernel/shells")
        .and_then(Value::as_u64)
        .unwrap_or(0) as usize;
    if !closed || shells != 1 || (0..3).any(|i| max[i] <= min[i]) {
        return Err(bad(
            "selection",
            "analytic fillet requires a closed body with finite bounds",
        ));
    }
    Ok(BodyMetrics {
        volume,
        closed,
        shell_count: shells,
        bounds: (min, max),
    })
}

fn local_corner_cut(
    body: &str,
    plan: &solidcraft_trivariate::stock_corner::StockCornerPlan,
    before: BodyMetrics,
    field: &str,
) -> Result<String, PluginError> {
    use solidcraft_trivariate::stock_corner::SketchSegment;
    let sketch_raw = commands::run(
        "sketch.create",
        &json!({
            "name":"SCT Local Corner Treatment",
            "plane":{"origin":plan.origin,"x_dir":plan.x_dir,"y_dir":plan.y_dir}
        })
        .to_string(),
    )?;
    let sketch: Value = serde_json::from_str(&sketch_raw).map_err(|e| bad(field, e))?;
    let sketch_id = sketch
        .get("sketch")
        .and_then(Value::as_u64)
        .ok_or_else(|| bad(field, "host did not return the temporary sketch feature id"))?;
    let built = (|| -> Result<String, PluginError> {
        let mut i = 0;
        while i < plan.profile_segments.len() {
            match &plan.profile_segments[i] {
                SketchSegment::Line { start, end } => {
                    let mut points = vec![*start, *end];
                    i += 1;
                    while let Some(SketchSegment::Line { start, end }) =
                        plan.profile_segments.get(i)
                    {
                        let previous = points.last().copied().expect("at least one point");
                        if previous != *start {
                            return Err(bad(
                                "selection",
                                "core profile has disconnected line segments",
                            ));
                        }
                        points.push(*end);
                        i += 1;
                    }
                    commands::run(
                        "sketch.line",
                        &json!({"points":points,"closed":false}).to_string(),
                    )?;
                }
                SketchSegment::Arc {
                    center,
                    start,
                    end,
                    clockwise,
                } => {
                    commands::run(
                        "sketch.arc.center_point",
                        &json!({"center":center,"start":start,"end":end,"clockwise":clockwise})
                            .to_string(),
                    )?;
                    i += 1;
                }
            }
        }
        commands::run("sketch.finish", "{}")?;
        commands::run("solid.extrude", &json!({"sketch":sketch_id,"distance":plan.tool_length,"direction":"positive","operation":"cut","targets":[body]}).to_string())
    })();
    let reply = match built {
        Ok(reply) => reply,
        Err(error) => {
            let _ = commands::run("timeline.delete", &json!({"feature":sketch_id}).to_string());
            return Err(error);
        }
    };
    let after = match body_metrics(body) {
        Ok(m) => m,
        Err(error) => {
            let _ = commands::run("timeline.delete", &json!({"feature":sketch_id}).to_string());
            return Err(error);
        }
    };
    if let Err(reason) = solidcraft_trivariate::stock_corner::validate_result(before, after, plan) {
        let _ = commands::run("timeline.delete", &json!({"feature":sketch_id}).to_string());
        return Err(bad("selection", reason));
    }
    Ok(reply)
}

fn face_boundary_edges(
    selected: &[Value],
    faces: &[FaceRow],
    edges: &[EdgeRow],
    body: &str,
) -> Result<Vec<u64>, PluginError> {
    let face_by_index: HashMap<u64, &FaceRow> = faces.iter().map(|f| (f.index, f)).collect();
    let name_to_face: HashMap<&str, &FaceRow> = faces
        .iter()
        .filter_map(|f| f.name.as_deref().map(|n| (n, f)))
        .collect();
    let mut group = BTreeSet::<String>::new();
    let mut todo = Vec::new();
    for item in selected {
        let index = item.get("index").and_then(Value::as_u64).ok_or_else(|| {
            bad(
                "selection",
                "each selected face needs its model.faces index",
            )
        })?;
        let face = face_by_index.get(&index).ok_or_else(|| {
            bad(
                "selection",
                format!("face index {index} is not present on `{body}`"),
            )
        })?;
        let name = face.name.as_deref().ok_or_else(|| {
            bad(
                "selection",
                format!("face {index} on `{body}` has no persistent name"),
            )
        })?;
        let p = point_witness(item)?;
        verify_face_witness(body, name, p)?;
        if group.insert(name.to_string()) {
            todo.push(name.to_string());
        }
    }
    let local_scale = edges.iter().map(|e| e.length).fold(0.0, f64::max);
    let tol = (local_scale * 1e-8).max(1e-9);
    while let Some(name) = todo.pop() {
        let Some(face) = name_to_face.get(name.as_str()) else {
            continue;
        };
        for edge in edges {
            let Some(edge_name) = edge.name.as_deref() else {
                continue;
            };
            if !edge_has_face(edge_name, &name) {
                continue;
            }
            // Coplanar faces sharing a named edge represent one selected planar region; omit
            // their internal partition seam from the native fillet/chamfer edge set.
            for other in faces {
                let Some(other_name) = other.name.as_deref() else {
                    continue;
                };
                if other_name != name
                    && edge_has_face(edge_name, other_name)
                    && coplanar(face, other, tol)
                    && group.insert(other_name.to_string())
                {
                    todo.push(other_name.to_string());
                }
            }
        }
    }
    let selected_faces = faces
        .iter()
        .filter(|f| f.name.as_ref().is_some_and(|n| group.contains(n)))
        .count();
    if selected_faces == faces.len() {
        return Err(bad(
            "selection",
            "select a face boundary, not the whole body",
        ));
    }
    let mut result = BTreeSet::new();
    for e in edges {
        let Some(name) = e.name.as_deref() else {
            continue;
        };
        let sides = edge_face_sides(name);
        let in_group = sides.iter().filter(|n| group.contains(**n)).count();
        if sides.len() == 2 && in_group == 1 {
            result.insert(e.index);
        }
    }
    if result.is_empty() {
        return Err(bad(
            "selection",
            "selected face has no resolvable boundary edges",
        ));
    }
    if result.len() == edges.len() {
        return Err(bad(
            "selection",
            "selected faces resolve to the whole body boundary",
        ));
    }
    Ok(result.into_iter().collect())
}

/// Convert host-owned selection and topology records into the geometry core's input.
/// Missing B-Rep witnesses make the local sweep ineligible, so the caller can preserve the
/// host's ordinary feature path.
fn stock_corner_plan(
    body: &str,
    selected_index: u64,
    edges: &[EdgeRow],
    metrics: BodyMetrics,
    amount: f64,
    chamfer: bool,
) -> Result<solidcraft_trivariate::stock_corner::StockCornerPlan, &'static str> {
    use solidcraft_trivariate::stock_corner::{Edge, Face, Input, Treatment, analyze};

    let selected = edges
        .iter()
        .find(|edge| edge.index == selected_index)
        .ok_or("selected edge is absent")?;
    let (Some(start), Some(end), Some(name)) =
        (selected.start, selected.end, selected.name.as_deref())
    else {
        return Err("selected edge lacks its B-Rep endpoints or adjacent face names");
    };
    let sides = edge_face_sides(name);
    if sides.len() != 2 || sides[0] == sides[1] {
        return Err("selected edge does not identify two distinct adjacent faces");
    }
    let face_rows = faces(body).map_err(|_| "could not read adjacent faces")?;
    let host_face = |side: &str| {
        face_rows
            .iter()
            .find(|face| face.name.as_deref() == Some(side))
    };
    let to_face = |face: &FaceRow| -> Result<Face, &'static str> {
        Ok(Face {
            normal: face.normal.ok_or("adjacent face has no plane normal")?,
            point: face.point.ok_or("adjacent face has no plane point")?,
        })
    };
    let adjacent_faces = [
        to_face(host_face(sides[0]).ok_or("first adjacent face is absent")?)?,
        to_face(host_face(sides[1]).ok_or("second adjacent face is absent")?)?,
    ];
    let mut support_trims = Vec::new();
    for edge in edges.iter().filter(|edge| edge.index != selected_index) {
        let Some(edge_name) = edge.name.as_deref() else {
            continue;
        };
        if !sides.iter().any(|side| edge_has_face(edge_name, side)) {
            continue;
        }
        let (Some(start), Some(end)) = (edge.start, edge.end) else {
            return Err("a support trim is missing endpoints");
        };
        support_trims.push(Edge {
            length: edge.length,
            start,
            end,
            mid: edge.mid,
        });
    }
    let input = Input {
        body: metrics,
        selected: Edge {
            length: selected.length,
            start,
            end,
            mid: selected.mid,
        },
        adjacent_faces,
        support_trims,
    };
    let treatment = if chamfer {
        Treatment::EqualDistanceChamfer(amount)
    } else {
        Treatment::Radius(amount)
    };
    analyze(&input, treatment)
}

pub(super) fn apply(entry: &str, params: &Value) -> Result<String, PluginError> {
    let field = if entry == "fillet" {
        "radius_mm"
    } else {
        "distance_mm"
    };
    let amount = params
        .get(field)
        .map(|v| {
            v.as_f64()
                .ok_or_else(|| bad(field, "must be a positive finite number of millimetres"))
        })
        .transpose()?
        .unwrap_or(0.25);
    if !amount.is_finite() || amount <= 0.0 {
        return Err(bad(
            field,
            "must be a positive finite number of millimetres",
        ));
    }
    let selection = params
        .get("selection")
        .and_then(Value::as_array)
        .ok_or_else(|| bad("selection", "expected a nonempty selection array"))?;
    let body = one_body(selection)?;
    body_exists(&body)?;
    let edges = edges(&body)?;
    let selected_edges: Vec<u64> = if selection.iter().any(|s| s["type"] == "face") {
        if selection.iter().any(|s| s["type"] == "edge") {
            return Err(bad(
                "selection",
                "select either edges or faces in one operation",
            ));
        }
        face_boundary_edges(selection, &faces(&body)?, &edges, &body)?
    } else {
        selection
            .iter()
            .map(|item| resolve_edge(item, &edges, &body))
            .collect::<Result<_, _>>()?
    };
    let commands_edges: Vec<Value> = selected_edges
        .iter()
        .map(|index| json!({"body":body,"index":index}))
        .collect();
    let command = if entry == "fillet" {
        "solid.fillet"
    } else {
        "solid.chamfer"
    };
    let tangent_chain = params
        .get("tangent_chain")
        .and_then(Value::as_bool)
        .unwrap_or(solidcraft_trivariate::stock_corner::DEFAULT_TANGENT_CHAIN);
    let mut args = json!({"body":body,"edges":commands_edges,"tangent_chain":tangent_chain});
    if entry == "fillet" {
        args["radius"] = json!(amount);
        args["fit"] = json!("exact");
    } else {
        args["distance"] = json!(amount);
    }
    if selected_edges.len() == 1
        && selection.len() == 1
        && selection[0].get("type").and_then(Value::as_str) == Some("edge")
    {
        // The core analyzes the actual adjacent planes and emits the local section. Use the
        // swept profile only when those planes contain a nearby trim; regular cases continue
        // through the host's native edge operation.
        if let Ok(before) = body_metrics(&body) {
            if let Ok(plan) = stock_corner_plan(
                &body,
                selected_edges[0],
                &edges,
                before,
                amount,
                entry == "chamfer",
            ) {
                if plan.nearby_support_trim {
                    return local_corner_cut(&body, &plan, before, field);
                }
            }
        }
    }
    commands::run(command, &args.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejects_empty_body_mixed_owner_and_malformed_selections() {
        assert!(
            one_body(&[])
                .unwrap_err()
                .to_string()
                .contains("select at least one")
        );
        let whole = json!([{"type":"body","body":"B"}]);
        assert!(
            one_body(whole.as_array().unwrap())
                .unwrap_err()
                .to_string()
                .contains("not a whole body")
        );
        let mixed = json!([{"type":"edge","body":"A"},{"type":"face","body":"B"}]);
        assert!(
            one_body(mixed.as_array().unwrap())
                .unwrap_err()
                .to_string()
                .contains("one body")
        );
    }
    #[test]
    fn point_witness_must_match_sparse_edge_index() {
        let edges = [EdgeRow {
            index: 115,
            name: Some("A|B".into()),
            length: 2.0,
            start: Some([0.0, 0.0, 0.0]),
            end: Some([2.0, 0.0, 0.0]),
            mid: Some([1.0, 0.0, 0.0]),
        }];
        let item = json!({"type":"edge","body":"B","index":115,"point":[0.5,0.0,0.0]});
        assert_eq!(resolve_edge(&item, &edges, "B").unwrap(), 115);
        let bad = json!({"type":"edge","body":"B","index":115,"point":[0.5,2.0,0.0]});
        assert!(resolve_edge(&bad, &edges, "B").is_err());
    }

    #[test]
    fn face_pair_names_keep_duplicate_sides_and_numbered_pieces() {
        assert_eq!(edge_face_sides("F1|F2~12"), vec!["F1", "F2"]);
        assert_eq!(edge_face_sides("F1|F1"), vec!["F1", "F1"]);
        assert!(coplanar(
            &FaceRow {
                index: 0,
                name: None,
                point: Some([0.0, 0.0, 0.0]),
                normal: Some([0.0, 0.0, 1.0])
            },
            &FaceRow {
                index: 1,
                name: None,
                point: Some([10000.0, 0.0, 0.0]),
                normal: Some([0.0, 0.0, 1.0])
            },
            1e-8
        ));
        assert!(!coplanar(
            &FaceRow {
                index: 0,
                name: None,
                point: Some([0.0, 0.0, 0.0]),
                normal: Some([0.0, 0.0, 1.0])
            },
            &FaceRow {
                index: 1,
                name: None,
                point: Some([0.0, 0.0, 1e-5]),
                normal: Some([0.0, 0.0, -1.0])
            },
            1e-8
        ));
    }
}
