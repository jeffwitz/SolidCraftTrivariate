//! Selection-scoped host B-Rep chamfer and fillet dispatch.
//!
//! The host owns B-Rep topology. We resolve its sparse `model.edges` / `model.faces` indices
//! against the named body, then ask SolidCraft's native solid operation to edit those exact
//! edges. Face picks expand only to the selected planar patch's perimeter; coplanar adjacent
//! faces are treated as one patch so an artificial partition seam is not selected. No tangent
//! chain propagation is requested. This does not edit trivariate recipe geometry.

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

#[derive(Clone, Debug)]
struct BodyMetrics {
    volume: f64,
    closed: bool,
    bounds: ([f64; 3], [f64; 3]),
}

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
        bounds: (min, max),
    })
}

#[derive(Clone, Copy, Debug)]
struct LocalCorner {
    origin: [f64; 3],
    x_dir: [f64; 3],
    y_dir: [f64; 3],
    edge_length: f64,
    tool_length: f64,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn unit_axis(v: [f64; 3], tolerance: f64) -> Option<usize> {
    let axis = (0..3).find(|&i| v[i].abs() > 1.0 - tolerance)?;
    ((0..3)
        .filter(|&i| i != axis)
        .all(|i| v[i].abs() <= tolerance))
    .then_some(axis)
}

/// Prove the picked edge is a full-span, convex, axis-aligned stock corner. The tool only
/// sweeps a quarter-circle corner prism over that full edge; curved or partial edges refuse.
fn prove_local_corner(
    body: &str,
    item: &Value,
    radius: f64,
    edges: &[EdgeRow],
    metrics: &BodyMetrics,
) -> Result<LocalCorner, PluginError> {
    if item.get("type").and_then(Value::as_str) != Some("edge") {
        return Err(bad(
            "selection",
            "analytic fallback requires one selected edge",
        ));
    }
    let index = item
        .get("index")
        .and_then(Value::as_u64)
        .ok_or_else(|| bad("selection", "selected edge has no index"))?;
    let edge = edges
        .iter()
        .find(|e| e.index == index)
        .ok_or_else(|| bad("selection", "selected edge vanished"))?;
    let (Some(a), Some(b), Some(name)) = (edge.start, edge.end, edge.name.as_deref()) else {
        return Err(bad(
            "selection",
            "analytic fallback requires a straight named edge with endpoints",
        ));
    };
    let (lo, hi) = metrics.bounds;
    let size = std::array::from_fn::<_, 3, _>(|i| hi[i] - lo[i]);
    let scale = size.into_iter().fold(1.0_f64, f64::max);
    let tol = scale * 2.0e-8;
    let chord_vec = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let chord = distance(a, b);
    if !chord.is_finite() || chord <= tol {
        return Err(bad(
            "selection",
            "analytic fallback requires a nonzero straight edge",
        ));
    }
    let Some(axis) = unit_axis(chord_vec.map(|v| v / chord), 1e-9) else {
        return Err(bad(
            "selection",
            "analytic fallback requires an axis-aligned straight edge",
        ));
    };
    if chord <= tol || (edge.length - chord).abs() > tol || (chord - size[axis]).abs() > tol {
        return Err(bad(
            "selection",
            "analytic fallback requires an edge spanning the full stock length",
        ));
    }
    let start = if a[axis] < b[axis] { a } else { b };
    let end = if a[axis] < b[axis] { b } else { a };
    if (start[axis] - lo[axis]).abs() > tol || (end[axis] - hi[axis]).abs() > tol {
        return Err(bad(
            "selection",
            "analytic fallback edge must span the complete stock bound",
        ));
    }
    let transverse: Vec<usize> = (0..3).filter(|&i| i != axis).collect();
    for &i in &transverse {
        if (start[i] - end[i]).abs() > tol
            || ((start[i] - lo[i]).abs() > tol && (start[i] - hi[i]).abs() > tol)
        {
            return Err(bad(
                "selection",
                "analytic fallback requires an outer stock corner",
            ));
        }
    }
    if radius >= size[transverse[0]].min(size[transverse[1]]) - tol {
        return Err(bad(
            "radius_mm",
            "radius must fit within both stock dimensions",
        ));
    }
    let face_rows = faces(body)?;
    let sides = edge_face_sides(name);
    if sides.len() != 2 || sides[0] == sides[1] {
        return Err(bad(
            "selection",
            "analytic fallback requires two distinct planar exterior faces",
        ));
    }
    let face_for = |side: &str| face_rows.iter().find(|f| f.name.as_deref() == Some(side));
    let f0 = face_for(sides[0])
        .ok_or_else(|| bad("selection", "could not resolve first adjacent face"))?;
    let f1 = face_for(sides[1])
        .ok_or_else(|| bad("selection", "could not resolve second adjacent face"))?;
    let (Some(n0), Some(p0), Some(n1), Some(p1)) = (f0.normal, f0.point, f1.normal, f1.point)
    else {
        return Err(bad(
            "selection",
            "analytic fallback requires planar adjacent faces",
        ));
    };
    let Some(ax0) = unit_axis(n0, 1e-9) else {
        return Err(bad("selection", "adjacent face is not axis-aligned"));
    };
    let Some(ax1) = unit_axis(n1, 1e-9) else {
        return Err(bad("selection", "adjacent face is not axis-aligned"));
    };
    if ax0 == axis || ax1 == axis || ax0 == ax1 || dot(n0, n1).abs() > 1e-9 {
        return Err(bad(
            "selection",
            "adjacent faces are not the two perpendicular exterior sides",
        ));
    }
    for (normal, point, side_axis) in [(n0, p0, ax0), (n1, p1, ax1)] {
        let bound = if normal[side_axis] > 0.0 {
            hi[side_axis]
        } else {
            lo[side_axis]
        };
        if (point[side_axis] - bound).abs() > tol || (start[side_axis] - bound).abs() > tol {
            return Err(bad("selection", "edge is not on the exterior stock planes"));
        }
    }
    let mut e_dir = [0.0; 3];
    e_dir[axis] = 1.0;
    let (x_dir, opposite) = if dot(cross(e_dir, n0), [-n1[0], -n1[1], -n1[2]]) > 1.0 - 1e-9 {
        (n0, n1)
    } else if dot(cross(e_dir, n1), [-n0[0], -n0[1], -n0[2]]) > 1.0 - 1e-9 {
        (n1, n0)
    } else {
        return Err(bad(
            "selection",
            "adjacent face orientation does not prove a convex outer corner",
        ));
    };
    let y_dir = cross(e_dir, x_dir);
    debug_assert!(dot(y_dir, [-opposite[0], -opposite[1], -opposite[2]]) > 1.0 - 1e-9);
    let inset = 1.0_f64.min(size[axis] / 4.0);
    let origin = [
        start[0] - e_dir[0] * inset,
        start[1] - e_dir[1] * inset,
        start[2] - e_dir[2] * inset,
    ];
    Ok(LocalCorner {
        origin,
        x_dir,
        y_dir,
        edge_length: size[axis],
        tool_length: size[axis] + 2.0 * inset,
    })
}

fn is_full_box_perimeter(edge: &EdgeRow, bounds: ([f64; 3], [f64; 3]), tol: f64) -> bool {
    let (Some(a), Some(b), Some(mid)) = (edge.start, edge.end, edge.mid) else {
        return false;
    };
    let chord = distance(a, b);
    if (edge.length - chord).abs() > tol || chord <= tol {
        return false;
    }
    let delta = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let Some(axis) = unit_axis(delta.map(|v| v / chord), 1e-9) else {
        return false;
    };
    let (lo, hi) = bounds;
    (a[axis].min(b[axis]) - lo[axis]).abs() <= tol
        && (a[axis].max(b[axis]) - hi[axis]).abs() <= tol
        && (mid[axis] - (lo[axis] + hi[axis]) * 0.5).abs() <= tol
        && (0..3).filter(|&i| i != axis).all(|i| {
            (a[i] - b[i]).abs() <= tol
                && ((a[i] - lo[i]).abs() <= tol || (a[i] - hi[i]).abs() <= tol)
        })
}

/// A nearby non-perimeter trim on either supporting plane means a native constant-radius
/// blend would have to trim another loop. Use the local analytic sweep before asking the host
/// fillet command, whose failed native feature can be retained by some host transactions.
fn has_nearby_support_trim(
    body: &str,
    selected: &EdgeRow,
    radius: f64,
    edges: &[EdgeRow],
    metrics: &BodyMetrics,
) -> Result<bool, PluginError> {
    let selected_name = selected
        .name
        .as_deref()
        .ok_or_else(|| bad("selection", "selected edge has no adjacent-face names"))?;
    let side_names = edge_face_sides(selected_name);
    let face_rows = faces(body)?;
    let scale = (0..3)
        .map(|i| metrics.bounds.1[i] - metrics.bounds.0[i])
        .fold(1.0_f64, f64::max);
    let tol = scale * 2.0e-8;
    let selected_start = selected
        .start
        .ok_or_else(|| bad("selection", "selected edge has no start point"))?;
    let selected_end = selected
        .end
        .ok_or_else(|| bad("selection", "selected edge has no end point"))?;
    for side in side_names {
        let Some(face) = face_rows.iter().find(|f| f.name.as_deref() == Some(side)) else {
            continue;
        };
        let (Some(normal), Some(point)) = (face.normal, face.point) else {
            continue;
        };
        for edge in edges {
            if edge.index == selected.index
                || !edge.name.as_deref().is_some_and(|n| edge_has_face(n, side))
                || is_full_box_perimeter(edge, metrics.bounds, tol)
            {
                continue;
            }
            let samples: Vec<[f64; 3]> = [edge.start, edge.end, edge.mid]
                .into_iter()
                .flatten()
                .collect();
            if samples.is_empty() {
                continue;
            }
            let lies_on_support = samples.iter().all(|p| {
                let d = [p[0] - point[0], p[1] - point[1], p[2] - point[2]];
                dot(d, normal).abs() <= tol
            });
            // Every point of this trim is within its total arc length of any one sample, so
            // this lower bound cannot miss an unsampled crossing between start/mid/end.
            let near_edge = samples.iter().any(|p| {
                segment_distance(*p, selected_start, selected_end) - edge.length <= radius + tol
            });
            if lies_on_support && near_edge {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn analytic_corner_cut(
    body: &str,
    radius: f64,
    corner: LocalCorner,
    before: &BodyMetrics,
) -> Result<String, PluginError> {
    // A constant-radius quarter-circle tool, extruded just past the two end planes.
    let sketch_raw = commands::run(
        "sketch.create",
        &json!({
            "name":"SCT Local Analytic Fillet",
            "plane":{"origin":corner.origin,"x_dir":corner.x_dir,"y_dir":corner.y_dir}
        })
        .to_string(),
    )?;
    let sketch: Value = serde_json::from_str(&sketch_raw).map_err(|e| bad("fillet", e))?;
    let sketch_id = sketch
        .get("sketch")
        .and_then(Value::as_u64)
        .ok_or_else(|| {
            bad(
                "fillet",
                "host did not return the temporary sketch feature id",
            )
        })?;
    let r = radius;
    let path = json!({"points":[[-r,0.0],[-r,-1.0],[1.0,-1.0],[1.0,r],[0.0,r]],"closed":false});
    let built = (|| -> Result<String, PluginError> {
        commands::run("sketch.line", &path.to_string())?;
        commands::run(
            "sketch.arc.center_point",
            &json!({"center":[-r,r],"start":[0.0,r],"end":[-r,0.0],"clockwise":true}).to_string(),
        )?;
        commands::run("sketch.finish", "{}")?;
        commands::run("solid.extrude", &json!({"sketch":sketch_id,"distance":corner.tool_length,"direction":"positive","operation":"cut","targets":[body]}).to_string())
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
    let bbox_same = (0..3).all(|i| {
        (before.bounds.0[i] - after.bounds.0[i]).abs() <= 2e-8
            && (before.bounds.1[i] - after.bounds.1[i]).abs() <= 2e-8
    });
    let removed = before.volume - after.volume;
    let upper = radius * radius * (1.0 - std::f64::consts::FRAC_PI_4) * corner.edge_length
        + before.volume * 2e-8;
    if !after.closed || !bbox_same || removed <= before.volume * 1e-12 || removed > upper {
        let _ = commands::run("timeline.delete", &json!({"feature":sketch_id}).to_string());
        return Err(bad(
            "fillet",
            "analytic local cut failed the one-solid, unchanged-bounds, or volume-envelope check",
        ));
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
    let mut args = json!({"body":body,"edges":commands_edges,"tangent_chain":false});
    if entry == "fillet" {
        args["radius"] = json!(amount);
        args["fit"] = json!("exact");
    } else {
        args["distance"] = json!(amount);
    }
    if entry == "fillet" && selected_edges.len() == 1 {
        // Do not probe with a native feature when a nearby trim loop is already known: the
        // host may retain failed feature diagnostics in the active transaction. Preserve native
        // support for regular edges and fall back only for the rigorously bounded corner case.
        if let Ok(before) = body_metrics(&body) {
            if let Ok(corner) = prove_local_corner(&body, &selection[0], amount, &edges, &before) {
                if let Some(selected) = edges.iter().find(|e| e.index == selected_edges[0]) {
                    if has_nearby_support_trim(&body, selected, amount, &edges, &before)
                        .unwrap_or(false)
                    {
                        return analytic_corner_cut(&body, amount, corner, &before);
                    }
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
