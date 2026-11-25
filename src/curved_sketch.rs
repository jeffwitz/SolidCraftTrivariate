//! Exact polynomial curves for rectangular sketch cuts on guided partitions.
//!
//! The scaffold constructor determines the topology. This adapter accepts only
//! whole rectangular feature sides, so every curved trace is represented once
//! and copied identically to all incident blocks. Through features use a 2/3
//! corner contraction and half-span midpoint handles (C1 joins). Blind cuts
//! use 4/5 contraction and quarter-span handles (G0 joins).

use crate::bracket::BracketBlock;
use crate::bracket_json::{decode_blocks, encode_blocks};
use crate::json::{object, Json};

const FAMILY: &str = "curved_sketch_layers_v1";
// Represents thirds, fifths, eighths, and the exact Coons-center division by 4.
const SCALE: i64 = 480;

fn n(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}
fn text(value: &str) -> Json {
    Json::String(value.to_owned())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Point {
    x: i64,
    y: i64,
    z: i64,
}

#[derive(Clone)]
struct Feature {
    id: String,
    vertices: Vec<[i64; 2]>,
    bounds: [i64; 4], // xmin, xmax, ymin, ymax
    shrink: (i64, i64),
    handle: (i64, i64),
    through: bool,
}

#[derive(Clone, Copy)]
struct EdgeCurve {
    a: [i64; 2],
    b: [i64; 2],
    outward: i64,
    transverse_span: i64,
    sag: (i64, i64),
}

fn xy(value: &Json) -> Result<[i64; 2], &'static str> {
    let [x, y] = value.items()? else {
        return Err("curved_sketch_expected_xy_point");
    };
    Ok([x.number()?, y.number()?])
}

fn bounds(vertices: &[[i64; 2]]) -> [i64; 4] {
    [
        vertices.iter().map(|p| p[0]).min().unwrap(),
        vertices.iter().map(|p| p[0]).max().unwrap(),
        vertices.iter().map(|p| p[1]).min().unwrap(),
        vertices.iter().map(|p| p[1]).max().unwrap(),
    ]
}

fn on_segment(p: [i64; 2], a: [i64; 2], b: [i64; 2]) -> bool {
    if a[0] == b[0] {
        p[0] == a[0] && a[1].min(b[1]) <= p[1] && p[1] <= a[1].max(b[1])
    } else if a[1] == b[1] {
        p[1] == a[1] && a[0].min(b[0]) <= p[0] && p[0] <= a[0].max(b[0])
    } else {
        false
    }
}

fn features(scaffold: &Json) -> Result<(Vec<Feature>, Vec<EdgeCurve>), &'static str> {
    let height = scaffold.field("height_mm")?.number::<i64>()?;
    let mut features = Vec::new();
    let mut edges = Vec::new();
    let outer = scaffold.field("outer")?;
    let outer_vertices = outer
        .field("vertices")?
        .items()?
        .iter()
        .map(xy)
        .collect::<Result<Vec<_>, _>>()?;
    if outer_vertices.len() < 4 {
        return Err("curved_sketch_bad_outer_loop");
    }
    let ob = bounds(&outer_vertices);
    for i in 0..outer_vertices.len() {
        let a = outer_vertices[i];
        let b = outer_vertices[(i + 1) % outer_vertices.len()];
        if a[1] == ob[2] && b[1] == ob[2] {
            edges.push(EdgeCurve {
                a,
                b,
                outward: -1,
                transverse_span: ob[3] - ob[2],
                sag: (1, 8),
            });
        }
    }
    for cut in scaffold.field("cuts")?.items()? {
        let vertices = cut
            .field("vertices")?
            .items()?
            .iter()
            .map(xy)
            .collect::<Result<Vec<_>, _>>()?;
        if vertices.len() != 4 {
            return Err("curved_sketch_requires_rectangular_features");
        }
        let b = bounds(&vertices);
        if b[0] == b[1] || b[2] == b[3] {
            return Err("curved_sketch_degenerate_feature");
        }
        for i in 0..4 {
            let a = vertices[i];
            let c = vertices[(i + 1) % 4];
            if a[0] != c[0] && a[1] != c[1] {
                return Err("curved_sketch_requires_axis_aligned_features");
            }
        }
        let interval = cut.field("z_interval_mm")?.items()?;
        if interval.len() != 2 {
            return Err("curved_sketch_bad_feature_interval");
        }
        let through = interval[0].number::<i64>()? == 0 && interval[1].number::<i64>()? == height;
        let shrink = if through { (2, 3) } else { (4, 5) };
        let handle = if through { (1, 2) } else { (1, 4) };
        let id = cut.field("id")?.text()?.to_owned();
        features.push(Feature {
            id,
            vertices: vertices.clone(),
            bounds: b,
            shrink,
            handle,
            through,
        });
        for i in 0..4 {
            let a = vertices[i];
            let c = vertices[(i + 1) % 4];
            let (outward, transverse_span) = if a[1] == c[1] {
                (if 2 * a[1] > b[2] + b[3] { 1 } else { -1 }, b[3] - b[2])
            } else {
                (if 2 * a[0] > b[0] + b[1] { 1 } else { -1 }, b[1] - b[0])
            };
            let sag = if through { (1, 3) } else { (1, 5) };
            edges.push(EdgeCurve {
                a,
                b: c,
                outward,
                transverse_span,
                sag,
            });
        }
    }
    for i in 0..features.len() {
        for j in 0..i {
            let a = features[i].bounds;
            let b = features[j].bounds;
            if a[0] < b[1] && b[0] < a[1] && a[2] < b[3] && b[2] < a[3] {
                return Err("curved_sketch_intersecting_features_out_of_scope");
            }
        }
    }
    Ok((features, edges))
}

fn scale_ratio(value: i64, num: i64, den: i64) -> Result<i64, &'static str> {
    let n = value
        .checked_mul(num)
        .and_then(|x| x.checked_mul(SCALE))
        .ok_or("curved_sketch_integer_overflow")?;
    if n % den != 0 {
        return Err("curved_sketch_common_scale_not_exact");
    }
    Ok(n / den)
}

fn remap(p: [i64; 2], z: i64, features: &[Feature]) -> Result<Point, &'static str> {
    for feature in features {
        let on_boundary = (0..feature.vertices.len()).any(|i| {
            on_segment(
                p,
                feature.vertices[i],
                feature.vertices[(i + 1) % feature.vertices.len()],
            )
        });
        if !on_boundary {
            continue;
        }
        let [xmin, xmax, ymin, ymax] = feature.bounds;
        // p' = center + shrink * (p - center), kept as an exact integer numerator.
        let tx = (xmin + xmax) * feature.shrink.1 + (2 * p[0] - xmin - xmax) * feature.shrink.0;
        let ty = (ymin + ymax) * feature.shrink.1 + (2 * p[1] - ymin - ymax) * feature.shrink.0;
        let den = 2 * feature.shrink.1;
        let x = scale_ratio(tx, 1, den)?;
        let y = scale_ratio(ty, 1, den)?;
        return Ok(Point { x, y, z: z * SCALE });
    }
    Ok(Point {
        x: p[0] * SCALE,
        y: p[1] * SCALE,
        z: z * SCALE,
    })
}

fn matching_curve(
    a: [i64; 2],
    b: [i64; 2],
    edges: &[EdgeCurve],
) -> Result<Option<usize>, &'static str> {
    for (i, edge) in edges.iter().enumerate() {
        if on_segment(a, edge.a, edge.b) && on_segment(b, edge.a, edge.b) {
            // Curved feature sides must be whole scaffold edges. This keeps
            // their intended parameterization and corner continuity explicit.
            if !((a == edge.a && b == edge.b) || (a == edge.b && b == edge.a)) {
                return Err("curved_sketch_feature_side_subdivided");
            }
            return Ok(Some(i));
        }
    }
    Ok(None)
}

fn edge_control(
    a: Point,
    b: Point,
    source_a: [i64; 2],
    source_b: [i64; 2],
    edges: &[EdgeCurve],
    seen: &mut [bool],
) -> Result<[i64; 4], &'static str> {
    let mut x = (a.x + b.x) / 2;
    let mut y = (a.y + b.y) / 2;
    if (a.x + b.x) % 2 != 0 || (a.y + b.y) % 2 != 0 {
        return Err("curved_sketch_edge_midpoint_not_exact");
    }
    if let Some(i) = matching_curve(source_a, source_b, edges)? {
        seen[i] = true;
        let edge = &edges[i];
        let numerator = edge.transverse_span * edge.sag.0 * SCALE;
        let denominator = edge.sag.1;
        if numerator % denominator != 0 {
            return Err("curved_sketch_handle_not_exact_on_common_scale");
        }
        let offset = edge.outward * numerator / denominator;
        if edge.a[0] == edge.b[0] {
            x += offset;
        } else {
            y += offset;
        }
    }
    Ok([x, y, (a.z + b.z) / 2, SCALE])
}

fn map_controls(
    block: BracketBlock,
    features: &[Feature],
    edges: &[EdgeCurve],
    seen: &mut [bool],
    top_shift: [i64; 2],
    height: i64,
) -> Result<BracketBlock, &'static str> {
    if block.degrees != [1, 1, 1]
        || block.controls.len() != 8
        || block.controls.iter().any(|p| p[3] != 1)
    {
        return Err("curved_sketch_requires_integral_trilinear_scaffold");
    }
    let raw = block.controls;
    let corner = |i: usize, j: usize, k: usize| raw[(i * 2 + j) * 2 + k];
    let z0 = corner(0, 0, 0)[2];
    let z1 = corner(0, 0, 1)[2];
    let scaled_shift = top_shift.map(|delta| delta * z0);
    if scaled_shift.iter().any(|x| x % height != 0) {
        return Err("curved_sketch_nonintegral_layer_offset");
    }
    let offset = scaled_shift.map(|x| x / height);
    let raw_xy = [
        corner(0, 0, 0),
        corner(1, 0, 0),
        corner(0, 1, 0),
        corner(1, 1, 0),
    ]
    .map(|p| [p[0] - offset[0], p[1] - offset[1]]);
    if [corner(1, 0, 0)[2], corner(0, 1, 0)[2], corner(1, 1, 0)[2]]
        .iter()
        .any(|&z| z != z0)
        || [corner(1, 0, 1)[2], corner(0, 1, 1)[2], corner(1, 1, 1)[2]]
            .iter()
            .any(|&z| z != z1)
    {
        return Err("curved_sketch_requires_planar_axial_layers");
    }
    let shift = [
        corner(0, 0, 1)[0] - corner(0, 0, 0)[0],
        corner(0, 0, 1)[1] - corner(0, 0, 0)[1],
    ];
    for (i, j) in [(1, 0), (0, 1), (1, 1)] {
        let a = corner(i, j, 0);
        let b = corner(i, j, 1);
        if [b[0] - a[0], b[1] - a[1]] != shift {
            return Err("curved_sketch_requires_uniform_axial_shift");
        }
    }
    let p = raw_xy
        .map(|q| remap(q, z0, features))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
    let [p00, p10, p01, p11] =
        <[Point; 4]>::try_from(p).map_err(|_| "curved_sketch_bad_corner_count")?;
    let mids = [
        edge_control(p00, p10, raw_xy[0], raw_xy[1], edges, seen)?,
        edge_control(p01, p11, raw_xy[2], raw_xy[3], edges, seen)?,
        edge_control(p00, p01, raw_xy[0], raw_xy[2], edges, seen)?,
        edge_control(p10, p11, raw_xy[1], raw_xy[3], edges, seen)?,
    ];
    let corners = [p00, p10, p01, p11].map(|p| [p.x, p.y, p.z, SCALE]);
    let sum = |axis: usize| mids.iter().map(|p| p[axis]).sum::<i64>();
    let mut center = [0_i64; 4];
    for axis in 0..3 {
        let numerator = 2 * sum(axis) - corners.iter().map(|p| p[axis]).sum::<i64>();
        if numerator % 4 != 0 {
            return Err("curved_sketch_coons_center_not_exact");
        }
        center[axis] = numerator / 4;
    }
    center[3] = SCALE;
    let mut controls = Vec::with_capacity(18);
    for i in 0..3 {
        for j in 0..3 {
            let q = match (i, j) {
                (0, 0) => corners[0],
                (2, 0) => corners[1],
                (0, 2) => corners[2],
                (2, 2) => corners[3],
                (1, 0) => mids[0],
                (1, 2) => mids[1],
                (0, 1) => mids[2],
                (2, 1) => mids[3],
                (1, 1) => center,
                _ => unreachable!(),
            };
            let q = [
                q[0] + offset[0] * SCALE,
                q[1] + offset[1] * SCALE,
                q[2],
                SCALE,
            ];
            controls.push(q);
            controls.push([
                q[0] + shift[0] * SCALE,
                q[1] + shift[1] * SCALE,
                z1 * SCALE,
                SCALE,
            ]);
        }
    }
    Ok(BracketBlock {
        id: block.id,
        degrees: [2, 2, 1],
        controls,
    })
}

/// Construct a rationally exact curved model from a sketch-layer scaffold.
/// Input: `{ "family":"curved_sketch_layers_v1", "scaffold": <recipe> }`.
pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&["family", "scaffold"])?;
    if recipe.field("family")?.text()? != FAMILY {
        return Err("curved_sketch_unknown_family");
    }
    let scaffold = recipe.field("scaffold")?;
    if scaffold.field("family")?.text()? != "sketch_layer_partition_v1" {
        return Err("curved_sketch_requires_guided_layer_scaffold");
    }
    let base = crate::interchange::build(scaffold)?;
    let (features, edges) = features(scaffold)?;
    let height = scaffold.field("height_mm")?.number::<i64>()?;
    let shift = scaffold
        .field("top_shift_mm")
        .ok()
        .map(xy)
        .transpose()?
        .unwrap_or([0, 0]);
    let mut seen = vec![false; edges.len()];
    let blocks = decode_blocks(&base)?
        .into_iter()
        .map(|block| map_controls(block, &features, &edges, &mut seen, shift, height))
        .collect::<Result<Vec<_>, _>>()?;
    if seen.contains(&false) {
        return Err("curved_sketch_feature_side_missing_from_scaffold");
    }
    let topology = base.field("topology")?.clone();
    let design_faces = base.field("design_faces")?.clone();
    let scaffold_plan = base.field("partition_plan").ok().cloned();
    let mut model = object([
        ("schema", text(FAMILY)),
        ("scaffold", scaffold.clone()),
        ("blocks", encode_blocks(blocks)),
        ("topology", topology),
        ("design_faces", design_faces),
        (
            "degree_convention",
            text("tensor_bernstein_xy_quadratic_z_linear"),
        ),
        (
            "curve_construction",
            object([
                (
                    "method",
                    text("contracted_rectangle_quadratic_midpoint_controls"),
                ),
                ("common_integer_scale", n(SCALE)),
                ("through_corner_scale", Json::Array(vec![n(2), n(3)])),
                (
                    "through_handle_fraction_of_reduced_span",
                    Json::Array(vec![n(1), n(2)]),
                ),
                ("blind_corner_scale", Json::Array(vec![n(4), n(5)])),
                (
                    "blind_handle_fraction_of_reduced_span",
                    Json::Array(vec![n(1), n(4)]),
                ),
                ("through_join", text("C1")),
                ("blind_join", text("G0")),
            ]),
        ),
        (
            "feature_curve_provenance",
            Json::Array(
                features
                    .iter()
                    .map(|f| {
                        object([
                            ("id", text(&f.id)),
                            ("kind", text(if f.through { "through" } else { "blind" })),
                            (
                                "source_loop",
                                Json::Array(
                                    f.vertices
                                        .iter()
                                        .map(|p| Json::Array(p.iter().copied().map(n).collect()))
                                        .collect(),
                                ),
                            ),
                            ("corner_scale_numerator", n(f.shrink.0)),
                            ("corner_scale_denominator", n(f.shrink.1)),
                            ("handle_numerator", n(f.handle.0)),
                            ("handle_denominator", n(f.handle.1)),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    if let (Some(plan), Json::Object(fields)) = (scaffold_plan, &mut model) {
        fields.insert("scaffold_partition_plan".into(), plan);
    }
    Ok(model)
}
