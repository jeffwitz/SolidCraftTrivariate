//! Produce a native SolidCraft construction sketch from the exact planar partition.
use crate::constraint_system::Ratio;
use crate::json::{object, Json};
use crate::sketch_layers::{partition, SketchLayers, SketchLoop};
use std::collections::{BTreeMap, BTreeSet};

fn num(n: impl ToString) -> Json {
    Json::Number(n.to_string())
}
fn string(s: impl ToString) -> Json {
    Json::String(s.to_string())
}
type Point = [Ratio; 2];

fn internal_edges(job: &SketchLayers) -> Result<(Vec<[Point; 2]>, usize), &'static str> {
    let result = partition(job)?;
    let mut incidence = BTreeMap::new();
    for block in &result.blocks {
        if block.degrees != [1, 1, 1] {
            return Err("unsupported curved sketch preview cell");
        }
        let mut corners = Vec::new();
        for index in [0, 4, 6, 2] {
            let [x, y, _, w] = block.controls[index].map(i128::from);
            corners.push([Ratio::new(x, w)?, Ratio::new(y, w)?]);
        }
        for i in 0..4 {
            let mut edge = [corners[i], corners[(i + 1) % 4]];
            edge.sort();
            *incidence.entry(edge).or_insert(0) += 1;
        }
    }
    if incidence.values().any(|&n| n > 2) {
        return Err("nonmanifold sketch preview edge");
    }
    Ok((
        incidence
            .into_iter()
            .filter_map(|(e, n)| (n == 2).then_some(e))
            .collect(),
        result.blocks.len(),
    ))
}

pub fn snapshot_script(script: &Json) -> Result<&Json, &'static str> {
    Ok(
        if script.field("script").is_ok() || script.field("schema").is_ok() {
            if script
                .field("schema")
                .ok()
                .map(Json::text)
                .transpose()?
                .is_some_and(|s| s != "solidcraft_trivariate_timeline_v1")
                || !script.field("omitted")?.items()?.is_empty()
            {
                return Err("unsupported or incomplete SolidCraft snapshot");
            }
            script.field("script")?
        } else {
            script
        },
    )
}

pub(crate) fn source_feature(script: &Json, source_id: u32) -> Result<(&Json, u32), &'static str> {
    let mut ids = BTreeSet::new();
    let mut source = None;
    for cmd in script.field("commands")?.items()? {
        if cmd.field("command")?.text()? != "timeline.add" {
            continue;
        }
        let feature = cmd.field("params")?.field("feature")?;
        let id = feature.field("id")?.number::<u32>()?;
        if !ids.insert(id) {
            return Err("duplicate source feature id");
        }
        if id == source_id {
            source = Some(feature);
        }
    }
    let source = source.ok_or("missing source sketch")?;
    if source.field("type")?.text()? != "sketch" {
        return Err("expected a sketch feature");
    }
    let target_id = ids
        .last()
        .and_then(|id| id.checked_add(1))
        .ok_or("feature id overflow")?;
    Ok((source, target_id))
}

pub fn descriptor(input: &Json, source_id: u32) -> Result<Json, &'static str> {
    let script = snapshot_script(input)?;
    let (source, target_id) = source_feature(script, source_id)?;
    let mut loops = crate::solidcraft_polygon::preview_contours(source)?
        .into_iter()
        .enumerate()
        .map(|(i, vertices)| SketchLoop {
            id: format!("loop{i}"),
            vertices,
        });
    let job = SketchLayers {
        outer: loops.next().ok_or("empty profile")?,
        cuts: loops.map(|l| (l, [0, 1])).collect(),
        guides: Vec::new(),
        height: 1,
        force_triangulation: false,
        operations: None,
        shift: [0, 0],
        strategy: Default::default(),
    };
    let (edges, blocks) = internal_edges(&job)?;
    let mut points = vec![object([
        ("id", string("origin")),
        ("pos", Json::Array(vec![num("0.0"), num("0.0")])),
        ("fixed", Json::Bool(true)),
    ])];
    let mut curves = Vec::new();
    let mut exact = Vec::new();
    for (i, edge) in edges.iter().enumerate() {
        let id = format!("l{}", i + 1);
        let start = points.len();
        for (p, suffix) in edge.iter().zip(["start", "end"]) {
            let pos = p
                .iter()
                .map(|q| {
                    Json::Number(format!(
                        "{:?}",
                        q.numerator() as f64 / q.denominator() as f64
                    ))
                })
                .collect();
            points.push(object([
                ("id", string(format!("{id}.{suffix}"))),
                ("pos", Json::Array(pos)),
                ("fixed", Json::Bool(true)),
            ]));
        }
        curves.push(object([
            ("id", string(id)),
            ("type", string("line")),
            ("a", num(start)),
            ("b", num(start + 1)),
            ("construction", Json::Bool(true)),
        ]));
        exact.push(Json::Array(
            edge.iter()
                .map(|p| {
                    Json::Array(
                        p.iter()
                            .map(|q| Json::Array(vec![num(q.numerator()), num(q.denominator())]))
                            .collect(),
                    )
                })
                .collect(),
        ));
    }
    let sketch = object([
        ("points", Json::Array(points)),
        (
            "counters",
            object([("k", num(0)), ("l", num(curves.len()))]),
        ),
        ("constraints", Json::Array(vec![])),
        ("curves", Json::Array(curves)),
    ]);
    let feature = object([
        ("type", string("sketch")),
        ("name", string(format!("SCT Partition {source_id}"))),
        ("plane", source.field("plane")?.clone()),
        ("id", num(target_id)),
        ("suppressed", Json::Bool(false)),
        ("sketch", sketch),
    ]);
    Ok(object([
        ("schema", string("solidcraft_sketch_preview_v1")),
        ("source_feature_id", num(source_id)),
        ("source_feature", source.clone()),
        ("source_script", script.clone()),
        ("feature", feature),
        ("internal_edges", Json::Array(exact)),
        ("block_count", num(blocks)),
    ]))
}
