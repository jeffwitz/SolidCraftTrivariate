//! Adapt the audited FieldCraft sketch envelope to world-space bore treatments.
use crate::{
    bracket::BracketBlock,
    bracket_edge_edit::{edit_compact_bracket, BoreEdge, EdgeTreatment, SpanMap},
    bracket_json::{decode_blocks, encode_blocks},
    interchange,
    json::{object, parse, Json},
};

fn section_features(source: &Json) -> Result<([BoreEdge; 2], i32), &'static str> {
    let transform = source.field("transform")?;
    if transform.field("matrix")? != &parse("[[1,0,0],[0,1,0],[0,0,1]]")?
        || transform.field("translation")? != &parse("[0,0,0]")?
    {
        return Err("expected an axis aligned FieldCraft source at the stock origin");
    }
    let denominator = transform.field("denominator")?.number::<i32>()?;
    if denominator <= 0 {
        return Err("nonpositive model scale");
    }
    // Work in the integer sketch units. A quarter of a physical mm is
    // denominator quarter-units here, including non-quarter-mm hole radii.
    let quarters = |value: &Json| {
        value
            .number::<i32>()?
            .checked_mul(4)
            .ok_or("bore coordinate overflow")
    };
    let features = source.field("source_model")?.field("bores")?.items()?;
    if features.len() != 2 {
        return Err("expected one bore per wing");
    }
    let bore = |name, axis, chart_prefix| {
        let feature = features
            .iter()
            .find(|f| f.field("id").and_then(Json::text) == Ok(name))
            .ok_or("missing wing bore")?;
        if feature.field("axis")?.text()? != ["x", "y", "z"][axis] {
            return Err("unsupported bore direction");
        }
        let [x, y, z] = feature.field("center_mm")?.items()? else {
            return Err("invalid bore center");
        };
        Ok(BoreEdge {
            chart_prefix,
            axis,
            center_quarters: [quarters(x)?, quarters(y)?, quarters(z)?],
            radius_quarters: quarters(feature.field("radius_mm")?)?,
        })
    };
    Ok((
        [
            bore("east", 1, "wing_east_")?,
            bore("north", 0, "wing_north_")?,
        ],
        denominator,
    ))
}

struct WorldEdit {
    logical_blocks: Vec<String>,
    elements: Vec<SpanMap>,
    lower: f64,
    has_singular_contact: bool,
}

fn edit(source: &Json, state: &str) -> Result<WorldEdit, &'static str> {
    if source.field("schema")?.text()? != "affine_image_v1"
        || interchange::audit(source)?.field("status")?.text()? != "valid"
    {
        return Err("expected an audited FieldCraft bracket");
    }
    let (bores, scale) = section_features(source)?;
    let treatment = match state {
        "sharp" => EdgeTreatment::Sharp,
        "chamfer" => EdgeTreatment::Chamfer {
            width_quarters: scale,
        },
        "fillet" => EdgeTreatment::Fillet {
            radius_quarters: scale,
        },
        _ => return Err("expected sharp, chamfer or fillet"),
    };
    let parent = decode_blocks(source)?;
    let sketch = decode_blocks(source.field("source_model")?)?;
    let candidate = edit_compact_bracket(&sketch, bores, treatment)?;
    let mut spans = candidate.elements;
    for span in &mut spans {
        for h in &mut span.controls {
            h[3] = h[3]
                .checked_mul(i64::from(scale))
                .ok_or("world weight overflow")?;
        }
        for z in &mut span.station_quarters {
            if *z % scale != 0 {
                return Err("unrepresentable physical axial station");
            }
            *z /= scale;
        }
    }
    let lower = if candidate.singular_boundary_qualified {
        0.0
    } else {
        (candidate.determinant_lower_bound_mm3 / f64::from(scale).powi(3)).next_down()
    };
    let candidate = WorldEdit {
        logical_blocks: candidate.logical_blocks,
        elements: spans,
        lower,
        has_singular_contact: candidate.singular_boundary_qualified,
    };
    crate::fixed_span_topology::verify(
        &parent,
        &candidate.elements.iter().map(cell).collect::<Vec<_>>(),
    )?;
    Ok(candidate)
}

fn cell(e: &crate::bracket_edge_edit::SpanMap) -> BracketBlock {
    BracketBlock {
        id: format!("{}#{}", e.logical_id, e.span_index),
        degrees: e.degrees,
        controls: e.controls.clone(),
    }
}

pub fn construct(source: &Json, state: &str) -> Result<Json, &'static str> {
    let candidate = edit(source, state)?;
    let mut volumes = Vec::new();
    for id in candidate.logical_blocks {
        let elements: Vec<_> = candidate
            .elements
            .iter()
            .filter(|e| e.logical_id == id)
            .collect();
        let breaks = elements
            .iter()
            .map(|e| {
                Json::Array(
                    e.parameter_interval_32nds
                        .into_iter()
                        .map(|n| Json::Number(n.to_string()))
                        .collect(),
                )
            })
            .collect();
        let direction = if id.starts_with("wing_") { "w" } else { "v" };
        let spans = elements.into_iter().map(cell).collect();
        volumes.push(object([
            ("id", Json::String(id)),
            ("span_axis", Json::String(direction.into())),
            ("parameter_intervals_32nds", Json::Array(breaks)),
            ("spans", encode_blocks(spans)),
        ]));
    }
    Ok(object([
        (
            "schema",
            Json::String("fieldcraft_bore_rim_spans_v1".into()),
        ),
        ("state", Json::String(state.into())),
        ("volumes", Json::Array(volumes)),
        (
            "determinant_lower_bound_mm3",
            Json::Number(candidate.lower.to_string()),
        ),
        (
            "singular_cap_contact",
            Json::Bool(candidate.has_singular_contact),
        ),
        ("parent_topology", source.field("topology")?.clone()),
        ("parent_design_faces", source.field("design_faces")?.clone()),
    ]))
}

pub fn export_step(source: &Json, state: &str) -> Result<String, &'static str> {
    let candidate = edit(source, state)?;
    let cells = candidate.elements.iter().map(cell).collect();
    crate::step_boundary::export_span_boundary(cells, source, true)
}
