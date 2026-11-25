//! Import supported SolidCraft XZ/YZ sketch extrusions through an exact G19 image.
use crate::{
    interchange,
    json::{object, Json},
    solidcraft_bridge,
};

fn as_object_mut(
    value: &mut Json,
) -> Result<&mut std::collections::BTreeMap<String, Json>, &'static str> {
    match value {
        Json::Object(fields) => Ok(fields),
        _ => Err("invalid SolidCraft history object"),
    }
}

fn set_sketch_plane(script: &mut Json) -> Result<(), &'static str> {
    let root = as_object_mut(script)?;
    let commands = root
        .get_mut("commands")
        .ok_or("missing SolidCraft timeline")?;
    let Json::Array(commands) = commands else {
        return Err("invalid SolidCraft timeline");
    };
    let first = commands.first_mut().ok_or("empty SolidCraft timeline")?;
    let first = as_object_mut(first)?;
    let params = first
        .get_mut("params")
        .ok_or("missing SolidCraft feature")?;
    let params = as_object_mut(params)?;
    let feature = params
        .get_mut("feature")
        .ok_or("missing SolidCraft feature")?;
    let feature = as_object_mut(feature)?;
    let plane = feature
        .get_mut("plane")
        .ok_or("missing SolidCraft sketch plane")?;
    let plane = as_object_mut(plane)?;
    let name = plane
        .get_mut("name")
        .ok_or("missing SolidCraft plane name")?;
    *name = Json::String("XY".into());
    let sketch = as_object_mut(
        feature
            .get_mut("sketch")
            .ok_or("missing SolidCraft sketch")?,
    )?;
    let counters = as_object_mut(
        sketch
            .get_mut("counters")
            .ok_or("missing SolidCraft sketch counters")?,
    )?;
    // Older histories omit the zero constraint counter. SC5's strict decoder
    // accepts the equivalent canonical spelling with k:0.
    counters
        .entry("k".into())
        .or_insert_with(|| Json::Number("0".into()));
    Ok(())
}

fn source_recipe(model: &Json) -> Result<Json, &'static str> {
    let family = model.field("schema")?.text()?;
    let fields: &[&str] = match family {
        "polygon_prism_v1" => &["vertices", "height_mm", "top_shift_mm"],
        "polygonal_perforated_prism_v1" => &["outer", "holes", "height_mm", "top_shift_mm"],
        "polygon_circular_bores_v1" => &["outer", "bores", "height_mm", "top_shift_mm"],
        _ => return Err("unsupported SolidCraft oriented-plane source family"),
    };
    let mut values = vec![("family", Json::String(family.into()))];
    for &key in fields {
        values.push((key, model.field(key)?.clone()));
    }
    Ok(Json::Object(
        values
            .into_iter()
            .map(|(key, value)| (key.into(), value))
            .collect(),
    ))
}

fn transform(plane: &str) -> Result<Json, &'static str> {
    let matrix: [[i32; 3]; 3] = match plane {
        // SolidCraft's origin XZ basis is (+X,-Z), normal +Y.
        "XZ" => [[1, 0, 0], [0, 0, 1], [0, -1, 0]],
        // Its origin YZ basis is (+Y,+Z), normal +X.
        "YZ" => [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
        _ => return Err("unsupported SolidCraft oriented plane"),
    };
    Ok(object([
        (
            "matrix",
            Json::Array(
                matrix
                    .iter()
                    .map(|row| {
                        Json::Array(row.iter().map(|x| Json::Number(x.to_string())).collect())
                    })
                    .collect(),
            ),
        ),
        (
            "translation",
            Json::Array(vec![Json::Number("0".into()); 3]),
        ),
        ("denominator", Json::Number("1".into())),
    ]))
}

/// Import a strict SC5/SC6 history on origin XZ or YZ and apply its verified
/// right-handed G19 world transform. The input history is not modified.
pub fn import_oriented_plane_timeline(script: &Json) -> Result<Json, &'static str> {
    script.keys(&["about", "commands"])?;
    let commands = script.field("commands")?.items()?;
    if commands.len() != 2 {
        return Err("SC7 accepts exactly the SC5/SC6 sketch and extrusion pair");
    }
    let sketch_feature = commands[0].field("params")?.field("feature")?;
    if sketch_feature.field("type")?.text()? != "sketch" {
        return Err("SC7 expected a sketch as the first timeline feature");
    }
    let plane = sketch_feature.field("plane")?;
    if plane.field("type")?.text()? != "origin" {
        return Err("SC7 supports origin planes only");
    }
    plane.keys(&["name", "type"])?;
    let plane_name = plane.field("name")?.text()?;
    let affine = transform(plane_name)?;
    let mut local_history = script.clone();
    set_sketch_plane(&mut local_history)?;
    let local_model = solidcraft_bridge::import_script(&local_history)?;
    let recipe = object([
        ("family", Json::String("affine_image_v1".into())),
        ("source", source_recipe(&local_model)?),
        ("transform", affine),
    ]);
    let candidate = interchange::build(&recipe)?;
    let report = interchange::audit(&candidate)?;
    if report.field("status")?.text()? != "valid" {
        return Err("SC7 affine image failed source-aware audit");
    }
    Ok(candidate)
}
