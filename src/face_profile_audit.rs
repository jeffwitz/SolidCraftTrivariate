//! Reproduce face-profile charts directly from their persisted local geometry.
use crate::{
    face_profile_extrusion::{self, FaceProfileExtrusion},
    json::{Json, object},
};

fn coords<const N: usize>(j: &Json) -> Result<[f64; N], &'static str> {
    j.items()?
        .iter()
        .map(Json::number)
        .collect::<Result<Vec<f64>, _>>()?
        .try_into()
        .map_err(|_| "wrong coordinate dimension")
}

pub(crate) fn build_recipe(j: &Json) -> Result<Json, &'static str> {
    j.keys(&[
        "family",
        "operation",
        "bounds_mm",
        "thickness_mm",
        "depth_mm",
        "tolerance_mm",
        "frame",
        "spans",
    ])?;
    if j.field("family")?.text()? != "face_profile_extrusion_v1" {
        return Err("wrong face-profile recipe family");
    }
    let bounds = j
        .field("bounds_mm")?
        .items()?
        .iter()
        .map(coords)
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| "expected two face intervals")?;
    let frame = j
        .field("frame")?
        .items()?
        .iter()
        .map(coords)
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| "expected origin and three frame axes")?;
    let spans = j
        .field("spans")?
        .items()?
        .iter()
        .map(|s| s.items()?.iter().map(coords).collect())
        .collect::<Result<Vec<_>, _>>()?;
    face_profile_extrusion::build(&FaceProfileExtrusion {
        bounds,
        frame,
        spans,
        thickness: j.field("thickness_mm")?.number()?,
        depth: j.field("depth_mm")?.number()?,
        operation: j.field("operation")?.text()?.to_owned(),
        tolerance: j.field("tolerance_mm")?.number()?,
    })
}

pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    let rebuilt = build_recipe(model.field("source_recipe")?)?;
    if rebuilt != *model {
        return Err("face profile metadata or geometry differs from reconstructed model");
    }
    Ok(object([
        ("status", Json::String("valid".into())),
        (
            "method",
            Json::String("certified_face_profile_extrusion".into()),
        ),
        (
            "determinant_lower_bound_mm3",
            model.field("determinant_lower_bound_mm3")?.clone(),
        ),
        (
            "blocks",
            Json::Number(model.field("blocks")?.items()?.len().to_string()),
        ),
    ]))
}
