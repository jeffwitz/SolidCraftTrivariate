//! Certified persistence of solids built from explicit planar sketch plans.
use crate::{
    bracket_audit::BracketCertificate,
    bracket_json::encode,
    g2_cross_plate,
    incidence::enumerate_faces,
    json::{object, Json},
};
use std::collections::BTreeSet;

pub(crate) const FAMILY: &str = "g2_cross_plate_v1";
pub(crate) const PRIMITIVE: &str = "primitive_sweep_v1";
pub(crate) const PRIMITIVE_METHOD: &str = "native_primitive_sweep_template";
pub(crate) const METHOD: &str = "periodic_cubic_transverse_bore_partition";

fn candidate(recipe: &Json) -> Result<(Json, f64), &'static str> {
    let family = recipe.field("family")?.text()?;
    let (plan, (blocks, mut faces)) = match family {
        FAMILY => {
            let (parameters, plan) = crate::sketch_partition_plan::plan(recipe)?;
            (plan, g2_cross_plate::construct_partition(&parameters)?)
        }
        PRIMITIVE => (
            crate::primitive_template_json::plan(recipe)?,
            crate::primitive_template_json::construct(recipe)?,
        ),
        _ => return Err("unknown swept sketch family"),
    };
    let mut lower = f64::INFINITY;
    for block in &blocks {
        let cert = crate::spatial_certificate::certify_spatial(
            block.degrees,
            &block.controls,
            Default::default(),
        );
        let bound = cert
            .determinant_lower_bound_mm3
            .ok_or("uncertified transverse plate block")?;
        lower = lower.min(bound);
    }
    if !lower.is_finite() || lower <= 0.0 {
        return Err("unresolved periodic plate regularity");
    }
    let incidence = enumerate_faces(&blocks)?;
    faces.sort_by(|a, b| a.0.cmp(&b.0));
    let labels: BTreeSet<_> = faces.iter().map(|(side, _)| side).collect();
    if labels.len() != faces.len() || labels != incidence.exterior.iter().collect() {
        return Err("incomplete periodic plate boundary provenance");
    }
    let cert = BracketCertificate {
        incidence,
        design_faces: faces,
        determinant_lower_bound_mm3: lower,
        radial_error_upper_mm: 0.0,
    };
    let mut model = encode(family, blocks, &cert);
    if let Json::Object(fields) = &mut model {
        fields.insert("source_recipe".into(), recipe.clone());
        fields.insert("sketch_plan".into(), plan);
    }
    Ok((model, lower))
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    candidate(recipe).map(|(model, _)| model)
}

pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    let (expected, lower) = candidate(model.field("source_recipe")?)?;
    if *model != expected {
        return Err("periodic_plate_stored_partition_mismatch");
    }
    let count =
        |value: &Json| Ok::<_, &'static str>(Json::Number(value.items()?.len().to_string()));
    let topology = model.field("topology")?;
    let source = model.field("source_recipe")?;
    let (method, bores) = if source.field("family")?.text()? == FAMILY {
        (METHOD, 2)
    } else {
        (
            PRIMITIVE_METHOD,
            usize::from(source.field("operation")?.field("kind")?.text()? == "revolve"),
        )
    };
    Ok(object([
        ("status", Json::String("valid".into())),
        ("method", Json::String(method.into())),
        ("blocks", count(model.field("blocks")?)?),
        ("shared_faces", count(topology.field("interfaces")?)?),
        ("exterior_sides", count(topology.field("exterior")?)?),
        (
            "determinant_lower_bound_mm3",
            Json::Number(lower.to_string()),
        ),
        ("through_bores", Json::Number(bores.to_string())),
        ("cavities", Json::Number("0".into())),
        ("hausdorff_error_upper_bound_mm", Json::Number("0".into())),
    ]))
}
