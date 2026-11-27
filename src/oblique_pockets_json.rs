//! Exchange and closed-domain certification of the bounded oblique pocket experiment.
use crate::{
    bracket_audit::BracketCertificate,
    bracket_json::{design, encode_blocks, topology},
    incidence::enumerate_faces,
    json::{object, parse, Json},
    oblique_pockets::build_oblique_blind_pockets,
    spatial_certificate::{certify_spatial, CertificateBudget},
};

pub(crate) const FAMILY: &str = "oblique_blind_pockets_v1";
pub(crate) const METHOD: &str = "certified_oblique_blind_pockets";
const FIXTURE: &str = include_str!("../tests/references/oblique_blind_pockets.recipe.json");

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    if recipe != &parse(FIXTURE)? {
        return Err("oblique pocket profiles or sweep parameters outside the qualified fixture");
    }
    let blocks = build_oblique_blind_pockets()?;
    let incidence = enumerate_faces(&blocks)?;
    let mut lower = f64::INFINITY;
    for block in &blocks {
        let bound = certify_spatial(block.degrees, &block.controls, CertificateBudget::default())
            .determinant_lower_bound_mm3
            .ok_or("oblique pocket chart lacks a strict Jacobian certificate")?;
        lower = lower.min(bound);
    }
    let design_faces = incidence
        .exterior
        .iter()
        .map(|side| {
            (
                side.clone(),
                format!(
                    "{}:{}{}",
                    side.block,
                    ["u", "v", "w"][side.axis],
                    u8::from(side.upper)
                ),
            )
        })
        .collect();
    let cert = BracketCertificate {
        incidence,
        determinant_lower_bound_mm3: lower,
        radial_error_upper_mm: 0.,
        design_faces,
    };
    Ok(object([
        ("schema", Json::String(FAMILY.into())),
        ("recipe", recipe.clone()),
        ("blocks", encode_blocks(blocks)),
        ("topology", topology(&cert.incidence)),
        ("design_faces", design(&cert)),
        (
            "determinant_lower_bound_mm3",
            Json::Number(lower.to_string()),
        ),
    ]))
}

pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    if build(model.field("recipe")?)? != *model {
        return Err("oblique pocket model differs from its certified construction");
    }
    let count =
        |j: &Json| -> Result<Json, &'static str> { Ok(Json::Number(j.items()?.len().to_string())) };
    let topology = model.field("topology")?;
    Ok(object([
        ("method", Json::String(METHOD.into())),
        ("status", Json::String("valid".into())),
        ("blocks", count(model.field("blocks")?)?),
        ("shared_faces", count(topology.field("interfaces")?)?),
        ("exterior_sides", count(topology.field("exterior")?)?),
        (
            "determinant_lower_bound_mm3",
            model.field("determinant_lower_bound_mm3")?.clone(),
        ),
    ]))
}
