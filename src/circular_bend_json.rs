//! Persisted circular bends and their independently reconstructed audit.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    circular_bend::{build_circular_bend, BendRecipe, BendSection},
    circular_bend_audit::audit_circular_bend,
    json::{object, Json},
};

pub(crate) const FAMILY: &str = "circular_bend_v1";
fn n(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn parse(x: &Json) -> Result<BendRecipe, &'static str> {
    Ok(BendRecipe {
        section: BendSection::parse(x.field("section")?.text()?)?,
        major_radius_mm: x.field("major_radius_mm")?.number()?,
        quadrants: x.field("quadrants")?.number()?,
    })
}
pub(crate) fn build(input: &Json) -> Result<Json, &'static str> {
    input.keys(&["family", "section", "major_radius_mm", "quadrants"])?;
    let recipe = parse(input)?;
    let blocks = build_circular_bend(recipe)?;
    let report = audit_circular_bend(recipe, &blocks)?;
    let mut model = encode(FAMILY, blocks, &report.certificate);
    if let Json::Object(fields) = &mut model {
        for key in ["section", "major_radius_mm", "quadrants"] {
            fields.insert(key.into(), input.field(key)?.clone());
        }
    }
    Ok(model)
}
pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    model.keys(&[
        "schema",
        "section",
        "major_radius_mm",
        "quadrants",
        "blocks",
        "design_faces",
        "topology",
    ])?;
    let recipe = parse(model)?;
    let blocks = decode_blocks(model)?;
    let c = audit_circular_bend(recipe, &blocks)?;
    if model.field("topology")? != &topology(&c.certificate.incidence) {
        return Err("circular_bend_topology_mismatch");
    }
    if model.field("design_faces")? != &design(&c.certificate) {
        return Err("circular_bend_design_face_mismatch");
    }
    Ok(object([
        ("status", Json::String("valid".into())),
        ("method", Json::String("certified_circular_bend".into())),
        ("blocks", n(blocks.len())),
        ("shared_faces", n(c.certificate.incidence.interfaces.len())),
        ("exterior_sides", n(c.certificate.incidence.exterior.len())),
        ("periodic_angle", Json::Bool(recipe.quadrants == 4)),
        ("cavities", n(c.cavities)),
        (
            "section_fidelity_upper_mm",
            n(c.certificate.radial_error_upper_mm),
        ),
        (
            "determinant_lower_bound_mm3",
            n(c.certificate.determinant_lower_bound_mm3),
        ),
        ("ideal_volume_mm3", n(c.ideal_volume_mm3)),
    ]))
}
