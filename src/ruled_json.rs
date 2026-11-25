//! Interchange for the globally invertible, height-dependent similarity sweep.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    json::{object, Json},
    ruled_sweep::{build_ruled_sweep, SweepRecipe, SweepSection},
    ruled_sweep_audit::audit_ruled_sweep,
};

const FAMILY: &str = "ruled_section_sweep_v1";
const PARAMETERS: [&str; 4] = [
    "section",
    "height_mm",
    "top_multiplier_eighths",
    "top_translation_eighths_mm",
];
fn number(v: impl ToString) -> Json {
    Json::Number(v.to_string())
}
fn parse(input: &Json) -> Result<SweepRecipe, &'static str> {
    let pair = |key| -> Result<[i32; 2], &'static str> {
        let [x, y] = input.field(key)?.items()? else {
            return Err("expected two integer sweep coefficients");
        };
        Ok([x.number()?, y.number()?])
    };
    Ok(SweepRecipe {
        section: SweepSection::parse(input.field("section")?.text()?)?,
        height_mm: input.field("height_mm")?.number()?,
        top_multiplier_eighths: pair("top_multiplier_eighths")?,
        top_translation_eighths_mm: pair("top_translation_eighths_mm")?,
    })
}
pub(super) fn build(recipe: &Json) -> Result<Json, &'static str> {
    let mut keys = vec!["family"];
    keys.extend(PARAMETERS);
    recipe.keys(&keys)?;
    let parameters = parse(recipe)?;
    let blocks = build_ruled_sweep(parameters)?;
    let cert = audit_ruled_sweep(parameters, &blocks)?;
    let mut model = encode(FAMILY, blocks, &cert.certificate);
    if let Json::Object(fields) = &mut model {
        for key in PARAMETERS {
            fields.insert(key.into(), recipe.field(key)?.clone());
        }
    }
    Ok(model)
}
pub(super) fn audit(model: &Json) -> Result<Json, &'static str> {
    let mut keys = vec!["schema", "blocks", "topology", "design_faces"];
    keys.extend(PARAMETERS);
    model.keys(&keys)?;
    let parameters = parse(model)?;
    let blocks = decode_blocks(model)?;
    let result = audit_ruled_sweep(parameters, &blocks)?;
    let cert = &result.certificate;
    if &topology(&cert.incidence) != model.field("topology")? {
        return Err("incidence_mismatch");
    }
    if &design(cert) != model.field("design_faces")? {
        return Err("design_faces");
    }
    Ok(object([
        ("status", Json::String("valid".into())),
        (
            "method",
            Json::String("certified_ruled_similarity_sweep".into()),
        ),
        ("blocks", number(blocks.len())),
        ("shared_faces", number(cert.incidence.interfaces.len())),
        ("exterior_sides", number(cert.incidence.exterior.len())),
        (
            "determinant_lower_bound_mm3",
            number(cert.determinant_lower_bound_mm3),
        ),
        (
            "radial_error_upper_bound_mm",
            number(cert.radial_error_upper_mm),
        ),
        (
            "minimum_squared_scale_ratio",
            Json::Array(result.minimum_squared_scale.map(number).into()),
        ),
        (
            "sweep_cap_features",
            object([
                ("bottom", Json::String("z0".into())),
                ("top", Json::String("z8".into())),
            ]),
        ),
        ("ideal_volume_mm3", number(result.ideal_volume_mm3)),
    ]))
}
