//! Versioned interchange of integer-homogeneous annular sectors.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    json::{object, Json},
    polygon_json::pair,
    rational_sector::{build_rational_sector, RationalSector},
    rational_sector_audit::audit_rational_sector,
};

pub(crate) const FAMILY: &str = "rational_sector_prism_v1";
pub(crate) const METHOD: &str = "unisolvent_homogeneous_polar_sector";
const PARAMETERS: [&str; 3] = ["half_angle_tangent", "radii_mm", "height_mm"];

pub(crate) fn parse(value: &Json) -> Result<RationalSector, &'static str> {
    Ok(RationalSector {
        half_angle_tangent: pair(value.field("half_angle_tangent")?)?,
        radii_mm: pair(value.field("radii_mm")?)?,
        height_mm: value.field("height_mm")?.number()?,
    })
}

pub(crate) fn build(value: &Json) -> Result<Json, &'static str> {
    value.keys(&[&["family"][..], &PARAMETERS].concat())?;
    let r = parse(value)?;
    let blocks = build_rational_sector(&r)?;
    let certificate = audit_rational_sector(&r, &blocks)?;
    let mut model = encode(FAMILY, blocks, &certificate);
    if let Json::Object(fields) = &mut model {
        for k in PARAMETERS {
            fields.insert(k.to_string(), value.field(k)?.clone());
        }
    }
    Ok(model)
}

pub(crate) fn audit(value: &Json) -> Result<Json, &'static str> {
    value.keys(
        &[
            &["schema", "blocks", "topology", "design_faces"][..],
            &PARAMETERS,
        ]
        .concat(),
    )?;
    let recipe = parse(value)?;
    let blocks = decode_blocks(value)?;
    let c = audit_rational_sector(&recipe, &blocks)?;
    if value.field("topology")? != &topology(&c.incidence)
        || value.field("design_faces")? != &design(&c)
    {
        return Err("sector_incidence_or_provenance_mismatch");
    }
    let n = |v: i64| Json::Number(v.to_string());
    let [ri, ro] = recipe.radii_mm.map(i64::from);
    let [num, den] = recipe.half_angle_tangent.map(i64::from);
    let coefficient = i64::from(recipe.height_mm) * (ro * ro - ri * ri);
    let rate = (coefficient * den * den) as f64 / (num * num + den * den) as f64;
    Ok(object([
        ("status", Json::String("valid".into())),
        ("method", Json::String(METHOD.into())),
        ("blocks", n(blocks.len() as i64)),
        ("shared_faces", n(c.incidence.interfaces.len() as i64)),
        ("exterior_sides", n(c.incidence.exterior.len() as i64)),
        (
            "determinant_lower_bound_mm3",
            Json::Number(c.determinant_lower_bound_mm3.to_string()),
        ),
        ("hausdorff_error_upper_bound_mm", n(0)),
        (
            "volume_exact",
            object([
                ("atan_coefficient_mm3", n(coefficient)),
                ("argument", value.field("half_angle_tangent")?.clone()),
            ]),
        ),
        ("volume_rate_dt_mm3", Json::Number(rate.to_string())),
    ]))
}
