//! Versioned exchange for polygonal prisms with conical countersinks.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    counterbores_json::parse_bore_profiles,
    countersinks::build_countersinks,
    countersinks_audit::audit_countersinks,
    json::{object, Json},
};

pub(crate) const FAMILY: &str = "polygon_countersinks_v1";
pub(crate) const METHOD: &str = "certified_polygon_countersinks";
const PARAMETERS: [&str; 4] = ["outer", "bores", "height_mm", "top_shift_mm"];

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&[&["family"][..], &PARAMETERS].concat())?;
    let r = parse_bore_profiles(recipe, "countersink_radius_mm")?;
    let blocks = build_countersinks(&r)?;
    let certified = audit_countersinks(&r, &blocks)?;
    let mut model = encode(FAMILY, blocks, &certified.certificate);
    if let Json::Object(fields) = &mut model {
        for k in PARAMETERS {
            fields.insert(k.into(), recipe.field(k)?.clone());
        }
    }
    Ok(model)
}

pub(crate) fn audit(data: &Json) -> Result<Json, &'static str> {
    data.keys(
        &[
            &["schema", "blocks", "topology", "design_faces"][..],
            &PARAMETERS,
        ]
        .concat(),
    )?;
    let r = parse_bore_profiles(data, "countersink_radius_mm")?;
    let blocks = decode_blocks(data)?;
    let cert = audit_countersinks(&r, &blocks)?;
    let c = &cert.certificate;
    if data.field("topology")? != &topology(&c.incidence)
        || data.field("design_faces")? != &design(c)
    {
        return Err("countersink_topology_or_provenance_mismatch");
    }
    let n = |v: i64| Json::Number(v.to_string());
    Ok(object([
        ("status", Json::String("valid".into())),
        ("method", Json::String(METHOD.into())),
        ("blocks", n(blocks.len() as i64)),
        ("shared_faces", n(c.incidence.interfaces.len() as i64)),
        ("exterior_sides", n(c.incidence.exterior.len() as i64)),
        ("through_bores", n(r.bores.len() as i64)),
        ("cavities", n(0)),
        (
            "determinant_lower_bound_mm3",
            Json::Number(c.determinant_lower_bound_mm3.to_string()),
        ),
        ("hausdorff_error_upper_bound_mm", n(0)),
        (
            "volume_exact",
            object([
                (
                    "twice_rational_part_mm3",
                    n(cert.outer_area_twice_mm2 * i64::from(r.height_mm)),
                ),
                (
                    "pi_coefficient_mm3",
                    Json::Array(vec![n(-cert.removed_pi_numerator), n(3)]),
                ),
            ]),
        ),
    ]))
}
