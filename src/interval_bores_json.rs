//! Versioned exchange for polygonal extrusions with bounded axial bore intervals.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    circular_polygon_json::bore_values,
    interval_bores::{build_interval_bore_prism, IntervalBore, IntervalBorePrism},
    interval_bores_audit::audit_interval_bore_prism,
    json::{object, Json},
    polygon_json::{pair, profile},
};

pub(crate) const FAMILY: &str = "polygon_interval_bores_v1";
pub(crate) const METHOD: &str = "certified_polygon_interval_bores";
const PARAMETERS: [&str; 4] = ["outer", "bores", "height_mm", "top_shift_mm"];

pub(crate) fn parse(data: &Json) -> Result<IntervalBorePrism, &'static str> {
    let list = data.field("bores")?.items()?;
    if !(1..=4).contains(&list.len()) {
        return Err("interval_bore_count");
    }
    let mut bores = Vec::with_capacity(list.len());
    for v in list {
        v.keys(&["id", "center_mm", "radius_mm", "guard_mm", "z_interval_mm"])?;
        bores.push(IntervalBore {
            bore: bore_values(v)?,
            z_interval_mm: pair(v.field("z_interval_mm")?)?,
        });
    }
    Ok(IntervalBorePrism {
        outer: profile(data.field("outer")?)?,
        bores,
        height_mm: data.field("height_mm")?.number()?,
        top_shift_mm: pair(data.field("top_shift_mm")?)?,
    })
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&[&["family"][..], &PARAMETERS].concat())?;
    let r = parse(recipe)?;
    let blocks = build_interval_bore_prism(&r)?;
    let certified = audit_interval_bore_prism(&r, &blocks)?;
    encoded_model(recipe, blocks, &certified.certificate)
}

pub(crate) fn encoded_model(
    recipe: &Json,
    blocks: Vec<crate::bracket::BracketBlock>,
    certificate: &crate::bracket_audit::BracketCertificate,
) -> Result<Json, &'static str> {
    let mut model = encode(FAMILY, blocks, certificate);
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
    let r = parse(data)?;
    let blocks = decode_blocks(data)?;
    let cert = audit_interval_bore_prism(&r, &blocks)?;
    let c = &cert.certificate;
    if data.field("topology")? != &topology(&c.incidence)
        || data.field("design_faces")? != &design(c)
    {
        return Err("interval_bore_topology_or_provenance_mismatch");
    }
    let n = |v: i64| Json::Number(v.to_string());
    Ok(object([
        ("status", Json::String("valid".into())),
        ("method", Json::String(METHOD.into())),
        ("blocks", n(blocks.len() as i64)),
        ("shared_faces", n(c.incidence.interfaces.len() as i64)),
        ("exterior_sides", n(c.incidence.exterior.len() as i64)),
        ("through_bores", n(cert.through_bores as i64)),
        ("cavities", n(cert.enclosed_voids as i64)),
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
                    n(-cert.removed_pi_volume_coefficient_mm3),
                ),
            ]),
        ),
    ]))
}
