//! Versioned exchange for polygonal extrusions with stepped through bores.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    counterbores::{build_counterbores, Counterbore, CounterboreOpening, CounterborePrism},
    counterbores_audit::audit_counterbores,
    json::{object, Json},
    polygon_json::{pair, profile},
};

pub(crate) const FAMILY: &str = "polygon_counterbores_v1";
pub(crate) const METHOD: &str = "certified_polygon_counterbores";
const PARAMETERS: [&str; 4] = ["outer", "bores", "height_mm", "top_shift_mm"];

pub(crate) fn parse(data: &Json) -> Result<CounterborePrism, &'static str> {
    let list = data.field("bores")?.items()?;
    if !(1..=4).contains(&list.len()) {
        return Err("counterbore_count");
    }
    let mut bores = Vec::with_capacity(list.len());
    for v in list {
        v.keys(&[
            "id",
            "center_mm",
            "radius_mm",
            "guard_mm",
            "counterbore_radius_mm",
            "depth_mm",
            "opening",
        ])?;
        let opening = match v.field("opening")?.text()? {
            "start" => CounterboreOpening::Start,
            "end" => CounterboreOpening::End,
            _ => return Err("unsupported_counterbore_opening"),
        };
        bores.push(Counterbore {
            id: v.field("id")?.text()?.to_owned(),
            center_mm: pair(v.field("center_mm")?)?,
            radius_mm: v.field("radius_mm")?.number()?,
            guard_mm: v.field("guard_mm")?.number()?,
            counterbore_radius_mm: v.field("counterbore_radius_mm")?.number()?,
            depth_mm: v.field("depth_mm")?.number()?,
            opening,
        });
    }
    Ok(CounterborePrism {
        outer: profile(data.field("outer")?)?,
        bores,
        height_mm: data.field("height_mm")?.number()?,
        top_shift_mm: pair(data.field("top_shift_mm")?)?,
    })
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&[&["family"][..], &PARAMETERS].concat())?;
    let r = parse(recipe)?;
    let blocks = build_counterbores(&r)?;
    let certified = audit_counterbores(&r, &blocks)?;
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
    let r = parse(data)?;
    let blocks = decode_blocks(data)?;
    let cert = audit_counterbores(&r, &blocks)?;
    let c = &cert.certificate;
    if data.field("topology")? != &topology(&c.incidence)
        || data.field("design_faces")? != &design(c)
    {
        return Err("counterbore_topology_or_provenance_mismatch");
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
