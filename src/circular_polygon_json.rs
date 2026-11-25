//! Exchange format for extruded polygons with exact circular through holes.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    circular_polygon::{build_circular_polygon_prism, CircularPolygonBore, CircularPolygonPrism},
    circular_polygon_audit::audit_circular_polygon_prism,
    json::{object, Json},
    polygon_json::{pair, profile},
};

pub(crate) const FAMILY: &str = "polygon_circular_bores_v1";
pub(crate) const METHOD: &str = "certified_polygon_circular_bores";
const PARAMETERS: [&str; 4] = ["outer", "bores", "height_mm", "top_shift_mm"];

fn parse(data: &Json) -> Result<CircularPolygonPrism, &'static str> {
    let list = data.field("bores")?.items()?;
    if !(1..=4).contains(&list.len()) {
        return Err("circular_polygon_bore_count");
    }
    let mut bores = Vec::with_capacity(list.len());
    for v in list {
        v.keys(&["id", "center_mm", "radius_mm", "guard_mm"])?;
        bores.push(CircularPolygonBore {
            id: v.field("id")?.text()?.to_owned(),
            center_mm: pair(v.field("center_mm")?)?,
            radius_mm: v.field("radius_mm")?.number()?,
            guard_mm: v.field("guard_mm")?.number()?,
        });
    }
    Ok(CircularPolygonPrism {
        outer: profile(data.field("outer")?)?,
        bores,
        height_mm: data.field("height_mm")?.number()?,
        top_shift_mm: pair(data.field("top_shift_mm")?)?,
    })
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&[&["family"][..], &PARAMETERS].concat())?;
    let r = parse(recipe)?;
    let blocks = build_circular_polygon_prism(&r)?;
    let certified = audit_circular_polygon_prism(&r, &blocks)?;
    let mut model = encode(FAMILY, blocks, &certified.certificate);
    if let Json::Object(fields) = &mut model {
        for key in PARAMETERS {
            fields.insert(key.into(), recipe.field(key)?.clone());
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
    let cert = audit_circular_polygon_prism(&r, &blocks)?;
    let c = &cert.certificate;
    if data.field("topology")? != &topology(&c.incidence)
        || data.field("design_faces")? != &design(c)
    {
        return Err("circular_polygon_topology_or_provenance_mismatch");
    }
    let n = |v: i64| Json::Number(v.to_string());
    Ok(object([
        ("status", Json::String("valid".into())),
        ("method", Json::String(METHOD.into())),
        ("blocks", n(blocks.len() as i64)),
        ("shared_faces", n(c.incidence.interfaces.len() as i64)),
        ("exterior_sides", n(c.incidence.exterior.len() as i64)),
        ("through_bores", n(r.bores.len() as i64)),
        (
            "determinant_lower_bound_mm3",
            Json::Number(c.determinant_lower_bound_mm3.to_string()),
        ),
        ("hausdorff_error_upper_bound_mm", n(0)),
        (
            "projected_area_exact",
            object([
                ("outer_area_twice_mm2", n(cert.outer_area_twice_mm2)),
                ("pi_coefficient_mm2", n(-cert.sum_squared_radii_mm2)),
            ]),
        ),
        ("extrusion_height_mm", n(r.height_mm.into())),
    ]))
}
