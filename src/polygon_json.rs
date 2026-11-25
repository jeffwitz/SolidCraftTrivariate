//! Versioned polygon extrusion exchange.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    json::{object, Json},
    polygon_prism::{build_polygon_prism, PolygonPrism, Vertex},
    polygon_prism_audit::audit_polygon_prism,
};
pub(crate) const FAMILY: &str = "polygon_prism_v1";
fn n(v: impl ToString) -> Json {
    Json::Number(v.to_string())
}
fn pair(v: &Json) -> Result<[i32; 2], &'static str> {
    let [a, b] = v.items()? else {
        return Err("expected two polygon coordinates");
    };
    Ok([a.number()?, b.number()?])
}
fn parse(j: &Json) -> Result<PolygonPrism, &'static str> {
    let items = j.field("vertices")?.items()?;
    if !(3..=24).contains(&items.len()) {
        return Err("polygon_prism_range");
    }
    let mut vertices = Vec::with_capacity(items.len());
    for v in items {
        v.keys(&["id", "xy_mm"])?;
        vertices.push(Vertex {
            id: v.field("id")?.text()?.to_owned(),
            xy_mm: pair(v.field("xy_mm")?)?,
        });
    }
    Ok(PolygonPrism {
        vertices,
        height_mm: j.field("height_mm")?.number()?,
        top_shift_mm: pair(j.field("top_shift_mm")?)?,
    })
}
pub(crate) fn build(j: &Json) -> Result<Json, &'static str> {
    j.keys(&["family", "vertices", "height_mm", "top_shift_mm"])?;
    let r = parse(j)?;
    let blocks = build_polygon_prism(&r)?;
    let audit = audit_polygon_prism(&r, &blocks)?;
    let mut model = encode(FAMILY, blocks, &audit.certificate);
    if let Json::Object(fields) = &mut model {
        for key in ["vertices", "height_mm", "top_shift_mm"] {
            fields.insert(key.into(), j.field(key)?.clone());
        }
    }
    Ok(model)
}
pub(crate) fn audit(j: &Json) -> Result<Json, &'static str> {
    j.keys(&[
        "schema",
        "vertices",
        "height_mm",
        "top_shift_mm",
        "blocks",
        "topology",
        "design_faces",
    ])?;
    let r = parse(j)?;
    let blocks = decode_blocks(j)?;
    let result = audit_polygon_prism(&r, &blocks)?;
    let c = &result.certificate;
    if j.field("topology")? != &topology(&c.incidence) || j.field("design_faces")? != &design(c) {
        return Err("polygon_topology_or_provenance_mismatch");
    }
    Ok(object([
        ("status", Json::String("valid".into())),
        ("method", Json::String("certified_polygon_prism".into())),
        ("blocks", n(blocks.len())),
        ("shared_faces", n(c.incidence.interfaces.len())),
        ("exterior_sides", n(c.incidence.exterior.len())),
        (
            "determinant_lower_bound_mm3",
            n(c.determinant_lower_bound_mm3),
        ),
        (
            "projected_area_twice_mm2",
            n(result.projected_area_twice_mm2),
        ),
        ("volume_sixths_mm3", n(result.volume_sixths_mm3)),
    ]))
}
