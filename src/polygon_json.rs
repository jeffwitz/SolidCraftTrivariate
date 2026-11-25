//! Versioned polygon extrusion exchange.
use crate::{
    bracket::BracketBlock,
    bracket_json::{decode_blocks, design, encode, topology},
    json::{object, Json},
    perforated_prism::{build_perforated_prism, LoopProfile, PerforatedPrism},
    perforated_prism_audit::audit_perforated_prism,
    polygon_prism::{build_polygon_prism, PolygonPrism, Vertex},
    polygon_prism_audit::{audit_polygon_prism, PolygonCertificate},
};
pub(crate) const FAMILY: &str = "polygon_prism_v1";
pub(crate) const PERFORATED: &str = "polygonal_perforated_prism_v1";
fn n(v: impl ToString) -> Json {
    Json::Number(v.to_string())
}
pub(super) fn pair(v: &Json) -> Result<[i32; 2], &'static str> {
    let [a, b] = v.items()? else {
        return Err("expected two polygon coordinates");
    };
    Ok([a.number()?, b.number()?])
}
fn vertices(j: &Json, max: usize) -> Result<Vec<Vertex>, &'static str> {
    let items = j.items()?;
    if !(3..=max).contains(&items.len()) {
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
    Ok(vertices)
}
pub(super) fn profile(j: &Json) -> Result<LoopProfile, &'static str> {
    j.keys(&["id", "vertices"])?;
    Ok(LoopProfile {
        id: j.field("id")?.text()?.into(),
        vertices: vertices(j.field("vertices")?, 32)?,
    })
}
enum Prism {
    Simple(PolygonPrism),
    Perforated(PerforatedPrism),
}
impl Prism {
    fn parse(j: &Json, family: &str) -> Result<Self, &'static str> {
        let height_mm = j.field("height_mm")?.number()?;
        let top_shift_mm = pair(j.field("top_shift_mm")?)?;
        match family {
            FAMILY => Ok(Self::Simple(PolygonPrism {
                vertices: vertices(j.field("vertices")?, 24)?,
                height_mm,
                top_shift_mm,
            })),
            PERFORATED => {
                let holes = j.field("holes")?.items()?;
                if holes.len() > 4 {
                    return Err("too many polygon holes");
                }
                Ok(Self::Perforated(PerforatedPrism {
                    outer: profile(j.field("outer")?)?,
                    holes: holes.iter().map(profile).collect::<Result<_, _>>()?,
                    height_mm,
                    top_shift_mm,
                }))
            }
            _ => Err("unknown polygon extrusion family"),
        }
    }
    fn build(&self) -> Result<Vec<BracketBlock>, &'static str> {
        match self {
            Self::Simple(r) => build_polygon_prism(r),
            Self::Perforated(r) => build_perforated_prism(r),
        }
    }
    fn audit(&self, blocks: &[BracketBlock]) -> Result<PolygonCertificate, &'static str> {
        match self {
            Self::Simple(r) => audit_polygon_prism(r, blocks),
            Self::Perforated(r) => audit_perforated_prism(r, blocks),
        }
    }
}
fn keys(family: &str) -> Result<&'static [&'static str], &'static str> {
    match family {
        FAMILY => Ok(&["vertices", "height_mm", "top_shift_mm"]),
        PERFORATED => Ok(&["outer", "holes", "height_mm", "top_shift_mm"]),
        _ => Err("unknown polygon extrusion family"),
    }
}
pub(crate) fn build(j: &Json) -> Result<Json, &'static str> {
    let family = j.field("family")?.text()?;
    let fields = keys(family)?;
    j.keys(&[&["family"][..], fields].concat())?;
    let r = Prism::parse(j, family)?;
    let blocks = r.build()?;
    let audit = r.audit(&blocks)?;
    let mut model = encode(family, blocks, &audit.certificate);
    if let Json::Object(map) = &mut model {
        for &key in fields {
            map.insert(key.into(), j.field(key)?.clone());
        }
    }
    Ok(model)
}
pub(crate) fn audit(j: &Json) -> Result<Json, &'static str> {
    let family = j.field("schema")?.text()?;
    j.keys(
        &[
            &["schema", "blocks", "topology", "design_faces"][..],
            keys(family)?,
        ]
        .concat(),
    )?;
    let r = Prism::parse(j, family)?;
    let blocks = decode_blocks(j)?;
    let result = r.audit(&blocks)?;
    let method = match r {
        Prism::Simple(_) => "certified_polygon_prism",
        Prism::Perforated(_) => "certified_perforated_polygon_prism",
    };
    let c = &result.certificate;
    if j.field("topology")? != &topology(&c.incidence) || j.field("design_faces")? != &design(c) {
        return Err("polygon_topology_or_provenance_mismatch");
    }
    Ok(object([
        ("status", Json::String("valid".into())),
        ("method", Json::String(method.into())),
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
