//! Stable JSON exchange for the fixed bracket. Exact controls stay integral.

use super::bracket::{bracket_blocks, BracketBlock};
use super::bracket_audit::{audit_bracket, BracketCertificate};
use super::incidence::BlockSide;
use super::json::{object, Json};

fn text(value: &str) -> Json {
    Json::String(value.into())
}
fn number(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}
fn side_name(side: &BlockSide) -> String {
    format!("{}{}", ["u", "v", "w"][side.axis], u8::from(side.upper))
}
fn side_json(side: &BlockSide) -> Json {
    Json::Array(vec![text(&side.block), text(&side_name(side))])
}

fn topology(cert: &BracketCertificate) -> Json {
    let interfaces = cert
        .incidence
        .interfaces
        .iter()
        .map(|face| {
            let o = face.orientation;
            object([
                ("first", side_json(&face.first)),
                ("second", side_json(&face.second)),
                (
                    "orientation",
                    object([
                        ("swap", Json::Bool(o.swap)),
                        ("flip_first", Json::Bool(o.flip_first)),
                        ("flip_second", Json::Bool(o.flip_second)),
                    ]),
                ),
            ])
        })
        .collect();
    object([
        ("interfaces", Json::Array(interfaces)),
        (
            "exterior",
            Json::Array(cert.incidence.exterior.iter().map(side_json).collect()),
        ),
    ])
}

fn design(cert: &BracketCertificate) -> Json {
    Json::Array(
        cert.design_faces
            .iter()
            .map(|(s, feature)| {
                object([
                    ("block", text(&s.block)),
                    ("side", text(&side_name(s))),
                    ("feature", text(feature)),
                ])
            })
            .collect(),
    )
}

pub(crate) fn build() -> Result<Json, &'static str> {
    let blocks = bracket_blocks()?;
    let cert = audit_bracket(&blocks)?;
    let blocks = blocks
        .into_iter()
        .map(|b| {
            object([
                ("id", text(&b.id)),
                (
                    "degrees",
                    Json::Array(b.degrees.into_iter().map(number).collect()),
                ),
                (
                    "controls",
                    Json::Array(
                        b.controls
                            .into_iter()
                            .map(|p| Json::Array(p.into_iter().map(number).collect()))
                            .collect(),
                    ),
                ),
            ])
        })
        .collect();
    Ok(object([
        ("schema", text("filleted_bracket_v1")),
        ("blocks", Json::Array(blocks)),
        ("topology", topology(&cert)),
        ("design_faces", design(&cert)),
    ]))
}

pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    model.keys(&["schema", "blocks", "topology", "design_faces"])?;
    let mut blocks = Vec::new();
    for b in model.field("blocks")?.items()? {
        b.keys(&["id", "degrees", "controls"])?;
        let degrees = b
            .field("degrees")?
            .items()?
            .iter()
            .map(Json::number::<usize>)
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .map_err(|_| "invalid degree count")?;
        let mut controls = Vec::new();
        for p in b.field("controls")?.items()? {
            controls.push(
                p.items()?
                    .iter()
                    .map(Json::number::<i64>)
                    .collect::<Result<Vec<_>, _>>()?
                    .try_into()
                    .map_err(|_| "invalid homogeneous control length")?,
            );
        }
        blocks.push(BracketBlock {
            id: b.field("id")?.text()?.into(),
            degrees,
            controls,
        });
    }
    let certificate = audit_bracket(&blocks)?;
    if model.field("topology")? != &topology(&certificate) {
        return Err("incidence_mismatch");
    }
    if model.field("design_faces")? != &design(&certificate) {
        return Err("design_faces");
    }
    Ok(object([
        ("status", text("valid")),
        ("method", text("certified_fixed_extrusion_partition")),
        ("blocks", number(blocks.len())),
        (
            "shared_faces",
            number(certificate.incidence.interfaces.len()),
        ),
        (
            "exterior_sides",
            number(certificate.incidence.exterior.len()),
        ),
        (
            "determinant_lower_bound_mm3",
            number(certificate.determinant_lower_bound_mm3),
        ),
        (
            "radial_error_upper_bound_mm",
            number(certificate.radial_error_upper_mm),
        ),
    ]))
}
