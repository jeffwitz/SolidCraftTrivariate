//! Stable JSON exchange for audited fixed-template block complexes. Exact controls stay integral.

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

fn certify(
    family: &str,
    blocks: &[BracketBlock],
) -> Result<(BracketCertificate, Option<f64>), &'static str> {
    if family == "filleted_bracket_v1" {
        return Ok((audit_bracket(blocks)?, None));
    }
    let certificate = crate::section_audit::audit_section(family, blocks)?;
    Ok((
        BracketCertificate {
            incidence: certificate.incidence,
            determinant_lower_bound_mm3: certificate.determinant_lower_bound_mm3,
            radial_error_upper_mm: certificate.radial_error_upper_mm,
            design_faces: certificate.design_faces,
        },
        Some(certificate.volume_mm3),
    ))
}

pub(crate) fn build(family: &str) -> Result<Json, &'static str> {
    let blocks = if family == "filleted_bracket_v1" {
        bracket_blocks()?
    } else {
        crate::sections::section_blocks(family)?
    };
    let (cert, _) = certify(family, &blocks)?;
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
        ("schema", text(family)),
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
    let (certificate, ideal_volume) = certify(model.field("schema")?.text()?, &blocks)?;
    if model.field("topology")? != &topology(&certificate) {
        return Err("incidence_mismatch");
    }
    if model.field("design_faces")? != &design(&certificate) {
        return Err("design_faces");
    }
    let mut report = object([
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
    ]);
    if let (Json::Object(fields), Some(volume)) = (&mut report, ideal_volume) {
        fields.insert("ideal_volume_mm3".into(), number(volume));
    }
    Ok(report)
}
