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
    Ok(encode(family, blocks, &cert))
}

fn encode(family: &str, blocks: Vec<BracketBlock>, cert: &BracketCertificate) -> Json {
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
    object([
        ("schema", text(family)),
        ("blocks", Json::Array(blocks)),
        ("topology", topology(cert)),
        ("design_faces", design(cert)),
    ])
}

fn bore_parameters(input: &Json) -> Result<([i32; 2], i32), &'static str> {
    let [cx, cy] = input.field("center_half_mm")?.items()? else {
        return Err("invalid bore center");
    };
    Ok((
        [cx.number()?, cy.number()?],
        input.field("radius_half_mm")?.number()?,
    ))
}

pub(crate) fn build_bore(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&["family", "center_half_mm", "radius_half_mm"])?;
    let (center, radius) = bore_parameters(recipe)?;
    let candidate = crate::movable_bore::build_movable_bore(center[0], center[1], radius)
        .map_err(|_| "invalid movable bore parameters")?;
    let blocks: Vec<_> = candidate.blocks.into_iter().map(|b| b.block).collect();
    let cert = crate::movable_bore_audit::audit_bore(center, radius, &blocks)?;
    let mut model = encode("movable_bore_v1", blocks, &cert);
    if let Json::Object(fields) = &mut model {
        fields.insert(
            "center_half_mm".into(),
            Json::Array(center.into_iter().map(number).collect()),
        );
        fields.insert("radius_half_mm".into(), number(radius));
    }
    Ok(model)
}

pub(crate) fn build_depth(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&["family", "depth_half_mm"])?;
    let depth = recipe.field("depth_half_mm")?.number::<i32>()?;
    let candidate =
        crate::pocket_depth::build_pocket_depth(depth).map_err(|_| "invalid pocket depth")?;
    let blocks: Vec<_> = candidate.blocks.into_iter().map(|b| b.block).collect();
    let cert = crate::pocket_depth_audit::audit_depth(depth, &blocks)?;
    let mut model = encode("pocket_depth_v1", blocks, &cert);
    if let Json::Object(fields) = &mut model {
        fields.insert("depth_half_mm".into(), number(depth));
    }
    Ok(model)
}

pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    let family = model.field("schema")?.text()?;
    if family == "pocket_depth_v1" {
        model.keys(&[
            "schema",
            "blocks",
            "topology",
            "design_faces",
            "depth_half_mm",
        ])?;
    } else if family == "movable_bore_v1" {
        model.keys(&[
            "schema",
            "blocks",
            "topology",
            "design_faces",
            "center_half_mm",
            "radius_half_mm",
        ])?;
    } else {
        model.keys(&["schema", "blocks", "topology", "design_faces"])?;
    }
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
    let (certificate, ideal_volume) = if family == "movable_bore_v1" {
        let (center, ticks) = bore_parameters(model)?;
        let radius = f64::from(ticks) * 0.5;
        (
            crate::movable_bore_audit::audit_bore(center, ticks, &blocks)?,
            Some(8. * (400. - std::f64::consts::PI * radius * radius)),
        )
    } else if family == "pocket_depth_v1" {
        let ticks = model.field("depth_half_mm")?.number::<i32>()?;
        let certificate = crate::pocket_depth_audit::audit_depth(ticks, &blocks)?;
        let depth = 0.5 * f64::from(ticks.min(16));
        (
            certificate,
            Some(3200. - 16. * std::f64::consts::PI * depth),
        )
    } else {
        certify(family, &blocks)?
    };
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
