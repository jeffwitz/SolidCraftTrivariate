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
pub(super) fn side_json(side: &BlockSide) -> Json {
    Json::Array(vec![text(&side.block), text(&side_name(side))])
}

pub(super) fn topology(incidence: &crate::incidence::Incidence) -> Json {
    let interfaces = incidence
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
            Json::Array(incidence.exterior.iter().map(side_json).collect()),
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

pub(crate) const PARAMETRIC_FAMILIES: [&str; 9] = [
    "movable_bore_v1",
    "pocket_depth_v1",
    "exact_bore_v1",
    "row_bores_v1",
    "overlap_bore_v1",
    "crescent_pocket_v1",
    "parameterized_bracket_v1",
    "travelling_bore_v1",
    "travelling_bore_bounded_v1",
];

fn parameters(family: &str) -> &'static [&'static str] {
    match family {
        "travelling_bore_v1" | "travelling_bore_bounded_v1" => &["center_x_um"],
        "parameterized_bracket_v1" => &["radii_sixteenths_mm"],
        "movable_bore_v1" => &["center_half_mm", "radius_half_mm"],
        "row_bores_v1" => &["half_widths_half_mm", "centers_half_mm", "radii_half_mm"],
        "exact_bore_v1" => &["center_um", "radius_um"],
        "pocket_depth_v1" | "crescent_pocket_v1" => &["depth_half_mm"],
        _ => &[],
    }
}

pub(crate) fn bracket_radii(input: &Json) -> Result<[u8; 2], &'static str> {
    let [x, y] = input.field("radii_sixteenths_mm")?.items()? else {
        return Err("expected two bracket radius ticks");
    };
    Ok([x.number()?, y.number()?])
}

fn bore_parameters(input: &Json, fields: &[&str]) -> Result<([i32; 2], i32), &'static str> {
    let [cx, cy] = input.field(fields[0])?.items()? else {
        return Err("invalid bore center");
    };
    Ok((
        [cx.number()?, cy.number()?],
        input.field(fields[1])?.number()?,
    ))
}

struct RowParameters {
    extents: [i32; 2],
    centers: Vec<[i32; 2]>,
    radii: Vec<i32>,
}

fn row_parameters(input: &Json) -> Result<RowParameters, &'static str> {
    fn pair(value: &Json) -> Result<[i32; 2], &'static str> {
        let [x, y] = value.items()? else {
            return Err("expected two half-mm coordinates");
        };
        Ok([x.number()?, y.number()?])
    }
    Ok(RowParameters {
        extents: pair(input.field("half_widths_half_mm")?)?,
        centers: input
            .field("centers_half_mm")?
            .items()?
            .iter()
            .map(pair)
            .collect::<Result<_, _>>()?,
        radii: input
            .field("radii_half_mm")?
            .items()?
            .iter()
            .map(Json::number)
            .collect::<Result<_, _>>()?,
    })
}

fn bounded_template(
    input: &Json,
    blocks: &[BracketBlock],
) -> Result<
    (
        crate::travelling_bore_bounded_audit::BoundedTravellingBoreCertificate,
        f64,
    ),
    &'static str,
> {
    let center = input.field("center_x_um")?.number::<i32>()?;
    let cert = crate::travelling_bore_bounded_audit::audit_travelling_bore_bounded(center, blocks)?;
    let ideal = crate::travelling_bore_events::analyze_travelling_bore(i64::from(center))?;
    Ok((cert, ideal.volume_mm3))
}

fn audit_template(
    input: &Json,
    family: &str,
    blocks: &[BracketBlock],
) -> Result<(BracketCertificate, Option<f64>), &'static str> {
    match family {
        "travelling_bore_bounded_v1" => {
            let (cert, volume) = bounded_template(input, blocks)?;
            Ok((cert.certificate, Some(volume)))
        }
        "travelling_bore_v1" => {
            let center = input.field("center_x_um")?.number::<i32>()?;
            let cert = crate::travelling_bore_audit::audit_travelling_bore(center, blocks)?;
            let physical =
                crate::travelling_bore_events::analyze_travelling_bore(i64::from(center))?;
            Ok((cert, Some(physical.volume_mm3)))
        }
        "parameterized_bracket_v1" => {
            let [rx, ry] = bracket_radii(input)?;
            let result = crate::bracket_parameters_audit::audit_bracket_parameters(rx, ry, blocks)?;
            Ok((result.certificate, Some(result.volume_mm3)))
        }
        "crescent_pocket_v1" => {
            let ticks = input.field("depth_half_mm")?.number::<i32>()?;
            let cert = crate::crescent_audit::audit_crescent_pocket(ticks, blocks)?;
            let pi = std::f64::consts::PI;
            let ideal =
                8. * (400. - 4. * pi) - f64::from(ticks) / 2. * (4. * pi / 3. + 2. * 3_f64.sqrt());
            Ok((cert, Some(ideal)))
        }
        "overlap_bore_v1" => {
            let cert = crate::overlap_bore_audit::audit_overlap(blocks)?;
            let ideal = 8. * (400. - 16. * std::f64::consts::PI / 3. - 2. * 3_f64.sqrt());
            Ok((cert, Some(ideal)))
        }
        "row_bores_v1" => {
            let row = row_parameters(input)?;
            let cert =
                crate::row_bores_audit::audit_row(row.extents, &row.centers, &row.radii, blocks)?;
            let disk_area = std::f64::consts::PI
                * row
                    .radii
                    .iter()
                    .map(|r| (f64::from(*r) / 2.).powi(2))
                    .sum::<f64>();
            let area = f64::from(row.extents[0]) * f64::from(row.extents[1]);
            Ok((cert, Some(8. * (area - disk_area))))
        }
        "movable_bore_v1" | "exact_bore_v1" => {
            let (center, ticks) = bore_parameters(input, parameters(family))?;
            let (cert, radius) = if family == "exact_bore_v1" {
                (
                    crate::exact_bore_audit::audit_exact_bore(center, ticks, blocks)?,
                    f64::from(ticks) / 1000.,
                )
            } else {
                (
                    crate::movable_bore_audit::audit_bore(center, ticks, blocks)?,
                    f64::from(ticks) / 2.,
                )
            };
            Ok((
                cert,
                Some(8. * (400. - std::f64::consts::PI * radius * radius)),
            ))
        }
        "pocket_depth_v1" => {
            let ticks = input.field("depth_half_mm")?.number::<i32>()?;
            let cert = crate::pocket_depth_audit::audit_depth(ticks, blocks)?;
            Ok((
                cert,
                Some(3200. - 8. * std::f64::consts::PI * f64::from(ticks.min(16))),
            ))
        }
        _ => certify(family, blocks),
    }
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    let family = recipe.field("family")?.text()?;
    let fields = parameters(family);
    let mut keys = vec!["family"];
    keys.extend(fields);
    recipe.keys(&keys)?;
    let blocks = match family {
        "filleted_bracket_v1" => bracket_blocks()?,
        "travelling_bore_bounded_v1" => {
            crate::travelling_bore_bounded::travelling_bore_bounded_blocks(
                recipe.field("center_x_um")?.number()?,
            )?
        }
        "travelling_bore_v1" => {
            crate::travelling_bore::travelling_bore_blocks(recipe.field("center_x_um")?.number()?)?
        }
        "parameterized_bracket_v1" => {
            let [rx, ry] = bracket_radii(recipe)?;
            crate::bracket_parameters::build_bracket_parameters(rx, ry)?.blocks
        }
        "crescent_pocket_v1" => crate::crescent_pocket::crescent_pocket_blocks(
            recipe.field("depth_half_mm")?.number()?,
        )?,
        "overlap_bore_v1" => crate::overlap_bore::overlap_blocks()?,
        "row_bores_v1" => {
            use crate::row_bores::{build_row_bores, RowBoresError::*};
            let row = row_parameters(recipe)?;
            build_row_bores(row.extents[0], row.extents[1], &row.centers, &row.radii)
                .map_err(|error| match error {
                    InvalidInput => "row_bore_invalid_input",
                    NonRepresentableBand => "row_bore_nonrepresentable_band",
                    DiskOutsideBand { .. } => "row_bore_outside_band",
                    DiskTangentToBoundary { .. } => "row_bore_boundary_contact",
                    DisksTangent { .. } => "row_bore_pair_tangent",
                    DisksIntersect { .. } => "row_bore_pair_overlap",
                    NumericalRange => "row_bore_numerical_range",
                })?
                .blocks
        }
        "movable_bore_v1" => {
            let ([x, y], r) = bore_parameters(recipe, fields)?;
            crate::movable_bore::build_movable_bore(x, y, r)
                .map_err(|_| "invalid movable bore parameters")?
                .blocks
                .into_iter()
                .map(|b| b.block)
                .collect()
        }
        "exact_bore_v1" => {
            let ([x, y], r) = bore_parameters(recipe, fields)?;
            crate::exact_bore::build_exact_bore(x, y, r)
                .map_err(|_| "invalid exact bore parameters")?
                .blocks
        }
        "pocket_depth_v1" => {
            crate::pocket_depth::build_pocket_depth(recipe.field("depth_half_mm")?.number()?)
                .map_err(|_| "invalid pocket depth")?
                .blocks
                .into_iter()
                .map(|b| b.block)
                .collect()
        }
        _ => crate::sections::section_blocks(family)?,
    };
    let (cert, _) = audit_template(recipe, family, &blocks)?;
    let mut model = encode(family, blocks, &cert);
    if let Json::Object(object) = &mut model {
        for &key in fields {
            object.insert(key.into(), recipe.field(key)?.clone());
        }
    }
    Ok(model)
}

pub(super) fn encode_blocks(blocks: Vec<BracketBlock>) -> Json {
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
    Json::Array(blocks)
}

fn encode(family: &str, blocks: Vec<BracketBlock>, cert: &BracketCertificate) -> Json {
    object([
        ("schema", text(family)),
        ("blocks", encode_blocks(blocks)),
        ("topology", topology(&cert.incidence)),
        ("design_faces", design(cert)),
    ])
}

pub(super) fn decode_blocks(model: &Json) -> Result<Vec<BracketBlock>, &'static str> {
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
    Ok(blocks)
}

pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    let family = model.field("schema")?.text()?;
    let mut keys = vec!["schema", "blocks", "topology", "design_faces"];
    keys.extend(parameters(family));
    model.keys(&keys)?;
    let blocks = decode_blocks(model)?;
    let (certificate, ideal_volume, bounded) = if family == "travelling_bore_bounded_v1" {
        let (cert, volume) = bounded_template(model, &blocks)?;
        let limits = (
            cert.radial_error_upper_mm,
            cert.hausdorff_error_upper_mm,
            cert.fidelity_budget_mm,
        );
        (cert.certificate, Some(volume), Some(limits))
    } else {
        let (cert, volume) = audit_template(model, family, &blocks)?;
        (cert, volume, None)
    };
    if model.field("topology")? != &topology(&certificate.incidence) {
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
    if let (Json::Object(fields), Some((radial, hausdorff, budget))) = (&mut report, bounded) {
        fields.insert("radial_error_upper_bound_mm".into(), number(radial));
        fields.insert("hausdorff_error_upper_bound_mm".into(), number(hausdorff));
        fields.insert("fidelity_budget_mm".into(), number(budget));
        fields.insert("fidelity_target".into(), text("ideal_clipped_bore_arcs"));
        fields.insert("representation".into(), text("quantized_rational_arcs"));
    }
    Ok(report)
}
