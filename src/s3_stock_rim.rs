//! Exact section edits on the horizontal outer perimeters of an audited S3 model.
use crate::{
    bracket::BracketBlock,
    bracket_json::{decode_blocks, encode_blocks},
    interchange,
    json::{object, Json},
};
use std::collections::BTreeSet;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum RimState {
    Sharp,
    Chamfer,
    Fillet,
}

impl RimState {
    pub fn name(self) -> &'static str {
        match self {
            Self::Sharp => "sharp",
            Self::Chamfer => "chamfer",
            Self::Fillet => "fillet",
        }
    }
}

fn text(s: &str) -> Json {
    Json::String(s.to_owned())
}

fn edited_spans(
    b: &BracketBlock,
    outer: bool,
    state: RimState,
) -> Result<Vec<BracketBlock>, &'static str> {
    if b.degrees != [3, 1, 1] || b.controls.len() != 16 {
        return Err("expected S3 tensor chart");
    }
    let height = |h: [i64; 4]| {
        if h[3] <= 0 || h[2] % h[3] != 0 {
            return Err("expected integral chart end heights");
        }
        Ok(h[2] / h[3])
    };
    let lo = height(b.controls[0])?;
    let hi = height(b.controls[1])?;
    if lo < 0 || hi > 12 || lo >= hi || (lo == 0 && hi == 12) {
        return Err("unsupported S3 layer interval");
    }
    for p in b.controls.as_chunks::<2>().0 {
        if p[0][0..2] != p[1][0..2]
            || p[0][3] != p[1][3]
            || height(p[0])? != lo
            || height(p[1])? != hi
        {
            return Err("expected vertical S3 extrusion");
        }
    }
    if lo != 0 && hi != 12 {
        return Ok(vec![b.clone()]);
    }
    let bottom = lo == 0;
    let cap = if bottom { lo } else { hi };
    let sign = if bottom { 1 } else { -1 };
    let mut band = BracketBlock {
        id: b.id.clone(),
        degrees: [3, 1, 2],
        controls: Vec::new(),
    };
    let mut body = BracketBlock {
        id: b.id.clone(),
        degrees: b.degrees,
        controls: Vec::new(),
    };
    for (row, pair) in b.controls.as_chunks::<2>().0.iter().enumerate() {
        let wall = outer && row % 2 == 0;
        let [x, y, _, w] = pair[0].map(i128::from);
        let pole = |offset_eighths: i128, z_eighths: i64, weight: i128| {
            // On the stock rectangle this homothety insets all four sides
            // by offset_eighths/8 mm. Apply it only to the outer wall row.
            let offset = if wall { offset_eighths } else { 0 };
            let h = [
                400 * x + offset * (50 * w - x),
                400 * y + offset * (50 * w - 2 * y),
                50 * i128::from(z_eighths) * w,
                400 * w,
            ];
            let mut p = [0; 4];
            for (dst, value) in p.iter_mut().zip(h) {
                *dst = i64::try_from(value * weight).map_err(|_| "rim pole overflow")?;
            }
            Ok(p)
        };
        let mut section = [
            pole(0, cap * 8 + 2 * sign, 1)?,
            pole(
                i128::from(state == RimState::Chamfer),
                cap * 8
                    + if wall && state == RimState::Fillet {
                        0
                    } else {
                        sign
                    },
                1,
            )?,
            pole(if state == RimState::Sharp { 0 } else { 2 }, cap * 8, 2)?,
        ];
        if bottom {
            section.reverse();
        }
        band.controls.extend(section);
        for z in [
            lo * 8 + if bottom { 2 } else { 0 },
            hi * 8 - if bottom { 0 } else { 2 },
        ] {
            body.controls.push(pole(0, z, 1)?);
        }
    }
    Ok(if bottom {
        vec![band, body]
    } else {
        vec![body, band]
    })
}

/// Return one logical volume per parent block. Each volume lists its exact
/// Bezier spans in ascending axial order; band spans use a local unit parameter.
pub fn construct(model: &Json, state: RimState) -> Result<Json, &'static str> {
    if model.field("schema")?.text()? != "s3_face_pocket_v1" {
        return Err("expected an S3 pocket model");
    }
    let bounds = model
        .field("source_recipe")?
        .field("parent_model")?
        .field("source_recipe")?
        .field("bounds_xyz_mm")?;
    if *bounds != crate::json::parse("[[0,100],[0,50],[0,12]]")? {
        return Err("unsupported S3 stock bounds");
    }
    if interchange::audit(model)?.field("status")?.text()? != "valid" {
        return Err("invalid S3 source");
    }
    let mut outer = BTreeSet::new();
    for face in model.field("design_faces")?.items()? {
        if matches!(
            face.field("feature")?.text()?,
            "stock/xmin" | "stock/xmax" | "stock/ymin" | "stock/ymax"
        ) {
            if face.field("side")?.text()? != "v0" {
                return Err("unexpected stock side");
            }
            outer.insert(face.field("block")?.text()?);
        }
    }
    let parent = decode_blocks(model)?;
    let mut cells = Vec::new();
    let mut volumes = Vec::new();
    let mut contact_count = 0usize;
    let mut regular_bound = f64::INFINITY;
    for block in &parent {
        let mut spans = edited_spans(block, outer.contains(block.id.as_str()), state)?;
        for (i, s) in spans.iter_mut().enumerate() {
            s.id = format!("{}#{i}", block.id);
        }
        let station = |p: [i64; 4]| {
            let numerator = i128::from(p[2]) * 4;
            let w = i128::from(p[3]);
            if numerator % w != 0 {
                return Err("off-grid axial break");
            }
            Ok(Json::Number((numerator / w).to_string()))
        };
        let mut breaks = vec![station(spans[0].controls[0])?];
        for span in &spans {
            breaks.push(station(span.controls[span.degrees[2]])?);
        }
        for span in &spans {
            use crate::spatial_certificate::{
                certify_cap_contact, certify_spatial, CertificateStatus,
            };
            if state == RimState::Fillet
                && span.degrees[2] == 2
                && outer.contains(block.id.as_str())
            {
                let top = span.controls[0][2] != 0;
                certify_cap_contact(span.degrees, &span.controls, false, top, 2)?;
                contact_count += 1;
            } else {
                let proof = certify_spatial(span.degrees, &span.controls, Default::default());
                if proof.status != CertificateStatus::Valid {
                    return Err("unqualified stock-rim volume map");
                }
                regular_bound = regular_bound.min(
                    proof
                        .determinant_lower_bound_mm3
                        .ok_or("missing determinant bound")?,
                );
            }
        }
        cells.extend_from_slice(&spans);
        volumes.push(object([
            ("id", text(&block.id)),
            ("axial_breaks_quarter_mm", Json::Array(breaks)),
            ("spans", encode_blocks(spans)),
        ]));
    }
    crate::fixed_span_topology::verify(&parent, &cells)?;
    Ok(object([
        ("schema", text("s3_stock_rim_spans_v1")),
        ("state", text(state.name())),
        (
            "contact_singular_spans",
            Json::Number(contact_count.to_string()),
        ),
        (
            "determinant_lower_bound_mm3",
            Json::Number(if contact_count == 0 {
                regular_bound.to_string()
            } else {
                "0".into()
            }),
        ),
        ("parent_design_faces", model.field("design_faces")?.clone()),
        ("volumes", Json::Array(volumes)),
        ("topology", model.field("topology")?.clone()),
    ]))
}

/// Construct the edited complex and export its sewn rational boundary.
pub fn export_step(model: &Json, state: RimState) -> Result<String, &'static str> {
    let grouped = construct(model, state)?;
    let mut blocks = Vec::new();
    for volume in grouped.field("volumes")?.items()? {
        blocks.extend_from_slice(volume.field("spans")?.items()?);
    }
    let cells = decode_blocks(&object([("blocks", Json::Array(blocks))]))?;
    crate::step_boundary::export_span_boundary(cells, model, false)
}
