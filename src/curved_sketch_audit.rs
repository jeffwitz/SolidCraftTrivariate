//! Exact replay of a curved sketch followed by closed-cube regularity and trace checks.
use crate::bracket_json::{decode_blocks, topology};
use crate::incidence::enumerate_faces;
use crate::json::{object, Json};
use crate::spatial_certificate::{certify_spatial, CertificateBudget, CertificateStatus};

pub(crate) const FAMILY: &str = "curved_sketch_layers_v1";
pub(crate) const METHOD: &str = "curved_sketch_trace_and_bernstein";

fn number(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}

fn certify(model: &Json) -> Result<Json, &'static str> {
    let blocks = decode_blocks(model)?;
    let faces = enumerate_faces(&blocks)?;
    if &topology(&faces) != model.field("topology")? {
        return Err("curved_sketch_trace_mismatch");
    }
    let mut minimum = f64::INFINITY;
    for block in &blocks {
        let certificate =
            certify_spatial(block.degrees, &block.controls, CertificateBudget::default());
        match certificate.status {
            CertificateStatus::Invalid => return Err("curved_sketch_nonpositive_jacobian"),
            CertificateStatus::Indeterminate => {
                return Err("curved_sketch_indeterminate_bernstein")
            }
            CertificateStatus::Valid => {
                minimum = minimum.min(
                    certificate
                        .determinant_lower_bound_mm3
                        .ok_or("missing_jacobian_bound")?,
                );
            }
        }
    }
    if !minimum.is_finite() || minimum <= 0.0 {
        return Err("missing_jacobian_bound");
    }
    // Inputs to this family are disjoint, strictly internal feature footprints.
    let source = model.field("scaffold")?;
    let height = source.field("height_mm")?.number::<i32>()?;
    let mut cavities = 0;
    for cut in source.field("cuts")?.items()? {
        let z = cut.field("z_interval_mm")?.items()?;
        if let [start, end] = z {
            cavities += usize::from(start.number::<i32>()? > 0 && end.number::<i32>()? < height);
        } else {
            return Err("invalid_sketch_cut_interval");
        }
    }
    Ok(object([
        ("status", Json::String("valid".into())),
        ("method", Json::String(METHOD.into())),
        ("blocks", number(blocks.len())),
        ("shared_faces", number(faces.interfaces.len())),
        ("exterior_sides", number(faces.exterior.len())),
        ("cavities", number(cavities)),
        ("determinant_lower_bound_mm3", number(minimum)),
    ]))
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    let model = crate::curved_sketch::build(recipe)?;
    certify(&model)?;
    Ok(model)
}

pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    let recipe = object([
        ("family", Json::String(FAMILY.into())),
        ("scaffold", model.field("scaffold")?.clone()),
    ]);
    if crate::curved_sketch::build(&recipe)? != *model {
        return Err("curved_sketch_replay_mismatch");
    }
    certify(model)
}
