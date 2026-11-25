//! Audited persistence of automatically partitioned multi-contour NURBS sketches.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    bracket_json::encode,
    incidence::enumerate_faces,
    json::{object, Json},
    nurbs_extrusion::determinant_bound_at,
    nurbs_multiloop::construct,
};
use std::collections::BTreeSet;

pub(crate) const FAMILY: &str = "nurbs_multiloop_v1";
pub(crate) const EDGE_FAMILY: &str = "nurbs_edge_loops_v1";
pub(crate) const METHOD: &str = "ordered_monotone_nurbs_contours_and_axial_layers";
const PARAMETERS: [&str; 3] = ["bounds_xy_mm", "height_mm", "loops"];
const EDGE_PARAMETERS: [&str; 4] = ["bounds_xy_mm", "height_mm", "edges", "cuts"];

fn number(v: impl ToString) -> Json {
    Json::Number(v.to_string())
}

fn bound(b: &BracketBlock) -> Result<f64, &'static str> {
    let p = b.controls.first().ok_or("empty multiloop block")?;
    if p[3] <= 0 || p[2] % p[3] != 0 {
        return Err("invalid multiloop layer altitude");
    }
    determinant_bound_at(b, p[2] / p[3])
}

fn candidate(recipe: &Json) -> Result<(Json, f64), &'static str> {
    let family = recipe.field("family")?.text()?;
    let normalized;
    let contours = match family {
        FAMILY => recipe,
        EDGE_FAMILY => {
            normalized = crate::nurbs_edge_loops::normalize(recipe)?;
            &normalized
        }
        _ => return Err("unknown multiloop recipe"),
    };
    let (blocks, mut faces) = construct(contours)?;
    let incidence = enumerate_faces(&blocks)?;
    faces.sort_by(|a, b| a.0.cmp(&b.0));
    let labelled: BTreeSet<_> = faces.iter().map(|(s, _)| s).collect();
    if faces.len() != labelled.len() || labelled != incidence.exterior.iter().collect() {
        return Err("incomplete multiloop boundary provenance");
    }
    let lower = blocks
        .iter()
        .map(bound)
        .try_fold(f64::INFINITY, |lowest, b| b.map(|v| lowest.min(v)))?;
    if !lower.is_finite() || lower <= 0. {
        return Err("unresolved multiloop determinant bound");
    }
    let cert = BracketCertificate {
        incidence,
        design_faces: faces,
        determinant_lower_bound_mm3: lower,
        radial_error_upper_mm: 0.,
    };
    let mut model = encode(family, blocks, &cert);
    if let Json::Object(fields) = &mut model {
        for key in PARAMETERS {
            fields.insert(key.into(), contours.field(key)?.clone());
        }
        if family == EDGE_FAMILY {
            for key in ["edges", "cuts"] {
                fields.insert(key.into(), recipe.field(key)?.clone());
            }
        }
    }
    Ok((model, lower))
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    candidate(recipe).map(|(model, _)| model)
}

pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    let family = model.field("schema")?.text()?;
    let parameters = match family {
        FAMILY => PARAMETERS.as_slice(),
        EDGE_FAMILY => EDGE_PARAMETERS.as_slice(),
        _ => return Err("unknown multiloop schema"),
    };
    let mut recipe = object([("family", Json::String(family.into()))]);
    if let Json::Object(fields) = &mut recipe {
        for &key in parameters {
            fields.insert(key.into(), model.field(key)?.clone());
        }
    }
    let (expected, lower) = candidate(&recipe)?;
    if *model != expected {
        return Err("multiloop_stored_partition_mismatch");
    }
    let height = model.field("height_mm")?.number::<i64>()?;
    let mut cavities = 0;
    let mut through = 0;
    for contour in model.field("loops")?.items()? {
        let [a, b] = contour.field("z_interval_mm")?.items()? else {
            return Err("invalid multiloop axial interval");
        };
        let (lo, hi) = (a.number::<i64>()?, b.number::<i64>()?);
        cavities += usize::from(lo > 0 && hi < height);
        through += usize::from(lo == 0 && hi == height);
    }
    let topology = model.field("topology")?;
    Ok(object([
        ("method", Json::String(METHOD.into())),
        ("status", Json::String("valid".into())),
        ("blocks", number(model.field("blocks")?.items()?.len())),
        (
            "shared_faces",
            number(topology.field("interfaces")?.items()?.len()),
        ),
        (
            "exterior_sides",
            number(topology.field("exterior")?.items()?.len()),
        ),
        ("determinant_lower_bound_mm3", number(lower)),
        ("hausdorff_error_upper_bound_mm", number(0)),
        ("cavities", number(cavities)),
        ("through_bores", number(through)),
    ]))
}
