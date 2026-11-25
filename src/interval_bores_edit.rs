//! Strict versioned interchange for an atomic sequence of interval edits.
use crate::{
    interval_bores_json::{self, FAMILY},
    interval_bores_session::{IntervalBoreSession, IntervalEditReport},
    json::{object, Json},
    polygon_json::pair,
};
use std::collections::BTreeSet;

fn number(n: impl ToString) -> Json {
    Json::Number(n.to_string())
}
fn names(items: Vec<String>) -> Json {
    Json::Array(items.into_iter().map(Json::String).collect())
}

fn report(id: &str, r: IntervalEditReport) -> Json {
    let events = r
        .events
        .into_iter()
        .map(|e| {
            object([
                ("moving_cap", Json::String(e.moving_cap)),
                ("fixed_plane", Json::String(e.fixed_plane)),
                ("height_mm", number(e.height_mm)),
                ("time", Json::Array(e.time.map(number).to_vec())),
            ])
        })
        .collect();
    object([
        ("bore", Json::String(id.into())),
        ("old_class", Json::String(r.old_class.into())),
        ("new_class", Json::String(r.new_class.into())),
        ("events", Json::Array(events)),
        ("created_features", names(r.created_features)),
        ("lost_features", names(r.lost_features)),
        ("remapped_features", names(r.remapped_features)),
        ("reused_layers", number(r.reused_layers)),
        ("built_layers", number(r.built_layers)),
        ("reused_blocks", number(r.reused_blocks)),
        ("built_blocks", number(r.built_blocks)),
        ("volume_change_pi_mm3", number(r.volume_change_pi_mm3)),
        ("construction_s", number(r.build_s)),
        ("audit_s", number(r.audit_s)),
    ])
}

/// Audit the incoming model and apply every edit before returning a publishable model.
pub fn edit(model: &Json, edits: &Json) -> Result<(Json, Json), String> {
    if model.field("schema")?.text()? != FAMILY {
        return Err("interval edits require an interval bore model".into());
    }
    interval_bores_json::audit(model)?;
    edits.keys(&["schema", "steps"])?;
    if edits.field("schema")?.text()? != "interval_bore_edits_v1" {
        return Err("unsupported interval edit schema".into());
    }
    let steps = edits.field("steps")?.items()?;
    if steps.is_empty() {
        return Err("empty interval edit sequence".into());
    }
    let mut session = IntervalBoreSession::new(interval_bores_json::parse(model)?)?;
    let mut reports = Vec::with_capacity(steps.len());
    for step in steps {
        step.keys(&["bore", "z_interval_mm", "require_features"])?;
        let id = step.field("bore")?.text()?;
        let interval = pair(step.field("z_interval_mm")?)?;
        let required: Vec<String> = step
            .field("require_features")?
            .items()?
            .iter()
            .map(|v| v.text().map(str::to_owned))
            .collect::<Result<_, _>>()?;
        if required.iter().collect::<BTreeSet<_>>().len() != required.len() {
            return Err("duplicate required interval edit feature".into());
        }
        reports.push(report(id, session.apply(id, interval, &required)?));
    }
    let state = session.into_current();
    let mut parameters = object([
        ("outer", model.field("outer")?.clone()),
        ("bores", model.field("bores")?.clone()),
        ("height_mm", model.field("height_mm")?.clone()),
        ("top_shift_mm", model.field("top_shift_mm")?.clone()),
    ]);
    if let Json::Object(fields) = &mut parameters {
        let Some(Json::Array(bores)) = fields.get_mut("bores") else {
            return Err("invalid interval bore list".into());
        };
        for (data, bore) in bores.iter_mut().zip(&state.recipe.bores) {
            if let Json::Object(fields) = data {
                fields.insert(
                    "z_interval_mm".into(),
                    Json::Array(bore.z_interval_mm.map(number).to_vec()),
                );
            }
        }
    }
    let candidate = interval_bores_json::encoded_model(
        &parameters,
        state.blocks,
        &state.certificate.certificate,
    )?;
    Ok((
        candidate,
        object([
            ("status", Json::String("updated".into())),
            ("steps", Json::Array(reports)),
        ]),
    ))
}
