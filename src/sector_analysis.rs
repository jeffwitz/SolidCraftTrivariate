//! A stable two-chart analysis space for a sector whose geometric charts can change.
use crate::{
    bracket::BracketBlock,
    incidence::BlockSide,
    json::{object, Json},
};
use std::collections::BTreeMap;

pub(crate) fn unpack<'a>(model: &Json, job: &'a Json) -> Result<(&'a Json, bool), &'static str> {
    if job.field("schema")?.text()? != "solidcraft_sector_analysis_job_v1" {
        return Ok((job, false));
    }
    job.keys(&["schema", "analysis_partition", "job"])?;
    if model.field("schema")?.text()? != crate::rational_sector_json::FAMILY {
        return Err("fixed angular partition requires a rational sector");
    }
    if job.field("analysis_partition")?.text()? != "angular_halves" {
        return Err("unsupported sector analysis partition");
    }
    let inner = job.field("job")?;
    if inner.field("schema")?.text()? != "solidcraft_mechanical_job_v1" {
        return Err("sector analysis wrapper requires a v1 mechanical job");
    }
    Ok((inner, true))
}

/// Remap physical sides only; the new interior angular interface has no design face.
pub(crate) fn fix_partition(
    blocks: &mut Vec<BracketBlock>,
    groups: &mut BTreeMap<BlockSide, String>,
) -> Result<(), &'static str> {
    if blocks.len() == 2 {
        return Ok(());
    }
    let [single] = blocks.as_slice() else {
        return Err("unexpected source angular partition");
    };
    let children = crate::rational_sector::split_for_analysis(single)?;
    let mut next = BTreeMap::new();
    for (side, tag) in groups.iter() {
        for (index, child) in children.iter().enumerate() {
            if side.axis != 0 || usize::from(side.upper) == index {
                next.insert(
                    BlockSide {
                        block: child.id.clone(),
                        axis: side.axis,
                        upper: side.upper,
                    },
                    tag.clone(),
                );
            }
        }
    }
    *blocks = children;
    *groups = next;
    Ok(())
}

pub(crate) fn annotate(
    result: &mut Json,
    model: &Json,
    wrapper: &Json,
) -> Result<(), &'static str> {
    let Json::Object(fields) = result else {
        return Err("invalid mechanical result");
    };
    if fields.contains_key("job") {
        fields.insert("job".into(), wrapper.clone());
    }
    fields.insert(
        "analysis_partition".into(),
        object([
            ("policy", Json::String("angular_halves".into())),
            (
                "source_charts",
                Json::Number(model.field("blocks")?.items()?.len().to_string()),
            ),
            ("analysis_root_charts", Json::Number("2".into())),
        ]),
    );
    Ok(())
}
