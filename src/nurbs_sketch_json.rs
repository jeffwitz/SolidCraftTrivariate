//! Native exchange for two synchronized, exact NURBS sketch boundaries.
use crate::{
    bracket_json::{encode_blocks, topology},
    incidence::enumerate_faces,
    json::{object, Json},
    nurbs_curve::{extract_bezier, BezierSpan},
    nurbs_extrusion::extrude_strips,
};
pub(crate) const FAMILY: &str = "nurbs_sketch_strip_v1";
pub(crate) const METHOD: &str = "exact_monotone_nurbs_strip_extrusion";

fn n(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}

pub(crate) fn spans(curve: &Json) -> Result<Vec<BezierSpan>, &'static str> {
    curve.keys(&["degree", "knots", "controls"])?;
    let degree = curve.field("degree")?.number()?;
    let knots = curve
        .field("knots")?
        .items()?
        .iter()
        .map(Json::number)
        .collect::<Result<Vec<i64>, _>>()?;
    let mut controls = Vec::new();
    for point in curve.field("controls")?.items()? {
        let [x, y, w] = point.items()? else {
            return Err("expected planar homogeneous triple");
        };
        controls.push([x.number()?, y.number()?, w.number()?]);
    }
    extract_bezier(degree, &knots, &controls)
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&["family", "lower", "upper", "height_mm"])?;
    let lower = spans(recipe.field("lower")?)?;
    let upper = spans(recipe.field("upper")?)?;
    let blocks = extrude_strips(&lower, &upper, recipe.field("height_mm")?.number()?)?;
    let links = enumerate_faces(&blocks)?;
    if links.interfaces.len() + 1 != blocks.len() {
        return Err("nurbs_strip_nonconforming_chain");
    }
    let faces = links
        .exterior
        .iter()
        .map(|s| {
            let (side, feature) = match (s.axis, s.upper) {
                (0, false) => ("u0", "left"),
                (0, true) => ("u1", "right"),
                (1, false) => ("v0", "lower"),
                (1, true) => ("v1", "upper"),
                (2, false) => ("w0", "start"),
                _ => ("w1", "end"),
            };
            object([
                ("block", Json::String(s.block.clone())),
                ("side", Json::String(side.into())),
                ("feature", Json::String(feature.into())),
            ])
        })
        .collect();
    let first = lower[0].controls[0];
    let last = *lower.last().unwrap().controls.last().unwrap();
    let decreasing =
        i128::from(first[0]) * i128::from(last[2]) > i128::from(last[0]) * i128::from(first[2]);
    let mut intervals: Vec<_> = lower.iter().map(|s| s.interval).collect();
    if decreasing {
        intervals.reverse();
        for interval in &mut intervals {
            interval.swap(0, 1);
        }
    }
    Ok(object([
        ("schema", Json::String(FAMILY.into())),
        ("lower", recipe.field("lower")?.clone()),
        ("upper", recipe.field("upper")?.clone()),
        ("height_mm", recipe.field("height_mm")?.clone()),
        (
            "span_intervals",
            Json::Array(
                intervals
                    .into_iter()
                    .map(|interval| Json::Array(interval.map(n).into()))
                    .collect(),
            ),
        ),
        ("blocks", encode_blocks(blocks)),
        ("topology", topology(&links)),
        ("design_faces", Json::Array(faces)),
    ]))
}

pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    let recipe = object([
        ("family", Json::String(FAMILY.into())),
        ("lower", model.field("lower")?.clone()),
        ("upper", model.field("upper")?.clone()),
        ("height_mm", model.field("height_mm")?.clone()),
    ]);
    if &build(&recipe)? != model {
        return Err("nurbs_strip_exact_reconstruction_mismatch");
    }
    let bound = crate::bracket_json::decode_blocks(model)?
        .iter()
        .map(crate::nurbs_extrusion::determinant_bound)
        .try_fold(f64::INFINITY, |lo, bound| bound.map(|b| lo.min(b)))?;
    let topology = model.field("topology")?;
    Ok(object([
        ("method", Json::String(METHOD.into())),
        ("status", Json::String("valid".into())),
        ("blocks", n(model.field("blocks")?.items()?.len())),
        (
            "shared_faces",
            n(topology.field("interfaces")?.items()?.len()),
        ),
        (
            "exterior_sides",
            n(topology.field("exterior")?.items()?.len()),
        ),
        ("cavities", n(0)),
        ("determinant_lower_bound_mm3", n(bound)),
        ("hausdorff_error_upper_bound_mm", n(0)),
        (
            "injectivity",
            Json::String("strict_ordered_abscissa_positive_strip_gap_positive_height".into()),
        ),
    ]))
}
