//! Deterministic exact NURBS-cover planning for audited gallery volume maps.
//!
//! This planner selects compatible source-cell chains and pairs. The output
//! layout is always rebuilt by `nurbs_tensor_assembly::assemble_cover`; the
//! planner never accepts or copies proposed control nets.

use crate::json::{object, Json};
use std::collections::{BTreeMap, BTreeSet};

type Result<T> = std::result::Result<T, &'static str>;

#[derive(Clone)]
struct SpanPlan {
    source: String,
    axis_map: [usize; 3],
    flips: [bool; 3],
    interval: [usize; 2],
}

#[derive(Clone)]
struct Group {
    id: String,
    spans: Vec<SpanPlan>,
}

fn num(n: usize) -> Json {
    Json::Number(n.to_string())
}

fn text(s: impl Into<String>) -> Json {
    Json::String(s.into())
}

fn array(values: impl IntoIterator<Item = Json>) -> Json {
    Json::Array(values.into_iter().collect())
}

fn triple(value: &Json) -> Result<[usize; 3]> {
    value
        .items()?
        .iter()
        .map(Json::number)
        .collect::<Result<Vec<_>>>()?
        .try_into()
        .map_err(|_| "expected three degrees")
}

fn side(value: &Json) -> Result<(String, usize, usize)> {
    let parts = value.items()?;
    if parts.len() != 2 {
        return Err("invalid source side");
    }
    let name = parts[0].text()?.to_owned();
    let (axis, endpoint) = side_name(parts[1].text()?)?;
    Ok((name, axis, endpoint))
}

fn side_name(name: &str) -> Result<(usize, usize)> {
    let bytes = name.as_bytes();
    if bytes.len() != 2 || !b"uvw".contains(&bytes[0]) || !b"01".contains(&bytes[1]) {
        return Err("invalid source side");
    }
    Ok((
        b"uvw"
            .iter()
            .position(|x| *x == bytes[0])
            .ok_or("invalid side axis")?,
        (bytes[1] - b'0') as usize,
    ))
}

fn bools(value: [bool; 3]) -> Json {
    array(value.into_iter().map(Json::Bool))
}

fn indices(value: [usize; 3]) -> Json {
    array(value.into_iter().map(num))
}

fn interval(value: [usize; 2]) -> Json {
    array(value.into_iter().map(num))
}

fn oriented_span(
    span: &SpanPlan,
    target: [usize; 3],
    source: &BTreeMap<String, [usize; 3]>,
) -> Result<Json> {
    let source_degree = source.get(&span.source).ok_or("unknown source block")?;
    let elevation = std::array::from_fn(|axis| {
        target[axis].saturating_sub(source_degree[span.axis_map[axis]])
    });
    Ok(object([
        ("source_block", text(span.source.clone())),
        ("axis_map", indices(span.axis_map)),
        ("flips", bools(span.flips)),
        ("degree_elevation", indices(elevation)),
        ("parameter_interval", interval(span.interval)),
    ]))
}

fn permutation_parity(axes: [usize; 3], flips: [bool; 3]) -> bool {
    let inversions = (0..3)
        .flat_map(|i| (i + 1..3).map(move |j| (i, j)))
        .filter(|&(i, j)| axes[i] > axes[j])
        .count();
    (inversions + flips.into_iter().filter(|x| *x).count()) % 2 == 0
}

fn source_degrees(model: &Json) -> Result<BTreeMap<String, [usize; 3]>> {
    let mut output = BTreeMap::new();
    for block in model.field("blocks")?.items()? {
        let id = block.field("id")?.text()?.to_owned();
        let degree = triple(block.field("degrees")?)?;
        if output.insert(id, degree).is_some() {
            return Err("duplicate source block id");
        }
    }
    if output.is_empty() {
        return Err("model has no source blocks");
    }
    Ok(output)
}

fn single_group(id: &str) -> Group {
    Group {
        id: id.to_owned(),
        spans: vec![SpanPlan {
            source: id.to_owned(),
            axis_map: [0, 1, 2],
            flips: [false; 3],
            interval: [0, 1],
        }],
    }
}

fn chain_groups(model: &Json, degrees: &BTreeMap<String, [usize; 3]>) -> Result<Vec<Group>> {
    let topology = model.field("topology")?;
    let interfaces = topology.field("interfaces")?.items()?;
    let mut pair_count = BTreeMap::<(String, String), usize>::new();
    for iface in interfaces {
        let (a, _, _) = side(iface.field("first")?)?;
        let (b, _, _) = side(iface.field("second")?)?;
        let key = if a <= b { (a, b) } else { (b, a) };
        *pair_count.entry(key).or_default() += 1;
    }
    let mut forward = BTreeMap::<String, String>::new();
    let mut predecessor = BTreeMap::<String, String>::new();
    for iface in interfaces {
        let orientation = iface.field("orientation")?;
        if ["swap", "flip_first", "flip_second"].iter().any(|key| {
            orientation
                .field(key)
                .ok()
                .and_then(|x| x.bool().ok())
                .unwrap_or(false)
        }) {
            continue;
        }
        let (a, axis_a, end_a) = side(iface.field("first")?)?;
        let (b, axis_b, end_b) = side(iface.field("second")?)?;
        let (lo, hi) = if axis_a == 2 && end_a == 1 && axis_b == 2 && end_b == 0 {
            (a, b)
        } else if axis_a == 2 && end_a == 0 && axis_b == 2 && end_b == 1 {
            (b, a)
        } else {
            continue;
        };
        let pair_key = if lo <= hi {
            (lo.clone(), hi.clone())
        } else {
            (hi.clone(), lo.clone())
        };
        if pair_count.get(&pair_key).copied().unwrap_or(0) != 1 {
            continue;
        }
        let Some(dl) = degrees.get(&lo) else { continue };
        let Some(dh) = degrees.get(&hi) else { continue };
        if dl != dh {
            continue;
        }
        if forward.contains_key(&lo) || predecessor.contains_key(&hi) {
            continue;
        }
        forward.insert(lo.clone(), hi.clone());
        predecessor.insert(hi, lo);
    }
    let mut result = Vec::new();
    let starts: Vec<_> = forward
        .keys()
        .filter(|id| !predecessor.contains_key(*id))
        .cloned()
        .collect();
    for start in starts {
        let mut ids = vec![start.clone()];
        let mut current = start;
        let mut seen = BTreeSet::from([current.clone()]);
        while let Some(next) = forward.get(&current) {
            if !seen.insert(next.clone()) {
                break;
            }
            ids.push(next.clone());
            current = next.clone();
        }
        if ids.len() < 2 {
            continue;
        }
        let spans = ids
            .iter()
            .enumerate()
            .map(|(i, source)| SpanPlan {
                source: source.clone(),
                axis_map: [2, 0, 1],
                flips: [false; 3],
                interval: [i, i + 1],
            })
            .collect();
        result.push(Group {
            id: format!("axial:{}", ids.join("+")),
            spans,
        });
    }
    Ok(result)
}

fn pair_group(
    a: &str,
    b: &str,
    iface: &Json,
    tangent_swap_a: bool,
    tangent_flips_a: [bool; 2],
) -> Option<Group> {
    let (_, normal_a, end_a) = side(iface.field("first").ok()?).ok()?;
    let (_, normal_b, end_b) = side(iface.field("second").ok()?).ok()?;
    let orientation = iface.field("orientation").ok();
    let flag = |key: &str| {
        orientation
            .and_then(|o| o.field(key).ok())
            .and_then(|v| v.bool().ok())
            .unwrap_or(false)
    };
    let swap = flag("swap");
    let face_flips = [flag("flip_first"), flag("flip_second")];
    let mut tangents_a: Vec<_> = (0..3).filter(|x| *x != normal_a).collect();
    if tangent_swap_a {
        tangents_a.swap(0, 1);
    }
    let tangents_b: Vec<_> = tangents_a
        .iter()
        .map(|axis| {
            let index_a = (0..3).filter(|x| *x != normal_a).position(|x| x == *axis)?;
            let index_b = if swap { 1 - index_a } else { index_a };
            (0..3).filter(|x| *x != normal_b).nth(index_b)
        })
        .collect::<Option<Vec<_>>>()?;
    let axes_a = [normal_a, tangents_a[0], tangents_a[1]];
    let axes_b = [normal_b, tangents_b[0], tangents_b[1]];
    let flips_a = [end_a == 0, tangent_flips_a[0], tangent_flips_a[1]];
    let flips_b = [
        end_b == 1,
        tangent_flips_a[0]
            ^ face_flips
                [usize::from(tangents_a[0] == (0..3).filter(|x| *x != normal_a).nth(1).unwrap())],
        tangent_flips_a[1]
            ^ face_flips
                [usize::from(tangents_a[1] == (0..3).filter(|x| *x != normal_a).nth(1).unwrap())],
    ];
    if !permutation_parity(axes_a, flips_a) || !permutation_parity(axes_b, flips_b) {
        return None;
    }
    Some(Group {
        id: format!("band:{a}+{b}"),
        spans: vec![
            SpanPlan {
                source: a.to_owned(),
                axis_map: axes_a,
                flips: flips_a,
                interval: [0, 1],
            },
            SpanPlan {
                source: b.to_owned(),
                axis_map: axes_b,
                flips: flips_b,
                interval: [1, 2],
            },
        ],
    })
}

fn group_degrees(group: &Group, source: &BTreeMap<String, [usize; 3]>) -> Result<[usize; 3]> {
    let mut out = [1; 3];
    for span in &group.spans {
        let degree = source.get(&span.source).ok_or("unknown source block")?;
        for axis in 0..3 {
            out[axis] = out[axis].max(degree[span.axis_map[axis]]);
        }
    }
    Ok(out)
}

fn subface(span: &SpanPlan, patch_id: &str, side_name: &str) -> Result<Json> {
    let bytes = side_name.as_bytes();
    if bytes.len() != 2 || !b"uvw".contains(&bytes[0]) || !b"01".contains(&bytes[1]) {
        return Err("invalid interface side");
    }
    let source_axis = b"uvw".iter().position(|x| *x == bytes[0]).unwrap();
    let patch_axis = span
        .axis_map
        .iter()
        .position(|x| *x == source_axis)
        .ok_or("invalid axis map")?;
    let endpoint = (bytes[1] - b'0') as usize ^ usize::from(span.flips[patch_axis]);
    Ok(object([
        ("patch_id", text(patch_id.to_owned())),
        ("source_block", text(span.source.clone())),
        ("side", text(side_name.to_owned())),
        (
            "patch_side",
            text(format!("{}{endpoint}", b"uvw"[patch_axis] as char)),
        ),
        ("axis_map", indices(span.axis_map)),
        ("flips", bools(span.flips)),
        ("parameter_interval", interval(span.interval)),
    ]))
}

fn cover_layout(
    model: &Json,
    groups: &[Group],
    degrees: &BTreeMap<String, [usize; 3]>,
) -> Result<Json> {
    let mut owner = BTreeMap::<String, (&str, &SpanPlan)>::new();
    let mut patch_values = Vec::new();
    for group in groups {
        let mut source_ids = Vec::new();
        for span in &group.spans {
            source_ids.push(text(span.source.clone()));
            if owner
                .insert(span.source.clone(), (group.id.as_str(), span))
                .is_some()
            {
                return Err("source block assigned more than once");
            }
        }
        let degree = group_degrees(group, degrees)?;
        let spans = group
            .spans
            .iter()
            .map(|span| oriented_span(span, degree, degrees))
            .collect::<Result<Vec<_>>>()?;
        let patch_kind = if group.spans.len() == 1 {
            "single_bezier"
        } else if group.id.starts_with("axial:") {
            "multi_span_nurbs"
        } else {
            "two_span_nurbs"
        };
        let mut patch = object([
            ("id", text(group.id.clone())),
            ("source_blocks", array(source_ids.clone())),
            ("cells", array(source_ids)),
            ("patch_kind", text(patch_kind)),
            ("degrees", indices(degree)),
            ("spans", array(spans)),
        ]);
        if let Json::Object(fields) = &mut patch {
            if group.spans.len() > 1 {
                fields.insert("exact_piece_extraction_expected".into(), Json::Bool(true));
            }
        }
        patch_values.push(patch);
    }
    if owner.len() != degrees.len() {
        return Err("patch layout does not cover every source block");
    }
    let mut external = Vec::new();
    let mut self_seams = Vec::new();
    for iface in model.field("topology")?.field("interfaces")?.items()? {
        let (a, _, _) = side(iface.field("first")?)?;
        let (b, _, _) = side(iface.field("second")?)?;
        let (patch_a, span_a) = *owner.get(&a).ok_or("missing source span")?;
        let (patch_b, span_b) = *owner.get(&b).ok_or("missing source span")?;
        let link = object([
            (
                "first",
                subface(span_a, patch_a, iface.field("first")?.items()?[1].text()?)?,
            ),
            (
                "second",
                subface(span_b, patch_b, iface.field("second")?.items()?[1].text()?)?,
            ),
            ("orientation", iface.field("orientation")?.clone()),
            ("status", text("exact_match")),
        ]);
        if patch_a == patch_b {
            let (id1, _, _) = side(iface.field("first")?)?;
            let (id2, _, _) = side(iface.field("second")?)?;
            let side1 = iface.field("first")?.items()?[1].text()?;
            let side2 = iface.field("second")?.items()?[1].text()?;
            if !same_grid_face(owner[&id1].1, side1, owner[&id2].1, side2)? {
                self_seams.push(link);
            }
        } else {
            external.push(link);
        }
    }
    let mut boundaries = Vec::new();
    for face in model.field("topology")?.field("exterior")?.items()? {
        let (id, _, _) = side(face)?;
        let span = owner.get(&id).ok_or("missing exterior source span")?;
        boundaries.push(subface(span.1, span.0, face.items()?[1].text()?)?);
    }
    let model_schema = model.field("schema")?.clone();
    let mut fields = vec![
        ("schema", text("gallery_exact_patch_map_v1")),
        ("method", text("deterministic_exact_pair_and_axial_chain_planning")),
        ("scope_note", text("Exact extracted cells assembled into disjoint C0 tensor-product NURBS covers. The planner is deterministic but does not claim a globally minimum patch count, global injectivity, or material coverage beyond the audited source model.")),
        ("source_schema", model_schema),
        ("source_block_count", num(degrees.len())),
        ("patches", array(patch_values)),
        ("interfaces", array(external)),
        ("boundary_faces", array(boundaries)),
        ("self_seams", array(self_seams)),
    ];
    for key in [
        "source_model",
        "source_recipe",
        "case_id",
        "source_model_sha256",
        "source_recipe_sha256",
    ] {
        if let Ok(value) = model.field(key) {
            fields.push((key, value.clone()));
        }
    }
    Ok(Json::Object(
        fields
            .into_iter()
            .map(|(key, value)| (key.to_owned(), value))
            .collect(),
    ))
}

fn same_grid_face(a: &SpanPlan, side_a: &str, b: &SpanPlan, side_b: &str) -> Result<bool> {
    let grid = |span: &SpanPlan, name: &str| -> Result<[[usize; 2]; 3]> {
        let (axis, source_endpoint) = side_name(name)?;
        let patch_axis = span
            .axis_map
            .iter()
            .position(|x| *x == axis)
            .ok_or("invalid map")?;
        let endpoint = source_endpoint ^ usize::from(span.flips[patch_axis]);
        let mut face = [[0, 1]; 3];
        face[0] = span.interval;
        face[patch_axis] = [span.interval[endpoint]; 2];
        Ok(face)
    };
    let x = grid(a, side_a)?;
    let y = grid(b, side_b)?;
    Ok(x == y)
}

/// Plan a deterministic cover and use the native exact assembler as its feasibility check.
pub fn plan(model: &Json) -> Result<Json> {
    crate::interchange::audit(model)?;
    if model.field("schema").ok().and_then(|x| x.text().ok())
        == Some("box_cylinder_boolean_suite_v1")
        && model.field("operation").ok().and_then(|x| x.text().ok()) == Some("cylinder")
    {
        // The tall cylinder has a separate periodic-cover contract. Reuse its
        // explicitly constructed two-patch cover rather than pretending that
        // a greedy pair cover solves the periodic closure problem.
        return crate::cylinder_tensor_covers::build_cover(model, 2)
            .map_err(|_| "custom tall-cylinder cover rejected");
    }
    let degrees = source_degrees(model)?;
    let mut groups = chain_groups(model, &degrees)?;
    let consumed: BTreeSet<_> = groups
        .iter()
        .flat_map(|g| g.spans.iter().map(|s| s.source.clone()))
        .collect();
    for id in degrees.keys().filter(|id| !consumed.contains(*id)) {
        groups.push(single_group(id));
    }
    groups.sort_by(|a, b| a.id.cmp(&b.id));
    // The source topology is the authority for the initial one-span cover.
    let baseline = cover_layout(model, &groups, &degrees)?;
    crate::nurbs_tensor_assembly::assemble_cover(model, &baseline)?;

    let interfaces = model.field("topology")?.field("interfaces")?.items()?;
    let mut counts = BTreeMap::<(String, String), usize>::new();
    for iface in interfaces {
        let (a, _, _) = side(iface.field("first")?)?;
        let (b, _, _) = side(iface.field("second")?)?;
        let key = if a <= b { (a, b) } else { (b, a) };
        *counts.entry(key).or_default() += 1;
    }
    let chain_members: BTreeSet<_> = groups
        .iter()
        .filter(|g| g.spans.len() > 1)
        .flat_map(|g| g.spans.iter().map(|s| s.source.clone()))
        .collect();
    for iface in interfaces {
        let (a, _, _) = side(iface.field("first")?)?;
        let (b, _, _) = side(iface.field("second")?)?;
        if a == b || chain_members.contains(&a) || chain_members.contains(&b) {
            continue;
        }
        let key = if a <= b {
            (a.clone(), b.clone())
        } else {
            (b.clone(), a.clone())
        };
        if counts.get(&key).copied().unwrap_or(0) != 1 {
            continue;
        }
        let Some(ga) = groups
            .iter()
            .position(|g| g.spans.len() == 1 && g.spans[0].source == a)
        else {
            continue;
        };
        let Some(gb) = groups
            .iter()
            .position(|g| g.spans.len() == 1 && g.spans[0].source == b)
        else {
            continue;
        };
        let mut accepted = None;
        'orient: for swap_a in [false, true] {
            for mask in 0..4 {
                let group = pair_group(&a, &b, iface, swap_a, [mask & 1 != 0, mask & 2 != 0]);
                let Some(group) = group else { continue };
                let mut candidate = groups.clone();
                let high = ga.max(gb);
                let low = ga.min(gb);
                candidate.remove(high);
                candidate.remove(low);
                candidate.push(group.clone());
                candidate.sort_by(|x, y| x.id.cmp(&y.id));
                if let Ok(layout) = cover_layout(model, &candidate, &degrees) {
                    if crate::nurbs_tensor_assembly::assemble_cover(model, &layout).is_ok() {
                        accepted = Some(candidate);
                        break 'orient;
                    }
                }
            }
        }
        if let Some(candidate) = accepted {
            groups = candidate;
        }
    }
    groups.sort_by(|a, b| a.id.cmp(&b.id));
    let layout = cover_layout(model, &groups, &degrees)?;
    crate::nurbs_tensor_assembly::assemble_cover(model, &layout)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(text: &str) -> Json {
        crate::json::parse(text).unwrap()
    }

    fn assert_complete_cover(model: &Json, cover: &Json) {
        assert_eq!(
            cover.field("source_block_count").unwrap().number::<usize>(),
            Ok(model.field("blocks").unwrap().items().unwrap().len())
        );
        let source: BTreeSet<_> = model
            .field("blocks")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .map(|b| b.field("id").unwrap().text().unwrap().to_owned())
            .collect();
        let mut covered = BTreeSet::new();
        for patch in cover.field("patches").unwrap().items().unwrap() {
            for cell in patch.field("cells").unwrap().items().unwrap() {
                assert!(covered.insert(cell.text().unwrap().to_owned()));
            }
        }
        assert_eq!(covered, source);
        assert_eq!(
            cover
                .field("boundary_faces")
                .unwrap()
                .items()
                .unwrap()
                .len(),
            model
                .field("topology")
                .unwrap()
                .field("exterior")
                .unwrap()
                .items()
                .unwrap()
                .len()
        );
    }

    #[test]
    fn plans_curved_gallery_cover_with_exact_native_assembly() {
        let model = read(include_str!(
            "../docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"
        ));
        let cover = plan(&model).unwrap();
        assert_complete_cover(&model, &cover);
        assert_eq!(cover.field("patches").unwrap().items().unwrap().len(), 5);
    }

    #[test]
    fn plans_tall_boolean_axial_chains_and_residual_pairs() {
        let model = read(include_str!(
            "../docs/site/box-cylinder/tall/union.model.json"
        ));
        let cover = plan(&model).unwrap();
        assert_complete_cover(&model, &cover);
        assert_eq!(cover.field("patches").unwrap().items().unwrap().len(), 7);
        assert_eq!(
            cover
                .field("patches")
                .unwrap()
                .items()
                .unwrap()
                .iter()
                .filter(|patch| patch.field("cells").unwrap().items().unwrap().len() == 3)
                .count(),
            5
        );
    }

    #[test]
    fn dispatches_periodic_cylinder_to_custom_cover() {
        let model = read(include_str!(
            "../docs/site/box-cylinder/tall/cylinder.model.json"
        ));
        let cover = plan(&model).unwrap();
        assert_complete_cover(&model, &cover);
        assert_eq!(cover.field("patch_count").unwrap().number::<usize>(), Ok(2));
        assert_eq!(cover.field("self_seams").unwrap().items().unwrap().len(), 3);
    }
}
