//! Exact removal of degree-one knots from tensor-product homogeneous control nets.
//!
//! Removal is accepted only when reinserting every removed knot reconstructs every
//! homogeneous component exactly using bounded i128 integer arithmetic.

use crate::json::{self, Json};
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq)]
pub struct ReducedLine {
    pub knots: Vec<i64>,
    pub controls: Vec<[i64; 4]>,
    pub removed_knots: Vec<i64>,
}

/// Remove every interior degree-one knot whose reinsertion is exact.
/// Nonremovable interior knots are retained; malformed vectors are rejected.
pub fn remove_degree_one_knots(
    knots: &[i64],
    controls: &[[i64; 4]],
) -> Result<ReducedLine, &'static str> {
    let keep = retained_planes(knots, controls, controls.len(), 1)?;
    let (knots, removed_knots) = filter_knots(knots, &keep);
    let controls = controls
        .iter()
        .zip(keep)
        .filter_map(|(p, k)| k.then_some(*p))
        .collect();
    Ok(ReducedLine {
        knots,
        controls,
        removed_knots,
    })
}

// The derivative of a linear spline is constant on each open interval. A plane
// can disappear exactly when the adjacent slopes agree on every tensor line.
// Equal slopes are transitive, so the mask needs just one pass over the input.
fn retained_planes(
    knots: &[i64],
    controls: &[[i64; 4]],
    count: usize,
    stride: usize,
) -> Result<Vec<bool>, &'static str> {
    if count < 2
        || count.checked_add(2) != Some(knots.len())
        || knots[0] != knots[1]
        || knots[count] != knots[count + 1]
        || knots[1..=count].windows(2).any(|w| w[0] >= w[1])
    {
        return Err("expected a clamped degree-one vector with distinct interior knots");
    }
    if controls.iter().any(|p| p[3] <= 0) {
        return Err("nonpositive homogeneous weight");
    }
    let mut keep = vec![false; count];
    keep[0] = true;
    keep[count - 1] = true;
    for (i, middle) in controls.iter().enumerate() {
        let plane = i / stride % count;
        if keep[plane] {
            continue;
        }
        let left = i128::from(knots[plane + 1]) - i128::from(knots[plane]);
        let right = i128::from(knots[plane + 2]) - i128::from(knots[plane + 1]);
        for c in 0..4 {
            // Both sides are bounded by the largest i64 coordinate times
            // the complete knot interval, which fits in signed i128.
            let before = i128::from(controls[i - stride][c]);
            let after = i128::from(controls[i + stride][c]);
            if right * before + left * after != (left + right) * i128::from(middle[c]) {
                keep[plane] = true;
                break;
            }
        }
    }
    Ok(keep)
}

fn filter_knots(knots: &[i64], keep: &[bool]) -> (Vec<i64>, Vec<i64>) {
    let mut kept = vec![knots[0]];
    let mut removed = Vec::new();
    for (&k, &retain) in knots[1..knots.len() - 1].iter().zip(keep) {
        if retain {
            kept.push(k);
        } else {
            removed.push(k);
        }
    }
    kept.push(knots[knots.len() - 1]);
    (kept, removed)
}

fn as_i64(v: &Json) -> Result<i64, &'static str> {
    v.number()
}
fn is_axis_permutation(map: &[i64]) -> bool {
    map.len() == 3
        && [0, 1, 2]
            .iter()
            .all(|axis| map.iter().filter(|&&x| x == *axis).count() == 1)
}
fn array_i64(v: &Json) -> Result<Vec<i64>, &'static str> {
    v.items()?.iter().map(as_i64).collect()
}
fn arr(values: impl IntoIterator<Item = Json>) -> Json {
    Json::Array(values.into_iter().collect())
}
fn num(n: i64) -> Json {
    Json::Number(n.to_string())
}
fn s(x: &str) -> Json {
    Json::String(x.to_owned())
}
fn obj(m: BTreeMap<String, Json>) -> Json {
    Json::Object(m)
}
fn field_clone(v: &Json, k: &str) -> Result<Json, &'static str> {
    Ok(v.field(k)?.clone())
}

/// Convert a published exact patch-map document to compact axial geometry.
/// The extrusion direction is the patch axis corresponding to source parameter w.
pub fn reduce_patch_document(input: &Json) -> Result<Json, &'static str> {
    let patches = input.field("patches")?.items()?;
    let mut compact = Vec::new();
    let mut original_knots = Vec::new();
    let mut provenance_patches = Vec::new();
    let mut reports = Vec::new();
    for patch in patches {
        let id = patch.field("id")?.text()?;
        let degrees = array_i64(patch.field("degrees")?)?;
        if degrees.len() != 3 {
            return Err("patch must have three parametric degrees");
        }
        let shape = array_i64(patch.field("control_shape")?)?;
        let knot_json = patch.field("knots")?.items()?;
        let knots: Vec<Vec<i64>> = knot_json.iter().map(array_i64).collect::<Result<_, _>>()?;
        if shape.len() != 3 || degrees.iter().any(|&d| !(1..=3).contains(&d)) || knots.len() != 3 {
            return Err("patch requires three dimensions and degrees in 1..=3");
        }
        let dims_checked: Vec<usize> = shape
            .iter()
            .map(|&x| usize::try_from(x).map_err(|_| "invalid control shape"))
            .collect::<Result<_, _>>()?;
        if dims_checked.contains(&0) {
            return Err("empty patch control dimension");
        }
        let spans = patch.field("spans")?.items()?;
        let first_map = array_i64(
            spans
                .first()
                .ok_or("patch has no source spans")?
                .field("axis_map")?,
        )?;
        if !is_axis_permutation(&first_map) {
            return Err("source axis map must be a permutation");
        }
        let axis = first_map
            .iter()
            .position(|&source| source == 2)
            .ok_or("source physical Z axis is unmapped")?;
        if degrees[axis] != 1 {
            return Err("physical Z patch direction must have degree one");
        }
        for span in spans {
            let map = array_i64(span.field("axis_map")?)?;
            if !is_axis_permutation(&map)
                || map.iter().position(|&source| source == 2) != Some(axis)
            {
                return Err("source Z axis mapping changes within patch spans");
            }
        }
        for a in 0..3 {
            let degree = degrees[a] as usize;
            if dims_checked[a] < degree + 1 {
                return Err("control count must be at least degree plus one");
            }
            if knots[a].windows(2).any(|w| w[0] > w[1])
                || usize::try_from(shape[a])
                    .ok()
                    .and_then(|n| n.checked_add(degrees[a] as usize + 1))
                    != Some(knots[a].len())
            {
                return Err("invalid knot ordering or knot/control count");
            }
            if knots[a][..=degree].iter().any(|&u| u != knots[a][0])
                || knots[a][knots[a].len() - degree - 1..]
                    .iter()
                    .any(|&u| u != *knots[a].last().unwrap())
            {
                return Err("knot vectors must be clamped");
            }
            if knots[a][degree] >= knots[a][knots[a].len() - degree - 1] {
                return Err("knot vector domain endpoints must be distinct");
            }
            let mut i = degree + 1;
            while i + degree + 1 < knots[a].len() {
                let value = knots[a][i];
                let mut multiplicity = 1;
                while i + multiplicity < knots[a].len() && knots[a][i + multiplicity] == value {
                    multiplicity += 1;
                }
                if multiplicity > degree {
                    return Err("interior knot multiplicity exceeds degree");
                }
                i += multiplicity;
            }
        }
        let raw_controls = patch.field("controls_homogeneous")?.items()?;
        let mut controls = Vec::new();
        for c in raw_controls {
            let a = array_i64(c)?;
            if a.len() != 4 {
                return Err("expected four homogeneous components");
            }
            if a[3] <= 0 {
                return Err("homogeneous weights must be positive");
            }
            controls.push([a[0], a[1], a[2], a[3]]);
        }
        let expected: usize = dims_checked.iter().try_fold(1usize, |n, &x| {
            n.checked_mul(x).ok_or("control shape overflow")
        })?;
        if controls.len() != expected {
            return Err("control shape does not match homogeneous controls");
        }
        let dims = dims_checked;
        let stride = dims[axis + 1..].iter().product::<usize>();
        let keep = retained_planes(&knots[axis], &controls, dims[axis], stride)?;
        let (reduced_knots, removed) = filter_knots(&knots[axis], &keep);
        let mut knot_out: Vec<Json> = knots
            .iter()
            .map(|k| arr(k.iter().copied().map(num)))
            .collect();
        knot_out[axis] = arr(reduced_knots.iter().copied().map(num));
        let mut out = BTreeMap::new();
        for key in [
            "id",
            "degrees",
            "control_shape",
            "homogeneous_common_scale_denominator",
            "patch_kind",
        ] {
            if let Ok(v) = patch.field(key) {
                out.insert(key.to_owned(), v.clone());
            }
        }
        let target_count = reduced_knots.len() - 2;
        let mut shape_new = dims.clone();
        shape_new[axis] = target_count;
        let shape_out: Vec<Json> = shape_new.iter().map(|&x| num(x as i64)).collect();
        let repacked = controls
            .iter()
            .enumerate()
            .filter_map(|(i, p)| keep[i / stride % dims[axis]].then_some(p));
        out.insert("control_shape".into(), arr(shape_out));
        out.insert("knots".into(), arr(knot_out));
        out.insert(
            "controls_homogeneous".into(),
            arr(repacked.map(|c| arr(c.iter().copied().map(num)))),
        );
        compact.push(obj(out));
        original_knots.push(json::object([
            ("patch_id", s(id)),
            ("knots", patch.field("knots")?.clone()),
        ]));
        let mut provenance = BTreeMap::new();
        if let Json::Object(m) = patch {
            for (k, v) in m {
                if k != "controls_homogeneous" {
                    provenance.insert(k.clone(), v.clone());
                }
            }
        }
        provenance_patches.push(obj(provenance));
        reports.push(json::object([
            ("patch_id", s(id)),
            ("physical_axis", num(axis as i64)),
            ("removed_knots", arr(removed.into_iter().map(num))),
        ]));
    }
    let mut metadata = BTreeMap::new();
    if let Json::Object(m) = input {
        for (k, v) in m {
            if k != "patches" {
                metadata.insert(k.clone(), v.clone());
            }
        }
    }
    let mut out = BTreeMap::new();
    out.insert("schema".into(), s("gallery_axial_knot_reduction_v1"));
    out.insert("source_schema".into(), field_clone(input, "schema")?);
    for key in [
        "source_model",
        "source_model_sha256",
        "source_recipe",
        "source_recipe_sha256",
    ] {
        if let Ok(v) = input.field(key) {
            out.insert(key.into(), v.clone());
        }
    }
    out.insert("analysis_knots".into(), arr(original_knots));
    out.insert("patches".into(), arr(compact));
    out.insert("metadata".into(), obj(metadata));
    out.insert(
        "provenance".into(),
        json::object([("patches", arr(provenance_patches))]),
    );
    out.insert("removals".into(), arr(reports));
    Ok(obj(out))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_uniform_and_nonuniform_removal() {
        let u = remove_degree_one_knots(
            &[0, 0, 1, 2, 3, 3],
            &[[0, 0, 0, 1], [1, 0, 0, 1], [2, 0, 0, 1], [3, 0, 0, 1]],
        )
        .unwrap();
        assert_eq!(u.knots, [0, 0, 3, 3]);
        assert_eq!(u.controls.len(), 2);
        let n = remove_degree_one_knots(
            &[0, 0, 1, 3, 5, 5],
            &[[0, 0, 0, 1], [1, 0, 0, 1], [3, 0, 0, 1], [5, 0, 0, 1]],
        )
        .unwrap();
        assert_eq!(n.knots, [0, 0, 5, 5]);
    }
    #[test]
    fn mixed_nonremovable_knot_is_retained() {
        let r = remove_degree_one_knots(
            &[0, 0, 1, 2, 3, 4, 4],
            &[
                [0, 0, 0, 1],
                [1, 0, 0, 1],
                [2, 0, 0, 1],
                [4, 0, 0, 1],
                [4, 0, 0, 1],
            ],
        )
        .unwrap();
        assert_eq!(r.knots, [0, 0, 2, 3, 4, 4]);
        assert_eq!(r.removed_knots, [1]);
    }
    #[test]
    fn rejects_malformed_and_homogeneous_mismatch() {
        assert!(remove_degree_one_knots(&[0, 0, 1, 2], &[[0, 0, 0, 1], [1, 0, 0, 1]]).is_err());
        let r = remove_degree_one_knots(
            &[0, 0, 1, 2, 2],
            &[[0, 0, 0, 1], [1, 0, 0, 1], [2, 0, 0, 2]],
        )
        .unwrap();
        assert_eq!(r.knots, [0, 0, 1, 2, 2]);
    }

    #[test]
    fn adapter_preserves_metadata_and_uses_mapped_physical_axis() {
        let mut net = Vec::new();
        for x in 0..4 {
            for y in 0..2 {
                for z in 0..2 {
                    net.push(format!("[{x},{y},{z},1]"));
                }
            }
        }
        let source_text = r#"{
          "schema":"gallery_exact_patch_map_v1","source_model":"m","source_model_sha256":"abc",
          "source_recipe":"r","source_recipe_sha256":"def","interfaces":[{"id":"i"}],"self_seams":[],"boundary_faces":[],
          "patches":[{"id":"p","degrees":[1,1,1],"control_shape":[4,2,2],
            "knots":[[0,0,1,2,3,3],[0,0,1,1],[0,0,1,1]],"controls_homogeneous":[__NET__],
            "spans":[{"axis_map":[2,0,1],"source_block":"b"}],"cells":["b"],"source_blocks":["b"]}]
        }"#.replace("__NET__",&net.join(","));
        let source = json::parse(&source_text).unwrap();
        let out = reduce_patch_document(&source).unwrap();
        assert_eq!(
            out.field("schema").unwrap().text().unwrap(),
            "gallery_axial_knot_reduction_v1"
        );
        let p = &out.field("patches").unwrap().items().unwrap()[0];
        assert_eq!(
            array_i64(p.field("control_shape").unwrap()).unwrap(),
            [2, 2, 2]
        );
        assert_eq!(
            array_i64(&p.field("knots").unwrap().items().unwrap()[0]).unwrap(),
            [0, 0, 3, 3]
        );
        assert_eq!(
            out.field("metadata").unwrap().field("interfaces").unwrap(),
            source.field("interfaces").unwrap()
        );
        assert_eq!(
            out.field("analysis_knots").unwrap().items().unwrap()[0]
                .field("knots")
                .unwrap(),
            source.field("patches").unwrap().items().unwrap()[0]
                .field("knots")
                .unwrap()
        );
        assert!(out
            .field("provenance")
            .unwrap()
            .field("patches")
            .unwrap()
            .items()
            .unwrap()[0]
            .field("controls_homogeneous")
            .is_err());
    }

    #[test]
    fn adapter_uses_global_plane_mask_when_one_line_is_not_affine() {
        let mut net = Vec::new();
        let second_line = [0, 0, 2, 4, 4];
        for (x, value) in second_line.into_iter().enumerate() {
            for y in 0..2 {
                for z in 0..2 {
                    let ordinate = if y == 0 { 0 } else { value };
                    net.push(format!("[{x},{ordinate},{z},1]"));
                }
            }
        }
        let source = json::parse(&format!(
            r#"{{
          "schema":"gallery_exact_patch_map_v1","interfaces":[],"self_seams":[],"boundary_faces":[],
          "patches":[{{"id":"p","degrees":[1,1,1],"control_shape":[5,2,2],
            "knots":[[0,0,1,2,3,4,4],[0,0,1,1],[0,0,1,1]],
            "controls_homogeneous":[{}],"spans":[{{"axis_map":[2,0,1]}}]}}]
        }}"#,
            net.join(",")
        ))
        .unwrap();
        let out = reduce_patch_document(&source).unwrap();
        let patch = &out.field("patches").unwrap().items().unwrap()[0];
        assert_eq!(
            array_i64(&patch.field("knots").unwrap().items().unwrap()[0]).unwrap(),
            [0, 0, 1, 3, 4, 4]
        );
        assert_eq!(
            out.field("removals").unwrap().items().unwrap()[0]
                .field("removed_knots")
                .unwrap()
                .items()
                .unwrap()
                .iter()
                .map(as_i64)
                .collect::<Result<Vec<_>, _>>()
                .unwrap(),
            [2]
        );
        assert_eq!(
            array_i64(patch.field("control_shape").unwrap()).unwrap(),
            [4, 2, 2]
        );
    }
}
