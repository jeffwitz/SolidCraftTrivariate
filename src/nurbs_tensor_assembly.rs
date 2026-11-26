//! Exact regular tensor assembly of rational Bézier cells on a supplied C0 knot grid.
use crate::bracket::BracketBlock;
use crate::constraint_system::Ratio;
use crate::json::Json;
use std::collections::{BTreeMap, BTreeSet};

type Result<T> = std::result::Result<T, &'static str>;
pub type Homogeneous = [Ratio; 4];

/// Axes and reversals map patch coordinates to the source cell's coordinates.
#[derive(Clone)]
pub struct Placement<'a> {
    pub source: &'a BracketBlock,
    pub cell: [usize; 3],
    pub axes: [usize; 3],
    pub flips: [bool; 3],
}

pub struct C0Tensor {
    degrees: [usize; 3],
    span_counts: [usize; 3],
    controls: Vec<Homogeneous>,
}

fn index(shape: [usize; 3], [i, j, k]: [usize; 3]) -> usize {
    (i * shape[1] + j) * shape[2] + k
}

fn coordinates([a, b, c]: [usize; 3]) -> impl Iterator<Item = [usize; 3]> {
    (0..a).flat_map(move |i| (0..b).flat_map(move |j| (0..c).map(move |k| [i, j, k])))
}

fn shape_for(p: [usize; 3], counts: [usize; 3]) -> Result<[usize; 3]> {
    let mut shape = [0; 3];
    for i in 0..3 {
        if !(1..=3).contains(&p[i]) || counts[i] == 0 {
            return Err("invalid degree or span count");
        }
        shape[i] = p[i]
            .checked_mul(counts[i])
            .and_then(|n| n.checked_add(1))
            .ok_or("tensor shape overflow")?;
    }
    Ok(shape)
}

fn orient_and_elevate(cell: &Placement, target: [usize; 3]) -> Result<Vec<Homogeneous>> {
    let source = cell.source;
    let mut sorted = cell.axes;
    sorted.sort();
    if sorted != [0, 1, 2] {
        return Err("invalid axis permutation");
    }
    let parity = (0..3)
        .flat_map(|i| (i + 1..3).map(move |j| (i, j)))
        .filter(|&(i, j)| cell.axes[i] > cell.axes[j])
        .count()
        + cell.flips.iter().filter(|&&v| v).count();
    if parity % 2 != 0 {
        return Err("orientation reversal");
    }
    let src_shape = shape_for(source.degrees, [1; 3])?;
    if src_shape.iter().product::<usize>() != source.controls.len() {
        return Err("invalid source controls");
    }
    let mut shape = cell.axes.map(|a| src_shape[a]);
    let mut net = Vec::with_capacity(source.controls.len());
    for ijk in coordinates(shape) {
        let mut src = [0; 3];
        for axis in 0..3 {
            src[cell.axes[axis]] = if cell.flips[axis] {
                shape[axis] - 1 - ijk[axis]
            } else {
                ijk[axis]
            };
        }
        let point = source.controls[index(src_shape, src)];
        if point[3] <= 0 {
            return Err("nonpositive source weight");
        }
        net.push(point.map(|x| Ratio::new(x.into(), 1).expect("i64 integer")));
    }
    for axis in 0..3 {
        if shape[axis] > target[axis] + 1 {
            return Err("degree reduction is not an elevation");
        }
        while shape[axis] < target[axis] + 1 {
            let old = shape;
            let new_degree = old[axis];
            shape[axis] += 1;
            let mut raised = Vec::with_capacity(shape.iter().product());
            for at in coordinates(shape) {
                let j = at[axis];
                let mut pred = at;
                pred[axis] = j.saturating_sub(1);
                let mut next = at;
                next[axis] = j.min(new_degree - 1);
                let x = net[index(old, pred)];
                let y = net[index(old, next)];
                let a = Ratio::new(j as i128, new_degree as i128)?;
                let b = Ratio::one().checked_sub(a)?;
                let mut point = [Ratio::zero(); 4];
                for d in 0..4 {
                    point[d] = a.checked_mul(x[d])?.checked_add(b.checked_mul(y[d])?)?;
                }
                raised.push(point);
            }
            net = raised;
        }
    }
    Ok(net)
}

impl C0Tensor {
    /// Assemble oriented cells. Shared controls must agree up to one positive scale per span.
    pub fn assemble(
        degrees: [usize; 3],
        counts: [usize; 3],
        mut cells: Vec<Placement>,
    ) -> Result<Self> {
        let shape = shape_for(degrees, counts)?;
        if counts.iter().try_fold(1usize, |n, &x| n.checked_mul(x)) != Some(cells.len()) {
            return Err("incomplete tensor grid");
        }
        cells.sort_by_key(|c| c.cell);
        let size = shape
            .iter()
            .try_fold(1usize, |n, &x| n.checked_mul(x))
            .ok_or("tensor shape overflow")?;
        let mut controls: Vec<Option<Homogeneous>> = vec![None; size];
        for (expected, placement) in coordinates(counts).zip(&cells) {
            if placement.cell != expected {
                return Err("duplicate or missing tensor cell");
            }
            let local = orient_and_elevate(placement, degrees)?;
            let slots: Vec<_> = coordinates(degrees.map(|d| d + 1))
                .map(|at| {
                    index(
                        shape,
                        std::array::from_fn(|a| degrees[a] * expected[a] + at[a]),
                    )
                })
                .collect();
            let mut scale = Ratio::one();
            for (&slot, point) in slots.iter().zip(&local) {
                if let Some(existing) = controls[slot] {
                    scale = existing[3]
                        .checked_mul(Ratio::new(point[3].denominator(), point[3].numerator())?)?;
                    break;
                }
            }
            for (slot, mut point) in slots.into_iter().zip(local) {
                for x in &mut point {
                    *x = x.checked_mul(scale)?;
                }
                if controls[slot].is_some_and(|old| old != point) {
                    return Err("incompatible homogeneous trace");
                }
                controls[slot] = Some(point);
            }
        }
        Ok(Self {
            degrees,
            span_counts: counts,
            controls: controls
                .into_iter()
                .collect::<Option<_>>()
                .ok_or("uncovered control")?,
        })
    }

    /// Exact Bézier extraction on a full-multiplicity knot grid, in patch orientation.
    pub fn extract(&self, cell: [usize; 3]) -> Result<Vec<Homogeneous>> {
        if (0..3).any(|i| cell[i] >= self.span_counts[i]) {
            return Err("cell outside knot domain");
        }
        let shape = shape_for(self.degrees, self.span_counts)?;
        Ok(coordinates(self.degrees.map(|d| d + 1))
            .map(|at| {
                self.controls[index(
                    shape,
                    std::array::from_fn(|i| self.degrees[i] * cell[i] + at[i]),
                )]
            })
            .collect())
    }
}

fn triple(j: &Json) -> Result<[usize; 3]> {
    j.items()?
        .iter()
        .map(Json::number)
        .collect::<Result<Vec<_>>>()?
        .try_into()
        .map_err(|_| "expected three indices")
}

fn integer_controls(net: &[Homogeneous]) -> Result<(i128, Json)> {
    let mut denominator = 1i128;
    for v in net.iter().flatten() {
        let (mut a, mut b) = (denominator, v.denominator());
        while b != 0 {
            (a, b) = (b, a % b);
        }
        denominator = (denominator / a)
            .checked_mul(v.denominator())
            .ok_or("common denominator overflow")?;
    }
    let mut controls = Vec::with_capacity(net.len());
    for point in net {
        let mut row = Vec::with_capacity(4);
        for x in point {
            let v = x
                .numerator()
                .checked_mul(denominator / x.denominator())
                .ok_or("scaled control overflow")?;
            row.push(Json::Number(v.to_string()));
        }
        controls.push(Json::Array(row));
    }
    Ok((denominator, Json::Array(controls)))
}

// A face has one constant coordinate and two unit intervals on the patch grid.
fn grid_face(span: &Json, side: &str) -> Result<[[usize; 2]; 3]> {
    let [axis_char, end] = side.as_bytes() else {
        return Err("invalid face side");
    };
    let src_axis = b"uvw"
        .iter()
        .position(|c| c == axis_char)
        .ok_or("invalid face axis")?;
    let axis = triple(span.field("axis_map")?)?
        .iter()
        .position(|&v| v == src_axis)
        .ok_or("missing axis")?;
    let flag = span
        .field("flips")?
        .items()?
        .get(axis)
        .ok_or("missing flip")?;
    let (b'0' | b'1', Json::Bool(flip)) = (end, flag) else {
        return Err("invalid endpoint or flip");
    };
    let interval = |r: &Json| -> Result<[usize; 2]> {
        let [a, b] = r.items()? else {
            return Err("invalid interval");
        };
        Ok([a.number()?, b.number()?])
    };
    let mut grid = if let Ok(v) = span.field("parameter_intervals") {
        let [i, j, k] = v.items()? else {
            return Err("invalid span location");
        };
        [interval(i)?, interval(j)?, interval(k)?]
    } else {
        [interval(span.field("parameter_interval")?)?, [0, 1], [0, 1]]
    };
    grid[axis] = [grid[axis][usize::from((*end == b'1') ^ flip)]; 2];
    Ok(grid)
}

/// Compare source faces and orientations with the subface records carried by a layout.
fn validate_topology(model: &Json, layout: &Json) -> Result<()> {
    if layout.field("schema")?.text()? != "gallery_exact_patch_map_v1"
        || layout.field("source_schema")? != model.field("schema")?
    {
        return Err("patch layout schema mismatch");
    }
    let mut owners = BTreeMap::new();
    let mut patch_ids = BTreeSet::new();
    for patch in layout.field("patches")?.items()? {
        let owner = patch.field("id")?.text()?;
        if !patch_ids.insert(owner) {
            return Err("duplicate patch id");
        }
        let mut ids = BTreeSet::new();
        for span in patch.field("spans")?.items()? {
            let id = span.field("source_block")?.text()?;
            if owners.insert(id, (owner, span)).is_some() {
                return Err("duplicate source block");
            }
            ids.insert(id);
        }
        for field in ["source_blocks", "cells"] {
            let entries = patch.field(field)?.items()?;
            if entries.len() != ids.len()
                || entries
                    .iter()
                    .map(Json::text)
                    .collect::<Result<BTreeSet<_>>>()?
                    != ids
            {
                return Err("patch sources disagree with its spans");
            }
        }
    }
    let key = |record: &Json| -> Result<(String, String)> {
        let id = record.field("source_block")?.text()?;
        let side = record.field("side")?.text()?;
        let &(owner, span) = owners.get(id).ok_or("unknown subface block")?;
        let [axis_char, end] = side.as_bytes() else {
            return Err("invalid subface side");
        };
        let source_axis = b"uvw"
            .iter()
            .position(|x| x == axis_char)
            .ok_or("invalid subface axis")?;
        let axis = triple(span.field("axis_map")?)?
            .iter()
            .position(|&a| a == source_axis)
            .ok_or("missing subface axis")?;
        let flip = span
            .field("flips")?
            .items()?
            .get(axis)
            .ok_or("missing subface flip")?;
        let mapped_end = match (end, flip) {
            (b'0', Json::Bool(false)) | (b'1', Json::Bool(true)) => 0,
            (b'1', Json::Bool(false)) | (b'0', Json::Bool(true)) => 1,
            _ => return Err("invalid subface endpoint"),
        };
        if owner != record.field("patch_id")?.text()?
            || record.field("patch_side")?.text()?
                != format!("{}{mapped_end}", b"uvw"[axis] as char)
        {
            return Err("incorrect subface owner or side");
        }
        for field in [
            "flips",
            "axis_map",
            "parameter_intervals",
            "parameter_interval",
        ] {
            if let Ok(value) = record.field(field) {
                if span.field(field)? != value {
                    return Err("incorrect subface parameter mapping");
                }
            }
        }
        if record.field("parameter_interval").is_err()
            && record.field("parameter_intervals").is_err()
        {
            return Err("missing subface location");
        }
        Ok((id.to_owned(), side.to_owned()))
    };
    let source_key = |face: &Json| -> Result<(String, String)> {
        let [id, side] = face.items()? else {
            return Err("invalid source face");
        };
        Ok((id.text()?.to_owned(), side.text()?.to_owned()))
    };
    let topology = model.field("topology")?;
    let mut exterior = topology
        .field("exterior")?
        .items()?
        .iter()
        .map(source_key)
        .collect::<Result<BTreeSet<_>>>()?;
    for face in layout.field("boundary_faces")?.items()? {
        if !exterior.remove(&key(face)?) {
            return Err("unknown or duplicate exterior subface");
        }
    }
    if !exterior.is_empty() {
        return Err("missing exterior subface");
    }
    let mut external = BTreeSet::new();
    let mut internal = BTreeSet::new();
    for link in topology.field("interfaces")?.items()? {
        let a = source_key(link.field("first")?)?;
        let b = source_key(link.field("second")?)?;
        let &(owner_a, span_a) = owners.get(a.0.as_str()).ok_or("missing source span")?;
        let &(owner_b, span_b) = owners.get(b.0.as_str()).ok_or("missing source span")?;
        let target = if owner_a == owner_b {
            if grid_face(span_a, &a.1)? == grid_face(span_b, &b.1)? {
                continue;
            }
            &mut internal
        } else {
            &mut external
        };
        target.insert((a, b, link.field("orientation")?.to_string()));
    }
    for (field, remaining) in [("interfaces", &mut external), ("self_seams", &mut internal)] {
        if field == "self_seams" && layout.field(field).is_err() {
            continue;
        }
        for link in layout.field(field)?.items()? {
            let pair = (
                key(link.field("first")?)?,
                key(link.field("second")?)?,
                link.field("orientation")?.to_string(),
            );
            if !remaining.remove(&pair) {
                return Err("invalid or repeated interface");
            }
        }
    }
    if !internal.is_empty() {
        return Err("missing periodic seam");
    }
    if !external.is_empty() {
        return Err("missing interpatch interface");
    }
    Ok(())
}

/// Rebuild every patch net from source geometry and orientation/placement data.
/// Knots express a unit C0 tensor grid. Input patch control nets are never read.
pub fn assemble_cover(model: &Json, layout: &Json) -> Result<Json> {
    validate_topology(model, layout)?;
    let blocks = crate::bracket_json::decode_blocks(model)?;
    let sources: BTreeMap<_, _> = blocks.iter().map(|b| (b.id.as_str(), b)).collect();
    let mut used = BTreeSet::new();
    let mut output = layout.clone();
    let Json::Object(root) = &mut output else {
        return Err("expected layout object");
    };
    let Some(Json::Array(patches)) = root.get_mut("patches") else {
        return Err("missing patches");
    };
    for patch in patches {
        let degrees = triple(patch.field("degrees")?)?;
        let mut cells = Vec::new();
        let mut counts = [0; 3];
        for span in patch.field("spans")?.items()? {
            let id = span.field("source_block")?.text()?;
            if !used.insert(id.to_owned()) {
                return Err("source cell used twice");
            }
            let source = *sources.get(id).ok_or("unknown source cell")?;
            let pair = |v: &Json| -> Result<usize> {
                let ab = v.items()?;
                if ab.len() != 2 {
                    return Err("invalid parameter interval");
                }
                let a = ab[0].number::<usize>()?;
                if a.checked_add(1) != Some(ab[1].number()?) {
                    return Err("expected unit knot interval");
                }
                Ok(a)
            };
            let cell = if let Ok(value) = span.field("parameter_intervals") {
                value
                    .items()?
                    .iter()
                    .map(pair)
                    .collect::<Result<Vec<_>>>()?
                    .try_into()
                    .map_err(|_| "expected three intervals")?
            } else {
                [pair(span.field("parameter_interval")?)?, 0, 0]
            };
            for i in 0..3 {
                counts[i] = counts[i].max(cell[i] + 1);
            }
            let flips = span
                .field("flips")?
                .items()?
                .iter()
                .map(|f| match f {
                    Json::Bool(v) => Ok(*v),
                    _ => Err("expected flip flag"),
                })
                .collect::<Result<Vec<_>>>()?
                .try_into()
                .map_err(|_| "expected three flips")?;
            cells.push(Placement {
                source,
                axes: triple(span.field("axis_map")?)?,
                flips,
                cell,
            });
        }
        let net = C0Tensor::assemble(degrees, counts, cells)?;
        let (factor, controls) = integer_controls(&net.controls)?;
        let Json::Object(fields) = patch else {
            return Err("invalid patch");
        };
        fields.insert("controls_homogeneous".into(), controls);
        fields.insert(
            "homogeneous_common_scale_denominator".into(),
            Json::Number(factor.to_string()),
        );
        let nums = |v: Vec<usize>| {
            Json::Array(v.into_iter().map(|x| Json::Number(x.to_string())).collect())
        };
        fields.insert(
            "control_shape".into(),
            nums(shape_for(degrees, counts)?.to_vec()),
        );
        let mut knots = Vec::new();
        for (d, n) in degrees.into_iter().zip(counts) {
            let mut k = vec![0; d + 1];
            for j in 1..n {
                k.extend(vec![j; d]);
            }
            k.extend(vec![n; d + 1]);
            knots.push(nums(k));
        }
        fields.insert("knots".into(), Json::Array(knots));
    }
    if used.len() != blocks.len() {
        return Err("unassigned source cells");
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn slab(x: i64, degree: usize, weight: i64) -> BracketBlock {
        let mut controls = Vec::new();
        for [i, j, k] in coordinates([2, degree + 1, 2]) {
            controls.push([
                (x + i as i64) * weight,
                j as i64 * 2 * weight / degree as i64,
                k as i64 * weight,
                weight,
            ]);
        }
        BracketBlock {
            id: format!("slab_{x}"),
            degrees: [1, degree, 1],
            controls,
        }
    }
    #[test]
    fn scaled_and_elevated_join_extracts_on_both_sides() {
        let left = slab(0, 1, 1);
        let right = slab(1, 2, 7);
        let cell = |source, pos| Placement {
            source,
            cell: [pos, 0, 0],
            axes: [0, 1, 2],
            flips: [false; 3],
        };
        let net = C0Tensor::assemble([1, 2, 1], [2, 1, 1], vec![cell(&right, 1), cell(&left, 0)])
            .unwrap();
        for n in 0..2 {
            let actual = net.extract([n, 0, 0]).unwrap();
            for (p, [i, j, k]) in actual.iter().zip(coordinates([2, 3, 2])) {
                assert_eq!(
                    *p,
                    [(n + i) as i128, j as i128, k as i128, 1].map(|v| Ratio::new(v, 1).unwrap())
                );
            }
        }
        assert!(net.extract([2, 0, 0]).is_err());
        let mut damaged = right.clone();
        damaged.controls[2][0] += 1;
        assert!(C0Tensor::assemble(
            [1, 2, 1],
            [2, 1, 1],
            vec![cell(&left, 0), cell(&damaged, 1)]
        )
        .is_err());
    }
}
