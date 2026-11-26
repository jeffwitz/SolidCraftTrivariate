//! Exact closure of isoparametric cuts over oriented tensor-product interfaces.
//!
//! A cut is represented by a rational parameter on one block axis.  On every
//! incident complete face, its parameter is transported through the declared
//! axis permutation and reversal.  The module refuses a crossing that is not
//! declared isoparametric, and checks the full child count before allocation.

use crate::bracket::BracketBlock;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

type HomogeneousPoint = [i128; 4];
type SplitNet = (Vec<HomogeneousPoint>, Vec<HomogeneousPoint>);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CutParameter {
    pub numerator: i128,
    pub denominator: i128,
}

impl CutParameter {
    pub fn new(numerator: i128, denominator: i128) -> Result<Self, &'static str> {
        if denominator <= 0 || numerator <= 0 || numerator >= denominator {
            return Err("mapped_cut_parameter_outside_open_unit_interval");
        }
        let g = gcd(numerator as u128, denominator as u128) as i128;
        Ok(Self {
            numerator: numerator / g,
            denominator: denominator / g,
        })
    }

    fn reversed(self) -> Result<Self, &'static str> {
        Self::new(
            self.denominator
                .checked_sub(self.numerator)
                .ok_or("mapped_cut_integer_overflow")?,
            self.denominator,
        )
    }
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a.max(1)
}

/// A complete-face incidence. `a_to_b[t]` gives the destination block axis
/// corresponding to the `t`th ascending tangential axis of face A, and whether
/// that coordinate is reversed. The face normal axes must be omitted.
///
/// This topology-level primitive trusts `isoparametric_admitted`; production
/// callers should populate it from exact incidence and recertify the assembled
/// children and interfaces after closure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaceMap {
    pub a_block: usize,
    pub a_normal_axis: usize,
    pub a_upper: bool,
    pub b_block: usize,
    pub b_normal_axis: usize,
    pub b_upper: bool,
    pub a_to_b: [(usize, bool); 2],
    /// Set only after exact homogeneous traces and orientation were admitted.
    pub isoparametric_admitted: bool,
    pub provenance: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CutRequest {
    pub block: usize,
    pub axis: usize,
    pub parameter: CutParameter,
    pub provenance: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PropagatedCut {
    pub block: usize,
    pub axis: usize,
    pub parameter: CutParameter,
    pub provenance: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PropagatedPlan {
    pub blocks: Vec<BracketBlock>,
    pub cuts: Vec<PropagatedCut>,
    pub source_to_children: Vec<(String, Vec<String>)>,
    pub boundary_children: Vec<BoundaryChild>,
    pub requested_cut_count: usize,
    pub propagated_cut_count: usize,
    pub child_count: usize,
    pub budget: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundaryChild {
    pub source_block: usize,
    pub child_block: usize,
    pub axis: usize,
    pub upper: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Key {
    block: usize,
    axis: usize,
    parameter: CutParameter,
}

fn checked_mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("mapped_cut_integer_overflow")
}

fn blend(a: [i128; 4], b: [i128; 4], n: i128, d: i128) -> Result<[i128; 4], &'static str> {
    let mut out = [0; 4];
    for c in 0..4 {
        out[c] = checked_mul(d - n, a[c])?
            .checked_add(checked_mul(n, b[c])?)
            .ok_or("mapped_cut_integer_overflow")?;
    }
    Ok(out)
}

/// Split one homogeneous Bernstein net at `t`, preserving a common integer
/// scale on each child. This is exact homogeneous de Casteljau subdivision.
fn split_curve(net: &[[i64; 4]], t: CutParameter) -> Result<SplitNet, &'static str> {
    let degree = net.len().checked_sub(1).ok_or("mapped_cut_empty_curve")?;
    if degree == 0 || degree > 8 || net.iter().any(|p| p[3] <= 0) {
        return Err("mapped_cut_unsupported_curve_degree_or_weight");
    }
    let d = t.denominator;
    let p = t.numerator;
    let mut rows = vec![net.iter().map(|x| x.map(i128::from)).collect::<Vec<_>>()];
    for level in 1..=degree {
        let previous = &rows[level - 1];
        let mut row = Vec::with_capacity(degree + 1 - level);
        for i in 0..(degree + 1 - level) {
            row.push(blend(previous[i], previous[i + 1], p, d)?);
        }
        rows.push(row);
    }
    let mut left = Vec::with_capacity(degree + 1);
    let mut right = Vec::with_capacity(degree + 1);
    for k in 0..=degree {
        left.push(scale_point(rows[k][0], pow(d, degree - k)?)?);
        right.push(scale_point(rows[degree - k][k], pow(d, k)?)?);
    }
    Ok((left, right))
}

fn pow(x: i128, n: usize) -> Result<i128, &'static str> {
    let mut out = 1_i128;
    for _ in 0..n {
        out = checked_mul(out, x)?;
    }
    Ok(out)
}

fn scale_point(p: [i128; 4], s: i128) -> Result<[i128; 4], &'static str> {
    let mut q = [0; 4];
    for c in 0..4 {
        q[c] = checked_mul(p[c], s)?;
    }
    Ok(q)
}

fn to_i64_net(mut net: Vec<[i128; 4]>) -> Result<Vec<[i64; 4]>, &'static str> {
    let mut g = 0_u128;
    for p in &net {
        for x in p {
            g = gcd(g, x.unsigned_abs());
        }
    }
    let g = i128::try_from(g.max(1)).map_err(|_| "mapped_cut_integer_overflow")?;
    net.drain(..)
        .map(|p| {
            Ok([
                i64::try_from(p[0] / g).map_err(|_| "mapped_cut_child_control_exceeds_i64")?,
                i64::try_from(p[1] / g).map_err(|_| "mapped_cut_child_control_exceeds_i64")?,
                i64::try_from(p[2] / g).map_err(|_| "mapped_cut_child_control_exceeds_i64")?,
                i64::try_from(p[3] / g).map_err(|_| "mapped_cut_child_control_exceeds_i64")?,
            ])
        })
        .collect()
}

fn shape(block: &BracketBlock) -> Result<[usize; 3], &'static str> {
    if block.degrees.iter().any(|d| *d == 0 || *d > 8) {
        return Err("mapped_cut_invalid_block_shape");
    }
    let s = block
        .degrees
        .map(|d| d.checked_add(1).ok_or("mapped_cut_invalid_block_shape"));
    let s = [s[0]?, s[1]?, s[2]?];
    if block.controls.len() != s.iter().product::<usize>() {
        return Err("mapped_cut_invalid_block_shape");
    }
    if block.controls.iter().any(|p| p[3] <= 0) {
        return Err("mapped_cut_nonpositive_weight");
    }
    Ok(s)
}

fn idx(s: [usize; 3], p: [usize; 3]) -> usize {
    (p[0] * s[1] + p[1]) * s[2] + p[2]
}

fn split_block(
    block: &BracketBlock,
    axis: usize,
    t: CutParameter,
) -> Result<(BracketBlock, BracketBlock), &'static str> {
    let s = shape(block)?;
    if axis > 2 {
        return Err("mapped_cut_invalid_axis");
    }
    let mut left = vec![[0_i128; 4]; block.controls.len()];
    let mut right = vec![[0_i128; 4]; block.controls.len()];
    let other: Vec<usize> = (0..3).filter(|&a| a != axis).collect();
    for i in 0..s[other[0]] {
        for j in 0..s[other[1]] {
            let mut line = Vec::with_capacity(s[axis]);
            for k in 0..s[axis] {
                let mut p = [0; 3];
                p[other[0]] = i;
                p[other[1]] = j;
                p[axis] = k;
                line.push(block.controls[idx(s, p)]);
            }
            let (l, r) = split_curve(&line, t)?;
            for k in 0..s[axis] {
                let mut p = [0; 3];
                p[other[0]] = i;
                p[other[1]] = j;
                p[axis] = k;
                left[idx(s, p)] = l[k];
                right[idx(s, p)] = r[k];
            }
        }
    }
    let child = |suffix: &str, controls| BracketBlock {
        id: format!("{}{}", block.id, suffix),
        degrees: block.degrees,
        controls,
    };
    Ok((
        child("_cut_lo", to_i64_net(left)?),
        child("_cut_hi", to_i64_net(right)?),
    ))
}

fn sort_cuts(cuts: &mut Vec<CutParameter>) -> Result<(), &'static str> {
    cuts.sort_by(|a, b| {
        crate::constraint_system::Ratio::new(a.numerator, a.denominator)
            .expect("validated cut parameter")
            .cmp_value(
                &crate::constraint_system::Ratio::new(b.numerator, b.denominator)
                    .expect("validated cut parameter"),
            )
    });
    cuts.dedup();
    Ok(())
}

fn validate_map(m: &FaceMap, blocks: usize) -> Result<(), &'static str> {
    if m.a_block >= blocks
        || m.b_block >= blocks
        || m.a_block == m.b_block
        || m.a_normal_axis > 2
        || m.b_normal_axis > 2
    {
        return Err("mapped_cut_invalid_face_map");
    }
    let a_tan: Vec<_> = (0..3).filter(|&x| x != m.a_normal_axis).collect();
    let b_tan: BTreeSet<_> = (0..3).filter(|&x| x != m.b_normal_axis).collect();
    let mapped: BTreeSet<_> = m.a_to_b.iter().map(|x| x.0).collect();
    if mapped != b_tan
        || m.a_to_b
            .iter()
            .enumerate()
            .any(|(i, x)| x.0 == m.b_normal_axis || a_tan[i] == m.a_normal_axis)
    {
        return Err("mapped_cut_invalid_axis_permutation");
    }
    Ok(())
}

/// Close cuts over the supplied complete-face graph, then subdivide all
/// affected blocks. The caller supplies only faces whose traces were already
/// admitted by exact homogeneous comparison. Positive Jacobian certificates
/// transfer to each child by positive affine reparameterization.
pub fn propagate_exact_cuts(
    blocks: &[BracketBlock],
    faces: &[FaceMap],
    requests: &[CutRequest],
    budget: usize,
) -> Result<PropagatedPlan, &'static str> {
    if blocks.is_empty() || requests.is_empty() || budget == 0 {
        return Err("mapped_cut_empty_input_or_budget");
    }
    for m in faces {
        validate_map(m, blocks.len())?;
    }
    for block in blocks {
        shape(block)?;
    }
    let mut queue = VecDeque::new();
    let mut seen = BTreeSet::new();
    let mut provenance: BTreeMap<Key, BTreeSet<String>> = BTreeMap::new();
    for r in requests {
        if r.block >= blocks.len() || r.axis > 2 {
            return Err("mapped_cut_request_out_of_range");
        }
        shape(&blocks[r.block])?;
        if CutParameter::new(r.parameter.numerator, r.parameter.denominator)? != r.parameter {
            return Err("mapped_cut_noncanonical_parameter");
        }
        queue.push_back((
            Key {
                block: r.block,
                axis: r.axis,
                parameter: r.parameter,
            },
            r.provenance.clone(),
        ));
    }
    while let Some((key, source)) = queue.pop_front() {
        provenance.entry(key).or_default().insert(source.clone());
        if !seen.insert((key, source.clone())) {
            continue;
        }
        for face in faces {
            let (source_block, source_normal, destination_block, destination_normal, mapping) =
                if face.a_block == key.block {
                    (
                        face.a_block,
                        face.a_normal_axis,
                        face.b_block,
                        face.b_normal_axis,
                        face.a_to_b,
                    )
                } else if face.b_block == key.block {
                    let mut inverse = [(0, false); 2];
                    for (a_index, &(b_axis, reversed)) in face.a_to_b.iter().enumerate() {
                        let b_index = (0..3)
                            .filter(|&a| a != face.b_normal_axis)
                            .position(|a| a == b_axis)
                            .ok_or("mapped_cut_invalid_axis_permutation")?;
                        inverse[b_index] = (x_axis_for(face.a_normal_axis, a_index), reversed);
                    }
                    (
                        face.b_block,
                        face.b_normal_axis,
                        face.a_block,
                        face.a_normal_axis,
                        inverse,
                    )
                } else {
                    continue;
                };
            let _ = (source_block, destination_block);
            if key.axis == source_normal {
                continue;
            }
            if !face.isoparametric_admitted {
                return Err("mapped_cut_incompatible_trace");
            }
            let source_tangent_index = (0..3)
                .filter(|&a| a != source_normal)
                .position(|a| a == key.axis)
                .ok_or("mapped_cut_invalid_axis_permutation")?;
            let (dest_axis, reversed) = mapping[source_tangent_index];
            if dest_axis == destination_normal {
                return Err("mapped_cut_invalid_axis_permutation");
            }
            let parameter = if reversed {
                key.parameter.reversed()?
            } else {
                key.parameter
            };
            queue.push_back((
                Key {
                    block: destination_block,
                    axis: dest_axis,
                    parameter,
                },
                source.clone(),
            ));
        }
    }
    let mut per_block: Vec<[Vec<CutParameter>; 3]> = (0..blocks.len())
        .map(|_| std::array::from_fn(|_| Vec::new()))
        .collect();
    for key in provenance.keys() {
        per_block[key.block][key.axis].push(key.parameter);
    }
    let mut predicted = 0_usize;
    for cuts in &mut per_block {
        for axis_cuts in cuts.iter_mut() {
            sort_cuts(axis_cuts)?;
        }
        let n = cuts
            .iter()
            .try_fold(1_usize, |n, c| n.checked_mul(c.len() + 1))
            .ok_or("mapped_cut_block_budget")?;
        predicted = predicted.checked_add(n).ok_or("mapped_cut_block_budget")?;
    }
    if predicted > budget {
        return Err("mapped_cut_block_budget");
    }

    let propagated_cut_count = provenance.len();
    let mut output = Vec::with_capacity(predicted);
    let mut source_to_children = Vec::with_capacity(blocks.len());
    let mut boundary_children = Vec::new();
    for (source_index, block) in blocks.iter().enumerate() {
        let mut pieces = vec![(block.clone(), [[true, true]; 3])];
        for (axis, axis_cuts) in per_block[source_index].iter().enumerate() {
            let mut cuts = axis_cuts.clone();
            sort_cuts(&mut cuts)?;
            if cuts.is_empty() {
                continue;
            }
            let mut next = Vec::new();
            for (mut piece, mut piece_boundary) in pieces {
                let mut previous = CutParameter {
                    numerator: 0,
                    denominator: 1,
                };
                for cut in &cuts {
                    let local_num = cut
                        .numerator
                        .checked_mul(previous.denominator)
                        .and_then(|x| {
                            x.checked_sub(previous.numerator.checked_mul(cut.denominator)?)
                        })
                        .ok_or("mapped_cut_integer_overflow")?;
                    let local_den = cut
                        .denominator
                        .checked_mul(
                            previous
                                .denominator
                                .checked_sub(previous.numerator)
                                .ok_or("mapped_cut_integer_overflow")?,
                        )
                        .ok_or("mapped_cut_integer_overflow")?;
                    let local = CutParameter::new(local_num, local_den)?;
                    let prefix = piece.id.clone();
                    let (mut left, right) = split_block(&piece, axis, local)?;
                    left.id = format!(
                        "{}_a{}_{}of{}",
                        prefix, axis, cut.numerator, cut.denominator
                    );
                    let mut left_boundary = piece_boundary;
                    left_boundary[axis][1] = false;
                    next.push((left, left_boundary));
                    piece_boundary[axis][0] = false;
                    piece = right;
                    piece.id = format!(
                        "{}_a{}_after{}of{}",
                        prefix, axis, cut.numerator, cut.denominator
                    );
                    previous = *cut;
                }
                next.push((piece, piece_boundary));
            }
            pieces = next;
        }
        let ids = pieces.iter().map(|(b, _)| b.id.clone()).collect::<Vec<_>>();
        for (child, boundary) in pieces {
            let child_index = output.len();
            for (axis, sides) in boundary.iter().enumerate() {
                for upper in [false, true] {
                    if sides[usize::from(upper)] {
                        boundary_children.push(BoundaryChild {
                            source_block: source_index,
                            child_block: child_index,
                            axis,
                            upper,
                        });
                    }
                }
            }
            output.push(child);
        }
        source_to_children.push((block.id.clone(), ids));
    }
    let cuts = provenance
        .into_iter()
        .map(|(key, values)| PropagatedCut {
            block: key.block,
            axis: key.axis,
            parameter: key.parameter,
            provenance: values.into_iter().collect(),
        })
        .collect::<Vec<_>>();
    Ok(PropagatedPlan {
        blocks: output,
        cuts,
        source_to_children,
        boundary_children,
        requested_cut_count: requests.len(),
        propagated_cut_count,
        child_count: predicted,
        budget,
    })
}

fn x_axis_for(normal: usize, tangent_index: usize) -> usize {
    (0..3)
        .filter(|&a| a != normal)
        .nth(tangent_index)
        .unwrap_or(0)
}
