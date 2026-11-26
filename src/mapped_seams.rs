//! Discover hanging rectangular traces and pull their edge coordinates back
//! to the adjacent tensor chart. Other face shapes remain subject to the final
//! exact whole-face incidence and exterior-label checks.
use crate::bracket::BracketBlock;
use crate::constraint_system::Ratio;
use crate::incidence::{enumerate_faces, BlockSide};
use crate::mapped_cut_propagation::{CutParameter, CutRequest, FaceMap};
use std::collections::{BTreeMap, BTreeSet};

fn tangents(normal: usize) -> [usize; 2] {
    match normal {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}
struct Rectangle<'a> {
    side: &'a BlockSide,
    block: usize,
    normal: usize,
    sign: i8,
    plane: Ratio,
    world: [usize; 2],
    starts: [Ratio; 2],
    ends: [Ratio; 2],
}
impl Rectangle<'_> {
    fn interval(&self, axis: usize) -> [Ratio; 2] {
        let i = usize::from(self.world[1] == axis);
        let (a, b) = (self.starts[i], self.ends[i]);
        if a.cmp_value(&b).is_lt() {
            [a, b]
        } else {
            [b, a]
        }
    }
}
fn rectangle<'a>(
    b: &BracketBlock,
    side: &'a BlockSide,
    index: usize,
) -> Result<Option<Rectangle<'a>>, &'static str> {
    let [s, t] = tangents(side.axis);
    let get = |i: usize, j: usize| {
        let mut v = [0; 3];
        v[side.axis] = if side.upper { b.degrees[side.axis] } else { 0 };
        v[s] = i;
        v[t] = j;
        b.controls[(v[0] * (b.degrees[1] + 1) + v[1]) * (b.degrees[2] + 1) + v[2]]
    };
    let p = get(0, 0);
    let endpoints = [get(b.degrees[s], 0), get(0, b.degrees[t])];
    if p[3] <= 0 || endpoints.iter().any(|q| q[3] != p[3]) {
        return Ok(None);
    }
    let mut world = [0; 2];
    let mut inc = [0i128; 2];
    for a in 0..2 {
        let changed: Vec<_> = (0..3).filter(|&k| p[k] != endpoints[a][k]).collect();
        if changed.len() != 1 {
            return Ok(None);
        }
        world[a] = changed[0];
        inc[a] = i128::from(endpoints[a][world[a]]) - i128::from(p[world[a]]);
    }
    if world[0] == world[1] {
        return Ok(None);
    }
    let n = 3 - world[0] - world[1];
    for i in 0..=b.degrees[s] {
        for j in 0..=b.degrees[t] {
            let q = get(i, j);
            if q[3] != p[3] || q[n] != p[n] {
                return Ok(None);
            }
            for (a, k, deg) in [(0, i, b.degrees[s]), (1, j, b.degrees[t])] {
                let delta = i128::from(q[world[a]]) - i128::from(p[world[a]]);
                if deg as i128 * delta != k as i128 * inc[a] {
                    return Ok(None);
                }
            }
        }
    }
    let cyclic = (world[0] + 1) % 3 == world[1];
    let sign = if cyclic ^ (side.axis == 1) ^ !side.upper ^ (inc[0] < 0) ^ (inc[1] < 0) {
        1
    } else {
        -1
    };
    let frac = |v: i64| Ratio::new(v.into(), p[3].into());
    Ok(Some(Rectangle {
        side,
        block: index,
        normal: n,
        sign,
        plane: frac(p[n])?,
        world,
        starts: [frac(p[world[0]])?, frac(p[world[1]])?],
        ends: [frac(endpoints[0][world[0]])?, frac(endpoints[1][world[1]])?],
    }))
}
fn request(
    receiver: &Rectangle,
    other: &Rectangle,
    seen: &mut BTreeSet<(usize, usize, CutParameter, String)>,
) -> Result<(), &'static str> {
    for (i, world) in receiver.world.into_iter().enumerate() {
        let [lo, hi] = receiver.interval(world);
        for x in other.interval(world) {
            if !(lo.cmp_value(&x).is_lt() && x.cmp_value(&hi).is_lt()) {
                continue;
            }
            let length = receiver.ends[i].checked_sub(receiver.starts[i])?;
            let u = x
                .checked_sub(receiver.starts[i])?
                .checked_mul(Ratio::new(length.denominator(), length.numerator())?)?;
            let side = other.side;
            seen.insert((
                receiver.block,
                tangents(receiver.side.axis)[i],
                CutParameter::new(u.numerator(), u.denominator())?,
                format!(
                    "face/{}/{}/{}",
                    side.block,
                    side.axis,
                    usize::from(side.upper)
                ),
            ));
        }
    }
    Ok(())
}
pub(crate) fn discover(
    blocks: &[BracketBlock],
) -> Result<(Vec<FaceMap>, Vec<CutRequest>), &'static str> {
    let incidence = enumerate_faces(blocks)?;
    let index: BTreeMap<_, _> = blocks
        .iter()
        .enumerate()
        .map(|(i, b)| (b.id.as_str(), i))
        .collect();
    let mut maps = vec![];
    for interface in incidence.interfaces {
        let a = interface.first;
        let b = interface.second;
        let d = interface.orientation;
        let to = tangents(b.axis);
        maps.push(FaceMap {
            a_block: index[a.block.as_str()],
            a_normal_axis: a.axis,
            a_upper: a.upper,
            b_block: index[b.block.as_str()],
            b_normal_axis: b.axis,
            b_upper: b.upper,
            a_to_b: [
                (to[usize::from(d.swap)], d.flip_first),
                (to[usize::from(!d.swap)], d.flip_second),
            ],
            isoparametric_admitted: true,
            provenance: format!(
                "{}:{}:{}/{}:{}:{}",
                a.block, a.axis, a.upper, b.block, b.axis, b.upper
            ),
        });
    }
    let mut planes: BTreeMap<(usize, Ratio), Vec<Rectangle>> = BTreeMap::new();
    let mut seen = BTreeSet::new();
    for face in &incidence.exterior {
        let id = index[face.block.as_str()];
        if let Some(r) = rectangle(&blocks[id], face, id)? {
            let group = planes.entry((r.normal, r.plane)).or_default();
            for q in group.iter().filter(|q| q.sign != r.sign && q.block != id) {
                if r.world.iter().all(|&w| {
                    let a = r.interval(w);
                    let b = q.interval(w);
                    a[0].cmp_value(&b[1]).is_lt() && b[0].cmp_value(&a[1]).is_lt()
                }) {
                    request(&r, q, &mut seen)?;
                    request(q, &r, &mut seen)?;
                }
            }
            group.push(r);
        }
    }
    Ok((
        maps,
        seen.into_iter()
            .map(|(block, axis, parameter, provenance)| CutRequest {
                block,
                axis,
                parameter,
                provenance,
            })
            .collect(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compact_hosts::plain_box;
    use crate::mapped_cut_propagation::propagate_exact_cuts;

    #[test]
    fn hanging_thirds_close_through_a_previously_conforming_neighbor() {
        let blocks = vec![
            plain_box("back".into(), [[-2, 0], [0, 3], [0, 2]]),
            plain_box("front".into(), [[0, 2], [0, 3], [0, 2]]),
            plain_box("low".into(), [[2, 4], [0, 1], [0, 2]]),
            plain_box("high".into(), [[2, 4], [1, 3], [0, 2]]),
        ];
        let (faces, cuts) = discover(&blocks).unwrap();
        assert_eq!(cuts.len(), 2);
        assert!(cuts.iter().all(|s| s.block == 1
            && s.axis == 1
            && s.parameter == CutParameter::new(1, 3).unwrap()));
        let fine = propagate_exact_cuts(&blocks, &faces, &cuts, 6).unwrap();
        assert_eq!(fine.blocks.len(), 6);
        assert!(discover(&fine.blocks).unwrap().1.is_empty());
        assert_eq!(enumerate_faces(&fine.blocks).unwrap().interfaces.len(), 7);
    }
}
