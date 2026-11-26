//! Exact assembly of certified trilinear sketch cells into straight-edged strips.

use std::collections::{BTreeMap, BTreeSet};

/// A strip represented by exact homogeneous section control points.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Patch {
    /// Four corners per station, ordered `[radial * 2 + z]`.
    pub sections: Vec<[[i64; 4]; 4]>,
    /// Indices into the input cell array, in strip order.
    pub cells: Vec<usize>,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Q(i128, i128);

impl Q {
    fn new(n: i128, d: i128) -> Self {
        assert_ne!(d, 0, "positive affine HEX8 precondition: nonzero weights");
        let (mut n, mut d) = if d < 0 { (-n, -d) } else { (n, d) };
        let g = gcd(n, d);
        n /= g;
        d /= g;
        Self(n, d)
    }
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.max(1)
}

type P = [Q; 3];
type S = [P; 4];

fn idx(i: usize, j: usize, k: usize) -> usize {
    4 * i + 2 * j + k
}

fn line_key(a: [i64; 4], b: [i64; 4]) -> Option<[i128; 3]> {
    let a = [a[0] as i128, a[1] as i128, a[3] as i128];
    let b = [b[0] as i128, b[1] as i128, b[3] as i128];
    let det = |x: i128, y: i128, z: i128, w: i128| x.checked_mul(y)?.checked_sub(z.checked_mul(w)?);
    let l = [
        det(a[1], b[2], a[2], b[1])?,
        det(a[2], b[0], a[0], b[2])?,
        det(a[0], b[1], a[1], b[0])?,
    ];
    if l == [0; 3] || l.contains(&i128::MIN) {
        return None;
    }
    let g = gcd(gcd(l[0], l[1]), l[2]);
    let mut l = [l[0] / g, l[1] / g, l[2] / g];
    if l.iter().find(|&&v| v != 0).is_some_and(|&v| v < 0) {
        l.iter_mut().for_each(|v| *v = -*v);
    }
    Some(l)
}

fn section_for_side(
    points: &[P; 8],
    along: usize,
    radial: usize,
    station: usize,
    side: usize,
    along_flip: bool,
) -> S {
    let mut out = [[Q(0, 1); 3]; 4];
    for r in 0..2 {
        for z in 0..2 {
            let mut a = [0; 3];
            a[along] = station ^ usize::from(along_flip);
            a[radial] = if r == 0 { side } else { side ^ 1 };
            a[2] = z;
            out[2 * r + z] = points[idx(a[0], a[1], a[2])];
        }
    }
    out
}

fn common_scale(sections: &[S]) -> Option<Vec<[[i64; 4]; 4]>> {
    let mut scale = 1_i128;
    for s in sections {
        for p in s {
            for q in p {
                scale = scale.checked_mul(q.1 / gcd(scale, q.1))?;
            }
        }
    }
    sections
        .iter()
        .map(|s| {
            let mut out = [[0_i64; 4]; 4];
            for (i, p) in s.iter().enumerate() {
                for axis in 0..3 {
                    out[i][axis] = i64::try_from(p[axis].0.checked_mul(scale / p[axis].1)?).ok()?;
                }
                out[i][3] = i64::try_from(scale).ok()?;
            }
            Some(out)
        })
        .collect()
}

#[derive(Clone)]
struct Span {
    id: usize,
    start: S,
    end: S,
}

/// Join exact cell sections into disjoint strips along shared straight stock edges.
///
/// Preconditions: each input is a positive affine HEX8, ordered `4*i+2*j+k`
/// (z fastest), with one positive, uniform homogeneous weight per cell.
/// Rational point equality is used, so uniformly rescaled cells interoperate.
pub fn assemble(input: &[[[i64; 4]; 8]]) -> Vec<Patch> {
    let mut groups: BTreeMap<[i128; 3], Vec<Span>> = BTreeMap::new();
    for (id, raw) in input.iter().enumerate() {
        let weight = raw[0][3];
        if weight <= 0 || raw.iter().any(|point| point[3] != weight) {
            continue;
        }
        let points = std::array::from_fn(|i| {
            std::array::from_fn(|axis| Q::new(i128::from(raw[i][axis]), i128::from(weight)))
        });
        // Four positive XY rotations. For each boundary side, orient the strip
        // tangent so (along, inward, +z) agrees with positive cell orientation.
        for (along, radial) in [(0, 1), (1, 0)] {
            for side in 0..2 {
                let raw0 = section_for_side(&points, along, radial, 0, side, false);
                let raw1 = section_for_side(&points, along, radial, 1, side, false);
                let fixed = side;
                let mut aidx = [0; 3];
                let mut bidx = [0; 3];
                aidx[radial] = fixed;
                bidx[along] = 1;
                bidx[radial] = fixed;
                let Some(line) =
                    line_key(raw[idx(aidx[0], aidx[1], 0)], raw[idx(bidx[0], bidx[1], 0)])
                else {
                    continue;
                };
                // This is the positive XY rotation for a positive input map.
                let reverse = (along == 1) ^ (side == 1);
                let start = if reverse { raw1 } else { raw0 };
                let end = if reverse { raw0 } else { raw1 };
                groups
                    .entry(line)
                    .or_default()
                    .push(Span { id, start, end });
            }
        }
    }
    let mut proposed: Vec<(Vec<usize>, Vec<S>)> = Vec::new();
    for spans in groups.values() {
        let mut from: BTreeMap<S, Vec<usize>> = BTreeMap::new();
        for (i, s) in spans.iter().enumerate() {
            from.entry(s.start).or_default().push(i);
        }
        let mut next = vec![None; spans.len()];
        let mut has_prev = vec![false; spans.len()];
        for (i, s) in spans.iter().enumerate() {
            if let Some(options) = from.get(&s.end) {
                if let Some(&j) = options.iter().find(|&&j| j != i && spans[j].id != s.id) {
                    next[i] = Some(j);
                    has_prev[j] = true;
                }
            }
        }
        let mut starts: Vec<_> = (0..spans.len()).filter(|&i| !has_prev[i]).collect();
        starts.sort_by_key(|&i| spans[i].id);
        let mut visited = BTreeSet::new();
        for start in starts.into_iter().chain(0..spans.len()) {
            if visited.contains(&start) {
                continue;
            }
            let mut ids = Vec::new();
            let mut seq = vec![spans[start].start];
            let mut at = start;
            loop {
                if !visited.insert(at) {
                    break;
                }
                ids.push(spans[at].id);
                seq.push(spans[at].end);
                match next[at] {
                    Some(j) if !visited.contains(&j) => at = j,
                    _ => break,
                }
            }
            if !ids.is_empty() && seq.first() != seq.last() {
                proposed.push((ids, seq));
            }
        }
    }
    proposed.sort_by(|a, b| b.0.len().cmp(&a.0.len()).then_with(|| a.0[0].cmp(&b.0[0])));
    let mut used = BTreeSet::new();
    let mut result = Vec::new();
    for (ids, seq) in proposed {
        if ids.iter().any(|i| used.contains(i)) {
            continue;
        }
        if let Some(sections) = common_scale(&seq) {
            used.extend(ids.iter().copied());
            result.push(Patch {
                sections,
                cells: ids,
            });
        }
    }
    // Source cells not selected for a multi-cell edge remain explicit singleton patches.
    for (id, source) in input.iter().enumerate() {
        if used.contains(&id) {
            continue;
        }
        // A singleton preserves the original homogeneous coefficients.
        let make = |station: usize| std::array::from_fn(|n| source[idx(station, n / 2, n % 2)]);
        used.insert(id);
        result.push(Patch {
            sections: vec![make(0), make(1)],
            cells: vec![id],
        });
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    fn hex(x0: i64, x1: i64, y0: i64, y1: i64, z0: i64, z1: i64, scale: i64) -> [[i64; 4]; 8] {
        let mut c = [[0; 4]; 8];
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    c[idx(i, j, k)] = [
                        [x0, x1][i] * scale,
                        [y0, y1][j] * scale,
                        [z0, z1][k] * scale,
                        scale,
                    ];
                }
            }
        }
        c
    }
    #[test]
    fn joins_exact_sections_across_different_scales() {
        let cells = [hex(0, 2, 0, 1, 0, 3, 1), hex(2, 5, 0, 1, 0, 3, 7)];
        let patches = assemble(&cells);
        let p = patches.iter().find(|p| p.cells.len() == 2).unwrap();
        assert_eq!(p.sections.len(), 3);
        assert_eq!(p.cells, vec![0, 1]);
        assert_eq!(p.sections[1][0], [2, 0, 0, 1]);
    }
    #[test]
    fn variable_weights_preserve_the_exact_singleton_net() {
        let mut cell = hex(0, 1, 0, 1, 0, 1, 1);
        cell[7][3] = 2;
        let patches = assemble(&[cell]);
        assert_eq!(patches[0].sections.concat(), cell);
    }

    #[test]
    fn nonmatching_cells_remain_singletons() {
        let cells = [hex(0, 1, 0, 1, 0, 1, 1), hex(2, 3, 0, 1, 0, 1, 1)];
        assert!(assemble(&cells).iter().all(|p| p.cells.len() == 1));
    }
}
