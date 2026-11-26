//! Optional Steiner vertices on host sides. A dynamic program traverses the
//! host and hole cycles together. Adjacent hosts reuse the first side trace.
use crate::{constraint_system::Ratio, sketch_convex_rings::quality};
use std::collections::{BTreeMap, BTreeSet};
type Point = [Ratio; 2];
type Quad = [Point; 4];
type Edge = (Point, Point);
fn key(a: Point, b: Point) -> Edge {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}
fn xy(p: Point) -> [f64; 2] {
    p.map(|v| v.numerator() as f64 / v.denominator() as f64)
}

fn proposals(a: Point, b: Point, hole: &[Point]) -> Option<Vec<Point>> {
    let [x, y] = [a, b].map(xy);
    let d = [y[0] - x[0], y[1] - x[1]];
    let length = d[0] * d[0] + d[1] * d[1];
    let mut slots = BTreeSet::from([4, 8, 12]);
    for &h in hole {
        let v = xy(h);
        let t = ((v[0] - x[0]) * d[0] + (v[1] - x[1]) * d[1]) / length;
        let s = (16. * t).round() as i32;
        if (1..16).contains(&s) {
            slots.insert(s);
        }
    }
    slots
        .into_iter()
        .map(|k| {
            let t = Ratio::new(i128::from(k), 16).ok()?;
            let coord = |j: usize| {
                b[j].checked_sub(a[j])
                    .ok()?
                    .checked_mul(t)
                    .ok()?
                    .checked_add(a[j])
                    .ok()
            };
            Some([coord(0)?, coord(1)?])
        })
        .collect()
}

#[derive(Clone, Copy)]
struct Step {
    end: usize,
    advance: usize,
    quad: Quad,
    condition: f64,
}
#[derive(Clone, Copy)]
struct State {
    count: usize,
    worst: f64,
    prev: Option<(usize, usize, usize)>,
}

fn walk(points: &[Point], required: &[bool], hole: &[Point]) -> Option<Vec<Quad>> {
    let (m, n) = (points.len(), hole.len());
    if m * n > 4096 || m * m * n > 250_000 {
        return None;
    }
    let mut steps = vec![Vec::new(); m * n];
    let next = |i: usize| {
        let mut ends = Vec::new();
        for k in i + 1..=i + m {
            ends.push(k);
            if required[k % m] {
                break;
            }
        }
        ends
    };
    for i in 0..m {
        for j in 0..n {
            let mut add = |end: usize, advance: usize, quad: Quad| {
                let condition = quality(quad);
                if condition.is_finite() {
                    steps[i * n + j].push(Step {
                        end: end - i,
                        advance,
                        quad,
                        condition,
                    })
                }
            };
            for k in next(i) {
                add(k, 1, [points[i], points[k % m], hole[(j + 1) % n], hole[j]]);
                if required[k % m] {
                    for t in next(k) {
                        add(t, 0, [points[i], points[k % m], points[t % m], hole[j]])
                    }
                }
            }
        }
    }
    // Fix the first hole vertex and try its spoke at each host vertex.
    // The directed cycles close only after both boundaries have been consumed.
    let mut best: Option<(usize, f64, Vec<Quad>)> = None;
    for start in 0..m {
        let mut states = vec![None; (m + 1) * (n + 1)];
        states[0] = Some(State {
            count: 0,
            worst: 0.,
            prev: None,
        });
        for i in 0..m {
            for j in 0..=n {
                let idx = i * (n + 1) + j;
                let Some(state) = states[idx] else { continue };
                for (s, t) in steps[((start + i) % m) * n + j % n].iter().enumerate() {
                    let end = i + t.end;
                    let h = j + t.advance;
                    if end > m || h > n {
                        continue;
                    }
                    let count = state.count + 1;
                    let worst = state.worst.max(t.condition);
                    let dest = &mut states[end * (n + 1) + h];
                    if dest.is_none_or(|old| {
                        count < old.count || count == old.count && worst < old.worst
                    }) {
                        *dest = Some(State {
                            count,
                            worst,
                            prev: Some((i, j, s)),
                        });
                    }
                }
            }
        }
        if let Some(last) = states[m * (n + 1) + n] {
            if best
                .as_ref()
                .is_some_and(|(c, w, _)| (*c, *w) <= (last.count, last.worst))
            {
                continue;
            }
            let (mut i, mut j) = (m, n);
            let mut path = Vec::new();
            while let Some((p, h, s)) = states[i * (n + 1) + j]?.prev {
                path.push(steps[((start + p) % m) * n + h % n][s].quad);
                (i, j) = (p, h);
            }
            path.reverse();
            best = Some((last.count, last.worst, path));
        }
    }
    best.map(|(_, _, q)| q)
}

pub(super) fn rings(
    hosts: &[Vec<Point>],
    holes: &[Vec<Point>],
    order: &[usize],
) -> Option<Vec<Quad>> {
    let mut traces = BTreeMap::<Edge, Vec<Point>>::new();
    let mut cells = Vec::new();
    for &id in order {
        let host = &hosts[id];
        let mut points = Vec::new();
        let mut required = Vec::new();
        let mut sides = Vec::new();
        for i in 0..host.len() {
            let a = host[i];
            let b = host[(i + 1) % host.len()];
            let edge = key(a, b);
            points.push(a);
            required.push(true);
            let previous = traces.get(&edge);
            let mut extra = match previous {
                Some(trace) => trace.clone(),
                None => proposals(a, b, &holes[id])?,
            };
            if previous.is_some() && a != edge.0 {
                extra.reverse()
            }
            let start = points.len();
            required.extend(std::iter::repeat_n(previous.is_some(), extra.len()));
            points.extend(extra);
            sides.push((edge, a, start..points.len()));
        }
        let ring = walk(&points, &required, &holes[id])?;
        let used: BTreeSet<_> = ring.iter().flatten().copied().collect();
        for (edge, a, range) in sides {
            traces.entry(edge).or_insert_with(|| {
                let mut selected: Vec<_> = points[range]
                    .iter()
                    .copied()
                    .filter(|p| used.contains(p))
                    .collect();
                if edge.0 != a {
                    selected.reverse()
                }
                selected
            });
        }
        cells.extend(ring);
    }
    Some(cells)
}
