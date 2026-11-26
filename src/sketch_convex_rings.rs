//! Exact quad strips between convex host polygons and convex cut boundaries.
//!
//! The dynamic program walks both CCW boundaries monotonically. A step either
//! advances one host and one cut vertex, or advances two host vertices while
//! holding the cut vertex. Inputs are paired host/cut polygons; shared host
//! edges are kept conforming by never inserting a private point into them.

use crate::constraint_system::Ratio;
use std::collections::BTreeMap;

type Point = [Ratio; 2];
type Quad = [Point; 4];

fn cross(a: Point, b: Point, c: Point) -> Result<Ratio, &'static str> {
    b[0].checked_sub(a[0])?
        .checked_mul(c[1].checked_sub(a[1])?)?
        .checked_sub(
            b[1].checked_sub(a[1])?
                .checked_mul(c[0].checked_sub(a[0])?)?,
        )
}

fn sign(x: Ratio) -> i8 {
    x.numerator().signum() as i8
}

fn on_segment(a: Point, b: Point, p: Point) -> Result<bool, &'static str> {
    Ok(sign(cross(a, b, p)?) == 0
        && p[0].cmp_value(if a[0].cmp_value(&b[0]).is_le() {
            &a[0]
        } else {
            &b[0]
        }) != std::cmp::Ordering::Less
        && p[0].cmp_value(if a[0].cmp_value(&b[0]).is_ge() {
            &a[0]
        } else {
            &b[0]
        }) != std::cmp::Ordering::Greater
        && p[1].cmp_value(if a[1].cmp_value(&b[1]).is_le() {
            &a[1]
        } else {
            &b[1]
        }) != std::cmp::Ordering::Less
        && p[1].cmp_value(if a[1].cmp_value(&b[1]).is_ge() {
            &a[1]
        } else {
            &b[1]
        }) != std::cmp::Ordering::Greater)
}

fn proper_intersection(a: Point, b: Point, c: Point, d: Point) -> Result<bool, &'static str> {
    let x = sign(cross(a, b, c)?);
    let y = sign(cross(a, b, d)?);
    let z = sign(cross(c, d, a)?);
    let w = sign(cross(c, d, b)?);
    Ok(x * y < 0 && z * w < 0)
}

fn area2(poly: &[Point]) -> Result<Ratio, &'static str> {
    let mut a = Ratio::zero();
    for i in 0..poly.len() {
        let p = poly[i];
        let q = poly[(i + 1) % poly.len()];
        a = a.checked_add(
            p[0].checked_mul(q[1])?
                .checked_sub(p[1].checked_mul(q[0])?)?,
        )?;
    }
    Ok(a)
}

fn valid_convex(poly: &[Point]) -> Result<bool, &'static str> {
    if poly.len() < 3 || sign(area2(poly)?) <= 0 {
        return Ok(false);
    }
    for i in 0..poly.len() {
        if sign(cross(
            poly[i],
            poly[(i + 1) % poly.len()],
            poly[(i + 2) % poly.len()],
        )?) < 0
        {
            return Ok(false);
        }
    }
    Ok(true)
}

fn strictly_inside_convex(p: Point, poly: &[Point]) -> Result<bool, &'static str> {
    for i in 0..poly.len() {
        if sign(cross(poly[i], poly[(i + 1) % poly.len()], p)?) <= 0 {
            return Ok(false);
        }
    }
    Ok(true)
}

fn admissible(q: Quad, hole: &[Point]) -> Result<bool, &'static str> {
    if !valid_convex(&q)? {
        return Ok(false);
    }
    // A ring cell may share a cut edge, but must not cover any open part of
    // the cut. Convexity makes vertex containment and proper edge crossings
    // sufficient for this exact boundary test.
    for &p in hole {
        if strictly_inside_convex(p, &q)? {
            return Ok(false);
        }
    }
    for i in 0..4 {
        for j in 0..hole.len() {
            if proper_intersection(q[i], q[(i + 1) % 4], hole[j], hole[(j + 1) % hole.len()])? {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

pub(super) fn quality(q: Quad) -> f64 {
    let point = |p: Point| {
        [
            p[0].numerator() as f64 / p[0].denominator() as f64,
            p[1].numerator() as f64 / p[1].denominator() as f64,
        ]
    };
    let v = q.map(point);
    let mut worst = 0.0_f64;
    for u in [0.0, 0.5, 1.0] {
        for w in [0.0, 0.5, 1.0] {
            let du = [
                (1.0 - w) * (v[1][0] - v[0][0]) + w * (v[2][0] - v[3][0]),
                (1.0 - w) * (v[1][1] - v[0][1]) + w * (v[2][1] - v[3][1]),
            ];
            let dv = [
                (1.0 - u) * (v[3][0] - v[0][0]) + u * (v[2][0] - v[1][0]),
                (1.0 - u) * (v[3][1] - v[0][1]) + u * (v[2][1] - v[1][1]),
            ];
            let det = du[0] * dv[1] - du[1] * dv[0];
            if det <= 0.0 {
                return f64::INFINITY;
            }
            let fro = du[0] * du[0] + du[1] * du[1] + dv[0] * dv[0] + dv[1] * dv[1];
            let disc = (fro * fro - 4.0 * det * det).max(0.0).sqrt();
            worst = worst.max((fro + disc) / (2.0 * det));
        }
    }
    worst
}

fn ring_dp(host: &[Point], hole: &[Point]) -> Result<Option<Vec<Quad>>, &'static str> {
    let (m, n) = (host.len(), hole.len());
    if m + n < 6
        || !(m + n).is_multiple_of(2)
        || m > 96
        || n > 96
        || m.saturating_mul(m).saturating_mul(n) > 1_000_000
    {
        return Ok(None);
    }
    let mut best: Option<(f64, f64, Vec<Quad>)> = None;
    for start in 0..m {
        let mut dp: BTreeMap<(usize, usize), (f64, f64, Vec<Quad>)> = BTreeMap::new();
        dp.insert((0, 0), (0.0, 0.0, Vec::new()));
        for i in 0..=m {
            for j in 0..=n {
                let Some((mx, sm, path)) = dp.get(&(i, j)).cloned() else {
                    continue;
                };
                if i == m && j == n {
                    continue;
                }
                let p = |k: usize| host[(start + k) % m];
                let h = |k: usize| hole[k % n];
                let mut candidates = Vec::new();
                if i < m && j < n {
                    candidates.push(((i + 1, j + 1), [p(i), p(i + 1), h(j + 1), h(j)]));
                }
                if i + 2 <= m {
                    candidates.push(((i + 2, j), [p(i), p(i + 1), p(i + 2), h(j)]));
                }
                for (state, q) in candidates {
                    if !admissible(q, hole)? {
                        continue;
                    }
                    let k = quality(q);
                    if !k.is_finite() {
                        continue;
                    }
                    let value = (mx.max(k), sm + k);
                    let replace = dp.get(&state).is_none_or(|old| {
                        value
                            .0
                            .total_cmp(&old.0)
                            .then_with(|| value.1.total_cmp(&old.1))
                            .is_lt()
                    });
                    if replace {
                        let mut next = path.clone();
                        next.push(q);
                        dp.insert(state, (value.0, value.1, next));
                    }
                }
            }
        }
        if let Some(candidate) = dp.remove(&(m, n)) {
            let replace = best.as_ref().is_none_or(|old| {
                candidate
                    .0
                    .total_cmp(&old.0)
                    .then_with(|| candidate.1.total_cmp(&old.1))
                    .is_lt()
            });
            if replace {
                best = Some(candidate);
            }
        }
    }
    Ok(best.map(|(_, _, path)| path))
}

fn point_cmp(a: Point, b: Point) -> std::cmp::Ordering {
    a[0].cmp_value(&b[0]).then_with(|| a[1].cmp_value(&b[1]))
}

fn edge_key(a: Point, b: Point) -> (Point, Point) {
    if point_cmp(a, b).is_le() {
        (a, b)
    } else {
        (b, a)
    }
}

/// Construct conforming convex quad rings for paired CCW convex hosts and cuts.
/// Every input vertex is retained. Host vertices lying on another host edge
/// are propagated before balancing; any added midpoint is propagated to every
/// owner of its edge. Returns ring cells only; blind-cut caps are caller-owned.
pub fn plan_convex_rings(hosts: &[Vec<Point>], holes: &[Vec<Point>]) -> Option<Vec<Quad>> {
    if hosts.len() != holes.len() || hosts.is_empty() {
        return None;
    }
    let mut refined = hosts.to_vec();
    // Normalize shared host boundaries first, including pre-existing points
    // present on only one side of a shared edge.
    let all_vertices: Vec<Point> = refined.iter().flatten().copied().collect();
    for poly in &mut refined {
        let mut normalized = Vec::new();
        for i in 0..poly.len() {
            let a = poly[i];
            let b = poly[(i + 1) % poly.len()];
            normalized.push(a);
            let mut interior = Vec::new();
            for &p in &all_vertices {
                if p != a && p != b && on_segment(a, b, p).ok()? {
                    interior.push(p);
                }
            }
            interior.sort_by(|x, y| {
                let axis = if a[0] != b[0] { 0 } else { 1 };
                let order = x[axis].cmp_value(&y[axis]);
                if a[axis].cmp_value(&b[axis]).is_le() {
                    order
                } else {
                    order.reverse()
                }
            });
            interior.dedup();
            normalized.extend(interior);
        }
        *poly = normalized;
    }
    for host in &refined {
        if host.len() < 3 {
            return None;
        }
    }
    for hole in holes {
        if hole.len() < 3 || !valid_convex(hole).ok()? {
            return None;
        }
    }
    // Split an exact edge in every host that shares it. This handles parity
    // and the m<n feasibility bound without creating T-junctions on chords.
    let limit = holes
        .iter()
        .map(Vec::len)
        .sum::<usize>()
        .saturating_mul(3)
        .max(8);
    for _ in 0..limit {
        let uses = {
            let mut map = BTreeMap::<(Point, Point), Vec<usize>>::new();
            for (h, poly) in refined.iter().enumerate() {
                for i in 0..poly.len() {
                    map.entry(edge_key(poly[i], poly[(i + 1) % poly.len()]))
                        .or_default()
                        .push(h);
                }
            }
            map
        };
        let invalid: Vec<usize> = (0..refined.len())
            .filter(|&h| {
                refined[h].len() < holes[h].len()
                    || !(refined[h].len() + holes[h].len()).is_multiple_of(2)
            })
            .collect();
        if invalid.is_empty() {
            break;
        }
        let mut choice: Option<((Point, Point), Vec<usize>, Ratio)> = None;
        for &h in &invalid {
            for i in 0..refined[h].len() {
                let a = refined[h][i];
                let b = refined[h][(i + 1) % refined[h].len()];
                let key = edge_key(a, b);
                let owners = &uses[&key];
                if owners.iter().any(|&owner| !invalid.contains(&owner)) {
                    continue;
                }
                let dx = b[0].checked_sub(a[0]).ok()?;
                let dy = b[1].checked_sub(a[1]).ok()?;
                let length = dx
                    .checked_mul(dx)
                    .ok()?
                    .checked_add(dy.checked_mul(dy).ok()?)
                    .ok()?;
                if choice
                    .as_ref()
                    .is_none_or(|(_, _, best)| length.cmp_value(best).is_gt())
                {
                    choice = Some((key, owners.clone(), length));
                }
            }
        }
        let (key, owners, _) = choice?;
        let half = Ratio::new(1, 2).ok()?;
        let mid = [
            key.0[0]
                .checked_add(key.1[0])
                .ok()?
                .checked_mul(half)
                .ok()?,
            key.0[1]
                .checked_add(key.1[1])
                .ok()?
                .checked_mul(half)
                .ok()?,
        ];
        for owner in owners {
            let poly = &mut refined[owner];
            let edge =
                (0..poly.len()).find(|&i| edge_key(poly[i], poly[(i + 1) % poly.len()]) == key)?;
            poly.insert(edge + 1, mid);
        }
    }
    if (0..refined.len()).any(|h| {
        refined[h].len() < holes[h].len() || !(refined[h].len() + holes[h].len()).is_multiple_of(2)
    }) {
        return None;
    }
    // If a balanced boundary has no feasible strip, try a bounded collection
    // of exact rational Steiner pairs on host edges. Paired insertions preserve
    // the count parity; shared edges receive identical points on every owner.
    let fractions = [(1_i128, 3_i128), (2, 3), (1, 4), (3, 4), (1, 2)];
    let mut rings: Vec<Option<Vec<Quad>>> = vec![None; refined.len()];
    for _ in 0..refined.len().saturating_mul(3).max(1) {
        let mut failed = None;
        for h in 0..refined.len() {
            if !valid_convex(&refined[h]).ok()? || !valid_convex(&holes[h]).ok()? {
                return None;
            }
            for &vertex in &holes[h] {
                if !strictly_inside_convex(vertex, &refined[h]).ok()? {
                    return None;
                }
            }
            rings[h] = ring_dp(&refined[h], &holes[h]).ok()?;
            if rings[h].is_none() {
                failed = Some(h);
                break;
            }
        }
        let Some(h) = failed else {
            return Some(rings.into_iter().flatten().flatten().collect());
        };
        let mut chosen = None;
        let mut edges: Vec<(Point, Point, Ratio)> = Vec::new();
        for i in 0..refined[h].len() {
            let a = refined[h][i];
            let b = refined[h][(i + 1) % refined[h].len()];
            let dx = b[0].checked_sub(a[0]).ok()?;
            let dy = b[1].checked_sub(a[1]).ok()?;
            edges.push((
                a,
                b,
                dx.checked_mul(dx)
                    .ok()?
                    .checked_add(dy.checked_mul(dy).ok()?)
                    .ok()?,
            ));
        }
        edges.sort_by(|a, b| b.2.cmp_value(&a.2));
        'edge: for (a, b, _) in edges {
            let key = edge_key(a, b);
            let owners: Vec<usize> = (0..refined.len())
                .filter(|&owner| {
                    (0..refined[owner].len()).any(|i| {
                        edge_key(
                            refined[owner][i],
                            refined[owner][(i + 1) % refined[owner].len()],
                        ) == key
                    })
                })
                .collect();
            for i in 0..fractions.len() {
                for j in (i + 1)..fractions.len() {
                    let t = Ratio::new(fractions[i].0, fractions[i].1).ok()?;
                    let u = Ratio::new(fractions[j].0, fractions[j].1).ok()?;
                    let point = |s: Ratio| -> Option<Point> {
                        Some([
                            a[0].checked_add(b[0].checked_sub(a[0]).ok()?.checked_mul(s).ok()?)
                                .ok()?,
                            a[1].checked_add(b[1].checked_sub(a[1]).ok()?.checked_mul(s).ok()?)
                                .ok()?,
                        ])
                    };
                    let mut trial = refined.clone();
                    for &owner in &owners {
                        let poly = &mut trial[owner];
                        let e = (0..poly.len())
                            .find(|&k| edge_key(poly[k], poly[(k + 1) % poly.len()]) == key)?;
                        let p = point(t)?;
                        let q = point(u)?;
                        if point_cmp(poly[e], poly[(e + 1) % poly.len()]).is_le() {
                            poly.splice(e + 1..e + 1, [p, q]);
                        } else {
                            poly.splice(e + 1..e + 1, [q, p]);
                        }
                    }
                    if ring_dp(&trial[h], &holes[h]).ok()?.is_some() {
                        chosen = Some(trial);
                        break 'edge;
                    }
                }
            }
        }
        refined = chosen?;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    fn p(x: i128, y: i128) -> Point {
        [Ratio::new(x, 1).unwrap(), Ratio::new(y, 1).unwrap()]
    }

    #[test]
    fn rectangle_host_is_balanced_for_octagonal_cut() {
        let host = vec![p(0, 0), p(12, 0), p(12, 12), p(0, 12)];
        let hole = vec![
            p(6, 2),
            p(9, 3),
            p(10, 6),
            p(9, 9),
            p(6, 10),
            p(3, 9),
            p(2, 6),
            p(3, 3),
        ];
        let cells = plan_convex_rings(&[host], &[hole]).expect("convex annular quad strip");
        assert!(!cells.is_empty());
        assert!(cells.iter().all(|q| valid_convex(q).unwrap()));
        let expected = area2(&[p(0, 0), p(12, 0), p(12, 12), p(0, 12)])
            .unwrap()
            .checked_sub(
                area2(&[
                    p(6, 2),
                    p(9, 3),
                    p(10, 6),
                    p(9, 9),
                    p(6, 10),
                    p(3, 9),
                    p(2, 6),
                    p(3, 3),
                ])
                .unwrap(),
            )
            .unwrap();
        let actual = cells
            .iter()
            .try_fold(Ratio::zero(), |sum, q| sum.checked_add(area2(q)?))
            .unwrap();
        assert_eq!(actual, expected);
    }

    #[test]
    fn odd_cut_gets_host_midpoint_for_blind_cap_conformity() {
        let host = vec![p(0, 0), p(10, 0), p(10, 10), p(0, 10)];
        let triangle = vec![p(4, 3), p(7, 5), p(4, 7)];
        let cells =
            plan_convex_rings(&[host], &[triangle]).expect("odd cut balanced by host edge split");
        assert!(!cells.is_empty());
        assert!(cells.iter().all(|q| valid_convex(q).unwrap()));
    }
}
