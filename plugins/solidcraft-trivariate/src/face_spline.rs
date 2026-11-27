//! Extract planar SolidCraft sketch curves as homogeneous Bézier spans.
//!
//! Fit splines follow the host's natural cubic, chord-parameter interpolation. Control
//! splines use their clamped knots and positive rational weights. The returned triples are
//! `(x*w, y*w, w)` in sketch millimetres.

use serde_json::Value;
use std::collections::HashMap;

type P = [f64; 2];

fn number(v: &Value, what: &str) -> Result<f64, String> {
    let x = v
        .as_f64()
        .ok_or_else(|| format!("{what} must be numeric"))?;
    if !x.is_finite() {
        return Err(format!("{what} must be finite"));
    }
    Ok(x)
}

fn point(v: &Value) -> Result<P, String> {
    let a = v.as_array().ok_or("point position must be an array")?;
    if a.len() < 2 {
        return Err("point position needs two coordinates".into());
    }
    Ok([number(&a[0], "point x")?, number(&a[1], "point y")?])
}

fn h(p: P) -> [f64; 3] {
    [p[0], p[1], 1.0]
}
fn add(a: P, b: P) -> P {
    [a[0] + b[0], a[1] + b[1]]
}
fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1]]
}
fn mul(a: P, s: f64) -> P {
    [a[0] * s, a[1] * s]
}
fn div(a: P, s: f64) -> P {
    [a[0] / s, a[1] / s]
}
fn norm(a: P) -> f64 {
    a[0].hypot(a[1])
}

/// Natural cubic fit points, exactly the host's chord-length Hermite construction.
fn fit_spans(p: &[P]) -> Result<Vec<Vec<[f64; 3]>>, String> {
    if p.len() < 2 {
        return Err("fit spline needs at least two points".into());
    }
    if p.len() == 2 {
        return Ok(vec![vec![h(p[0]), h(p[1])]]);
    }
    let hs: Vec<f64> = p
        .windows(2)
        .map(|w| norm(sub(w[1], w[0])).max(1e-9))
        .collect();
    let m = p.len() - 2;
    let mut diag = vec![0.0; m];
    let mut up = vec![0.0; m];
    let mut low = vec![0.0; m];
    let mut rhs = vec![[0.0, 0.0]; m];
    for i in 0..m {
        let (a, b) = (hs[i], hs[i + 1]);
        diag[i] = 2.0 * (a + b);
        low[i] = a;
        up[i] = b;
        rhs[i] = mul(
            sub(div(sub(p[i + 2], p[i + 1]), b), div(sub(p[i + 1], p[i]), a)),
            6.0,
        );
    }
    for i in 1..m {
        let w = low[i] / diag[i - 1];
        diag[i] -= w * up[i - 1];
        rhs[i] = sub(rhs[i], mul(rhs[i - 1], w));
    }
    let mut sol = vec![[0.0, 0.0]; m];
    for i in (0..m).rev() {
        let next = if i + 1 < m {
            mul(sol[i + 1], up[i])
        } else {
            [0.0, 0.0]
        };
        sol[i] = div(sub(rhs[i], next), diag[i]);
    }
    let second = |i: usize| {
        if i == 0 || i == p.len() - 1 {
            [0.0, 0.0]
        } else {
            sol[i - 1]
        }
    };
    Ok((0..p.len() - 1)
        .map(|i| {
            let (a, b, hh) = (p[i], p[i + 1], hs[i]);
            let slope = div(sub(b, a), hh);
            let d0 = sub(
                slope,
                mul(add(mul(second(i), 2.0), second(i + 1)), hh / 6.0),
            );
            let d1 = add(
                slope,
                mul(add(second(i), mul(second(i + 1), 2.0)), hh / 6.0),
            );
            vec![
                h(a),
                h(add(a, mul(d0, hh / 3.0))),
                h(sub(b, mul(d1, hh / 3.0))),
                h(b),
            ]
        })
        .collect())
}

fn periodic_fit_spans(p: &[P]) -> Result<Vec<Vec<[f64; 3]>>, String> {
    if p.len() < 4 || p.first() != p.last() {
        return Err(
            "periodic fit spline needs at least three distinct samples and an exact closing point"
                .into(),
        );
    }
    let n = p.len() - 1;
    let hs: Vec<f64> = p.windows(2).map(|w| norm(sub(w[1], w[0]))).collect();
    if hs.iter().any(|h| !h.is_finite() || *h <= 1e-9) {
        return Err("periodic fit spline has a degenerate chord".into());
    }
    // Sherman–Morrison reduces the cyclic interpolation system to two
    // tridiagonal solves. Eliminate both coordinates and the correction
    // vector together; memory and work stay linear in the number of points.
    let mut diagonal = vec![0.; n];
    let mut rhs = vec![[0.; 3]; n];
    let corner = hs[n - 1];
    let gamma = -2. * (corner + hs[0]);
    for i in 0..n {
        let prev = (i + n - 1) % n;
        let (hp, hn) = (hs[prev], hs[i]);
        diagonal[i] = 2. * (hp + hn);
        let d = mul(
            sub(div(sub(p[i + 1], p[i]), hn), div(sub(p[i], p[prev]), hp)),
            6.,
        );
        rhs[i] = [d[0], d[1], 0.];
    }
    diagonal[0] -= gamma;
    diagonal[n - 1] -= corner * corner / gamma;
    rhs[0][2] = gamma;
    rhs[n - 1][2] = corner;
    for i in 1..n {
        let ratio = hs[i - 1] / diagonal[i - 1];
        diagonal[i] -= ratio * hs[i - 1];
        for k in 0..3 {
            rhs[i][k] -= ratio * rhs[i - 1][k];
        }
    }
    for i in (0..n).rev() {
        for k in 0..3 {
            if i + 1 < n {
                rhs[i][k] -= hs[i] * rhs[i + 1][k];
            }
            rhs[i][k] /= diagonal[i];
        }
    }
    let end = corner / gamma;
    let denominator = 1. + rhs[0][2] + end * rhs[n - 1][2];
    let fix: [f64; 2] = std::array::from_fn(|k| (rhs[0][k] + end * rhs[n - 1][k]) / denominator);
    let second = |i: usize| {
        [
            rhs[i][0] - fix[0] * rhs[i][2],
            rhs[i][1] - fix[1] * rhs[i][2],
        ]
    };
    let mut out = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let (start, end, hh) = (p[i], p[i + 1], hs[i]);
        let chord = sub(end, start);
        let t0 = sub(
            chord,
            mul(add(mul(second(i), 2.0), second(j)), hh * hh / 6.0),
        );
        let t1 = add(
            chord,
            mul(add(second(i), mul(second(j), 2.0)), hh * hh / 6.0),
        );
        out.push(vec![
            h(start),
            h(add(start, div(t0, 3.0))),
            h(sub(end, div(t1, 3.0))),
            h(end),
        ]);
    }
    let first = out
        .first()
        .and_then(|s| s.first())
        .copied()
        .ok_or("periodic spline has no spans")?;
    let last = out
        .last_mut()
        .and_then(|s| s.last_mut())
        .ok_or("periodic spline has no spans")?;
    *last = first;
    Ok(out)
}

fn insert_knot_once(
    p: usize,
    u: f64,
    knots: &[f64],
    c: &[[f64; 3]],
) -> Result<(Vec<f64>, Vec<[f64; 3]>), String> {
    let n = c.len() - 1;
    if knots.len() != c.len() + p + 1 {
        return Err("control spline knot count is inconsistent".into());
    }
    let mut k = p;
    while k + 1 < knots.len() && knots[k + 1] <= u {
        k += 1;
    }
    let s = knots.iter().filter(|x| **x == u).count();
    if s >= p {
        return Err("control spline knot multiplicity exceeds degree".into());
    }
    let mut q = vec![[0.0; 3]; c.len() + 1];
    for i in 0..=k - p {
        q[i] = c[i];
    }
    for i in (k - s)..=n {
        q[i + 1] = c[i];
    }
    for i in (k - p + 1)..=(k - s) {
        let den = knots[i + p] - knots[i];
        if den == 0.0 {
            return Err("degenerate control spline knot interval".into());
        }
        let a = (u - knots[i]) / den;
        for d in 0..3 {
            q[i][d] = a * c[i][d] + (1.0 - a) * c[i - 1][d];
        }
    }
    let mut out = knots.to_vec();
    out.insert(k + 1, u);
    Ok((out, q))
}

/// Split a clamped NURBS into rational Bezier pieces by exact knot insertion.
fn control_spans(
    degree: usize,
    mut knots: Vec<f64>,
    points: &[P],
    weights: &[f64],
) -> Result<Vec<Vec<[f64; 3]>>, String> {
    if degree == 0 || degree > 9 || points.len() <= degree {
        return Err("unsupported control spline degree/count".into());
    }
    if knots.len() != points.len() + degree + 1 || knots.windows(2).any(|w| w[1] < w[0]) {
        return Err("invalid control spline knots".into());
    }
    if !weights.is_empty() && weights.len() != points.len() {
        return Err("control spline needs one weight per point".into());
    }
    let weights: Vec<f64> = if weights.is_empty() {
        vec![1.0; points.len()]
    } else {
        weights.to_vec()
    };
    if weights.iter().any(|w| !w.is_finite() || *w <= 0.0) {
        return Err("control spline weights must be finite and positive".into());
    }
    let (lo, hi) = (knots[degree], knots[points.len()]);
    if !lo.is_finite()
        || !hi.is_finite()
        || hi <= lo
        || knots[..=degree].iter().any(|x| *x != lo)
        || knots[points.len()..].iter().any(|x| *x != hi)
    {
        return Err("control spline must have a finite clamped domain".into());
    }
    let mut c: Vec<[f64; 3]> = points
        .iter()
        .zip(weights)
        .map(|(p, w)| [p[0] * w, p[1] * w, w])
        .collect();
    extract_knot_spans(degree, &mut knots, &mut c, lo, hi)
}

fn extract_knot_spans(
    degree: usize,
    knots: &mut Vec<f64>,
    c: &mut Vec<[f64; 3]>,
    lo: f64,
    hi: f64,
) -> Result<Vec<Vec<[f64; 3]>>, String> {
    let internals: Vec<f64> = knots
        .iter()
        .copied()
        .filter(|u| *u >= lo && *u <= hi)
        .collect();
    for u in internals {
        while knots.iter().filter(|x| **x == u).count() < degree {
            let (next_knots, next_c) = insert_knot_once(degree, u, knots, c)?;
            *knots = next_knots;
            *c = next_c;
        }
    }
    let mut out = Vec::new();
    for i in degree..c.len() {
        if knots[i + 1] > knots[i] && knots[i] >= lo && knots[i + 1] <= hi {
            out.push(c[i - degree..=i].to_vec());
        }
    }
    if out.is_empty() {
        return Err("control spline has no nonempty spans".into());
    }
    Ok(out)
}

/// Extract all non-construction curves in the sketch, preserving each curve's native span order.
pub(super) fn spans(sketch: &Value) -> Result<Vec<Vec<[f64; 3]>>, String> {
    let root = sketch.get("sketch").unwrap_or(sketch);
    let pts = root
        .get("points")
        .and_then(Value::as_array)
        .ok_or("sketch points are missing")?;
    let mut by_id = HashMap::<String, P>::new();
    let mut by_index = Vec::new();
    for (i, v) in pts.iter().enumerate() {
        let p = point(v.get("pos").unwrap_or(v))?;
        by_index.push(p);
        if let Some(id) = v.get("id").and_then(Value::as_str) {
            by_id.insert(id.to_string(), p);
        }
        by_id.insert(i.to_string(), p);
    }
    let getp = |v: &Value| -> Result<P, String> {
        if let Some(id) = v.as_str() {
            return by_id
                .get(id)
                .copied()
                .ok_or_else(|| format!("unknown sketch point `{id}`"));
        }
        let i = v
            .as_u64()
            .ok_or("curve point reference must be an index or id")? as usize;
        by_index
            .get(i)
            .copied()
            .ok_or("curve point reference is out of range".into())
    };
    let curves = root
        .get("curves")
        .and_then(Value::as_array)
        .ok_or("sketch curves are missing")?;
    let mut out = Vec::new();
    for curve in curves {
        if curve.get("construction").and_then(Value::as_bool) == Some(true) {
            continue;
        }
        // Auto-projected face edges are linked references, not authored contour curves.
        if curve.get("link").and_then(Value::as_str).is_some() {
            continue;
        }
        match curve
            .get("type")
            .and_then(Value::as_str)
            .ok_or("curve type is missing")?
        {
            "line" => {
                let a = getp(curve.get("a").ok_or("line start is missing")?)?;
                let b = getp(curve.get("b").ok_or("line end is missing")?)?;
                out.push(vec![h(a), h(b)]);
            }
            "spline" => {
                let refs = curve
                    .get("pts")
                    .and_then(Value::as_array)
                    .ok_or("spline points are missing")?;
                let p = refs.iter().map(getp).collect::<Result<Vec<_>, _>>()?;
                if curve.get("control").and_then(Value::as_bool) == Some(true) {
                    let degree = curve.get("degree").and_then(Value::as_u64).unwrap_or(3) as usize;
                    let nurbs = curve.get("nurbs");
                    let periodic = curve.get("periodic").and_then(Value::as_bool) == Some(true);
                    if periodic && nurbs.is_none() {
                        if degree < 3 || p.len() < degree + 1 {
                            return Err("periodic control spline requires degree 3 or greater and degree+1 controls".into());
                        }
                        let n = p.len();
                        let mut wrapped = p.clone();
                        wrapped.extend(p.iter().take(degree).copied());
                        let knots: Vec<f64> = (0..=n + 2 * degree).map(|i| i as f64).collect();
                        let mut c: Vec<[f64; 3]> = wrapped.iter().map(|q| h(*q)).collect();
                        out.extend(extract_knot_spans(
                            degree,
                            &mut knots.clone(),
                            &mut c,
                            degree as f64,
                            (n + degree) as f64,
                        )?);
                        continue;
                    }
                    let knots = nurbs
                        .and_then(|n| n.get("knots"))
                        .and_then(Value::as_array)
                        .map(|a| {
                            a.iter()
                                .map(|v| number(v, "knot"))
                                .collect::<Result<Vec<_>, _>>()
                        })
                        .transpose()?
                        .unwrap_or_else(|| uniform_knots(p.len(), degree));
                    let weights = nurbs
                        .and_then(|n| n.get("weights"))
                        .and_then(Value::as_array)
                        .map(|a| {
                            a.iter()
                                .map(|v| number(v, "weight"))
                                .collect::<Result<Vec<_>, _>>()
                        })
                        .transpose()?
                        .unwrap_or_default();
                    out.extend(control_spans(degree, knots, &p, &weights)?);
                } else {
                    let periodic = curve.get("periodic").and_then(Value::as_bool) == Some(true);
                    let mut pieces = if periodic {
                        periodic_fit_spans(&p)?
                    } else {
                        fit_spans(&p)?
                    };
                    if let Some(handles) = curve
                        .get("handles")
                        .and_then(Value::as_array)
                        .filter(|a| !a.is_empty())
                    {
                        let hp = handles.iter().map(getp).collect::<Result<Vec<_>, _>>()?;
                        let count = if periodic {
                            p.len().saturating_sub(1)
                        } else {
                            p.len()
                        };
                        if hp.len() != count {
                            return Err(
                                "fit spline tangent handle count does not match its points".into(),
                            );
                        }
                        for i in 0..pieces.len() {
                            let j = if i + 1 == count { 0 } else { i + 1 };
                            let chord = norm(sub(p[i + 1], p[i])).max(1e-9);
                            let unit = |k: usize| {
                                let before = if k > 0 {
                                    Some(norm(sub(p[k], p[k - 1])))
                                } else if periodic {
                                    Some(norm(sub(p[count], p[count - 1])))
                                } else {
                                    None
                                };
                                let after = if k + 1 < p.len() {
                                    Some(norm(sub(p[k + 1], p[k])))
                                } else {
                                    None
                                };
                                match (before, after) {
                                    (Some(a), Some(b)) => (a + b) / 6.0,
                                    (Some(a), None) | (None, Some(a)) => a / 3.0,
                                    _ => 1.0 / 3.0,
                                }
                                .max(1e-9)
                            };
                            let d0 = if periodic || i + 1 < p.len() {
                                div(sub(hp[i], p[i]), unit(i))
                            } else {
                                div(sub(p[i], hp[i]), unit(i))
                            };
                            let d1 = if periodic || j + 1 < p.len() {
                                div(sub(hp[j], p[j]), unit(j))
                            } else {
                                div(sub(p[j], hp[j]), unit(j))
                            };
                            pieces[i][1] = h(add(p[i], mul(d0, chord / 3.0)));
                            pieces[i][2] = h(sub(p[i + 1], mul(d1, chord / 3.0)));
                        }
                    }
                    out.extend(pieces);
                }
            }
            other => return Err(format!("unsupported sketch contour curve `{other}`")),
        }
    }
    if out.is_empty() {
        return Err("sketch has no supported contour curves".into());
    }
    Ok(out)
}

fn uniform_knots(n: usize, degree: usize) -> Vec<f64> {
    let d = degree.min(n.saturating_sub(1));
    let spans = n - d;
    let mut k = vec![0.0; d + 1];
    k.extend((1..spans).map(|i| i as f64 / spans as f64));
    k.extend(std::iter::repeat(1.0).take(d + 1));
    k
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn extracts_line_and_host_fit_spans() {
        let s = serde_json::json!({"points":[{"id":"a","pos":[0,0]},{"id":"b","pos":[1,0]},{"id":"c","pos":[2,1]},{"id":"d","pos":[3,0]}],"curves":[{"type":"line","a":"a","b":"b"},{"type":"spline","pts":["b","c","d"],"control":false}]});
        let b = spans(&s).unwrap();
        assert_eq!(b.len(), 3);
        assert_eq!(b[0].len(), 2);
        assert_eq!(b[1].len(), 4);
        assert_eq!(b[2].len(), 4);
        assert_eq!(b[1][0], [1.0, 0.0, 1.0]);
        assert_eq!(b[2][3], [3.0, 0.0, 1.0]);
    }
    #[test]
    fn extracts_rational_control_bezier_and_uniform_curve() {
        let s = serde_json::json!({"points":[{"id":"a","pos":[0,0]},{"id":"b","pos":[1,2]},{"id":"c","pos":[2,0]}],"curves":[{"type":"spline","control":true,"degree":2,"pts":["a","b","c"],"nurbs":{"knots":[0,0,0,1,1,1],"weights":[1,0.5,1]}}]});
        let b = spans(&s).unwrap();
        assert_eq!(b.len(), 1);
        assert_eq!(b[0].len(), 3);
        assert_eq!(b[0][1], [0.5, 1.0, 0.5]);
    }
    #[test]
    fn rational_internal_knot_spans_match_independent_deboor_points() {
        let controls = [
            [0.0, 0.0, 1.0],
            [2.0, 1.0, 0.5],
            [3.0, 4.0, 1.0],
            [8.0, 2.0, 2.0],
        ];
        let knots = [0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0];
        let projected = [[0.0, 0.0], [4.0, 2.0], [3.0, 4.0], [4.0, 1.0]];
        let pieces = control_spans(2, knots.to_vec(), &projected, &[1.0, 0.5, 1.0, 2.0]).unwrap();
        assert_eq!(pieces.len(), 2);
        fn de_boor(degree: usize, knots: &[f64], cp: &[[f64; 3]], u: f64) -> [f64; 2] {
            let n = cp.len();
            let mut k = degree;
            while k + 1 < n && knots[k + 1] <= u {
                k += 1;
            }
            let mut q = cp[k - degree..=k].to_vec();
            for r in 1..=degree {
                for j in (r..=degree).rev() {
                    let i = k - degree + j;
                    let a = (u - knots[i]) / (knots[i + degree + 1 - r] - knots[i]);
                    for c in 0..3 {
                        q[j][c] = (1.0 - a) * q[j - 1][c] + a * q[j][c];
                    }
                }
            }
            [q[degree][0] / q[degree][2], q[degree][1] / q[degree][2]]
        }
        fn bezier(cp: &[[f64; 3]], t: f64) -> [f64; 2] {
            let mut q = cp.to_vec();
            for r in 1..q.len() {
                for j in 0..q.len() - r {
                    for c in 0..3 {
                        q[j][c] = (1.0 - t) * q[j][c] + t * q[j + 1][c];
                    }
                }
            }
            [q[0][0] / q[0][2], q[0][1] / q[0][2]]
        }
        for (span, (lo, hi)) in pieces.iter().zip([(0.0, 1.0), (1.0, 2.0)]) {
            for t in [0.0, 0.13, 0.5, 0.87, 1.0] {
                let a = bezier(span, t);
                let b = de_boor(2, &knots, &controls, lo + (hi - lo) * t);
                assert!(
                    (a[0] - b[0]).abs() < 2e-13 && (a[1] - b[1]).abs() < 2e-13,
                    "{a:?} != {b:?}"
                );
            }
        }
    }
    #[test]
    fn rejects_bad_weight() {
        let s = serde_json::json!({"points":[{"id":"a","pos":[0,0]},{"id":"b","pos":[1,1]},{"id":"c","pos":[2,0]}],"curves":[{"type":"spline","control":true,"degree":2,"pts":["a","b","c"],"nurbs":{"knots":[0,0,0,1,1,1],"weights":[1,0,1]}}]});
        assert!(spans(&s).is_err());
    }
    #[test]
    fn extracts_native_closed_fit_control_and_ignores_projected_reference() {
        let mut pts = vec![serde_json::json!({"id":"origin","pos":[0,0]})];
        for (i, p) in [
            [1, 1],
            [2, 1],
            [3, 1],
            [4, 1],
            [10, 10],
            [15, 7],
            [20, 10],
            [15, 13],
            [35, 10],
            [40, 7],
            [45, 10],
            [40, 13],
        ]
        .iter()
        .enumerate()
        {
            pts.push(serde_json::json!({"id":format!("p{}",i+1),"pos":p}));
        }
        let s = serde_json::json!({"sketch":{"points":pts,"curves":[
            {"type":"line","id":"projection","link":"j1","a":1,"b":2},
            {"type":"spline","id":"s1","pts":[5,6,7,8,5],"periodic":true,"degree":3},
            {"type":"spline","id":"s2","pts":[9,10,11,12],"periodic":true,"control":true,"degree":3}
        ]}});
        let b = spans(&s).unwrap();
        assert_eq!(b.len(), 8);
        assert!(b.iter().all(|span| span.len() == 4));
        for i in 0..4 {
            assert_eq!(b[i].last().unwrap(), b[(i + 1) % 4].first().unwrap());
        }
        for i in 4..8 {
            assert_eq!(b[i].last().unwrap(), b[4 + (i + 1) % 4].first().unwrap());
        }
    }
    #[test]
    fn extracts_open_degree_five_host_control_spline() {
        let pts: Vec<_> = (0..6)
            .map(|i| serde_json::json!({"id":format!("p{}",i),"pos":[i,i*i]}))
            .collect();
        let s = serde_json::json!({"points":pts,"curves":[{"type":"spline","control":true,"degree":5,"pts":["p0","p1","p2","p3","p4","p5"]}]});
        let b = spans(&s).unwrap();
        assert_eq!(b.len(), 1);
        assert_eq!(b[0].len(), 6);
    }

    #[test]
    fn extracts_closed_degree_five_periodic_control_spline() {
        let pts: Vec<_> = (0..6)
            .map(|i| serde_json::json!({"id":format!("p{}",i),"pos":[i,(i*i)%7]}))
            .collect();
        let refs: Vec<_> = (0..6).map(|i| serde_json::json!(format!("p{i}"))).collect();
        let s = serde_json::json!({"points":pts,"curves":[{"type":"spline","control":true,"periodic":true,"degree":5,"pts":refs}]});
        let b = spans(&s).unwrap();
        assert_eq!(b.len(), 6);
        assert!(b.iter().all(|span| span.len() == 6));
        for i in 0..b.len() {
            let a = b[i].last().unwrap();
            let c = b[(i + 1) % b.len()].first().unwrap();
            assert!((a[0] / a[2] - c[0] / c[2]).abs() < 1e-13);
            assert!((a[1] / a[2] - c[1] / c[2]).abs() < 1e-13);
        }
    }
}
