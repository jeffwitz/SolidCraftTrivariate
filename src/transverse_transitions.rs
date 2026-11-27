//! Exact first-order transverse transition classification for extruded blocks.
//!
//! This is a census diagnostic. It proves an affine identity between first jets
//! on a face when successful; it does not prove a global reparameterization.

use crate::constraint_system::Ratio;
use crate::json::{object, Json};

type Q = Ratio;

trait CheckedRatioExt {
    fn add(self, other: Self) -> Option<Self>
    where
        Self: Sized;
    fn sub(self, other: Self) -> Option<Self>
    where
        Self: Sized;
    fn mul(self, other: Self) -> Option<Self>
    where
        Self: Sized;
    fn div(self, other: Self) -> Option<Self>
    where
        Self: Sized;
}

impl CheckedRatioExt for Ratio {
    fn add(self, other: Self) -> Option<Self> {
        self.checked_add(other).ok()
    }
    fn sub(self, other: Self) -> Option<Self> {
        self.checked_sub(other).ok()
    }
    fn mul(self, other: Self) -> Option<Self> {
        self.checked_mul(other).ok()
    }
    fn div(self, other: Self) -> Option<Self> {
        self.checked_div(other).ok()
    }
}

fn qnew(numerator: i128, denominator: i128) -> Option<Q> {
    Q::new(numerator, denominator).ok()
}

fn qtext(value: Q) -> String {
    if value.denominator() == 1 {
        value.numerator().to_string()
    } else {
        format!("{}/{}", value.numerator(), value.denominator())
    }
}

type P = Vec<Q>;
fn trim(mut p: P) -> P {
    while p.last().is_some_and(|value| value.numerator() == 0) {
        p.pop();
    }
    p
}
fn padd(a: &P, b: &P, scale: i128) -> Option<P> {
    let mut c = vec![qnew(0, 1)?; a.len().max(b.len())];
    for (i, x) in a.iter().enumerate() {
        c[i] = c[i].add(*x)?
    }
    for (i, x) in b.iter().enumerate() {
        c[i] = c[i].add(x.mul(qnew(scale, 1)?)?)?
    }
    Some(trim(c))
}
fn pmul(a: &P, b: &P) -> Option<P> {
    if a.is_empty() || b.is_empty() {
        return Some(vec![]);
    }
    let mut c = vec![qnew(0, 1)?; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            c[i + j] = c[i + j].add(x.mul(*y)?)?
        }
    }
    Some(trim(c))
}
fn deriv(a: &P) -> Option<P> {
    let mut out = Vec::new();
    for (i, x) in a.iter().enumerate().skip(1) {
        out.push(x.mul(qnew(i as i128, 1)?)?)
    }
    Some(trim(out))
}
fn cross(a: &[P; 2], b: &[P; 2]) -> Option<P> {
    padd(&pmul(&a[0], &b[1])?, &pmul(&a[1], &b[0])?, -1)
}
fn bernstein(row: &[Q]) -> Option<P> {
    if row.is_empty() {
        return Some(vec![]);
    }
    let d = row.len() - 1;
    let mut out = vec![qnew(0, 1)?; d + 1];
    for (i, v) in row.iter().enumerate() {
        let choose = binom(d, i)?;
        for j in 0..=d - i {
            let term = v.mul(qnew(
                choose
                    .checked_mul(binom(d - i, j)?)?
                    .checked_mul(if j % 2 == 0 { 1 } else { -1 })?,
                1,
            )?)?;
            out[i + j] = out[i + j].add(term)?;
        }
    }
    Some(trim(out))
}
fn binom(n: usize, k: usize) -> Option<i128> {
    let k = k.min(n - k);
    let mut r = 1i128;
    for i in 0..k {
        r = r
            .checked_mul((n - i) as i128)?
            .checked_div((i + 1) as i128)?;
    }
    Some(r)
}
fn quotient_affine(a: &P, b: &P) -> Option<Option<(Q, Q)>> {
    if b.is_empty() || b[0].numerator() == 0 {
        return None;
    }
    let c0 = if a.is_empty() {
        qnew(0, 1)?
    } else {
        a[0].div(b[0])?
    };
    let a1 = if a.len() > 1 { a[1] } else { qnew(0, 1)? };
    let b1 = if b.len() > 1 { b[1] } else { qnew(0, 1)? };
    let c1 = a1.sub(c0.mul(b1)?)?.div(b[0])?;
    let product = pmul(&vec![c0, c1], b)?;
    let residual = padd(a, &product, -1)?;
    Some(if residual.is_empty() {
        Some((c0, c1))
    } else {
        None
    })
}

fn parse_face(v: &Json) -> Result<(&str, &str), &'static str> {
    let a = v.items()?;
    if a.len() != 2 {
        return Err("expected block id and side");
    };
    Ok((a[0].text()?, a[1].text()?))
}
fn block_jet(block: &Json, side: &str, reverse: bool) -> Result<([P; 2], [P; 2], P), &'static str> {
    let deg = block.field("degrees")?.items()?;
    if deg.len() != 3 {
        return Err("expected three degrees");
    };
    let du = deg[0].number::<usize>()?;
    let dv = deg[1].number::<usize>()?;
    let dw = deg[2].number::<usize>()?;
    if du == 0 || dv == 0 || du > 64 || dv > 64 || dw != 1 {
        return Err("unsupported transition degree");
    };
    let controls = block.field("controls")?.items()?;
    if controls.len() != 2 * (du + 1) * (dv + 1)
        || controls
            .iter()
            .any(|point| point.items().map_or(true, |values| values.len() != 4))
    {
        return Err("invalid transverse control net");
    }
    let common_weight = controls
        .first()
        .and_then(|c| c.items().ok())
        .and_then(|c| c.get(3))
        .and_then(|w| w.number::<i128>().ok())
        .filter(|w| *w > 0)
        .filter(|w| {
            controls.iter().all(|c| {
                c.items()
                    .ok()
                    .and_then(|v| v.get(3))
                    .and_then(|w| w.number::<i128>().ok())
                    == Some(*w)
            })
        });
    let mut grid = vec![vec![Vec::<Q>::new(); dv + 1]; du + 1];
    for i in 0..=du {
        for j in 0..=dv {
            let v = &controls[2 * (i * (dv + 1) + j)];
            let mut c = Vec::new();
            for x in v.items()? {
                let value = x.number::<i128>()?;
                c.push(if let Some(w) = common_weight {
                    qnew(value, w).ok_or("rational overflow")?
                } else {
                    qnew(value, 1).ok_or("rational overflow")?
                })
            }
            grid[i][j] = c;
        }
    }
    let axis = match side.as_bytes() {
        [b'u', b'0'] => 0,
        [b'u', b'1'] => 0,
        [b'v', b'0'] => 1,
        [b'v', b'1'] => 1,
        _ => return Err("unsupported side"),
    };
    let end = side.as_bytes()[1] == b'1';
    let d = if axis == 0 { du } else { dv };
    let e = if end { d } else { 0 };
    let n = if end {
        d.saturating_sub(1)
    } else {
        (d + 1).min(d)
    };
    let len = if axis == 0 { dv + 1 } else { du + 1 };
    let mut h: [Vec<Q>; 4] = std::array::from_fn(|_| Vec::new());
    let mut near: [Vec<Q>; 4] = std::array::from_fn(|_| Vec::new());
    for t in 0..len {
        let (i, j, ni, nj) = if axis == 0 {
            (e, t, n, t)
        } else {
            (t, e, t, n)
        };
        for c in 0..4 {
            h[c].push(grid[i][j][c]);
            near[c].push(grid[ni][nj][c]);
        }
    }
    // Parameter reversal changes every polynomial's argument. Bernstein reversal
    // is exact by reversing the control row before conversion.
    if reverse {
        for c in 0..4 {
            h[c].reverse();
            near[c].reverse();
        }
    }
    let hp: [P; 4] = std::array::from_fn(|c| bernstein(&h[c]).unwrap_or_default());
    let np: [P; 4] = std::array::from_fn(|c| bernstein(&near[c]).unwrap_or_default());
    let weight = hp[3].clone();
    let wt = deriv(&weight).ok_or("rational overflow")?;
    let weight_normal = padd(&np[3], &hp[3], -1)
        .ok_or("rational overflow")?
        .into_iter()
        .map(|x| x.mul(qnew(d as i128, 1).unwrap()).unwrap())
        .collect::<P>();
    let mut derivs: [[P; 2]; 2] = std::array::from_fn(|_| std::array::from_fn(|_| Vec::new()));
    for c in 0..2 {
        let tangent = deriv(&hp[c]).ok_or("rational overflow")?;
        let normal = padd(&np[c], &hp[c], -1).ok_or("rational overflow")?;
        let normal = normal
            .into_iter()
            .map(|x| x.mul(qnew(d as i128, 1).unwrap()).unwrap())
            .collect::<P>();
        let dt = padd(
            &pmul(&tangent, &weight).ok_or("rational overflow")?,
            &pmul(&hp[c], &wt).ok_or("rational overflow")?,
            -1,
        )
        .ok_or("rational overflow")?;
        let dn = padd(
            &pmul(&normal, &weight).ok_or("rational overflow")?,
            &pmul(&hp[c], &weight_normal).ok_or("rational overflow")?,
            -1,
        )
        .ok_or("rational overflow")?;
        derivs[0][c] = dt;
        derivs[1][c] = dn;
    }
    let weight2 = pmul(&weight, &weight).ok_or("rational overflow")?;
    Ok((derivs[0].clone(), derivs[1].clone(), weight2))
}

type AffinePair = (Q, Q);
type DirectionalResult =
    Result<(&'static str, Option<AffinePair>, Option<AffinePair>), &'static str>;

fn directional(a: &Json, sa: &str, b: &Json, sb: &str, reverse: bool) -> DirectionalResult {
    let (ta, na, wa) = block_jet(a, sa, false)?;
    let (_, nb, wb) = block_jet(b, sb, reverse)?;
    let ca = cross(&na, &ta).ok_or("rational overflow")?;
    let cb = cross(&nb, &ta).ok_or("rational overflow")?;
    let den = pmul(&ca, &wb).ok_or("rational overflow")?;
    let nalpha = pmul(&cb, &wa).ok_or("rational overflow")?;
    let beta_num =
        pmul(&cross(&na, &nb).ok_or("rational overflow")?, &wa).ok_or("rational overflow")?;
    let Some(alpha) = quotient_affine(&nalpha, &den).ok_or("rational overflow")? else {
        return Ok(("outside_affine_transition_test", None, None));
    };
    let Some(beta) = quotient_affine(&beta_num, &den).ok_or("rational overflow")? else {
        return Ok(("outside_affine_transition_test", None, None));
    };
    if alpha
        .0
        .mul(alpha.0.add(alpha.1).ok_or("rational overflow")?)
        .ok_or("rational overflow")?
        .numerator()
        <= 0
    {
        return Ok(("transition_singular_on_face", None, None));
    }
    let kind = if alpha.1.numerator() != 0 || beta.1.numerator() != 0 {
        "affine_transition"
    } else {
        "constant_transition"
    };
    Ok((kind, Some(alpha), Some(beta)))
}
fn coeffs(alpha: Option<(Q, Q)>, beta: Option<(Q, Q)>) -> Json {
    let pair = |p: Option<(Q, Q)>| {
        p.map_or(Json::Array(vec![]), |(a, b)| {
            Json::Array(vec![Json::String(qtext(a)), Json::String(qtext(b))])
        })
    };
    if alpha.is_none() || beta.is_none() {
        object([])
    } else {
        object([("alpha", pair(alpha)), ("beta", pair(beta))])
    }
}

/// Classify every interface in one model, returning the same per-link schema as
/// `scripts/census_transverse_transitions.py`.
pub fn classify(model: &Json) -> Result<Json, &'static str> {
    let blocks = model.field("blocks")?.items()?;
    let mut by_id = std::collections::BTreeMap::new();
    for b in blocks {
        by_id.insert(b.field("id")?.text()?, b);
    }
    let mut records = Vec::new();
    for link in model.field("topology")?.field("interfaces")?.items()? {
        let first = link.field("first")?;
        let second = link.field("second")?;
        let (a, sa) = parse_face(first)?;
        let (b, sb) = parse_face(second)?;
        let orient = link.field("orientation")?;
        let swap = orient.field("swap")?.bool()?;
        let ff = orient.field("flip_first")?.bool()?;
        let fs = orient.field("flip_second")?.bool()?;
        let side_axial = sa.starts_with('w') || sb.starts_with('w');
        let aa = by_id.get(a).ok_or("unknown block")?;
        let bb = by_id.get(b).ok_or("unknown block")?;
        let mut coeff = object([]);
        let mut reverse_coeff = object([]);
        let mut reversed = None;
        let status;
        if side_axial {
            status = "axial_or_nonlateral"
        } else {
            let shift = |x: &Json| -> Result<Option<[Q; 3]>, &'static str> {
                let dg = x.field("degrees")?.items()?;
                if dg[2].number::<usize>()? != 1 {
                    return Ok(None);
                }
                let c = x.field("controls")?.items()?;
                let mut got = None;
                let (pairs, _) = c.as_chunks::<2>();
                for pair in pairs {
                    let u = pair[0].items()?;
                    let v = pair[1].items()?;
                    let mut s = [Q::zero(); 3];
                    for i in 0..3 {
                        let wu = u[3].number::<i128>()?;
                        let wv = v[3].number::<i128>()?;
                        if wu != wv || wu <= 0 {
                            return Ok(None);
                        }
                        s[i] = qnew(
                            v[i].number::<i128>()?
                                .checked_sub(u[i].number::<i128>()?)
                                .ok_or("rational overflow")?,
                            wu,
                        )
                        .ok_or("rational overflow")?;
                    }
                    if got.is_some_and(|old| old != s) {
                        return Ok(None);
                    }
                    got = Some(s);
                }
                Ok(got)
            };
            let x = shift(aa)?;
            let y = shift(bb)?;
            if swap || fs || x.is_none() || x != y {
                status = "unsupported_transition_frame"
            } else {
                let (s, c1, c2) = directional(aa, sa, bb, sb, ff)?;
                status = s;
                coeff = coeffs(c1, c2);
                let (s, c1, c2) = directional(bb, sb, aa, sa, ff)?;
                reversed = Some(s);
                reverse_coeff = coeffs(c1, c2);
            }
        }
        let face = |f: &Json| Json::Array(f.items().unwrap().to_vec());
        let mut fields = std::collections::BTreeMap::new();
        fields.insert("faces".into(), Json::Array(vec![face(first), face(second)]));
        fields.insert("status".into(), Json::String(status.into()));
        if let Json::Object(values) = coeff {
            fields.extend(values);
        }
        if let Some(reversed) = reversed {
            fields.insert("reversed_status".into(), Json::String(reversed.into()));
            fields.insert("reversed_coefficients".into(), reverse_coeff);
        }
        records.push(Json::Object(fields));
    }
    Ok(Json::Array(records))
}
