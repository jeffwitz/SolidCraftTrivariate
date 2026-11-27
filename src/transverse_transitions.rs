//! Exact first-order transverse transition classification for extruded blocks.
//!
//! This is a census diagnostic. It proves an affine identity between first jets
//! on a face when successful; it does not prove a global reparameterization.

use crate::constraint_system::Ratio;
use crate::json::{Json, object};

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

fn pdiv_rem(a: &P, b: &P) -> Option<(P, P)> {
    if b.is_empty() {
        return None;
    }
    let mut r = a.clone();
    let mut q = vec![qnew(0, 1)?; a.len().saturating_sub(b.len()) + 1];
    while !r.is_empty() && r.len() >= b.len() {
        let shift = r.len() - b.len();
        let factor = r.last()?.div(*b.last()?)?;
        q[shift] = factor;
        let mut term = vec![qnew(0, 1)?; shift];
        term.extend(
            b.iter()
                .map(|x| x.mul(factor))
                .collect::<Option<Vec<_>>>()?,
        );
        r = padd(&r, &term, -1)?;
    }
    Some((trim(q), r))
}

fn pgcd(a: &P, b: &P) -> Option<P> {
    let (mut x, mut y) = (trim(a.clone()), trim(b.clone()));
    while !y.is_empty() {
        let (_, r) = pdiv_rem(&x, &y)?;
        x = y;
        y = r;
    }
    if x.is_empty() {
        return Some(x);
    }
    let lead = *x.last()?;
    x.iter().map(|v| v.div(lead)).collect()
}

fn split_beta(a: &P, b: &P, beta: &P) -> Option<(P, P)> {
    // Solve beta = a * beta_second + b * beta_first for linear factors.
    let z = qnew(0, 1)?;
    let at = |p: &P, i: usize| p.get(i).copied().unwrap_or(z);
    let mut m = [
        [at(b, 0), z, at(a, 0), z, at(beta, 0)],
        [at(b, 1), at(b, 0), at(a, 1), at(a, 0), at(beta, 1)],
        [z, at(b, 1), z, at(a, 1), at(beta, 2)],
    ];
    let mut row = 0;
    for col in 0..4 {
        let Some(pivot) = (row..3).find(|r| m[*r][col].numerator() != 0) else {
            continue;
        };
        m.swap(row, pivot);
        let d = m[row][col];
        for v in &mut m[row][col..] {
            *v = v.div(d)?;
        }
        let normalized = m[row];
        for (r, other) in m.iter_mut().enumerate() {
            if r == row {
                continue;
            }
            let f = other[col];
            for (v, pivot) in other[col..].iter_mut().zip(&normalized[col..]) {
                *v = v.sub(f.mul(*pivot)?)?;
            }
        }
        row += 1;
        if row == 3 {
            break;
        }
    }
    for trailing in &m[row..] {
        if trailing[..4].iter().all(|v| v.numerator() == 0) && trailing[4].numerator() != 0 {
            return None;
        }
    }
    let mut x = [z; 4];
    for r in (0..row).rev() {
        let col = (0..4).find(|c| m[r][*c].numerator() != 0)?;
        let rhs = (col + 1..4).try_fold(m[r][4], |v, c| v.sub(m[r][c].mul(x[c])?))?;
        x[col] = rhs.div(m[r][col])?;
    }
    Some((vec![x[0], x[1]], vec![x[2], x[3]]))
}

fn normalized_control_points(block: &Json) -> Result<Option<Vec<[Q; 3]>>, &'static str> {
    let controls = block.field("controls")?.items()?;
    let Some(weight) = controls
        .first()
        .and_then(|p| p.items().ok())
        .and_then(|p| p.get(3))
        .and_then(|w| w.number::<i128>().ok())
    else {
        return Ok(None);
    };
    if weight <= 0 {
        return Ok(None);
    }
    let mut points = Vec::with_capacity(controls.len());
    for control in controls {
        let Ok(values) = control.items() else {
            return Ok(None);
        };
        if values.len() != 4
            || values[3].number::<i128>().ok() != Some(weight)
            || values[3].number::<i128>().ok().is_none_or(|w| w <= 0)
        {
            return Ok(None);
        }
        let mut xyz = [Q::zero(); 3];
        for k in 0..3 {
            let Ok(value) = values[k].number::<i128>() else {
                return Ok(None);
            };
            xyz[k] = qnew(value, weight).ok_or("rational overflow")?;
        }
        points.push(xyz);
    }
    Ok(Some(points))
}

fn has_constant_positive_weight(block: &Json) -> Result<bool, &'static str> {
    let controls = block.field("controls")?.items()?;
    let Some(weight) = controls
        .first()
        .and_then(|p| p.items().ok())
        .and_then(|p| p.get(3))
        .and_then(|w| w.number::<i128>().ok())
    else {
        return Ok(false);
    };
    Ok(weight > 0
        && controls.iter().all(|p| {
            p.items()
                .ok()
                .is_some_and(|v| v.len() == 4 && v[3].number::<i128>().ok() == Some(weight))
        }))
}

fn has_planar_extrusion_support(block: &Json) -> Result<bool, &'static str> {
    let degrees = block.field("degrees")?.items()?;
    if degrees.len() != 3 || degrees.iter().any(|d| d.number::<usize>().ok() != Some(1)) {
        return Ok(false);
    }
    let Some(points) = normalized_control_points(block)? else {
        return Ok(false);
    };
    if points.len() != 8 {
        return Ok(false);
    }
    let mut base_z = None;
    let mut height = None;
    for i in 0..4 {
        let low = points[2 * i];
        let high = points[2 * i + 1];
        if base_z.is_some_and(|z| z != low[2]) {
            return Ok(false);
        }
        base_z = Some(low[2]);
        let delta = [
            high[0].sub(low[0]).ok_or("rational overflow")?,
            high[1].sub(low[1]).ok_or("rational overflow")?,
            high[2].sub(low[2]).ok_or("rational overflow")?,
        ];
        if delta[0].numerator() != 0 || delta[1].numerator() != 0 || delta[2].numerator() <= 0 {
            return Ok(false);
        }
        if height.is_some_and(|h| h != delta[2]) {
            return Ok(false);
        }
        height = Some(delta[2]);
    }
    Ok(true)
}

fn oriented_endpoints(
    block: &Json,
    side: &str,
    reverse: bool,
) -> Result<Option<[[[Q; 3]; 2]; 2]>, &'static str> {
    let degrees = block.field("degrees")?.items()?;
    if degrees.len() != 3 || degrees.iter().any(|d| d.number::<usize>().ok() != Some(1)) {
        return Ok(None);
    }
    let Some(points) = normalized_control_points(block)? else {
        return Ok(None);
    };
    if points.len() != 8 {
        return Ok(None);
    }
    let (axis, end) = match side.as_bytes() {
        [b'u', b'0'] => (0, 0),
        [b'u', b'1'] => (0, 1),
        [b'v', b'0'] => (1, 0),
        [b'v', b'1'] => (1, 1),
        _ => return Ok(None),
    };
    let mut result = [[[Q::zero(); 3]; 2]; 2];
    for (z, z_row) in result.iter_mut().enumerate() {
        for (t0, endpoint) in z_row.iter_mut().enumerate() {
            let t = if reverse { 1 - t0 } else { t0 };
            let (i, j) = if axis == 0 { (end, t) } else { (t, end) };
            *endpoint = points[2 * (2 * i + j) + z];
        }
    }
    Ok(Some(result))
}

fn as_g1_result(
    a: &Json,
    sa: &str,
    b: &Json,
    sb: &str,
    reverse: bool,
) -> Result<Json, &'static str> {
    let unsupported = || object([("status", Json::String("unsupported_family".into()))]);
    let da = a.field("degrees")?.items()?;
    let db = b.field("degrees")?.items()?;
    if da.len() != 3
        || db.len() != 3
        || da
            .iter()
            .chain(db)
            .any(|d| d.number::<usize>().ok() != Some(1))
    {
        return Ok(unsupported());
    }
    if !has_constant_positive_weight(a)? || !has_constant_positive_weight(b)? {
        return Ok(object([(
            "status",
            Json::String("unsupported_projective_weights".into()),
        )]));
    }
    if !has_planar_extrusion_support(a)? || !has_planar_extrusion_support(b)? {
        return Ok(object([(
            "status",
            Json::String("unsupported_extrusion_support".into()),
        )]));
    }
    let Some(trace_a) = oriented_endpoints(a, sa, false)? else {
        return Ok(object([(
            "status",
            Json::String("unsupported_trace".into()),
        )]));
    };
    let Some(trace_b) = oriented_endpoints(b, sb, reverse)? else {
        return Ok(object([(
            "status",
            Json::String("unsupported_trace".into()),
        )]));
    };
    if trace_a != trace_b {
        return Ok(object([(
            "status",
            Json::String("unsupported_trace".into()),
        )]));
    }
    let (ta, na, _) = block_jet(a, sa, false)?;
    let (tb, nb, _) = block_jet(b, sb, reverse)?;
    if tb != ta {
        return Ok(object([(
            "status",
            Json::String("unsupported_trace".into()),
        )]));
    }
    // The determinant syzygy constructs gluing data:
    // det(n_b,t)n_a + det(t,n_a)n_b + det(n_a,n_b)t = 0.
    let mut alpha_a = cross(&nb, &ta).ok_or("rational overflow")?;
    let mut alpha_b = cross(&ta, &na).ok_or("rational overflow")?;
    let mut beta = cross(&na, &nb).ok_or("rational overflow")?;
    let common = pgcd(&pgcd(&alpha_a, &alpha_b).ok_or("rational overflow")?, &beta)
        .ok_or("rational overflow")?;
    if common.is_empty() {
        return Ok(object([("status", Json::String("indeterminate".into()))]));
    }
    for p in [&mut alpha_a, &mut alpha_b, &mut beta] {
        let (q, r) = pdiv_rem(p, &common).ok_or("rational overflow")?;
        if !r.is_empty() {
            return Err("inexact gluing factor division");
        }
        *p = q;
    }
    let coeff = |p: &P| -> Json {
        Json::Array(if p.is_empty() {
            vec![Json::String("0".into())]
        } else {
            p.iter().map(|x| Json::String(qtext(*x))).collect()
        })
    };
    // AS-G1 permits quadratic beta; the two alpha gluing functions are linear.
    let linear = alpha_a.len() <= 2 && alpha_b.len() <= 2 && beta.len() <= 3;
    let relatively_prime = pgcd(&alpha_a, &alpha_b).is_some_and(|p| p.len() == 1);
    let Some((beta_first, beta_second)) = split_beta(&alpha_a, &alpha_b, &beta) else {
        return Ok(object([
            ("status", Json::String("not_as_g1".into())),
            ("reason", Json::String("no_linear_beta_split".into())),
        ]));
    };
    let split_check = padd(
        &padd(
            &pmul(&alpha_a, &beta_second).ok_or("rational overflow")?,
            &pmul(&alpha_b, &beta_first).ok_or("rational overflow")?,
            1,
        )
        .ok_or("rational overflow")?,
        &beta,
        -1,
    )
    .ok_or("rational overflow")?;
    let identity_check = (0..2).all(|c| {
        let left = pmul(&alpha_a, &na[c])
            .and_then(|x| pmul(&alpha_b, &nb[c]).and_then(|y| padd(&x, &y, 1)))
            .and_then(|x| pmul(&beta, &ta[c]).and_then(|y| padd(&x, &y, 1)));
        left.is_some_and(|p| p.is_empty())
    });
    let endpoint_signs = |p: &P| -> Option<(i8, i8)> {
        let x = p.first()?.numerator();
        let y = p
            .iter()
            .try_fold(qnew(0, 1)?, |sum, c| sum.add(*c))?
            .numerator();
        if x == 0 || y == 0 || x.signum() != y.signum() {
            None
        } else {
            Some((x.signum() as i8, y.signum() as i8))
        }
    };
    // Normals are directed inward on each chart; compatible inward gluing has
    // equal-sign alpha factors because the two physical inward vectors oppose.
    let orientable = endpoint_signs(&alpha_a)
        .zip(endpoint_signs(&alpha_b))
        .is_some_and(|(x, y)| x.0 == y.0 && x.1 == y.1);
    let no_zero = |p: &P| -> Option<bool> {
        let (s0, s1) = endpoint_signs(p)?;
        if p.len() <= 2 {
            return Some(s0 == s1);
        }
        if p.len() != 3 {
            return None;
        }
        let (a, b) = (p[2], p[1]);
        if a.numerator() == 0 {
            return Some(s0 == s1);
        }
        let vertex = b.mul(qnew(-1, 2)?)?.div(a)?;
        if vertex.numerator() <= 0 || vertex.sub(qnew(1, 1)?)?.numerator() >= 0 {
            return Some(s0 == s1);
        }
        let value = p[0].sub(b.mul(b)?.div(a.mul(qnew(4, 1)?)?)?)?;
        Some(if s0 > 0 {
            value.numerator() > 0
        } else {
            value.numerator() < 0
        })
    };
    let regular = no_zero(&cross(&ta, &na).ok_or("rational overflow")?).unwrap_or(false)
        && no_zero(&cross(&ta, &nb).ok_or("rational overflow")?).unwrap_or(false);
    let admitted = linear
        && relatively_prime
        && orientable
        && regular
        && split_check.is_empty()
        && identity_check;
    Ok(object([
        (
            "status",
            Json::String(if admitted { "as_g1" } else { "not_as_g1" }.into()),
        ),
        (
            "gluing_factors",
            object([
                ("alpha_first", coeff(&alpha_a)),
                ("alpha_second", coeff(&alpha_b)),
                ("beta_first", coeff(&beta_first)),
                ("beta_second", coeff(&beta_second)),
                ("beta", coeff(&beta)),
            ]),
        ),
        ("common_factor", coeff(&common)),
        ("relatively_prime", Json::Bool(relatively_prime)),
        ("orientable", Json::Bool(orientable)),
        ("regular_on_closed_edge", Json::Bool(regular)),
        ("trace_tangent_match", Json::Bool(true)),
        ("beta_split_verified", Json::Bool(split_check.is_empty())),
        ("syzygy_verified", Json::Bool(identity_check)),
        (
            "proof",
            Json::String(
                "exact determinant syzygy; polynomial gcd division has zero remainder".into(),
            ),
        ),
    ]))
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
    let convert = |rows: [Vec<Q>; 4]| -> Result<[P; 4], &'static str> {
        rows.into_iter()
            .map(|row| bernstein(&row).ok_or("rational overflow"))
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .map_err(|_| "invalid transverse row count")
    };
    let hp = convert(h)?;
    let np = convert(near)?;
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
        let mut as_g1 = object([("status", Json::String("unsupported_frame".into()))]);
        let status;
        if side_axial {
            status = "axial_or_nonlateral";
            as_g1 = object([("status", Json::String("nonlateral_interface".into()))]);
        } else {
            let shift = |x: &Json| -> Result<Option<[Q; 3]>, &'static str> {
                let Ok(dg) = x.field("degrees").and_then(Json::items) else {
                    return Ok(None);
                };
                if dg.len() != 3 {
                    return Ok(None);
                }
                let Ok(du) = dg[0].number::<usize>() else {
                    return Ok(None);
                };
                let Ok(dv) = dg[1].number::<usize>() else {
                    return Ok(None);
                };
                if du == 0
                    || dv == 0
                    || du > 64
                    || dv > 64
                    || dg[2].number::<usize>().ok() != Some(1)
                {
                    return Ok(None);
                }
                let Some(expected) = (du + 1).checked_mul(dv + 1).and_then(|m| m.checked_mul(2))
                else {
                    return Ok(None);
                };
                let Ok(c) = x.field("controls").and_then(Json::items) else {
                    return Ok(None);
                };
                if c.len() != expected || c.iter().any(|p| p.items().map_or(true, |v| v.len() != 4))
                {
                    return Ok(None);
                }
                let mut got = None;
                let (pairs, remainder) = c.as_chunks::<2>();
                if !remainder.is_empty() {
                    return Ok(None);
                }
                for pair in pairs {
                    let u = pair[0].items()?;
                    let v = pair[1].items()?;
                    let mut s = [Q::zero(); 3];
                    for i in 0..3 {
                        let Ok(wu) = u[3].number::<i128>() else {
                            return Ok(None);
                        };
                        let Ok(wv) = v[3].number::<i128>() else {
                            return Ok(None);
                        };
                        if wu != wv || wu <= 0 {
                            return Ok(None);
                        }
                        let Ok(lower) = u[i].number::<i128>() else {
                            return Ok(None);
                        };
                        let Ok(upper) = v[i].number::<i128>() else {
                            return Ok(None);
                        };
                        s[i] = qnew(upper.checked_sub(lower).ok_or("rational overflow")?, wu)
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
                as_g1 = as_g1_result(aa, sa, bb, sb, ff)?;
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
        fields.insert("as_g1".into(), as_g1);
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
