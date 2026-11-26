//! Certified piecewise cubic approximation of the named S3 offset rail.
//!
//! The certificate bounds the interpolation error to the exact constant-distance
//! rail using interval automatic differentiation on each dyadic parameter span.

const TOL: f64 = 1.0e-3;
const RADIUS: f64 = 0.25;

#[derive(Clone, Copy, Debug)]
struct I {
    lo: f64,
    hi: f64,
}
impl I {
    fn new(lo: f64, hi: f64) -> Self {
        Self {
            lo: down(lo),
            hi: up(hi),
        }
    }
    fn add(self, b: Self) -> Self {
        Self::new(self.lo + b.lo, self.hi + b.hi)
    }
    fn neg(self) -> Self {
        Self::new(-self.hi, -self.lo)
    }
    fn sub(self, b: Self) -> Self {
        self.add(b.neg())
    }
    fn mul(self, b: Self) -> Self {
        let v = [
            self.lo * b.lo,
            self.lo * b.hi,
            self.hi * b.lo,
            self.hi * b.hi,
        ];
        Self::new(
            v.iter().copied().fold(f64::INFINITY, f64::min),
            v.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        )
    }
    fn reciprocal(self) -> Option<Self> {
        if self.lo <= 0.0 {
            None
        } else {
            Some(Self::new(1.0 / self.hi, 1.0 / self.lo))
        }
    }
    fn sqrt(self) -> Option<Self> {
        if self.lo <= 0.0 {
            None
        } else {
            Some(Self::new(self.lo.sqrt(), self.hi.sqrt()))
        }
    }
    fn absmax(self) -> f64 {
        self.lo.abs().max(self.hi.abs())
    }
}
fn down(x: f64) -> f64 {
    x - (x.abs() + 1.0) * 8.0 * f64::EPSILON
}
fn up(x: f64) -> f64 {
    x + (x.abs() + 1.0) * 8.0 * f64::EPSILON
}

// Taylor jets carry interval enclosures of derivatives divided by factorial.
#[derive(Clone, Copy)]
struct K {
    c: [I; 5],
}
impl K {
    fn c(x: f64) -> Self {
        let mut c = [I::new(0., 0.); 5];
        c[0] = I::new(x, x);
        Self { c }
    }
    fn variable(a: f64, b: f64) -> Self {
        let mut c = [I::new(0., 0.); 5];
        c[0] = I::new(a, b);
        c[1] = I::new(1., 1.);
        Self { c }
    }
    fn add(self, b: Self) -> Self {
        Self {
            c: std::array::from_fn(|i| self.c[i].add(b.c[i])),
        }
    }
    fn neg(self) -> Self {
        Self {
            c: self.c.map(I::neg),
        }
    }
    fn sub(self, b: Self) -> Self {
        self.add(b.neg())
    }
    fn mul(self, b: Self) -> Self {
        Self {
            c: std::array::from_fn(|n| {
                (0..=n).fold(I::new(0., 0.), |s, k| s.add(self.c[k].mul(b.c[n - k])))
            }),
        }
    }
    fn inv(self) -> Option<Self> {
        let c0 = self.c[0].reciprocal()?;
        let mut c = [I::new(0., 0.); 5];
        c[0] = c0;
        for n in 1..=4 {
            let sum = (1..=n).fold(I::new(0., 0.), |s, k| s.add(self.c[k].mul(c[n - k])));
            c[n] = sum.neg().mul(c0);
        }
        Some(Self { c })
    }
    fn sqrt(self) -> Option<Self> {
        let mut c = [I::new(0., 0.); 5];
        c[0] = self.c[0].sqrt()?;
        let denom = c[0].mul(I::new(2., 2.)).reciprocal()?;
        for n in 1..=4 {
            let sum = (1..n).fold(I::new(0., 0.), |s, k| s.add(c[k].mul(c[n - k])));
            c[n] = self.c[n].sub(sum).mul(denom);
        }
        Some(Self { c })
    }
}
fn poly4(c: &[f64], t: K) -> K {
    c.iter().rev().fold(K::c(0.), |a, &x| a.mul(t).add(K::c(x)))
}
fn offset_jet(t0: f64, t1: f64, side: f64) -> Result<[K; 2], &'static str> {
    let t = K::variable(t0, t1);
    let (x, y, w) = (poly4(X, t), poly4(Y, t), poly4(W, t));
    let xp = poly4(&[X[1], 2. * X[2], 3. * X[3]], t);
    let yp = poly4(&[Y[1], 2. * Y[2], 3. * Y[3]], t);
    let wp = poly4(&[W[1], 2. * W[2], 3. * W[3]], t);
    let dx = xp.mul(w).sub(x.mul(wp));
    let dy = yp.mul(w).sub(y.mul(wp));
    let speed = dx
        .mul(dx)
        .add(dy.mul(dy))
        .sqrt()
        .ok_or("speed interval includes zero")?;
    let invw = w.inv().ok_or("weight interval includes zero")?;
    let invs = speed.inv().ok_or("speed interval includes zero")?;
    Ok([
        x.mul(invw).sub(dy.mul(invs).mul(K::c(side * RADIUS))),
        y.mul(invw).add(dx.mul(invs).mul(K::c(side * RADIUS))),
    ])
}
fn fourth_derivative_bound(a: f64, b: f64, side: f64) -> Result<f64, &'static str> {
    // Subdivide only the certificate domain; this never adds a geometric span.
    let mut bound: f64 = 0.0;
    let mut domains = Vec::new();
    for i in 0..16 {
        domains.push((
            a + (b - a) * (i as f64 / 16.0),
            a + (b - a) * ((i + 1) as f64 / 16.0),
            0u8,
        ));
    }
    while let Some((lo, hi, depth)) = domains.pop() {
        match offset_jet(lo, hi, side) {
            Ok(jets) => {
                bound = bound.max(24.0 * jets[0].c[4].absmax().hypot(jets[1].c[4].absmax()))
            }
            Err("weight interval includes zero" | "speed interval includes zero") if depth < 20 => {
                let m = (lo + hi) * 0.5;
                domains.push((m, hi, depth + 1));
                domains.push((lo, m, depth + 1));
            }
            Err(e) => return Err(e),
        }
    }
    Ok(up(bound))
}

// Homogeneous source span from the exact feasibility witness.
const X: &[f64] = &[824., -30., -516., 322.];
const Y: &[f64] = &[324., 30., -216., 134.];
const W: &[f64] = &[10., 0., -6., 4.];

fn curve(t: f64, side: f64) -> Result<[f64; 2], &'static str> {
    // The terminal tangent is exactly horizontal, so its contact point is rational.
    if t == 1.0 {
        return Ok([75.0, 34.0 - side * RADIUS]);
    }
    let r = offset_jet(t, t, side)?;
    Ok([r[0].c[0].lo, r[1].c[0].lo])
}

#[derive(Clone, Debug)]
pub struct RailSpan {
    pub interval: [f64; 2],
    pub controls: [[f64; 3]; 4],
    pub coefficient_error_bound_mm: f64,
    pub error_bound_mm: f64,
}
#[derive(Clone, Debug)]
pub struct OffsetRail {
    pub source_block: &'static str,
    pub side: i8,
    pub radius_mm: f64,
    pub spans: Vec<RailSpan>,
    pub continuity: &'static str,
}

/// Construct one signed contact rail over the complete named source span.
pub fn build(side: i8) -> Result<OffsetRail, &'static str> {
    if ![-1, 1].contains(&side) {
        return Err("rail side must be -1 or 1");
    }
    let mut spans = Vec::new();
    let mut stack = vec![(0.0, 1.0, 0u8)];
    while let Some((a, b, depth)) = stack.pop() {
        let h = b - a;
        let fourth = match fourth_derivative_bound(a, b, side as f64) {
            Ok(v) => v,
            Err("weight interval includes zero" | "speed interval includes zero") if depth < 30 => {
                let m = (a + b) * 0.5;
                stack.push((m, b, depth + 1));
                stack.push((a, m, depth + 1));
                continue;
            }
            Err(e) => return Err(e),
        };
        // Cubic Hermite remainder: max |(t-a)^2(t-b)^2| / 4! = h^4/384.
        let interpolation_error = up(fourth * h.powi(4) / 384.0);
        if interpolation_error > TOL {
            if depth >= 30 {
                return Err("dyadic subdivision budget exhausted");
            }
            let m = (a + b) * 0.5;
            stack.push((m, b, depth + 1));
            stack.push((a, m, depth + 1));
            continue;
        }
        let p = curve(a, side as f64)?;
        let q = curve(b, side as f64)?;
        let start = offset_jet(a, a, side as f64)?;
        let d0 = [
            start[0].c[1].lo + start[0].c[1].hi,
            start[1].c[1].lo + start[1].c[1].hi,
        ]
        .map(|x| x * 0.5);
        let end = offset_jet(b, b, side as f64)?;
        let d1 = [
            end[0].c[1].lo + end[0].c[1].hi,
            end[1].c[1].lo + end[1].c[1].hi,
        ]
        .map(|x| x * 0.5);
        let controls = [
            [p[0], p[1], 1.],
            [p[0] + h * d0[0] / 3., p[1] + h * d0[1] / 3., 1.],
            [q[0] - h * d1[0] / 3., q[1] - h * d1[1] / 3., 1.],
            [q[0], q[1], 1.],
        ];
        // Bound each stored control against its ideal Hermite coefficient enclosure.
        let ideal = [
            [start[0].c[0], start[1].c[0]],
            [
                start[0].c[0].add(start[0].c[1].mul(I::new(h / 3., h / 3.))),
                start[1].c[0].add(start[1].c[1].mul(I::new(h / 3., h / 3.))),
            ],
            [
                end[0].c[0].sub(end[0].c[1].mul(I::new(h / 3., h / 3.))),
                end[1].c[0].sub(end[1].c[1].mul(I::new(h / 3., h / 3.))),
            ],
            [end[0].c[0], end[1].c[0]],
        ];
        let coefficient_error = controls
            .iter()
            .zip(ideal)
            .map(|(p, q)| {
                q[0].sub(I::new(p[0], p[0]))
                    .absmax()
                    .hypot(q[1].sub(I::new(p[1], p[1])).absmax())
            })
            .fold(0.0, f64::max);
        let err = up(fourth * h.powi(4) / 384.0 + coefficient_error);
        if err > TOL {
            if depth >= 30 {
                return Err("dyadic subdivision budget exhausted");
            }
            let m = (a + b) * 0.5;
            stack.push((m, b, depth + 1));
            stack.push((a, m, depth + 1));
            continue;
        }
        spans.push(RailSpan {
            interval: [a, b],
            controls,
            coefficient_error_bound_mm: up(coefficient_error),
            error_bound_mm: err,
        });
    }
    spans.sort_by(|a, b| a.interval[0].total_cmp(&b.interval[0]));
    Ok(OffsetRail {
        source_block: "s3_i_u_zbase_g2_s0_p0_z2_src0",
        side,
        radius_mm: RADIUS,
        spans,
        continuity: "C0",
    })
}

pub fn source_endpoint_offset(side: i8, end: bool) -> Result<[f64; 2], &'static str> {
    curve(if end { 1. } else { 0. }, side as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn both_sides_cover_the_whole_named_span_with_one_micron_bounds() {
        for side in [-1, 1] {
            let r = build(side).unwrap();
            assert_eq!(r.spans.first().unwrap().interval[0], 0.);
            assert_eq!(r.spans.last().unwrap().interval[1], 1.);
            assert!(r.spans.iter().all(|s| s.error_bound_mm <= TOL));
            for w in r.spans.windows(2) {
                assert_eq!(w[0].interval[1], w[1].interval[0]);
                assert_eq!(w[0].controls[3], w[1].controls[0]);
            }
        }
    }
    #[test]
    fn rejects_invalid_orientation() {
        assert!(build(0).is_err());
    }
    #[test]
    fn endpoint_rails_match_the_symbolic_speed_values() {
        let a = curve(0., 1.).unwrap();
        let b = curve(1., 1.).unwrap();
        assert!(a.iter().all(|x| x.is_finite()));
        assert_eq!(b, [75.0, 33.75]);
        assert!((a[0] - (82.4 - 0.25 / 2.0_f64.sqrt())).abs() < 1e-10);
    }
}
