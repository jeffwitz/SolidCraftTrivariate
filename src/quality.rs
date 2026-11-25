//! Pointwise shape diagnostics and conservative full-block distortion bounds.

use crate::bracket::BracketBlock;

#[derive(Debug, Clone, Copy)]
pub struct QualitySample {
    /// Descending physical stretches per unit parameter.
    pub singular_values: [f64; 3],
    pub condition_number: f64,
    pub scaled_determinant: f64,
}

fn norm(x: [f64; 3]) -> f64 {
    x[0].hypot(x[1]).hypot(x[2])
}
fn dot(x: [f64; 3], y: [f64; 3]) -> f64 {
    x[0] * y[0] + x[1] * y[1] + x[2] * y[2]
}

/// One-sided Jacobi works on the columns themselves, avoiding the loss of a
/// small singular value from explicitly forming the Gram matrix.
pub fn at_jacobian(columns: [[f64; 3]; 3]) -> Result<QualitySample, &'static str> {
    let lengths = columns.map(norm);
    if lengths.iter().any(|x| !x.is_finite() || *x <= 0.) {
        return Err("invalid Jacobian column norm");
    }
    let [a, b, c]: [[f64; 3]; 3] = std::array::from_fn(|i| columns[i].map(|x| x / lengths[i]));
    let det = dot(
        a,
        [
            b[1] * c[2] - b[2] * c[1],
            b[2] * c[0] - b[0] * c[2],
            b[0] * c[1] - b[1] * c[0],
        ],
    );
    let mut u = columns;
    for _ in 0..32 {
        let mut changed = false;
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            let p = u[i];
            let q = u[j];
            let np = norm(p);
            let nq = norm(q);
            if np == 0. || nq == 0. {
                return Err("singular Jacobian");
            }
            let corr = dot(p.map(|v| v / np), q.map(|v| v / nq));
            if corr.abs() <= 8. * f64::EPSILON {
                continue;
            }
            let scale = np.max(nq);
            let l = np / scale;
            let r = nq / scale;
            let gamma = corr * l * r;
            if gamma == 0. {
                return Err("Jacobian singular-value range indeterminate");
            }
            let delta = (r * r - l * l) * 0.5;
            let t = if delta == 0. {
                1.
            } else {
                gamma * delta.signum() / (delta.abs() + delta.hypot(gamma))
            };
            let cosine = 1. / 1f64.hypot(t);
            let sine = t * cosine;
            u[i] = std::array::from_fn(|k| cosine * p[k] - sine * q[k]);
            u[j] = std::array::from_fn(|k| sine * p[k] + cosine * q[k]);
            changed = true;
        }
        if !changed {
            let mut s = u.map(norm);
            s.sort_by(|a, b| b.total_cmp(a));
            let cond = s[0] / s[2];
            if !cond.is_finite() {
                return Err("nonfinite Jacobian condition");
            }
            return Ok(QualitySample {
                singular_values: s,
                condition_number: cond,
                scaled_determinant: det.clamp(-1., 1.),
            });
        }
    }
    Err("Jacobian singular-value iteration did not converge")
}

pub struct QualityBounds {
    pub column_norm_upper: [f64; 3],
    pub maximum_stretch_upper: f64,
    pub minimum_stretch_lower: f64,
    pub condition_upper: f64,
    pub scaled_determinant_lower: f64,
}

fn up_add(a: f64, b: f64) -> f64 {
    (a + b).next_up()
}
fn up_mul(a: f64, b: f64) -> f64 {
    (a * b).next_up()
}
fn up_div(a: f64, b: f64) -> f64 {
    (a / b).next_up()
}
fn upper_magnitude(x: i128) -> f64 {
    (x.unsigned_abs() as f64).next_up()
}
fn vector_bound(x: [i128; 3]) -> f64 {
    x.into_iter().map(upper_magnitude).fold(0., up_add)
}

/// Use a positive lower bound on the *same block's* closed-domain determinant
/// from an independent certificate. The column bounds follow from differentiating
/// X = H/W and applying the Bernstein convex hull property to its two factors.
pub fn bound_block(
    block: &BracketBlock,
    determinant_lower: f64,
) -> Result<QualityBounds, &'static str> {
    if !matches!(block.degrees, [1, 1, 1] | [2, 1, 1])
        || !determinant_lower.is_finite()
        || determinant_lower <= 0.
    {
        return Err("unsupported quality-bound input");
    }
    let dims = block.degrees.map(|d| d + 1);
    if block.controls.len() != dims.iter().product() || block.controls.iter().any(|p| p[3] <= 0) {
        return Err("invalid rational block for quality bound");
    }
    let origin = block.controls[0].map(i128::from);
    let w0_lower = (origin[3] as f64).next_down();
    let w_min = block.controls.iter().map(|p| p[3]).min().unwrap();
    let mut shifted = Vec::with_capacity(block.controls.len());
    let mut radius: f64 = 0.;
    for p in &block.controls {
        let p = p.map(i128::from);
        // w0*(H_i-P0*w_i) in integer arithmetic, so translations cancel exactly.
        let h = std::array::from_fn(|i| p[i] * origin[3] - origin[i] * p[3]);
        let coord = up_div(up_div(vector_bound(h), w0_lower), (p[3] as f64).next_down());
        radius = radius.max(coord);
        shifted.push(h);
    }
    let index = |v: [usize; 3]| (v[0] * dims[1] + v[1]) * dims[2] + v[2];
    let mut bounds = [0.; 3];
    for (axis, m) in bounds.iter_mut().enumerate() {
        let mut delta_h: f64 = 0.;
        let mut delta_w: f64 = 0.;
        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    let mut q = [i, j, k];
                    if q[axis] + 1 == dims[axis] {
                        continue;
                    }
                    let a = index(q);
                    q[axis] += 1;
                    let b = index(q);
                    let mut diff = [0; 3];
                    for (c, value) in diff.iter_mut().enumerate() {
                        *value = shifted[b][c]
                            .checked_sub(shifted[a][c])
                            .ok_or("quality-bound integer range")?;
                    }
                    delta_h = delta_h.max(up_div(vector_bound(diff), w0_lower));
                    delta_w = delta_w.max(upper_magnitude(
                        i128::from(block.controls[b][3]) - i128::from(block.controls[a][3]),
                    ));
                }
            }
        }
        *m = up_div(
            up_mul(
                block.degrees[axis] as f64,
                up_add(delta_h, up_mul(radius, delta_w)),
            ),
            (w_min as f64).next_down(),
        );
    }
    let frob2 = bounds.iter().map(|&b| up_mul(b, b)).fold(0., up_add);
    let maximum = frob2.sqrt().next_up();
    let minimum = ((2. * determinant_lower) / frob2).next_down();
    let condition = up_div(maximum, minimum);
    let scaled = (determinant_lower / bounds.iter().copied().fold(1., up_mul))
        .next_down()
        .min(1.);
    if !condition.is_finite()
        || minimum <= 0.
        || !minimum.is_finite()
        || !scaled.is_finite()
        || scaled <= 0.
    {
        return Err("quality-bound numerical range");
    }
    Ok(QualityBounds {
        column_norm_upper: bounds,
        maximum_stretch_upper: maximum,
        minimum_stretch_lower: minimum,
        condition_upper: condition,
        scaled_determinant_lower: scaled,
    })
}
