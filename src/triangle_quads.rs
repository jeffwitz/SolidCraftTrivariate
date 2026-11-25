//! Exact conversion of a quadratic rational triangle to three biquadratic
//! tensor patches. Each corner owns the quadrilateral from its two incident
//! edge midpoints to the barycenter in the reference triangle.
fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}
const CONTROL: [[usize; 3]; 3] = [[0, 3, 5], [3, 1, 4], [5, 4, 2]];

/// Triangle coefficients in A,B,C,AB,BC,CA order. Planar homogeneous triples
/// (wx,wy,w) must share a gauge. Returns A,B,C patches with index u*3+v.
/// Positivity of weights is checked here; physical regularity is a separate audit.
pub fn quadrangulate(controls: [[i64; 3]; 6]) -> Result<[[[i64; 3]; 9]; 3], &'static str> {
    if controls.iter().any(|p| p[2] <= 0) {
        return Err("nonpositive triangle weight");
    }
    let mut patches = [[[0; 3]; 9]; 3];
    for (corner, patch) in patches.iter_mut().enumerate() {
        let next = (corner + 1) % 3;
        let prev = (corner + 2) % 3;
        // Bilinear barycentric restriction. All four corners use denominator six.
        let mut bary = [[0_i128; 3]; 4];
        bary[0][corner] = 6;
        bary[1][corner] = 3;
        bary[1][prev] = 3;
        bary[2][corner] = 3;
        bary[2][next] = 3;
        bary[3] = [2; 3];
        let mut net = [[0_i128; 3]; 9];
        for (i, a) in bary.iter().enumerate() {
            for (j, b) in bary.iter().enumerate() {
                let (u, v) = (i / 2 + j / 2, i % 2 + j % 2);
                // Product of Bernstein degree-one factors, cleared by a
                // common factor four instead of dividing the middle rows.
                let factor = (if u == 1 { 1 } else { 2 }) * (if v == 1 { 1 } else { 2 });
                for (s, &as_) in a.iter().enumerate() {
                    for (t, &bt) in b.iter().enumerate() {
                        let h = controls[CONTROL[s][t]];
                        for (out, x) in net[u * 3 + v].iter_mut().zip(h) {
                            *out += factor * as_ * bt * i128::from(x);
                        }
                    }
                }
            }
        }
        // Per coefficient the multiplier sum is 144, well inside i128 even
        // for an extremal i64 input. One common reduction preserves the map.
        let common = net.iter().flatten().copied().fold(0, gcd);
        for (dst, point) in patch.iter_mut().zip(net) {
            for (out, x) in dst.iter_mut().zip(point) {
                *out = i64::try_from(x / common)
                    .map_err(|_| "quadratic triangle restriction exceeds i64")?;
            }
        }
    }
    Ok(patches)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reproduces_a_weighted_triangle_and_its_parameter_restrictions() {
        let h = [
            [0, 0, 2],
            [8, 0, 4],
            [0, 6, 3],
            [3, -1, 3],
            [4, 4, 4],
            [-1, 2, 2],
        ];
        let patches = quadrangulate(h).unwrap();
        for (i, patch) in patches.iter().enumerate() {
            for u in [0., 0.2, 0.5, 1.] {
                for v in [0., 0.37, 1.] {
                    let mut l = [0.; 3];
                    l[i] = (1. - u) * (1. - v) + (u * (1. - v) + v * (1. - u)) / 2. + u * v / 3.;
                    l[(i + 1) % 3] = u * (1. - v) / 2. + u * v / 3.;
                    l[(i + 2) % 3] = v * (1. - u) / 2. + u * v / 3.;
                    let tri: [f64; 3] = std::array::from_fn(|d| {
                        (0..3)
                            .flat_map(|a| {
                                (0..3).map(move |b| l[a] * l[b] * h[CONTROL[a][b]][d] as f64)
                            })
                            .sum()
                    });
                    let basis = |x: f64| [(1. - x).powi(2), 2. * x * (1. - x), x * x];
                    let (bu, bv) = (basis(u), basis(v));
                    let quad: [f64; 3] = std::array::from_fn(|d| {
                        (0..9)
                            .map(|n| bu[n / 3] * bv[n % 3] * patch[n][d] as f64)
                            .sum()
                    });
                    for d in 0..2 {
                        assert!((tri[d] / tri[2] - quad[d] / quad[2]).abs() < 2e-14);
                    }
                }
            }
        }
    }
}
