//! Apply a common rational rigid frame to each homogeneous control net.
use crate::face_profile_extrusion::SectionPatch;
const SCALE: i128 = 1_000_000_000;
fn coefficient(x: f64) -> Result<i128, &'static str> {
    if !x.is_finite() || x.abs() > 1e7 {
        return Err("face frame exceeds reconstruction range");
    }
    Ok((x * SCALE as f64).round() as i128)
}
fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}
pub(crate) fn lift(
    p: &SectionPatch,
    depth: f64,
    frame: [[f64; 3]; 4],
) -> Result<Vec<[i64; 4]>, &'static str> {
    let mut net = Vec::with_capacity(p.controls.len() * 2);
    for &[x, y, w] in &p.controls {
        for z in [0., depth] {
            let mut h = [0i128; 4];
            for (d, c) in h.iter_mut().enumerate().take(3) {
                *c = coefficient(frame[1][d])? * i128::from(x)
                    + coefficient(frame[2][d])? * i128::from(y)
                    + coefficient(frame[0][d] + z * frame[3][d])? * i128::from(w);
            }
            h[3] = SCALE * i128::from(w);
            net.push(h);
        }
    }
    let factor = net.iter().flatten().copied().fold(0, gcd).max(1);
    net.into_iter()
        .map(|p| {
            let mut control = [0; 4];
            for (i, n) in p.iter().enumerate() {
                control[i] = i64::try_from(n / factor)
                    .map_err(|_| "face frame exceeds exact homogeneous range")?
            }
            Ok(control)
        })
        .collect()
}
