//! Moving bore with quantized rational boundary arcs. Coordinates on straight
//! clipping lines stay exact. Fidelity to the ideal circle is an audit result.
use crate::{
    bracket::BracketBlock,
    travelling_bore::{build_travel, isqrt},
};

const TICKS_PER_MM: i128 = 1_000_000_000_000;
const TICKS_PER_UM: i128 = TICKS_PER_MM / 1_000;

/// Return the nearest femtometre to the positive ordinate of an ideal cut.
pub(crate) fn circle_height_fm(distance_um: i128) -> i128 {
    let square = (2_000_i128.pow(2) - distance_um.pow(2)) * TICKS_PER_UM.pow(2);
    let lower = isqrt(square);
    lower + i128::from(square - lower * lower > lower)
}

/// Triangle controls in a single scale 2Q. Straight sides use the exact
/// homogeneous midpoint. For a circular side between unit directions a,b,
/// the central homogeneous numerator relative to the center is
/// R*(a+b)/|a+b|, and its weight is cos(delta_angle/2)=|a+b|/2.
pub(crate) fn triangle_net(
    points: [[i128; 2]; 3],
    circular: [bool; 3],
    center_um: i128,
) -> Result<[[i64; 3]; 6], &'static str> {
    let mut net = [[0_i128; 3]; 6];
    for (dst, [x, y]) in net.iter_mut().zip(points) {
        *dst = [2 * x, 2 * y, 2 * TICKS_PER_MM];
    }
    let center = center_um * TICKS_PER_UM;
    let radius = 2 * TICKS_PER_MM;
    for (edge, (a, b)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {
        if circular[a] && circular[b] {
            let unit = |i: usize| {
                let x = (points[i][0] - center) as f64;
                let y = points[i][1] as f64;
                let length = x.hypot(y);
                [x / length, y / length]
            };
            let (va, vb) = (unit(a), unit(b));
            let (x, y) = (va[0] + vb[0], va[1] + vb[1]);
            let length = x.hypot(y);
            let w = length / 2.;
            let round = |v: f64| -> Result<i128, &'static str> {
                if v.is_finite() && v.abs() < (1_u64 << 53) as f64 {
                    Ok(v.round() as i128)
                } else {
                    Err("bounded arc coefficient outside exact rounding range")
                }
            };
            // A cardinal endpoint has an exact horizontal or vertical tangent.
            // Preserve it when quantizing the middle control, including at
            // centers with a fractional millimetre x-coordinate.
            let x_cardinal = [a, b].into_iter().find(|&i| points[i][1] == 0);
            let y_cardinal = [a, b].into_iter().find(|&i| points[i][0] == center);
            let quantum = if let Some(i) = x_cardinal {
                let mut a = points[i][0].abs();
                let mut b = TICKS_PER_MM;
                while b != 0 {
                    (a, b) = (b, a % b);
                }
                TICKS_PER_MM / a
            } else {
                1
            };
            let iw = quantum * round(TICKS_PER_MM as f64 * w / quantum as f64)?;
            let mut wx = round(
                iw as f64 * (center as f64 + radius as f64 * x / (length * w))
                    / TICKS_PER_MM as f64,
            )?;
            let mut wy = round(iw as f64 * radius as f64 * y / (length * w) / TICKS_PER_MM as f64)?;
            if let Some(i) = x_cardinal {
                wx = points[i][0] * iw / TICKS_PER_MM;
            }
            if let Some(i) = y_cardinal {
                wy = points[i][1] * iw / TICKS_PER_MM;
            }
            net[3 + edge] = [2 * wx, 2 * wy, 2 * iw];
        } else {
            net[3 + edge] = std::array::from_fn(|axis| (net[a][axis] + net[b][axis]) / 2);
        }
    }
    let mut controls = [[0_i64; 3]; 6];
    for (output, input) in controls.iter_mut().flatten().zip(net.into_iter().flatten()) {
        *output = i64::try_from(input).map_err(|_| "bounded triangle outside i64")?;
    }
    Ok(controls)
}

/// Build quantized rational geometry at an integer-micrometre center.
pub fn travelling_bore_bounded_blocks(center_um: i32) -> Result<Vec<BracketBlock>, &'static str> {
    build_travel(center_um, true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cuts_are_rounded_by_integer_comparison() {
        for distance in 0..=2000 {
            let y = circle_height_fm(distance);
            let exact4 = 4 * (2000_i128.pow(2) - distance.pow(2)) * TICKS_PER_UM.pow(2);
            assert!((2 * y - 1).max(0).pow(2) <= exact4 && exact4 <= (2 * y + 1).pow(2));
        }
    }

    #[test]
    fn supports_irrational_and_one_micrometre_approaches() {
        for c in [
            -5000, -4001, -3999, -2001, -1999, -100, -1, 1, 7999, 8001, 9999, 10001, 11900, 11999,
            12001, 14000,
        ] {
            let blocks = travelling_bore_bounded_blocks(c).unwrap_or_else(|e| panic!("{c}: {e}"));
            let incidence = crate::incidence::enumerate_faces(&blocks).unwrap();
            assert!(!incidence.interfaces.is_empty());
        }
        for c in [-5001, -4000, 0, 8000, 12000, 14001] {
            assert!(travelling_bore_bounded_blocks(c).is_err());
        }
    }
}
