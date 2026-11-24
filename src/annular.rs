//! Fixed two-block rational quarter-annulus construction.

/// One homogeneous tensor-product Bézier block in millimetres.
#[derive(Clone, Debug, PartialEq)]
pub struct AnnularBlock {
    pub id: &'static str,
    pub degrees: [usize; 3],
    /// Homogeneous `(wx, wy, wz, w)` controls; the `w` index varies fastest.
    pub controls: Vec<[f64; 4]>,
}

/// Two quarter-annulus blocks joined across their horizontal interface.
#[derive(Clone, Debug, PartialEq)]
pub struct AnnularPair {
    pub heights_mm: [u32; 2],
    pub blocks: Vec<AnnularBlock>,
    pub interface: [&'static str; 3],
}

const MAX_HEIGHT_MM: u32 = 1000;
const DEGREES: [usize; 3] = [2, 1, 1];

/// Construct two stacked rational quarter-annuli with radii 1 and 2 mm.
///
/// The angular direction increases in `u`, the radial direction decreases
/// from radius 2 to radius 1 in `v`, and altitude increases in `w`. The two
/// supplied positive heights are bounded to 1000 mm each. Controls are
/// flattened in `(u, v, w)` order, with `w` varying fastest.
pub fn build_annular(heights_mm: [u32; 2]) -> Result<AnnularPair, &'static str> {
    if heights_mm
        .iter()
        .any(|height| *height == 0 || *height > MAX_HEIGHT_MM)
    {
        return Err("heights must be positive and at most 1000 mm");
    }

    let arc_weights = [1.0_f64, 0.5_f64.sqrt(), 1.0_f64];
    let mut blocks = Vec::with_capacity(2);
    let mut z_start = 0.0_f64;

    for (index, height_mm) in heights_mm.iter().copied().enumerate() {
        let z_end = z_start + f64::from(height_mm);
        let mut controls = Vec::with_capacity(12);

        for (angle_index, weight) in arc_weights.iter().copied().enumerate() {
            for radius in [2.0_f64, 1.0_f64] {
                let (x, y) = match angle_index {
                    0 => (radius, 0.0),
                    1 => (radius, radius),
                    _ => (0.0, radius),
                };
                for z in [z_start, z_end] {
                    controls.push([weight * x, weight * y, weight * z, weight]);
                }
            }
        }

        blocks.push(AnnularBlock {
            id: if index == 0 { "lower" } else { "upper" },
            degrees: DEGREES,
            controls,
        });
        z_start = z_end;
    }

    Ok(AnnularPair {
        heights_mm,
        blocks,
        interface: ["lower.w1", "upper.w0", "identity_uv"],
    })
}
