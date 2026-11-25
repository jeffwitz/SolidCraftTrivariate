//! Exact rational material charts for intersecting X/Z through-bores of equal radius.
//!
//! In each octant the two cylinder inequalities reduce to
//! min(|x|,|z|)^2 + y^2 >= r^2. Sweep the meridional section in the
//! dominant coordinate from the symmetry plane to the outside of the stock.
//! Coordinates in this description are relative to the axis intersection.

use crate::bracket::BracketBlock;

struct Octant {
    center: [i64; 3],
    sign: [i64; 3],
    extent: [i64; 3],
    radius: i64,
    guard: i64,
}

impl Octant {
    fn world(&self, local: [i64; 3], weight: i64) -> [i64; 4] {
        let p: [i64; 3] =
            std::array::from_fn(|i| self.center[i] * weight + self.sign[i] * local[i]);
        [p[0], p[1], p[2], weight]
    }

    fn sector(&self, dominant: usize, arc: i64, id: String) -> BracketBlock {
        let mut controls = Vec::with_capacity(12);
        let handedness = if dominant == 0 { 1 } else { -1 };
        let forward = handedness * self.sign.iter().product::<i64>() > 0;
        // Stereographic q runs through [0, 1/2] then [1/2, 1].
        // Use the blossom at (p/2,t/2) to form exact quadratic controls.
        for (i, (p, t)) in [(arc, arc), (arc, arc + 1), (arc + 1, arc + 1)]
            .into_iter()
            .enumerate()
        {
            let i = i as i64;
            let circle = [
                self.radius * (4 - p * t),
                self.radius * 2 * (p + t),
                4 + p * t,
            ];
            // The outer edge is linear with constant weight. Its midpoint
            // control preserves the degree-two trace on the neighboring box.
            let outer = if arc == 0 {
                [4 * self.guard, 2 * i * self.extent[1], 4]
            } else {
                [2 * (2 - i) * self.guard, 4 * self.extent[1], 4]
            };
            for [minor, height, weight] in [circle, outer] {
                for upper in [false, true] {
                    let mut point = [minor, height, minor];
                    if upper == forward {
                        point[dominant] = self.extent[dominant] * weight;
                    }
                    controls.push(self.world(point, weight));
                }
            }
        }
        BracketBlock {
            id,
            degrees: [2, 1, 1],
            controls,
        }
    }

    fn corner(&self, id: String) -> BracketBlock {
        let intervals: [[i64; 2]; 3] = std::array::from_fn(|i| {
            let start = if i == 1 { 0 } else { self.guard };
            let mut ends = [start, self.extent[i]].map(|v| self.center[i] + self.sign[i] * v);
            ends.sort();
            ends
        });
        let mut controls = Vec::with_capacity(8);
        for x in intervals[0] {
            for y in intervals[1] {
                for z in intervals[2] {
                    controls.push([x, y, z, 1]);
                }
            }
        }
        BracketBlock {
            id,
            degrees: [1; 3],
            controls,
        }
    }
}

/// Partition the material of an axis-aligned stock after both through-cuts.
/// Each sign octant has four swept circular sectors and one outer brick.
/// The guard is an internal partition coordinate; changing it preserves the part.
pub fn steinmetz_box(
    bounds: [[i32; 2]; 3],
    center: [i32; 3],
    radius: i32,
    guard: i32,
    id: &str,
) -> Result<Vec<BracketBlock>, &'static str> {
    if radius <= 0 || guard <= radius || id.is_empty() {
        return Err("invalid_steinmetz_parameters");
    }
    let center = center.map(i64::from);
    let radius = i64::from(radius);
    let guard = i64::from(guard);
    let mut blocks = Vec::with_capacity(40);
    for octant in 0..8 {
        let sign = std::array::from_fn(|i| if octant & (1 << i) == 0 { -1 } else { 1 });
        let extent = std::array::from_fn(|i| {
            let end = bounds[i][usize::from(sign[i] > 0)];
            sign[i] * (i64::from(end) - center[i])
        });
        if extent[0] <= guard || extent[2] <= guard || extent[1] <= radius {
            return Err("steinmetz_partition_clearance");
        }
        let chart = Octant {
            center,
            sign,
            extent,
            radius,
            guard,
        };
        for dominant in [0, 2] {
            for arc in [0, 1] {
                blocks.push(chart.sector(
                    dominant,
                    arc,
                    format!("{id}_o{octant}_axis{dominant}_arc{arc}"),
                ));
            }
        }
        blocks.push(chart.corner(format!("{id}_o{octant}_corner")));
    }
    Ok(blocks)
}
