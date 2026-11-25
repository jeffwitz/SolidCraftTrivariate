//! Four-block rational chart for a rectangular prism with an axial round bore.
//!
//! Each chart spans one side of the rectangle at `v = 0` and one quarter of
//! the circular bore at `v = 1`. The angular join points are the exact
//! rational points `(3, 4) / 5` and its quarter turns. The alternating
//! endpoint weights make the radial traces identical at the four joins.

use crate::bracket::BracketBlock;
use std::convert::TryFrom;

type HPoint = [i128; 3]; // Two homogeneous section coordinates and a weight.

fn checked_i64(value: i128) -> Result<i64, &'static str> {
    i64::try_from(value).map_err(|_| "compact wing control exceeds i64")
}

fn side_corners(bounds: [[i32; 2]; 3], plane: [usize; 2], side: usize) -> ([i128; 2], [i128; 2]) {
    let [a, b] = plane;
    let (lo_a, hi_a) = (i128::from(bounds[a][0]), i128::from(bounds[a][1]));
    let (lo_b, hi_b) = (i128::from(bounds[b][0]), i128::from(bounds[b][1]));
    // For the counter-clockwise ring the arc sectors are north, west, south,
    // east; their endpoints follow the outer perimeter counter-clockwise.
    let ccw = [
        ([hi_a, hi_b], [lo_a, hi_b]),
        ([lo_a, hi_b], [lo_a, lo_b]),
        ([lo_a, lo_b], [hi_a, lo_b]),
        ([hi_a, lo_b], [hi_a, hi_b]),
    ];
    ccw[side]
}

/// Build the four exact rational charts of a rectangular prism with one
/// circular through-bore. `axis` selects the extrusion coordinate; the other
/// two coordinates form the circular section. The bore must lie strictly
/// inside that rectangle. The returned maps have degrees `[2, 1, 1]`, with
/// `u` circumferential, `v=0` on the rectangle and `v=1` on the circle, and
/// `w` along the extrusion axis.
pub(crate) fn build(
    axis: usize,
    bounds: [[i32; 2]; 3],
    center: [i32; 3],
    radius: i32,
    id: &str,
) -> Result<Vec<BracketBlock>, &'static str> {
    if axis > 2 || radius <= 0 || id.is_empty() {
        return Err("compact_wing_invalid_parameters");
    }
    if bounds.iter().any(|pair| pair[0] >= pair[1]) {
        return Err("compact_wing_invalid_bounds");
    }
    let plane = [(axis + 1) % 3, (axis + 2) % 3];
    let (a, b) = (plane[0], plane[1]);
    let (ca, cb, r) = (
        i128::from(center[a]),
        i128::from(center[b]),
        i128::from(radius),
    );
    if ca - r <= i128::from(bounds[a][0])
        || ca + r >= i128::from(bounds[a][1])
        || cb - r <= i128::from(bounds[b][0])
        || cb + r >= i128::from(bounds[b][1])
    {
        return Err("compact_wing_bore_must_be_strictly_inside_bounds");
    }

    // The cyclic axis order makes (plane[0], plane[1], axis) positive, so the
    // same counter-clockwise section orientation works for every extrusion.
    let radial = [
        [3_i128 * r, 4_i128 * r],
        [-4_i128 * r, 3_i128 * r],
        [-3_i128 * r, -4_i128 * r],
        [4_i128 * r, -3_i128 * r],
    ];
    let axial = [i128::from(bounds[axis][0]), i128::from(bounds[axis][1])];
    let mut blocks = Vec::with_capacity(4);
    for sector in 0..4 {
        let (outer_start, outer_end) = side_corners(bounds, plane, sector);
        let start = radial[sector];
        let end = radial[(sector + 1) % 4];

        // A quarter-circle conic has tangent-intersection control C=P+Q.
        // With endpoint weights (10,20) or (20,10), its middle weight is 10.
        // This avoids irrational sqrt(2) entirely.
        let endpoint_weights = if sector % 2 == 0 {
            [10_i128, 20_i128]
        } else {
            [20, 10]
        };
        let inner = [
            [
                endpoint_weights[0] * ca + endpoint_weights[0] * start[0] / 5,
                endpoint_weights[0] * cb + endpoint_weights[0] * start[1] / 5,
                endpoint_weights[0],
            ],
            [
                10 * ca + 2 * (start[0] + end[0]),
                10 * cb + 2 * (start[1] + end[1]),
                10,
            ],
            [
                endpoint_weights[1] * ca + endpoint_weights[1] * end[0] / 5,
                endpoint_weights[1] * cb + endpoint_weights[1] * end[1] / 5,
                endpoint_weights[1],
            ],
        ];

        let outer = [
            [10 * outer_start[0], 10 * outer_start[1], 10],
            [
                5 * (outer_start[0] + outer_end[0]),
                5 * (outer_start[1] + outer_end[1]),
                10,
            ],
            [10 * outer_end[0], 10 * outer_end[1], 10],
        ];

        let mut controls = Vec::with_capacity(12);
        for u in 0..3 {
            for v in 0..2 {
                let h: HPoint = if v == 0 { outer[u] } else { inner[u] };
                for z in axial {
                    let mut xyz_w = [0_i128; 4];
                    xyz_w[a] = h[0];
                    xyz_w[b] = h[1];
                    xyz_w[axis] = h[2] * z;
                    xyz_w[3] = h[2];
                    controls.push([
                        checked_i64(xyz_w[0])?,
                        checked_i64(xyz_w[1])?,
                        checked_i64(xyz_w[2])?,
                        checked_i64(xyz_w[3])?,
                    ]);
                }
            }
        }
        let side = ["north", "west", "south", "east"][sector];
        blocks.push(BracketBlock {
            id: format!("{id}_{side}"),
            degrees: [2, 1, 1],
            controls,
        });
    }
    Ok(blocks)
}
