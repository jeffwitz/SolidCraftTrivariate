//! Exact rational annular charts for a transverse bore in a square guard.
use crate::bracket::BracketBlock;
use std::collections::BTreeSet;

fn add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("sketch_bore_integer_overflow")
}

fn mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("sketch_bore_integer_overflow")
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn block(id: String, controls: Vec<[i128; 4]>) -> Result<BracketBlock, &'static str> {
    if controls.len() != 12 || controls.iter().any(|p| p[3] <= 0) {
        return Err("sketch_bore_invalid_net");
    }
    let common = controls
        .iter()
        .flatten()
        .fold(0_u128, |g, x| gcd(g, x.unsigned_abs()));
    if common == 0 {
        return Err("sketch_bore_zero_net");
    }
    let common = i128::try_from(common).map_err(|_| "sketch_bore_integer_overflow")?;
    let controls = controls
        .into_iter()
        .map(|p| {
            let p = p.map(|v| v / common);
            Ok([
                i64::try_from(p[0]).map_err(|_| "sketch_bore_control_exceeds_i64")?,
                i64::try_from(p[1]).map_err(|_| "sketch_bore_control_exceeds_i64")?,
                i64::try_from(p[2]).map_err(|_| "sketch_bore_control_exceeds_i64")?,
                i64::try_from(p[3]).map_err(|_| "sketch_bore_control_exceeds_i64")?,
            ])
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(BracketBlock {
        id,
        degrees: [2, 1, 1],
        controls,
    })
}

fn restrict(p: [[i128; 3]; 3], half: usize) -> [[i128; 3]; 3] {
    let mut result = [[0_i128; 3]; 3];
    for coordinate in 0..3 {
        let [a, b, c] = [p[0][coordinate], p[1][coordinate], p[2][coordinate]];
        let values = if half == 0 {
            [4 * a, 2 * (a + b), a + 2 * b + c]
        } else {
            [a + 2 * b + c, 2 * (b + c), 4 * c]
        };
        for i in 0..3 {
            result[i][coordinate] = values[i];
        }
    }
    result
}

fn source_index(i: usize, quadrant: usize) -> usize {
    if quadrant % 2 == 1 {
        2 - i
    } else {
        i
    }
}

fn rotate(point: [i128; 3], quadrant: usize) -> [i128; 3] {
    let [x, y, w] = point;
    match quadrant {
        0 => [x, y, w],
        1 => [-x, y, w],
        2 => [-x, -y, w],
        _ => [x, -y, w],
    }
}

fn outer_segment(half_width: i128, quadrant: usize, half: usize) -> [[i128; 2]; 3] {
    let g = half_width;
    let (start, corner, end) = match quadrant {
        0 => ([g, 0], [g, g], [0, g]),
        1 => ([0, g], [-g, g], [-g, 0]),
        2 => ([-g, 0], [-g, -g], [0, -g]),
        _ => ([0, -g], [g, -g], [g, 0]),
    };
    let (a, b) = if half == 0 {
        (start, corner)
    } else {
        (corner, end)
    };
    [
        [2 * a[0], 2 * a[1]],
        [a[0] + b[0], a[1] + b[1]],
        [2 * b[0], 2 * b[1]],
    ]
}

/// Build eight exact half-quadrant charts around a bore in an axis-aligned square guard.
/// `interval` gives absolute world coordinates along `axis`; transverse coordinates
/// are centered at `center[(axis+1)%3]` and `center[(axis+2)%3]`.
pub fn bore_charts(
    axis: usize,
    center: [i32; 3],
    interval: [i32; 2],
    half_width: i32,
    radius: i32,
    id: &str,
) -> Result<Vec<BracketBlock>, &'static str> {
    if axis > 2 || interval[0] >= interval[1] || radius <= 0 || half_width <= radius {
        return Err("sketch_bore_invalid_domain");
    }
    if id.is_empty() || id.len() > 24 || !id.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_')
    {
        return Err("sketch_bore_invalid_id");
    }
    let axes = [(axis + 1) % 3, (axis + 2) % 3];
    let c0 = i128::from(center[axes[0]]);
    let c1 = i128::from(center[axes[1]]);
    let scale = 2_i128;
    let r = i128::from(radius);
    let arc = [[2 * r, 0, 2], [2 * r, 2 * r, 2], [0, 4 * r, 4]];
    let names = ["east", "north", "west", "south"];
    let mut blocks = Vec::with_capacity(8);
    let mut ids = BTreeSet::new();
    for quadrant in 0..4 {
        let mut full = [[0_i128; 3]; 3];
        for u in 0..3 {
            full[u] = rotate(arc[source_index(u, quadrant)], quadrant);
        }
        for half in 0..2 {
            let inner = restrict(full, half);
            let outer = outer_segment(i128::from(half_width), quadrant, half);
            let mut controls = Vec::with_capacity(12);
            for u in 0..3 {
                for v in 0..2 {
                    for axial in interval {
                        let (x, y, weight) = if v == 0 {
                            (outer[u][0], outer[u][1], scale)
                        } else {
                            (
                                mul(inner[u][0], scale)?,
                                mul(inner[u][1], scale)?,
                                mul(inner[u][2], scale)?,
                            )
                        };
                        let mut world = [0_i128; 3];
                        world[axes[0]] = add(mul(c0, weight)?, x)?;
                        world[axes[1]] = add(mul(c1, weight)?, y)?;
                        world[axis] = mul(i128::from(axial), weight)?;
                        controls.push([world[0], world[1], world[2], weight]);
                    }
                }
            }
            let block_id = format!(
                "bore_{id}_sector_{}_half_{half}_axis_{axis}_{}_{}",
                names[quadrant], interval[0], interval[1]
            );
            if !ids.insert(block_id.clone()) {
                return Err("sketch_bore_duplicate_block_id");
            }
            blocks.push(block(block_id, controls)?);
        }
    }
    Ok(blocks)
}
