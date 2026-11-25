//! Exact circular through-bores in a polygonal sheared prism.
use crate::{
    bracket::BracketBlock,
    perforated_prism::{build_perforated_prism, LoopProfile, PerforatedPrism},
    polygon_prism::Vertex,
};
use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CircularPolygonBore {
    pub id: String,
    pub center_mm: [i32; 2],
    pub radius_mm: i32,
    pub guard_mm: i32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CircularPolygonPrism {
    pub outer: LoopProfile,
    pub bores: Vec<CircularPolygonBore>,
    pub height_mm: i32,
    pub top_shift_mm: [i32; 2],
}

const RING_GAUGE: i128 = 8;

fn valid_bore_id(id: &str) -> bool {
    !id.is_empty() && id.len() <= 18 && id.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_')
}

fn checked_control(value: i128) -> Result<i64, &'static str> {
    i64::try_from(value).map_err(|_| "circular polygon control exceeds i64")
}

fn transformed_relative(point: [i128; 3], quadrant: usize) -> [i128; 3] {
    let [x, y, w] = point;
    match quadrant {
        0 => [x, y, w],
        1 => [-x, y, w],
        2 => [-x, -y, w],
        _ => [x, -y, w],
    }
}

fn source_index(i: usize, quadrant: usize) -> usize {
    if quadrant % 2 == 1 {
        2 - i
    } else {
        i
    }
}

fn quarter_net(
    bore: &CircularPolygonBore,
    height: i32,
    shift: [i32; 2],
    quadrant: usize,
    half: usize,
) -> Result<BracketBlock, &'static str> {
    if quadrant > 3 || half > 1 {
        return Err("invalid circular polygon sector");
    }
    // The uncut sector is represented at scale two so the diamond-edge
    // midpoint is integral. All other coordinates and weights are integral.
    let r = i128::from(bore.radius_mm);
    let g = i128::from(bore.guard_mm);
    let center = bore.center_mm.map(i128::from);
    let outer = [[2 * g, 0, 2], [g, g, 2], [0, 2 * g, 2]];
    let circle = [[2 * r, 0, 2], [2 * r, 2 * r, 2], [0, 4 * r, 4]];
    let mut base = [[[0_i128; 4]; 2]; 3]; // [u][v][x,y,z,w], before extrusion layer
    for (i, row) in base.iter_mut().enumerate() {
        let si = source_index(i, quadrant);
        for (v, source) in [(0, outer[si]), (1, circle[si])] {
            let [rx, ry, w] = transformed_relative(source, quadrant);
            row[v] = [center[0] * w + rx, center[1] * w + ry, 0, w];
        }
    }
    // Add the two extrusion layers and their common in-plane translation.
    let mut layered = [[[[0_i128; 4]; 2]; 2]; 3];
    for (i, row) in layered.iter_mut().enumerate() {
        for (v, layers) in row.iter_mut().enumerate() {
            let [x, y, _, w] = base[i][v];
            for (layer, z) in [0_i128, i128::from(height)].into_iter().enumerate() {
                layers[layer] = [
                    x + i128::from(shift[0]) * w * i128::from(layer > 0),
                    y + i128::from(shift[1]) * w * i128::from(layer > 0),
                    z * w,
                    w,
                ];
            }
        }
    }
    // The complete profile coefficients are 2*H. Restrict exactly to either
    // half interval and serialize all homogeneous controls at gauge eight.
    let mut controls = Vec::with_capacity(12);
    let restricted: [[[[i128; 4]; 2]; 2]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|v| {
            std::array::from_fn(|w| {
                let b0 = layered[0][v][w];
                let b1 = layered[1][v][w];
                let b2 = layered[2][v][w];
                match (half, i) {
                    (0, 0) => b0.map(|x| x * 4),
                    (0, 1) => std::array::from_fn(|k| 2 * (b0[k] + b1[k])),
                    (0, 2) => std::array::from_fn(|k| b0[k] + 2 * b1[k] + b2[k]),
                    (1, 0) => std::array::from_fn(|k| b0[k] + 2 * b1[k] + b2[k]),
                    (1, 1) => std::array::from_fn(|k| 2 * (b1[k] + b2[k])),
                    (1, 2) => b2.map(|x| x * 4),
                    _ => unreachable!("validated circular polygon sector indices"),
                }
            })
        })
    });
    for u in &restricted {
        for v in u {
            for numerator in v {
                controls.push([
                    checked_control(numerator[0])?,
                    checked_control(numerator[1])?,
                    checked_control(numerator[2])?,
                    checked_control(numerator[3])?,
                ]);
            }
        }
    }
    // The construction above uses the base scale two and the exact split
    // factors. RING_GAUGE documents the physical homogeneous scale.
    debug_assert_eq!(RING_GAUGE, 8);
    if controls.iter().any(|p| p[3] <= 0) {
        return Err("nonpositive circular polygon weight");
    }
    let sector = ["east", "north", "west", "south"][quadrant];
    Ok(BracketBlock {
        id: format!(
            "bore_{}_sector_{}_half_{}_z0_{}",
            bore.id, sector, half, height
        ),
        degrees: [2, 1, 1],
        controls,
    })
}

/// Construct an audited-domain candidate made of a G18 exterior and exact
/// rational annuli around one to four circular through-bores.
pub fn build_circular_polygon_prism(
    recipe: &CircularPolygonPrism,
) -> Result<Vec<BracketBlock>, &'static str> {
    if !(1..=4).contains(&recipe.bores.len())
        || recipe.outer.vertices.len() + 4 * recipe.bores.len() > 32
    {
        return Err("circular polygon bore count or vertex budget out of range");
    }
    let mut bore_ids = BTreeSet::new();
    let mut all_ids = BTreeSet::new();
    if recipe.outer.id.is_empty() || !all_ids.insert(recipe.outer.id.clone()) {
        return Err("invalid circular polygon outer loop id");
    }
    for vertex in &recipe.outer.vertices {
        if vertex.id.is_empty() || !all_ids.insert(vertex.id.clone()) {
            return Err("duplicate circular polygon vertex identity");
        }
    }
    let mut bores = recipe.bores.iter().collect::<Vec<_>>();
    bores.sort_by(|a, b| a.id.cmp(&b.id));
    let mut guard_loops = Vec::with_capacity(bores.len());
    for bore in &bores {
        if !valid_bore_id(&bore.id) || !bore_ids.insert(bore.id.clone()) {
            return Err("invalid or duplicate circular polygon bore id");
        }
        if !all_ids.insert(bore.id.clone()) {
            return Err("circular polygon bore id collides with a profile identity");
        }
        if bore.radius_mm < 1 || i64::from(bore.guard_mm) < 2 * i64::from(bore.radius_mm) {
            return Err("circular polygon guard must be at least twice the positive radius");
        }
        let loop_id = format!("guard_{}", bore.id);
        if !all_ids.insert(loop_id.clone()) {
            return Err("circular polygon guard loop id collision");
        }
        let c = bore.center_mm.map(i64::from);
        let g = i64::from(bore.guard_mm);
        let points = [
            [c[0] + g, c[1]],
            [c[0], c[1] + g],
            [c[0] - g, c[1]],
            [c[0], c[1] - g],
        ];
        let suffixes = ["e", "n", "w", "s"];
        let mut vertices = Vec::with_capacity(4);
        for (point, suffix) in points.into_iter().zip(suffixes) {
            if point.iter().any(|x| !(-64..=64).contains(x)) {
                return Err("circular polygon guard vertex outside native domain");
            }
            let id = format!("g_{}_{}", bore.id, suffix);
            if !all_ids.insert(id.clone()) {
                return Err("circular polygon guard vertex id collision");
            }
            vertices.push(Vertex {
                id,
                xy_mm: [point[0] as i32, point[1] as i32],
            });
        }
        guard_loops.push(LoopProfile {
            id: loop_id,
            vertices,
        });
    }
    let mut blocks = build_perforated_prism(&PerforatedPrism {
        outer: recipe.outer.clone(),
        holes: guard_loops,
        height_mm: recipe.height_mm,
        top_shift_mm: recipe.top_shift_mm,
    })?;
    for bore in bores {
        for quadrant in 0..4 {
            for half in 0..2 {
                blocks.push(quarter_net(
                    bore,
                    recipe.height_mm,
                    recipe.top_shift_mm,
                    quadrant,
                    half,
                )?);
            }
        }
    }
    blocks.sort_by(|a, b| a.id.cmp(&b.id));
    Ok(blocks)
}
