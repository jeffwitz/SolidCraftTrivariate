//! Exact straight-sided quadrangulation of a simple named polygon and a sheared extrusion.
use crate::bracket::BracketBlock;
use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Vertex {
    pub id: String,
    pub xy_mm: [i32; 2],
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolygonPrism {
    pub vertices: Vec<Vertex>,
    pub height_mm: i32,
    pub top_shift_mm: [i32; 2],
}

fn orient(a: [i32; 2], b: [i32; 2], c: [i32; 2]) -> i64 {
    let cross = |u: usize, v: usize| i64::from(b[u] - a[u]) * i64::from(c[v] - a[v]);
    cross(0, 1) - cross(1, 0)
}
fn on_segment(a: [i32; 2], b: [i32; 2], p: [i32; 2]) -> bool {
    orient(a, b, p) == 0 && (0..2).all(|d| (a[d].min(b[d])..=a[d].max(b[d])).contains(&p[d]))
}
fn touch(a: [i32; 2], b: [i32; 2], c: [i32; 2], d: [i32; 2]) -> bool {
    let (ac, ad, ca, cb) = (
        orient(a, b, c),
        orient(a, b, d),
        orient(c, d, a),
        orient(c, d, b),
    );
    (ac.signum() * ad.signum() < 0 && ca.signum() * cb.signum() < 0)
        || on_segment(a, b, c)
        || on_segment(a, b, d)
        || on_segment(c, d, a)
        || on_segment(c, d, b)
}

fn simple_contour(recipe: &PolygonPrism) -> Result<Vec<&Vertex>, &'static str> {
    let mut polygon = recipe.vertices.iter().collect::<Vec<_>>();
    let n = polygon.len();
    if !(3..=24).contains(&n)
        || !(1..=64).contains(&recipe.height_mm)
        || recipe.top_shift_mm.iter().any(|v| !(-32..=32).contains(v))
    {
        return Err("polygon_prism_range");
    }
    let mut ids = BTreeSet::new();
    let mut points = BTreeSet::new();
    for v in &polygon {
        if v.id.is_empty()
            || v.id.len() > 24
            || !v.id.bytes().all(|c| c.is_ascii_alphanumeric() || c == b'_')
            || !ids.insert(&v.id)
        {
            return Err("polygon_vertex_id");
        }
        if v.xy_mm.iter().any(|x| !(-64..=64).contains(x)) || !points.insert(v.xy_mm) {
            return Err("polygon_vertex_coordinates");
        }
    }
    let mut area = 0;
    for i in 0..n {
        let a = polygon[i].xy_mm;
        let b = polygon[(i + 1) % n].xy_mm;
        let c = polygon[(i + 2) % n].xy_mm;
        if orient(a, b, c) == 0 {
            return Err("polygon_collinear_adjacent_edges");
        }
        area += orient([0, 0], a, b);
        for j in i + 1..n {
            if j == i + 1 || (i == 0 && j == n - 1) {
                continue;
            }
            if touch(a, b, polygon[j].xy_mm, polygon[(j + 1) % n].xy_mm) {
                return Err("polygon_nonadjacent_edge_contact");
            }
        }
    }
    if area == 0 {
        return Err("polygon_zero_area");
    }
    if area < 0 {
        polygon.reverse();
    }
    let first = polygon
        .iter()
        .enumerate()
        .min_by_key(|(_, v)| &v.id)
        .unwrap()
        .0;
    polygon.rotate_left(first);
    Ok(polygon)
}

fn ear(polygon: &[&Vertex], i: usize) -> bool {
    let n = polygon.len();
    let (prev, next) = ((i + n - 1) % n, (i + 1) % n);
    let (a, b, c) = (polygon[prev].xy_mm, polygon[i].xy_mm, polygon[next].xy_mm);
    orient(a, b, c) > 0
        && polygon.iter().enumerate().all(|(j, p)| {
            [prev, i, next].contains(&j)
                || orient(a, b, p.xy_mm) < 0
                || orient(b, c, p.xy_mm) < 0
                || orient(c, a, p.xy_mm) < 0
        })
}

pub(crate) fn append_triangle(
    mut t: [&Vertex; 3],
    recipe: &PolygonPrism,
    blocks: &mut Vec<BracketBlock>,
) {
    let first = (0..3).min_by_key(|&i| &t[i].id).unwrap();
    t.rotate_left(first);
    let centroid: [i64; 2] =
        std::array::from_fn(|d| 2 * t.iter().map(|v| i64::from(v.xy_mm[d])).sum::<i64>());
    for i in 0..3 {
        let a = t[i];
        let next = t[(i + 1) % 3];
        let prev = t[(i + 2) % 3];
        let midpoint = |b: &Vertex| std::array::from_fn(|d| 3 * i64::from(a.xy_mm[d] + b.xy_mm[d]));
        let mut controls = Vec::with_capacity(8);
        for [x, y] in [
            a.xy_mm.map(|v| 6 * i64::from(v)),
            midpoint(prev),
            midpoint(next),
            centroid,
        ] {
            controls.push([x, y, 0, 6]);
            controls.push([
                x + 6 * i64::from(recipe.top_shift_mm[0]),
                y + 6 * i64::from(recipe.top_shift_mm[1]),
                6 * i64::from(recipe.height_mm),
                6,
            ]);
        }
        blocks.push(BracketBlock {
            id: format!("tri/{}/{}/{}/{}", t[0].id, t[1].id, t[2].id, a.id),
            degrees: [1, 1, 1],
            controls,
        });
    }
}

/// Triangulate by exact ears in a deterministic contour order, then split each triangle into three quads.
pub fn build_polygon_prism(recipe: &PolygonPrism) -> Result<Vec<BracketBlock>, &'static str> {
    let mut contour = simple_contour(recipe)?;
    let mut blocks = Vec::with_capacity(3 * (contour.len() - 2));
    while contour.len() > 3 {
        let n = contour.len();
        let i = (0..n)
            .filter(|&i| ear(&contour, i))
            .min_by_key(|&i| &contour[i].id)
            .ok_or("polygon_no_valid_ear")?;
        append_triangle(
            [contour[(i + n - 1) % n], contour[i], contour[(i + 1) % n]],
            recipe,
            &mut blocks,
        );
        contour.remove(i);
    }
    append_triangle([contour[0], contour[1], contour[2]], recipe, &mut blocks);
    blocks.sort_by(|a, b| a.id.cmp(&b.id));
    Ok(blocks)
}
