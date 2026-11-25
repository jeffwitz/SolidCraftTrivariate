//! Exact constrained triangulation of polygonal profiles with strict holes.
use crate::bracket::BracketBlock;
use crate::polygon_prism::{append_triangle, PolygonPrism, Vertex};
use std::collections::{BTreeSet, HashSet};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LoopProfile {
    pub id: String,
    pub vertices: Vec<Vertex>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PerforatedPrism {
    pub outer: LoopProfile,
    pub holes: Vec<LoopProfile>,
    pub height_mm: i32,
    pub top_shift_mm: [i32; 2],
}

fn orient(a: [i32; 2], b: [i32; 2], c: [i32; 2]) -> i64 {
    i64::from(b[0] - a[0]) * i64::from(c[1] - a[1])
        - i64::from(b[1] - a[1]) * i64::from(c[0] - a[0])
}
fn on_segment(a: [i32; 2], b: [i32; 2], p: [i32; 2]) -> bool {
    orient(a, b, p) == 0 && (0..2).all(|k| a[k].min(b[k]) <= p[k] && p[k] <= a[k].max(b[k]))
}
fn intersects(a: [i32; 2], b: [i32; 2], c: [i32; 2], d: [i32; 2]) -> bool {
    let (x, y, z, w) = (
        orient(a, b, c),
        orient(a, b, d),
        orient(c, d, a),
        orient(c, d, b),
    );
    (x.signum() * y.signum() < 0 && z.signum() * w.signum() < 0)
        || on_segment(a, b, c)
        || on_segment(a, b, d)
        || on_segment(c, d, a)
        || on_segment(c, d, b)
}
fn point_in_loop(p2: [i64; 2], scale: i64, loop_: &[Vertex]) -> bool {
    let mut inside = false;
    for i in 0..loop_.len() {
        let a = loop_[i].xy_mm.map(i64::from);
        let b = loop_[(i + 1) % loop_.len()].xy_mm.map(i64::from);
        let (px, py) = (p2[0], p2[1]);
        let cross = (b[0] - a[0]) * (py - scale * a[1]) - (b[1] - a[1]) * (px - scale * a[0]);
        if cross == 0
            && (scale * a[0]).min(scale * b[0]) <= px
            && px <= (scale * a[0]).max(scale * b[0])
            && (scale * a[1]).min(scale * b[1]) <= py
            && py <= (scale * a[1]).max(scale * b[1])
        {
            return false;
        }
        if (scale * a[1] > py) != (scale * b[1] > py) {
            let lhs = (px - scale * a[0]) * (b[1] - a[1]);
            let rhs = (py - scale * a[1]) * (b[0] - a[0]);
            if (b[1] > a[1]) == (lhs < rhs) {
                inside = !inside;
            }
        }
    }
    inside
}
fn material(p2: [i64; 2], scale: i64, loops: &[LoopProfile]) -> bool {
    point_in_loop(p2, scale, &loops[0].vertices)
        && loops[1..]
            .iter()
            .all(|l| !point_in_loop(p2, scale, &l.vertices))
}
fn validate_loop(
    lp: &LoopProfile,
    ids: &mut BTreeSet<String>,
    vertices: &mut Vec<Vertex>,
) -> Result<(), &'static str> {
    if lp.id.is_empty()
        || lp.id.len() > 24
        || !lp
            .id
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || c == b'_')
        || lp.vertices.len() < 3
    {
        return Err("perforated_loop");
    }
    let n = lp.vertices.len();
    for (i, v) in lp.vertices.iter().enumerate() {
        if v.id.is_empty()
            || v.id.len() > 24
            || !v.id.bytes().all(|c| c.is_ascii_alphanumeric() || c == b'_')
            || !ids.insert(v.id.clone())
        {
            return Err("perforated_vertex_id");
        }
        if v.xy_mm.iter().any(|x| !(-64..=64).contains(x)) {
            return Err("perforated_coordinate_range");
        }
        if orient(
            lp.vertices[i].xy_mm,
            lp.vertices[(i + 1) % n].xy_mm,
            lp.vertices[(i + 2) % n].xy_mm,
        ) == 0
        {
            return Err("perforated_collinear_adjacent_edges");
        }
        for j in i + 1..n {
            if j == i + 1 || (i == 0 && j == n - 1) {
                continue;
            }
            if intersects(
                v.xy_mm,
                lp.vertices[(i + 1) % n].xy_mm,
                lp.vertices[j].xy_mm,
                lp.vertices[(j + 1) % n].xy_mm,
            ) {
                return Err("perforated_loop_contact");
            }
        }
        vertices.push(v.clone());
    }
    Ok(())
}

/// Build a deterministic constrained triangulation, split each triangle into three quads,
/// then extrude with the same exact homogeneous convention as `polygon_prism`.
pub fn build_perforated_prism(recipe: &PerforatedPrism) -> Result<Vec<BracketBlock>, &'static str> {
    if recipe.holes.len() > 4
        || !(1..=64).contains(&recipe.height_mm)
        || recipe.top_shift_mm.iter().any(|x| !(-32..=32).contains(x))
    {
        return Err("perforated_range");
    }
    let mut loops = vec![recipe.outer.clone()];
    loops.extend(recipe.holes.clone());
    let mut ids = BTreeSet::new();
    let mut all = Vec::new();
    let mut loop_ids = BTreeSet::new();
    for lp in &loops {
        if !loop_ids.insert(lp.id.clone()) {
            return Err("perforated_loop_id");
        }
        validate_loop(lp, &mut ids, &mut all)?;
    }
    if !(3..=32).contains(&all.len()) {
        return Err("perforated_vertex_count");
    }
    // Normalize orientations: outer CCW, holes CW. Validate strict containment and no contact.
    for (li, lp) in loops.iter_mut().enumerate() {
        let area: i64 = (0..lp.vertices.len())
            .map(|i| {
                let a = lp.vertices[i].xy_mm;
                let b = lp.vertices[(i + 1) % lp.vertices.len()].xy_mm;
                i64::from(a[0]) * i64::from(b[1]) - i64::from(a[1]) * i64::from(b[0])
            })
            .sum();
        if area == 0 {
            return Err("perforated_zero_area");
        }
        if (area > 0) != (li == 0) {
            lp.vertices.reverse();
        }
    }
    all = loops
        .iter()
        .flat_map(|lp| lp.vertices.iter().cloned())
        .collect();
    let mut edges = Vec::<(usize, usize)>::new();
    let mut edge_set = HashSet::<(usize, usize)>::new();
    let mut offset = 0;
    for lp in &loops {
        for i in 0..lp.vertices.len() {
            let e = (offset + i, offset + (i + 1) % lp.vertices.len());
            edges.push(e);
            edge_set.insert((e.0.min(e.1), e.0.max(e.1)));
        }
        offset += lp.vertices.len();
    }
    for h in 1..loops.len() {
        let p = loops[h].vertices[0].xy_mm.map(i64::from).map(|x| 2 * x);
        if !point_in_loop(p, 2, &loops[0].vertices) {
            return Err("perforated_hole_outside");
        }
        for (other, loop_profile) in loops.iter().enumerate().skip(1) {
            if other != h && point_in_loop(p, 2, &loop_profile.vertices) {
                return Err("perforated_nested_holes");
            }
        }
    }
    for i in 0..edges.len() {
        for j in i + 1..edges.len() {
            let (a, b) = edges[i];
            let (c, d) = edges[j];
            let adjacent = a == c || a == d || b == c || b == d;
            if !adjacent && intersects(all[a].xy_mm, all[b].xy_mm, all[c].xy_mm, all[d].xy_mm) {
                return Err("perforated_boundary_contact");
            }
            if adjacent {
                let shared = if a == c || a == d { a } else { b };
                let aa = if a == shared { b } else { a };
                let bb = if c == shared || d == shared {
                    if c == shared {
                        d
                    } else {
                        c
                    }
                } else {
                    c
                };
                if orient(all[aa].xy_mm, all[shared].xy_mm, all[bb].xy_mm) == 0 {
                    return Err("perforated_collinear_adjacent_edges");
                }
            }
        }
    }
    // Candidate diagonals are globally ordered; insert a segment iff it is visible in material and uncrossed.
    let n = all.len();
    let mut candidates = Vec::new();
    for a in 0..n {
        for b in a + 1..n {
            if edge_set.contains(&(a, b)) {
                continue;
            }
            let p = all[a].xy_mm;
            let q = all[b].xy_mm;
            let mid = [i64::from(p[0] + q[0]), i64::from(p[1] + q[1])];
            if !material(mid, 2, &loops) {
                continue;
            }
            let mut visible = true;
            for &(c, d) in &edges {
                if c == a || d == a || c == b || d == b {
                    continue;
                }
                if intersects(p, q, all[c].xy_mm, all[d].xy_mm) {
                    visible = false;
                    break;
                }
            }
            if visible {
                let (id_a, id_b) = if all[a].id < all[b].id {
                    (all[a].id.clone(), all[b].id.clone())
                } else {
                    (all[b].id.clone(), all[a].id.clone())
                };
                candidates.push((
                    (i64::from(p[0] - q[0])).pow(2) + (i64::from(p[1] - q[1])).pow(2),
                    id_a,
                    id_b,
                    a,
                    b,
                ));
            }
        }
    }
    candidates.sort();
    for (_, _, _, a, b) in candidates {
        let p = all[a].xy_mm;
        let q = all[b].xy_mm;
        let mut ok = true;
        for &(c, d) in &edges {
            if c == a || d == a || c == b || d == b {
                continue;
            }
            if intersects(p, q, all[c].xy_mm, all[d].xy_mm) {
                ok = false;
                break;
            }
        }
        if ok {
            edges.push((a, b));
            edge_set.insert((a, b));
        }
    }
    let mut triangles = Vec::<[usize; 3]>::new();
    for a in 0..n {
        for b in a + 1..n {
            for c in b + 1..n {
                if !edge_set.contains(&(a, b))
                    || !edge_set.contains(&(b, c))
                    || !edge_set.contains(&(a, c))
                {
                    continue;
                }
                let (x, mut y, mut z) = (a, b, c);
                let area = orient(all[x].xy_mm, all[y].xy_mm, all[z].xy_mm);
                if area == 0 {
                    continue;
                }
                if area < 0 {
                    std::mem::swap(&mut y, &mut z);
                }
                let centroid = [
                    i64::from(all[x].xy_mm[0] + all[y].xy_mm[0] + all[z].xy_mm[0]),
                    i64::from(all[x].xy_mm[1] + all[y].xy_mm[1] + all[z].xy_mm[1]),
                ];
                if !material(centroid, 3, &loops) {
                    continue;
                }
                let has_inside = (0..n).any(|p| {
                    p != x
                        && p != y
                        && p != z
                        && orient(all[x].xy_mm, all[y].xy_mm, all[p].xy_mm) >= 0
                        && orient(all[y].xy_mm, all[z].xy_mm, all[p].xy_mm) >= 0
                        && orient(all[z].xy_mm, all[x].xy_mm, all[p].xy_mm) >= 0
                });
                if !has_inside {
                    triangles.push([x, y, z]);
                }
            }
        }
    }
    let expected = n + 2 * recipe.holes.len() - 2;
    if triangles.len() != expected {
        return Err("perforated_triangulation_incomplete");
    }
    let prism = PolygonPrism {
        vertices: all.clone(),
        height_mm: recipe.height_mm,
        top_shift_mm: recipe.top_shift_mm,
    };
    let mut blocks = Vec::with_capacity(3 * triangles.len());
    for t in triangles {
        append_triangle([&all[t[0]], &all[t[1]], &all[t[2]]], &prism, &mut blocks);
    }
    blocks.sort_by(|a, b| a.id.cmp(&b.id));
    Ok(blocks)
}
