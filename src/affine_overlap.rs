//! Exact pairwise relation checks for integer affine parallelepipeds.
//!
//! Vertices use the tensor bit order `v0 + u*e0 + v*e1 + w*e2` with index
//! `u + 2*v + 4*w`. Coordinates are integer micrometres bounded to one metre.
//! This module handles only pair relations and declared whole-feature contact;
//! it does not prove global coverage of a multi-block complex.

use std::collections::BTreeSet;

const COORD_LIMIT_UM: i64 = 1_000_000;

/// One nondegenerate, directly oriented affine parallelepiped.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AffineParallelepiped {
    pub id: String,
    pub vertices_um: [[i64; 3]; 8],
    /// Stable global topological vertex identities; shared contact vertices
    /// carry the same identity in both blocks.
    pub vertex_ids: [u64; 8],
}

/// Exact pairwise relation. Boundary dimension is zero, one, or two.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairRelation {
    StrictlySeparated,
    BoundaryContact { dimension: usize },
    InteriorOverlap,
    FirstContainsSecond,
    SecondContainsFirst,
    Coincident,
    Indeterminate,
}

/// The claimed shared topological feature, expressed by global vertex IDs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContactFeature {
    Vertex(u64),
    Edge([u64; 2]),
    Face([u64; 4]),
}

/// Result of checking the physical intersection against a declared identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactDecision {
    Accepted,
    Rejected,
    Indeterminate,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContactAudit {
    pub decision: ContactDecision,
    pub relation: PairRelation,
    pub reason: &'static str,
}

type V = [i128; 3];

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct P2 {
    x: i128,
    y: i128,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct RP2 {
    x: i128,
    y: i128,
    d: i128,
}

/// Classify the two closed solids using exact separating-axis predicates.
pub fn classify_pair(
    first: &AffineParallelepiped,
    second: &AffineParallelepiped,
) -> Result<PairRelation, &'static str> {
    let a = validate(first)?;
    let b = validate(second)?;
    let axes = sat_axes(&a, &b)?;
    let mut touching = Vec::new();
    for axis in axes {
        if axis == [0, 0, 0] {
            continue;
        }
        let (amin, amax) = projection(&a, axis)?;
        let (bmin, bmax) = projection(&b, axis)?;
        if amax < bmin || bmax < amin {
            return Ok(PairRelation::StrictlySeparated);
        }
        if amax == bmin {
            touching.push((axis, true));
        } else if bmax == amin {
            touching.push((axis, false));
        }
    }
    if !touching.is_empty() {
        let mut dimension = None;
        for (axis, a_max) in touching {
            let Ok((amin, amax)) = projection(&a, axis) else {
                return Ok(PairRelation::Indeterminate);
            };
            let Ok((bmin, bmax)) = projection(&b, axis) else {
                return Ok(PairRelation::Indeterminate);
            };
            let (sa, sb) = if a_max { (amax, bmin) } else { (amin, bmax) };
            let a_support = support_points(&a, axis, sa);
            let b_support = support_points(&b, axis, sb);
            let Ok(pa) = project_support(&a_support, axis) else {
                return Ok(PairRelation::Indeterminate);
            };
            let Ok(pb) = project_support(&b_support, axis) else {
                return Ok(PairRelation::Indeterminate);
            };
            let Ok(Some(found)) = intersection_dimension(&pa, &pb) else {
                return Ok(PairRelation::Indeterminate);
            };
            if dimension.is_some_and(|previous| previous != found) {
                return Ok(PairRelation::Indeterminate);
            }
            dimension = Some(found);
        }
        return Ok(dimension.map_or(PairRelation::Indeterminate, |dimension| {
            PairRelation::BoundaryContact { dimension }
        }));
    }

    let mut a_in_b = true;
    let mut b_in_a = true;
    for p in &a {
        match point_in(&b, *p) {
            Ok(inside) => a_in_b &= inside,
            Err(_) => return Ok(PairRelation::Indeterminate),
        }
    }
    for p in &b {
        match point_in(&a, *p) {
            Ok(inside) => b_in_a &= inside,
            Err(_) => return Ok(PairRelation::Indeterminate),
        }
    }
    Ok(match (a_in_b, b_in_a) {
        (true, true) => PairRelation::Coincident,
        (true, false) => PairRelation::SecondContainsFirst,
        (false, true) => PairRelation::FirstContainsSecond,
        (false, false) => PairRelation::InteriorOverlap,
    })
}

/// Verify that the pair's exact geometric intersection is precisely the
/// declared shared vertex, edge, or whole face.
pub fn audit_declared_contact(
    first: &AffineParallelepiped,
    second: &AffineParallelepiped,
    feature: &ContactFeature,
) -> Result<ContactAudit, &'static str> {
    let relation = classify_pair(first, second)?;
    let accepted = match feature {
        ContactFeature::Vertex(id) => {
            let (ia, ib) = shared_indices(first, second, &[*id]);
            if ia.is_none() || ib.is_none() {
                return Ok(ContactAudit {
                    decision: ContactDecision::Rejected,
                    relation,
                    reason: "declared vertex identity is not shared at equal coordinates",
                });
            }
            matches!(relation, PairRelation::BoundaryContact { dimension: 0 })
        }
        ContactFeature::Edge(ids) => {
            let (ia, ib) = shared_indices(first, second, ids);
            if ia.is_none()
                || ib.is_none()
                || !is_edge(first, ia.as_ref().unwrap())
                || !is_edge(second, ib.as_ref().unwrap())
                || !distinct_points(first, ia.as_ref().unwrap())
            {
                return Ok(ContactAudit {
                    decision: ContactDecision::Rejected,
                    relation,
                    reason: "declared edge identities do not form the same nondegenerate edge",
                });
            }
            matches!(relation, PairRelation::BoundaryContact { dimension: 1 })
        }
        ContactFeature::Face(ids) => {
            let (ia, ib) = shared_indices(first, second, ids);
            if ia.is_none()
                || ib.is_none()
                || !is_face(first, ia.as_ref().unwrap())
                || !is_face(second, ib.as_ref().unwrap())
                || !noncollinear_face(first, ia.as_ref().unwrap())
            {
                return Ok(ContactAudit {
                    decision: ContactDecision::Rejected,
                    relation,
                    reason:
                        "declared face identities do not form the same nondegenerate whole face",
                });
            }
            matches!(relation, PairRelation::BoundaryContact { dimension: 2 })
        }
    };
    if accepted {
        Ok(ContactAudit {
            decision: ContactDecision::Accepted,
            relation,
            reason: "the exact pair intersection equals the declared whole feature",
        })
    } else if relation == PairRelation::Indeterminate {
        Ok(ContactAudit {
            decision: ContactDecision::Indeterminate,
            relation,
            reason: "exact contact dimension could not be established",
        })
    } else {
        Ok(ContactAudit {
            decision: ContactDecision::Rejected,
            relation,
            reason: "geometric intersection is not exactly the declared feature",
        })
    }
}

fn validate(block: &AffineParallelepiped) -> Result<[V; 8], &'static str> {
    if block.id.is_empty()
        || block
            .vertices_um
            .iter()
            .flatten()
            .any(|x| x.unsigned_abs() > COORD_LIMIT_UM as u64)
    {
        return Err("empty block ID or coordinate outside +/-1 m domain");
    }
    let ids: BTreeSet<_> = block.vertex_ids.iter().copied().collect();
    if ids.len() != 8 {
        return Err("vertex identities must be unique within a block");
    }
    let v0 = block.vertices_um[0].map(i128::from);
    let edges = [
        sub(block.vertices_um[1].map(i128::from), v0),
        sub(block.vertices_um[2].map(i128::from), v0),
        sub(block.vertices_um[4].map(i128::from), v0),
    ];
    let det = dot(edges[0], cross(edges[1], edges[2])?);
    if det <= 0 {
        return Err("parallelepiped must be nondegenerate and directly oriented");
    }
    for i in 0..8 {
        let expected = add(
            v0,
            add_scaled(
                edges[0],
                (i & 1) as i128,
                add_scaled(
                    edges[1],
                    ((i >> 1) & 1) as i128,
                    scale(edges[2], ((i >> 2) & 1) as i128),
                ),
            ),
        );
        if block.vertices_um[i].map(i128::from) != expected {
            return Err("vertices do not follow the affine tensor bit ordering");
        }
    }
    Ok(block.vertices_um.map(|p| p.map(i128::from)))
}

fn sat_axes(a: &[V; 8], b: &[V; 8]) -> Result<Vec<V>, &'static str> {
    let ea = edges(a);
    let eb = edges(b);
    let mut axes = Vec::with_capacity(15);
    for e in [&ea, &eb] {
        axes.push(cross(e[0], e[1])?);
        axes.push(cross(e[1], e[2])?);
        axes.push(cross(e[2], e[0])?);
    }
    for x in ea {
        for y in eb {
            axes.push(cross(x, y)?);
        }
    }
    Ok(axes)
}

fn edges(v: &[V; 8]) -> [V; 3] {
    [sub(v[1], v[0]), sub(v[2], v[0]), sub(v[4], v[0])]
}

fn projection(v: &[V; 8], axis: V) -> Result<(i128, i128), &'static str> {
    let mut iter = v.iter().map(|p| checked_dot(*p, axis));
    let first = iter.next().ok_or("empty vertices")??;
    let mut min = first;
    let mut max = first;
    for p in iter {
        let p = p?;
        min = min.min(p);
        max = max.max(p);
    }
    Ok((min, max))
}

fn point_in(block: &[V; 8], p: V) -> Result<bool, &'static str> {
    let e = edges(block);
    let det = dot(e[0], cross(e[1], e[2])?);
    let q = sub(p, block[0]);
    let n0 = dot(q, cross(e[1], e[2])?);
    let n1 = dot(e[0], cross(q, e[2])?);
    let n2 = dot(e[0], cross(e[1], q)?);
    Ok([n0, n1, n2].iter().all(|n| *n >= 0 && *n <= det))
}

fn support_points(vertices: &[V; 8], axis: V, projection_value: i128) -> Vec<V> {
    vertices
        .iter()
        .copied()
        .filter(|p| dot(*p, axis) == projection_value)
        .collect()
}

fn project_support(points: &[V], normal: V) -> Result<Vec<P2>, &'static str> {
    let drop = (0..3)
        .max_by_key(|i| normal[*i].unsigned_abs())
        .ok_or("missing normal")?;
    if normal[drop] == 0 {
        return Err("zero normal");
    }
    Ok(points
        .iter()
        .map(|p| match drop {
            0 => P2 { x: p[1], y: p[2] },
            1 => P2 { x: p[0], y: p[2] },
            _ => P2 { x: p[0], y: p[1] },
        })
        .collect())
}

fn intersection_dimension(a: &[P2], b: &[P2]) -> Result<Option<usize>, &'static str> {
    let a = convex_hull(a.to_vec());
    let b = convex_hull(b.to_vec());
    if a.is_empty() || b.is_empty() {
        return Ok(None);
    }
    let mut points = Vec::new();
    for p in &a {
        if point_in_hull(*p, &b) {
            points.push(RP2::integer(*p));
        }
    }
    for p in &b {
        if point_in_hull(*p, &a) {
            points.push(RP2::integer(*p));
        }
    }
    for (a0, a1) in hull_segments(&a) {
        for (b0, b1) in hull_segments(&b) {
            segment_crossings(a0, a1, b0, b1, &mut points)?;
        }
    }
    points.sort();
    points.dedup();
    if points.is_empty() {
        return Ok(None);
    }
    let mut uncertain = false;
    for i in 0..points.len() {
        for j in i + 1..points.len() {
            for k in j + 1..points.len() {
                match orientation_rational(points[i], points[j], points[k]) {
                    Ok(value) if value != 0 => return Ok(Some(2)),
                    Ok(_) => {}
                    Err(_) => uncertain = true,
                }
            }
        }
    }
    if uncertain {
        return Err("intersection dimension exceeds exact i128 range");
    }
    Ok(Some(if points.len() > 1 { 1 } else { 0 }))
}

impl RP2 {
    fn integer(p: P2) -> Self {
        Self {
            x: p.x,
            y: p.y,
            d: 1,
        }
    }
    fn new(mut x: i128, mut y: i128, mut d: i128) -> Result<Self, &'static str> {
        if d == 0 {
            return Err("zero rational denominator");
        }
        if d < 0 {
            x = x.checked_neg().ok_or("rational overflow")?;
            y = y.checked_neg().ok_or("rational overflow")?;
            d = d.checked_neg().ok_or("rational overflow")?;
        }
        let g = gcd(gcd(x.unsigned_abs(), y.unsigned_abs()), d as u128) as i128;
        Ok(Self {
            x: x / g,
            y: y / g,
            d: d / g,
        })
    }
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

fn convex_hull(mut p: Vec<P2>) -> Vec<P2> {
    p.sort();
    p.dedup();
    if p.len() <= 2 {
        return p;
    }
    let mut lo = Vec::new();
    for q in &p {
        while lo.len() >= 2 && orient2(lo[lo.len() - 2], lo[lo.len() - 1], *q) <= 0 {
            lo.pop();
        }
        lo.push(*q);
    }
    let mut hi = Vec::new();
    for q in p.iter().rev() {
        while hi.len() >= 2 && orient2(hi[hi.len() - 2], hi[hi.len() - 1], *q) <= 0 {
            hi.pop();
        }
        hi.push(*q);
    }
    lo.pop();
    hi.pop();
    lo.extend(hi);
    lo
}

fn orient2(a: P2, b: P2, c: P2) -> i128 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}
fn hull_segments(p: &[P2]) -> Vec<(P2, P2)> {
    match p.len() {
        0 | 1 => Vec::new(),
        2 => vec![(p[0], p[1])],
        _ => (0..p.len()).map(|i| (p[i], p[(i + 1) % p.len()])).collect(),
    }
}
fn point_in_hull(p: P2, h: &[P2]) -> bool {
    match h.len() {
        0 => false,
        1 => p == h[0],
        2 => orient2(h[0], h[1], p) == 0 && on_segment(p, h[0], h[1]),
        _ => (0..h.len()).all(|i| orient2(h[i], h[(i + 1) % h.len()], p) >= 0),
    }
}
fn on_segment(p: P2, a: P2, b: P2) -> bool {
    p.x >= a.x.min(b.x) && p.x <= a.x.max(b.x) && p.y >= a.y.min(b.y) && p.y <= a.y.max(b.y)
}
fn cross2(a: P2, b: P2) -> i128 {
    a.x * b.y - a.y * b.x
}

fn segment_crossings(a: P2, b: P2, c: P2, d: P2, out: &mut Vec<RP2>) -> Result<(), &'static str> {
    let r = P2 {
        x: b.x - a.x,
        y: b.y - a.y,
    };
    let s = P2 {
        x: d.x - c.x,
        y: d.y - c.y,
    };
    let q = P2 {
        x: c.x - a.x,
        y: c.y - a.y,
    };
    let den = cross2(r, s);
    let tn = cross2(q, s);
    let un = cross2(q, r);
    if den == 0 {
        if cross2(q, r) != 0 {
            return Ok(());
        }
        for p in [a, b, c, d] {
            if on_segment(p, a, b) && on_segment(p, c, d) {
                out.push(RP2::integer(p));
            }
        }
        return Ok(());
    }
    if !in_unit(tn, den) || !in_unit(un, den) {
        return Ok(());
    }
    let x =
        a.x.checked_mul(den)
            .and_then(|v| r.x.checked_mul(tn).and_then(|w| v.checked_add(w)))
            .ok_or("intersection overflow")?;
    let y =
        a.y.checked_mul(den)
            .and_then(|v| r.y.checked_mul(tn).and_then(|w| v.checked_add(w)))
            .ok_or("intersection overflow")?;
    out.push(RP2::new(x, y, den)?);
    Ok(())
}
fn in_unit(n: i128, d: i128) -> bool {
    if d > 0 {
        n >= 0 && n <= d
    } else {
        n <= 0 && n >= d
    }
}

fn orientation_rational(a: RP2, b: RP2, c: RP2) -> Result<i128, &'static str> {
    let ax =
        b.x.checked_mul(a.d)
            .and_then(|v| a.x.checked_mul(b.d).and_then(|w| v.checked_sub(w)))
            .ok_or("orientation overflow")?;
    let ay =
        b.y.checked_mul(a.d)
            .and_then(|v| a.y.checked_mul(b.d).and_then(|w| v.checked_sub(w)))
            .ok_or("orientation overflow")?;
    let bx =
        c.x.checked_mul(a.d)
            .and_then(|v| a.x.checked_mul(c.d).and_then(|w| v.checked_sub(w)))
            .ok_or("orientation overflow")?;
    let by =
        c.y.checked_mul(a.d)
            .and_then(|v| a.y.checked_mul(c.d).and_then(|w| v.checked_sub(w)))
            .ok_or("orientation overflow")?;
    ax.checked_mul(by)
        .and_then(|u| ay.checked_mul(bx).and_then(|v| u.checked_sub(v)))
        .ok_or("orientation overflow")
}

fn shared_indices(
    a: &AffineParallelepiped,
    b: &AffineParallelepiped,
    ids: &[u64],
) -> (Option<Vec<usize>>, Option<Vec<usize>>) {
    let lookup = |block: &AffineParallelepiped| {
        ids.iter()
            .map(|id| block.vertex_ids.iter().position(|x| x == id))
            .collect::<Option<Vec<_>>>()
    };
    let ia = lookup(a);
    let ib = lookup(b);
    let equal = ia
        .as_ref()
        .zip(ib.as_ref())
        .map(|(x, y)| {
            x.iter()
                .zip(y)
                .all(|(i, j)| a.vertices_um[*i] == b.vertices_um[*j])
        })
        .unwrap_or(false);
    if equal {
        (ia, ib)
    } else {
        (None, None)
    }
}

fn is_edge(block: &AffineParallelepiped, indices: &[usize]) -> bool {
    indices.len() == 2
        && (indices[0] ^ indices[1]).count_ones() == 1
        && block.vertices_um[indices[0]] != block.vertices_um[indices[1]]
}

fn is_face(block: &AffineParallelepiped, indices: &[usize]) -> bool {
    if indices.len() != 4 {
        return false;
    }
    (0..3).any(|bit| {
        let fixed = (indices[0] >> bit) & 1;
        indices.iter().all(|i| ((i >> bit) & 1) == fixed)
            && indices.iter().copied().collect::<BTreeSet<_>>().len() == 4
            && indices.iter().all(|i| *i < 8)
            && indices
                .iter()
                .map(|i| block.vertices_um[*i])
                .collect::<BTreeSet<_>>()
                .len()
                == 4
            && indices.iter().all(|i| ((i >> bit) & 1) == fixed)
    })
}

fn distinct_points(block: &AffineParallelepiped, indices: &[usize]) -> bool {
    indices.len() == 2 && block.vertices_um[indices[0]] != block.vertices_um[indices[1]]
}

fn noncollinear_face(block: &AffineParallelepiped, indices: &[usize]) -> bool {
    let p = indices
        .iter()
        .map(|i| block.vertices_um[*i].map(i128::from))
        .collect::<Vec<_>>();
    p.len() == 4
        && cross(sub(p[1], p[0]), sub(p[2], p[0]))
            .map(|n| n != [0, 0, 0])
            .unwrap_or(false)
}

fn add(a: V, b: V) -> V {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn sub(a: V, b: V) -> V {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn scale(a: V, s: i128) -> V {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn add_scaled(a: V, s: i128, b: V) -> V {
    add(scale(a, s), b)
}
fn dot(a: V, b: V) -> i128 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn checked_dot(a: V, b: V) -> Result<i128, &'static str> {
    a[0].checked_mul(b[0])
        .and_then(|x| a[1].checked_mul(b[1]).and_then(|y| x.checked_add(y)))
        .and_then(|x| a[2].checked_mul(b[2]).and_then(|z| x.checked_add(z)))
        .ok_or("projection overflow")
}
fn cross(a: V, b: V) -> Result<V, &'static str> {
    let c = [
        a[1].checked_mul(b[2])
            .and_then(|x| a[2].checked_mul(b[1]).and_then(|y| x.checked_sub(y))),
        a[2].checked_mul(b[0])
            .and_then(|x| a[0].checked_mul(b[2]).and_then(|y| x.checked_sub(y))),
        a[0].checked_mul(b[1])
            .and_then(|x| a[1].checked_mul(b[0]).and_then(|y| x.checked_sub(y))),
    ];
    Ok([
        c[0].ok_or("cross product overflow")?,
        c[1].ok_or("cross product overflow")?,
        c[2].ok_or("cross product overflow")?,
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cube(
        id: &str,
        origin: [i64; 3],
        edges: [[i64; 3]; 3],
        ids: [u64; 8],
    ) -> AffineParallelepiped {
        let mut vertices_um = [[0; 3]; 8];
        for (i, vertex) in vertices_um.iter_mut().enumerate() {
            for (j, coordinate) in vertex.iter_mut().enumerate() {
                *coordinate = origin[j]
                    + ((i & 1) as i64) * edges[0][j]
                    + (((i >> 1) & 1) as i64) * edges[1][j]
                    + (((i >> 2) & 1) as i64) * edges[2][j];
            }
        }
        AffineParallelepiped {
            id: id.into(),
            vertices_um,
            vertex_ids: ids,
        }
    }

    fn ids(start: u64) -> [u64; 8] {
        std::array::from_fn(|i| start + i as u64)
    }

    #[test]
    fn classifies_separation_overlap_containment_and_coincidence() {
        let unit = [[1000, 0, 0], [0, 1000, 0], [0, 0, 1000]];
        let a = cube("a", [0, 0, 0], unit, ids(0));
        assert_eq!(
            classify_pair(&a, &cube("b", [2000, 0, 0], unit, ids(10))).unwrap(),
            PairRelation::StrictlySeparated
        );
        assert_eq!(
            classify_pair(&a, &cube("b", [500, 0, 0], unit, ids(10))).unwrap(),
            PairRelation::InteriorOverlap
        );
        assert_eq!(
            classify_pair(
                &a,
                &cube(
                    "b",
                    [200, 200, 200],
                    [[200, 0, 0], [0, 200, 0], [0, 0, 200]],
                    ids(10)
                )
            )
            .unwrap(),
            PairRelation::FirstContainsSecond
        );
        assert_eq!(
            classify_pair(&a, &cube("b", [0, 0, 0], unit, ids(10))).unwrap(),
            PairRelation::Coincident
        );
    }

    #[test]
    fn accepts_only_intersection_equal_to_declared_feature() {
        let unit = [[1000, 0, 0], [0, 1000, 0], [0, 0, 1000]];
        let mut ai = ids(100);
        let mut bi = ids(200);
        ai[7] = 7;
        bi[0] = 7;
        let a = cube("a", [0, 0, 0], unit, ai);
        let b = cube("b", [1000, 1000, 1000], unit, bi);
        assert_eq!(
            audit_declared_contact(&a, &b, &ContactFeature::Vertex(7))
                .unwrap()
                .decision,
            ContactDecision::Accepted
        );

        // The blocks meet along an edge. Declaring only one endpoint is not enough.
        let mut ai = ids(300);
        let mut bi = ids(400);
        ai[3] = 9;
        ai[7] = 10;
        bi[0] = 9;
        bi[4] = 10;
        let edge_a = cube("edge-a", [0, 0, 0], unit, ai);
        let edge_b = cube("edge-b", [1000, 1000, 0], unit, bi);
        assert_eq!(
            classify_pair(&edge_a, &edge_b).unwrap(),
            PairRelation::BoundaryContact { dimension: 1 }
        );
        assert_eq!(
            audit_declared_contact(&edge_a, &edge_b, &ContactFeature::Vertex(9))
                .unwrap()
                .decision,
            ContactDecision::Rejected
        );
        assert_eq!(
            audit_declared_contact(&edge_a, &edge_b, &ContactFeature::Edge([9, 10]))
                .unwrap()
                .decision,
            ContactDecision::Accepted
        );

        let mut ai = ids(500);
        let mut bi = ids(600);
        for (j, index) in [1, 3, 5, 7].into_iter().enumerate() {
            ai[index] = 30 + j as u64;
        }
        for (j, index) in [0, 2, 4, 6].into_iter().enumerate() {
            bi[index] = 30 + j as u64;
        }
        let face_a = cube("face-a", [0, 0, 0], unit, ai);
        let face_b = cube("face-b", [1000, 0, 0], unit, bi);
        assert_eq!(
            audit_declared_contact(&face_a, &face_b, &ContactFeature::Face([30, 31, 32, 33]))
                .unwrap()
                .decision,
            ContactDecision::Accepted
        );

        let partial_face =
            audit_declared_contact(&face_a, &face_b, &ContactFeature::Face([30, 31, 32, 999]))
                .unwrap();
        assert_eq!(partial_face.decision, ContactDecision::Rejected);
    }

    #[test]
    fn support_polygon_intersection_finds_oblique_vertex_and_edge_contacts() {
        let cube_edges = [[4000, 0, 0], [0, 4000, 0], [0, 0, 4000]];
        let face = cube("face", [0, 0, 0], cube_edges, ids(0));
        // The oblique block touches the relative interior of z=4 at one vertex.
        let oblique_vertex = cube(
            "oblique-vertex",
            [1000, 1000, 4000],
            [[1000, 0, 1000], [0, 1000, 1000], [-1000, -1000, 1000]],
            ids(20),
        );
        assert_eq!(
            classify_pair(&face, &oblique_vertex).unwrap(),
            PairRelation::BoundaryContact { dimension: 0 }
        );

        // One oblique edge lies in the face interior; the two solids meet on that segment.
        let oblique_edge = cube(
            "oblique-edge",
            [1000, 1000, 4000],
            [[2000, 1000, 0], [-1000, 0, 1000], [0, 0, 1000]],
            ids(40),
        );
        assert_eq!(
            classify_pair(&face, &oblique_edge).unwrap(),
            PairRelation::BoundaryContact { dimension: 1 }
        );
    }

    #[test]
    fn refuses_bad_orientation_and_identity_or_returns_indeterminate_for_unknown_contact() {
        let unit = [[1000, 0, 0], [0, 1000, 0], [0, 0, 1000]];
        let mut reversed = cube("bad", [0, 0, 0], [unit[1], unit[0], unit[2]], ids(0));
        reversed.vertices_um.swap(1, 2);
        assert!(classify_pair(&reversed, &cube("good", [0, 0, 0], unit, ids(10))).is_err());
        let mut extreme = cube("extreme", [0, 0, 0], unit, ids(20));
        extreme.vertices_um[0][0] = i64::MIN;
        assert!(classify_pair(&extreme, &cube("good", [0, 0, 0], unit, ids(30))).is_err());
        let a = cube("a", [0, 0, 0], unit, ids(100));
        let b = cube("b", [1000, 1000, 1000], unit, ids(200));
        let out = audit_declared_contact(&a, &b, &ContactFeature::Vertex(999)).unwrap();
        assert_eq!(out.decision, ContactDecision::Rejected);
        assert_ne!(out.relation, PairRelation::Indeterminate);
    }
}
