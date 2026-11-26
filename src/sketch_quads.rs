//! Exact deterministic quadrangulation of simple counter-clockwise rational faces.
//!
//! Each face is ear-clipped into positive rational triangles. Every triangle
//! is replaced by three convex quads using edge midpoints and its centroid.
//! Consequently every input boundary edge is split at its exact midpoint in
//! every adjacent face, independent of that face's triangulation.

use crate::{
    constraint_system::Ratio,
    sketch_arrangement::{arrange, Arrangement, Role, Segment},
    sketch_layers::{edges, SketchLayers},
};
use std::cmp::Ordering;

/// Provenance of a quadrilateral edge lying on an input face edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BoundaryEdgeSource {
    /// Index `i` refers to input edge `face[i] -> face[(i+1)%n]`.
    pub edge_index: usize,
    /// Zero is the first half from the edge start to its midpoint; one is the second half.
    pub half: u8,
}

/// Quadrilateral coordinates with per-edge source provenance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuadWithProvenance {
    pub vertices: [[Ratio; 2]; 4],
    /// `Some` only on the two half-edges subdividing an original face edge.
    pub edge_sources: [Option<BoundaryEdgeSource>; 4],
}

fn compare_positive(mut a: u128, mut b: u128, mut c: u128, mut d: u128) -> Ordering {
    let mut reversed = false;
    loop {
        let qa = a / b;
        let qc = c / d;
        if qa != qc {
            let order = qa.cmp(&qc);
            return if reversed { order.reverse() } else { order };
        }
        let ra = a % b;
        let rc = c % d;
        match (ra == 0, rc == 0) {
            (true, true) => return Ordering::Equal,
            (true, false) => {
                return if reversed {
                    Ordering::Greater
                } else {
                    Ordering::Less
                }
            }
            (false, true) => {
                return if reversed {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            }
            (false, false) => {
                (a, b, c, d) = (b, ra, d, rc);
                reversed = !reversed;
            }
        }
    }
}

/// Compare canonical ratios exactly without cross multiplication overflow.
fn cmp_ratio(a: Ratio, b: Ratio) -> Ordering {
    match (a.numerator().signum(), b.numerator().signum()) {
        (x, y) if x != y => x.cmp(&y),
        (0, 0) => Ordering::Equal,
        (1, 1) => compare_positive(
            a.numerator() as u128,
            a.denominator() as u128,
            b.numerator() as u128,
            b.denominator() as u128,
        ),
        (-1, -1) => compare_positive(
            a.numerator().unsigned_abs(),
            a.denominator() as u128,
            b.numerator().unsigned_abs(),
            b.denominator() as u128,
        )
        .reverse(),
        _ => unreachable!("sign comparison covers all canonical ratios"),
    }
}

fn cmp_point(a: [Ratio; 2], b: [Ratio; 2]) -> Ordering {
    cmp_ratio(a[0], b[0]).then_with(|| cmp_ratio(a[1], b[1]))
}

fn sub(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_sub(b)
}

pub(crate) fn cross(a: [Ratio; 2], b: [Ratio; 2], c: [Ratio; 2]) -> Result<Ratio, &'static str> {
    let abx = sub(b[0], a[0])?;
    let aby = sub(b[1], a[1])?;
    let acx = sub(c[0], a[0])?;
    let acy = sub(c[1], a[1])?;
    abx.checked_mul(acy)?.checked_sub(aby.checked_mul(acx)?)
}

fn sign(x: Ratio) -> i8 {
    x.numerator().signum() as i8
}

fn between(a: Ratio, x: Ratio, b: Ratio) -> bool {
    cmp_ratio(a, x) != Ordering::Greater && cmp_ratio(x, b) != Ordering::Greater
        || cmp_ratio(b, x) != Ordering::Greater && cmp_ratio(x, a) != Ordering::Greater
}

fn on_segment(a: [Ratio; 2], b: [Ratio; 2], p: [Ratio; 2]) -> bool {
    between(a[0], p[0], b[0]) && between(a[1], p[1], b[1])
}

fn segments_intersect(
    a: [Ratio; 2],
    b: [Ratio; 2],
    c: [Ratio; 2],
    d: [Ratio; 2],
) -> Result<bool, &'static str> {
    let (o1, o2, o3, o4) = (
        sign(cross(a, b, c)?),
        sign(cross(a, b, d)?),
        sign(cross(c, d, a)?),
        sign(cross(c, d, b)?),
    );
    if o1 * o2 < 0 && o3 * o4 < 0 {
        return Ok(true);
    }
    Ok((o1 == 0 && on_segment(a, b, c))
        || (o2 == 0 && on_segment(a, b, d))
        || (o3 == 0 && on_segment(c, d, a))
        || (o4 == 0 && on_segment(c, d, b)))
}

fn signed_area2(face: &[[Ratio; 2]]) -> Result<Ratio, &'static str> {
    let mut area = Ratio::zero();
    for i in 0..face.len() {
        let a = face[i];
        let b = face[(i + 1) % face.len()];
        area = area.checked_add(
            a[0].checked_mul(b[1])?
                .checked_sub(a[1].checked_mul(b[0])?)?,
        )?;
    }
    Ok(area)
}

fn validate_face(face: &[[Ratio; 2]]) -> Result<(), &'static str> {
    if face.len() < 3 {
        return Err("a face needs at least three vertices");
    }
    for i in 0..face.len() {
        if face[i] == face[(i + 1) % face.len()] {
            return Err("a face has a zero-length edge");
        }
        let previous = face[(i + face.len() - 1) % face.len()];
        let current = face[i];
        let next = face[(i + 1) % face.len()];
        if sign(cross(previous, current, next)?) == 0 && !on_segment(previous, next, current) {
            return Err("adjacent collinear edges overlap");
        }
    }
    let area = signed_area2(face)?;
    if sign(area) <= 0 {
        return Err("face vertices must have positive counter-clockwise area");
    }
    for i in 0..face.len() {
        let i_next = (i + 1) % face.len();
        for j in (i + 1)..face.len() {
            let j_next = (j + 1) % face.len();
            if i == j || i_next == j || j_next == i {
                continue;
            }
            if segments_intersect(face[i], face[i_next], face[j], face[j_next])? {
                return Err("face boundary is not simple");
            }
        }
    }
    Ok(())
}

fn point_in_ccw_triangle(
    point: [Ratio; 2],
    triangle: [[Ratio; 2]; 3],
) -> Result<bool, &'static str> {
    Ok(sign(cross(triangle[0], triangle[1], point)?) >= 0
        && sign(cross(triangle[1], triangle[2], point)?) >= 0
        && sign(cross(triangle[2], triangle[0], point)?) >= 0)
}

fn triangulate(face: &[[Ratio; 2]]) -> Result<Vec<[usize; 3]>, &'static str> {
    let mut remaining = (0..face.len()).collect::<Vec<_>>();
    let mut triangles = Vec::with_capacity(face.len() - 2);
    while remaining.len() > 3 {
        let mut ear = None;
        for position in 0..remaining.len() {
            let previous = remaining[(position + remaining.len() - 1) % remaining.len()];
            let current = remaining[position];
            let next = remaining[(position + 1) % remaining.len()];
            let triangle = [face[previous], face[current], face[next]];
            if sign(cross(triangle[0], triangle[1], triangle[2])?) <= 0 {
                continue;
            }
            let mut blocked = false;
            for &candidate in &remaining {
                if candidate != previous
                    && candidate != current
                    && candidate != next
                    && point_in_ccw_triangle(face[candidate], triangle)?
                {
                    blocked = true;
                    break;
                }
            }
            if blocked {
                continue;
            }
            // The ear diagonal may not touch or cross another current boundary edge.
            for edge in 0..remaining.len() {
                let a = remaining[edge];
                let b = remaining[(edge + 1) % remaining.len()];
                if a == previous || b == previous || a == next || b == next {
                    continue;
                }
                if segments_intersect(face[previous], face[next], face[a], face[b])? {
                    blocked = true;
                    break;
                }
            }
            if !blocked {
                ear = Some((position, [previous, current, next]));
                break;
            }
        }
        let (position, triangle) = ear.ok_or("exact ear triangulation found no admissible ear")?;
        triangles.push(triangle);
        remaining.remove(position);
    }
    let last = [remaining[0], remaining[1], remaining[2]];
    if sign(cross(face[last[0]], face[last[1]], face[last[2]])?) <= 0 {
        return Err("triangulation ended with a degenerate face");
    }
    triangles.push(last);
    Ok(triangles)
}

pub(crate) fn midpoint(a: [Ratio; 2], b: [Ratio; 2]) -> Result<[Ratio; 2], &'static str> {
    let half = Ratio::new(1, 2)?;
    Ok([
        a[0].checked_add(b[0])?.checked_mul(half)?,
        a[1].checked_add(b[1])?.checked_mul(half)?,
    ])
}

pub(crate) fn centroid(triangle: [[Ratio; 2]; 3]) -> Result<[Ratio; 2], &'static str> {
    let third = Ratio::new(1, 3)?;
    Ok([
        triangle[0][0]
            .checked_add(triangle[1][0])?
            .checked_add(triangle[2][0])?
            .checked_mul(third)?,
        triangle[0][1]
            .checked_add(triangle[1][1])?
            .checked_add(triangle[2][1])?
            .checked_mul(third)?,
    ])
}

fn triangle_edge_source(
    face_len: usize,
    start: usize,
    end: usize,
    half: u8,
) -> Option<BoundaryEdgeSource> {
    if (start + 1) % face_len == end {
        Some(BoundaryEdgeSource {
            edge_index: start,
            half,
        })
    } else if (end + 1) % face_len == start {
        Some(BoundaryEdgeSource {
            edge_index: end,
            half,
        })
    } else {
        None
    }
}

fn quad_with_edges(
    vertices: [[Ratio; 2]; 4],
    sources: [Option<BoundaryEdgeSource>; 4],
) -> Result<QuadWithProvenance, &'static str> {
    for i in 0..4 {
        if sign(cross(
            vertices[i],
            vertices[(i + 1) % 4],
            vertices[(i + 2) % 4],
        )?) <= 0
        {
            return Err("generated quadrilateral is not strictly convex and counter-clockwise");
        }
    }
    Ok(QuadWithProvenance {
        vertices,
        edge_sources: sources,
    })
}

/// Quadrangulate a simple CCW face and retain the source edge of each boundary half-edge.
pub fn quadrangulate_with_provenance(
    face: &[[Ratio; 2]],
) -> Result<Vec<QuadWithProvenance>, &'static str> {
    validate_face(face)?;
    let triangles = triangulate(face)?;
    let mut output = Vec::with_capacity(3 * triangles.len());
    for [a, b, c] in triangles {
        let tri = [face[a], face[b], face[c]];
        let [mab, mbc, mca] = [
            midpoint(tri[0], tri[1])?,
            midpoint(tri[1], tri[2])?,
            midpoint(tri[2], tri[0])?,
        ];
        let center = centroid(tri)?;
        output.push(quad_with_edges(
            [tri[0], mab, center, mca],
            [
                triangle_edge_source(face.len(), a, b, 0),
                None,
                None,
                triangle_edge_source(face.len(), c, a, 1),
            ],
        )?);
        output.push(quad_with_edges(
            [tri[1], mbc, center, mab],
            [
                triangle_edge_source(face.len(), b, c, 0),
                None,
                None,
                triangle_edge_source(face.len(), a, b, 1),
            ],
        )?);
        output.push(quad_with_edges(
            [tri[2], mca, center, mbc],
            [
                triangle_edge_source(face.len(), c, a, 0),
                None,
                None,
                triangle_edge_source(face.len(), b, c, 1),
            ],
        )?);
    }
    let mut area = Ratio::zero();
    for quad in &output {
        area = area.checked_add(signed_area2(&quad.vertices)?)?;
    }
    if area != signed_area2(face)? {
        return Err("quadrangulation does not preserve the exact oriented area");
    }
    Ok(output)
}

/// Exact deterministic quadrangulation. Each input edge is split at its midpoint.
pub fn quadrangulate(face: &[[Ratio; 2]]) -> Result<Vec<[[Ratio; 2]; 4]>, &'static str> {
    Ok(quadrangulate_with_provenance(face)?
        .into_iter()
        .map(|quad| quad.vertices)
        .collect())
}

/// Sort exact rational points lexicographically without using `Ratio::Ord`.
pub fn compare_points_exact(a: [Ratio; 2], b: [Ratio; 2]) -> Ordering {
    cmp_point(a, b)
}

/// Encode sketch segments on the doubled coordinate grid used by the layer planner.
pub(crate) fn footprint_segments(job: &SketchLayers) -> Vec<Segment> {
    let mut segments = Vec::new();
    for contour in std::iter::once(&job.outer).chain(job.cuts.iter().map(|(ring, _)| ring)) {
        for (index, (a, b)) in edges(contour).enumerate() {
            segments.push(Segment {
                id: format!("boundary/{}/{index}", contour.id),
                a: a.map(|x| 2 * x),
                b: b.map(|x| 2 * x),
                role: Role::Boundary,
            });
        }
    }
    segments.extend(job.guides.iter().map(|g| Segment {
        id: format!("guide/{}", g.id),
        a: g.a.map(|x| 2 * x),
        b: g.b.map(|x| 2 * x),
        role: Role::Construction,
    }));
    segments
}

/// Generate guides in integer ticks, then express the exact arrangement in the
/// doubled millimeter coordinates consumed by the layer builder.
pub(crate) fn scaled_guided_arrangement(
    job: &SketchLayers,
    scale: i32,
    make_guides: impl FnOnce(&SketchLayers) -> Option<Vec<Segment>>,
) -> Result<Option<Arrangement>, &'static str> {
    if !(1..=8).contains(&scale) || !job.guides.is_empty() {
        return Ok(None);
    }
    let mut ticks = job.clone();
    for ring in std::iter::once(&mut ticks.outer).chain(ticks.cuts.iter_mut().map(|c| &mut c.0)) {
        for p in &mut ring.vertices {
            for coordinate in p {
                *coordinate = coordinate
                    .checked_mul(scale)
                    .ok_or("guide_coordinate_overflow")?;
            }
        }
    }
    let Some(guides) = make_guides(&ticks) else {
        return Ok(None);
    };
    ticks.guides = guides;
    let Some(mut graph) = guided_quad_arrangement(&ticks)? else {
        return Ok(None);
    };
    let factor = Ratio::new(1, i128::from(scale))?;
    for p in &mut graph.vertices {
        for coordinate in p {
            *coordinate = coordinate.checked_mul(factor)?;
        }
    }
    Ok(Some(graph))
}

/// Keep complete quadrilateral faces when the sketch itself supplies a conforming
/// subdivision of the stock. Dropped bridge edges and unconnected inner loops
/// must not turn into a spurious face or an unmatched half of a material interface.
pub(crate) fn guided_quad_arrangement(
    job: &SketchLayers,
) -> Result<Option<Arrangement>, &'static str> {
    let graph = arrange(&footprint_segments(job))?;
    let outer = format!("boundary/{}/", job.outer.id);
    let mut uses = vec![[0_usize; 2]; graph.edges.len()];
    for face in &graph.faces {
        if face.vertices.len() != 4 {
            return Ok(None);
        }
        for i in 0..4 {
            let point = |j: usize| graph.vertices[face.vertices[j % 4]];
            if cross(point(i), point(i + 1), point(i + 2))?.numerator() <= 0 {
                return Ok(None);
            }
        }
        for side in &face.edges {
            uses[side.edge][usize::from(side.forward)] += 1;
        }
    }
    for (edge, directions) in graph.edges.iter().zip(uses) {
        let stock_boundary = edge
            .provenance
            .iter()
            .any(|source| source.segment_id.starts_with(&outer));
        let conforming = if stock_boundary {
            directions == [1, 0] || directions == [0, 1]
        } else {
            directions == [1, 1]
        };
        if !conforming {
            return Ok(None);
        }
    }
    Ok(Some(graph))
}

/// Name a guided cell without reference to its rank in the full footprint.
/// In the rectangular pocket template, the wall and the depth event persist
/// when the pocket is moved or resized. Other plans keep unchanged cells by
/// their exact local coordinates, independent of the rest of the graph.
pub(crate) fn guided_cell_id(
    quad: &[[Ratio; 2]; 4],
    layer: [i32; 2],
    job: &SketchLayers,
    method: &str,
) -> String {
    if method == "guided_rectangular_blind_pocket" {
        let corners = &job.outer.vertices;
        let mut side = "footprint";
        for (axis, low_name, high_name) in [(0, "west", "east"), (1, "south", "north")] {
            let values = corners.iter().map(|p| p[axis]);
            for (limit, name) in [
                (values.clone().min().unwrap(), low_name),
                (values.max().unwrap(), high_name),
            ] {
                let value = crate::sketch_layers::ratio(limit);
                if quad.iter().filter(|p| p[axis] == value).count() == 2 {
                    side = name;
                }
            }
        }
        let floor = job.cuts[0].1[0];
        let slab = if layer[1] == floor { "below" } else { "beside" };
        return format!("pocket_{}_{side}_{slab}", job.cuts[0].0.id);
    }
    use std::fmt::Write;
    let mut corners = *quad;
    corners.sort_by(|a, b| cmp_point(*a, *b));
    let mut id = String::from("quad");
    for p in corners {
        for coordinate in p {
            write!(
                id,
                "_{}d{}",
                coordinate.numerator(),
                coordinate.denominator()
            )
            .unwrap();
        }
    }
    write!(id, "_z{}_{}", layer[0], layer[1]).unwrap();
    id
}
