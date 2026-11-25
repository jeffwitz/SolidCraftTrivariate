//! Exact planar arrangements of finite sketch segments.
use crate::constraint_system::Ratio;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Role {
    Boundary,
    Construction,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    pub id: String,
    pub a: [i32; 2],
    pub b: [i32; 2],
    pub role: Role,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Provenance {
    pub segment_id: String,
    pub role: Role,
    /// True when the source segment runs from `ArrangementEdge.vertices[0]`
    /// to `ArrangementEdge.vertices[1]`.
    pub forward: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArrangementEdge {
    pub vertices: [usize; 2],
    pub provenance: Vec<Provenance>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceClass {
    Material,
    Void,
    Exterior,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaceEdge {
    pub edge: usize,
    /// True when this face traverses the edge in its canonical vertex order.
    pub forward: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Face {
    /// Indices into `Arrangement.vertices`, ordered counter-clockwise.
    pub vertices: Vec<usize>,
    /// One oriented edge per consecutive vertex pair, including the closing edge.
    pub edges: Vec<FaceEdge>,
    pub class: FaceClass,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Arrangement {
    pub vertices: Vec<[Ratio; 2]>,
    pub edges: Vec<ArrangementEdge>,
    /// Positive bounded faces. The unbounded face has no finite polygon and is omitted.
    pub faces: Vec<Face>,
}

type Point = [Ratio; 2];

#[derive(Clone)]
struct AtomicEdge {
    a: usize,
    b: usize,
    provenance: Vec<Provenance>,
}

#[derive(Clone, Copy)]
struct HalfEdge {
    edge: usize,
    from: usize,
    to: usize,
    forward: bool,
}

fn r(value: i128) -> Ratio {
    Ratio::new(value, 1).expect("unit denominator")
}

fn add(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_add(b)
}

fn sub(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_sub(b)
}

fn mul(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_mul(b)
}

fn div(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    if b.numerator() == 0 {
        return Err("sketch_arrangement_zero_divisor");
    }
    let numerator = a
        .numerator()
        .checked_mul(b.denominator())
        .ok_or("sketch_arrangement_rational_overflow")?;
    let denominator = a
        .denominator()
        .checked_mul(b.numerator())
        .ok_or("sketch_arrangement_rational_overflow")?;
    Ratio::new(numerator, denominator).map_err(|_| "sketch_arrangement_rational_overflow")
}

fn numeric_cmp(a: Ratio, b: Ratio) -> Result<Ordering, &'static str> {
    let left = a
        .numerator()
        .checked_mul(b.denominator())
        .ok_or("sketch_arrangement_rational_overflow")?;
    let right = b
        .numerator()
        .checked_mul(a.denominator())
        .ok_or("sketch_arrangement_rational_overflow")?;
    Ok(left.cmp(&right))
}

fn cross(a: Point, b: Point, c: Point) -> Result<Ratio, &'static str> {
    let ab = [sub(b[0], a[0])?, sub(b[1], a[1])?];
    let ac = [sub(c[0], a[0])?, sub(c[1], a[1])?];
    sub(mul(ab[0], ac[1])?, mul(ab[1], ac[0])?)
}

fn parameter(a: [i32; 2], b: [i32; 2], p: Point) -> Result<Ratio, &'static str> {
    let dx = i64::from(b[0]) - i64::from(a[0]);
    let dy = i64::from(b[1]) - i64::from(a[1]);
    if dx != 0 {
        div(sub(p[0], r(i128::from(a[0])))?, r(i128::from(dx)))
    } else {
        div(sub(p[1], r(i128::from(a[1])))?, r(i128::from(dy)))
    }
}

fn point_at(a: [i32; 2], b: [i32; 2], t: Ratio) -> Result<Point, &'static str> {
    Ok([
        add(
            r(i128::from(a[0])),
            mul(r(i128::from(i64::from(b[0]) - i64::from(a[0]))), t)?,
        )?,
        add(
            r(i128::from(a[1])),
            mul(r(i128::from(i64::from(b[1]) - i64::from(a[1]))), t)?,
        )?,
    ])
}

fn on_segment(p: Point, a: [i32; 2], b: [i32; 2]) -> Result<bool, &'static str> {
    let pa = [r(i128::from(a[0])), r(i128::from(a[1]))];
    let pb = [r(i128::from(b[0])), r(i128::from(b[1]))];
    if cross(pa, pb, p)? != Ratio::zero() {
        return Ok(false);
    }
    let t = parameter(a, b, p)?;
    Ok(numeric_cmp(t, Ratio::zero())? != Ordering::Less
        && numeric_cmp(t, Ratio::one())? != Ordering::Greater)
}

fn add_split(splits: &mut [Vec<Ratio>], index: usize, value: Ratio) {
    if !splits[index].contains(&value) {
        splits[index].push(value);
    }
}

fn segment_intersections(
    left: &Segment,
    right: &Segment,
    i: usize,
    j: usize,
    splits: &mut [Vec<Ratio>],
) -> Result<(), &'static str> {
    let p = left.a.map(|x| r(i128::from(x)));
    let q = right.a.map(|x| r(i128::from(x)));
    let rvec = [
        r(i128::from(i64::from(left.b[0]) - i64::from(left.a[0]))),
        r(i128::from(i64::from(left.b[1]) - i64::from(left.a[1]))),
    ];
    let svec = [
        r(i128::from(i64::from(right.b[0]) - i64::from(right.a[0]))),
        r(i128::from(i64::from(right.b[1]) - i64::from(right.a[1]))),
    ];
    let qp = [sub(q[0], p[0])?, sub(q[1], p[1])?];
    let denominator = sub(mul(rvec[0], svec[1])?, mul(rvec[1], svec[0])?)?;
    if denominator != Ratio::zero() {
        let t = div(
            sub(mul(qp[0], svec[1])?, mul(qp[1], svec[0])?)?,
            denominator,
        )?;
        let u = div(
            sub(mul(qp[0], rvec[1])?, mul(qp[1], rvec[0])?)?,
            denominator,
        )?;
        if numeric_cmp(t, Ratio::zero())? != Ordering::Less
            && numeric_cmp(t, Ratio::one())? != Ordering::Greater
            && numeric_cmp(u, Ratio::zero())? != Ordering::Less
            && numeric_cmp(u, Ratio::one())? != Ordering::Greater
        {
            add_split(splits, i, t);
            add_split(splits, j, u);
        }
        return Ok(());
    }

    if cross(p, [add(p[0], rvec[0])?, add(p[1], rvec[1])?], q)? != Ratio::zero() {
        return Ok(());
    }
    // Collinear overlap is represented once geometrically and retains every
    // source label. Every endpoint on the other segment becomes a split node.
    for endpoint in [right.a, right.b] {
        let point = endpoint.map(|x| r(i128::from(x)));
        if on_segment(point, left.a, left.b)? {
            add_split(splits, i, parameter(left.a, left.b, point)?);
        }
    }
    for endpoint in [left.a, left.b] {
        let point = endpoint.map(|x| r(i128::from(x)));
        if on_segment(point, right.a, right.b)? {
            add_split(splits, j, parameter(right.a, right.b, point)?);
        }
    }
    Ok(())
}

fn angle_half(v: Point) -> usize {
    if v[1] > Ratio::zero() || (v[1] == Ratio::zero() && v[0] >= Ratio::zero()) {
        0
    } else {
        1
    }
}

fn angle_less(origin: Point, left: Point, right: Point) -> Result<bool, &'static str> {
    let l = [sub(left[0], origin[0])?, sub(left[1], origin[1])?];
    let rr = [sub(right[0], origin[0])?, sub(right[1], origin[1])?];
    let hl = angle_half(l);
    let hr = angle_half(rr);
    if hl != hr {
        return Ok(hl < hr);
    }
    let determinant = sub(mul(l[0], rr[1])?, mul(l[1], rr[0])?)?;
    if determinant != Ratio::zero() {
        return Ok(determinant > Ratio::zero());
    }
    let ld = add(mul(l[0], l[0])?, mul(l[1], l[1])?)?;
    let rd = add(mul(rr[0], rr[0])?, mul(rr[1], rr[1])?)?;
    Ok(numeric_cmp(ld, rd)? == Ordering::Less)
}

fn area(vertices: &[usize], points: &[Point]) -> Result<Ratio, &'static str> {
    let mut sum = Ratio::zero();
    for i in 0..vertices.len() {
        let a = points[vertices[i]];
        let b = points[vertices[(i + 1) % vertices.len()]];
        sum = add(sum, sub(mul(a[0], b[1])?, mul(a[1], b[0])?)?)?;
    }
    div(sum, r(2))
}

fn sample_in_face(vertices: &[usize], points: &[Point]) -> Result<Point, &'static str> {
    let mut edges = Vec::with_capacity(vertices.len());
    for i in 0..vertices.len() {
        edges.push((
            points[vertices[i]],
            points[vertices[(i + 1) % vertices.len()]],
        ));
    }
    for (a, b) in &edges {
        let dx = sub(b[0], a[0])?;
        let dy = sub(b[1], a[1])?;
        let abs_x = Ratio::new(
            i128::try_from(dx.numerator().unsigned_abs())
                .map_err(|_| "sketch_arrangement_rational_overflow")?,
            dx.denominator(),
        )?;
        let abs_y = Ratio::new(
            i128::try_from(dy.numerator().unsigned_abs())
                .map_err(|_| "sketch_arrangement_rational_overflow")?,
            dy.denominator(),
        )?;
        let magnitude = add(add(abs_x, abs_y)?, Ratio::one())?;
        let midpoint = [div(add(a[0], b[0])?, r(2))?, div(add(a[1], b[1])?, r(2))?];
        for power in 4..=120 {
            let denominator = 1_i128
                .checked_shl(power)
                .ok_or("sketch_arrangement_rational_overflow")?;
            let epsilon = div(div(Ratio::one(), magnitude)?, r(denominator))?;
            let candidate = [
                sub(midpoint[0], mul(dy, epsilon)?)?,
                add(midpoint[1], mul(dx, epsilon)?)?,
            ];
            let mut inside = false;
            for (left, right) in &edges {
                if (numeric_cmp(left[1], candidate[1])? == Ordering::Greater)
                    != (numeric_cmp(right[1], candidate[1])? == Ordering::Greater)
                {
                    let x = add(
                        left[0],
                        div(
                            mul(sub(candidate[1], left[1])?, sub(right[0], left[0])?)?,
                            sub(right[1], left[1])?,
                        )?,
                    )?;
                    if numeric_cmp(x, candidate[0])? == Ordering::Greater {
                        inside = !inside;
                    }
                }
            }
            if inside {
                return Ok(candidate);
            }
        }
    }
    Err("sketch_arrangement_face_has_no_exact_interior_sample")
}

fn boundary_crossings(
    point: Point,
    edges: &[AtomicEdge],
    points: &[Point],
) -> Result<usize, &'static str> {
    let mut crossings = 0;
    for edge in edges {
        if !edge.provenance.iter().any(|p| p.role == Role::Boundary) {
            continue;
        }
        let a = points[edge.a];
        let b = points[edge.b];
        if (numeric_cmp(a[1], point[1])? == Ordering::Greater)
            != (numeric_cmp(b[1], point[1])? == Ordering::Greater)
        {
            let x = add(
                a[0],
                div(
                    mul(sub(point[1], a[1])?, sub(b[0], a[0])?)?,
                    sub(b[1], a[1])?,
                )?,
            )?;
            if numeric_cmp(x, point[0])? == Ordering::Greater {
                crossings += 1;
            }
        }
    }
    Ok(crossings)
}

fn split_repeated_cycle(cycle: Vec<usize>, halves: &[HalfEdge]) -> Vec<Vec<usize>> {
    let mut positions = BTreeMap::new();
    for (i, &half) in cycle.iter().enumerate() {
        let vertex = halves[half].from;
        if let Some(&j) = positions.get(&vertex) {
            if i > j {
                let inner = cycle[j..i].to_vec();
                let mut outer = cycle[i..].to_vec();
                outer.extend_from_slice(&cycle[..j]);
                let mut result = split_repeated_cycle(inner, halves);
                result.extend(split_repeated_cycle(outer, halves));
                return result;
            }
        } else {
            positions.insert(vertex, i);
        }
    }
    vec![cycle]
}

/// Build an exact planar graph, split all intersections and return positive bounded faces.
pub fn arrange(segments: &[Segment]) -> Result<Arrangement, &'static str> {
    if segments.is_empty() {
        return Err("sketch_arrangement_empty");
    }
    let mut ids = BTreeSet::new();
    for segment in segments {
        if segment.id.is_empty() || !ids.insert(segment.id.as_str()) {
            return Err("sketch_arrangement_invalid_or_duplicate_segment_id");
        }
        if segment.a == segment.b {
            return Err("sketch_arrangement_zero_length_segment");
        }
    }

    let mut splits = vec![vec![Ratio::zero(), Ratio::one()]; segments.len()];
    for i in 0..segments.len() {
        for j in i + 1..segments.len() {
            segment_intersections(&segments[i], &segments[j], i, j, &mut splits)?;
        }
    }

    let mut point_map = BTreeMap::<Point, usize>::new();
    let mut points = Vec::<Point>::new();
    let mut edge_map = BTreeMap::<(usize, usize), Vec<Provenance>>::new();
    for (i, segment) in segments.iter().enumerate() {
        let mut ts = splits[i].clone();
        for a in 1..ts.len() {
            let value = ts[a];
            let mut j = a;
            while j > 0 && numeric_cmp(value, ts[j - 1])? == Ordering::Less {
                ts[j] = ts[j - 1];
                j -= 1;
            }
            ts[j] = value;
        }
        for pair in ts.windows(2) {
            if pair[0] == pair[1] {
                continue;
            }
            let pa = point_at(segment.a, segment.b, pair[0])?;
            let pb = point_at(segment.a, segment.b, pair[1])?;
            let va = *point_map.entry(pa).or_insert_with(|| {
                let id = points.len();
                points.push(pa);
                id
            });
            let vb = *point_map.entry(pb).or_insert_with(|| {
                let id = points.len();
                points.push(pb);
                id
            });
            if va == vb {
                return Err("sketch_arrangement_collapsed_atomic_edge");
            }
            let (key, forward) = if va < vb {
                ((va, vb), true)
            } else {
                ((vb, va), false)
            };
            edge_map.entry(key).or_default().push(Provenance {
                segment_id: segment.id.clone(),
                role: segment.role,
                forward,
            });
        }
    }
    let mut edges = Vec::with_capacity(edge_map.len());
    for ((a, b), mut provenance) in edge_map {
        provenance.sort_by(|x, y| {
            (&x.segment_id, x.role, x.forward).cmp(&(&y.segment_id, y.role, y.forward))
        });
        edges.push(AtomicEdge { a, b, provenance });
    }

    let mut boundary_degree = vec![0_usize; points.len()];
    for edge in &edges {
        if edge.provenance.iter().any(|p| p.role == Role::Boundary) {
            boundary_degree[edge.a] += 1;
            boundary_degree[edge.b] += 1;
        }
    }
    if boundary_degree.iter().any(|degree| degree % 2 != 0) {
        return Err("sketch_arrangement_open_boundary");
    }

    let mut public_edges = Vec::with_capacity(edges.len());
    let mut halves = Vec::with_capacity(edges.len() * 2);
    let mut outgoing = vec![Vec::<usize>::new(); points.len()];
    for (index, edge) in edges.iter().enumerate() {
        public_edges.push(ArrangementEdge {
            vertices: [edge.a, edge.b],
            provenance: edge.provenance.clone(),
        });
        let first = halves.len();
        halves.push(HalfEdge {
            edge: index,
            from: edge.a,
            to: edge.b,
            forward: true,
        });
        halves.push(HalfEdge {
            edge: index,
            from: edge.b,
            to: edge.a,
            forward: false,
        });
        outgoing[edge.a].push(first);
        outgoing[edge.b].push(first + 1);
    }
    for (vertex, star) in outgoing.iter_mut().enumerate() {
        for i in 1..star.len() {
            let value = star[i];
            let mut j = i;
            while j > 0
                && angle_less(
                    points[vertex],
                    points[halves[value].to],
                    points[halves[star[j - 1]].to],
                )?
            {
                star[j] = star[j - 1];
                j -= 1;
            }
            star[j] = value;
        }
    }

    let mut next = vec![usize::MAX; halves.len()];
    for (index, half) in halves.iter().enumerate() {
        let reverse = index ^ 1;
        let star = &outgoing[half.to];
        let position = star
            .iter()
            .position(|&candidate| candidate == reverse)
            .ok_or("sketch_arrangement_invalid_half_edge_rotation")?;
        next[index] = star[(position + star.len() - 1) % star.len()];
    }

    let mut visited = vec![false; halves.len()];
    let mut faces = Vec::new();
    for start in 0..halves.len() {
        if visited[start] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut current = start;
        loop {
            if visited[current] {
                if current != start {
                    return Err("sketch_arrangement_nonmanifold_face_walk");
                }
                break;
            }
            visited[current] = true;
            cycle.push(current);
            current = next[current];
        }
        for simple in split_repeated_cycle(cycle, &halves) {
            let vertices = simple.iter().map(|&h| halves[h].from).collect::<Vec<_>>();
            if vertices.len() < 3
                || vertices.iter().copied().collect::<BTreeSet<_>>().len() != vertices.len()
            {
                continue;
            }
            let signed_area = area(&vertices, &points)?;
            if signed_area <= Ratio::zero() {
                continue;
            }
            let sample = sample_in_face(&vertices, &points)?;
            let crossings = boundary_crossings(sample, &edges, &points)?;
            let class = match crossings {
                0 => FaceClass::Exterior,
                n if n % 2 == 1 => FaceClass::Material,
                _ => FaceClass::Void,
            };
            let boundary_edges = simple
                .iter()
                .map(|&h| FaceEdge {
                    edge: halves[h].edge,
                    forward: halves[h].forward,
                })
                .collect();
            faces.push(Face {
                vertices,
                edges: boundary_edges,
                class,
            });
        }
    }
    faces.sort_by(|a, b| {
        a.vertices
            .iter()
            .min()
            .cmp(&b.vertices.iter().min())
            .then_with(|| a.vertices.len().cmp(&b.vertices.len()))
    });
    if faces.is_empty() {
        return Err("sketch_arrangement_no_bounded_faces");
    }
    Ok(Arrangement {
        vertices: points,
        edges: public_edges,
        faces,
    })
}
