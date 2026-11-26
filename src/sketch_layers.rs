//! Feature-independent partition of a polygonal sketch, projected through axial events.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    constraint_system::Ratio,
    incidence::enumerate_faces,
    sketch_arrangement::{arrange, Role, Segment},
    sketch_quads::quadrangulate,
};
use std::collections::{BTreeMap, BTreeSet};

type Point = [Ratio; 2];
type Quad = [Point; 4];
#[derive(Clone)]
pub struct SketchLoop {
    pub id: String,
    pub vertices: Vec<[i32; 2]>,
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PlannerStrategy {
    #[default]
    Automatic,
    Triangulate,
    RingPlanarRerank,
    RingPlanarLookahead,
}
#[derive(Clone)]
pub struct SketchLayers {
    pub outer: SketchLoop,
    pub cuts: Vec<(SketchLoop, [i32; 2])>,
    pub guides: Vec<Segment>,
    pub height: i32,
    pub force_triangulation: bool,
    /// Ordered add (true) / remove (false) operations, if the recipe carries Boolean history.
    pub operations: Option<Vec<bool>>,
    pub shift: [i32; 2],
    /// Opt-in planner variant. `Automatic` preserves the historical policy.
    pub strategy: PlannerStrategy,
}
pub struct Partition {
    pub blocks: Vec<BracketBlock>,
    pub certificate: BracketCertificate,
    pub volume: Ratio,
    pub planar_faces: usize,
    pub cavities: usize,
    pub rectangular_grid: bool,
    pub layer_counts: Vec<([i32; 2], usize)>,
    pub generated_levels: Vec<i32>,
    pub method: &'static str,
    pub candidate_comparison: Option<(bool, usize, bool, usize)>,
    pub graph: crate::sketch_arrangement::Arrangement,
}

pub(crate) fn ratio(i: i32) -> Ratio {
    Ratio::new(i128::from(i), 1).expect("integer")
}
pub(crate) fn orient(a: [i32; 2], b: [i32; 2], p: [i32; 2]) -> i64 {
    i64::from(b[0] - a[0]) * i64::from(p[1] - a[1])
        - i64::from(b[1] - a[1]) * i64::from(p[0] - a[0])
}
pub(crate) fn edges(l: &SketchLoop) -> impl Iterator<Item = ([i32; 2], [i32; 2])> + '_ {
    l.vertices
        .iter()
        .copied()
        .zip(l.vertices.iter().copied().cycle().skip(1))
}
pub(crate) fn hits(a: [i32; 2], b: [i32; 2], c: [i32; 2], d: [i32; 2]) -> bool {
    (0..2).all(|k| a[k].min(b[k]) <= c[k].max(d[k]) && c[k].min(d[k]) <= a[k].max(b[k]))
        && orient(a, b, c).signum() * orient(a, b, d).signum() <= 0
        && orient(c, d, a).signum() * orient(c, d, b).signum() <= 0
}
pub(crate) fn area2(l: &SketchLoop) -> i64 {
    signed_area2(l).abs()
}
pub(crate) fn signed_area2(l: &SketchLoop) -> i64 {
    edges(l)
        .map(|(a, b)| i64::from(a[0]) * i64::from(b[1]) - i64::from(a[1]) * i64::from(b[0]))
        .sum::<i64>()
}
fn cross(a: Point, b: Point, p: Point) -> Result<Ratio, &'static str> {
    b[0].checked_sub(a[0])?
        .checked_mul(p[1].checked_sub(a[1])?)?
        .checked_sub(
            b[1].checked_sub(a[1])?
                .checked_mul(p[0].checked_sub(a[0])?)?,
        )
}
pub(crate) fn inside(l: &SketchLoop, p: Point) -> Result<bool, &'static str> {
    let mut yes = false;
    for (a, b) in edges(l) {
        let above = |y| p[1].checked_sub(ratio(y)).map(|v| v.numerator() < 0);
        if above(a[1])? != above(b[1])? {
            let s = cross(a.map(ratio), b.map(ratio), p)?.numerator();
            if (s > 0) == (b[1] > a[1]) {
                yes = !yes;
            }
        }
    }
    Ok(yes)
}
fn on_edge(l: &SketchLoop, p: Point) -> Result<bool, &'static str> {
    for (a, b) in edges(l) {
        if cross(a.map(ratio), b.map(ratio), p)? != Ratio::zero() {
            continue;
        }
        let mut in_box = true;
        for k in 0..2 {
            in_box &= p[k].checked_sub(ratio(a[k].min(b[k])))?.numerator() >= 0
                && p[k].checked_sub(ratio(a[k].max(b[k])))?.numerator() <= 0;
        }
        if in_box {
            return Ok(true);
        }
    }
    Ok(false)
}
pub(crate) fn valid_name(s: &str) -> bool {
    !s.is_empty()
        && s.len() <= 48
        && s.bytes()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, b'_' | b'-'))
}
pub(crate) fn valid_point(p: [i32; 2]) -> bool {
    p.iter().all(|v| (-4096..=4096).contains(v))
}
pub(crate) fn validate(job: &SketchLayers) -> Result<(), &'static str> {
    if job
        .operations
        .as_ref()
        .is_some_and(|ops| ops.len() != job.cuts.len())
    {
        return Err("invalid_boolean_sketch_operations");
    }
    if !(1..=4096).contains(&job.height)
        || !valid_point(job.shift)
        || job.cuts.len() > 12
        || job.guides.len() > 64
    {
        return Err("sketch_layer_budget_or_dimensions");
    }
    let loops: Vec<_> = std::iter::once(&job.outer)
        .chain(job.cuts.iter().map(|(l, _)| l))
        .collect();
    let mut ids = BTreeSet::new();
    for (i, l) in loops.iter().enumerate() {
        if !valid_name(&l.id)
            || !ids.insert(&l.id)
            || !(3..=64).contains(&l.vertices.len())
            || !l.vertices.iter().copied().all(valid_point)
        {
            return Err("invalid_sketch_loop");
        }
        if l.vertices.iter().collect::<BTreeSet<_>>().len() != l.vertices.len() || area2(l) == 0 {
            return Err("degenerate_sketch_loop");
        }
        for (j, (a, b)) in edges(l).enumerate() {
            for (k, (c, d)) in edges(l).enumerate().skip(j + 1) {
                if k == j + 1 || (j == 0 && k + 1 == l.vertices.len()) {
                    if orient(a, b, c) == 0 && orient(a, b, d) == 0 {
                        let (p, q, r) = if k == j + 1 { (a, b, d) } else { (c, a, b) };
                        if (0..2)
                            .map(|k| i64::from(q[k] - p[k]) * i64::from(r[k] - q[k]))
                            .sum::<i64>()
                            <= 0
                        {
                            return Err("sketch_loop_backtracks");
                        }
                    }
                } else if hits(a, b, c, d) {
                    return Err("sketch_loop_self_intersection");
                }
            }
        }
        for other in &loops[..i] {
            if edges(l).any(|(a, b)| edges(other).any(|(c, d)| hits(a, b, c, d))) {
                return Err("sketch_contours_touch_or_intersect");
            }
        }
        if i > 0 && !inside(&job.outer, l.vertices[0].map(ratio))? {
            return Err("cut_outside_outer_sketch");
        }
        if i > 1 && job.operations.is_none() {
            for other in &loops[1..i] {
                if inside(other, l.vertices[0].map(ratio))?
                    || inside(l, other.vertices[0].map(ratio))?
                {
                    return Err("nested_cut_requires_boolean_history");
                }
            }
        }
    }
    for (_, z) in &job.cuts {
        if z[0] < 0 || z[0] >= z[1] || z[1] > job.height {
            return Err("invalid_sketch_cut_interval");
        }
    }
    for g in &job.guides {
        if !valid_name(&g.id)
            || !ids.insert(&g.id)
            || !valid_point(g.a)
            || !valid_point(g.b)
            || g.a == g.b
        {
            return Err("invalid_sketch_construction_line");
        }
    }
    for (i, g) in job.guides.iter().enumerate() {
        for p in [g.a, g.b] {
            let attached = std::iter::once(&job.outer)
                .chain(job.cuts.iter().map(|(l, _)| l))
                .any(|l| edges(l).any(|(a, b)| hits(a, b, p, p)))
                || job
                    .guides
                    .iter()
                    .enumerate()
                    .any(|(k, h)| i != k && hits(h.a, h.b, p, p));
            if !attached && inside(&job.outer, p.map(ratio))? {
                return Err("dangling_sketch_construction_endpoint");
            }
        }
    }
    Ok(())
}
fn mean(points: &[Point]) -> Result<Point, &'static str> {
    let mut p = [Ratio::zero(); 2];
    for v in points {
        for k in 0..2 {
            p[k] = p[k].checked_add(v[k])?;
        }
    }
    let s = Ratio::new(1, points.len() as i128)?;
    Ok([p[0].checked_mul(s)?, p[1].checked_mul(s)?])
}
fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}
fn make_block(
    q: &Quad,
    z: [i32; 2],
    job: &SketchLayers,
    id: String,
) -> Result<BracketBlock, &'static str> {
    let mut pts = Vec::new();
    let mut w = 1i128;
    for i in [0, 3, 1, 2] {
        for height in z {
            let s = Ratio::new(i128::from(height), i128::from(job.height))?;
            let p = [
                q[i][0].checked_add(ratio(job.shift[0]).checked_mul(s)?)?,
                q[i][1].checked_add(ratio(job.shift[1]).checked_mul(s)?)?,
                ratio(height),
                Ratio::one(),
            ];
            for v in p {
                w = (w / gcd(w, v.denominator()))
                    .checked_mul(v.denominator())
                    .ok_or("sketch_homogeneous_scale_overflow")?;
            }
            pts.push(p);
        }
    }
    let mut controls = Vec::new();
    for p in pts {
        let mut h = [0i64; 4];
        for k in 0..4 {
            h[k] = i64::try_from(
                p[k].numerator()
                    .checked_mul(w / p[k].denominator())
                    .ok_or("sketch_coordinate_overflow")?,
            )
            .map_err(|_| "sketch_coordinate_range")?;
        }
        controls.push(h);
    }
    Ok(BracketBlock {
        id,
        degrees: [1, 1, 1],
        controls,
    })
}
fn graph(
    job: &SketchLayers,
) -> Result<(crate::sketch_arrangement::Arrangement, Vec<i32>, bool), &'static str> {
    let mut lines = crate::sketch_quads::footprint_segments(job);
    let mut levels = BTreeSet::new();
    for contour in std::iter::once(&job.outer).chain(job.cuts.iter().map(|(l, _)| l)) {
        let ys = contour.vertices.iter().map(|p| p[1]);
        levels.insert(ys.clone().min().unwrap() + ys.max().unwrap());
    }
    let mut bounds = [[i32::MAX, i32::MIN]; 2];
    for s in &lines {
        for p in [s.a, s.b] {
            for k in 0..2 {
                bounds[k][0] = bounds[k][0].min(p[k] - 2);
                bounds[k][1] = bounds[k][1].max(p[k] + 2);
            }
        }
    }
    let [[x0, x1], [y0, y1]] = bounds;
    let rectangular =
        !job.force_triangulation && lines.iter().all(|s| s.a[0] == s.b[0] || s.a[1] == s.b[1]);
    if rectangular {
        levels = lines.iter().flat_map(|s| [s.a[1], s.b[1]]).collect();
        let xs: BTreeSet<_> = lines.iter().flat_map(|s| [s.a[0], s.b[0]]).collect();
        for x in xs {
            lines.push(Segment {
                id: format!("grid/x/{x}"),
                a: [x, y0],
                b: [x, y1],
                role: Role::Construction,
            });
        }
    }
    let frame = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
    for i in 0..4 {
        lines.push(Segment {
            id: format!("frame/{i}"),
            a: frame[i],
            b: frame[(i + 1) % 4],
            role: Role::Construction,
        });
    }
    for &y in &levels {
        lines.push(Segment {
            id: format!("bridge/{y}"),
            a: [x0, y],
            b: [x1, y],
            role: Role::Construction,
        });
    }
    Ok((arrange(&lines)?, levels.into_iter().collect(), rectangular))
}

/// Resolve ordered material history on a single open slab.
fn material(job: &SketchLayers, membership: &[usize], layer: [i32; 2]) -> bool {
    membership
        .iter()
        .rev()
        .find(|&&i| job.cuts[i].1[0] <= layer[0] && layer[1] <= job.cuts[i].1[1])
        .is_none_or(|&i| job.operations.as_ref().is_some_and(|ops| ops[i]))
}

/// Exclusive areas of contour regions, independent of the quadrangulation.
fn zones(job: &SketchLayers) -> Result<Vec<(i64, Vec<usize>)>, &'static str> {
    let mut regions = vec![(area2(&job.outer), Vec::new())];
    for (i, (l, _)) in job.cuts.iter().enumerate() {
        let mut path = Vec::new();
        for (j, (other, _)) in job.cuts.iter().enumerate() {
            if i == j || inside(other, l.vertices[0].map(ratio))? {
                path.push(j);
            }
        }
        regions.push((area2(l), path));
    }
    for i in 0..job.cuts.len() {
        let parent = regions[i + 1]
            .1
            .iter()
            .filter(|&&k| k != i)
            .min_by_key(|&&k| area2(&job.cuts[k].0))
            .map_or(0, |&k| k + 1);
        regions[parent].0 -= area2(&job.cuts[i].0);
    }
    Ok(regions)
}

/// Flood both material and air in the common conforming footprint.
fn connected_voids(
    quads: &[(Quad, Vec<usize>)],
    layers: &[Vec<bool>],
) -> Result<usize, &'static str> {
    let n = quads.len();
    if n == 0 || layers.is_empty() {
        return Err("empty_sketch_material");
    }
    let mut neighbors = vec![Vec::new(); n];
    let mut edges = BTreeMap::new();
    for (i, (q, _)) in quads.iter().enumerate() {
        for k in 0..4 {
            let (a, b) = (q[k], q[(k + 1) % 4]);
            let edge = if a < b { (a, b) } else { (b, a) };
            if let Some(j) = edges.insert(edge, i) {
                neighbors[i].push(j);
                neighbors[j].push(i);
            }
        }
    }
    for pair in layers.windows(2) {
        for (i, row) in neighbors.iter().enumerate() {
            for &j in row {
                if pair[0][i] != pair[0][j] && pair[0][i] == pair[1][j] && pair[0][j] == pair[1][i]
                {
                    return Err("nonmanifold_sketch_layer_transition");
                }
            }
        }
    }
    let mut seen = vec![false; n * layers.len()];
    let mut solids = 0;
    let mut cavities = 0;
    for root in 0..seen.len() {
        if seen[root] {
            continue;
        }
        let state = layers[root / n][root % n];
        let mut exterior = false;
        let mut work = vec![root];
        seen[root] = true;
        while let Some(v) = work.pop() {
            let l = v / n;
            let i = v % n;
            exterior |= l == 0 || l + 1 == layers.len();
            for t in neighbors[i]
                .iter()
                .map(|&j| l * n + j)
                .chain((l > 0).then(|| v - n))
                .chain((l + 1 < layers.len()).then(|| v + n))
            {
                if !seen[t] && layers[t / n][t % n] == state {
                    seen[t] = true;
                    work.push(t);
                }
            }
        }
        if state {
            solids += 1;
        } else if !exterior {
            cavities += 1;
        }
    }
    if solids != 1 {
        return Err("sketch_disconnected_material_components");
    }
    Ok(cavities)
}

/// Build the four corner-to-corner separators for a rectangular blind pocket.
/// The resulting planar graph has four annular quads and one pocket-footprint
/// quad, so it can be propagated unchanged below the pocket and omitted above it.
fn rectangular_blind_pocket_guides(job: &SketchLayers) -> Option<SketchLayers> {
    if job.force_triangulation || !job.guides.is_empty() || job.cuts.len() != 1 {
        return None;
    }
    if job.operations.as_ref().is_some_and(|ops| ops != &[false]) {
        return None;
    }
    let (pocket, z) = &job.cuts[0];
    if z[0] <= 0 || z[1] != job.height {
        return None;
    }

    fn rectangle(loop_: &SketchLoop) -> Option<[[i32; 2]; 4]> {
        if loop_.vertices.len() != 4 {
            return None;
        }
        let xmin = loop_.vertices.iter().map(|p| p[0]).min()?;
        let xmax = loop_.vertices.iter().map(|p| p[0]).max()?;
        let ymin = loop_.vertices.iter().map(|p| p[1]).min()?;
        let ymax = loop_.vertices.iter().map(|p| p[1]).max()?;
        if xmin == xmax || ymin == ymax {
            return None;
        }
        let expected = [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]];
        let actual: BTreeSet<_> = loop_.vertices.iter().copied().collect();
        (actual.len() == 4 && expected.into_iter().all(|p| actual.contains(&p))).then_some(expected)
    }

    let outer = rectangle(&job.outer)?;
    let inner = rectangle(pocket)?;
    if !(outer[0][0] < inner[0][0]
        && inner[1][0] < outer[1][0]
        && outer[0][1] < inner[0][1]
        && inner[2][1] < outer[2][1])
    {
        return None;
    }

    let mut candidate = job.clone();
    for i in 0..4 {
        candidate.guides.push(Segment {
            id: format!("pocket-corner-{i}"),
            a: outer[i],
            b: inner[i],
            role: Role::Construction,
        });
    }
    Some(candidate)
}

fn sample_condition(p: &Partition) -> f64 {
    let mut worst: f64 = 1.0;
    for block in &p.blocks {
        let controls = block.controls.iter().map(|v| v.map(|c| c as f64)).collect();
        let Ok(map) = crate::RationalBezierMap::new(block.degrees, controls) else {
            return f64::INFINITY;
        };
        for u in [0., 0.5, 1.] {
            for v in [0., 0.5, 1.] {
                for w in [0., 0.5, 1.] {
                    let Ok(sample) = map.evaluate([u, v, w]) else {
                        return f64::INFINITY;
                    };
                    let Ok(q) = crate::quality::at_jacobian(sample.jacobian_columns) else {
                        return f64::INFINITY;
                    };
                    worst = worst.max(q.condition_number);
                }
            }
        }
    }
    worst
}

// Block count is primary. Treat roundoff-scale differences as ties so that
// symmetric candidates retain their deterministic proposal order.
fn prefer(candidate: &Partition, current: &Partition, quality_ties: bool) -> bool {
    match candidate.blocks.len().cmp(&current.blocks.len()) {
        std::cmp::Ordering::Less => true,
        std::cmp::Ordering::Greater => false,
        std::cmp::Ordering::Equal if quality_ties => {
            sample_condition(candidate) < sample_condition(current) * (1.0 - 1e-12)
        }
        _ => false,
    }
}

/// Candidate set used by the e5bd5c6 model format.
pub(crate) fn partition_rectangular_hosts(job: &SketchLayers) -> Result<Partition, &'static str> {
    validate(job)?;
    rectangular_host_candidates(job)
}

fn rectangular_host_candidates(job: &SketchLayers) -> Result<Partition, &'static str> {
    let mut best = partition_initial_candidates(job, true);
    for (axis, method) in [(0, "rectangular_hosts_x"), (1, "rectangular_hosts_y")] {
        let Some(guides) = crate::sketch_quad_separators::rectangular_hosts(job, axis) else {
            continue;
        };
        let mut proposed = job.clone();
        proposed.guides = guides;
        let Ok(Some(graph)) = crate::sketch_quads::guided_quad_arrangement(&proposed) else {
            continue;
        };
        let Ok(candidate) = partition_from_graph(job, Some((graph, method))) else {
            continue;
        };
        let replace = match &best {
            Ok(old) => prefer(&candidate, old, true),
            Err(_) => true,
        };
        if replace {
            best = Ok(candidate);
        }
    }
    best
}

pub(crate) fn partition_before_rings(job: &SketchLayers) -> Result<Partition, &'static str> {
    validate(job)?;
    pre_ring_candidates(job)
}

fn pre_ring_candidates(job: &SketchLayers) -> Result<Partition, &'static str> {
    let mut best = rectangular_host_candidates(job);
    for scale in [1, 2, 4, 8] {
        let graph = crate::sketch_quads::scaled_guided_arrangement(job, scale, |scaled| {
            crate::sketch_quad_separators::convex_host_guides(scaled)
        });
        if let Ok(Some(graph)) = graph {
            if let Ok(candidate) =
                partition_from_graph(job, Some((graph, "convex_host_separators")))
            {
                let replace = match &best {
                    Ok(old) => prefer(&candidate, old, true),
                    Err(_) => true,
                };
                if replace {
                    best = Ok(candidate);
                }
            }
        }
    }
    if let Ok(Some(graph)) = crate::sketch_radial::footprint(job) {
        if let Ok(candidate) = partition_from_graph(job, Some((graph, "radial_convex_footprint"))) {
            let replace = match &best {
                Ok(old) => prefer(&candidate, old, true),
                Err(_) => true,
            };
            if replace {
                best = Ok(candidate);
            }
        }
    }
    best
}

pub fn partition(job: &SketchLayers) -> Result<Partition, &'static str> {
    if job.strategy == PlannerStrategy::RingPlanarLookahead {
        return partition_ring_planar_lookahead(job);
    }
    if job.strategy == PlannerStrategy::RingPlanarRerank {
        return partition_ring_planar_rerank(job);
    }
    partition_automatic(job)
}

fn partition_automatic(job: &SketchLayers) -> Result<Partition, &'static str> {
    validate(job)?;
    let mut best = pre_ring_candidates(job);
    crate::sketch_ring_footprint::visit(job, |graph| {
        if let Ok(candidate) =
            partition_from_graph(job, Some((graph, "separated_convex_quad_rings")))
        {
            if best
                .as_ref()
                .map_or(true, |old| prefer(&candidate, old, true))
            {
                best = Ok(candidate);
            }
        }
    });
    best
}

/// Experimental opt-in selector that scores bounded host-separator/ring
/// candidates by their planar quadrilateral quality. Block count remains the
/// primary objective; all accepted candidates pass the same exact partition
/// construction and certificate path as [`partition`].
pub fn partition_ring_planar_lookahead(job: &SketchLayers) -> Result<Partition, &'static str> {
    partition_ring_planar_select(job, true)
}

/// Ablation that re-ranks only the default bounded separation proposals.
pub fn partition_ring_planar_rerank(job: &SketchLayers) -> Result<Partition, &'static str> {
    partition_ring_planar_select(job, false)
}

fn partition_ring_planar_select(
    job: &SketchLayers,
    expanded_endpoints: bool,
) -> Result<Partition, &'static str> {
    validate(job)?;
    let mut best = partition_automatic(job)?;
    let mut best_quality = planar_quality_score(job, &best)?;
    let mut consider_graph =
        |graph| consider_planar_candidate(job, &mut best, &mut best_quality, graph);
    if expanded_endpoints {
        crate::sketch_ring_footprint::visit_lookahead(job, &mut consider_graph);
    } else {
        crate::sketch_ring_footprint::visit(job, &mut consider_graph);
    }
    drop(consider_graph);
    if expanded_endpoints && best.rectangular_grid {
        let graph = best.graph.clone();
        if let Ok(slides) = slide_graphs(job, &graph, 128) {
            for slide in slides {
                consider_planar_candidate(job, &mut best, &mut best_quality, slide);
            }
        }
    }
    Ok(best)
}

fn consider_planar_candidate(
    job: &SketchLayers,
    best: &mut Partition,
    best_quality: &mut f64,
    graph: crate::sketch_arrangement::Arrangement,
) {
    let Ok(candidate) = partition_from_graph(job, Some((graph, "separated_convex_quad_rings")))
    else {
        return;
    };
    let Ok(candidate_quality) = planar_quality_score(job, &candidate) else {
        return;
    };
    let replace = match candidate.blocks.len().cmp(&best.blocks.len()) {
        std::cmp::Ordering::Less => true,
        std::cmp::Ordering::Greater => false,
        std::cmp::Ordering::Equal => candidate_quality < *best_quality * (1.0 - 1e-12),
    };
    if replace {
        *best = candidate;
        *best_quality = candidate_quality;
    }
}

fn point_on_segment(a: Point, b: Point, p: Point) -> Result<bool, &'static str> {
    if cross(a, b, p)? != Ratio::zero() {
        return Ok(false);
    }
    for k in 0..2 {
        let (lo, hi) = if a[k].cmp_value(&b[k]).is_le() {
            (a[k], b[k])
        } else {
            (b[k], a[k])
        };
        if p[k].cmp_value(&lo).is_lt() || p[k].cmp_value(&hi).is_gt() {
            return Ok(false);
        }
    }
    Ok(true)
}

fn source_segments(job: &SketchLayers) -> Vec<(Point, Point)> {
    let graph_point = |point: [i32; 2]| {
        point.map(|value| Ratio::new(i128::from(value) * 2, 1).expect("bounded sketch point"))
    };
    let mut segments = Vec::new();
    for contour in std::iter::once(&job.outer).chain(job.cuts.iter().map(|(loop_, _)| loop_)) {
        segments.extend(edges(contour).map(|(a, b)| (graph_point(a), graph_point(b))));
    }
    segments.extend(
        job.guides
            .iter()
            .map(|guide| (graph_point(guide.a), graph_point(guide.b))),
    );
    segments
}

fn slid_graph_is_valid(
    graph: &crate::sketch_arrangement::Arrangement,
) -> Result<bool, &'static str> {
    for face in &graph.faces {
        if face.vertices.len() != face.edges.len() || face.vertices.len() < 3 {
            return Ok(false);
        }
        let mut twice_area = Ratio::zero();
        for i in 0..face.vertices.len() {
            let a = graph.vertices[face.vertices[i]];
            let b = graph.vertices[face.vertices[(i + 1) % face.vertices.len()]];
            let edge = graph
                .edges
                .get(face.edges[i].edge)
                .ok_or("unknown_arrangement_edge")?;
            let pair = edge.vertices;
            let expected = [
                face.vertices[i],
                face.vertices[(i + 1) % face.vertices.len()],
            ];
            if pair != expected && pair != [expected[1], expected[0]] {
                return Ok(false);
            }
            twice_area = twice_area.checked_add(
                a[0].checked_mul(b[1])?
                    .checked_sub(a[1].checked_mul(b[0])?)?,
            )?;
        }
        if twice_area.numerator() <= 0 {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Make bounded exact alternatives by sliding unpinned straight-through
/// arrangement junctions. Every incident face shares the same moved vertex.
fn slide_graphs(
    job: &SketchLayers,
    graph: &crate::sketch_arrangement::Arrangement,
    budget: usize,
) -> Result<Vec<crate::sketch_arrangement::Arrangement>, &'static str> {
    let segments = source_segments(job);
    let graph_point = |point: [i32; 2]| {
        point.map(|value| Ratio::new(i128::from(value) * 2, 1).expect("bounded sketch point"))
    };
    let pinned: Vec<Point> = std::iter::once(&job.outer)
        .chain(job.cuts.iter().map(|(loop_, _)| loop_))
        .flat_map(|loop_| loop_.vertices.iter().map(|p| graph_point(*p)))
        .chain(
            job.guides
                .iter()
                .flat_map(|guide| [graph_point(guide.a), graph_point(guide.b)]),
        )
        .collect();
    let mut neighbors = vec![BTreeSet::new(); graph.vertices.len()];
    for edge in &graph.edges {
        let [a, b] = edge.vertices;
        neighbors[a].insert(b);
        neighbors[b].insert(a);
    }
    let fractions = [
        Ratio::new(1, 16)?,
        Ratio::new(1, 8)?,
        Ratio::new(1, 4)?,
        Ratio::new(1, 2)?,
    ];
    let mut variants = Vec::new();
    let mut seen = Vec::<(usize, Point)>::new();
    for (vertex, &point) in graph.vertices.iter().enumerate() {
        if pinned.contains(&point) {
            continue;
        }
        let incident: Vec<_> = neighbors[vertex].iter().copied().collect();
        for i in 0..incident.len() {
            for j in i + 1..incident.len() {
                let (a, b) = (graph.vertices[incident[i]], graph.vertices[incident[j]]);
                if a == b || point == a || point == b || !point_on_segment(a, b, point)? {
                    continue;
                }
                for target in [a, b] {
                    for fraction in fractions {
                        let moved = [
                            point[0].checked_add(
                                target[0].checked_sub(point[0])?.checked_mul(fraction)?,
                            )?,
                            point[1].checked_add(
                                target[1].checked_sub(point[1])?.checked_mul(fraction)?,
                            )?,
                        ];
                        if moved == point
                            || graph.vertices.contains(&moved)
                            || seen.contains(&(vertex, moved))
                        {
                            continue;
                        }
                        let mut preserves_sources = true;
                        for &(start, end) in &segments {
                            if point_on_segment(start, end, point)?
                                && !point_on_segment(start, end, moved)?
                            {
                                preserves_sources = false;
                                break;
                            }
                        }
                        if !preserves_sources {
                            continue;
                        }
                        seen.push((vertex, moved));
                        let mut candidate = graph.clone();
                        candidate.vertices[vertex] = moved;
                        if !slid_graph_is_valid(&candidate)? {
                            continue;
                        }
                        variants.push(candidate);
                        if variants.len() >= budget {
                            return Ok(variants);
                        }
                    }
                }
            }
        }
    }
    Ok(variants)
}

pub fn planar_quality_score(
    job: &SketchLayers,
    partition: &Partition,
) -> Result<f64, &'static str> {
    let half = Ratio::new(1, 2)?;
    let mut worst: f64 = 1.0;
    let mut count = 0usize;
    for face in &partition.graph.faces {
        let mut ring = face
            .vertices
            .iter()
            .map(|&i| partition.graph.vertices[i])
            .map(|p| Ok([p[0].checked_mul(half)?, p[1].checked_mul(half)?]))
            .collect::<Result<Vec<_>, &'static str>>()?;
        let cells = if partition.rectangular_grid {
            vec![ring
                .as_slice()
                .try_into()
                .map_err(|_| "nonquadrilateral_ring_face")?]
        } else {
            let first = ring
                .iter()
                .enumerate()
                .min_by_key(|(_, p)| *p)
                .ok_or("empty_ring_face")?
                .0;
            ring.rotate_left(first);
            quadrangulate(&ring)?
        };
        for q in cells {
            let center = mean(&q)?;
            if !inside(&job.outer, center)? {
                continue;
            }
            let mut through_void = false;
            for (cut, interval) in &job.cuts {
                if *interval == [0, job.height] && inside(cut, center)? {
                    through_void = true;
                    break;
                }
            }
            if through_void {
                continue;
            }
            worst = worst.max(crate::sketch_convex_rings::quality(q));
            count += 1;
        }
    }
    if count == 0 {
        return Err("empty_planar_quality_domain");
    }
    Ok(worst)
}

pub(crate) fn partition_classic(job: &SketchLayers) -> Result<Partition, &'static str> {
    validate(job)?;
    partition_initial_candidates(job, false)
}

fn partition_initial_candidates(
    job: &SketchLayers,
    quality_ties: bool,
) -> Result<Partition, &'static str> {
    let Some((candidate_job, candidate_method)) = (!job.force_triangulation)
        .then(|| guided_candidate_job(job))
        .flatten()
    else {
        return partition_from_graph(job, None);
    };
    let Some(candidate_graph) = crate::sketch_quads::guided_quad_arrangement(&candidate_job)
        .ok()
        .flatten()
    else {
        return partition_from_graph(job, None);
    };
    let baseline = partition_from_graph(job, None);
    let guided = partition_from_graph(job, Some((candidate_graph, candidate_method)));
    match (baseline, guided) {
        (Ok(base), Ok(mut candidate)) if prefer(&candidate, &base, quality_ties) => {
            candidate.candidate_comparison =
                Some((true, base.blocks.len(), true, candidate.blocks.len()));
            Ok(candidate)
        }
        (Ok(mut base), Ok(candidate)) => {
            base.candidate_comparison =
                Some((true, base.blocks.len(), true, candidate.blocks.len()));
            Ok(base)
        }
        (Ok(mut base), Err(_)) => {
            base.candidate_comparison = Some((true, base.blocks.len(), false, 0));
            Ok(base)
        }
        (Err(_), Ok(mut candidate)) => {
            candidate.candidate_comparison = Some((false, 0, true, candidate.blocks.len()));
            Ok(candidate)
        }
        (Err(error), Err(_)) => Err(error),
    }
}

/// Replay the pre-planner partition for stored models created before guided
/// quadrilateral footprints were introduced.
pub(crate) fn partition_legacy(job: &SketchLayers) -> Result<Partition, &'static str> {
    validate(job)?;
    partition_from_graph(job, None)
}

fn guided_candidate_job(job: &SketchLayers) -> Option<(SketchLayers, &'static str)> {
    if !job.guides.is_empty() {
        return Some((job.clone(), "guided_user_quad_footprint"));
    }
    if let Some(candidate) = rectangular_blind_pocket_guides(job) {
        return Some((candidate, "guided_rectangular_blind_pocket"));
    }
    if let Some(guides) = crate::sketch_quad_separators::local_quad_guides(job) {
        let mut candidate = job.clone();
        candidate.guides = guides;
        return Some((candidate, "guided_quad_footprint"));
    }
    job.cuts
        .is_empty()
        .then(|| (job.clone(), "guided_quad_footprint"))
}

fn partition_from_graph(
    job: &SketchLayers,
    guided: Option<(crate::sketch_arrangement::Arrangement, &'static str)>,
) -> Result<Partition, &'static str> {
    let guided_method = guided.as_ref().map(|(_, method)| *method);
    let (graph, levels, rectangular_grid) = if let Some((graph, _)) = guided {
        // Every face in this candidate has already been checked as a strict
        // convex quadrilateral, so keep its four source sides intact.
        (graph, Vec::new(), true)
    } else {
        graph(job)?
    };
    let mut quads = Vec::new();
    let half = Ratio::new(1, 2)?;
    for face in &graph.faces {
        let mut ring = face
            .vertices
            .iter()
            .map(|&i| graph.vertices[i])
            .map(|p| Ok([p[0].checked_mul(half)?, p[1].checked_mul(half)?]))
            .collect::<Result<Vec<_>, &'static str>>()?;
        let first = ring
            .iter()
            .enumerate()
            .min_by_key(|(_, p)| *p)
            .ok_or("empty_sketch_face")?
            .0;
        ring.rotate_left(first);
        let cells = if rectangular_grid {
            vec![ring
                .as_slice()
                .try_into()
                .map_err(|_| "nonrectangular_grid_cell")?]
        } else {
            quadrangulate(&ring)?
        };
        for mut q in cells {
            let c = mean(&q)?;
            if !inside(&job.outer, c)? {
                continue;
            }
            let start = q.iter().enumerate().min_by_key(|(_, p)| *p).unwrap().0;
            q.rotate_left(start);
            let mut cut = Vec::new();
            for (i, (l, _)) in job.cuts.iter().enumerate() {
                if inside(l, c)? {
                    cut.push(i);
                }
            }
            quads.push((q, cut));
        }
    }
    quads.sort_by_key(|(q, _)| *q);
    let zs: Vec<_> = job
        .cuts
        .iter()
        .flat_map(|(_, z)| *z)
        .chain([0, job.height])
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let regions = zones(job)?;
    let mut states = Vec::new();
    let mut count = Vec::new();
    let mut blocks = Vec::new();
    let mut locations = BTreeMap::new();
    let mut volume = Ratio::zero();
    let mut bound = f64::INFINITY;
    for pair in zs.windows(2) {
        let layer = [pair[0], pair[1]];
        let dz = layer[1] - layer[0];
        let begin = blocks.len();
        let mut layer_area = Ratio::zero();
        let state: Vec<_> = quads
            .iter()
            .map(|(_, path)| material(job, path, layer))
            .collect();
        for (i, (q, _)) in quads.iter().enumerate() {
            if !state[i] {
                continue;
            }
            for k in 0..4 {
                let d = cross(q[k], q[(k + 1) % 4], q[(k + 3) % 4])?.checked_mul(ratio(dz))?;
                if d.numerator() <= 0 {
                    return Err("sketch_quad_nonpositive_jacobian");
                }
                // Both integer conversions can round; enclose numerator and denominator separately.
                bound = bound.min(
                    ((d.numerator() as f64).next_down() / (d.denominator() as f64).next_up())
                        .next_down(),
                );
            }
            let a = cross(q[0], q[1], q[2])?
                .checked_add(cross(q[0], q[2], q[3])?)?
                .checked_mul(half)?;
            layer_area = layer_area.checked_add(a)?;
            let id = if let Some(method) = guided_method {
                crate::sketch_quads::guided_cell_id(q, layer, job, method)
            } else {
                format!("face{i}_z{}_{}", layer[0], layer[1])
            };
            locations.insert(id.clone(), (i, layer));
            blocks.push(make_block(q, layer, job, id)?);
            if blocks.len() > 10000 {
                return Err("sketch_partition_block_budget");
            }
        }
        let target = regions
            .iter()
            .filter(|(_, path)| material(job, path, layer))
            .map(|(a, _)| *a)
            .sum::<i64>();
        states.push(state);
        if layer_area != Ratio::new(i128::from(target), 2)? {
            return Err("sketch_layer_area_mismatch");
        }
        volume = volume.checked_add(layer_area.checked_mul(ratio(dz))?)?;
        count.push((layer, blocks.len() - begin));
    }
    let cavities = connected_voids(&quads, &states)?;
    let prefix = if job.operations.is_some() {
        "feature"
    } else {
        "cut"
    };
    let incidence = enumerate_faces(&blocks)?;
    let mut tags = Vec::new();
    for side in &incidence.exterior {
        let &(i, z) = locations.get(&side.block).ok_or("unknown_sketch_block")?;
        let (q, cut) = &quads[i];
        let tag = if side.axis == 2 {
            let h = z[usize::from(side.upper)];
            if h == 0 {
                "start".into()
            } else if h == job.height {
                "end".into()
            } else {
                let j = cut
                    .iter()
                    .rev()
                    .find(|&&j| job.cuts[j].1.contains(&h))
                    .ok_or("unclassified_internal_sketch_cap")?;
                let (l, interval) = &job.cuts[*j];
                let end = if h == interval[0] { "lower" } else { "upper" };
                format!("{prefix}/{}/{end}", l.id)
            }
        } else {
            let (a, b) = match (side.axis, side.upper) {
                (0, false) => (0, 3),
                (0, true) => (1, 2),
                (1, false) => (0, 1),
                (1, true) => (3, 2),
                _ => return Err("unknown_sketch_side"),
            };
            let midpoint = mean(&[q[a], q[b]])?;
            if on_edge(&job.outer, midpoint)? {
                format!("outer/{}", job.outer.id)
            } else {
                let mut label = None;
                for (l, interval) in &job.cuts {
                    if interval[0] <= z[0] && z[1] <= interval[1] && on_edge(l, midpoint)? {
                        label = Some(format!("{prefix}/{}/wall", l.id));
                        break;
                    }
                }
                label.ok_or("sketch_exposed_internal_separatrix")?
            }
        };
        tags.push((side.clone(), tag));
    }
    let planar_faces = graph.faces.len();
    Ok(Partition {
        blocks,
        certificate: BracketCertificate {
            incidence,
            design_faces: tags,
            radial_error_upper_mm: 0.,
            determinant_lower_bound_mm3: bound,
        },
        volume,
        planar_faces,
        cavities,
        layer_counts: count,
        generated_levels: levels,
        method: guided_method.unwrap_or(if rectangular_grid {
            "orthogonal_feature_footprint"
        } else {
            "triangulated_feature_footprint"
        }),
        candidate_comparison: None,
        rectangular_grid,
        graph,
    })
}
