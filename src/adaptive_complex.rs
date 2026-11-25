//! Orientable cubical chain complex for an adaptive multipatch Q2 partition.

use crate::{
    adaptive_partition::{AdaptivePartition, FaceRect, LeafKey, LeafSide},
    adaptive_space::AdaptiveSpace,
    incidence::{BlockSide, Orientation},
    json::{object, Json},
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ChartPoint {
    pub root: String,
    pub coordinates: [u32; 3],
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SignedId {
    pub id: usize,
    pub sign: i8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComplexVertex {
    pub id: usize,
    pub charts: Vec<ChartPoint>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ChartSegment {
    pub root: String,
    pub start: [u32; 3],
    pub end: [u32; 3],
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComplexEdge {
    pub id: usize,
    /// The first vertex has coefficient -1 and the second +1 in ∂₁.
    pub vertices: [usize; 2],
    pub provenance: Vec<ChartSegment>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FaceSource {
    Interface {
        first: LeafSide,
        second: LeafSide,
        first_rect: FaceRect,
        second_rect: FaceRect,
        orientation: Orientation,
    },
    Exterior {
        side: LeafSide,
        rect: FaceRect,
        features: Vec<String>,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComplexFace {
    pub id: usize,
    pub boundary: Vec<SignedId>,
    pub incident_cells: Vec<SignedId>,
    pub source: FaceSource,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComplexCell {
    pub id: usize,
    pub leaf: LeafKey,
    pub block_id: String,
    pub boundary: Vec<SignedId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrientedComplex {
    pub grid_denominator: u32,
    pub vertices: Vec<ComplexVertex>,
    pub edges: Vec<ComplexEdge>,
    pub faces: Vec<ComplexFace>,
    pub cells: Vec<ComplexCell>,
    pub boundary_faces: Vec<usize>,
}

#[derive(Debug, Clone)]
struct RawFace {
    first: LeafSide,
    first_rect: FaceRect,
    second: Option<(LeafSide, FaceRect, Orientation)>,
    source: FaceSource,
}

#[derive(Debug, Clone)]
struct RawEdge {
    provenance: BTreeSet<ChartSegment>,
}

#[derive(Debug, Clone)]
struct FaceCycle {
    source: FaceSource,
    boundary: Vec<EdgeUse>,
    incident_cells: Vec<SignedId>,
}

type EdgeUse = ((usize, usize), ChartSegment, i8);

struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new() -> Self {
        Self { parent: Vec::new() }
    }
    fn add(&mut self) -> usize {
        let id = self.parent.len();
        self.parent.push(id);
        id
    }
    fn find(&mut self, value: usize) -> usize {
        if self.parent[value] != value {
            self.parent[value] = self.find(self.parent[value]);
        }
        self.parent[value]
    }
    fn join(&mut self, left: usize, right: usize) {
        let (left, right) = (self.find(left), self.find(right));
        if left != right {
            let (low, high) = if left < right {
                (left, right)
            } else {
                (right, left)
            };
            self.parent[high] = low;
        }
    }
}

fn tangent_axes(axis: usize) -> [usize; 2] {
    match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}

fn side_point(side: &LeafSide, uv: [u32; 2], denominator: u32) -> ChartPoint {
    let [a, b] = tangent_axes(side.root_side.axis);
    let mut coordinates = [0; 3];
    let width = denominator >> side.leaf.level;
    coordinates[side.root_side.axis] =
        side.leaf.index[side.root_side.axis] * width + u32::from(side.root_side.upper) * width;
    coordinates[a] = uv[0];
    coordinates[b] = uv[1];
    ChartPoint {
        root: side.leaf.root.clone(),
        coordinates,
    }
}

fn corners(rect: FaceRect) -> [[u32; 2]; 4] {
    let [[u0, u1], [v0, v1]] = rect.ranges;
    [[u0, v0], [u1, v0], [u1, v1], [u0, v1]]
}

fn flip(value: u32, denominator: u32, reverse: bool) -> u32 {
    if reverse {
        denominator - value
    } else {
        value
    }
}

fn first_to_second(point: [u32; 2], denominator: u32, orientation: Orientation) -> [u32; 2] {
    let a = flip(point[0], denominator, orientation.flip_first);
    let b = flip(point[1], denominator, orientation.flip_second);
    if orientation.swap {
        [b, a]
    } else {
        [a, b]
    }
}

fn second_to_first(point: [u32; 2], denominator: u32, orientation: Orientation) -> [u32; 2] {
    if orientation.swap {
        [
            flip(point[1], denominator, orientation.flip_first),
            flip(point[0], denominator, orientation.flip_second),
        ]
    } else {
        [
            flip(point[0], denominator, orientation.flip_first),
            flip(point[1], denominator, orientation.flip_second),
        ]
    }
}

fn side_orientation(side: &BlockSide) -> i8 {
    let axis_sign = if side.axis.is_multiple_of(2) { 1 } else { -1 };
    axis_sign * if side.upper { 1 } else { -1 }
}

fn add_chart_point(
    point: ChartPoint,
    indices: &mut BTreeMap<ChartPoint, usize>,
    union: &mut UnionFind,
) -> usize {
    if let Some(&index) = indices.get(&point) {
        index
    } else {
        let index = union.add();
        indices.insert(point, index);
        index
    }
}

fn make_raw_faces(
    adaptive: &AdaptiveSpace,
    partition: &AdaptivePartition,
) -> Result<Vec<RawFace>, &'static str> {
    let report = partition.snapshot()?;
    let mut feature_by_leaf_side = BTreeMap::<LeafSide, Vec<String>>::new();
    for (feature, sides) in &adaptive.boundary_patches {
        for side in sides {
            feature_by_leaf_side
                .entry(side.clone())
                .or_default()
                .push(feature.clone());
        }
    }
    for features in feature_by_leaf_side.values_mut() {
        features.sort();
        features.dedup();
    }
    let mut result = Vec::with_capacity(report.interfaces.len() + report.exterior.len());
    for patch in report.interfaces {
        result.push(RawFace {
            first: patch.first.clone(),
            first_rect: patch.first_subrect,
            second: Some((
                patch.second.clone(),
                patch.second_subrect,
                patch.orientation,
            )),
            source: FaceSource::Interface {
                first: patch.first,
                second: patch.second,
                first_rect: patch.first_subrect,
                second_rect: patch.second_subrect,
                orientation: patch.orientation,
            },
        });
    }
    for patch in report.exterior {
        let features = feature_by_leaf_side
            .get(&patch.side)
            .cloned()
            .unwrap_or_default();
        result.push(RawFace {
            first: patch.side.clone(),
            first_rect: patch.root_rect,
            second: None,
            source: FaceSource::Exterior {
                side: patch.side,
                rect: patch.root_rect,
                features,
            },
        });
    }
    Ok(result)
}

fn edge_vertices(
    start: [u32; 2],
    end: [u32; 2],
    side: &LeafSide,
    points: &BTreeMap<ChartPoint, usize>,
    global: &BTreeMap<ChartPoint, usize>,
    denominator: u32,
) -> Result<Vec<(usize, [u32; 3])>, &'static str> {
    let varying = if start[0] != end[0] { 0 } else { 1 };
    if (start[0] == end[0]) == (start[1] == end[1]) {
        return Err("degenerate or non-axis-aligned face edge");
    }
    let fixed = 1 - varying;
    let low = start[varying].min(end[varying]);
    let high = start[varying].max(end[varying]);
    let mut candidates = points
        .keys()
        .filter_map(|point| {
            if point.root != side.leaf.root {
                return None;
            }
            let [a, b] = tangent_axes(side.root_side.axis);
            let uv = [point.coordinates[a], point.coordinates[b]];
            let width = denominator >> side.leaf.level;
            let plane = side.leaf.index[side.root_side.axis] * width
                + u32::from(side.root_side.upper) * width;
            let normal_matches = point.coordinates[side.root_side.axis] == plane;
            (normal_matches && uv[fixed] == start[fixed] && (low..=high).contains(&uv[varying]))
                .then(|| (uv[varying], global[point], point.coordinates))
        })
        .collect::<Vec<_>>();
    candidates.sort_by_key(|entry| entry.0);
    candidates.dedup_by_key(|entry| entry.0);
    if start[varying] > end[varying] {
        candidates.reverse();
    }
    let ids = candidates
        .into_iter()
        .map(|entry| (entry.1, entry.2))
        .collect::<Vec<_>>();
    if ids.len() < 2 || ids.first().map(|entry| entry.0) == ids.last().map(|entry| entry.0) {
        return Err("open or collapsed face boundary segment");
    }
    Ok(ids)
}

fn face_cycle(
    side: &LeafSide,
    rect: FaceRect,
    points: &BTreeMap<ChartPoint, usize>,
    global: &BTreeMap<ChartPoint, usize>,
    edge_ids: &mut BTreeMap<(usize, usize), RawEdge>,
    denominator: u32,
) -> Result<Vec<EdgeUse>, &'static str> {
    let direction = side_orientation(&side.root_side);
    let mut boundary = corners(rect);
    if direction < 0 {
        boundary = [boundary[0], boundary[3], boundary[2], boundary[1]];
    }
    let mut cycle = Vec::new();
    for index in 0..4 {
        let start = boundary[index];
        let end = boundary[(index + 1) % 4];
        let chain = edge_vertices(start, end, side, points, global, denominator)?;
        for pair in chain.windows(2) {
            if pair[0].0 == pair[1].0 {
                return Err("collapsed oriented edge");
            }
            let key = (pair[0].0.min(pair[1].0), pair[0].0.max(pair[1].0));
            let sign = if pair[0].0 < pair[1].0 { 1 } else { -1 };
            let segment = ChartSegment {
                root: side.leaf.root.clone(),
                start: if pair[0].1 <= pair[1].1 {
                    pair[0].1
                } else {
                    pair[1].1
                },
                end: if pair[0].1 <= pair[1].1 {
                    pair[1].1
                } else {
                    pair[0].1
                },
            };
            edge_ids
                .entry(key)
                .or_insert_with(|| RawEdge {
                    provenance: BTreeSet::new(),
                })
                .provenance
                .insert(segment.clone());
            cycle.push((key, segment, sign));
        }
    }
    Ok(cycle)
}

fn cell_id_map(
    report_leaves: &[LeafKey],
    adaptive: &AdaptiveSpace,
) -> Result<(Vec<ComplexCell>, BTreeMap<LeafKey, usize>), &'static str> {
    let mut cells = Vec::with_capacity(report_leaves.len());
    let mut ids = BTreeMap::new();
    for leaf in report_leaves {
        let id = cells.len();
        let block_id = adaptive
            .leaf_blocks
            .get(leaf)
            .ok_or("leaf complex cell has no analysis block")?
            .clone();
        if ids.insert(leaf.clone(), id).is_some() {
            return Err("duplicate adaptive leaf in cell complex");
        }
        cells.push(ComplexCell {
            id,
            leaf: leaf.clone(),
            block_id,
            boundary: Vec::new(),
        });
    }
    Ok((cells, ids))
}

fn validate_chains(complex: &OrientedComplex) -> Result<(), &'static str> {
    for face in &complex.faces {
        let mut vertex_sum = BTreeMap::<usize, i32>::new();
        for incidence in &face.boundary {
            let edge = complex
                .edges
                .get(incidence.id)
                .ok_or("face references missing edge")?;
            *vertex_sum.entry(edge.vertices[0]).or_default() -= i32::from(incidence.sign);
            *vertex_sum.entry(edge.vertices[1]).or_default() += i32::from(incidence.sign);
        }
        if vertex_sum.values().any(|&value| value != 0) {
            return Err("boundary of a 2-cell is not a closed cycle");
        }
    }
    for cell in &complex.cells {
        let mut edge_sum = BTreeMap::<usize, i32>::new();
        for incidence in &cell.boundary {
            let face = complex
                .faces
                .get(incidence.id)
                .ok_or("cell references missing face")?;
            for edge in &face.boundary {
                *edge_sum.entry(edge.id).or_default() +=
                    i32::from(incidence.sign) * i32::from(edge.sign);
            }
        }
        if edge_sum.values().any(|&value| value != 0) {
            return Err("boundary of a 3-cell is not a closed 2-chain");
        }
    }
    let mut face_usage = BTreeMap::<usize, Vec<i8>>::new();
    for cell in &complex.cells {
        for incidence in &cell.boundary {
            face_usage
                .entry(incidence.id)
                .or_default()
                .push(incidence.sign);
        }
    }
    for face in &complex.faces {
        let uses = face_usage.get(&face.id).ok_or("oriented face is unused")?;
        let expected = if complex.boundary_faces.contains(&face.id) {
            1
        } else {
            2
        };
        if uses.len() != expected
            || (expected == 1 && uses[0] != 1)
            || (expected == 2
                && (uses[0] == uses[1] || uses.iter().map(|x| i16::from(*x)).sum::<i16>() != 0))
        {
            return Err("oriented face has invalid cell incidence");
        }
    }
    let mut boundary_edge_usage = BTreeMap::<usize, Vec<i8>>::new();
    for &face_id in &complex.boundary_faces {
        for edge in &complex.faces[face_id].boundary {
            boundary_edge_usage
                .entry(edge.id)
                .or_default()
                .push(edge.sign);
        }
    }
    if boundary_edge_usage.values().any(|uses| {
        uses.len() != 2 || uses[0] == uses[1] || i16::from(uses[0]) + i16::from(uses[1]) != 0
    }) {
        return Err("exterior boundary is not a closed oriented surface");
    }
    Ok(())
}

fn segment_provenance_connected(
    provenance: &BTreeSet<ChartSegment>,
    glued_segments: &BTreeSet<(ChartSegment, ChartSegment)>,
) -> bool {
    if provenance.is_empty() {
        return false;
    }
    let mut reached = BTreeSet::from([provenance.iter().next().unwrap().clone()]);
    loop {
        let before = reached.len();
        for (left, right) in glued_segments {
            if provenance.contains(left) && provenance.contains(right) {
                if reached.contains(left) {
                    reached.insert(right.clone());
                }
                if reached.contains(right) {
                    reached.insert(left.clone());
                }
            }
        }
        if reached.len() == before {
            return reached.len() == provenance.len();
        }
    }
}

/// Construct and verify the signed cellular boundary operators through C₃.
pub fn build_oriented_complex(
    adaptive: &AdaptiveSpace,
    partition: &AdaptivePartition,
) -> Result<OrientedComplex, &'static str> {
    let report = partition.snapshot()?;
    if adaptive.space.blocks.len() != report.leaves.len()
        || adaptive.leaf_blocks.len() != report.leaves.len()
    {
        return Err("adaptive cells and leaf partition disagree");
    }
    let raw_faces = make_raw_faces(adaptive, partition)?;
    let denominator = partition.grid_denominator();
    let mut chart_indices = BTreeMap::<ChartPoint, usize>::new();
    let mut union = UnionFind::new();
    for face in &raw_faces {
        for uv in corners(face.first_rect) {
            add_chart_point(
                side_point(&face.first, uv, denominator),
                &mut chart_indices,
                &mut union,
            );
        }
        if let Some((second, rect, _)) = &face.second {
            for uv in corners(*rect) {
                add_chart_point(
                    side_point(second, uv, denominator),
                    &mut chart_indices,
                    &mut union,
                );
            }
        }
    }
    for face in &raw_faces {
        let Some((second, second_rect, orientation)) = &face.second else {
            continue;
        };
        let [a, b] = tangent_axes(face.first.root_side.axis);
        let width = denominator >> face.first.leaf.level;
        let normal = face.first.leaf.index[face.first.root_side.axis] * width
            + u32::from(face.first.root_side.upper) * width;
        let first_nodes = chart_indices
            .keys()
            .filter(|point| {
                if point.root != face.first.leaf.root
                    || point.coordinates[face.first.root_side.axis] != normal
                {
                    return false;
                }
                let uv = [point.coordinates[a], point.coordinates[b]];
                (0..2).all(|axis| {
                    face.first_rect.ranges[axis][0] <= uv[axis]
                        && uv[axis] <= face.first_rect.ranges[axis][1]
                })
            })
            .cloned()
            .collect::<Vec<_>>();
        let [sa, sb] = tangent_axes(second.root_side.axis);
        let second_width = denominator >> second.leaf.level;
        let second_normal = second.leaf.index[second.root_side.axis] * second_width
            + u32::from(second.root_side.upper) * second_width;
        for point in first_nodes {
            let uv = [point.coordinates[a], point.coordinates[b]];
            let mapped = first_to_second(uv, denominator, *orientation);
            if !(0..2).all(|axis| {
                second_rect.ranges[axis][0] <= mapped[axis]
                    && mapped[axis] <= second_rect.ranges[axis][1]
            }) {
                continue;
            }
            let first_id = chart_indices[&point];
            let mut coordinates = [0; 3];
            coordinates[second.root_side.axis] = second_normal;
            coordinates[sa] = mapped[0];
            coordinates[sb] = mapped[1];
            let second_point = ChartPoint {
                root: second.leaf.root.clone(),
                coordinates,
            };
            let second_id = add_chart_point(second_point, &mut chart_indices, &mut union);
            union.join(first_id, second_id);
        }
        let second_nodes = chart_indices
            .keys()
            .filter(|point| {
                if point.root != second.leaf.root
                    || point.coordinates[second.root_side.axis] != second_normal
                {
                    return false;
                }
                let uv = [point.coordinates[sa], point.coordinates[sb]];
                (0..2).all(|axis| {
                    second_rect.ranges[axis][0] <= uv[axis]
                        && uv[axis] <= second_rect.ranges[axis][1]
                })
            })
            .cloned()
            .collect::<Vec<_>>();
        for point in second_nodes {
            let uv = [point.coordinates[sa], point.coordinates[sb]];
            let mapped = second_to_first(uv, denominator, *orientation);
            if !(0..2).all(|axis| {
                face.first_rect.ranges[axis][0] <= mapped[axis]
                    && mapped[axis] <= face.first_rect.ranges[axis][1]
            }) {
                continue;
            }
            let second_id = chart_indices[&point];
            let mut coordinates = [0; 3];
            coordinates[face.first.root_side.axis] = normal;
            coordinates[a] = mapped[0];
            coordinates[b] = mapped[1];
            let first_point = ChartPoint {
                root: face.first.leaf.root.clone(),
                coordinates,
            };
            let first_id = add_chart_point(first_point, &mut chart_indices, &mut union);
            union.join(first_id, second_id);
        }
    }

    let mut minimum_chart = BTreeMap::<usize, ChartPoint>::new();
    let mut representative_by_chart = BTreeMap::new();
    for (point, &index) in &chart_indices {
        let representative = union.find(index);
        representative_by_chart.insert(point.clone(), representative);
        minimum_chart
            .entry(representative)
            .and_modify(|old| {
                if point < old {
                    *old = point.clone();
                }
            })
            .or_insert_with(|| point.clone());
    }
    let mut ordered_representatives = minimum_chart.into_iter().collect::<Vec<_>>();
    ordered_representatives.sort_by(|a, b| a.1.cmp(&b.1));
    let representative_ids = ordered_representatives
        .iter()
        .enumerate()
        .map(|(id, (representative, _))| (*representative, id))
        .collect::<BTreeMap<_, _>>();
    let global_vertex_by_chart = representative_by_chart
        .iter()
        .map(|(point, representative)| (point.clone(), representative_ids[representative]))
        .collect::<BTreeMap<_, _>>();
    let mut charts_by_vertex = vec![Vec::new(); representative_ids.len()];
    for (point, &id) in &global_vertex_by_chart {
        charts_by_vertex[id].push(point.clone());
    }
    let vertices = charts_by_vertex
        .into_iter()
        .enumerate()
        .map(|(id, mut charts)| {
            charts.sort();
            charts.dedup();
            ComplexVertex { id, charts }
        })
        .collect::<Vec<_>>();

    let mut raw_edges = BTreeMap::<(usize, usize), RawEdge>::new();
    let mut glued_segments = BTreeSet::<(ChartSegment, ChartSegment)>::new();
    let mut face_cycles = Vec::with_capacity(raw_faces.len());
    let mut boundary_face_ids = Vec::new();
    let (_, cell_ids) = cell_id_map(&report.leaves, adaptive)?;
    for face in raw_faces {
        let first_cycle = face_cycle(
            &face.first,
            face.first_rect,
            &chart_indices,
            &global_vertex_by_chart,
            &mut raw_edges,
            denominator,
        )?;
        let first_leaf = cell_ids[&face.first.leaf];
        let incident_cells = if let Some((second, second_rect, _)) = &face.second {
            let second_cycle = face_cycle(
                second,
                *second_rect,
                &chart_indices,
                &global_vertex_by_chart,
                &mut raw_edges,
                denominator,
            )?;
            let first_chain = first_cycle
                .iter()
                .map(|(edge, segment, sign)| (*edge, (segment.clone(), *sign)))
                .collect::<BTreeMap<_, _>>();
            let second_chain = second_cycle
                .iter()
                .map(|(edge, segment, sign)| (*edge, (segment.clone(), *sign)))
                .collect::<BTreeMap<_, _>>();
            if first_chain.len() != first_cycle.len()
                || second_chain.len() != second_cycle.len()
                || first_chain.len() != second_chain.len()
                || first_chain.iter().any(|(edge, (_, sign))| {
                    second_chain
                        .get(edge)
                        .is_none_or(|(_, other_sign)| *other_sign != -*sign)
                })
            {
                return Err("D4 interface face cycles are not oppositely oriented");
            }
            for (edge, (segment, _)) in &first_chain {
                let (other, _) = &second_chain[edge];
                let pair = if segment < other {
                    (segment.clone(), other.clone())
                } else {
                    (other.clone(), segment.clone())
                };
                glued_segments.insert(pair);
            }
            vec![
                SignedId {
                    id: first_leaf,
                    sign: 1,
                },
                SignedId {
                    id: cell_ids[&second.leaf],
                    sign: -1,
                },
            ]
        } else {
            boundary_face_ids.push(face_cycles.len());
            vec![SignedId {
                id: first_leaf,
                sign: 1,
            }]
        };
        face_cycles.push(FaceCycle {
            source: face.source,
            boundary: first_cycle,
            incident_cells,
        });
    }

    for edge in raw_edges.values() {
        if edge.provenance.len() < 2 {
            continue;
        }
        if !segment_provenance_connected(&edge.provenance, &glued_segments) {
            return Err("ambiguous multiedge: distinct chart segments share endpoints without D4 identification");
        }
    }

    let edges = raw_edges
        .into_iter()
        .enumerate()
        .map(|(id, ((first, second), edge))| ComplexEdge {
            id,
            vertices: [first, second],
            provenance: edge.provenance.into_iter().collect(),
        })
        .collect::<Vec<_>>();
    let edge_ids = edges
        .iter()
        .map(|edge| ((edge.vertices[0], edge.vertices[1]), edge.id))
        .collect::<BTreeMap<_, _>>();
    let faces = face_cycles
        .into_iter()
        .enumerate()
        .map(|(id, cycle)| ComplexFace {
            id,
            boundary: cycle
                .boundary
                .into_iter()
                .map(|(edge, _, sign)| SignedId {
                    id: edge_ids[&edge],
                    sign,
                })
                .collect(),
            incident_cells: cycle.incident_cells,
            source: cycle.source,
        })
        .collect::<Vec<_>>();
    let mut cells = cell_id_map(&report.leaves, adaptive)?.0;
    for face in &faces {
        for incidence in &face.incident_cells {
            cells[incidence.id].boundary.push(SignedId {
                id: face.id,
                sign: incidence.sign,
            });
        }
    }
    for cell in &mut cells {
        cell.boundary.sort_by_key(|face| face.id);
    }
    let complex = OrientedComplex {
        grid_denominator: denominator,
        vertices,
        edges,
        faces,
        cells,
        boundary_faces: boundary_face_ids,
    };
    validate_chains(&complex)?;
    Ok(complex)
}

fn text(value: &str) -> Json {
    Json::String(value.to_owned())
}
fn number(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}
fn signed_ids(values: &[SignedId]) -> Json {
    Json::Array(
        values
            .iter()
            .map(|value| Json::Array(vec![number(value.id), number(value.sign)]))
            .collect(),
    )
}
fn leaf_json(leaf: &LeafKey) -> Json {
    object([
        ("root", text(&leaf.root)),
        ("level", number(leaf.level)),
        (
            "index",
            Json::Array(leaf.index.iter().copied().map(number).collect()),
        ),
    ])
}
fn side_json(side: &LeafSide) -> Json {
    object([
        ("root_block", text(&side.root_side.block)),
        ("axis", number(side.root_side.axis)),
        ("upper", Json::Bool(side.root_side.upper)),
        ("leaf", leaf_json(&side.leaf)),
    ])
}
fn rect_json(rect: FaceRect) -> Json {
    Json::Array(
        rect.ranges
            .iter()
            .map(|range| Json::Array(range.iter().copied().map(number).collect()))
            .collect(),
    )
}

impl OrientedComplex {
    pub fn to_json(&self) -> Json {
        let vertices = self
            .vertices
            .iter()
            .map(|vertex| {
                object([
                    ("id", number(vertex.id)),
                    (
                        "charts",
                        Json::Array(
                            vertex
                                .charts
                                .iter()
                                .map(|point| {
                                    object([
                                        ("root", text(&point.root)),
                                        (
                                            "coordinates",
                                            Json::Array(
                                                point
                                                    .coordinates
                                                    .iter()
                                                    .copied()
                                                    .map(number)
                                                    .collect(),
                                            ),
                                        ),
                                    ])
                                })
                                .collect(),
                        ),
                    ),
                ])
            })
            .collect();
        let edges = self
            .edges
            .iter()
            .map(|edge| {
                object([
                    ("id", number(edge.id)),
                    (
                        "vertices",
                        Json::Array(edge.vertices.iter().copied().map(number).collect()),
                    ),
                    (
                        "boundary",
                        Json::Array(vec![
                            Json::Array(vec![number(edge.vertices[0]), number(-1)]),
                            Json::Array(vec![number(edge.vertices[1]), number(1)]),
                        ]),
                    ),
                    (
                        "provenance",
                        Json::Array(
                            edge.provenance
                                .iter()
                                .map(|segment| {
                                    object([
                                        ("root", text(&segment.root)),
                                        (
                                            "start",
                                            Json::Array(
                                                segment.start.iter().copied().map(number).collect(),
                                            ),
                                        ),
                                        (
                                            "end",
                                            Json::Array(
                                                segment.end.iter().copied().map(number).collect(),
                                            ),
                                        ),
                                    ])
                                })
                                .collect(),
                        ),
                    ),
                ])
            })
            .collect();
        let faces = self
            .faces
            .iter()
            .map(|face| {
                let (kind, source) = match &face.source {
                    FaceSource::Interface {
                        first,
                        second,
                        first_rect,
                        second_rect,
                        orientation,
                    } => (
                        "interface",
                        object([
                            ("first", side_json(first)),
                            ("second", side_json(second)),
                            ("first_rect", rect_json(*first_rect)),
                            ("second_rect", rect_json(*second_rect)),
                            (
                                "orientation",
                                object([
                                    ("swap", Json::Bool(orientation.swap)),
                                    ("flip_first", Json::Bool(orientation.flip_first)),
                                    ("flip_second", Json::Bool(orientation.flip_second)),
                                ]),
                            ),
                        ]),
                    ),
                    FaceSource::Exterior {
                        side,
                        rect,
                        features,
                    } => (
                        "exterior",
                        object([
                            ("side", side_json(side)),
                            ("rect", rect_json(*rect)),
                            (
                                "features",
                                Json::Array(features.iter().map(|feature| text(feature)).collect()),
                            ),
                        ]),
                    ),
                };
                object([
                    ("id", number(face.id)),
                    ("kind", text(kind)),
                    ("source", source),
                    ("boundary", signed_ids(&face.boundary)),
                    ("incident_cells", signed_ids(&face.incident_cells)),
                ])
            })
            .collect();
        let cells = self
            .cells
            .iter()
            .map(|cell| {
                object([
                    ("id", number(cell.id)),
                    ("leaf", leaf_json(&cell.leaf)),
                    ("block_id", text(&cell.block_id)),
                    ("boundary", signed_ids(&cell.boundary)),
                ])
            })
            .collect();
        object([
            ("schema", text("adaptive_oriented_cell_complex_v1")),
            ("grid_denominator", number(self.grid_denominator)),
            ("orientation_convention", text("block_chart_u_v_w_positive; each exterior face uses its outward chart orientation")),
            ("vertices", Json::Array(vertices)),
            ("edges", Json::Array(edges)),
            ("faces", Json::Array(faces)),
            ("cells", Json::Array(cells)),
            ("boundary_faces", Json::Array(self.boundary_faces.iter().copied().map(number).collect())),
            ("boundary_verified", Json::Bool(true)),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::{segment_provenance_connected, ChartSegment};
    use std::collections::BTreeSet;

    fn segment(root: &str) -> ChartSegment {
        ChartSegment {
            root: root.to_owned(),
            start: [0, 0, 0],
            end: [0, 2, 0],
        }
    }

    #[test]
    fn equal_endpoints_require_an_exact_d4_segment_link() {
        let left = segment("left_chart");
        let right = segment("right_chart");
        let provenance = BTreeSet::from([left.clone(), right.clone()]);
        assert!(!segment_provenance_connected(&provenance, &BTreeSet::new()));
        assert!(segment_provenance_connected(
            &provenance,
            &BTreeSet::from([(left, right)])
        ));
    }
}
