//! Exact dyadic octree planning for conforming multipatch complexes.
//!
//! This module plans topology only. It does not infer adjacency from physical
//! coordinates and does not construct or certify refined geometry.

use crate::incidence::{BlockSide, Incidence, Orientation};
use std::collections::BTreeSet;

const MAX_DEPTH: u8 = 20;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LeafKey {
    pub root: String,
    pub level: u8,
    pub index: [u32; 3],
}

impl LeafKey {
    pub fn root(root: impl Into<String>) -> Self {
        Self {
            root: root.into(),
            level: 0,
            index: [0; 3],
        }
    }
    fn children(&self) -> [Self; 8] {
        std::array::from_fn(|child| Self {
            root: self.root.clone(),
            level: self.level + 1,
            index: std::array::from_fn(|axis| 2 * self.index[axis] + ((child >> axis) & 1) as u32),
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct LeafSide {
    pub leaf: LeafKey,
    /// Semantic side of the unsplit root. This carries boundary-feature
    /// provenance through all subsequent dyadic refinements.
    pub root_side: BlockSide,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct FaceRect {
    /// Half-open intervals in an integer grid whose denominator is
    /// `AdaptivePartition::grid_denominator()`. Axis order is increasing root
    /// coordinate after omitting the face normal.
    pub ranges: [[u32; 2]; 2],
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FacePatch {
    pub first: LeafSide,
    pub second: LeafSide,
    /// The same overlap expressed in the first/root-master face coordinates.
    pub master_rect: FaceRect,
    /// Overlap rectangles in each root's own face parameter coordinates.
    pub first_subrect: FaceRect,
    pub second_subrect: FaceRect,
    /// D4 orientation from the second trace into the first trace.
    pub orientation: Orientation,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExteriorPatch {
    pub side: LeafSide,
    /// Rectangle in the semantic root side's parameter coordinates.
    pub root_rect: FaceRect,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RefinementReport {
    pub requested: Vec<LeafKey>,
    pub induced: Vec<LeafKey>,
    pub leaves: Vec<LeafKey>,
    pub interfaces: Vec<FacePatch>,
    pub exterior: Vec<ExteriorPatch>,
}

#[derive(Debug, Clone)]
pub struct AdaptivePartition {
    roots: Vec<String>,
    incidence: Incidence,
    leaves: BTreeSet<LeafKey>,
    max_depth: u8,
    max_cells: usize,
}

fn side_key(side: &BlockSide) -> (&str, usize, bool) {
    (&side.block, side.axis, side.upper)
}

fn tangent_axes(axis: usize) -> [usize; 2] {
    match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}

fn valid_side(side: &BlockSide, roots: &BTreeSet<String>) -> bool {
    side.axis < 3 && roots.contains(&side.block)
}

fn integer_scale(max_depth: u8) -> u32 {
    1_u32 << max_depth
}

fn leaf_rect(leaf: &LeafKey, axis: usize, max_depth: u8) -> FaceRect {
    let scale = 1_u32 << (max_depth - leaf.level);
    let [a, b] = tangent_axes(axis);
    FaceRect {
        ranges: [
            [leaf.index[a] * scale, (leaf.index[a] + 1) * scale],
            [leaf.index[b] * scale, (leaf.index[b] + 1) * scale],
        ],
    }
}

fn intersect_rect(a: FaceRect, b: FaceRect) -> Option<FaceRect> {
    let ranges = std::array::from_fn(|axis| {
        [
            a.ranges[axis][0].max(b.ranges[axis][0]),
            a.ranges[axis][1].min(b.ranges[axis][1]),
        ]
    });
    (ranges.iter().all(|r| r[0] < r[1])).then_some(FaceRect { ranges })
}

fn transform_point_to_second(uv: [u32; 2], side: u32, orientation: Orientation) -> [u32; 2] {
    let a = if orientation.flip_first {
        side - uv[0]
    } else {
        uv[0]
    };
    let b = if orientation.flip_second {
        side - uv[1]
    } else {
        uv[1]
    };
    if orientation.swap { [b, a] } else { [a, b] }
}

fn transform_rect_second_to_first(rect: FaceRect, side: u32, o: Orientation) -> FaceRect {
    let mapped = |axis: usize| -> [u32; 2] {
        let input_axis = if o.swap { 1 - axis } else { axis };
        let flip = if o.swap {
            if axis == 0 {
                o.flip_first
            } else {
                o.flip_second
            }
        } else if axis == 0 {
            o.flip_first
        } else {
            o.flip_second
        };
        let [lo, hi] = rect.ranges[input_axis];
        if flip {
            [side - hi, side - lo]
        } else {
            [lo, hi]
        }
    };
    FaceRect {
        ranges: [mapped(0), mapped(1)],
    }
}

fn transform_rect_first_to_second(rect: FaceRect, side: u32, o: Orientation) -> FaceRect {
    let mapped = |axis: usize| -> [u32; 2] {
        let input_axis = if o.swap { 1 - axis } else { axis };
        let flip = if o.swap {
            if axis == 0 {
                o.flip_second
            } else {
                o.flip_first
            }
        } else if axis == 0 {
            o.flip_first
        } else {
            o.flip_second
        };
        let [lo, hi] = rect.ranges[input_axis];
        if flip {
            [side - hi, side - lo]
        } else {
            [lo, hi]
        }
    };
    FaceRect {
        ranges: [mapped(0), mapped(1)],
    }
}

fn locate_leaf(
    leaves: &BTreeSet<LeafKey>,
    root: &str,
    point: [u32; 3],
    point_depth: u8,
    max_depth: u8,
) -> Option<LeafKey> {
    for level in (0..=max_depth.min(point_depth)).rev() {
        let shift = point_depth - level;
        let index = point.map(|coord| coord >> shift);
        let candidate = LeafKey {
            root: root.to_owned(),
            level,
            index,
        };
        if leaves.contains(&candidate) {
            return Some(candidate);
        }
    }
    None
}

fn query_quadrants(leaf: &LeafKey, axis: usize, upper: bool) -> Vec<([u32; 3], u8)> {
    let level = leaf.level;
    let point_depth = level + 2;
    let mut result = Vec::with_capacity(4);
    let tangents = tangent_axes(axis);
    for bit0 in 0..2_u32 {
        for bit1 in 0..2_u32 {
            let mut point = [0_u32; 3];
            let normal_boundary = leaf.index[axis] + u32::from(upper);
            let normal_scaled = normal_boundary << (point_depth - level);
            point[axis] = if upper {
                normal_scaled + 1
            } else {
                normal_scaled - 1
            };
            for (which, bit) in [bit0, bit1].into_iter().enumerate() {
                let t = tangents[which];
                point[t] = (leaf.index[t] << 2) + 2 * bit + 1;
            }
            result.push((point, point_depth));
        }
    }
    result
}

fn external_query_points(leaf: &LeafKey, axis: usize) -> Vec<([u32; 2], u8)> {
    let level = leaf.level;
    let point_depth = level + 2;
    let tangents = tangent_axes(axis);
    let mut points = Vec::with_capacity(4);
    for bit0 in 0..2_u32 {
        for bit1 in 0..2_u32 {
            points.push((
                [
                    (leaf.index[tangents[0]] << 2) + 2 * bit0 + 1,
                    (leaf.index[tangents[1]] << 2) + 2 * bit1 + 1,
                ],
                point_depth,
            ));
        }
    }
    points
}

fn root_face_neighbor(
    leaves: &BTreeSet<LeafKey>,
    side: &BlockSide,
    point: [u32; 2],
    point_depth: u8,
    max_depth: u8,
) -> Option<LeafKey> {
    let mut local = [0_u32; 3];
    let [a, b] = tangent_axes(side.axis);
    local[a] = point[0];
    local[b] = point[1];
    let denom_depth = point_depth;
    let denom = 1_u32 << denom_depth;
    local[side.axis] = if side.upper { denom - 1 } else { 1 };
    locate_leaf(leaves, &side.block, local, denom_depth, max_depth)
}

fn face_patch(
    first: &LeafKey,
    first_side: &BlockSide,
    second: &LeafKey,
    second_side: &BlockSide,
    orientation: Orientation,
    max_depth: u8,
) -> Result<FacePatch, &'static str> {
    let first_rect = leaf_rect(first, first_side.axis, max_depth);
    let second_rect = leaf_rect(second, second_side.axis, max_depth);
    let scale = integer_scale(max_depth);
    let second_in_first = transform_rect_second_to_first(second_rect, scale, orientation);
    let master_rect = intersect_rect(first_rect, second_in_first)
        .ok_or("adjacent leaf faces have no positive-area overlap")?;
    let first_subrect = master_rect;
    let second_subrect = transform_rect_first_to_second(master_rect, scale, orientation);
    Ok(FacePatch {
        first: LeafSide {
            leaf: first.clone(),
            root_side: first_side.clone(),
        },
        second: LeafSide {
            leaf: second.clone(),
            root_side: second_side.clone(),
        },
        master_rect,
        first_subrect,
        second_subrect,
        orientation,
    })
}

impl AdaptivePartition {
    /// Validate a complete root incidence and initialize one leaf per root.
    pub fn new(
        roots: &[String],
        incidence: &Incidence,
        max_depth: u8,
        max_cells: usize,
    ) -> Result<Self, &'static str> {
        if roots.is_empty() || max_depth > MAX_DEPTH || max_cells < roots.len() {
            return Err("invalid adaptive partition limits or empty roots");
        }
        let root_set = roots.iter().cloned().collect::<BTreeSet<_>>();
        if root_set.len() != roots.len() || roots.iter().any(|r| r.is_empty()) {
            return Err("root identities must be unique and nonempty");
        }
        let mut sides = BTreeSet::new();
        for side in &incidence.exterior {
            if !valid_side(side, &root_set) || !sides.insert(side.clone()) {
                return Err("invalid or duplicate exterior root side");
            }
        }
        for interface in &incidence.interfaces {
            if !valid_side(&interface.first, &root_set)
                || !valid_side(&interface.second, &root_set)
                || interface.first == interface.second
                || !sides.insert(interface.first.clone())
                || !sides.insert(interface.second.clone())
            {
                return Err("invalid or multiply assigned root interface side");
            }
        }
        if sides.len() != 6 * root_set.len()
            || root_set.iter().any(|root| {
                (0..3).any(|axis| {
                    [false, true].into_iter().any(|upper| {
                        !sides.contains(&BlockSide {
                            block: root.clone(),
                            axis,
                            upper,
                        })
                    })
                })
            })
        {
            return Err("every root side must occur exactly once in the initial incidence");
        }
        let mut normalized = incidence.clone();
        normalized.interfaces.sort_by(|a, b| {
            side_key(&a.first)
                .cmp(&side_key(&b.first))
                .then_with(|| side_key(&a.second).cmp(&side_key(&b.second)))
        });
        normalized.exterior.sort();
        let leaves = root_set.iter().cloned().map(LeafKey::root).collect();
        Ok(Self {
            roots: root_set.into_iter().collect(),
            incidence: normalized,
            leaves,
            max_depth,
            max_cells,
        })
    }

    pub fn leaves(&self) -> &BTreeSet<LeafKey> {
        &self.leaves
    }
    pub fn roots(&self) -> &[String] {
        &self.roots
    }
    pub fn max_depth(&self) -> u8 {
        self.max_depth
    }
    pub fn max_cells(&self) -> usize {
        self.max_cells
    }
    pub fn grid_denominator(&self) -> u32 {
        integer_scale(self.max_depth)
    }

    /// Return the initial or current interface plan, including whole root
    /// interfaces before any refinement is requested.
    pub fn interfaces(&self) -> Result<Vec<FacePatch>, &'static str> {
        let mut patches = Vec::new();
        // Internal same-root faces are found from the leaf on the lower side.
        for leaf in &self.leaves {
            for axis in 0..3 {
                if leaf.index[axis] + 1 >= (1_u32 << leaf.level) {
                    continue;
                }
                for (point, depth) in query_quadrants(leaf, axis, true) {
                    if let Some(other) =
                        locate_leaf(&self.leaves, &leaf.root, point, depth, self.max_depth)
                    {
                        if other != *leaf {
                            let side_first = BlockSide {
                                block: leaf.root.clone(),
                                axis,
                                upper: true,
                            };
                            let side_second = BlockSide {
                                block: leaf.root.clone(),
                                axis,
                                upper: false,
                            };
                            let patch = face_patch(
                                leaf,
                                &side_first,
                                &other,
                                &side_second,
                                Orientation {
                                    swap: false,
                                    flip_first: false,
                                    flip_second: false,
                                },
                                self.max_depth,
                            )?;
                            patches.push(patch);
                        }
                    }
                }
            }
        }

        // Incidence between distinct roots is always expressed in the first
        // root's face frame, then mapped through the declared D4 orientation.
        for interface in &self.incidence.interfaces {
            for leaf in self
                .leaves
                .iter()
                .filter(|l| l.root == interface.first.block)
            {
                let boundary_index = if interface.first.upper {
                    (1_u32 << leaf.level) - 1
                } else {
                    0
                };
                if leaf.index[interface.first.axis] != boundary_index {
                    continue;
                }
                for (point, depth) in external_query_points(leaf, interface.first.axis) {
                    let scale = 1_u32 << depth;
                    let mapped = transform_point_to_second(point, scale, interface.orientation);
                    if let Some(other) = root_face_neighbor(
                        &self.leaves,
                        &interface.second,
                        mapped,
                        depth,
                        self.max_depth,
                    ) {
                        let patch = face_patch(
                            leaf,
                            &interface.first,
                            &other,
                            &interface.second,
                            interface.orientation,
                            self.max_depth,
                        )?;
                        patches.push(patch);
                    }
                }
            }
        }
        patches.sort_by(|a, b| {
            a.first
                .cmp(&b.first)
                .then_with(|| a.second.cmp(&b.second))
                .then_with(|| a.master_rect.cmp(&b.master_rect))
        });
        patches.dedup_by(|a, b| {
            a.first == b.first && a.second == b.second && a.master_rect == b.master_rect
        });
        Ok(patches)
    }

    pub fn exterior_patches(&self) -> Vec<ExteriorPatch> {
        let mut patches = Vec::new();
        for root_side in &self.incidence.exterior {
            for leaf in self.leaves.iter().filter(|l| l.root == root_side.block) {
                let boundary_index = if root_side.upper {
                    (1_u32 << leaf.level) - 1
                } else {
                    0
                };
                if leaf.index[root_side.axis] == boundary_index {
                    patches.push(ExteriorPatch {
                        side: LeafSide {
                            leaf: leaf.clone(),
                            root_side: root_side.clone(),
                        },
                        root_rect: leaf_rect(leaf, root_side.axis, self.max_depth),
                    });
                }
            }
        }
        patches.sort_by(|a, b| a.side.cmp(&b.side));
        patches
    }

    pub fn snapshot(&self) -> Result<RefinementReport, &'static str> {
        Ok(RefinementReport {
            requested: Vec::new(),
            induced: Vec::new(),
            leaves: self.leaves.iter().cloned().collect(),
            interfaces: self.interfaces()?,
            exterior: self.exterior_patches(),
        })
    }

    /// Refine requested active leaves and propagate face 2:1 balance.
    /// Mutation is committed only after every requested split, propagated
    /// split, budget check, and interface reconstruction succeeds.
    pub fn refine(&mut self, requests: &[LeafKey]) -> Result<RefinementReport, &'static str> {
        if requests.is_empty() {
            return Err("refinement request must be nonempty");
        }
        let requested = requests.iter().cloned().collect::<BTreeSet<_>>();
        if requested.len() != requests.len()
            || requested.iter().any(|leaf| !self.leaves.contains(leaf))
        {
            return Err("requested leaf is duplicated or inactive");
        }
        let mut leaves = self.leaves.clone();
        let mut induced = BTreeSet::new();
        let mut pending = requested.iter().cloned().collect::<Vec<_>>();
        while let Some(leaf) = pending.pop() {
            if !leaves.contains(&leaf) {
                continue;
            }
            if leaf.level >= self.max_depth {
                return Err("maximum refinement depth exceeded");
            }
            leaves.remove(&leaf);
            leaves.extend(leaf.children());
            if leaves.len() > self.max_cells {
                return Err("maximum cell budget exceeded");
            }
        }

        loop {
            let old = AdaptivePartition {
                roots: self.roots.clone(),
                incidence: self.incidence.clone(),
                leaves: leaves.clone(),
                max_depth: self.max_depth,
                max_cells: self.max_cells,
            };
            let interfaces = old.interfaces()?;
            let mut need_split = BTreeSet::new();
            for patch in interfaces {
                let la = patch.first.leaf.level;
                let lb = patch.second.leaf.level;
                if la.abs_diff(lb) > 1 {
                    need_split.insert(if la < lb {
                        patch.first.leaf
                    } else {
                        patch.second.leaf
                    });
                }
            }
            if need_split.is_empty() {
                break;
            }
            for leaf in need_split {
                if leaf.level >= self.max_depth {
                    return Err("maximum refinement depth exceeded during balance propagation");
                }
                if !leaves.remove(&leaf) {
                    continue;
                }
                if !requested.contains(&leaf) {
                    induced.insert(leaf.clone());
                }
                leaves.extend(leaf.children());
                if leaves.len() > self.max_cells {
                    return Err("maximum cell budget exceeded during balance propagation");
                }
            }
        }

        let candidate = Self {
            roots: self.roots.clone(),
            incidence: self.incidence.clone(),
            leaves,
            max_depth: self.max_depth,
            max_cells: self.max_cells,
        };
        let report = RefinementReport {
            requested: requested.into_iter().collect(),
            induced: induced.into_iter().collect(),
            leaves: candidate.leaves.iter().cloned().collect(),
            interfaces: candidate.interfaces()?,
            exterior: candidate.exterior_patches(),
        };
        self.leaves = candidate.leaves;
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn side(root: &str, axis: usize, upper: bool) -> BlockSide {
        BlockSide {
            block: root.into(),
            axis,
            upper,
        }
    }
    fn one_interface(orientation: Orientation) -> Incidence {
        Incidence {
            interfaces: vec![crate::incidence::Interface {
                first: side("a", 0, true),
                second: side("b", 0, false),
                orientation,
            }],
            exterior: [
                side("a", 0, false),
                side("a", 1, false),
                side("a", 1, true),
                side("a", 2, false),
                side("a", 2, true),
                side("b", 0, true),
                side("b", 1, false),
                side("b", 1, true),
                side("b", 2, false),
                side("b", 2, true),
            ]
            .into(),
        }
    }
    fn partition(depth: u8, cells: usize, orientation: Orientation) -> AdaptivePartition {
        AdaptivePartition::new(
            &["a".into(), "b".into()],
            &one_interface(orientation),
            depth,
            cells,
        )
        .unwrap()
    }
    fn id(o: Orientation) -> Orientation {
        o
    }

    #[test]
    fn initial_state_exposes_whole_interfaces_and_semantic_boundary() {
        let p = partition(
            4,
            100,
            id(Orientation {
                swap: false,
                flip_first: false,
                flip_second: false,
            }),
        );
        let initial = p.snapshot().unwrap();
        assert_eq!(initial.leaves.len(), 2);
        assert_eq!(initial.interfaces.len(), 1);
        assert_eq!(initial.exterior.len(), 10);
        assert_eq!(initial.interfaces[0].master_rect.ranges, [[0, 16], [0, 16]]);
        assert_eq!(p.grid_denominator(), 16);
    }

    #[test]
    fn one_side_refinement_produces_four_d4_face_patches() {
        let orientations = [false, true]
            .into_iter()
            .flat_map(|swap| {
                [false, true].into_iter().flat_map(move |flip_first| {
                    [false, true]
                        .into_iter()
                        .map(move |flip_second| Orientation {
                            swap,
                            flip_first,
                            flip_second,
                        })
                })
            })
            .collect::<Vec<_>>();
        for orientation in orientations {
            let mut p = partition(4, 100, orientation);
            let report = p.refine(&[LeafKey::root("a")]).unwrap();
            assert_eq!(report.requested, vec![LeafKey::root("a")]);
            assert!(report.induced.is_empty());
            assert_eq!(report.leaves.len(), 9);
            assert_eq!(report.interfaces.len(), 16);
            let cross = report
                .interfaces
                .iter()
                .filter(|face| face.first.leaf.root != face.second.leaf.root)
                .collect::<Vec<_>>();
            assert_eq!(cross.len(), 4);
            assert_eq!(report.interfaces.len() - cross.len(), 12);
            assert_eq!(report.exterior.len(), 25);
            assert!(cross.iter().all(|face| {
                face.orientation == orientation
                    && face.second_subrect
                        == transform_rect_first_to_second(
                            face.master_rect,
                            p.grid_denominator(),
                            orientation,
                        )
            }));
            let area = cross
                .iter()
                .map(|face| {
                    u64::from(face.master_rect.ranges[0][1] - face.master_rect.ranges[0][0])
                        * u64::from(face.master_rect.ranges[1][1] - face.master_rect.ranges[1][0])
                })
                .sum::<u64>();
            assert_eq!(area, u64::from(p.grid_denominator()).pow(2));
            for (i, first) in cross.iter().enumerate() {
                for second in &cross[i + 1..] {
                    assert!(intersect_rect(first.master_rect, second.master_rect).is_none());
                }
            }
        }
    }

    #[test]
    fn cross_root_two_to_one_balance_is_induced_and_atomic() {
        let mut p = partition(
            4,
            200,
            Orientation {
                swap: false,
                flip_first: false,
                flip_second: false,
            },
        );
        p.refine(&[LeafKey::root("a")]).unwrap();
        let face_child = LeafKey {
            root: "a".into(),
            level: 1,
            index: [1, 0, 0],
        };
        let report = p.refine(std::slice::from_ref(&face_child)).unwrap();
        assert!(report.induced.contains(&LeafKey::root("b")));
        assert_eq!(report.leaves.len(), 23);
        assert!(
            report
                .interfaces
                .iter()
                .all(|f| f.first.leaf.level.abs_diff(f.second.leaf.level) <= 1)
        );

        let mut bounded = partition(
            4,
            9,
            Orientation {
                swap: false,
                flip_first: false,
                flip_second: false,
            },
        );
        let before = bounded.snapshot().unwrap();
        assert!(bounded.refine(&[LeafKey::root("a")]).is_ok());
        let before_failure = bounded.snapshot().unwrap();
        assert!(
            bounded
                .refine(&[LeafKey {
                    root: "a".into(),
                    level: 1,
                    index: [0, 0, 0]
                }])
                .is_err()
        );
        assert_eq!(bounded.snapshot().unwrap(), before_failure);
        assert_eq!(before.leaves.len(), 2);
    }

    #[test]
    fn permutation_and_invalid_incidence_are_deterministic() {
        let orient = Orientation {
            swap: true,
            flip_first: false,
            flip_second: true,
        };
        let incidence = one_interface(orient);
        let forward = partition(3, 100, orient).snapshot().unwrap();
        let mut permuted = incidence.clone();
        permuted.exterior.reverse();
        permuted.interfaces.reverse();
        let other = AdaptivePartition::new(&["b".into(), "a".into()], &permuted, 3, 100)
            .unwrap()
            .snapshot()
            .unwrap();
        assert_eq!(forward, other);

        let mut first = partition(3, 100, orient);
        let mut second = partition(3, 100, orient);
        let a = LeafKey::root("a");
        let b = LeafKey::root("b");
        let report_a = first.refine(&[a.clone(), b.clone()]).unwrap();
        let report_b = second.refine(&[b, a]).unwrap();
        assert_eq!(report_a, report_b);

        let mut invalid = incidence;
        invalid.exterior.push(side("a", 0, true));
        assert!(AdaptivePartition::new(&["a".into(), "b".into()], &invalid, 3, 100).is_err());
    }

    #[test]
    fn e6_roots_retain_the_certified_whole_face_incidence() {
        let blocks = crate::bracket::bracket_blocks().unwrap();
        let incidence = crate::incidence::enumerate_faces(&blocks).unwrap();
        let roots = blocks.iter().map(|b| b.id.clone()).collect::<Vec<_>>();
        let partition = AdaptivePartition::new(&roots, &incidence, 8, 10_000).unwrap();
        let snapshot = partition.snapshot().unwrap();
        assert_eq!(snapshot.leaves.len(), 36);
        assert_eq!(snapshot.interfaces.len(), 55);
        assert_eq!(snapshot.exterior.len(), 106);
    }
}
