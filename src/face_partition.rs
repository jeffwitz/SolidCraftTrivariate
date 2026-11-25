//! Exact dyadic coverage of a validated master face by interface or boundary patches.

use crate::face_constraints::{nested_constraints, FaceNet, Row, Span};
use crate::incidence::{BlockSide, Orientation};
use std::collections::BTreeSet;

/// Stable identity of one whole topological subface.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct TopologicalFaceId(String);

impl TopologicalFaceId {
    /// Create a nonempty topological face identity.
    pub fn new(value: impl Into<String>) -> Result<Self, &'static str> {
        let value = value.into();
        if value.is_empty() {
            return Err("topological face identity must be nonempty");
        }
        Ok(Self(value))
    }

    /// Borrow the stable identity string.
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// A master block side whose exact homogeneous net has been validated.
#[derive(Debug, Clone)]
pub struct MasterSide {
    side: BlockSide,
    net: FaceNet,
}

impl MasterSide {
    /// Validate and retain a master side and its homogeneous Bézier face net.
    pub fn validate(side: BlockSide, net: FaceNet) -> Result<Self, &'static str> {
        validate_side(&side)?;
        validate_net(&net)?;
        Ok(Self { side, net })
    }

    /// Return the block side designated as master.
    pub fn side(&self) -> &BlockSide {
        &self.side
    }

    /// Return the validated homogeneous master face net.
    pub fn net(&self) -> &FaceNet {
        &self.net
    }
}

/// The use assigned to a dyadic rectangle of the master side.
#[derive(Debug, Clone)]
pub enum FaceUse {
    /// A design boundary with the named feature.
    Boundary { feature: String },
    /// An interface to one uniquely identified slave block side.
    Interface {
        slave: BlockSide,
        net: FaceNet,
        orientation: Orientation,
    },
}

/// One half-or-whole dyadic rectangle on the master square.
#[derive(Debug, Clone)]
pub struct FacePatch {
    pub id: TopologicalFaceId,
    pub spans: [Span; 2],
    pub usage: FaceUse,
}

/// The classified contents of one returned topological subface.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PartitionUse {
    Boundary {
        feature: String,
    },
    Interface {
        slave: BlockSide,
        orientation: Orientation,
        constraints: Vec<Row>,
    },
}

/// A validated face assignment and its exact field transfer, when applicable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PartitionCell {
    pub id: TopologicalFaceId,
    pub spans: [Span; 2],
    pub usage: PartitionUse,
}

/// Exact local partition of one master side.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FacePartition {
    pub master: BlockSide,
    /// Cells are ordered by their lower-left dyadic coordinate, then face ID.
    pub cells: Vec<PartitionCell>,
    /// Atomic quadrant IDs in `(u-lower/v-lower, u-lower/v-upper,
    /// u-upper/v-lower, u-upper/v-upper)` order.
    pub atomic_cells: [TopologicalFaceId; 4],
    /// Partial boundary patches require basis support adaptation before BCs.
    /// No whole-face DOFs are blocked by this result.
    pub requires_boundary_basis_adaptation: bool,
}

fn validate_side(side: &BlockSide) -> Result<(), &'static str> {
    if side.block.is_empty() || side.axis > 2 {
        return Err("invalid block side identity");
    }
    Ok(())
}

fn validate_net(net: &FaceNet) -> Result<(), &'static str> {
    if net.degrees.iter().any(|degree| !(1..=2).contains(degree))
        || net.controls.len() != (net.degrees[0] + 1) * (net.degrees[1] + 1)
        || net.controls.iter().any(|control| control[3] <= 0)
    {
        return Err("invalid homogeneous face net");
    }
    Ok(())
}

fn span_cells(span: Span) -> &'static [usize] {
    match span {
        Span::Full => &[0, 1],
        Span::Lower => &[0],
        Span::Upper => &[1],
    }
}

fn span_range(span: Span) -> (usize, usize) {
    match span {
        Span::Full => (0, 2),
        Span::Lower => (0, 1),
        Span::Upper => (1, 2),
    }
}

fn patch_sort_key(patch: &PartitionCell) -> ((usize, usize, usize, usize), String) {
    let (u0, u1) = span_range(patch.spans[0]);
    let (v0, v1) = span_range(patch.spans[1]);
    ((u0, v0, u1, v1), patch.id.0.clone())
}

/// Validate exact dyadic coverage and build T1 constraints on every interface.
///
/// This proves only coverage of the supplied master parameter square. It does
/// not prove that the physical neighboring blocks form a globally conforming
/// or nonoverlapping volume mesh. A partial boundary is identified explicitly
/// so a caller can adapt the field basis before applying boundary conditions.
pub fn partition_master(
    master: &MasterSide,
    patches: &[FacePatch],
) -> Result<FacePartition, &'static str> {
    if patches.is_empty() {
        return Err("master face has no coverage patches");
    }
    let mut topological_ids = BTreeSet::new();
    let mut slave_sides = BTreeSet::new();
    let mut atom_owner: [Option<usize>; 4] = [None, None, None, None];
    let mut cells = Vec::with_capacity(patches.len());
    let mut partial_boundary = false;

    for (patch_index, patch) in patches.iter().enumerate() {
        if !topological_ids.insert(patch.id.as_str().to_owned()) {
            return Err("duplicate topological face identity");
        }
        let full = patch.spans == [Span::Full, Span::Full];
        let usage = match &patch.usage {
            FaceUse::Boundary { feature } => {
                if feature.is_empty() {
                    return Err("boundary feature identity must be nonempty");
                }
                if !full {
                    partial_boundary = true;
                }
                PartitionUse::Boundary {
                    feature: feature.clone(),
                }
            }
            FaceUse::Interface {
                slave,
                net,
                orientation,
            } => {
                validate_side(slave)?;
                if slave == &master.side {
                    return Err("master side cannot be its own slave");
                }
                if !slave_sides.insert(slave.clone()) {
                    return Err("duplicate slave side identity");
                }
                let constraints = nested_constraints(&master.net, net, patch.spans, *orientation)?;
                PartitionUse::Interface {
                    slave: slave.clone(),
                    orientation: *orientation,
                    constraints,
                }
            }
        };
        for &u in span_cells(patch.spans[0]) {
            for &v in span_cells(patch.spans[1]) {
                let atom = 2 * u + v;
                if atom_owner[atom].replace(patch_index).is_some() {
                    return Err("overlapping master face coverage");
                }
            }
        }
        cells.push(PartitionCell {
            id: patch.id.clone(),
            spans: patch.spans,
            usage,
        });
    }

    if atom_owner.iter().any(Option::is_none) {
        return Err("master face coverage has a gap");
    }
    let atomic_cells = std::array::from_fn(|atom| {
        cells[atom_owner[atom].expect("coverage was checked")]
            .id
            .clone()
    });
    cells.sort_by_key(patch_sort_key);
    Ok(FacePartition {
        master: master.side.clone(),
        cells,
        atomic_cells,
        requires_boundary_basis_adaptation: partial_boundary,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn master_side() -> MasterSide {
        MasterSide::validate(
            BlockSide {
                block: "master".into(),
                axis: 0,
                upper: true,
            },
            FaceNet {
                degrees: [1, 1],
                controls: vec![[0, 0, 0, 1], [0, 1, 0, 1], [1, 0, 0, 1], [1, 1, 0, 1]],
            },
        )
        .unwrap()
    }

    fn span_start(span: Span) -> i64 {
        match span {
            Span::Full | Span::Lower => 0,
            Span::Upper => 1,
        }
    }

    fn slave_net(spans: [Span; 2]) -> FaceNet {
        let u0 = span_start(spans[0]);
        let v0 = span_start(spans[1]);
        let u1 = u0 + if spans[0] == Span::Full { 2 } else { 1 };
        let v1 = v0 + if spans[1] == Span::Full { 2 } else { 1 };
        FaceNet {
            degrees: [1, 1],
            controls: vec![
                [u0, v0, 0, 2],
                [u0, v1, 0, 2],
                [u1, v0, 0, 2],
                [u1, v1, 0, 2],
            ],
        }
    }

    fn interface_patch(id: &str, spans: [Span; 2], slave_id: &str) -> FacePatch {
        FacePatch {
            id: TopologicalFaceId::new(id).unwrap(),
            spans,
            usage: FaceUse::Interface {
                slave: BlockSide {
                    block: slave_id.into(),
                    axis: 1,
                    upper: false,
                },
                net: slave_net(spans),
                orientation: Orientation {
                    swap: false,
                    flip_first: false,
                    flip_second: false,
                },
            },
        }
    }

    fn four_interfaces() -> Vec<FacePatch> {
        let mut patches = Vec::new();
        for u in [Span::Lower, Span::Upper] {
            for v in [Span::Lower, Span::Upper] {
                let suffix = format!(
                    "{}{}",
                    if u == Span::Lower { 0 } else { 1 },
                    if v == Span::Lower { 0 } else { 1 }
                );
                patches.push(interface_patch(
                    &format!("tf{suffix}"),
                    [u, v],
                    &format!("slave{suffix}"),
                ));
            }
        }
        patches
    }

    #[test]
    fn four_atomic_interfaces_cover_without_gaps_and_are_order_independent() {
        let master = master_side();
        assert_eq!(master.side().block, "master");
        assert_eq!(master.net().degrees, [1, 1]);
        let patches = four_interfaces();
        let forward = partition_master(&master, &patches).unwrap();
        let reversed: Vec<_> = patches.into_iter().rev().collect();
        let backward = partition_master(&master, &reversed).unwrap();
        assert_eq!(forward.atomic_cells, backward.atomic_cells);
        assert_eq!(forward.cells, backward.cells);
        assert_eq!(forward.cells.len(), 4);
        assert!(!forward.requires_boundary_basis_adaptation);
        assert!(forward.cells.iter().all(|cell| matches!(&cell.usage,
            PartitionUse::Interface { constraints, .. } if constraints.len() == 4)));
    }

    #[test]
    fn gap_overlap_duplicate_face_and_slave_are_rejected() {
        let master = master_side();
        let mut gap = four_interfaces();
        gap.pop();
        assert_eq!(
            partition_master(&master, &gap),
            Err("master face coverage has a gap")
        );

        let mut overlap = four_interfaces();
        overlap.push(interface_patch(
            "extra",
            [Span::Lower, Span::Lower],
            "slave_extra",
        ));
        assert_eq!(
            partition_master(&master, &overlap),
            Err("overlapping master face coverage")
        );

        let mut duplicate_face = four_interfaces();
        duplicate_face[1].id = duplicate_face[0].id.clone();
        assert_eq!(
            partition_master(&master, &duplicate_face),
            Err("duplicate topological face identity")
        );

        let mut duplicate_slave = four_interfaces();
        let first_slave = match &duplicate_slave[0].usage {
            FaceUse::Interface { slave, .. } => slave.clone(),
            FaceUse::Boundary { .. } => unreachable!(),
        };
        if let FaceUse::Interface { slave, .. } = &mut duplicate_slave[1].usage {
            *slave = first_slave;
        }
        assert_eq!(
            partition_master(&master, &duplicate_slave),
            Err("duplicate slave side identity")
        );
    }

    #[test]
    fn partial_boundary_is_reported_for_basis_adaptation_only() {
        let master = master_side();
        let patches = vec![
            FacePatch {
                id: TopologicalFaceId::new("interface_lower").unwrap(),
                spans: [Span::Lower, Span::Full],
                usage: FaceUse::Interface {
                    slave: BlockSide {
                        block: "slave_lower".into(),
                        axis: 1,
                        upper: false,
                    },
                    net: slave_net([Span::Lower, Span::Full]),
                    orientation: Orientation {
                        swap: false,
                        flip_first: false,
                        flip_second: false,
                    },
                },
            },
            FacePatch {
                id: TopologicalFaceId::new("wall_upper").unwrap(),
                spans: [Span::Upper, Span::Full],
                usage: FaceUse::Boundary {
                    feature: "design_wall".into(),
                },
            },
        ];
        let partition = partition_master(&master, &patches).unwrap();
        assert!(partition.requires_boundary_basis_adaptation);
        assert_eq!(
            partition.atomic_cells,
            [
                TopologicalFaceId::new("interface_lower").unwrap(),
                TopologicalFaceId::new("interface_lower").unwrap(),
                TopologicalFaceId::new("wall_upper").unwrap(),
                TopologicalFaceId::new("wall_upper").unwrap(),
            ]
        );
        assert!(partition.cells.iter().all(|cell| match &cell.usage {
            PartitionUse::Boundary { .. } => true,
            PartitionUse::Interface { constraints, .. } => !constraints.is_empty(),
        }));
    }

    #[test]
    fn altered_slave_trace_is_rejected_by_t1() {
        let master = master_side();
        let mut patch = interface_patch("tf", [Span::Lower, Span::Lower], "slave");
        if let FaceUse::Interface { net, .. } = &mut patch.usage {
            net.controls[3][0] += 1;
        }
        assert_eq!(
            partition_master(&master, &[patch]),
            Err("incompatible restricted homogeneous trace")
        );
    }
}
