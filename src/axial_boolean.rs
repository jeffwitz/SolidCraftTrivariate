//! Regularized booleans for products of one fixed G2 section and axial intervals.
//!
//! The operation is performed on open axial cells between sorted interval
//! bounds. Selected cells are retained as separate layers, including at a
//! coplanar union contact; their common caps are represented as internal faces
//! in the returned incidence. Cross-section intersections are outside scope.

use crate::bracket::BracketBlock;
use crate::incidence::{enumerate_faces, BlockSide, Incidence};
use crate::section_audit::audit_section;
use crate::sections::section_blocks;
use std::collections::BTreeMap;
use std::convert::TryFrom;

const AXIAL_LIMIT_MM: i32 = 1_000;

/// Boolean operation on two products of the same section and z interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxialBooleanOperation {
    Union,
    Intersection,
    Difference,
}

/// A G2 section family. Only the rectangle and pierced square are currently supported.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxialSectionFamily {
    ReferenceRectangleV1,
    CoreDiskV1,
    FullAnnulusV1,
    PiercedSquareV1,
    RoundedCornerV1,
    LayeredPocketV1,
}

impl AxialSectionFamily {
    fn name(self) -> &'static str {
        match self {
            Self::ReferenceRectangleV1 => "reference_rectangle_v1",
            Self::CoreDiskV1 => "core_disk_v1",
            Self::FullAnnulusV1 => "full_annulus_v1",
            Self::PiercedSquareV1 => "pierced_square_v1",
            Self::RoundedCornerV1 => "rounded_corner_v1",
            Self::LayeredPocketV1 => "layered_pocket_v1",
        }
    }

    fn supported(self) -> bool {
        matches!(self, Self::ReferenceRectangleV1 | Self::PiercedSquareV1)
    }
}

/// Axis along which the two product intervals are defined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxialAxis {
    X,
    Y,
    Z,
}

/// One solid defined as a fixed G2 section extruded over an integer-mm interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AxialOperand {
    pub section: AxialSectionFamily,
    pub axis: AxialAxis,
    pub interval_mm: [i32; 2],
}

/// Identity of an input operand contributing a face.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum OperandLabel {
    A,
    B,
}

/// One selected open axial cell, with its closure retained in the output model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectedAxialInterval {
    pub bounds_mm: [i32; 2],
    pub operands: Vec<OperandLabel>,
}

/// Original face attribute and orientation for one contributor to a result face.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaceProvenance {
    pub operand: OperandLabel,
    /// G2 source tag, such as `x10`, `circle_r4`, `z0`, or `z8`.
    pub source_face: String,
    /// True only for an oriented cap inherited from the subtracted B tool.
    pub orientation_reversed: bool,
}

/// A labelled exterior face of the result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AxialBooleanFace {
    pub side: BlockSide,
    pub provenance: Vec<FaceProvenance>,
}

/// Complete layered model for one regularized axial boolean.
#[derive(Debug, Clone, PartialEq)]
pub struct AxialBooleanModel {
    pub operation: AxialBooleanOperation,
    pub section: AxialSectionFamily,
    pub axis: AxialAxis,
    pub selected_intervals: Vec<SelectedAxialInterval>,
    pub component_count: usize,
    pub blocks: Vec<BracketBlock>,
    pub incidence: Incidence,
    pub faces: Vec<AxialBooleanFace>,
    /// Ideal target measure in mm³, based on the section's analytic G2 measure.
    pub volume_mm3: f64,
}

/// Specific refusal classes for the currently qualified axial-product domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxialBooleanError {
    UnsupportedSection,
    DifferentSections,
    IncompatibleAxes,
    UnsupportedAxis,
    InvertedInterval,
    NullInterval,
    IntervalOutOfRange,
    UnsupportedSectionModel,
    NumericalRange,
    Incidence(&'static str),
    ProvenanceMismatch,
}

#[derive(Clone, Copy)]
struct AxialCell {
    bounds: [i32; 2],
    in_a: bool,
    in_b: bool,
}

#[derive(Clone)]
struct OutputBlockMeta {
    source_id: String,
    interval: [i32; 2],
    in_a: bool,
    in_b: bool,
}

fn validate_interval(interval: [i32; 2]) -> Result<(), AxialBooleanError> {
    if interval[0] < -AXIAL_LIMIT_MM || interval[1] > AXIAL_LIMIT_MM {
        return Err(AxialBooleanError::IntervalOutOfRange);
    }
    if interval[0] > interval[1] {
        return Err(AxialBooleanError::InvertedInterval);
    }
    if interval[0] == interval[1] {
        return Err(AxialBooleanError::NullInterval);
    }
    Ok(())
}

fn selected(operation: AxialBooleanOperation, in_a: bool, in_b: bool) -> bool {
    match operation {
        AxialBooleanOperation::Union => in_a || in_b,
        AxialBooleanOperation::Intersection => in_a && in_b,
        AxialBooleanOperation::Difference => in_a && !in_b,
    }
}

fn make_cells(
    operation: AxialBooleanOperation,
    interval_a: [i32; 2],
    interval_b: [i32; 2],
) -> Vec<AxialCell> {
    let mut bounds = vec![interval_a[0], interval_a[1], interval_b[0], interval_b[1]];
    bounds.sort_unstable();
    bounds.dedup();
    bounds
        .windows(2)
        .map(|w| {
            let [lo, hi] = [w[0], w[1]];
            let in_a = lo >= interval_a[0] && hi <= interval_a[1];
            let in_b = lo >= interval_b[0] && hi <= interval_b[1];
            AxialCell {
                bounds: [lo, hi],
                in_a: in_a && selected(operation, in_a, in_b),
                in_b: in_b && selected(operation, in_a, in_b),
            }
        })
        .collect()
}

fn reslice_block(
    source: &BracketBlock,
    id: String,
    interval: [i32; 2],
) -> Result<BracketBlock, AxialBooleanError> {
    let mut controls = Vec::with_capacity(source.controls.len());
    for [x, y, z, weight] in &source.controls {
        let w = i128::from(*weight);
        let z = i128::from(*z);
        let output_z = if z == 0 {
            i128::from(interval[0]) * w
        } else if z == 8 * w {
            i128::from(interval[1]) * w
        } else {
            return Err(AxialBooleanError::UnsupportedSectionModel);
        };
        controls.push([
            *x,
            *y,
            i64::try_from(output_z).map_err(|_| AxialBooleanError::NumericalRange)?,
            *weight,
        ]);
    }
    Ok(BracketBlock {
        id,
        degrees: source.degrees,
        controls,
    })
}

fn source_faces(
    source_blocks: &[BracketBlock],
    section: AxialSectionFamily,
) -> Result<BTreeMap<(String, usize, bool), String>, AxialBooleanError> {
    let certificate = audit_section(section.name(), source_blocks)
        .map_err(|_| AxialBooleanError::UnsupportedSectionModel)?;
    Ok(certificate
        .design_faces
        .into_iter()
        .map(|(side, label)| ((side.block, side.axis, side.upper), label))
        .collect())
}

fn cell_operands(cell: AxialCell) -> Vec<OperandLabel> {
    let mut result = Vec::with_capacity(2);
    if cell.in_a {
        result.push(OperandLabel::A);
    }
    if cell.in_b {
        result.push(OperandLabel::B);
    }
    result
}

fn cap_sources(
    operation: AxialBooleanOperation,
    output_upper: bool,
    z: i32,
    interval_a: [i32; 2],
    interval_b: [i32; 2],
) -> Vec<FaceProvenance> {
    let output_normal = if output_upper { 1 } else { -1 };
    let mut result = Vec::new();
    for (operand, interval) in [(OperandLabel::A, interval_a), (OperandLabel::B, interval_b)] {
        let reversed = operation == AxialBooleanOperation::Difference && operand == OperandLabel::B;
        let endpoints = [(interval[0], -1, "z0"), (interval[1], 1, "z8")];
        for (coordinate, mut source_normal, source_face) in endpoints {
            if coordinate != z {
                continue;
            }
            if reversed {
                source_normal = -source_normal;
            }
            if source_normal == output_normal {
                result.push(FaceProvenance {
                    operand,
                    source_face: source_face.to_string(),
                    orientation_reversed: reversed,
                });
            }
        }
    }
    result
}

/// Build a regularized boolean of two solids in the same G2 cross-section.
///
/// Both operands must use the z axis and the same supported section family.
/// Axial bounds are ordered integer millimeters in `[-1000,1000]`. Selected
/// subintervals are retained as distinct layers rather than fused, so every
/// internal cap match remains visible in the exact face incidence.
pub fn build_axial_boolean(
    operation: AxialBooleanOperation,
    a: AxialOperand,
    b: AxialOperand,
) -> Result<AxialBooleanModel, AxialBooleanError> {
    if a.axis != b.axis {
        return Err(AxialBooleanError::IncompatibleAxes);
    }
    if a.axis != AxialAxis::Z {
        return Err(AxialBooleanError::UnsupportedAxis);
    }
    if a.section != b.section {
        return Err(AxialBooleanError::DifferentSections);
    }
    if !a.section.supported() {
        return Err(AxialBooleanError::UnsupportedSection);
    }
    validate_interval(a.interval_mm)?;
    validate_interval(b.interval_mm)?;

    // Audit the source recipe before transforming any cells.
    let source_blocks =
        section_blocks(a.section.name()).map_err(|_| AxialBooleanError::UnsupportedSectionModel)?;
    let source_labels = source_faces(&source_blocks, a.section)?;
    let cells = make_cells(operation, a.interval_mm, b.interval_mm);
    let mut selected_intervals = Vec::new();
    let mut blocks = Vec::new();
    let mut metadata = BTreeMap::<String, OutputBlockMeta>::new();
    let mut component_count = 0usize;
    let mut was_selected = false;
    let mut selected_length = 0_i128;

    for (index, cell) in cells.into_iter().enumerate() {
        let is_selected = cell.in_a || cell.in_b;
        if !is_selected {
            was_selected = false;
            continue;
        }
        if !was_selected {
            component_count += 1;
        }
        was_selected = true;
        selected_length = selected_length
            .checked_add(i128::from(cell.bounds[1] - cell.bounds[0]))
            .ok_or(AxialBooleanError::NumericalRange)?;
        selected_intervals.push(SelectedAxialInterval {
            bounds_mm: cell.bounds,
            operands: cell_operands(cell),
        });

        for source in &source_blocks {
            let id = format!("axial_boolean_cell_{index:02}_{}", source.id);
            let output = reslice_block(source, id.clone(), cell.bounds)?;
            metadata.insert(
                id,
                OutputBlockMeta {
                    source_id: source.id.clone(),
                    interval: cell.bounds,
                    in_a: cell.in_a,
                    in_b: cell.in_b,
                },
            );
            blocks.push(output);
        }
    }

    let incidence = enumerate_faces(&blocks).map_err(AxialBooleanError::Incidence)?;
    let mut faces = Vec::with_capacity(incidence.exterior.len());
    for side in &incidence.exterior {
        let meta = metadata
            .get(&side.block)
            .ok_or(AxialBooleanError::ProvenanceMismatch)?;
        let mut provenance = if side.axis == 2 {
            let z = if side.upper {
                meta.interval[1]
            } else {
                meta.interval[0]
            };
            cap_sources(operation, side.upper, z, a.interval_mm, b.interval_mm)
        } else {
            let label = source_labels
                .get(&(meta.source_id.clone(), side.axis, side.upper))
                .ok_or(AxialBooleanError::ProvenanceMismatch)?;
            let mut labels = Vec::with_capacity(2);
            if meta.in_a {
                labels.push(FaceProvenance {
                    operand: OperandLabel::A,
                    source_face: label.clone(),
                    orientation_reversed: false,
                });
            }
            if meta.in_b {
                labels.push(FaceProvenance {
                    operand: OperandLabel::B,
                    source_face: label.clone(),
                    orientation_reversed: false,
                });
            }
            labels
        };
        provenance.sort_by_key(|entry| (entry.operand, entry.source_face.clone()));
        if provenance.is_empty() {
            return Err(AxialBooleanError::ProvenanceMismatch);
        }
        faces.push(AxialBooleanFace {
            side: side.clone(),
            provenance,
        });
    }
    faces.sort_by(|left, right| left.side.cmp(&right.side));

    let area = match a.section {
        AxialSectionFamily::ReferenceRectangleV1 => 400.0,
        AxialSectionFamily::PiercedSquareV1 => 400.0 - 16.0 * std::f64::consts::PI,
        _ => return Err(AxialBooleanError::UnsupportedSection),
    };
    let volume_mm3 = area * selected_length as f64;
    if !volume_mm3.is_finite() || volume_mm3 < 0.0 {
        return Err(AxialBooleanError::NumericalRange);
    }

    Ok(AxialBooleanModel {
        operation,
        section: a.section,
        axis: a.axis,
        selected_intervals,
        component_count,
        blocks,
        incidence,
        faces,
        volume_mm3,
    })
}
