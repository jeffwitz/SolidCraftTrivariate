//! Variable-depth square pocket with explicit selection events.
//!
//! Depth uses half-millimeter integer ticks. The fixed plate is 20 x 20 x 8 mm;
//! depth saturates at breakthrough, so overdrilling does not extend the domain.
//! Homogeneous nets use one common factor of two to encode every half-millimeter
//! z coordinate exactly. This constructor does not itself qualify or adopt the
//! resulting geometry.

use crate::bracket::{BracketBlock, HOLE_EAST_PROFILE};
use crate::incidence::BlockSide;
use std::convert::TryFrom;

const SECTORS: [&str; 4] = ["east", "north", "west", "south"];
const MAX_DEPTH_TICKS: i32 = 18;
const BREAKTHROUGH_TICKS: i32 = 16;

#[derive(Debug, Clone, Copy)]
struct LayerSpan {
    z_ticks: [i32; 2],
    moving: [bool; 2],
    jets_available: bool,
}

/// Why a variable-depth pocket request was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PocketDepthError {
    /// Depth must be a finite integer tick in `[0, 18]`.
    ParameterOutOfDomain,
    /// An exact intermediate or persisted coefficient exceeds its integer type.
    NumericalRange,
}

/// Physical domain topology selected by the requested depth.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PocketRegime {
    /// No material is removed (`d = 0`).
    Solid,
    /// A blind bore has a floor strictly inside the plate (`0 < d < 8 mm`).
    BlindPocket,
    /// The bore reaches the lower plate face (`d >= 8 mm`).
    ThroughBore,
}

/// A change in the set of selected geometric faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PocketSelectionEvent {
    None,
    ZeroDepth,
    Breakthrough,
    OverdrillSaturated,
}

/// Design-face selections available in the current model state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectionState {
    pub available: Vec<&'static str>,
}

/// A stable logical design-face selection on one block side.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DesignFaceSelection {
    pub side: BlockSide,
    pub label: &'static str,
}

/// A local leaf and its control jet, when the current regime has a stable jet.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PocketDepthBlock {
    pub block: BracketBlock,
    /// Homogeneous control derivative per physical mm of depth.
    /// `None` at the zero-depth and breakthrough selection events.
    pub dcontrols_dd: Option<Vec<[i64; 4]>>,
}

/// Constructed domain signature and its explicit selection state.
#[derive(Debug, Clone, PartialEq)]
pub struct PocketDepth {
    /// Requested depth in half-millimeter ticks, including overdrill.
    pub depth_ticks: i32,
    pub effective_depth_mm: f64,
    pub floor_z_mm: Option<f64>,
    pub regime: PocketRegime,
    pub event: PocketSelectionEvent,
    /// Distinguishes the full, blind-pocket and through-bore partitions.
    pub partition_signature: &'static str,
    pub blocks: Vec<PocketDepthBlock>,
    pub design_faces: Vec<DesignFaceSelection>,
    pub selection_state: SelectionState,
    /// Ideal target volume in mm³, using the exact circular-bore formula.
    pub volume_mm3: f64,
    /// Ideal volume derivative per physical mm; absent at selection events.
    pub dvolume_dd_mm2: Option<f64>,
}

fn checked_mul(a: i128, b: i128) -> Result<i128, PocketDepthError> {
    a.checked_mul(b).ok_or(PocketDepthError::NumericalRange)
}

fn to_i64(value: i128) -> Result<i64, PocketDepthError> {
    i64::try_from(value).map_err(|_| PocketDepthError::NumericalRange)
}

fn core_disk_profile() -> Result<[[i128; 3]; 6], PocketDepthError> {
    let mut profile = [[0_i128; 3]; 6];
    for angle in 0..3 {
        let circle = HOLE_EAST_PROFILE[angle * 2 + 1];
        let square = HOLE_EAST_PROFILE[angle * 2];
        if square[0] % 10 != 0 || square[1] % 10 != 0 {
            return Err(PocketDepthError::NumericalRange);
        }
        profile[angle * 2] = circle;
        profile[angle * 2 + 1] = [square[0] / 10, square[1] / 10, square[2]];
    }
    Ok(profile)
}

fn rotate_profile(
    profile: &[[i128; 3]; 6],
    quarter_turns: usize,
) -> Result<[[i128; 3]; 6], PocketDepthError> {
    let mut result = [[0_i128; 3]; 6];
    for (index, [mut x, mut y, weight]) in profile.iter().copied().enumerate() {
        for _ in 0..quarter_turns {
            (x, y) = (y.checked_neg().ok_or(PocketDepthError::NumericalRange)?, x);
        }
        result[index] = [x, y, weight];
    }
    Ok(result)
}

fn box_leaf(
    id: String,
    z0_ticks: i32,
    z1_ticks: i32,
    moving_z0: bool,
    moving_z1: bool,
    jets_available: bool,
) -> Result<PocketDepthBlock, PocketDepthError> {
    let mut controls = Vec::with_capacity(8);
    let mut jet = Vec::with_capacity(8);
    for x in [-1_i128, 1] {
        for y in [-1_i128, 1] {
            for (z, moving) in [(z0_ticks, moving_z0), (z1_ticks, moving_z1)] {
                controls.push([
                    to_i64(checked_mul(x, 2)?)?,
                    to_i64(checked_mul(y, 2)?)?,
                    i64::from(z),
                    2,
                ]);
                jet.push([0, 0, if moving { -2 } else { 0 }, 0]);
            }
        }
    }
    Ok(PocketDepthBlock {
        block: BracketBlock {
            id,
            degrees: [1, 1, 1],
            controls,
        },
        dcontrols_dd: jets_available.then_some(jet),
    })
}

fn sector_leaf(
    id: String,
    profile: &[[i128; 3]; 6],
    span: LayerSpan,
) -> Result<PocketDepthBlock, PocketDepthError> {
    let mut controls = Vec::with_capacity(12);
    let mut jet = Vec::with_capacity(12);
    for [x, y, profile_weight] in profile {
        let weight = checked_mul(*profile_weight, 2)?;
        for (z, moving) in span.z_ticks.into_iter().zip(span.moving) {
            controls.push([
                to_i64(checked_mul(*x, 2)?)?,
                to_i64(checked_mul(*y, 2)?)?,
                to_i64(checked_mul(i128::from(z), *profile_weight)?)?,
                to_i64(weight)?,
            ]);
            jet.push([
                0,
                0,
                to_i64(if moving {
                    checked_mul(-2, *profile_weight)?
                } else {
                    0
                })?,
                0,
            ]);
        }
    }
    Ok(PocketDepthBlock {
        block: BracketBlock {
            id,
            degrees: [2, 1, 1],
            controls,
        },
        dcontrols_dd: span.jets_available.then_some(jet),
    })
}

fn add_sector_layer(
    blocks: &mut Vec<PocketDepthBlock>,
    id_prefix: &str,
    layer: &str,
    profile: &[[i128; 3]; 6],
    span: LayerSpan,
) -> Result<(), PocketDepthError> {
    for (turns, sector) in SECTORS.iter().enumerate() {
        let rotated = rotate_profile(profile, turns)?;
        blocks.push(sector_leaf(
            format!("{id_prefix}_{layer}_{sector}"),
            &rotated,
            span,
        )?);
    }
    Ok(())
}

fn outer_sector_index(id: &str) -> Option<usize> {
    ["east", "north", "west", "south"]
        .iter()
        .position(|sector| id.contains(sector))
}

fn design_faces(blocks: &[PocketDepthBlock], regime: PocketRegime) -> Vec<DesignFaceSelection> {
    let mut result = Vec::new();
    for item in blocks {
        let id = item.block.id.as_str();
        let layer = if id.ends_with("_lower") || id.contains("_lower_") {
            "lower"
        } else if id.ends_with("_upper") || id.contains("_upper_") {
            "upper"
        } else {
            "core"
        };
        let mut push = |axis, upper, label| {
            result.push(DesignFaceSelection {
                side: BlockSide {
                    block: id.to_string(),
                    axis,
                    upper,
                },
                label,
            });
        };

        if id.starts_with("pocket_outer_sector_") {
            if let Some(index) = outer_sector_index(id) {
                let boundary = ["x10", "y10", "x-10", "y-10"][index];
                // v=0 is the fixed square side; v=1 is the circular wall.
                push(1, false, boundary);
                if regime == PocketRegime::ThroughBore || layer == "upper" {
                    push(1, true, "bore_wall");
                }
            }
        }

        if layer == "upper" {
            push(2, true, "top");
            if regime == PocketRegime::BlindPocket {
                // The underside of the upper outer-sector leaf is an interface.
            }
        } else if layer == "lower" {
            if regime == PocketRegime::Solid {
                push(2, true, "top");
            } else if id == "pocket_core_lower" || id.starts_with("pocket_disk_sector_") {
                push(2, true, "pocket_floor");
            }
        }

        if layer == "lower" || (regime == PocketRegime::ThroughBore && layer == "upper") {
            push(2, false, "bottom");
        }
    }
    result.sort_by(|a, b| a.side.cmp(&b.side));
    result
}

/// Construct the square plate with a blind pocket or a through bore.
///
/// `depth_ticks` is measured in half-millimeter units and is accepted in
/// `[0,18]`. The effective removal depth is clamped at 8 mm. Logical IDs are
/// stable within each leaf role; zero-height leaves are always omitted.
pub fn build_pocket_depth(depth_ticks: i32) -> Result<PocketDepth, PocketDepthError> {
    if !(0..=MAX_DEPTH_TICKS).contains(&depth_ticks) {
        return Err(PocketDepthError::ParameterOutOfDomain);
    }

    let regime = if depth_ticks == 0 {
        PocketRegime::Solid
    } else if depth_ticks < BREAKTHROUGH_TICKS {
        PocketRegime::BlindPocket
    } else {
        PocketRegime::ThroughBore
    };
    let event = match depth_ticks {
        0 => PocketSelectionEvent::ZeroDepth,
        BREAKTHROUGH_TICKS => PocketSelectionEvent::Breakthrough,
        value if value > BREAKTHROUGH_TICKS => PocketSelectionEvent::OverdrillSaturated,
        _ => PocketSelectionEvent::None,
    };
    let floor_ticks = (16 - depth_ticks).max(0);
    let jets_available = matches!(
        event,
        PocketSelectionEvent::None | PocketSelectionEvent::OverdrillSaturated
    );

    let mut blocks = Vec::new();
    let disk_profile = core_disk_profile()?;
    if regime != PocketRegime::ThroughBore {
        let (z1_ticks, moving_z1) = if regime == PocketRegime::Solid {
            (16, false)
        } else {
            (floor_ticks, true)
        };
        if z1_ticks > 0 {
            blocks.push(box_leaf(
                "pocket_core_lower".into(),
                0,
                z1_ticks,
                false,
                moving_z1,
                jets_available,
            )?);
            add_sector_layer(
                &mut blocks,
                "pocket_disk_sector",
                "lower",
                &disk_profile,
                LayerSpan {
                    z_ticks: [0, z1_ticks],
                    moving: [false, moving_z1],
                    jets_available,
                },
            )?;
            add_sector_layer(
                &mut blocks,
                "pocket_outer_sector",
                "lower",
                &HOLE_EAST_PROFILE,
                LayerSpan {
                    z_ticks: [0, z1_ticks],
                    moving: [false, moving_z1],
                    jets_available,
                },
            )?;
        }
    }

    if regime != PocketRegime::Solid {
        let (z0_ticks, moving_z0) = if regime == PocketRegime::ThroughBore {
            (0, false)
        } else {
            (floor_ticks, true)
        };
        if 16 > z0_ticks {
            add_sector_layer(
                &mut blocks,
                "pocket_outer_sector",
                "upper",
                &HOLE_EAST_PROFILE,
                LayerSpan {
                    z_ticks: [z0_ticks, 16],
                    moving: [moving_z0, false],
                    jets_available,
                },
            )?;
        }
    }

    let partition_signature = match regime {
        PocketRegime::Solid => "square_plate_v1",
        PocketRegime::BlindPocket => "blind_pocket_v1",
        PocketRegime::ThroughBore => "through_bore_v1",
    };

    let depth_mm = f64::from(depth_ticks) * 0.5;
    let effective_depth_mm = depth_mm.min(8.0);
    let volume_mm3 = 3200.0 - 16.0 * std::f64::consts::PI * effective_depth_mm;
    if !volume_mm3.is_finite() || volume_mm3 <= 0.0 {
        return Err(PocketDepthError::NumericalRange);
    }
    let dvolume_dd_mm2 = match event {
        PocketSelectionEvent::ZeroDepth | PocketSelectionEvent::Breakthrough => None,
        PocketSelectionEvent::OverdrillSaturated => Some(0.0),
        PocketSelectionEvent::None if regime == PocketRegime::BlindPocket => {
            Some(-16.0 * std::f64::consts::PI)
        }
        PocketSelectionEvent::None => None,
    };

    let mut available: Vec<&'static str> = design_faces(&blocks, regime)
        .into_iter()
        .map(|face| face.label)
        .collect();
    available.sort_unstable();
    available.dedup();

    Ok(PocketDepth {
        depth_ticks,
        effective_depth_mm,
        floor_z_mm: if regime == PocketRegime::BlindPocket {
            Some(8.0 - effective_depth_mm)
        } else {
            None
        },
        regime,
        event,
        partition_signature,
        design_faces: design_faces(&blocks, regime),
        blocks,
        selection_state: SelectionState { available },
        volume_mm3,
        dvolume_dd_mm2,
    })
}
