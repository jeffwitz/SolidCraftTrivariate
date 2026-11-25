//! Transactional edits and fixed-position face selections for coaxial programs.

use crate::bracket::BracketBlock;
use crate::bracket_audit::BracketCertificate;
use crate::coaxial_audit::{audit_coaxial_program, CoaxialCertificate};
use crate::coaxial_program::{build_coaxial_program, CoaxialProgram, CylinderOp};

const GEOMETRY_TOLERANCE_MM: f64 = 1.0e-8;
const PROBE_STEP_MM: f64 = 1.0e-7;
const TANGENTIAL_STEP_MM: f64 = 1.0e-5;
const NORMAL_TOLERANCE: f64 = 1.0e-7;

/// A physical selection is a fixed position and outward unit normal in mm.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceSelection {
    pub point_mm: [f64; 3],
    pub normal: [f64; 3],
    pub feature: String,
}

/// Outcome when the same physical selection is checked on a later model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SelectionStatus {
    Attached { support: String },
    Disappeared,
    AmbiguousBoundary,
    ProvenanceChanged { feature: String },
}

/// One selection reconciliation produced by an edit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectionEvent {
    pub index: usize,
    pub status: SelectionStatus,
}

/// A fully reconstructed and audited state published by a session.
#[derive(Debug, Clone, PartialEq)]
pub struct CoaxialSessionState {
    pub generation: u64,
    pub program: CoaxialProgram,
    pub blocks: Vec<BracketBlock>,
    pub certificate: BracketCertificate,
    pub ideal_volume_mm3: f64,
    pub cavities: usize,
    pub selections: Vec<FaceSelection>,
    pub selection_events: Vec<SelectionEvent>,
}

/// Mutable session with atomic model replacement.
#[derive(Debug, Default)]
pub struct CoaxialSession {
    state: Option<CoaxialSessionState>,
}

#[derive(Clone, Debug)]
struct BoundaryCandidate {
    support: String,
    normal: [f64; 3],
    feature: String,
}

fn finite3(v: [f64; 3]) -> bool {
    v.into_iter().all(f64::is_finite)
}

fn norm(v: [f64; 3]) -> f64 {
    v.into_iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn normalized(v: [f64; 3]) -> Option<[f64; 3]> {
    let n = norm(v);
    (n.is_finite() && n > 0.0).then(|| v.map(|x| x / n))
}

fn add_scaled(p: [f64; 3], d: [f64; 3], s: f64) -> [f64; 3] {
    std::array::from_fn(|i| p[i] + s * d[i])
}

/// Evaluate the ordered regularized occupancy and last operation that changed it.
fn material_at(program: &CoaxialProgram, p: [f64; 3]) -> (bool, usize) {
    let r = p[0].hypot(p[1]);
    let mut material = p[0].abs() < 10.0
        && p[1].abs() < 10.0
        && p[2] > 0.0
        && p[2] < f64::from(program.base_height_mm);
    let mut changed_by = 0;
    for (index, feature) in program.features.iter().enumerate() {
        if r < f64::from(feature.radius_mm)
            && p[2] > f64::from(feature.z_mm[0])
            && p[2] < f64::from(feature.z_mm[1])
        {
            let next = feature.op == CylinderOp::Add;
            if next != material {
                material = next;
                changed_by = index + 1;
            }
        }
    }
    (material, changed_by)
}

fn owner(program: &CoaxialProgram, indices: [usize; 2]) -> String {
    let i = indices[0].max(indices[1]);
    if i == 0 {
        "base".into()
    } else {
        program.features[i - 1].id.clone()
    }
}

fn add_candidate(
    candidates: &mut Vec<BoundaryCandidate>,
    support: String,
    normal: [f64; 3],
    before: (bool, usize),
    after: (bool, usize),
    program: &CoaxialProgram,
) {
    if before.0 && !after.0 {
        candidates.push(BoundaryCandidate {
            support,
            normal,
            feature: owner(program, [before.1, after.1]),
        });
    }
}

fn physical_boundaries(program: &CoaxialProgram, p: [f64; 3]) -> Vec<BoundaryCandidate> {
    let mut candidates = Vec::new();
    let e = PROBE_STEP_MM;
    for z in program
        .features
        .iter()
        .flat_map(|f| f.z_mm)
        .chain([0, program.base_height_mm])
        .collect::<std::collections::BTreeSet<_>>()
    {
        if (p[2] - f64::from(z)).abs() > GEOMETRY_TOLERANCE_MM {
            continue;
        }
        let radial = normalized([p[0], p[1], 0.0]).unwrap_or([1.0, 0.0, 0.0]);
        for tangent in [radial, radial.map(|x| -x), [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]] {
            let sample = add_scaled(p, tangent, TANGENTIAL_STEP_MM);
            for normal in [[0., 0., 1.], [0., 0., -1.]] {
                let inside = material_at(program, add_scaled(sample, normal, -e));
                let outside = material_at(program, add_scaled(sample, normal, e));
                add_candidate(
                    &mut candidates,
                    format!("z:{z}"),
                    normal,
                    inside,
                    outside,
                    program,
                );
            }
        }
    }

    for axis in 0..2 {
        for side in [-1_i32, 1_i32] {
            if (p[axis] - 10.0 * f64::from(side)).abs() > GEOMETRY_TOLERANCE_MM {
                continue;
            }
            let mut normal = [0.0; 3];
            normal[axis] = f64::from(side);
            let name = if axis == 0 { "x" } else { "y" };
            for tangent_axis in 0..3 {
                if tangent_axis == axis {
                    continue;
                }
                for sign in [-1.0, 1.0] {
                    let mut tangent = [0.0; 3];
                    tangent[tangent_axis] = sign;
                    let sample = add_scaled(p, tangent, TANGENTIAL_STEP_MM);
                    let inside = material_at(program, add_scaled(sample, normal, -e));
                    let outside = material_at(program, add_scaled(sample, normal, e));
                    add_candidate(
                        &mut candidates,
                        format!("{name}:{side}"),
                        normal,
                        inside,
                        outside,
                        program,
                    );
                }
            }
        }
    }

    let radius = p[0].hypot(p[1]);
    if radius > 0.0 {
        let radial = [p[0] / radius, p[1] / radius, 0.0];
        let radii: std::collections::BTreeSet<_> =
            program.features.iter().map(|f| f.radius_mm).collect();
        for r in radii {
            if (radius - f64::from(r)).abs() > GEOMETRY_TOLERANCE_MM {
                continue;
            }
            for dz in [-TANGENTIAL_STEP_MM, 0.0, TANGENTIAL_STEP_MM] {
                let sample = add_scaled(p, [0.0, 0.0, 1.0], dz);
                let inward = material_at(program, add_scaled(sample, radial, -e));
                let outward = material_at(program, add_scaled(sample, radial, e));
                if inward.0 && !outward.0 {
                    add_candidate(
                        &mut candidates,
                        format!("r:{r}"),
                        radial,
                        inward,
                        outward,
                        program,
                    );
                } else if outward.0 && !inward.0 {
                    let normal = radial.map(|x| -x);
                    add_candidate(
                        &mut candidates,
                        format!("r:{r}"),
                        normal,
                        outward,
                        inward,
                        program,
                    );
                }
            }
        }
    }
    candidates.sort_by(|a, b| a.support.cmp(&b.support));
    candidates
        .dedup_by(|a, b| a.support == b.support && a.feature == b.feature && a.normal == b.normal);
    candidates
}

fn evaluate_selection(
    program: &CoaxialProgram,
    selection: &FaceSelection,
) -> Result<SelectionStatus, &'static str> {
    if !finite3(selection.point_mm) || !finite3(selection.normal) || selection.feature.is_empty() {
        return Err("coaxial_selection_invalid");
    }
    let normal = normalized(selection.normal).ok_or("coaxial_selection_invalid")?;
    if (norm(selection.normal) - 1.0).abs() > NORMAL_TOLERANCE {
        return Err("coaxial_selection_normal_not_unit");
    }
    let mut candidates = physical_boundaries(program, selection.point_mm);
    match candidates.len() {
        0 => Ok(SelectionStatus::Disappeared),
        1 => {
            let candidate = candidates.pop().expect("one candidate");
            let error = (0..3)
                .map(|i| (candidate.normal[i] - normal[i]).powi(2))
                .sum::<f64>()
                .sqrt();
            if error > NORMAL_TOLERANCE {
                return Ok(SelectionStatus::Disappeared);
            }
            if candidate.feature == selection.feature {
                Ok(SelectionStatus::Attached {
                    support: candidate.support,
                })
            } else {
                Ok(SelectionStatus::ProvenanceChanged {
                    feature: candidate.feature,
                })
            }
        }
        _ => Ok(SelectionStatus::AmbiguousBoundary),
    }
}

fn build_state(
    program: CoaxialProgram,
    generation: u64,
) -> Result<CoaxialSessionState, &'static str> {
    let blocks = build_coaxial_program(&program)?;
    let CoaxialCertificate {
        certificate,
        ideal_volume_mm3,
        cavities,
    } = audit_coaxial_program(&program, &blocks)?;
    Ok(CoaxialSessionState {
        generation,
        program,
        blocks,
        certificate,
        ideal_volume_mm3,
        cavities,
        selections: Vec::new(),
        selection_events: Vec::new(),
    })
}

impl CoaxialSession {
    /// Create a session by building and auditing its initial state.
    pub fn new(program: CoaxialProgram) -> Result<Self, &'static str> {
        Ok(Self {
            state: Some(build_state(program, 0)?),
        })
    }

    pub fn state(&self) -> Option<&CoaxialSessionState> {
        self.state.as_ref()
    }

    /// Add a selection only when it identifies one current physical face.
    pub fn select(&mut self, selection: FaceSelection) -> Result<SelectionEvent, &'static str> {
        let state = self.state.as_mut().ok_or("coaxial_session_empty")?;
        match evaluate_selection(&state.program, &selection)? {
            status @ SelectionStatus::Attached { .. } => {
                let index = state.selections.len();
                let event = SelectionEvent { index, status };
                state.selections.push(selection);
                state.selection_events = vec![event.clone()];
                Ok(event)
            }
            SelectionStatus::Disappeared => Err("coaxial_selection_disappeared"),
            SelectionStatus::AmbiguousBoundary => Err("coaxial_selection_ambiguous"),
            SelectionStatus::ProvenanceChanged { .. } => Err("coaxial_selection_provenance"),
        }
    }

    /// Build and audit before replacing the visible state; failures preserve it.
    pub fn replace(
        &mut self,
        program: CoaxialProgram,
    ) -> Result<&CoaxialSessionState, &'static str> {
        let previous = self.state.as_ref().ok_or("coaxial_session_empty")?;
        let generation = previous
            .generation
            .checked_add(1)
            .ok_or("coaxial_session_generation_overflow")?;
        let mut candidate = build_state(program, generation)?;
        candidate.selections = previous.selections.clone();
        candidate.selection_events = candidate
            .selections
            .iter()
            .enumerate()
            .map(|(index, selection)| {
                let status = evaluate_selection(&candidate.program, selection)
                    .unwrap_or(SelectionStatus::Disappeared);
                SelectionEvent { index, status }
            })
            .collect();
        self.state = Some(candidate);
        Ok(self.state.as_ref().expect("published candidate"))
    }
}
