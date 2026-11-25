//! Transactional session for audited, variable-depth pocket models.
//!
//! Candidate geometry is audited before its analysis space and stiffness are
//! built. The current published state and both dependency caches are replaced
//! only after every candidate stage succeeds.

use crate::analysis_space::{from_certified_faces, AnalysisSpace};
use crate::bracket::BracketBlock;
use crate::bracket_audit::BracketCertificate;
use crate::elasticity::{assemble, gauss, Stiffness};
use crate::incidence::Incidence;
use crate::pocket_depth::{build_pocket_depth, PocketDepth, PocketDepthError};
use crate::pocket_depth_audit::audit_depth;
use std::sync::Arc;

/// Requested depth and analysis material/integration settings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PocketRequest {
    pub depth_ticks: i32,
    pub young_mpa: f64,
    pub poisson: f64,
    pub quadrature_order: usize,
}

/// Explicit selection changes between the prior and candidate states.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SelectionDelta {
    pub created: Vec<&'static str>,
    pub disappeared: Vec<&'static str>,
    /// Selections removed without mapping to a replacement selection.
    pub unassigned: Vec<&'static str>,
}

/// Which cached dependencies were actually reused by this update.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CacheReuse {
    pub topology: bool,
    pub certified_geometry: bool,
    pub analysis_space: bool,
    pub stiffness: bool,
}

/// Fully committed, audited session state.
pub struct PocketSessionState {
    pub request: PocketRequest,
    pub domain: PocketDepth,
    pub certificate: Arc<BracketCertificate>,
    pub analysis_space: Arc<AnalysisSpace>,
    pub stiffness: Arc<Stiffness>,
    pub selection_delta: SelectionDelta,
    pub cache_reuse: CacheReuse,
}

/// Failure in candidate construction or analysis. No state is published on error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PocketSessionError {
    InvalidDepth,
    InvalidMaterial,
    InvalidQuadratureOrder,
    Audit(&'static str),
    SelectionMismatch,
    AnalysisSpace(&'static str),
    Assembly(&'static str),
}

struct GeometryCache {
    blocks: Vec<BracketBlock>,
    certificate: Arc<BracketCertificate>,
    analysis_space: Arc<AnalysisSpace>,
}

struct StiffnessCache {
    blocks: Vec<BracketBlock>,
    young_bits: u64,
    poisson_bits: u64,
    order: usize,
    stiffness: Arc<Stiffness>,
}

struct TopologyCache {
    partition_signature: &'static str,
    incidence: Arc<Incidence>,
}

/// Mutable transactional session. Use `update` to change depth or analysis settings.
#[derive(Default)]
pub struct PocketSession {
    state: Option<PocketSessionState>,
    topology_cache: Option<TopologyCache>,
    geometry_cache: Option<GeometryCache>,
    stiffness_cache: Option<StiffnessCache>,
}

fn validate_request(request: PocketRequest) -> Result<(), PocketSessionError> {
    if !request.young_mpa.is_finite()
        || request.young_mpa <= 0.0
        || !request.poisson.is_finite()
        || !(-1.0..0.5).contains(&request.poisson)
    {
        return Err(PocketSessionError::InvalidMaterial);
    }
    gauss(request.quadrature_order).map_err(|_| PocketSessionError::InvalidQuadratureOrder)?;
    Ok(())
}

fn selection_delta(previous: Option<&PocketDepth>, candidate: &PocketDepth) -> SelectionDelta {
    let Some(previous) = previous else {
        return SelectionDelta {
            created: candidate.selection_state.available.clone(),
            disappeared: Vec::new(),
            unassigned: Vec::new(),
        };
    };
    let mut created: Vec<_> = candidate
        .selection_state
        .available
        .iter()
        .copied()
        .filter(|label| !previous.selection_state.available.contains(label))
        .collect();
    let mut disappeared: Vec<_> = previous
        .selection_state
        .available
        .iter()
        .copied()
        .filter(|label| !candidate.selection_state.available.contains(label))
        .collect();
    created.sort_unstable();
    created.dedup();
    disappeared.sort_unstable();
    disappeared.dedup();
    // A label absent from the new state's available selection set has no
    // replacement unless its exact stable label is present again.
    let unassigned = disappeared.clone();
    SelectionDelta {
        created,
        disappeared,
        unassigned,
    }
}

impl PocketSession {
    /// Create and commit an initial audited state.
    pub fn new(request: PocketRequest) -> Result<Self, PocketSessionError> {
        let mut session = Self::default();
        session.update(request)?;
        Ok(session)
    }

    /// Return the last committed state, if the session has been initialized.
    pub fn state(&self) -> Option<&PocketSessionState> {
        self.state.as_ref()
    }

    /// Build, audit, and assemble a candidate transactionally.
    pub fn update(
        &mut self,
        request: PocketRequest,
    ) -> Result<&PocketSessionState, PocketSessionError> {
        validate_request(request)?;
        let domain = build_pocket_depth(request.depth_ticks).map_err(|error| match error {
            PocketDepthError::ParameterOutOfDomain => PocketSessionError::InvalidDepth,
            PocketDepthError::NumericalRange => PocketSessionError::InvalidDepth,
        })?;
        let blocks: Vec<BracketBlock> = domain
            .blocks
            .iter()
            .map(|leaf| leaf.block.clone())
            .collect();

        let geometry_reused = self
            .geometry_cache
            .as_ref()
            .is_some_and(|cache| cache.blocks == blocks);
        let (certificate, incidence, analysis_space, topology_reused) = if geometry_reused {
            let cache = self
                .geometry_cache
                .as_ref()
                .expect("checked geometry cache");
            let topology_reused = self.topology_cache.as_ref().is_some_and(|topology| {
                topology.partition_signature == domain.partition_signature
                    && *topology.incidence == cache.certificate.incidence
            });
            let incidence = if topology_reused {
                Arc::clone(
                    &self
                        .topology_cache
                        .as_ref()
                        .expect("checked topology cache")
                        .incidence,
                )
            } else {
                Arc::new(cache.certificate.incidence.clone())
            };
            (
                Arc::clone(&cache.certificate),
                incidence,
                Arc::clone(&cache.analysis_space),
                topology_reused,
            )
        } else {
            // Audit is the publication gate: no analysis object exists before it.
            let certificate = Arc::new(
                audit_depth(request.depth_ticks, &blocks).map_err(PocketSessionError::Audit)?,
            );
            let mut constructed_faces: Vec<_> = domain
                .design_faces
                .iter()
                .map(|face| (face.side.clone(), face.label.to_string()))
                .collect();
            constructed_faces.sort_by(|left, right| left.0.cmp(&right.0));
            if constructed_faces != certificate.design_faces {
                return Err(PocketSessionError::SelectionMismatch);
            }
            let topology_reused = self.topology_cache.as_ref().is_some_and(|topology| {
                topology.partition_signature == domain.partition_signature
                    && *topology.incidence == certificate.incidence
            });
            let incidence = if topology_reused {
                Arc::clone(
                    &self
                        .topology_cache
                        .as_ref()
                        .expect("checked topology cache")
                        .incidence,
                )
            } else {
                Arc::new(certificate.incidence.clone())
            };
            let analysis_space = Arc::new(
                from_certified_faces(&blocks, &incidence, &certificate.design_faces)
                    .map_err(PocketSessionError::AnalysisSpace)?,
            );
            (certificate, incidence, analysis_space, topology_reused)
        };

        let young_bits = request.young_mpa.to_bits();
        let poisson_bits = request.poisson.to_bits();
        let stiffness_reused = self.stiffness_cache.as_ref().is_some_and(|cache| {
            cache.blocks == blocks
                && cache.young_bits == young_bits
                && cache.poisson_bits == poisson_bits
                && cache.order == request.quadrature_order
        });
        let stiffness = if stiffness_reused {
            Arc::clone(
                &self
                    .stiffness_cache
                    .as_ref()
                    .expect("checked stiffness cache")
                    .stiffness,
            )
        } else {
            Arc::new(
                assemble(
                    &analysis_space,
                    request.young_mpa,
                    request.poisson,
                    request.quadrature_order,
                )
                .map_err(PocketSessionError::Assembly)?,
            )
        };

        let delta = selection_delta(self.state.as_ref().map(|state| &state.domain), &domain);
        let partition_signature = domain.partition_signature;
        let state = PocketSessionState {
            request,
            domain,
            certificate: Arc::clone(&certificate),
            analysis_space: Arc::clone(&analysis_space),
            stiffness: Arc::clone(&stiffness),
            selection_delta: delta,
            cache_reuse: CacheReuse {
                topology: topology_reused,
                certified_geometry: geometry_reused,
                analysis_space: geometry_reused,
                stiffness: stiffness_reused,
            },
        };

        // Commit all dependencies and the public state only after every stage succeeded.
        self.geometry_cache = Some(GeometryCache {
            blocks: blocks.clone(),
            certificate,
            analysis_space,
        });
        self.topology_cache = Some(TopologyCache {
            partition_signature,
            incidence,
        });
        self.stiffness_cache = Some(StiffnessCache {
            blocks,
            young_bits,
            poisson_bits,
            order: request.quadrature_order,
            stiffness,
        });
        self.state = Some(state);
        Ok(self.state.as_ref().expect("committed state"))
    }
}

#[cfg(test)]
mod tests {
    use super::{PocketRequest, PocketSession, PocketSessionError};
    use std::sync::Arc;

    fn request(depth_ticks: i32) -> PocketRequest {
        PocketRequest {
            depth_ticks,
            young_mpa: 210_000.0,
            poisson: 0.3,
            quadrature_order: 1,
        }
    }

    fn same_stiffness(left: &crate::elasticity::Stiffness, right: &crate::elasticity::Stiffness) {
        assert_eq!(left.offsets, right.offsets);
        assert_eq!(left.columns, right.columns);
        assert_eq!(
            left.integrated_volume.to_bits(),
            right.integrated_volume.to_bits()
        );
        assert_eq!(left.values, right.values);
    }

    #[test]
    fn depth_sequence_reuses_only_invariant_dependencies() {
        let mut session = PocketSession::new(request(0)).unwrap();
        assert!(!session
            .state()
            .unwrap()
            .domain
            .selection_state
            .available
            .contains(&"pocket_floor"));
        for (depth, topology_reused) in [(2, false), (8, true), (15, true), (16, false)] {
            let state = session.update(request(depth)).unwrap();
            assert_eq!(state.cache_reuse.topology, topology_reused);
            assert!(!state.cache_reuse.certified_geometry);
            assert!(!state.cache_reuse.analysis_space);
            assert!(!state.cache_reuse.stiffness);
            if depth == 8 || depth == 15 {
                assert!(state.selection_delta.created.is_empty());
                assert!(state.selection_delta.disappeared.is_empty());
            }
        }

        let state = session.update(request(17)).unwrap();
        assert!(state.cache_reuse.topology);
        assert!(state.cache_reuse.certified_geometry);
        assert!(state.cache_reuse.analysis_space);
        assert!(state.cache_reuse.stiffness);

        let before_overdrill = session.state().unwrap();
        let certificate = Arc::clone(&before_overdrill.certificate);
        let space = Arc::clone(&before_overdrill.analysis_space);
        let stiffness = Arc::clone(&before_overdrill.stiffness);
        let state = session.update(request(18)).unwrap();
        assert!(state.cache_reuse.topology);
        assert!(state.cache_reuse.certified_geometry);
        assert!(state.cache_reuse.analysis_space);
        assert!(state.cache_reuse.stiffness);
        assert!(Arc::ptr_eq(&certificate, &state.certificate));
        assert!(Arc::ptr_eq(&space, &state.analysis_space));
        assert!(Arc::ptr_eq(&stiffness, &state.stiffness));
        assert!(state.selection_delta.created.is_empty());
        assert!(state.selection_delta.disappeared.is_empty());

        let state = session.update(request(8)).unwrap();
        assert!(!state.cache_reuse.certified_geometry);
        assert_eq!(state.selection_delta.created, vec!["pocket_floor"]);
        assert!(state.selection_delta.disappeared.is_empty());
        let cold = PocketSession::new(request(8)).unwrap();
        let cold_state = cold.state().unwrap();
        assert_eq!(state.domain.blocks, cold_state.domain.blocks);
        assert_eq!(*state.analysis_space, *cold_state.analysis_space);
        same_stiffness(&state.stiffness, &cold_state.stiffness);
    }

    #[test]
    fn material_and_order_updates_reuse_geometry_but_assemble_again() {
        let mut session = PocketSession::new(request(8)).unwrap();
        let first = session.state().unwrap();
        let certificate = Arc::clone(&first.certificate);
        let space = Arc::clone(&first.analysis_space);
        let stiffness = Arc::clone(&first.stiffness);

        let mut material = request(8);
        material.young_mpa = 70_000.0;
        let state = session.update(material).unwrap();
        assert!(state.cache_reuse.certified_geometry);
        assert!(state.cache_reuse.topology);
        assert!(state.cache_reuse.topology);
        assert!(state.cache_reuse.analysis_space);
        assert!(!state.cache_reuse.stiffness);
        assert!(Arc::ptr_eq(&certificate, &state.certificate));
        assert!(Arc::ptr_eq(&space, &state.analysis_space));
        assert!(!Arc::ptr_eq(&stiffness, &state.stiffness));

        let current_stiffness = Arc::clone(&state.stiffness);
        material.quadrature_order = 2;
        let state = session.update(material).unwrap();
        assert!(state.cache_reuse.certified_geometry);
        assert!(state.cache_reuse.analysis_space);
        assert!(!state.cache_reuse.stiffness);
        assert!(Arc::ptr_eq(&certificate, &state.certificate));
        assert!(Arc::ptr_eq(&space, &state.analysis_space));
        assert!(!Arc::ptr_eq(&current_stiffness, &state.stiffness));

        let mut combined = request(9);
        combined.young_mpa = 205_000.0;
        combined.poisson = 0.28;
        combined.quadrature_order = 3;
        let state = session.update(combined).unwrap();
        assert!(state.cache_reuse.topology);
        assert!(!state.cache_reuse.certified_geometry);
        assert!(!state.cache_reuse.analysis_space);
        assert!(!state.cache_reuse.stiffness);
        let cold = PocketSession::new(combined).unwrap();
        same_stiffness(
            &session.state().unwrap().stiffness,
            &cold.state().unwrap().stiffness,
        );
    }

    #[test]
    fn selection_delta_and_rejected_updates_are_transactional() {
        let mut session = PocketSession::new(request(15)).unwrap();
        let state = session.update(request(16)).unwrap();
        assert_eq!(state.selection_delta.disappeared, vec!["pocket_floor"]);
        assert_eq!(state.selection_delta.unassigned, vec!["pocket_floor"]);
        assert!(!state.selection_delta.created.contains(&"bottom"));

        let previous_request = session.state().unwrap().request;
        let previous_certificate = Arc::clone(&session.state().unwrap().certificate);
        let previous_space = Arc::clone(&session.state().unwrap().analysis_space);
        let previous_stiffness = Arc::clone(&session.state().unwrap().stiffness);
        assert_eq!(
            session.update(request(-1)).err(),
            Some(PocketSessionError::InvalidDepth)
        );
        let mut invalid_material = request(8);
        invalid_material.young_mpa = f64::NAN;
        assert_eq!(
            session.update(invalid_material).err(),
            Some(PocketSessionError::InvalidMaterial)
        );
        let state = session.state().unwrap();
        assert_eq!(state.request, previous_request);
        assert!(Arc::ptr_eq(&previous_certificate, &state.certificate));
        assert!(Arc::ptr_eq(&previous_space, &state.analysis_space));
        assert!(Arc::ptr_eq(&previous_stiffness, &state.stiffness));
    }
}
