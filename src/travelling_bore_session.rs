//! Transactional translation of the bore through the pocket and plate edge.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    incidence::BlockSide,
    travelling_bore::travelling_bore_blocks,
    travelling_bore_audit::audit_travelling_bore,
    travelling_bore_bounded::travelling_bore_bounded_blocks,
    travelling_bore_bounded_audit::audit_travelling_bore_bounded,
    travelling_bore_events::{analyze_travelling_bore, TravellingBoreReport},
};
use std::collections::{BTreeMap, BTreeSet};

/// Construction contract selected for every transaction in a session.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TravelGeometry {
    #[default]
    ExactMicrometreCuts,
    BoundedRationalArcs,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TravelFidelity {
    ExactCircularArcs,
    CertifiedApproximation {
        radial_upper_mm: f64,
        hausdorff_upper_mm: f64,
    },
}

#[derive(Debug)]
pub struct TravelState {
    pub center_x_um: i32,
    pub physical: TravellingBoreReport,
    pub fidelity: TravelFidelity,
    pub blocks: Vec<BracketBlock>,
    pub certificate: BracketCertificate,
    pub selections: BTreeMap<String, Vec<BlockSide>>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct TravelChange {
    /// Critical center positions along the directed straight translation.
    /// The departure point is excluded and the arrival point included.
    pub crossed_centers_um: Vec<i32>,
    pub created_features: Vec<String>,
    pub lost_features: Vec<String>,
    /// Surviving design identities whose supporting block sides have changed.
    pub remapped_features: Vec<String>,
}

#[derive(Default)]
pub struct TravellingBoreSession {
    geometry: TravelGeometry,
    current: Option<TravelState>,
}

impl TravellingBoreSession {
    pub fn with_geometry(geometry: TravelGeometry) -> Self {
        Self {
            geometry,
            current: None,
        }
    }

    pub fn current(&self) -> Option<&TravelState> {
        self.current.as_ref()
    }

    /// A missing required selection or an invalid new geometry aborts the
    /// transaction. Feature names refer to physical design identities.
    pub fn apply(
        &mut self,
        center_x_um: i32,
        required_features: &[&str],
    ) -> Result<TravelChange, &'static str> {
        let physical = analyze_travelling_bore(i64::from(center_x_um))?;
        let (blocks, certificate, fidelity) = match self.geometry {
            TravelGeometry::ExactMicrometreCuts => {
                let blocks = travelling_bore_blocks(center_x_um)?;
                let cert = audit_travelling_bore(center_x_um, &blocks)?;
                (blocks, cert, TravelFidelity::ExactCircularArcs)
            }
            TravelGeometry::BoundedRationalArcs => {
                let blocks = travelling_bore_bounded_blocks(center_x_um)?;
                let result = audit_travelling_bore_bounded(center_x_um, &blocks)?;
                let fidelity = TravelFidelity::CertifiedApproximation {
                    radial_upper_mm: result.radial_error_upper_mm,
                    hausdorff_upper_mm: result.hausdorff_error_upper_mm,
                };
                (blocks, result.certificate, fidelity)
            }
        };
        let mut selections: BTreeMap<String, Vec<BlockSide>> = BTreeMap::new();
        for (side, feature) in &certificate.design_faces {
            selections
                .entry(feature.clone())
                .or_default()
                .push(side.clone());
        }
        if required_features
            .iter()
            .any(|f| !selections.contains_key(*f))
        {
            return Err("travelling_bore_required_selection_absent");
        }
        let mut crossed_centers_um = Vec::new();
        let new_names: BTreeSet<_> = selections.keys().cloned().collect();
        let mut old_names = BTreeSet::new();
        let mut remapped_features = Vec::new();
        if let Some(old) = &self.current {
            old_names.extend(old.selections.keys().cloned());
            for feature in new_names.intersection(&old_names) {
                if old.selections[feature] != selections[feature] {
                    remapped_features.push(feature.clone());
                }
            }
            for candidate in [-4000, -2000, 0, 8000, 10000, 12000] {
                if (old.center_x_um < candidate && candidate <= center_x_um)
                    || (center_x_um <= candidate && candidate < old.center_x_um)
                {
                    crossed_centers_um.push(candidate);
                }
            }
            if center_x_um < old.center_x_um {
                crossed_centers_um.reverse();
            }
        }
        let change = TravelChange {
            crossed_centers_um,
            created_features: new_names.difference(&old_names).cloned().collect(),
            lost_features: old_names.difference(&new_names).cloned().collect(),
            remapped_features,
        };
        self.current = Some(TravelState {
            center_x_um,
            physical,
            fidelity,
            blocks,
            certificate,
            selections,
        });
        Ok(change)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tangency_and_unrepresentable_intersection_preserve_the_accepted_state() {
        let mut session = TravellingBoreSession::default();
        session.apply(3000, &["bore_wall"]).unwrap();
        let change = session.apply(-800, &["bore_wall"]).unwrap();
        assert_eq!(change.crossed_centers_um, [0]);
        let old = session.current().unwrap();
        let blocks = old.blocks.clone();
        let certificate = old.certificate.clone();
        let selections = old.selections.clone();
        for refused in [0, -100] {
            assert!(session.apply(refused, &["bore_wall"]).is_err());
            let state = session.current().unwrap();
            assert_eq!(state.center_x_um, -800);
            assert_eq!(state.blocks, blocks);
            assert_eq!(state.certificate, certificate);
            assert_eq!(state.selections, selections);
        }
    }

    #[test]
    fn bounded_edits_preserve_fidelity_and_selections_across_rejected_events() {
        let mut session = TravellingBoreSession::with_geometry(TravelGeometry::BoundedRationalArcs);
        session.apply(-100, &["bore_wall"]).unwrap();
        assert!(session.apply(0, &[]).is_err());
        assert_eq!(session.current().unwrap().center_x_um, -100);
        session.apply(11900, &["bore_wall"]).unwrap();
        let saved = session.current().unwrap().fidelity;
        let TravelFidelity::CertifiedApproximation {
            hausdorff_upper_mm, ..
        } = saved
        else {
            panic!("missing bounded fidelity")
        };
        assert!((0.0..=1e-8).contains(&hausdorff_upper_mm));
        assert!(session.apply(13000, &["bore_wall"]).is_err());
        assert_eq!(session.current().unwrap().fidelity, saved);
        let mut cold = TravellingBoreSession::with_geometry(TravelGeometry::BoundedRationalArcs);
        cold.apply(11900, &[]).unwrap();
        let (a, b) = (session.current().unwrap(), cold.current().unwrap());
        assert_eq!(a.blocks, b.blocks);
        assert_eq!(a.selections, b.selections);
        assert_eq!(a.certificate, b.certificate);
        assert_eq!(a.fidelity, b.fidelity);
    }

    #[test]
    fn crossing_the_edge_reconciles_physical_selections_with_cold_replay() {
        let mut session = TravellingBoreSession::default();
        session.apply(11920, &["bore_wall"]).unwrap();
        assert_eq!(
            session.apply(13000, &["bore_wall"]),
            Err("travelling_bore_required_selection_absent")
        );
        assert_eq!(session.current().unwrap().center_x_um, 11920);
        let outside = session.apply(13000, &[]).unwrap();
        assert_eq!(outside.crossed_centers_um, [12000]);
        assert!(outside.lost_features.contains(&"bore_wall".to_string()));
        let reverse = session.apply(-800, &["bore_wall"]).unwrap();
        assert_eq!(reverse.crossed_centers_um, [12000, 10000, 8000, 0]);
        assert!(reverse.created_features.contains(&"bore_wall".to_string()));
        let mut cold = TravellingBoreSession::default();
        cold.apply(-800, &["bore_wall"]).unwrap();
        let (a, b) = (session.current().unwrap(), cold.current().unwrap());
        assert_eq!(a.blocks, b.blocks);
        assert_eq!(a.selections, b.selections);
        assert_eq!(a.certificate, b.certificate);
    }
}
