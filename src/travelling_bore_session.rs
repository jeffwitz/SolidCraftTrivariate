//! Transactional translation of the bore through the pocket and plate edge.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    incidence::BlockSide,
    travelling_bore::travelling_bore_blocks,
    travelling_bore_audit::audit_travelling_bore,
    travelling_bore_events::{analyze_travelling_bore, TravellingBoreReport},
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug)]
pub struct TravelState {
    pub center_x_um: i32,
    pub physical: TravellingBoreReport,
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
    current: Option<TravelState>,
}

impl TravellingBoreSession {
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
        let blocks = travelling_bore_blocks(center_x_um)?;
        let certificate = audit_travelling_bore(center_x_um, &blocks)?;
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
