//! Native geometry on the closed unit cube. Physical coordinates are in mm.

pub mod analysis_space;
pub mod annular;
pub mod annular_audit;
pub mod bezier;
pub mod bracket;
pub mod bracket_audit;
mod bracket_json;
pub mod incidence;
pub mod elasticity;
pub mod interchange;
pub mod json;
pub mod planar;
pub use bezier::{Evaluation, MapError, RationalBezierMap};
