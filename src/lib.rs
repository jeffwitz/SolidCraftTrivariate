//! Native geometry on the closed unit cube. Physical coordinates are in mm.

pub mod annular;
pub mod annular_audit;
pub mod bezier;
pub mod bracket;
pub mod incidence;
pub mod interchange;
pub mod json;
pub mod planar;
pub use bezier::{Evaluation, MapError, RationalBezierMap};
