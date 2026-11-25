//! Native geometry on the closed unit cube. Physical coordinates are in mm.

pub mod analysis_space;
pub mod annular;
pub mod annular_audit;
pub mod bezier;
pub mod bracket;
pub mod bracket_audit;
pub mod bracket_solve;
mod bracket_json;
pub mod incidence;
pub mod h_refine;
pub mod elasticity;
pub mod constraint_system;
pub mod face_constraints;
pub mod face_partition;
pub mod interchange;
pub mod json;
pub mod planar;
pub mod rational_sweep;
pub mod sections;
pub mod stress;
pub use bezier::{Evaluation, MapError, RationalBezierMap};
