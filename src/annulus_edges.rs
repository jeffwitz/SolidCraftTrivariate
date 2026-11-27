//! Edge edits on the reference annulus, in the original or prepared volume partition.
//! The original four-chart quadratic fillet has singular boundary rails.
use crate::json::{object, Json};
mod linear;
mod original;
mod prepared;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum EdgeProfile {
    Sharp,
    Chamfer,
    Fillet,
}
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum EdgeLayout {
    OriginalLinear,
    OriginalQuadratic,
    PreparedTen,
    PreparedTwenty,
}

/// Return all knots and homogeneous control points in millimetres.
pub fn build_edge_model(layout: EdgeLayout, profile: EdgeProfile) -> Result<Json, &'static str> {
    match layout {
        EdgeLayout::OriginalLinear if profile == EdgeProfile::Fillet => {
            Err("linear axial spans cannot represent the circular fillet")
        }
        EdgeLayout::OriginalLinear => linear::document(profile == EdgeProfile::Chamfer),
        EdgeLayout::OriginalQuadratic => original::document(profile),
        EdgeLayout::PreparedTen => Ok(prepared::document(profile, 2)),
        EdgeLayout::PreparedTwenty => Ok(prepared::document(profile, 4)),
    }
}
fn string(s: &str) -> Json {
    Json::String(s.into())
}
fn numbers<T: ToString>(values: impl IntoIterator<Item = T>) -> Json {
    Json::Array(
        values
            .into_iter()
            .map(|n| Json::Number(n.to_string()))
            .collect(),
    )
}
fn points(net: Vec<[i64; 4]>) -> Json {
    Json::Array(net.into_iter().map(numbers).collect())
}
fn axial_knots(degree: usize) -> Json {
    let mut knots = vec![0.; degree + 1];
    for interior in [0.25, 7.75] {
        knots.extend(std::iter::repeat_n(interior, degree));
    }
    knots.extend(std::iter::repeat_n(8., degree + 1));
    numbers(knots)
}
fn original_document<'a>(
    degrees: [usize; 3],
    knots: [Json; 3],
    blocks: impl Iterator<Item = (&'a String, Vec<[i64; 4]>)>,
) -> Json {
    object([
        ("degrees", numbers(degrees)),
        ("knots", Json::Array(knots.into())),
        ("fastest_index", string("axial")),
        (
            "blocks",
            Json::Array(
                blocks
                    .map(|(id, net)| object([("id", string(id)), ("xyzw", points(net))]))
                    .collect(),
            ),
        ),
    ])
}
