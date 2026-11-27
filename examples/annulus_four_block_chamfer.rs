use solidcraft_trivariate::annulus_edges::{build_edge_model, EdgeLayout, EdgeProfile};
#[path = "support/annulus_legacy_output.rs"]
mod legacy_output;

fn main() {
    let state = match std::env::args().nth(1).as_deref() {
        Some("sharp") => EdgeProfile::Sharp,
        Some("chamfer") | None => EdgeProfile::Chamfer,
        _ => panic!("expected sharp or chamfer"),
    };
    let model = build_edge_model(EdgeLayout::OriginalLinear, state).expect("valid annulus state");
    legacy_output::print_model(&model);
}
