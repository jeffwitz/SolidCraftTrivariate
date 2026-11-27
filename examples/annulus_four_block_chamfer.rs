use solidcraft_trivariate::annulus_edges::{build_edge_model, EdgeLayout, EdgeProfile};
fn main() {
    let state = match std::env::args().nth(1).as_deref() {
        Some("sharp") => EdgeProfile::Sharp,
        Some("chamfer") | None => EdgeProfile::Chamfer,
        _ => panic!("expected sharp or chamfer"),
    };
    println!(
        "{}",
        build_edge_model(EdgeLayout::OriginalLinear, state).expect("valid annulus state")
    );
}
