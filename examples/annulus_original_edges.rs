use solidcraft_trivariate::annulus_edges::{build_edge_model, EdgeLayout, EdgeProfile};
fn main() {
    let state = match std::env::args().nth(1).as_deref() {
        Some("sharp") => EdgeProfile::Sharp,
        Some("chamfer") => EdgeProfile::Chamfer,
        Some("fillet") => EdgeProfile::Fillet,
        _ => panic!("expected sharp, chamfer or fillet"),
    };
    println!(
        "{}",
        build_edge_model(EdgeLayout::OriginalQuadratic, state).expect("valid annulus state")
    );
}
