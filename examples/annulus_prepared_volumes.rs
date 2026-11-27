use solidcraft_trivariate::annulus_edges::{build_edge_model, EdgeLayout, EdgeProfile};
#[path = "support/annulus_legacy_output.rs"]
mod legacy_output;

fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let state = match args.first().map(String::as_str) {
        Some("sharp") => EdgeProfile::Sharp,
        Some("chamfer") => EdgeProfile::Chamfer,
        Some("fillet") | None => EdgeProfile::Fillet,
        _ => panic!("expected sharp, chamfer or fillet"),
    };
    let sectors = args.get(1).map_or(2, |s| s.parse::<usize>().unwrap());
    let model = build_edge_model(
        match sectors {
            2 => EdgeLayout::PreparedTen,
            4 => EdgeLayout::PreparedTwenty,
            _ => panic!("expected two or four sectors"),
        },
        state,
    )
    .expect("valid annulus state");
    legacy_output::print_model(&model);
}
