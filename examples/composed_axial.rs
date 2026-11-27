use solidcraft_trivariate::composed_axial::compose;

fn main() {
    let result = compose([0, 8], [4, 12], [5, 9]).expect("composed axial boolean");
    println!(
        "first intervals={:?}; second intervals={:?}; volume_mm3={:.12}; carried maps={}; face lineage={:?}",
        result.first.selected_intervals,
        result.second.selected_intervals,
        result.second.volume_mm3,
        result.carried_maps.len(),
        result.lineage
    );
}
