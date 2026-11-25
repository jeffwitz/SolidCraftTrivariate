//! Measure the conforming Q2 inventories and local geometric effect of a radius edit.
use solidcraft_trivariate::{
    analysis_space::row_analysis_space,
    json::{object, Json},
    row_bores::build_row_bores,
};
use std::time::Instant;

fn main() -> Result<(), &'static str> {
    let mut rows = Vec::new();
    for count in [1, 2, 4, 8, 16] {
        let extents = [10 * count, 20];
        let centers: Vec<_> = (0..count).map(|i| [20 * i + 10 - extents[0], 0]).collect();
        let mut radii = vec![4; count as usize];
        let make = |r: &[i32]| {
            build_row_bores(extents[0], extents[1], &centers, r).map_err(|_| "invalid row")
        };
        let time = Instant::now();
        let initial = make(&radii)?;
        let space = row_analysis_space(extents, &centers, &radii, &initial.blocks)?;
        let elapsed = time.elapsed().as_secs_f64();
        radii[count as usize / 2] += 1;
        let updated = make(&radii)?;
        let edited_space = row_analysis_space(extents, &centers, &radii, &updated.blocks)?;
        assert_eq!(edited_space.ndof, space.ndof);
        let changed = initial
            .blocks
            .iter()
            .zip(&updated.blocks)
            .filter(|(a, b)| a != b)
            .count();
        assert_eq!(changed, 4);
        let n = |v: usize| Json::Number(v.to_string());
        rows.push(object([
            ("bores", n(count as usize)),
            ("blocks", n(initial.blocks.len())),
            ("changed_block_nets", n(changed)),
            ("scalar_dofs", n(space.ndof)),
            (
                "geometry_and_audited_space_s",
                Json::Number(elapsed.to_string()),
            ),
        ]));
    }
    println!("{}", Json::Array(rows));
    Ok(())
}
