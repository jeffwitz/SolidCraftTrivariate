//! Run one reproducible displacement qualification at a specified h level.

use solidcraft_trivariate::{
    analysis_space::build_analysis_space,
    bracket::bracket_blocks,
    bracket_solve::solve_bracket,
    h_refine::refine_uniform,
    json::{object, Json},
};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let [levels, order] = args.as_slice() else {
        return Err("usage: bracket_convergence LEVEL GAUSS_ORDER".into());
    };
    let levels = levels.parse::<usize>()?;
    let order = order.parse::<usize>()?;
    let begin = Instant::now();
    let mut space = build_analysis_space(&bracket_blocks()?)?;
    for _ in 0..levels {
        space = refine_uniform(&space)?;
    }
    let prep = begin.elapsed().as_secs_f64();
    let start = Instant::now();
    let result = solve_bracket(&space, order)?;
    let n = |v: f64| Json::Number(v.to_string());
    let integer = |v: usize| Json::Number(v.to_string());
    let vector = |v: [f64; 3]| Json::Array(v.map(n).into());
    let fillet = solidcraft_trivariate::stress::fillet_limits(&space, &result.displacement_mm)?;
    println!(
        "{}",
        object([
            ("level", integer(levels)),
            ("gauss_order", integer(order)),
            ("blocks", integer(space.blocks.len())),
            ("scalar_dofs", integer(space.ndof)),
            ("vector_dofs", integer(3 * space.ndof)),
            ("prepare_s", n(prep)),
            ("assemble_solve_s", n(start.elapsed().as_secs_f64())),
            ("volume_mm3", n(result.volume_mm3)),
            ("loaded_area_mm2", n(result.load_area_x60_mm2)),
            ("pcg_iterations", integer(result.pcg_iterations)),
            ("relative_residual", n(result.relative_residual)),
            ("force_n", vector(result.applied_force_n)),
            ("reaction_n", vector(result.reaction_force_n)),
            ("strain_energy_n_mm", n(result.strain_energy_n_mm)),
            ("external_work_n_mm", n(result.external_work_n_mm)),
            ("mean_tip_y_mm", n(result.mean_y_displacement_mm)),
            ("point_tip_y_mm", n(result.point_y_displacement_mm)),
            (
                "fillet_limits_mpa",
                Json::Object(fillet.into_iter().map(|(id, v)| (id, n(v))).collect())
            ),
        ])
    );
    Ok(())
}
