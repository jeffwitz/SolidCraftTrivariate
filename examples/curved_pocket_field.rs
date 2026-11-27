//! Solve one fixed-height field on a curved pocket described by an audited model.
use solidcraft_trivariate::{json::parse, pocket_convergence};
use std::{env, error::Error, fs};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = env::args().skip(1).collect();
    let [model, degree, intervals, continuity, job] = args.as_slice() else {
        return Err(
            "usage: curved_pocket_field MODEL DEGREE INTERVALS C0|C1|radial|footprint JOB".into(),
        );
    };
    use pocket_convergence::SectionJoins;
    let joins = match continuity.as_str() {
        "C0" => SectionJoins::None,
        "radial" => SectionJoins::Radial,
        "footprint" => SectionJoins::Footprint,
        "C1" => SectionJoins::All,
        _ => return Err("unknown seam selection".into()),
    };
    let result = pocket_convergence::solve_model_selected(
        &parse(&fs::read_to_string(model)?)?,
        degree.parse()?,
        intervals.parse()?,
        joins,
        &parse(&fs::read_to_string(job)?)?,
    )?;
    println!("{result}");
    Ok(())
}
