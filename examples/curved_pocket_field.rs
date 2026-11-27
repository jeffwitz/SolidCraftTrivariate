//! Solve one fixed-height field on a curved pocket described by an audited model.
use solidcraft_trivariate::{json::parse, pocket_convergence};
use std::{env, error::Error, fs};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = env::args().skip(1).collect();
    if (2..=3).contains(&args.len()) {
        let (model_path, job_path, intervals) = match args.as_slice() {
            [model, job] => (model, job, 1),
            [model, mode_or_intervals, job] if mode_or_intervals == "auto" => (model, job, 1),
            [model, intervals, job] => (model, job, intervals.parse()?),
            _ => return Err("usage: curved_pocket_field MODEL JOB [INTERVALS] | MODEL auto JOB [INTERVALS] | MODEL DEGREE INTERVALS C0|C1|radial|footprint JOB".into()),
        };
        let model = parse(&fs::read_to_string(model_path)?)?;
        let job = parse(&fs::read_to_string(job_path)?)?;
        println!(
            "{}",
            pocket_convergence::solve_model_auto(&model, intervals, &job)?
        );
        return Ok(());
    }
    let [model, degree, intervals, continuity, job] = args.as_slice() else {
        return Err("usage: curved_pocket_field MODEL JOB | MODEL INTERVALS JOB | MODEL auto JOB | MODEL DEGREE INTERVALS C0|C1|radial|footprint JOB".into());
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
