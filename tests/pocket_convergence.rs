use solidcraft_trivariate::{json::Json, pocket_convergence};
use std::{fs, path::Path};

fn measure(p: usize, n: usize, c1: bool, q: usize) -> Json {
    let row = pocket_convergence::solve(p, n, c1, q).unwrap();
    let get = |key: &str| row.field(key).unwrap().number::<f64>().unwrap();
    assert!((get("volume_mm3") / 46800. - 1.).abs() < 1e-10);
    assert!(get("relative_residual") < 1e-9);
    assert!(get("virtual_work_relative_defect") < 1e-7);
    assert!(get("relative_energy_norm") > 0.);
    row
}

fn save(p: usize, n: usize, c1: bool, q: usize) {
    let row = measure(p, n, c1, q);
    let dir = Path::new("build/pocket-vector");
    fs::create_dir_all(dir).unwrap();
    fs::write(
        dir.join(format!("p{p}-n{n}-c{}-q{q}.json", usize::from(c1))),
        row.to_string(),
    )
    .unwrap();
}

fn affine_power_job(degree: usize, quadrature: usize) -> Json {
    solidcraft_trivariate::json::parse(&format!(
        "{{\"schema\":\"section_affine_power_mms_v1\",\"quadrature\":{quadrature},\"coefficients\":[1,2,0],\"offset\":10,\"scale\":150,\"degree\":{degree}}}"
    ))
    .unwrap()
}

fn save_control(name: &str, p: usize, c1: bool, degree: usize, q: usize) -> Json {
    let row = pocket_convergence::solve_published(p, 1, c1, &affine_power_job(degree, q)).unwrap();
    let dir = Path::new("build/curved-pocket-vector/controls");
    fs::create_dir_all(dir).unwrap();
    fs::write(dir.join(format!("{name}-q{q}.json")), row.to_string()).unwrap();
    row
}

fn save_curved(p: usize, n: usize, joins: pocket_convergence::SectionJoins, q: usize) {
    let row = pocket_convergence::solve_published(
        p,
        n,
        joins != pocket_convergence::SectionJoins::None,
        &affine_power_job(5, q),
    )
    .unwrap();
    let dir = Path::new("build/curved-pocket-vector");
    fs::create_dir_all(dir).unwrap();
    fs::write(
        dir.join(format!("p{p}-n{n}-{joins:?}-q{q}.json").to_lowercase()),
        row.to_string(),
    )
    .unwrap();
}

#[test]
fn pocket_quartic_jet_constraints_reduce_dimension() {
    let a = measure(4, 1, false, 8);
    let b = measure(4, 1, true, 8);
    let dim = |row: &Json| row.field("scalar_dofs").unwrap().number::<usize>().unwrap();
    assert!(dim(&a) > dim(&b));
    assert!(pocket_convergence::solve(4, 0, true, 8).is_err());
}

#[test]
#[ignore = "release vector h-refinement and quadrature replication"]
fn vector_section_sweep() {
    for n in [1, 2, 4, 8] {
        for p in [3, 4] {
            for c1 in [false, true] {
                save(p, n, c1, 8);
                if n >= 4 {
                    save(p, n, c1, 10);
                }
            }
        }
    }
}

#[test]
fn cubic_coarse_partition_matches_the_existing_native_extraction() {
    use solidcraft_trivariate::{inherited_space, interchange, json::parse, section_manufactured};
    let recipe = parse(include_str!("references/blind_bilinear_pocket.recipe.json")).unwrap();
    let model = interchange::build(&recipe).unwrap();
    let reference = inherited_space::build_plane_space(&model, 1, 1, false, &[], 3).unwrap();
    let job = parse(
        r#"{"schema":"section_quintic_mms_v1","quadrature":8,"x_origin_mm":-10,"x_scale_mm":70}"#,
    )
    .unwrap();
    let old = section_manufactured::solve(&model, &reference, &job).unwrap();
    let new = measure(3, 1, false, 8);
    for key in [
        "free_vector_dofs",
        "relative_energy_norm",
        "relative_displacement_L2",
        "energy_N_mm",
        "volume_mm3",
    ] {
        let x = old.field(key).unwrap().number::<f64>().unwrap();
        let y = new.field(key).unwrap().number::<f64>().unwrap();
        assert!((x / y - 1.).abs() < 1e-8, "{key}: {x} vs {y}");
    }
}

#[test]
fn published_biquadratic_pocket_geometry_is_used_and_affine_fields_are_reproduced() {
    for (p, degree) in [(3, 1), (4, 2)] {
        for c1 in [false, true] {
            let result =
                pocket_convergence::solve_published(p, 1, c1, &affine_power_job(degree, 8))
                    .unwrap();
            let get = |key: &str| result.field(key).unwrap().number::<f64>().unwrap();
            assert!((get("volume_mm3") / 48976. - 1.).abs() < 1e-10);
            assert!(
                get("relative_energy_norm") < 1e-8,
                "p={p} C1={c1}: {}",
                get("relative_energy_norm")
            );
            assert!(
                get("relative_displacement_L2") < 1e-8,
                "p={p} C1={c1}: {}",
                get("relative_displacement_L2")
            );
            assert_eq!(
                result
                    .field("chart_count")
                    .unwrap()
                    .number::<usize>()
                    .unwrap(),
                5
            );
        }
    }
}

#[test]
fn published_curved_field_control_loads_are_recorded_at_q10() {
    for (name, p, c1, degree) in [
        ("affine-p3-c0", 3, false, 1),
        ("affine-p3-c1", 3, true, 1),
        ("quadratic-p4-c0", 4, false, 2),
        ("quadratic-p4-c1", 4, true, 2),
    ] {
        let row = save_control(name, p, c1, degree, 10);
        assert!(
            row.field("relative_residual")
                .unwrap()
                .number::<f64>()
                .unwrap()
                < 1e-9
        );
        assert!(
            row.field("virtual_work_relative_defect")
                .unwrap()
                .number::<f64>()
                .unwrap()
                < 1e-7
        );
    }
}

#[test]
fn published_curved_radial_and_footprint_c1_are_separate_diagnostics() {
    use solidcraft_trivariate::json::parse;
    let model = parse(include_str!(
        "../docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"
    ))
    .unwrap();
    let job = affine_power_job(5, 10);
    let dir = Path::new("build/curved-pocket-vector/selective");
    fs::create_dir_all(dir).unwrap();
    for joins in [
        pocket_convergence::SectionJoins::Radial,
        pocket_convergence::SectionJoins::Footprint,
    ] {
        let row = pocket_convergence::solve_model_selected(&model, 3, 1, joins, &job).unwrap();
        fs::write(
            dir.join(format!("p3-n1-{joins:?}-q10.json").to_lowercase()),
            row.to_string(),
        )
        .unwrap();
    }
}

#[test]
fn published_curved_primary_p4_coarse_q12_is_stable() {
    let value = |row: &Json, key: &str| row.field(key).unwrap().number::<f64>().unwrap();
    let dir = Path::new("build/curved-pocket-vector");
    fs::create_dir_all(dir).unwrap();
    for c1 in [false, true] {
        let rows: Vec<_> = [8, 10, 12]
            .into_iter()
            .map(|q| {
                pocket_convergence::solve_published(4, 1, c1, &affine_power_job(5, q)).unwrap()
            })
            .collect();
        let e10 = value(&rows[1], "relative_energy_norm");
        let e12 = value(&rows[2], "relative_energy_norm");
        assert!((e10 - e12).abs() < 1e-7);
        assert!(
            (value(&rows[0], "energy_N_mm") / value(&rows[2], "energy_N_mm") - 1.).abs() < 1e-7
        );
        let name = if c1 { "all" } else { "none" };
        fs::write(
            dir.join(format!("p4-n1-{name}-primary-q12.json")),
            rows[2].to_string(),
        )
        .unwrap();
    }
}

#[test]
#[ignore = "full published biquadratic C0/C1 h and quadrature sweep"]
fn published_curved_vector_sweep() {
    for n in [1, 2, 4, 8] {
        for p in [3, 4] {
            for joins in [
                pocket_convergence::SectionJoins::None,
                pocket_convergence::SectionJoins::All,
            ] {
                for q in [8, 10] {
                    save_curved(p, n, joins, q);
                }
            }
        }
    }
    for joins in [
        pocket_convergence::SectionJoins::None,
        pocket_convergence::SectionJoins::All,
    ] {
        save_curved(4, 1, joins, 12);
    }
}
