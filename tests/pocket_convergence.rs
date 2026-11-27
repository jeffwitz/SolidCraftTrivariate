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
