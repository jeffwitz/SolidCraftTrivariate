use solidcraft_trivariate::{
    crossed_bending,
    g2_cross_plate::boolean,
    json::{parse, Json},
};
use std::{fs, process::Command};

fn model(rise_mm: i64) -> Json {
    let mut request = parse(include_str!("jobs/crossed_partial_rise3.recipe.json")).unwrap();
    let Json::Object(fields) = &mut request else {
        unreachable!()
    };
    fields.insert(
        "transverse_rise_mm".into(),
        Json::Number(rise_mm.to_string()),
    );
    boolean::report(&request, "union").unwrap()
}

fn oracle() -> Vec<Json> {
    parse(include_str!(
        "../docs/results/crossed-interface-iga/partial_bending_study.json"
    ))
    .unwrap()
    .items()
    .unwrap()
    .to_vec()
}

fn reference(rise: i64, nv: usize, nw: usize, quadrature: usize) -> Json {
    oracle()
        .into_iter()
        .find(|row| {
            row.field("rise_mm").unwrap().number::<i64>().unwrap() == rise
                && row.field("refine_v").unwrap().number::<usize>().unwrap() == nv
                && row.field("refine_w").unwrap().number::<usize>().unwrap() == nw
                && row
                    .field("quadrature_order")
                    .unwrap()
                    .number::<usize>()
                    .unwrap()
                    == quadrature
        })
        .expect("published independent S2 mechanics case")
}

fn close(actual: &Json, expected: &Json, relative: f64, absolute: f64) {
    let a = actual.number::<f64>().unwrap();
    let e = expected.number::<f64>().unwrap();
    assert!(
        (a - e).abs() <= absolute.max(relative * e.abs()),
        "{a} != {e}"
    );
}

fn check_case(model: &Json, nv: usize, nw: usize, q: usize) -> Json {
    let result = crossed_bending::analyze(model, nv, nw, q).unwrap();
    let expected = reference(3, nv, nw, q);
    for key in [
        "scalar_dofs_after_full_identifications",
        "partial_trace_rank",
        "free_vector_dofs",
    ] {
        assert_eq!(
            result.field(key).unwrap().number::<usize>().unwrap(),
            expected.field(key).unwrap().number::<usize>().unwrap(),
            "{key} at v={nv}, w={nw}, q={q}"
        );
    }
    close(
        result.field("loaded_end_mean_ux_mm").unwrap(),
        expected.field("loaded_end_mean_ux_mm").unwrap(),
        3e-8,
        2e-10,
    );
    close(
        result.field("volume_mm3").unwrap(),
        expected.field("volume_mm3").unwrap(),
        2e-10,
        2e-9,
    );
    close(
        result.field("force_N").unwrap(),
        expected.field("force_N").unwrap(),
        2e-10,
        2e-9,
    );
    assert!(
        result
            .field("affine_energy_error")
            .unwrap()
            .number::<f64>()
            .unwrap()
            < 1e-8
    );
    assert!(
        result
            .field("reduced_equilibrium_residual")
            .unwrap()
            .number::<f64>()
            .unwrap()
            < 4e-10
    );
    assert!(
        result
            .field("max_partial_face_geometry_mismatch_mm")
            .unwrap()
            .number::<f64>()
            .unwrap()
            < 1e-12
    );
    assert!(
        result
            .field("reduced_symmetry_relative")
            .unwrap()
            .number::<f64>()
            .unwrap()
            < 1e-12
    );
    result
}

fn check_case_for_rise(model: &Json, rise: i64, nv: usize, nw: usize, q: usize) -> Json {
    let result = crossed_bending::analyze(model, nv, nw, q).unwrap();
    let expected = reference(rise, nv, nw, q);
    for key in [
        "scalar_dofs_after_full_identifications",
        "partial_trace_rank",
        "free_vector_dofs",
    ] {
        assert_eq!(
            result.field(key).unwrap().number::<usize>().unwrap(),
            expected.field(key).unwrap().number::<usize>().unwrap(),
            "{key} at rise={rise}, v={nv}, w={nw}, q={q}"
        );
    }
    close(
        result.field("loaded_end_mean_ux_mm").unwrap(),
        expected.field("loaded_end_mean_ux_mm").unwrap(),
        3e-8,
        2e-10,
    );
    for key in ["volume_mm3", "force_N", "strain_energy_Nmm"] {
        let relative = if key == "strain_energy_Nmm" {
            2e-9
        } else {
            2e-10
        };
        close(
            result.field(key).unwrap(),
            expected.field(key).unwrap(),
            relative,
            2e-9,
        );
    }
    assert!(
        result
            .field("reduced_equilibrium_residual")
            .unwrap()
            .number::<f64>()
            .unwrap()
            < 4e-10
    );
    assert!(
        result
            .field("max_partial_face_geometry_mismatch_mm")
            .unwrap()
            .number::<f64>()
            .unwrap()
            < 1e-12
    );
    assert!(
        result
            .field("reduced_symmetry_relative")
            .unwrap()
            .number::<f64>()
            .unwrap()
            < 1e-12
    );
    result
}

#[test]
fn native_p322_reproduces_the_published_s2_refinement_and_quadrature_controls() {
    let model = model(3);
    let coarse = check_case(&model, 1, 1, 8);
    let fine_q8 = check_case(&model, 4, 1, 8);
    let fine_q14 = check_case(&model, 4, 1, 14);
    let u8 = fine_q8
        .field("loaded_end_mean_ux_mm")
        .unwrap()
        .number::<f64>()
        .unwrap();
    let u14 = fine_q14
        .field("loaded_end_mean_ux_mm")
        .unwrap()
        .number::<f64>()
        .unwrap();
    assert!((u8 - u14).abs() < 1e-8);
    assert!(
        coarse
            .field("free_vector_dofs")
            .unwrap()
            .number::<usize>()
            .unwrap()
            < fine_q8
                .field("free_vector_dofs")
                .unwrap()
                .number::<usize>()
                .unwrap()
    );
    assert_eq!(
        fine_q8
            .field("symmetry_check_method")
            .unwrap()
            .text()
            .unwrap(),
        "six_matrix_free_bilinear_probes"
    );
}

#[test]
fn native_p322_retains_the_two_published_w_refinement_regressions() {
    let model = model(3);
    for (nv, nw) in [(1, 2), (2, 2)] {
        check_case_for_rise(&model, 3, nv, nw, 8);
    }
}

#[test]
fn native_p322_refuses_other_solid_families_and_unsupported_refinements() {
    let mut other = model(3);
    let Json::Object(fields) = &mut other else {
        unreachable!()
    };
    fields.insert("operation".into(), Json::String("intersection".into()));
    assert!(crossed_bending::analyze(&other, 1, 1, 8).is_err());
    let valid = model(3);
    assert!(crossed_bending::analyze(&valid, 5, 1, 8).is_err());
    assert!(crossed_bending::analyze(&valid, 1, 1, 15).is_err());
}

#[test]
fn native_cli_replays_source_and_preserves_output_after_malformed_geometry() {
    let root = std::env::temp_dir().join(format!(
        "solidcraft-crossed-bending-cli-{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let input = root.join("model.json");
    let output = root.join("result.json");
    let valid = model(3);
    fs::write(&input, valid.to_string()).unwrap();
    let ok = Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
        .args([
            "crossed-bending",
            input.to_str().unwrap(),
            "1",
            "1",
            "8",
            output.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        ok.status.success(),
        "{}",
        String::from_utf8_lossy(&ok.stdout)
    );
    let cli_result = parse(&fs::read_to_string(&output).unwrap()).unwrap();
    let expected = reference(3, 1, 1, 8);
    close(
        cli_result.field("loaded_end_mean_ux_mm").unwrap(),
        expected.field("loaded_end_mean_ux_mm").unwrap(),
        3e-8,
        2e-10,
    );

    let mut malformed = valid;
    let Json::Object(fields) = &mut malformed else {
        unreachable!()
    };
    fields.insert("blocks".into(), Json::Array(Vec::new()));
    fs::write(&input, malformed.to_string()).unwrap();
    let sentinel = b"preserve this prior result\n";
    fs::write(&output, sentinel).unwrap();
    let failed = Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
        .args([
            "crossed-bending",
            input.to_str().unwrap(),
            "1",
            "1",
            "8",
            output.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert_eq!(fs::read(&output).unwrap(), sentinel);
    fs::remove_dir_all(root).unwrap();
}
