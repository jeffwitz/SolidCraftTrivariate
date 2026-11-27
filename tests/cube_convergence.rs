use solidcraft_trivariate::{
    cube_convergence,
    json::{object, Json},
};
use std::{fs, path::Path};

fn number(value: &Json) -> f64 {
    value.number().unwrap()
}
fn usize_number(value: &Json) -> usize {
    value.number().unwrap()
}

fn save_report(name: &str, value: &Json) {
    let path = Path::new("build/cube-convergence").join(name);
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    fs::write(path, value.to_string()).unwrap();
}

fn measure(row: &Json, key: &str) -> f64 {
    number(row.field("result").unwrap().field(key).unwrap())
}

fn rate(a: &Json, b: &Json, key: &str) -> f64 {
    (measure(a, key) / measure(b, key)).ln() / 2_f64.ln()
}

fn band_verdict(energy: f64, l2: f64) -> &'static str {
    if (2.4..=3.6).contains(&energy) && (3.2..=4.8).contains(&l2) {
        "pass"
    } else {
        "fail"
    }
}

#[test]
fn cube_extraction_smoke_has_complete_affine_parent() {
    let report = cube_convergence::run(&[1], 4).unwrap();
    let rows = report.field("observations").unwrap().items().unwrap();
    assert_eq!(rows.len(), 2);
    for row in rows {
        assert_eq!(usize_number(row.field("cell_count").unwrap()), 1);
        assert_eq!(usize_number(row.field("scalar_dofs").unwrap()), 64);
        assert!(number(row.field("max_geometry_position_error_mm").unwrap()) < 2e-12);
        assert!(number(row.field("max_jacobian_error").unwrap()) < 2e-12);
        assert!(number(row.field("max_affine_value_error_mm").unwrap()) < 2e-11);
        assert!(number(row.field("max_affine_gradient_error").unwrap()) < 2e-11);
        assert!(
            number(
                row.field("result")
                    .unwrap()
                    .field("relative_residual")
                    .unwrap()
            ) < 1e-9
        );
    }
}

#[test]
#[ignore = "preregistered q8/q12 convergence gate; run optimized with --ignored"]
fn uniform_cubic_cube_meets_preregistered_convergence_and_quadrature_bands() {
    let q8 = cube_convergence::run(&[1, 2, 4], 8).unwrap();
    save_report("q8-baseline.json", &q8);
    let q12 = cube_convergence::run(&[4], 12).unwrap();
    save_report("q12-fine.json", &q12);
    let rows = q8.field("observations").unwrap().items().unwrap();
    assert_eq!(rows.len(), 6);
    for row in rows {
        let cells = usize_number(row.field("cells_per_axis").unwrap());
        let ndof = usize_number(row.field("scalar_dofs").unwrap());
        let continuity = usize_number(row.field("interior_continuity").unwrap());
        let expected = if continuity == 2 {
            cells + 3
        } else {
            3 * cells + 1
        };
        assert_eq!(ndof, expected.pow(3));
        assert_eq!(usize_number(row.field("cell_count").unwrap()), cells.pow(3));
        assert!(number(row.field("min_sampled_determinant").unwrap()) > 0.);
        assert!(number(row.field("max_geometry_position_error_mm").unwrap()) < 2e-12);
        assert!(number(row.field("max_jacobian_error").unwrap()) < 2e-12);
        assert!(number(row.field("max_affine_value_error_mm").unwrap()) < 2e-11);
        assert!(number(row.field("max_affine_gradient_error").unwrap()) < 2e-11);
        let result = row.field("result").unwrap();
        assert!((number(result.field("volume_mm3").unwrap()) - 1.).abs() < 2e-12);
        assert!(number(result.field("relative_residual").unwrap()) < 1e-9);
        assert!(number(result.field("virtual_work_relative_defect").unwrap()) < 1e-8);
    }
    for continuity in [2usize, 0] {
        let subset: Vec<_> = rows
            .iter()
            .filter(|r| usize_number(r.field("interior_continuity").unwrap()) == continuity)
            .collect();
        let l2: Vec<_> = subset
            .iter()
            .map(|r| {
                number(
                    r.field("result")
                        .unwrap()
                        .field("relative_displacement_L2")
                        .unwrap(),
                )
            })
            .collect();
        let energy: Vec<_> = subset
            .iter()
            .map(|r| {
                number(
                    r.field("result")
                        .unwrap()
                        .field("relative_energy_norm")
                        .unwrap(),
                )
            })
            .collect();
        assert!(l2[0] > l2[1] && l2[1] > l2[2], "L2 must decrease: {l2:?}");
        assert!(
            energy[0] > energy[1] && energy[1] > energy[2],
            "energy error must decrease: {energy:?}"
        );
        let energy_order = (energy[1] / energy[2]).ln() / 2_f64.ln();
        let l2_order = (l2[1] / l2[2]).ln() / 2_f64.ln();
        assert!(
            (2.4..=3.6).contains(&energy_order),
            "energy order {energy_order}"
        );
        assert!((3.2..=4.8).contains(&l2_order), "L2 order {l2_order}");
    }
    let c2fine = &rows[2];
    let c0fine = &rows[5];
    assert!(
        usize_number(c2fine.field("scalar_dofs").unwrap())
            < usize_number(c0fine.field("scalar_dofs").unwrap())
    );
    let error_ratio = number(
        c2fine
            .field("result")
            .unwrap()
            .field("relative_energy_norm")
            .unwrap(),
    ) / number(
        c0fine
            .field("result")
            .unwrap()
            .field("relative_energy_norm")
            .unwrap(),
    );
    assert!(
        error_ratio <= 1.5,
        "C2/C0 fine energy-error ratio {error_ratio}"
    );

    let q12rows = q12.field("observations").unwrap().items().unwrap();
    assert_eq!(q12rows.len(), 2);
    for (a, b) in [(&rows[2], &q12rows[0]), (&rows[5], &q12rows[1])] {
        let qa = a.field("result").unwrap();
        let qb = b.field("result").unwrap();
        for key in ["relative_displacement_L2", "relative_energy_norm"] {
            let x = number(qa.field(key).unwrap());
            let y = number(qb.field(key).unwrap());
            assert!(
                (x - y).abs() / x.abs().max(1e-30) <= 2e-8,
                "{key}: q8={x}, q12={y}"
            );
        }
        assert!(
            (number(qa.field("volume_mm3").unwrap()) - number(qb.field("volume_mm3").unwrap()))
                .abs()
                <= 2e-12
        );
    }
}

#[test]
#[ignore = "optimized numerical replay and n=8 asymptotic supplement"]
fn cube_convergence_replay() {
    let q8 = cube_convergence::run(&[1, 2, 4], 8).unwrap();
    let q12 = cube_convergence::run(&[4], 12).unwrap();
    let fine = cube_convergence::run(&[4, 8], 8).unwrap();
    let rows = q8.field("observations").unwrap().items().unwrap();
    let q12rows = q12.field("observations").unwrap().items().unwrap();
    let finerows = fine.field("observations").unwrap().items().unwrap();
    assert_eq!(rows.len(), 6);
    assert_eq!(q12rows.len(), 2);
    assert_eq!(finerows.len(), 4);

    // Freeze every original row so the replay cannot silently substitute a
    // different manufactured problem, field partition, or reported error.
    let expected = [
        (64, 0.11525539675419867, 0.046643965123194145),
        (125, 0.016467473233362222, 0.003411203434748927),
        (343, 0.0031451460233471525, 0.00044404581226912827),
        (64, 0.11525539675419867, 0.046643965123194145),
        (343, 0.01161815128700155, 0.0021770625980608604),
        (2197, 0.0015007567420937593, 0.00013316929930979986),
    ];
    for (row, (dofs, energy, l2)) in rows.iter().zip(expected) {
        assert_eq!(usize_number(row.field("scalar_dofs").unwrap()), dofs);
        assert!((measure(row, "relative_energy_norm") - energy).abs() < 2e-10);
        assert!((measure(row, "relative_displacement_L2") - l2).abs() < 2e-10);
        let cells = usize_number(row.field("cells_per_axis").unwrap());
        assert_eq!(usize_number(row.field("cell_count").unwrap()), cells.pow(3));
        assert!(number(row.field("min_sampled_determinant").unwrap()) > 0.);
        assert!(number(row.field("max_geometry_position_error_mm").unwrap()) < 2e-12);
        assert!(number(row.field("max_jacobian_error").unwrap()) < 2e-12);
        assert!(number(row.field("max_affine_value_error_mm").unwrap()) < 2e-11);
        assert!(number(row.field("max_affine_gradient_error").unwrap()) < 2e-11);
        let result = row.field("result").unwrap();
        assert!((number(result.field("volume_mm3").unwrap()) - 1.).abs() < 2e-12);
        assert!(number(result.field("relative_residual").unwrap()) < 1e-9);
        assert!(number(result.field("virtual_work_relative_defect").unwrap()) < 1e-8);
    }

    let primary_c2_energy = rate(&rows[1], &rows[2], "relative_energy_norm");
    let primary_c2_l2 = rate(&rows[1], &rows[2], "relative_displacement_L2");
    let primary_c0_energy = rate(&rows[4], &rows[5], "relative_energy_norm");
    let primary_c0_l2 = rate(&rows[4], &rows[5], "relative_displacement_L2");
    assert!((primary_c2_energy - 2.388420299683446).abs() < 2e-8);
    assert!((primary_c2_l2 - 2.9415).abs() < 0.01);
    assert_eq!(band_verdict(primary_c2_energy, primary_c2_l2), "fail");
    assert!((primary_c0_energy - 2.9526184669557263).abs() < 2e-8);
    assert!((primary_c0_l2 - 4.031049461908917).abs() < 2e-8);
    assert_eq!(band_verdict(primary_c0_energy, primary_c0_l2), "pass");
    let primary_energy_ratio =
        measure(&rows[2], "relative_energy_norm") / measure(&rows[5], "relative_energy_norm");
    assert!((primary_energy_ratio - 2.0957).abs() < 0.002);

    let fine_c2_energy = rate(&finerows[0], &finerows[1], "relative_energy_norm");
    let fine_c2_l2 = rate(&finerows[0], &finerows[1], "relative_displacement_L2");
    let fine_c0_energy = rate(&finerows[2], &finerows[3], "relative_energy_norm");
    let fine_c0_l2 = rate(&finerows[2], &finerows[3], "relative_displacement_L2");
    for row in finerows {
        assert!(
            number(
                row.field("result")
                    .unwrap()
                    .field("relative_residual")
                    .unwrap()
            ) < 1e-9
        );
        assert!(
            (number(row.field("result").unwrap().field("volume_mm3").unwrap()) - 1.).abs() < 2e-12
        );
    }
    let supplement_c2_ratio = measure(&finerows[1], "relative_energy_norm")
        / measure(&finerows[3], "relative_energy_norm");
    assert!(
        supplement_c2_ratio > 1.5,
        "retain observed C2/C0 discrepancy"
    );
    let supplement_c2_energy_band = (2.4..=3.6).contains(&fine_c2_energy);
    let supplement_c2_l2_band = (3.2..=4.8).contains(&fine_c2_l2);
    let supplement_c0_energy_band = (2.4..=3.6).contains(&fine_c0_energy);
    let supplement_c0_l2_band = (3.2..=4.8).contains(&fine_c0_l2);

    let mut q12_defects = Vec::new();
    for (a, b) in [(&rows[2], &q12rows[0]), (&rows[5], &q12rows[1])] {
        for key in ["relative_displacement_L2", "relative_energy_norm"] {
            let x = measure(a, key);
            let y = measure(b, key);
            let defect = (x - y).abs() / x.abs().max(1e-30);
            assert!(defect <= 2e-8, "q8/q12 {key} relative defect {defect}");
            q12_defects.push(object_for_pair(key, defect));
        }
        assert!(
            (number(a.field("result").unwrap().field("volume_mm3").unwrap())
                - number(b.field("result").unwrap().field("volume_mm3").unwrap()))
            .abs()
                < 2e-12
        );
    }
    let log_covariance = |indices: &[usize], key: &str| {
        let xs: Vec<f64> = indices
            .iter()
            .map(|&i| (usize_number(rows[i].field("scalar_dofs").unwrap()) as f64).ln())
            .collect();
        let ys: Vec<f64> = indices
            .iter()
            .map(|&i| measure(&rows[i], key).ln())
            .collect();
        let xm = xs.iter().sum::<f64>() / xs.len() as f64;
        let ym = ys.iter().sum::<f64>() / ys.len() as f64;
        xs.iter()
            .zip(&ys)
            .map(|(x, y)| (x - xm) * (y - ym))
            .sum::<f64>()
            / (xs.len() - 1) as f64
    };
    let analytic = object([
        (
            "body_axial_resultant_N",
            Json::Number((7. * std::f64::consts::PI).to_string()),
        ),
        (
            "top_traction_axial_resultant_N",
            Json::Number((-3.5 * std::f64::consts::PI).to_string()),
        ),
        (
            "bottom_support_reaction_N",
            Json::Number((-3.5 * std::f64::consts::PI).to_string()),
        ),
        (
            "exact_strain_energy_N_mm",
            Json::Number((7. * std::f64::consts::PI.powi(2) / 8.).to_string()),
        ),
    ]);
    let report = object([
        ("schema", Json::String("cube_convergence_replay_v1".into())),
        ("primary_q8", q8.clone()),
        ("q12_finest_quadrature", q12.clone()),
        ("supplement_q8_levels_4_to_8", fine.clone()),
        ("analytic_diagnostics", analytic),
        (
            "primary_verdicts",
            object_for_rates(
                primary_c2_energy,
                primary_c2_l2,
                primary_c0_energy,
                primary_c0_l2,
            ),
        ),
        (
            "supplement_verdicts",
            object_for_rates(fine_c2_energy, fine_c2_l2, fine_c0_energy, fine_c0_l2),
        ),
        ("q8_q12_relative_defects", Json::Array(q12_defects)),
        (
            "acceptance_deviations",
            Json::Array(vec![
                Json::String("primary C2 rate band failed at n=2 to n=4".into()),
                Json::String("primary C2/C0 fine energy error ratio exceeded 1.5".into()),
                Json::String("initial strict gate rejected the frozen observations".into()),
            ]),
        ),
        (
            "primary_C2_to_C0_energy_error_ratio_n4",
            Json::Number(primary_energy_ratio.to_string()),
        ),
        (
            "supplement_rate_bands",
            object([
                ("C2_energy", Json::Bool(supplement_c2_energy_band)),
                ("C2_displacement_L2", Json::Bool(supplement_c2_l2_band)),
                ("C0_energy", Json::Bool(supplement_c0_energy_band)),
                ("C0_displacement_L2", Json::Bool(supplement_c0_l2_band)),
            ]),
        ),
        (
            "log_dof_error_sample_covariance",
            object([
                (
                    "C2_energy",
                    Json::Number(log_covariance(&[0, 1, 2], "relative_energy_norm").to_string()),
                ),
                (
                    "C2_displacement_L2",
                    Json::Number(
                        log_covariance(&[0, 1, 2], "relative_displacement_L2").to_string(),
                    ),
                ),
                (
                    "C0_energy",
                    Json::Number(log_covariance(&[3, 4, 5], "relative_energy_norm").to_string()),
                ),
                (
                    "C0_displacement_L2",
                    Json::Number(
                        log_covariance(&[3, 4, 5], "relative_displacement_L2").to_string(),
                    ),
                ),
            ]),
        ),
    ]);
    save_report("replay.json", &report);
    println!("Cube convergence replay: primary C2={primary_c2_energy:.6}/{primary_c2_l2:.6} ({}) C0={primary_c0_energy:.6}/{primary_c0_l2:.6} ({}); supplement C2={fine_c2_energy:.6}/{fine_c2_l2:.6} ({}) C0={fine_c0_energy:.6}/{fine_c0_l2:.6} ({})",
        band_verdict(primary_c2_energy,primary_c2_l2),band_verdict(primary_c0_energy,primary_c0_l2),
        band_verdict(fine_c2_energy,fine_c2_l2),band_verdict(fine_c0_energy,fine_c0_l2));
}

fn object_for_rates(c2e: f64, c2l: f64, c0e: f64, c0l: f64) -> Json {
    object([
        (
            "C2_maximum",
            object([
                ("energy_order", Json::Number(c2e.to_string())),
                ("displacement_L2_order", Json::Number(c2l.to_string())),
                ("verdict", Json::String(band_verdict(c2e, c2l).into())),
            ]),
        ),
        (
            "C0_bernstein",
            object([
                ("energy_order", Json::Number(c0e.to_string())),
                ("displacement_L2_order", Json::Number(c0l.to_string())),
                ("verdict", Json::String(band_verdict(c0e, c0l).into())),
            ]),
        ),
    ])
}

fn object_for_pair(key: &str, value: f64) -> Json {
    object([
        ("quantity", Json::String(key.into())),
        ("relative_defect", Json::Number(value.to_string())),
    ])
}
