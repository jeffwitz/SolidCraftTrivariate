use solidcraft_trivariate::{
    analysis_space::{AnalysisBlock, AnalysisSpace},
    elasticity::assemble as assemble_q2,
    extracted_space::{
        assemble_elements, compare_galerkin, solve_axial_manufactured, solve_inherited,
        ExtractedElement,
    },
    inherited_space::build_space,
    json::parse,
};
use std::collections::BTreeMap;

fn axis_identity(degree: usize) -> solidcraft_trivariate::extracted_space::AxisExtraction {
    let n = degree + 1;
    solidcraft_trivariate::extracted_space::AxisExtraction {
        degree,
        basis_count: n,
        row_offsets: (0..=n).collect(),
        columns: (0..n).collect(),
        coefficients: vec![solidcraft_trivariate::constraint_system::Ratio::one(); n],
    }
}

fn q2_cube_element() -> (ExtractedElement, AnalysisSpace) {
    let axes = [axis_identity(2), axis_identity(2), axis_identity(2)];
    let mut controls = Vec::new();
    let mut dofs = Vec::new();
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                dofs.push((i * 3 + j) * 3 + k);
                controls.push([i as f64, j as f64, k as f64, 1.]);
            }
        }
    }
    let block = AnalysisBlock {
        id: "q2_unit_cube".into(),
        controls: controls.iter().map(|p| p.map(|x| x as i128)).collect(),
        dofs: dofs.clone(),
    };
    let element = ExtractedElement {
        axes,
        dofs,
        field_weights: vec![1.; 27],
        geometry_controls: controls,
    };
    let space = AnalysisSpace {
        ndof: 27,
        blocks: vec![block],
        boundary_dofs: BTreeMap::new(),
    };
    (element, space)
}

fn cylinder() -> solidcraft_trivariate::json::Json {
    parse(include_str!(
        "../docs/site/box-cylinder/tall/cylinder.model.json"
    ))
    .unwrap()
}
fn boundary_faces(
    model: &solidcraft_trivariate::json::Json,
    space: &solidcraft_trivariate::inherited_space::InheritedSpace,
) -> Vec<(usize, usize, bool)> {
    let mut faces = std::collections::BTreeSet::new();
    for face in model.field("design_faces").unwrap().items().unwrap() {
        let id = face.field("block").unwrap().text().unwrap();
        let side = face.field("side").unwrap().text().unwrap();
        let axis = match side.chars().next().unwrap() {
            'u' => 0,
            'v' => 1,
            _ => 2,
        };
        let upper = side.ends_with('1');
        let cell = space.cells.iter().position(|c| c.id == id).unwrap();
        faces.insert((cell, axis, upper));
    }
    faces.into_iter().collect()
}

#[test]
fn q2_adapter_reproduces_the_legacy_sparse_stiffness() {
    let (element, space) = q2_cube_element();
    let old = assemble_q2(&space, 210_000., 0.3, 4).unwrap();
    let new = assemble_elements(&[element], 27, 210_000., 0.3, 4)
        .unwrap()
        .stiffness;
    assert_eq!(old.offsets, new.offsets);
    assert_eq!(old.columns, new.columns);
    for (a, b) in old
        .values
        .iter()
        .flatten()
        .flatten()
        .zip(new.values.iter().flatten().flatten())
    {
        assert!((a - b).abs() <= 2e-11 * a.abs().max(1.));
    }
    assert!((old.integrated_volume - 8.).abs() < 1e-12);
}

#[test]
fn cylinder_torsion_matches_the_independent_204_vector_reference() {
    let model = cylinder();
    let space = build_space(&model, 1, 0, false, &[]).unwrap();
    assert_eq!(space.ndof * 3, 204);
    let exact_torque = 0.9817477042468105;
    let exact_energy = 0.017671458676442584;
    let job = parse(include_str!(
        "references/inherited_cylinder_torsion.job.json"
    ))
    .unwrap();
    let public = solve_inherited(&model, &job).unwrap();
    assert_eq!(
        public
            .field("vector_dofs")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        204
    );
    let public_torque = public
        .field("applied_torque_n_mm")
        .unwrap()
        .number::<f64>()
        .unwrap();
    assert!((public_torque / exact_torque - 1.).abs() < 2e-10);
    assert!(
        (public
            .field("strain_energy_n_mm")
            .unwrap()
            .number::<f64>()
            .unwrap()
            / exact_energy
            - 1.)
            .abs()
            < 2e-10
    );
    assert!(
        public
            .field("torque_balance_relative")
            .unwrap()
            .number::<f64>()
            .unwrap()
            < 2e-10
    );
    assert!(
        public
            .field("relative_residual")
            .unwrap()
            .number::<f64>()
            .unwrap()
            < 1e-9
    );
    assert!(
        public
            .field("l2_displacement_error")
            .unwrap()
            .number::<f64>()
            .unwrap()
            < 1e-9
    );
    assert!(
        public
            .field("torque_relative_error")
            .unwrap()
            .number::<f64>()
            .unwrap()
            < 2e-10
    );
    let mut invalid = job.clone();
    if let solidcraft_trivariate::json::Json::Object(fields) = &mut invalid {
        fields.insert(
            "quadrature_oder".into(),
            solidcraft_trivariate::json::Json::Number("8".into()),
        );
    }
    assert!(solve_inherited(&model, &invalid).is_err());
    let pocket = parse(include_str!(
        "../docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"
    ))
    .unwrap();
    assert!(solve_inherited(&pocket, &job).is_err());
}

#[test]
fn higher_axial_continuity_preserves_galerkin_operator_and_improves_manufactured_error() {
    let model = cylinder();
    for (p, c) in [(1, 0), (2, 0), (2, 1), (3, 0), (3, 2)] {
        let s = build_space(&model, p, c, false, &[]).unwrap();
        for cell in &s.cells {
            let (g, b) = cell.sample([0.231, 0.437, 0.619]).unwrap();
            let sum = b.values.iter().sum::<f64>();
            let position = std::array::from_fn::<_, 3, _>(|axis| {
                b.dofs
                    .iter()
                    .enumerate()
                    .map(|(i, d)| b.values[i] * s.control_xyz[*d][axis])
                    .sum::<f64>()
            });
            assert!(
                (sum - 1.).abs() < 1e-10,
                "partition p={p} c={c} {} {sum}",
                cell.id
            );
            for axis in 0..3 {
                assert!(
                    (position[axis] - g.position[axis]).abs() < 1e-9,
                    "affine p={p} c={c} {} {:?} {:?}",
                    cell.id,
                    position,
                    g.position
                );
            }
            for axis in 0..3 {
                for dir in 0..3 {
                    let grad = b
                        .dofs
                        .iter()
                        .enumerate()
                        .map(|(i, d)| b.physical_gradients[i][dir] * s.control_xyz[*d][axis])
                        .sum::<f64>();
                    let expected = if axis == dir { 1. } else { 0. };
                    assert!(
                        (grad - expected).abs() < 1e-8,
                        "gradient p={p} c={c} {} ({axis},{dir}) {grad}",
                        cell.id
                    );
                }
            }
        }
    }
    let c0 = build_space(&model, 2, 0, false, &[]).unwrap();
    let c1 = build_space(&model, 2, 1, false, &[]).unwrap();
    let gal = compare_galerkin(&c0, &c1, 2.6, 0.3, 8, |x| {
        [0., 0., (std::f64::consts::PI * x[2] / 36.).sin()]
    })
    .unwrap();
    assert!(gal.same_quadrature_points);
    assert!(gal.constrained_scalar_dofs < gal.c0_scalar_dofs);
    assert!(gal.stiffness_relative_frobenius_defect < 3e-11, "{gal:?}");
    assert!(gal.force_relative_l2_defect < 3e-11, "{gal:?}");
    eprintln!("C0→C1 Galerkin {gal:?}");

    let mut errors = Vec::new();
    let mut l2s = Vec::new();
    let mut dofs = Vec::new();
    for (p, c) in [(1, 0), (2, 0), (2, 1), (3, 0), (3, 2)] {
        let space = build_space(&model, p, c, false, &[]).unwrap();
        let bounds = [-12., 24.];
        let result = solve_axial_manufactured(
            space.ndof,
            &space.control_xyz,
            bounds,
            space.cells.len(),
            210_000.,
            0.3,
            10,
            1e-10,
            4000,
            &boundary_faces(&model, &space),
            |cell, uvw| space.cells[cell].sample(uvw),
        )
        .unwrap();
        errors.push(result.energy_error_sqrt_n_mm);
        l2s.push(result.l2_displacement_error_mm_sqrt_mm3);
        dofs.push(space.ndof);
        assert!(result.relative_residual < 1e-8);
    }
    eprintln!("manufactured axial degree/continuity dofs={dofs:?}, energy={errors:?}, l2={l2s:?}");
    assert!(
        dofs[0] < dofs[2] && dofs[2] < dofs[4],
        "degree/continuity sequence scalar DOFs={dofs:?}"
    );
    assert!(
        errors[2] < errors[0] && errors[4] < errors[2],
        "degree/continuity sequence energy={errors:?}, l2={l2s:?}"
    );
    assert!(
        errors[1] < errors[0] && errors[3] < errors[1],
        "uniform-degree sequence energy={errors:?}, l2={l2s:?}"
    );
}

#[test]
fn blind_pocket_mms_reports_error_dofs_and_separate_costs() {
    let model = parse(include_str!(
        "../docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"
    ))
    .unwrap();
    let base = parse(include_str!(
        "references/inherited_blind_pocket_mms.job.json"
    ))
    .unwrap();
    let mut rows = Vec::new();
    for (degree, continuity) in [(1, 0), (2, 1), (3, 2)] {
        let mut job = base.clone();
        if let solidcraft_trivariate::json::Json::Object(fields) = &mut job {
            fields.insert(
                "degree".into(),
                solidcraft_trivariate::json::Json::Number(degree.to_string()),
            );
            fields.insert(
                "continuity".into(),
                solidcraft_trivariate::json::Json::Number(continuity.to_string()),
            );
        }
        let result = solve_inherited(&model, &job).unwrap();
        let scalar_dofs = result
            .field("scalar_dofs")
            .unwrap()
            .number::<usize>()
            .unwrap();
        let energy = result
            .field("energy_error")
            .unwrap()
            .number::<f64>()
            .unwrap();
        let l2 = result
            .field("l2_displacement_error")
            .unwrap()
            .number::<f64>()
            .unwrap();
        let assembly = result
            .field("assembly_ms")
            .unwrap()
            .number::<f64>()
            .unwrap();
        let solve = result.field("solve_ms").unwrap().number::<f64>().unwrap();
        let error_eval = result
            .field("error_evaluation_ms")
            .unwrap()
            .number::<f64>()
            .unwrap();
        let nnz = result
            .field("matrix_nonzeros")
            .unwrap()
            .number::<usize>()
            .unwrap();
        rows.push((
            degree,
            continuity,
            scalar_dofs,
            nnz,
            energy,
            l2,
            assembly,
            solve,
            error_eval,
        ));
    }
    eprintln!("blind pocket manufactured degree/continuity, scalar dofs, nnz, energy, L2, assembly/solve/error ms: {rows:?}");
    assert!(rows.windows(2).all(|pair| pair[1].2 > pair[0].2));
    assert!(
        rows.windows(2).all(|pair| pair[1].4 < pair[0].4),
        "pocket energy errors: {rows:?}"
    );
    assert!(
        rows.windows(2).all(|pair| pair[1].5 < pair[0].5),
        "pocket displacement errors: {rows:?}"
    );
    assert!(rows
        .iter()
        .all(|r| r.6.is_finite() && r.7.is_finite() && r.8.is_finite() && r.3 > 0));
}
