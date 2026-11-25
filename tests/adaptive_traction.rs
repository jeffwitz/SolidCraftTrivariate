use solidcraft_trivariate::{
    adaptive_traction::{integrate_adaptive_traction, SurfaceQuadrature},
    analysis_space::{AnalysisBlock, AnalysisSpace},
    boundary_traction::{integrate_traction, LoadPatch},
    elasticity::Vector,
    incidence::BlockSide,
};
use std::collections::{BTreeMap, BTreeSet};

fn space_from_controls(blocks: Vec<AnalysisBlock>) -> AnalysisSpace {
    let ndof = blocks
        .iter()
        .flat_map(|b| b.dofs.iter().copied())
        .max()
        .unwrap_or(0)
        + 1;
    AnalysisSpace {
        ndof,
        blocks,
        boundary_dofs: BTreeMap::new(),
    }
}

fn cube_space() -> AnalysisSpace {
    let mut controls = Vec::with_capacity(27);
    for i in 0..3_i128 {
        for j in 0..3_i128 {
            for k in 0..3_i128 {
                controls.push([i, j, k, 2]);
            }
        }
    }
    space_from_controls(vec![AnalysisBlock {
        id: "cube".into(),
        controls,
        dofs: (0..27).collect(),
    }])
}

fn circle_ring_space() -> AnalysisSpace {
    let radii = [4_i128, 3, 2]; // Gauge-two controls for physical radii 2, 1.5, 1.
    let weights = [2_i128, 2, 4];
    let mut blocks = Vec::new();
    for quadrant in 0..4 {
        let mut controls = Vec::with_capacity(27);
        for (i, &weight) in weights.iter().enumerate() {
            for radius_twice in radii {
                let [x, y] = match i {
                    0 => [radius_twice, 0],
                    1 => [radius_twice, radius_twice],
                    _ => [0, 2 * radius_twice],
                };
                let [x, y] = match quadrant {
                    0 => [x, y],
                    1 => [-y, x],
                    2 => [-x, -y],
                    _ => [y, -x],
                };
                for k in 0..3_i128 {
                    controls.push([x, y, k * weight, weight]);
                }
            }
        }
        let index = quadrant * 27;
        blocks.push(AnalysisBlock {
            id: format!("arc{quadrant}"),
            controls,
            dofs: (index..index + 27).collect(),
        });
    }
    space_from_controls(blocks)
}

fn skew_space() -> AnalysisSpace {
    let mut controls = Vec::with_capacity(27);
    for i in 0..3_i128 {
        for j in 0..3_i128 {
            for k in 0..3_i128 {
                controls.push([10 * i + 15 * j, 2 * j, 4 * i + 8 * k, 2]);
            }
        }
    }
    space_from_controls(vec![AnalysisBlock {
        id: "skew".into(),
        controls,
        dofs: (0..27).collect(),
    }])
}

fn side(block: &str, axis: usize, upper: bool) -> BlockSide {
    BlockSide {
        block: block.into(),
        axis,
        upper,
    }
}

fn config(tolerance: f64, depth: u8, panels: usize) -> SurfaceQuadrature {
    SurfaceQuadrature {
        relative_tolerance: tolerance,
        max_depth: depth,
        max_panels: panels,
    }
}

fn error<T>(result: Result<T, &'static str>) -> &'static str {
    match result {
        Err(reason) => reason,
        Ok(_) => panic!("expected adaptive surface integration to be refused"),
    }
}

#[test]
fn plane_and_extreme_oblique_faces_integrate_exactly_and_respect_budget() {
    let cube = cube_space();
    let top = side("cube", 2, true);
    let result = integrate_adaptive_traction(
        &cube,
        &BTreeSet::from([top.clone()]),
        &[LoadPatch::whole(top.clone())],
        3,
        config(1e-10, 3, 100),
        |_, _| [0., 0., 3.],
    )
    .unwrap();
    assert!((result.load.area_mm2 - 1.).abs() < 1e-13);
    assert!((result.load.resultant_n[2] - 3.0).abs() < 1e-14);
    assert_eq!(result.panels, 5);
    assert_eq!(result.deepest, 1);
    assert!(result.maximum_estimated_relative_error <= 1e-10);

    let skew = skew_space();
    let oblique = side("skew", 0, true);
    let pressure = integrate_adaptive_traction(
        &skew,
        &BTreeSet::from([oblique.clone()]),
        &[LoadPatch::whole(oblique.clone())],
        3,
        config(1e-8, 2, 100),
        |_, n| n.map(|x| 2.5 * x),
    )
    .unwrap();
    let exact_area = (16_f64.powi(2) + 120_f64.powi(2)).sqrt();
    assert!((pressure.load.area_mm2 - exact_area).abs() / exact_area < 1e-13);
    for (actual, normal) in
        pressure
            .load
            .resultant_n
            .iter()
            .zip([16. / exact_area, -120. / exact_area, 0.])
    {
        assert!((actual - 2.5 * exact_area * normal).abs() < 1e-11);
    }

    let too_small = integrate_adaptive_traction(
        &cube,
        &BTreeSet::from([top.clone()]),
        &[LoadPatch::whole(top)],
        2,
        config(1e-8, 2, 4),
        |_, _| [0., 0., 1.],
    );
    assert_eq!(
        error(too_small),
        "adaptive surface quadrature panel budget exceeded"
    );
}

#[test]
fn circular_pressure_refines_and_handles_zero_global_resultant_scale() {
    let space = circle_ring_space();
    let sides = (0..4)
        .map(|q| side(&format!("arc{q}"), 1, true))
        .collect::<Vec<_>>();
    let exterior = sides.iter().cloned().collect::<BTreeSet<_>>();
    let patches = sides
        .iter()
        .cloned()
        .map(LoadPatch::whole)
        .collect::<Vec<_>>();
    let adaptive = integrate_adaptive_traction(
        &space,
        &exterior,
        &patches,
        8,
        config(1e-8, 8, 20000),
        |_, normal| normal.map(|x| 3.0 * x),
    )
    .unwrap();
    assert!(adaptive.panels > 4);
    assert!(adaptive.deepest > 1);
    assert!(adaptive.maximum_estimated_relative_error <= 1e-8);
    let exact_area = 4.0 * std::f64::consts::PI;
    assert!((adaptive.load.area_mm2 - exact_area).abs() / exact_area < 1e-8);
    assert!(adaptive.load.resultant_n.iter().all(|x| x.abs() < 1e-9));
    assert!(adaptive.load.moment_n_mm.iter().all(|x| x.abs() < 1e-9));

    let reference_patches = sides
        .iter()
        .flat_map(|face| {
            (0..8).flat_map(move |i| {
                (0..8).map(move |j| LoadPatch {
                    side: face.clone(),
                    depth: 3,
                    intervals: [[i, i + 1], [j, j + 1]],
                })
            })
        })
        .collect::<Vec<_>>();
    let reference = integrate_traction(&space, &exterior, &reference_patches, 16, |_, normal| {
        normal.map(|x| 3.0 * x)
    })
    .unwrap();
    assert!((adaptive.load.area_mm2 - reference.area_mm2).abs() / reference.area_mm2 < 1e-8);
    let error = adaptive
        .load
        .nodal_force_n
        .iter()
        .zip(&reference.nodal_force_n)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt();
    let scale = reference
        .nodal_force_n
        .iter()
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt();
    assert!(error / scale < 1e-8);
    let mut affine_u = vec![0.0; 3 * space.ndof];
    for block in &space.blocks {
        for (&dof, h) in block.dofs.iter().zip(&block.controls) {
            let point = [
                h[0] as f64 / h[3] as f64,
                h[1] as f64 / h[3] as f64,
                h[2] as f64 / h[3] as f64,
            ];
            affine_u[3 * dof] = 0.1 * point[0];
            affine_u[3 * dof + 1] = -0.05 * point[1];
            affine_u[3 * dof + 2] = 0.03 * point[2];
        }
    }
    let adaptive_work = adaptive
        .load
        .nodal_force_n
        .iter()
        .zip(&affine_u)
        .map(|(f, u)| f * u)
        .sum::<f64>();
    let reference_work = reference
        .nodal_force_n
        .iter()
        .zip(&affine_u)
        .map(|(f, u)| f * u)
        .sum::<f64>();
    assert!((adaptive_work - reference_work).abs() / reference_work.abs().max(1.0) < 1e-8);

    let uniform_traction = integrate_adaptive_traction(
        &space,
        &exterior,
        &patches,
        8,
        config(1e-8, 8, 20000),
        |_, _| [0.5, -0.25, 1.0],
    )
    .unwrap();
    let traction_reference =
        integrate_traction(&space, &exterior, &reference_patches, 16, |_, _| {
            [0.5, -0.25, 1.0]
        })
        .unwrap();
    let traction_diff = uniform_traction
        .load
        .nodal_force_n
        .iter()
        .zip(&traction_reference.nodal_force_n)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt();
    let traction_scale = traction_reference
        .nodal_force_n
        .iter()
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt();
    assert!(traction_diff / traction_scale < 1e-8);

    let zero = integrate_adaptive_traction(
        &space,
        &exterior,
        &patches,
        8,
        config(1e-8, 2, 100),
        |_, _| [0., 0., 0.],
    )
    .unwrap();
    assert!(zero.maximum_estimated_relative_error <= 1e-8);
    assert_eq!(zero.load.resultant_n, [0.; 3]);
}

#[test]
fn configuration_depth_selection_and_face_overlap_are_explicitly_refused() {
    assert!(config(f64::NAN, 4, 100).validate().is_err());
    assert!(config(1e-8, 0, 100).validate().is_err());
    let space = cube_space();
    let top = side("cube", 2, true);
    let exterior = BTreeSet::from([top.clone()]);
    let oscillatory = |x: Vector, _| [0., 0., (100. * x[0]).sin()];
    let depth_failure = integrate_adaptive_traction(
        &space,
        &exterior,
        &[LoadPatch::whole(top.clone())],
        1,
        config(1e-14, 1, 100),
        oscillatory,
    );
    assert_eq!(
        error(depth_failure),
        "adaptive surface quadrature tolerance not met within depth budget"
    );

    let overlap = [LoadPatch::whole(top.clone()), LoadPatch::whole(top.clone())];
    let overlap_error = integrate_adaptive_traction(
        &space,
        &exterior,
        &overlap,
        2,
        config(1e-8, 2, 100),
        |_, _| [1., 0., 0.],
    );
    assert_eq!(error(overlap_error), "overlapping traction patches");
    let wrong = side("cube", 2, false);
    let bad_side = integrate_adaptive_traction(
        &space,
        &exterior,
        &[LoadPatch::whole(wrong)],
        2,
        config(1e-8, 2, 100),
        |_, _| [1., 0., 0.],
    );
    assert_eq!(
        error(bad_side),
        "traction requires a certified exterior side"
    );
}
