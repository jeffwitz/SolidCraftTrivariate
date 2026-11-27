use solidcraft_trivariate::{
    constraint_system::Ratio,
    json::{parse, Json},
    thb_strip::{build_thb_strip_from_recipe, build_thb_strip_from_recipe_cells},
};

fn recipe() -> solidcraft_trivariate::json::Json {
    parse(include_str!("references/nurbs_strip.recipe.json")).unwrap()
}

#[test]
fn exact_interior_thb_and_prepared_share_the_marked_partition() {
    let strip = build_thb_strip_from_recipe(&recipe()).unwrap();
    assert_eq!(
        (
            strip.source_parent_ndof,
            strip.initial_ndof,
            strip.prepared_ndof,
            strip.uniform_fine_ndof
        ),
        (20, 80, 539, 2299)
    );
    assert_eq!(
        (
            strip.coarse_active_count,
            strip.fine_active_count,
            strip.ndof
        ),
        (539, 5, 544)
    );
    assert_eq!(strip.transfer_rank, 539);
    assert_eq!(strip.uniform_fine_rank, 544);
    assert_eq!(strip.leaves.len(), 240);
    assert!(strip.common_partition_verified);
    assert_eq!(strip.prepared_shape, [11, 7, 7]);
    assert_eq!(
        strip.prepared_zero_bottom_lift().unwrap().free_roots.len(),
        462
    );
    assert_eq!(
        strip.leaves.iter().filter(|leaf| leaf.level == 0).count(),
        112
    );
    assert_eq!(
        strip.leaves.iter().filter(|leaf| leaf.level == 1).count(),
        128
    );
    let prepared_coefficients: Vec<_> = (0..strip.prepared_ndof)
        .map(|i| Ratio::new((i % 13) as i128, 7).unwrap())
        .collect();
    let thb_coefficients = strip.transfer_field(&prepared_coefficients).unwrap();
    for (i, _leaf) in strip.leaves.iter().enumerate() {
        let point = [0.371, 0.413, 0.619];
        let (geometry, prepared) = strip.sample_prepared(i, point).unwrap();
        let (geometry_thb, local) = strip.sample(i, point).unwrap();
        assert_eq!(geometry.position, geometry_thb.position);
        assert_eq!(geometry.jacobian_columns, geometry_thb.jacobian_columns);
        assert!((prepared.values.iter().sum::<f64>() - 1.).abs() < 2e-13);
        assert!((local.values.iter().sum::<f64>() - 1.).abs() < 2e-13);
        let prepared_value: f64 = prepared
            .dofs
            .iter()
            .zip(&prepared.values)
            .map(|(&dof, &value)| {
                value * prepared_coefficients[dof].numerator() as f64
                    / prepared_coefficients[dof].denominator() as f64
            })
            .sum();
        let local_value: f64 = local
            .dofs
            .iter()
            .zip(&local.values)
            .map(|(&dof, &value)| {
                value * thb_coefficients[dof].numerator() as f64
                    / thb_coefficients[dof].denominator() as f64
            })
            .sum();
        assert!((prepared_value - local_value).abs() < 2e-12);
        for axis in 0..3 {
            let prepared_gradient: f64 = prepared
                .dofs
                .iter()
                .zip(&prepared.physical_gradients)
                .map(|(&dof, gradient)| {
                    gradient[axis] * prepared_coefficients[dof].numerator() as f64
                        / prepared_coefficients[dof].denominator() as f64
                })
                .sum();
            let local_gradient: f64 = local
                .dofs
                .iter()
                .zip(&local.physical_gradients)
                .map(|(&dof, gradient)| {
                    gradient[axis] * thb_coefficients[dof].numerator() as f64
                        / thb_coefficients[dof].denominator() as f64
                })
                .sum();
            assert!((prepared_gradient - local_gradient).abs() < 2e-10);
        }
    }
}

#[test]
fn alternative_interior_box_is_supported_and_invalid_boxes_are_refused() {
    let recipe = recipe();
    let other = build_thb_strip_from_recipe_cells(&recipe, [1, 1, 1], [7, 3, 3]).unwrap();
    assert!(other.common_partition_verified);
    assert!(other.leaves.len() < 4096);
    assert!(build_thb_strip_from_recipe_cells(&recipe, [0, 1, 1], [2, 3, 3]).is_err());
    assert!(build_thb_strip_from_recipe_cells(&recipe, [1, 2, 1], [3, 3, 3]).is_err());
    assert!(build_thb_strip_from_recipe_cells(&recipe, [2, 1, 1], [6, 3, 8]).is_err());
}

#[test]
fn reversed_source_orientation_is_refused() {
    let mut recipe = recipe();
    let Json::Object(root) = &mut recipe else {
        panic!("fixture recipe must be an object")
    };
    for key in ["lower", "upper"] {
        let Some(Json::Object(curve)) = root.get_mut(key) else {
            panic!("fixture curve missing")
        };
        let Some(Json::Array(controls)) = curve.get_mut("controls") else {
            panic!("fixture controls missing")
        };
        for control in controls {
            let Json::Array(values) = control else {
                panic!("fixture control must be an array")
            };
            let x = values[0].number::<i64>().unwrap();
            let w = values[2].number::<i64>().unwrap();
            values[0] = Json::Number((100 * w - x).to_string());
        }
    }
    assert!(build_thb_strip_from_recipe(&recipe).is_err());
}

#[test]
fn nonuniform_source_mesh_is_refused_by_shared_admission() {
    let mut recipe = recipe();
    let Json::Object(root) = &mut recipe else {
        panic!("fixture recipe must be an object")
    };
    for key in ["lower", "upper"] {
        let Some(Json::Object(curve)) = root.get_mut(key) else {
            panic!("fixture curve missing")
        };
        let Some(Json::Array(knots)) = curve.get_mut("knots") else {
            panic!("fixture knots missing")
        };
        knots[4] = Json::Number("0".into());
    }
    assert!(build_thb_strip_from_recipe(&recipe).is_err());
}
