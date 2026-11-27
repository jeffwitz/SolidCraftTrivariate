use solidcraft_trivariate::{
    bracket::BracketBlock,
    constraint_system::Ratio,
    incidence::{BlockSide, Orientation},
    json::{object, Json},
    pair_parent::{recover_model_interfaces, recover_pair},
};

fn side(block: &BracketBlock, axis: usize, upper: bool) -> BlockSide {
    BlockSide {
        block: block.id.clone(),
        axis,
        upper,
    }
}

fn make_piece(id: &str, degrees: [usize; 3], join_axis: usize, values: &[i64]) -> BracketBlock {
    let mut controls = Vec::new();
    for i in 0..=degrees[0] {
        for j in 0..=degrees[1] {
            for k in 0..=degrees[2] {
                let p = [i, j, k];
                let normal = values[p[join_axis]];
                let x = if join_axis == 0 {
                    normal
                } else {
                    (i as i64) * 10
                };
                let y = if join_axis == 1 {
                    normal
                } else {
                    (j as i64) * 10
                };
                let z = if join_axis == 2 {
                    normal
                } else {
                    (k as i64) * 10
                };
                controls.push([x, y, z, 1]);
            }
        }
    }
    BracketBlock {
        id: id.into(),
        degrees,
        controls,
    }
}

fn scale(block: &BracketBlock, factor: i64, id: &str) -> BracketBlock {
    BracketBlock {
        id: id.into(),
        degrees: block.degrees,
        controls: block
            .controls
            .iter()
            .map(|p| p.map(|x| x * factor))
            .collect(),
    }
}

fn oriented_second(
    source: &BracketBlock,
    first_axis: usize,
    second_axis: usize,
    second_upper: bool,
    orientation: Orientation,
    id: &str,
) -> BracketBlock {
    let first_free: Vec<_> = (0..3).filter(|&a| a != first_axis).collect();
    let second_free: Vec<_> = (0..3).filter(|&a| a != second_axis).collect();
    let mut input_degrees = [0; 3];
    input_degrees[second_axis] = source.degrees[first_axis];
    for (out_tangent, &parent_axis) in first_free.iter().enumerate() {
        let in_tangent = if orientation.swap {
            1 - out_tangent
        } else {
            out_tangent
        };
        input_degrees[second_free[in_tangent]] = source.degrees[parent_axis];
    }
    let mut controls = vec![[0; 4]; source.controls.len()];
    for i in 0..=source.degrees[0] {
        for j in 0..=source.degrees[1] {
            for k in 0..=source.degrees[2] {
                let output = [i, j, k];
                let mut input = [0; 3];
                input[second_axis] = if second_upper {
                    source.degrees[first_axis] - output[first_axis]
                } else {
                    output[first_axis]
                };
                for (out_tangent, &parent_axis) in first_free.iter().enumerate() {
                    let in_tangent = if orientation.swap {
                        1 - out_tangent
                    } else {
                        out_tangent
                    };
                    let flip = if orientation.swap {
                        if out_tangent == 0 {
                            orientation.flip_second
                        } else {
                            orientation.flip_first
                        }
                    } else if out_tangent == 0 {
                        orientation.flip_first
                    } else {
                        orientation.flip_second
                    };
                    let source_axis = second_free[in_tangent];
                    input[source_axis] = if flip {
                        input_degrees[source_axis] - output[parent_axis]
                    } else {
                        output[parent_axis]
                    };
                }
                let src = (i * (source.degrees[1] + 1) + j) * (source.degrees[2] + 1) + k;
                let dst = (input[0] * (input_degrees[1] + 1) + input[1]) * (input_degrees[2] + 1)
                    + input[2];
                controls[dst] = source.controls[src];
            }
        }
    }
    BracketBlock {
        id: id.into(),
        degrees: input_degrees,
        controls,
    }
}

fn elevate_y_to_quadratic(block: &BracketBlock, id: &str) -> BracketBlock {
    assert_eq!(block.degrees[1], 1);
    let mut controls = Vec::new();
    for i in 0..=block.degrees[0] {
        for j in 0..=2 {
            for k in 0..=block.degrees[2] {
                let left = (i * 2) * (block.degrees[2] + 1) + k;
                let right = (i * 2 + 1) * (block.degrees[2] + 1) + k;
                let p0 = block.controls[left];
                let p1 = block.controls[right];
                let point = if j == 0 {
                    p0
                } else if j == 2 {
                    p1
                } else {
                    std::array::from_fn(|c| (p0[c] + p1[c]) / 2)
                };
                controls.push(point);
            }
        }
    }
    BracketBlock {
        id: id.into(),
        degrees: [block.degrees[0], 2, block.degrees[2]],
        controls,
    }
}

fn identity() -> Orientation {
    Orientation {
        swap: false,
        flip_first: false,
        flip_second: false,
    }
}

fn block_json(block: &BracketBlock) -> Json {
    object([
        ("id", Json::String(block.id.clone())),
        (
            "degrees",
            Json::Array(
                block
                    .degrees
                    .iter()
                    .map(|x| Json::Number(x.to_string()))
                    .collect(),
            ),
        ),
        (
            "controls",
            Json::Array(
                block
                    .controls
                    .iter()
                    .map(|p| Json::Array(p.iter().map(|x| Json::Number(x.to_string())).collect()))
                    .collect(),
            ),
        ),
    ])
}

#[test]
fn recovers_unequal_split_and_positive_projective_gauge() {
    let first = make_piece("a", [2, 1, 2], 0, &[0, 10, 20]);
    let second_base = make_piece("b0", [2, 1, 2], 0, &[20, 40, 60]);
    let second = scale(&second_base, 2, "b");
    let recovered = recover_pair(
        &first,
        side(&first, 0, true),
        &second,
        side(&second, 0, false),
        identity(),
    )
    .unwrap();
    assert_eq!(recovered.split_parameter, Ratio::new(1, 3).unwrap());
    assert_eq!(
        recovered.first_to_second_interval_ratio,
        Ratio::new(1, 2).unwrap()
    );
    assert_eq!(recovered.second_to_first_gauge, Ratio::new(1, 2).unwrap());
    assert_eq!(recovered.degrees, [2, 1, 2]);
    assert_eq!(recovered.controls[0][3], Ratio::one());
    assert!(recovered
        .document()
        .to_string()
        .contains("mathematical_parent_recovery_not_feature_provenance"));
}

#[test]
fn recovers_a_swapped_and_flipped_face_on_different_axes() {
    let degrees = [2, 2, 2];
    let first = make_piece("a", degrees, 1, &[0, 10, 20]);
    let canonical_second = make_piece("b0", degrees, 1, &[20, 40, 60]);
    let orientation = Orientation {
        swap: true,
        flip_first: true,
        flip_second: false,
    };
    let second = oriented_second(&canonical_second, 1, 2, false, orientation, "b");
    let recovered = recover_pair(
        &first,
        side(&first, 1, true),
        &second,
        side(&second, 2, false),
        orientation,
    )
    .unwrap();
    assert_eq!(recovered.split_parameter, Ratio::new(1, 3).unwrap());
    assert_eq!(recovered.degrees, degrees);
    assert_eq!(
        recovered.first_to_second_interval_ratio,
        Ratio::new(1, 2).unwrap()
    );
}

#[test]
fn reverses_both_input_normal_directions_and_recovers_inverse_ratio() {
    let right = make_piece("right", [2, 1, 2], 0, &[20, 40, 60]);
    let left = make_piece("left", [2, 1, 2], 0, &[0, 10, 20]);
    let recovered = recover_pair(
        &right,
        side(&right, 0, false),
        &left,
        side(&left, 0, true),
        identity(),
    )
    .unwrap();
    assert_eq!(recovered.split_parameter, Ratio::new(2, 3).unwrap());
    assert_eq!(
        recovered.first_to_second_interval_ratio,
        Ratio::new(2, 1).unwrap()
    );
}

#[test]
fn elevates_tangential_degree_before_comparing_full_traces() {
    let first = make_piece("a", [2, 1, 2], 0, &[0, 10, 20]);
    let second_linear = make_piece("b0", [2, 1, 2], 0, &[20, 40, 60]);
    let second = elevate_y_to_quadratic(&second_linear, "b");
    let recovered = recover_pair(
        &first,
        side(&first, 0, true),
        &second,
        side(&second, 0, false),
        identity(),
    )
    .unwrap();
    assert_eq!(recovered.degrees, [2, 2, 2]);
    assert_eq!(recovered.split_parameter, Ratio::new(1, 3).unwrap());
}

#[test]
fn rejects_a_pair_that_matches_the_first_jet_but_not_the_full_restriction() {
    let first = make_piece("a", [3, 1, 1], 0, &[0, 10, 20, 30]);
    let mut second = make_piece("b", [3, 1, 1], 0, &[30, 50, 70, 90]);
    // This changes a higher normal Bernstein coefficient but leaves the face
    // and first normal jet unchanged.
    second.controls[2 * 4][0] += 1;
    let refusal = recover_pair(
        &first,
        side(&first, 0, true),
        &second,
        side(&second, 0, false),
        identity(),
    )
    .unwrap_err();
    assert_eq!(
        refusal,
        "full restrictions do not reproduce both supplied homogeneous nets"
    );
}

#[test]
fn refuses_wrong_face_and_nonpositive_homogeneous_weights() {
    let first = make_piece("a", [2, 1, 2], 0, &[0, 10, 20]);
    let mut second = make_piece("b", [2, 1, 2], 0, &[21, 41, 61]);
    assert!(recover_pair(
        &first,
        side(&first, 0, true),
        &second,
        side(&second, 0, false),
        identity()
    )
    .is_err());
    second.controls[0][3] = 0;
    assert_eq!(
        recover_pair(
            &first,
            side(&first, 0, true),
            &second,
            side(&second, 0, false),
            identity()
        )
        .unwrap_err(),
        "pair parent requires positive homogeneous weights"
    );
}

#[test]
fn refuses_opposing_nonzero_first_jets_instead_of_marking_them_indeterminate() {
    let first = make_piece("a", [2, 1, 2], 0, &[0, 10, 20]);
    let second = make_piece("b", [2, 1, 2], 0, &[20, 10, 0]);
    assert_eq!(
        recover_pair(
            &first,
            side(&first, 0, true),
            &second,
            side(&second, 0, false),
            identity(),
        )
        .unwrap_err(),
        "pair parent first jets imply no positive interval ratio"
    );
}

#[test]
fn model_report_preserves_interface_order_and_embeds_parent_evidence() {
    let first = make_piece("a", [2, 1, 2], 0, &[0, 10, 20]);
    let second = make_piece("b", [2, 1, 2], 0, &[20, 40, 60]);
    let interface = object([
        (
            "first",
            Json::Array(vec![Json::String("a".into()), Json::String("u1".into())]),
        ),
        (
            "second",
            Json::Array(vec![Json::String("b".into()), Json::String("u0".into())]),
        ),
        (
            "orientation",
            object([
                ("swap", Json::Bool(false)),
                ("flip_first", Json::Bool(false)),
                ("flip_second", Json::Bool(false)),
            ]),
        ),
    ]);
    let topology = object([("interfaces", Json::Array(vec![interface]))]);
    let model = object([
        ("schema", Json::String("pair_parent_test_v1".into())),
        (
            "blocks",
            Json::Array(vec![block_json(&first), block_json(&second)]),
        ),
        ("topology", topology),
    ]);
    let report = recover_model_interfaces(&model).unwrap();
    assert_eq!(
        report
            .field("interface_count")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        1
    );
    assert_eq!(
        report
            .field("recovered_count")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        1
    );
    let row = &report.field("interfaces").unwrap().items().unwrap()[0];
    assert_eq!(
        row.field("interface_index")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        0
    );
    assert_eq!(row.field("status").unwrap().text().unwrap(), "recovered");
    assert_eq!(
        row.field("recovery")
            .unwrap()
            .field("schema")
            .unwrap()
            .text()
            .unwrap(),
        "pair_bezier_parent_v1"
    );
}
