use solidcraft_trivariate::bracket::BracketBlock;
use solidcraft_trivariate::host_composition::{
    compose_with_host, AffinePiece, HostCompositionError, HostNodeStatus, HostRepresentation,
    PiecewiseAffineHost, RationalParameter,
};

fn quadratic_block(id: &str, x_numerators: [i64; 3], weights: [i64; 3]) -> BracketBlock {
    let mut controls = Vec::with_capacity(12);
    for i in 0..3 {
        for j in 0..2 {
            for k in 0..2 {
                let weight = weights[i];
                controls.push([
                    x_numerators[i],
                    j as i64 * weight,
                    k as i64 * weight,
                    weight,
                ]);
            }
        }
    }
    BracketBlock {
        id: id.into(),
        degrees: [2, 1, 1],
        controls,
    }
}

fn linear_block(id: &str) -> BracketBlock {
    let mut controls = Vec::with_capacity(8);
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                controls.push([i as i64, j as i64, k as i64, 1]);
            }
        }
    }
    BracketBlock {
        id: id.into(),
        degrees: [1, 1, 1],
        controls,
    }
}

fn piece(x_slope: i64, x_intercept: i64) -> AffinePiece {
    AffinePiece {
        matrix: [[x_slope, 0, 0], [0, 2, 0], [0, 0, 3]],
        translation: [x_intercept, 1, -2],
    }
}

fn host_at_two() -> PiecewiseAffineHost {
    PiecewiseAffineHost {
        breaks: vec![2],
        pieces: vec![piece(2, 1), piece(3, -1)],
    }
}

fn bernstein(degree: usize, x: f64) -> Vec<f64> {
    match degree {
        1 => vec![1.0 - x, x],
        2 => vec![(1.0 - x).powi(2), 2.0 * x * (1.0 - x), x * x],
        _ => unreachable!(),
    }
}

fn bernstein_derivative(degree: usize, x: f64) -> Vec<f64> {
    match degree {
        1 => vec![-1.0, 1.0],
        2 => vec![-2.0 * (1.0 - x), 2.0 - 4.0 * x, 2.0 * x],
        _ => unreachable!(),
    }
}

/// Independent tensor-product rational evaluation and all three map columns.
fn evaluate(block: &BracketBlock, q: [f64; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let bases: [Vec<f64>; 3] = std::array::from_fn(|axis| bernstein(block.degrees[axis], q[axis]));
    let derivatives: [Vec<f64>; 3] =
        std::array::from_fn(|axis| bernstein_derivative(block.degrees[axis], q[axis]));
    let mut h = [0.0; 4];
    let mut dh = [[0.0; 4]; 3];
    let mut index = 0;
    for i in 0..=block.degrees[0] {
        for j in 0..=block.degrees[1] {
            for k in 0..=block.degrees[2] {
                let control = block.controls[index];
                index += 1;
                let factors = [bases[0][i], bases[1][j], bases[2][k]];
                let ids = [i, j, k];
                let factor = factors[0] * factors[1] * factors[2];
                for c in 0..4 {
                    h[c] += factor * control[c] as f64;
                }
                for axis in 0..3 {
                    let mut derivative_factor = derivatives[axis][ids[axis]];
                    for (other, component) in factors.iter().enumerate() {
                        if other != axis {
                            derivative_factor *= component;
                        }
                    }
                    for c in 0..4 {
                        dh[axis][c] += derivative_factor * control[c] as f64;
                    }
                }
            }
        }
    }
    let point = [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
    let jacobian = std::array::from_fn(|axis| {
        std::array::from_fn(|c| (dh[axis][c] - point[c] * dh[axis][3]) / h[3])
    });
    (point, jacobian)
}

fn apply_piece(piece: &AffinePiece, point: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|row| {
        piece.matrix[row]
            .iter()
            .zip(point)
            .map(|(&a, x)| a as f64 * x)
            .sum::<f64>()
            + piece.translation[row] as f64
    })
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn exact_ratio(parameter: RationalParameter, numerator: i128, denominator: i128) {
    assert_eq!(parameter.numerator, numerator);
    assert_eq!(parameter.denominator, denominator);
}

#[test]
fn affine_piecewise_composition_matches_direct_map_and_one_sided_jets() {
    // X-2W has Bernstein coefficients [-1,0,1], hence its unique root is
    // exactly u=1/2. W=[1,2,1] is nonuniform while X/W=[1,2,3] is monotone.
    let source = quadratic_block("rational_root", [1, 4, 3], [1, 2, 1]);
    let host = host_at_two();
    let composition = compose_with_host(&source, &host).unwrap();
    assert_eq!(
        composition.representation,
        HostRepresentation::SplitAtRationalNodes
    );
    assert_eq!(composition.blocks.len(), 2);
    assert_eq!(composition.node_events.len(), 1);
    assert_eq!(composition.node_events[0].plane_x, 2);
    match composition.node_events[0].status {
        HostNodeStatus::InternalRational(root) => exact_ratio(root, 1, 2),
        status => panic!("expected rational interior event, got {status:?}"),
    }

    let first = &composition.blocks[0];
    let second = &composition.blocks[1];
    assert_eq!(first.host_piece_index, 0);
    assert_eq!(second.host_piece_index, 1);
    exact_ratio(first.parameter_interval[0], 0, 1);
    exact_ratio(first.parameter_interval[1], 1, 2);
    exact_ratio(second.parameter_interval[0], 1, 2);
    exact_ratio(second.parameter_interval[1], 1, 1);

    let samples: [[f64; 3]; 3] = [[0.13, 0.37, 0.61], [0.5, 0.5, 0.5], [0.87, 0.19, 0.73]];
    for (chart, affine) in [(first, &host.pieces[0]), (second, &host.pieces[1])] {
        let [lo, hi] = chart
            .parameter_interval
            .map(|parameter| parameter.numerator as f64 / parameter.denominator as f64);
        let span = hi - lo;
        for local in samples {
            let source_u = lo + span * local[0];
            let (source_point, source_jacobian) = evaluate(&source, [source_u, local[1], local[2]]);
            let expected_point = apply_piece(affine, source_point);
            let (actual_point, actual_jacobian) = evaluate(&chart.block, local);
            for component in 0..3 {
                assert_close(actual_point[component], expected_point[component], 1e-12);
            }
            for axis in 0..3 {
                for (row, observed) in actual_jacobian[axis].iter().enumerate() {
                    let source_factor = if axis == 0 { span } else { 1.0 };
                    let expected = (0..3)
                        .map(|column| {
                            affine.matrix[row][column] as f64
                                * source_jacobian[axis][column]
                                * source_factor
                        })
                        .sum::<f64>();
                    assert_close(*observed, expected, 1e-11);
                }
            }
        }
    }

    // Both restrictions meet at the host node, while their one-sided x jets
    // reflect the distinct physical slopes.
    let (left_node, left_jet) = evaluate(&first.block, [1.0, 0.31, 0.62]);
    let (right_node, right_jet) = evaluate(&second.block, [0.0, 0.31, 0.62]);
    let (source_node, source_jet) = evaluate(&source, [0.5, 0.31, 0.62]);
    for piece in &host.pieces {
        let expected = apply_piece(piece, source_node);
        let actual = if piece.matrix[0][0] == host.pieces[0].matrix[0][0] {
            left_node
        } else {
            right_node
        };
        for component in 0..3 {
            assert_close(actual[component], expected[component], 1e-12);
        }
    }
    for component in 0..3 {
        assert_close(left_node[component], right_node[component], 1e-12);
    }
    for row in 0..3 {
        let left_expected = (0..3)
            .map(|column| host.pieces[0].matrix[row][column] as f64 * source_jet[0][column] * 0.5)
            .sum::<f64>();
        let right_expected = (0..3)
            .map(|column| host.pieces[1].matrix[row][column] as f64 * source_jet[0][column] * 0.5)
            .sum::<f64>();
        assert_close(left_jet[0][row], left_expected, 1e-11);
        assert_close(right_jet[0][row], right_expected, 1e-11);
    }
    assert!((right_jet[0][0] - left_jet[0][0]).abs() > 1e-3);
}

#[test]
fn root_free_linear_host_composition_is_exact_and_reports_no_crossing() {
    let source = linear_block("linear_source");
    let mut host = host_at_two();
    host.pieces[0] = piece(2, 1);
    let composition = compose_with_host(&source, &host).unwrap();
    assert_eq!(composition.representation, HostRepresentation::Unsplit);
    assert_eq!(composition.blocks.len(), 1);
    assert_eq!(composition.blocks[0].host_piece_index, 0);
    assert_eq!(composition.node_events.len(), 1);
    assert!(matches!(
        composition.node_events[0].status,
        HostNodeStatus::NotCrossed
    ));
    for q in [[0.12, 0.27, 0.83], [0.5, 0.5, 0.5], [0.9, 0.73, 0.14]] {
        let (source_point, _) = evaluate(&source, q);
        let expected = apply_piece(&host.pieces[0], source_point);
        let (actual, _) = evaluate(&composition.blocks[0].block, q);
        for component in 0..3 {
            assert_close(actual[component], expected[component], 1e-12);
        }
    }
}

#[test]
fn irrational_root_is_refused_for_monotone_nonuniform_quadratic() {
    // At x=1 the Bernstein polynomial X-W has coefficients [-1,0,3]
    // and a root (-1+sqrt(3))/2, which is irrational. X/W control positions
    // are [0,1,4] and increase strictly, isolating the irrational-root case.
    let source = quadratic_block("irrational_root", [0, 2, 4], [1, 2, 1]);
    let host = PiecewiseAffineHost {
        breaks: vec![1],
        pieces: vec![piece(2, 0), piece(3, -1)],
    };
    assert!(matches!(
        compose_with_host(&source, &host),
        Err(HostCompositionError::IrrationalPreimage)
    ));
}

#[test]
fn broken_and_degenerate_hosts_are_rejected() {
    let source = quadratic_block("rational_root", [1, 4, 3], [1, 2, 1]);
    let broken = PiecewiseAffineHost {
        breaks: vec![2],
        pieces: vec![piece(2, 1), piece(3, 0)],
    };
    assert!(matches!(
        compose_with_host(&source, &broken),
        Err(HostCompositionError::DiscontinuousHost)
    ));

    let degenerate = PiecewiseAffineHost {
        breaks: vec![2],
        pieces: vec![
            AffinePiece {
                matrix: [[2, 0, 0], [0, 2, 0], [0, 0, 0]],
                translation: [1, 0, 0],
            },
            piece(3, -1),
        ],
    };
    assert!(matches!(
        compose_with_host(&source, &degenerate),
        Err(HostCompositionError::NonPositiveDeterminant)
    ));
}
