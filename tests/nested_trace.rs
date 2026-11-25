use solidcraft_trivariate::face_constraints::{nested_constraints, FaceNet, Span};
use solidcraft_trivariate::incidence::Orientation;

fn span_matrix(degree: usize, output_degree: usize, span: Span) -> (i128, Vec<Vec<i128>>) {
    let (denominator, rows) = match (degree, span) {
        (1, Span::Full) => (1, vec![vec![1, 0], vec![0, 1]]),
        (1, Span::Lower) => (2, vec![vec![2, 0], vec![1, 1]]),
        (1, Span::Upper) => (2, vec![vec![1, 1], vec![0, 2]]),
        (2, Span::Full) => (1, vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]]),
        (2, Span::Lower) => (4, vec![vec![4, 0, 0], vec![2, 2, 0], vec![1, 2, 1]]),
        (2, Span::Upper) => (4, vec![vec![1, 2, 1], vec![0, 2, 2], vec![0, 0, 4]]),
        _ => panic!("reference supports degree one and two"),
    };
    assert!(output_degree >= degree && output_degree <= 2);
    if output_degree == degree {
        return (denominator, rows);
    }
    let elevated = vec![
        rows[0].iter().map(|x| 2 * x).collect(),
        rows[0].iter().zip(&rows[1]).map(|(a, b)| a + b).collect(),
        rows[1].iter().map(|x| 2 * x).collect(),
    ];
    (denominator * 2, elevated)
}

fn master_net(degrees: [usize; 2]) -> FaceNet {
    let mut controls = Vec::new();
    for i in 0..=degrees[0] {
        for j in 0..=degrees[1] {
            let weight = 2 + (3 * i + 2 * j) as i64;
            let x = 1 + (2 * i + j) as i64;
            let y = 2 + (i * i + 3 * j) as i64;
            let z = 4 + (2 * i + j + i * j) as i64;
            controls.push([weight * x, weight * y, weight * z, weight]);
        }
    }
    FaceNet { degrees, controls }
}

fn restrict_controls(
    master: &FaceNet,
    output_degrees: [usize; 2],
    spans: [Span; 2],
) -> Vec<[i128; 4]> {
    let [p, q] = master.degrees;
    let (d0, t0) = span_matrix(p, output_degrees[0], spans[0]);
    let (d1, t1) = span_matrix(q, output_degrees[1], spans[1]);
    let mut result = Vec::with_capacity((output_degrees[0] + 1) * (output_degrees[1] + 1));
    for row0 in &t0 {
        for row1 in &t1 {
            let mut point = [0_i128; 4];
            for (a, &u_factor) in row0.iter().enumerate() {
                for (b, &v_factor) in row1.iter().enumerate() {
                    let factor = u_factor * v_factor;
                    let source = master.controls[a * (q + 1) + b];
                    for component in 0..4 {
                        point[component] += factor * i128::from(source[component]);
                    }
                }
            }
            result.push(point);
        }
    }
    assert!(d0 > 0 && d1 > 0);
    result
}

fn slave_storage_from_aligned(
    aligned: &[[i128; 4]],
    aligned_degrees: [usize; 2],
    orientation: Orientation,
    scale: i64,
) -> FaceNet {
    let slave_degrees = if orientation.swap {
        [aligned_degrees[1], aligned_degrees[0]]
    } else {
        aligned_degrees
    };
    let mut controls = vec![[0_i64; 4]; (slave_degrees[0] + 1) * (slave_degrees[1] + 1)];
    for i in 0..=aligned_degrees[0] {
        for j in 0..=aligned_degrees[1] {
            let a = if orientation.flip_first {
                aligned_degrees[0] - i
            } else {
                i
            };
            let b = if orientation.flip_second {
                aligned_degrees[1] - j
            } else {
                j
            };
            let (slave_i, slave_j) = if orientation.swap { (b, a) } else { (a, b) };
            let point = aligned[i * (aligned_degrees[1] + 1) + j];
            let index = slave_i * (slave_degrees[1] + 1) + slave_j;
            controls[index] = point.map(|value| {
                i64::try_from(value * i128::from(scale)).expect("reference coefficient fits i64")
            });
        }
    }
    FaceNet {
        degrees: slave_degrees,
        controls,
    }
}

fn restricted_weighted_field(
    master: &FaceNet,
    output_degrees: [usize; 2],
    spans: [Span; 2],
    field: &[i64],
) -> (Vec<i128>, Vec<i128>) {
    let [p, q] = master.degrees;
    assert_eq!(field.len(), master.controls.len());
    let (_, t0) = span_matrix(p, output_degrees[0], spans[0]);
    let (_, t1) = span_matrix(q, output_degrees[1], spans[1]);
    let mut numerator = Vec::new();
    let mut denominator = Vec::new();
    for row0 in &t0 {
        for row1 in &t1 {
            let mut n = 0_i128;
            let mut d = 0_i128;
            for (a, &u_factor) in row0.iter().enumerate() {
                for (b, &v_factor) in row1.iter().enumerate() {
                    let factor = u_factor * v_factor;
                    let index = a * (q + 1) + b;
                    let weight = i128::from(master.controls[index][3]);
                    d += factor * weight;
                    n += factor * weight * i128::from(field[index]);
                }
            }
            numerator.push(n);
            denominator.push(d);
        }
    }
    (numerator, denominator)
}

fn inverse_orient_scalar(
    aligned: &[i128],
    aligned_degrees: [usize; 2],
    orientation: Orientation,
) -> Vec<i128> {
    let slave_degrees = if orientation.swap {
        [aligned_degrees[1], aligned_degrees[0]]
    } else {
        aligned_degrees
    };
    let mut result = vec![0; aligned.len()];
    for i in 0..=aligned_degrees[0] {
        for j in 0..=aligned_degrees[1] {
            let (slave_i, slave_j) = if orientation.swap {
                let a = if orientation.flip_first {
                    aligned_degrees[0] - i
                } else {
                    i
                };
                let b = if orientation.flip_second {
                    aligned_degrees[1] - j
                } else {
                    j
                };
                (b, a)
            } else {
                (
                    if orientation.flip_first {
                        aligned_degrees[0] - i
                    } else {
                        i
                    },
                    if orientation.flip_second {
                        aligned_degrees[1] - j
                    } else {
                        j
                    },
                )
            };
            result[slave_i * (slave_degrees[1] + 1) + slave_j] =
                aligned[i * (aligned_degrees[1] + 1) + j];
        }
    }
    result
}

fn orientation(swap: bool, flip_first: bool, flip_second: bool) -> Orientation {
    Orientation {
        swap,
        flip_first,
        flip_second,
    }
}

fn check_fields(master: &FaceNet, slave: &FaceNet, spans: [Span; 2], orient: Orientation) {
    let rows = nested_constraints(master, slave, spans, orient).unwrap();
    assert_eq!(rows.len(), slave.controls.len());
    let fields = [
        (0..master.controls.len())
            .map(|i| (i * i + 3 * i + 1) as i64)
            .collect::<Vec<_>>(),
        (0..master.controls.len())
            .map(|i| 7 * i as i64 - 4)
            .collect::<Vec<_>>(),
        (0..master.controls.len())
            .map(|i| ((i % 3) as i64 - 1) * 11)
            .collect::<Vec<_>>(),
    ];
    for field in fields {
        let aligned_degrees = if orient.swap {
            [slave.degrees[1], slave.degrees[0]]
        } else {
            slave.degrees
        };
        let (aligned_numerator, aligned_denominator) =
            restricted_weighted_field(master, aligned_degrees, spans, &field);
        let expected_numerator = inverse_orient_scalar(&aligned_numerator, aligned_degrees, orient);
        let expected_denominator =
            inverse_orient_scalar(&aligned_denominator, aligned_degrees, orient);
        for (slave_index, row) in rows.iter().enumerate() {
            assert!(row.denominator > 0);
            assert!(row.terms.iter().all(|(index, _)| *index < field.len()));
            let actual_numerator: i128 = row
                .terms
                .iter()
                .map(|(index, coefficient)| coefficient * i128::from(field[*index]))
                .sum();
            assert_eq!(
                actual_numerator * expected_denominator[slave_index],
                expected_numerator[slave_index] * row.denominator,
                "field restriction differs at stored slave control {slave_index}"
            );
        }
    }
}

#[test]
fn all_eight_orientations_return_constraints_in_original_slave_order() {
    for degrees in [[1, 1], [1, 2], [2, 1], [2, 2]] {
        let master = master_net(degrees);
        let spans = [Span::Lower, Span::Upper];
        let aligned = restrict_controls(&master, degrees, spans);
        for swap in [false, true] {
            for flip_first in [false, true] {
                for flip_second in [false, true] {
                    let orient = orientation(swap, flip_first, flip_second);
                    let slave = slave_storage_from_aligned(&aligned, degrees, orient, 3);
                    check_fields(&master, &slave, spans, orient);
                }
            }
        }
    }
}

#[test]
fn degree_elevation_restrictions_match_independent_weighted_fields() {
    for master_degrees in [[1, 1], [1, 2], [2, 1]] {
        let master = master_net(master_degrees);
        let spans = [Span::Lower, Span::Upper];
        let output_degrees = [2, 2];
        let aligned = restrict_controls(&master, output_degrees, spans);
        for swap in [false, true] {
            for flip_first in [false, true] {
                for flip_second in [false, true] {
                    let orient = orientation(swap, flip_first, flip_second);
                    let slave = slave_storage_from_aligned(&aligned, output_degrees, orient, 3);
                    check_fields(&master, &slave, spans, orient);
                }
            }
        }
    }
}

#[test]
fn full_lower_upper_spans_and_mixed_degree_restrictions_match_independent_fields() {
    let master = master_net([2, 1]);
    for first in [Span::Full, Span::Lower, Span::Upper] {
        for second in [Span::Full, Span::Lower, Span::Upper] {
            let spans = [first, second];
            let aligned = restrict_controls(&master, master.degrees, spans);
            let slave = slave_storage_from_aligned(
                &aligned,
                master.degrees,
                orientation(false, false, false),
                5,
            );
            check_fields(&master, &slave, spans, orientation(false, false, false));
        }
    }
}

#[test]
fn common_homogeneous_rescaling_preserves_constraint_values() {
    let master = master_net([2, 1]);
    let spans = [Span::Upper, Span::Lower];
    let aligned = restrict_controls(&master, master.degrees, spans);
    let orient = orientation(true, true, false);
    let slave = slave_storage_from_aligned(&aligned, master.degrees, orient, 3);
    let base = nested_constraints(&master, &slave, spans, orient).unwrap();

    let mut scaled_master = master.clone();
    for point in &mut scaled_master.controls {
        for component in point {
            *component *= 7;
        }
    }
    let mut scaled_slave = slave.clone();
    for point in &mut scaled_slave.controls {
        for component in point {
            *component *= 11;
        }
    }
    let scaled = nested_constraints(&scaled_master, &scaled_slave, spans, orient).unwrap();
    let field: Vec<i64> = (0..master.controls.len())
        .map(|i| i as i64 * 5 - 2)
        .collect();
    for (row_index, (left, right)) in base.iter().zip(&scaled).enumerate() {
        let evaluate = |row: &solidcraft_trivariate::face_constraints::Row| -> (i128, i128) {
            let numerator = row
                .terms
                .iter()
                .map(|(index, coefficient)| coefficient * i128::from(field[*index]))
                .sum();
            (numerator, row.denominator)
        };
        let (left_n, left_d) = evaluate(left);
        let (right_n, right_d) = evaluate(right);
        assert_eq!(
            left_n * right_d,
            right_n * left_d,
            "rescaled row {row_index}"
        );
    }
}

#[test]
fn altered_slave_weight_or_homogeneous_coefficient_is_rejected() {
    let master = master_net([2, 1]);
    let spans = [Span::Lower, Span::Upper];
    let aligned = restrict_controls(&master, master.degrees, spans);
    let orient = orientation(false, true, false);
    let slave = slave_storage_from_aligned(&aligned, master.degrees, orient, 2);

    let mut bad_weight = slave.clone();
    bad_weight.controls[2][3] += 1;
    assert!(nested_constraints(&master, &bad_weight, spans, orient).is_err());

    let mut bad_point = slave;
    bad_point.controls[4][0] += 1;
    assert!(nested_constraints(&master, &bad_point, spans, orient).is_err());
}
