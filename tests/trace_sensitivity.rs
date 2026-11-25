use solidcraft_trivariate::face_constraints::{FaceNet, Span};
use solidcraft_trivariate::incidence::Orientation;
use solidcraft_trivariate::trace_sensitivity::nested_weight_jets;

fn fixture() -> (FaceNet, FaceNet) {
    let master: Vec<_> = [[1, 0, 1], [1, 1, 1], [0, 2, 2]]
        .into_iter()
        .flat_map(|[x, y, w]| [0, 1].map(|z| [x, y, z * w, w]))
        .collect();
    let mut slave = vec![];
    for stencil in [[4, 0, 0], [2, 2, 0], [1, 2, 1]] {
        for z in 0..2 {
            slave.push(std::array::from_fn(|c| {
                (0..3).map(|i| stencil[i] * master[2 * i + z][c]).sum()
            }));
        }
    }
    (
        FaceNet {
            degrees: [2, 1],
            controls: master,
        },
        FaceNet {
            degrees: [2, 1],
            controls: slave,
        },
    )
}

#[test]
fn derivative_matches_analytic_trace_and_second_order_taylor_remainder() {
    let (master, slave) = fixture();
    let identity = Orientation {
        swap: false,
        flip_first: false,
        flip_second: false,
    };
    let spans = [Span::Lower, Span::Full];
    // Independent scalar coefficient at the midpoint of the original curve.
    let f = |t: f64| (2. + 3. * (1. + t * t)) / (4. + t * t);
    let t: f64 = 1.;
    let exact_derivative =
        2. * t * (3. * (4. + t * t) - (2. + 3. * (1. + t * t))) / (4. + t * t).powi(2);
    let dw = [0., 0., 0., 0., 2., 2.];
    let s_dw = [0., 0., 0., 0., 2., 2.];
    let rows = nested_weight_jets(&master, &dw, &slave, &s_dw, spans, identity).unwrap();
    let values = [0., 0., 1., 1., 3., 3.];
    let row = &rows[4];
    let value: f64 = row.iter().map(|c| c.value * values[c.master]).sum();
    let derivative: f64 = row.iter().map(|c| c.derivative * values[c.master]).sum();
    assert!((value - f(t)).abs() < 1e-14 && (derivative - exact_derivative).abs() < 1e-14);
    for r in &rows {
        assert!(r.iter().map(|c| c.derivative).sum::<f64>().abs() < 1e-15);
    }
    let remainder = |h: f64| (f(t + h) - value - h * derivative).abs();
    assert!(remainder(1e-3) / remainder(5e-4) > 3.5);
    let mut mismatched = s_dw;
    mismatched[1] += 1.;
    assert!(nested_weight_jets(&master, &dw, &slave, &mismatched, spans, identity).is_err());
    let gauge: Vec<_> = slave
        .controls
        .iter()
        .zip(s_dw)
        .map(|(p, d)| d + 2. * p[3] as f64)
        .collect();
    assert_eq!(
        nested_weight_jets(&master, &dw, &slave, &gauge, spans, identity).unwrap(),
        rows
    );
}

#[test]
fn common_weight_scale_has_zero_jet_under_all_face_orientations() {
    let (master, canonical) = fixture();
    let d: Vec<_> = master.controls.iter().map(|p| p[3] as f64).collect();
    for swap in [false, true] {
        for flip_first in [false, true] {
            for flip_second in [false, true] {
                let o = Orientation {
                    swap,
                    flip_first,
                    flip_second,
                };
                let mut slave = FaceNet {
                    degrees: if swap { [1, 2] } else { [2, 1] },
                    controls: vec![[0; 4]; 6],
                };
                for i in 0..3 {
                    for j in 0..2 {
                        let (a, b) = (
                            if flip_first { 2 - i } else { i },
                            if flip_second { 1 - j } else { j },
                        );
                        let local = if swap { 3 * b + a } else { 2 * a + b };
                        slave.controls[local] = canonical.controls[2 * i + j];
                    }
                }
                let dw: Vec<_> = slave.controls.iter().map(|p| 3. * p[3] as f64).collect();
                let rows =
                    nested_weight_jets(&master, &d, &slave, &dw, [Span::Lower, Span::Full], o)
                        .unwrap();
                assert!(rows
                    .into_iter()
                    .flatten()
                    .all(|c| c.derivative.abs() < 1e-15));
            }
        }
    }
}
