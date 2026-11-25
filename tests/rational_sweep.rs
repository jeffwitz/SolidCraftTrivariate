use solidcraft_trivariate::{rational_sweep, RationalBezierMap};

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.17e} != {b:.17e}, tol={tol:.3e}");
}

fn physical(t: f64, h: f64, [u, v, w]: [f64; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let angle = 2.0 * (t * u).atan();
    let radius = 2.0 - v;
    let (s, c) = angle.sin_cos();
    let dangle = 2.0 * t / (1.0 + (t * u).powi(2));
    (
        [radius * c, radius * s, h * w],
        [
            [-radius * s * dangle, radius * c * dangle, 0.0],
            [-c, -s, 0.0],
            [0.0, 0.0, h],
        ],
    )
}

#[test]
fn restricted_rational_maps_match_independent_physical_formula() {
    for t in [0.37, 1.2, 1.0 + 2.0_f64.sqrt() + 0.2, 3.8] {
        let h = 7.25;
        let sweep = rational_sweep::build_rational_sweep(t, h).unwrap();
        let mut count = 0;
        for block in &sweep.blocks {
            let map = RationalBezierMap::new(block.degrees, block.controls.clone()).unwrap();
            let [a, b] = block.parameter_interval;
            let delta = b - a;
            for s in [0.0, 0.173, 0.5, 0.827, 1.0] {
                for v in [0.0, 0.31, 1.0] {
                    for w in [0.0, 0.63, 1.0] {
                        let actual = map.evaluate([s, v, w]).unwrap();
                        let (position, gradients) = physical(t, h, [a + delta * s, v, w]);
                        for c in 0..3 {
                            close(actual.position[c], position[c], 2e-13);
                            close(actual.jacobian_columns[0][c], gradients[0][c] * delta, 2e-12);
                            close(actual.jacobian_columns[1][c], gradients[1][c], 2e-13);
                            close(actual.jacobian_columns[2][c], gradients[2][c], 2e-13);
                        }
                        count += 1;
                    }
                }
            }
        }
        assert_eq!(count, sweep.blocks.len() * 45);
    }
}

#[test]
fn analytic_jets_and_fixed_volume_taylor_checks() {
    let (t, h) = (1.37, 8.5);
    let sweep = rational_sweep::build_rational_sweep(t, h).unwrap();
    let eval = rational_sweep::evaluate_sweep(t, h, [0.37, 0.42, 0.61]).unwrap();
    let step = 1e-5;
    let plus = rational_sweep::evaluate_sweep(t + step, h, [0.37, 0.42, 0.61]).unwrap();
    let minus = rational_sweep::evaluate_sweep(t - step, h, [0.37, 0.42, 0.61]).unwrap();
    for c in 0..3 {
        close((plus.position_mm[c] - minus.position_mm[c]) / (2.0 * step), eval.dposition_dt_mm[c], 2e-9);
        close(eval.position_mm[c], physical(t, h, [0.37, 0.42, 0.61]).0[c], 2e-14);
    }
    close(sweep.volume_mm3, 3.0 * h * t.atan(), 2e-14);
    close(sweep.dvolume_dt_mm3, 3.0 * h / (1.0 + t * t), 2e-14);
    close(sweep.dvolume_dh_mm2, 3.0 * t.atan(), 2e-14);
    let vh_plus = rational_sweep::build_rational_sweep(t, h + step).unwrap().volume_mm3;
    let vh_minus = rational_sweep::build_rational_sweep(t, h - step).unwrap().volume_mm3;
    close((vh_plus - vh_minus) / (2.0 * step), sweep.dvolume_dh_mm2, 2e-10);
    let vt_plus = rational_sweep::build_rational_sweep(t + step, h).unwrap().volume_mm3;
    let vt_minus = rational_sweep::build_rational_sweep(t - step, h).unwrap().volume_mm3;
    close((vt_plus - vt_minus) / (2.0 * step), sweep.dvolume_dt_mm3, 2e-9);

    // The volume is linear in height, so its second-order Taylor remainder vanishes.
    let dh = 2e-3;
    let predicted = sweep.volume_mm3 + sweep.dvolume_dh_mm2 * dh;
    close(rational_sweep::build_rational_sweep(t, h + dh).unwrap().volume_mm3, predicted, 2e-14);
}

#[test]
fn coefficient_derivatives_converge_under_centered_differences() {
    let (t, h) = (1.31, 6.75);
    let base = rational_sweep::build_rational_sweep(t, h).unwrap();
    let jet = base.blocks[0].dcontrols_dt.as_ref().unwrap();
    let height_jet = &base.blocks[0].dcontrols_dh;
    let steps = [0.001, 0.0005, 0.00025];
    let mut errors = Vec::new();
    for step in steps {
        let p = rational_sweep::build_rational_sweep(t + step, h).unwrap();
        let m = rational_sweep::build_rational_sweep(t - step, h).unwrap();
        let mut max_error: f64 = 0.0;
        let hp = rational_sweep::build_rational_sweep(t, h + step).unwrap();
        let hm = rational_sweep::build_rational_sweep(t, h - step).unwrap();
        for index in [0, 4, 8, 11] {
            for (component, expected) in jet[index].iter().enumerate() {
                let fd = (p.blocks[0].controls[index][component] - m.blocks[0].controls[index][component]) / (2.0 * step);
                max_error = max_error.max((fd - expected).abs());
                let fd_height = (hp.blocks[0].controls[index][component] - hm.blocks[0].controls[index][component]) / (2.0 * step);
                max_error = max_error.max((fd_height - height_jet[index][component]).abs());
            }
        }
        errors.push(max_error);
    }
    assert!(errors.iter().all(|e| *e <= 1e-7), "centered derivative errors: {errors:?}");
    // At these steps truncation error is already below the binary64 subtraction
    // floor; require the measured plateau to remain within the declared bound.
}

#[test]
fn chart_event_is_bit_exact_and_physical_gradient_is_continuous() {
    let event = 1.0 + 2.0_f64.sqrt();
    let below = f64::from_bits(event.to_bits() - 1);
    let above = f64::from_bits(event.to_bits() + 1);
    let (one, exact, two) = (
        rational_sweep::build_rational_sweep(below, 4.0).unwrap(),
        rational_sweep::build_rational_sweep(event, 4.0).unwrap(),
        rational_sweep::build_rational_sweep(above, 4.0).unwrap(),
    );
    assert_eq!(one.blocks.len(), 1);
    assert_eq!(two.blocks.len(), 2);
    assert!(!one.representation_event && exact.representation_event && !two.representation_event);
    assert!(exact.blocks.iter().all(|b| b.dcontrols_dt.is_none()));
    assert_eq!(one.domain.family, exact.domain.family);
    assert_eq!(exact.domain.family, two.domain.family);
    assert_eq!(one.domain.start_angle_rad, exact.domain.start_angle_rad);
    assert_eq!(exact.domain.start_angle_rad, two.domain.start_angle_rad);
    close(one.domain.end_angle_rad, exact.domain.end_angle_rad, 2e-15);
    close(exact.domain.end_angle_rad, two.domain.end_angle_rad, 2e-15);
    for point in [[0.23, 0.4, 0.7], [0.76, 0.8, 0.2]] {
        let a = rational_sweep::evaluate_sweep(below, 4.0, point).unwrap();
        let b = rational_sweep::evaluate_sweep(event, 4.0, point).unwrap();
        let c = rational_sweep::evaluate_sweep(above, 4.0, point).unwrap();
        for i in 0..3 {
            close(a.dposition_dt_mm[i], b.dposition_dt_mm[i], 2e-15);
            close(b.dposition_dt_mm[i], c.dposition_dt_mm[i], 2e-15);
        }
    }
}

#[test]
fn invalid_parameters_are_refused() {
    use rational_sweep::SweepError::{NonFiniteParameter, ParameterOutOfDomain};
    for (t, h) in [(f64::NAN, 1.0), (1.0, f64::INFINITY), (f64::NEG_INFINITY, 1.0)] {
        assert_eq!(rational_sweep::build_rational_sweep(t, h), Err(NonFiniteParameter));
    }
    for (t, h) in [(0.0, 1.0), (4.01, 1.0), (1.0, 0.49), (1.0, 32.1)] {
        assert_eq!(rational_sweep::build_rational_sweep(t, h), Err(ParameterOutOfDomain));
    }
    assert_eq!(rational_sweep::evaluate_sweep(1.0, 1.0, [f64::NAN, 0.0, 0.0]), Err(NonFiniteParameter));
    assert_eq!(rational_sweep::evaluate_sweep(1.0, 1.0, [1.1, 0.0, 0.0]), Err(ParameterOutOfDomain));
}
