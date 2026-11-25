//! Signed lift jets, nonhomogeneous boundary data and an independently rebuilt objective.
use solidcraft_trivariate::constrained_elasticity::solve_constrained;
use solidcraft_trivariate::constraint_system::{resolve, Lift, Ratio};
use solidcraft_trivariate::elasticity::Stiffness;
use solidcraft_trivariate::face_constraints::Row;
use solidcraft_trivariate::mechanical_sensitivity::{
    differentiate, LiftRate, MechanicalDerivative,
};

fn matrix(t: f64, rate: bool) -> Stiffness {
    let mut columns = vec![];
    let mut values = vec![];
    for i in 0..4 {
        for j in 0..4 {
            columns.push(j);
            values.push(std::array::from_fn(|a| {
                std::array::from_fn(|b| {
                    let base = if i == j && a == b {
                        7. + (i + a) as f64
                    } else if i == j {
                        0.25
                    } else if a == b {
                        -0.5
                    } else {
                        -0.06
                    };
                    let slope = if i == j && a == b {
                        0.4 + 0.1 * a as f64
                    } else {
                        0.03
                    };
                    if rate {
                        slope
                    } else {
                        base + t * slope
                    }
                })
            }));
        }
    }
    Stiffness {
        offsets: vec![0, 4, 8, 12, 16],
        columns,
        values,
        integrated_volume: 1.,
    }
}

fn lifts(num: i128, den: i128) -> [Lift; 3] {
    // Master weights proportional to 3 and 1+t². A dependent row uses the
    // rational midpoint weights and a nonzero value is imposed on root 0.
    let w0 = 3 * den * den;
    let w1 = den * den + num * num;
    let equation = [(
        2,
        Row {
            terms: vec![(0, w0), (1, w1)],
            denominator: w0 + w1,
        },
    )];
    std::array::from_fn(|c| {
        let prescribed = Ratio::new((c as i128 + 1) * num * num, den * den).unwrap();
        resolve(4, &equation, &[(0, prescribed)]).unwrap()
    })
}

fn load(t: f64) -> Vec<f64> {
    (0..12).map(|i| (i as f64 - 3.) * (1. + 0.2 * t)).collect()
}
fn gradient() -> Vec<f64> {
    (0..12).map(|i| (0.7 * i as f64).cos()).collect()
}
fn product(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn evaluate(num: i128, den: i128) -> (Vec<f64>, f64) {
    let t = num as f64 / den as f64;
    let lifts = lifts(num, den);
    let q = solve_constrained(&matrix(t, false), lifts.each_ref(), &load(t), 1e-13, 100).unwrap();
    let objective = product(&gradient(), &q.displacement_mm) + t * t / 2.;
    (q.displacement_mm, objective)
}

#[test]
fn tangent_and_adjoint_match_full_rebuilds_with_moving_constraints_and_lifting() {
    let t = 1.;
    let l = lifts(1, 1);
    let rates: [Vec<LiftRate>; 3] = std::array::from_fn(|c| {
        let mut r = vec![LiftRate::default(); 4];
        let g = (c + 1) as f64;
        r[0].offset = 2. * g;
        // (1+t²)/(4+t²) evaluated at t=1 has derivative 6/25.
        r[2].terms = vec![(0, 6. / 25.)];
        r[2].offset = (3. / 5. * 2. - 6. / 25.) * g;
        r
    });
    let fprime: Vec<_> = (0..12).map(|i| (i as f64 - 3.) * 0.2).collect();
    let result = differentiate(
        MechanicalDerivative {
            stiffness: &matrix(t, false),
            stiffness_rate: &matrix(t, true),
            force: &load(t),
            force_rate: &fprime,
            lifts: l.each_ref(),
            lift_rates: rates.each_ref().map(Vec::as_slice),
            objective_gradient: &gradient(),
            objective_partial: t,
            fixed_regime: true,
        },
        1e-13,
        100,
    )
    .unwrap();
    assert!((result.objective_tangent - result.objective_adjoint).abs() < 1e-11);
    let (_, base) = evaluate(1, 1);
    let mut remainders = vec![];
    for den in [1000, 2000, 4000] {
        let (plus, jplus) = evaluate(den + 1, den);
        let (minus, jminus) = evaluate(den - 1, den);
        let h = 1. / den as f64;
        let error = ((jplus - jminus) / (2. * h) - result.objective_adjoint).abs();
        assert!(error < 2. * h * h, "error {error} at h={h}");
        for ((a, b), derivative) in plus.iter().zip(minus).zip(&result.physical_tangent) {
            assert!(((a - b) / (2. * h) - derivative).abs() < 2. * h * h);
        }
        remainders.push((jplus - base - h * result.objective_adjoint).abs());
    }
    assert!(remainders.windows(2).all(|pair| pair[0] / pair[1] > 3.5));
}
