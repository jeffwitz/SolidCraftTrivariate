use solidcraft_trivariate::{
    adaptive_partition::{AdaptivePartition, LeafKey},
    adaptive_space::build_adaptive_space,
    analysis_space::{AnalysisBlock, AnalysisSpace},
    constrained_elasticity::{solve_preconditioned, ReducedPreconditioner},
    constraint_system::{resolve, Ratio},
    elasticity::Stiffness,
    incidence::{BlockSide, Incidence},
    multilevel::{MountedSpace, MultilevelPreconditioner},
};
use std::collections::BTreeMap;

fn diagonal(n: usize) -> Stiffness {
    Stiffness {
        offsets: (0..=n).collect(),
        columns: (0..n).collect(),
        values: (0..n)
            .map(|i| {
                let a = 3. + (i % 5) as f64;
                [[a, 0.2, -0.1], [0.2, a + 2., 0.03], [-0.1, 0.03, a + 1.]]
            })
            .collect(),
        integrated_volume: 1.,
    }
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[test]
fn coarse_factor_solves_coupled_components_and_rejects_invalid_inputs() {
    let k = diagonal(7);
    let lift = resolve(7, &[], &[(1, Ratio::zero())]).unwrap();
    assert!(MultilevelPreconditioner::new(&k, &lift, 32).is_err());
    let mut invalid = diagonal(7);
    invalid.values[0][0][1] += 0.1;
    assert!(MultilevelPreconditioner::new(&invalid, &lift, 1 << 20).is_err());
    invalid = diagonal(7);
    invalid.values[0][0][0] = -1.;
    assert!(MultilevelPreconditioner::new(&invalid, &lift, 1 << 20).is_err());
    let mg = MultilevelPreconditioner::new(&k, &lift, 1 << 20).unwrap();
    assert!(mg.apply(&[0.; 3]).is_err());
    let force: Vec<_> = (0..21).map(|i| (i as f64 + 0.35).cos()).collect();
    let answer = solve_preconditioned(&k, [&lift; 3], &force, 1e-12, 2, &mg).unwrap();
    assert_eq!(answer.iterations, 1);
    assert!(answer.relative_residual < 1e-14);
}

#[test]
fn nested_cycle_is_positive_symmetric_and_preserves_essential_values() {
    let root = AnalysisSpace {
        ndof: 27,
        blocks: vec![AnalysisBlock {
            id: "cube".into(),
            controls: (0..27)
                .map(|i| [(i / 9) * 5, ((i / 3) % 3) * 5, (i % 3) * 5, 10])
                .collect(),
            dofs: (0..27).collect(),
        }],
        boundary_dofs: BTreeMap::new(),
    };
    let sides = (0..3)
        .flat_map(|axis| {
            [false, true].map(move |upper| BlockSide {
                block: "cube".into(),
                axis,
                upper,
            })
        })
        .collect();
    let p = AdaptivePartition::new(
        &["cube".into()],
        &Incidence {
            interfaces: vec![],
            exterior: sides,
        },
        3,
        100,
    )
    .unwrap();
    let a = build_adaptive_space(&root, &p).unwrap();
    let mut q = p.clone();
    q.refine(&[LeafKey::root("cube")]).unwrap();
    let b = build_adaptive_space(&root, &q).unwrap();
    let coarse = diagonal(a.space.ndof);
    let fine = diagonal(b.space.ndof);
    let mut mg = MultilevelPreconditioner::new(&coarse, &a.lift, 1 << 20).unwrap();
    let source = MountedSpace {
        adaptive: &a,
        partition: &p,
        lift: &a.lift,
    };
    let target = MountedSpace {
        adaptive: &b,
        partition: &q,
        lift: &b.lift,
    };
    let fixed = resolve(b.space.ndof, &[], &[(0, Ratio::zero())]).unwrap();
    assert!(mg
        .add_level(
            source,
            MountedSpace {
                lift: &fixed,
                ..target
            },
            &fine
        )
        .is_err());
    assert_eq!(mg.level_count(), 1);
    mg.add_level(source, target, &fine).unwrap();
    let x: Vec<_> = (0..fine.size()).map(|i| ((i + 7) as f64).sin()).collect();
    let y: Vec<_> = (0..fine.size())
        .map(|i| ((13 * i + 1) as f64).cos())
        .collect();
    let bx = mg.apply(&x).unwrap();
    let by = mg.apply(&y).unwrap();
    assert!(dot(&x, &bx) > 0. && dot(&y, &by) > 0.);
    assert!((dot(&x, &by) - dot(&y, &bx)).abs() < 1e-11);
    // Partition of unity for the scalar Bernstein subdivision on this unit cube.
    // A constant component vector is unchanged by P, hence <P1,r> = sum(r).
    let transfer =
        solidcraft_trivariate::adaptive_transfer::build_transfer(&a, &p, &b, &q).unwrap();
    let mut transpose_times_x = vec![0.; a.space.ndof];
    for (i, row) in transfer.rows.iter().enumerate() {
        for &(j, w) in row {
            transpose_times_x[j] += x[i] * w.numerator() as f64 / w.denominator() as f64;
        }
    }
    assert!(
        (transpose_times_x.iter().sum::<f64>() - x[..b.space.ndof].iter().sum::<f64>()).abs()
            < 1e-12
    );
    let solved = solve_preconditioned(&fine, [&b.lift; 3], &x, 1e-11, 200, &mg).unwrap();
    assert!(solved.relative_residual <= 1e-11);
    assert!(mg.add_level(source, target, &fine).is_err());
}
