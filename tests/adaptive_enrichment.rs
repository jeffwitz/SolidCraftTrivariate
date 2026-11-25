#[path = "../examples/support/enrichment.rs"]
mod enrichment;

use enrichment::{bulk_mark, difference_energy_by_parent, galerkin_identity_relative_error};
use solidcraft_trivariate::{
    adaptive_partition::{AdaptivePartition, LeafKey},
    adaptive_space::build_adaptive_space,
    analysis_space::{AnalysisBlock, AnalysisSpace},
    incidence::{BlockSide, Incidence},
};
use std::collections::BTreeMap;

fn cube() -> (AnalysisSpace, Incidence) {
    let mut controls = Vec::with_capacity(27);
    let mut dofs = Vec::with_capacity(27);
    for i in 0..3_i128 {
        for j in 0..3_i128 {
            for k in 0..3_i128 {
                controls.push([i, j, k, 2]);
                dofs.push(dofs.len());
            }
        }
    }
    let sides = (0..3)
        .flat_map(|axis| {
            [false, true].map(move |upper| BlockSide {
                block: "cube".to_owned(),
                axis,
                upper,
            })
        })
        .collect::<Vec<_>>();
    (
        AnalysisSpace {
            ndof: 27,
            blocks: vec![AnalysisBlock {
                id: "cube".to_owned(),
                controls,
                dofs,
            }],
            boundary_dofs: BTreeMap::new(),
        },
        Incidence {
            interfaces: Vec::new(),
            exterior: sides,
        },
    )
}

#[test]
fn local_difference_energy_and_dorfler_marking_match_affine_cube_oracle() {
    let (parent, incidence) = cube();
    let current_partition = AdaptivePartition::new(&["cube".to_owned()], &incidence, 1, 8).unwrap();
    let mut target_partition = current_partition.clone();
    let all_leaves = target_partition
        .leaves()
        .iter()
        .cloned()
        .collect::<Vec<_>>();
    target_partition.refine(&all_leaves).unwrap();
    let target = build_adaptive_space(&parent, &target_partition).unwrap();

    // The exact affine field u=(x,0,0) has constant epsilon_xx=1. With
    // E=2, nu=1/4, its strain energy density is lambda+2*mu=12/5.
    let mut delta = vec![0.0; 3 * target.space.ndof];
    for block in &target.space.blocks {
        for (&dof, control) in block.dofs.iter().zip(&block.controls) {
            delta[3 * dof] = control[0] as f64 / control[3] as f64;
        }
    }
    let energy = difference_energy_by_parent(
        &target,
        &target_partition,
        &current_partition,
        &delta,
        2.0,
        0.25,
        8,
    )
    .unwrap();
    assert_eq!(energy.len(), 1);
    assert!((energy.values().sum::<f64>() - 2.4).abs() < 1e-12);

    let marked = bulk_mark(&energy, 0.5).unwrap();
    assert_eq!(
        marked,
        vec![current_partition.leaves().iter().next().unwrap().clone()]
    );
    let tied = (0..4)
        .map(|index| {
            (
                LeafKey {
                    root: "cube".to_owned(),
                    level: 1,
                    index: [index, 0, 0],
                },
                1.0,
            )
        })
        .collect::<BTreeMap<_, _>>();
    assert_eq!(
        bulk_mark(&tied, 0.5).unwrap(),
        tied.keys().take(2).cloned().collect::<Vec<_>>()
    );
    assert!(bulk_mark(&BTreeMap::new(), 0.5).unwrap().is_empty());
    assert!(bulk_mark(&energy, f64::NAN).is_err());
}

#[test]
fn nested_galerkin_energy_identity_reports_relative_discrepancy() {
    assert!(galerkin_identity_relative_error(2.0, 3.0, 4.0).unwrap() < 1e-15);
    assert!(galerkin_identity_relative_error(2.0, 3.0, 4.01).unwrap() > 1e-3);
    assert!(galerkin_identity_relative_error(0.0, 1.0, 1.0).is_err());
}
