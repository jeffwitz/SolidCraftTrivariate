use solidcraft_trivariate::bracket::BracketBlock;
use solidcraft_trivariate::mapped_cut_propagation::{
    propagate_exact_cuts, CutParameter, CutRequest, FaceMap,
};

fn adjacent_blocks() -> Vec<BracketBlock> {
    let mut curved = Vec::new();
    for i in 0..=2_i64 {
        for j in 0..=1_i64 {
            for k in 0..=1_i64 {
                curved.push([i * 2, j * 2, k * 2, 2]);
            }
        }
    }
    let mut box_controls = Vec::new();
    for i in 0..=1_i64 {
        for j in 0..=1_i64 {
            for k in 0..=1_i64 {
                box_controls.push([2 + i, j, k, 1]);
            }
        }
    }
    vec![
        BracketBlock {
            id: "host".into(),
            degrees: [2, 1, 1],
            controls: curved,
        },
        BracketBlock {
            id: "neighbor".into(),
            degrees: [1, 1, 1],
            controls: box_controls,
        },
    ]
}

#[test]
fn rational_cut_propagates_through_swapped_reversed_trace_and_splits_exactly() {
    let blocks = adjacent_blocks();
    let face = FaceMap {
        a_block: 0,
        a_normal_axis: 0,
        a_upper: true,
        b_block: 1,
        b_normal_axis: 0,
        b_upper: false,
        a_to_b: [(2, true), (1, false)],
        isoparametric_admitted: true,
        provenance: "host/neighbor face".into(),
    };
    let plan = propagate_exact_cuts(
        &blocks,
        &[face],
        &[CutRequest {
            block: 0,
            axis: 1,
            parameter: CutParameter::new(1, 3).unwrap(),
            provenance: "bore separator".into(),
        }],
        4,
    )
    .unwrap();
    assert_eq!(plan.requested_cut_count, 1);
    assert_eq!(plan.propagated_cut_count, 2);
    assert_eq!(plan.child_count, 4);
    assert_eq!(plan.blocks.len(), 4);
    assert!(plan
        .cuts
        .iter()
        .any(|c| c.block == 0 && c.axis == 1 && c.parameter == CutParameter::new(1, 3).unwrap()));
    assert!(plan
        .cuts
        .iter()
        .any(|c| c.block == 1 && c.axis == 2 && c.parameter == CutParameter::new(2, 3).unwrap()));
    assert_eq!(
        plan.source_to_children
            .iter()
            .map(|x| x.1.len())
            .collect::<Vec<_>>(),
        [2, 2]
    );
    assert!(plan
        .blocks
        .iter()
        .all(|b| b.controls.iter().all(|p| p[3] > 0)));
}

#[test]
fn incompatible_trace_and_budget_are_explicit_refusals() {
    let blocks = adjacent_blocks();
    let mut face = FaceMap {
        a_block: 0,
        a_normal_axis: 0,
        a_upper: true,
        b_block: 1,
        b_normal_axis: 0,
        b_upper: false,
        a_to_b: [(1, false), (2, false)],
        isoparametric_admitted: false,
        provenance: "unmatched trace".into(),
    };
    let request = [CutRequest {
        block: 0,
        axis: 1,
        parameter: CutParameter::new(1, 2).unwrap(),
        provenance: "test".into(),
    }];
    assert_eq!(
        propagate_exact_cuts(&blocks, &[face.clone()], &request, 10),
        Err("mapped_cut_incompatible_trace")
    );
    face.isoparametric_admitted = true;
    assert_eq!(
        propagate_exact_cuts(&blocks, &[face], &request, 3),
        Err("mapped_cut_block_budget")
    );
}

#[test]
fn duplicate_requests_close_once_and_preserve_all_origins() {
    let blocks = adjacent_blocks();
    let face = FaceMap {
        a_block: 0,
        a_normal_axis: 0,
        a_upper: true,
        b_block: 1,
        b_normal_axis: 0,
        b_upper: false,
        a_to_b: [(1, false), (2, false)],
        isoparametric_admitted: true,
        provenance: "shared".into(),
    };
    let p = CutParameter::new(1, 2).unwrap();
    let requests = [
        CutRequest {
            block: 0,
            axis: 1,
            parameter: p,
            provenance: "bore-a".into(),
        },
        CutRequest {
            block: 0,
            axis: 1,
            parameter: CutParameter::new(2, 4).unwrap(),
            provenance: "bore-b".into(),
        },
    ];
    let plan = propagate_exact_cuts(&blocks, &[face], &requests, 4).unwrap();
    assert_eq!(plan.propagated_cut_count, 2);
    assert_eq!(plan.child_count, 4);
    assert_eq!(
        plan.cuts
            .iter()
            .find(|c| c.block == 0)
            .unwrap()
            .provenance
            .len(),
        2
    );
}
