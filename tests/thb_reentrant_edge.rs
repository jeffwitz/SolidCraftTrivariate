use solidcraft_trivariate::thb_reentrant_edge::{
    build_hierarchy, build_hierarchy_with_widths, build_uniform,
};
use std::collections::BTreeMap;

fn trace(h: &solidcraft_trivariate::thb_reentrant_edge::ReentrantHierarchy, p: usize, level: usize, uv: [f64; 2]) -> BTreeMap<usize, f64> {
    let mut result = BTreeMap::new();
    for (id, value, _) in h.sample_patch(p, level, uv).unwrap() {
        if value.abs() > 1e-14 {
            *result.entry(id).or_insert(0.0) += value;
        }
    }
    result
}

#[test]
fn level1_thb_smoke_has_exact_partition_transfer_and_matching_seams() {
    for level in 0..=1 {
        let h = build_hierarchy(level).unwrap();
        assert_eq!(h.maximum_level, level);
        assert!(h.audit_exact_partition().unwrap());
        assert!(h.audit_exact_transfer().unwrap(), "root transfer failed at level {level}");
        assert!(h.audit_leaf_cover().unwrap());
        assert!(h.active_levels_by_leaf().unwrap().iter().all(|levels| levels.len() <= 2));
        assert!(h.exact_seam_trace_pairs > 0);
        assert!(h.active_scalar_dofs <= h.uniform_fine_scalar_dofs);
        if level > 0 {
            assert!(h.leaves > build_hierarchy(level - 1).unwrap().leaves);
        }
        for t in [0.173, 0.417, 0.731] {
            let southwest_u0 = trace(&h, 0, level, [0.0, t]);
            let southeast_v0 = trace(&h, 1, level, [t, 0.0]);
            assert_eq!(southwest_u0, southeast_v0);
            let southwest_v0 = trace(&h, 0, level, [t, 0.0]);
            let northwest_u0 = trace(&h, 2, level, [0.0, t]);
            assert_eq!(southwest_v0, northwest_u0);
        }
    }
}

#[test]
#[ignore = "release qualification: exact multilevel transfer, level-3 admissibility and uniform extraction"]
fn level3_thb_transfer_and_admissibility_qualification() {
    for level in 0..=3 {
        let h = build_hierarchy(level).unwrap();
        assert!(h.audit_exact_partition().unwrap());
        assert!(h.audit_exact_transfer().unwrap());
        assert!(h.audit_leaf_cover().unwrap());
        if level == 3 {
            assert!(h.audit_nested_transfers().unwrap());
            let mut histogram = BTreeMap::new();
            for levels in h.active_levels_by_leaf().unwrap() {
                let span = levels.last().unwrap() - levels[0] + 1;
                *histogram.entry((levels[0], *levels.last().unwrap(), span, levels.len())).or_insert(0usize) += 1;
            }
            assert_eq!(histogram, BTreeMap::from([
                ((0, 0, 1, 1), 144),
                ((0, 1, 2, 2), 117),
                ((1, 1, 1, 1), 27),
                ((1, 2, 2, 2), 117),
                ((2, 2, 1, 1), 27),
                ((2, 3, 2, 2), 117),
                ((3, 3, 1, 1), 75),
            ]));
        }
    }
}

#[test]
#[ignore = "release qualification: wide boxes and level-3 Bernstein extraction"]
fn wider_nested_boxes_are_admitted_and_uniform_spaces_are_controls() {
    let narrow = build_hierarchy_with_widths(3, &[4, 4, 4]).unwrap();
    let wide = build_hierarchy_with_widths(3, &[6, 6, 6]).unwrap();
    let uniform = build_uniform(3).unwrap();
    assert!(wide.audit_exact_partition().unwrap());
    assert!(wide.audit_exact_transfer().unwrap());
    assert!(wide.audit_leaf_cover().unwrap());
    assert!(wide.active_levels_by_leaf().unwrap().iter().all(|levels| levels.len() <= 2));
    assert_eq!(wide.active_scalar_dofs, 1277);
    assert_eq!(wide.leaves, 1164);
    let wide_leaves = wide.integration_leaves().unwrap();
    assert_eq!(wide_leaves.len(), wide.leaves);
    assert!(wide_leaves.iter().all(|leaf| leaf.bernstein.len() == 16
        && leaf.bernstein.iter().all(|row| row.len() == leaf.dofs.len())
        && leaf.bernstein.iter().all(|row| (row.iter().sum::<f64>() - 1.0).abs() < 1e-12)));
    assert!(wide.active_scalar_dofs > narrow.active_scalar_dofs);
    assert!(wide.active_scalar_dofs < uniform.active_scalar_dofs);
    assert_eq!(uniform.active_scalar_dofs, uniform.uniform_fine_scalar_dofs);
    assert!(uniform.audit_exact_transfer().is_err());
    assert!(uniform.audit_nested_transfers().is_err());
    assert!(build_hierarchy_with_widths(3, &[3, 4, 4]).is_err());
    assert!(build_hierarchy_with_widths(3, &[4, 4]).is_err());
    assert!(build_hierarchy(4).is_err());
    assert!(build_hierarchy(usize::MAX).is_err());
    assert!(build_uniform(usize::MAX).is_err());
    let leaves = uniform.integration_leaves().unwrap();
    assert_eq!(leaves.len(), uniform.leaves);
    assert!(leaves.iter().all(|leaf| leaf.bernstein.len() == 16
        && leaf.bernstein.iter().all(|row| row.len() == leaf.dofs.len())
        && leaf.bernstein.iter().all(|row| (row.iter().sum::<f64>() - 1.0).abs() < 1e-12)));
}
