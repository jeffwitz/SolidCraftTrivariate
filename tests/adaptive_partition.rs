//! Independent unit-voxel oracle for the complete oriented face partition.
use solidcraft_trivariate::{
    adaptive_partition::{AdaptivePartition, LeafKey, LeafSide},
    bracket::bracket_blocks,
    bracket_audit::audit_bracket,
    incidence::{BlockSide, Incidence, Interface, Orientation},
};
use std::collections::{BTreeMap, BTreeSet};

type Point = ([u32; 2], LeafSide);
fn pair(a: Point, b: Point) -> [Point; 2] {
    if a < b {
        [a, b]
    } else {
        [b, a]
    }
}
fn trans([u, v]: [u32; 2], n: u32, o: Orientation) -> [u32; 2] {
    let x = if o.flip_first { n - 1 - u } else { u };
    let y = if o.flip_second { n - 1 - v } else { v };
    if o.swap {
        [y, x]
    } else {
        [x, y]
    }
}
fn leaf_side(leaf: &LeafKey, axis: usize, upper: bool) -> LeafSide {
    LeafSide {
        leaf: leaf.clone(),
        root_side: BlockSide {
            block: leaf.root.clone(),
            axis,
            upper,
        },
    }
}
fn verify(plan: &AdaptivePartition, inc: &Incidence) {
    let n = plan.grid_denominator();
    let mut cells = BTreeMap::new();
    let mut roots = BTreeSet::new();
    for leaf in plan.leaves() {
        roots.insert(&leaf.root);
        let w = n >> leaf.level;
        for i in 0..w {
            for j in 0..w {
                for k in 0..w {
                    let pos = std::array::from_fn(|a| leaf.index[a] * w + [i, j, k][a]);
                    assert!(cells.insert((&leaf.root, pos), leaf).is_none());
                }
            }
        }
    }
    assert_eq!(cells.len(), roots.len() * (n as usize).pow(3));
    let mut expected = BTreeSet::new();
    for (&(root, xyz), &a) in &cells {
        for d in 0..3 {
            let mut to = xyz;
            to[d] += 1;
            if to[d] == n {
                continue;
            }
            let b = cells[&(root, to)];
            if a != b {
                let uv: Vec<_> = (0..3).filter(|&i| i != d).map(|i| xyz[i]).collect();
                let uv = [uv[0], uv[1]];
                expected.insert(pair(
                    (uv, leaf_side(a, d, true)),
                    (uv, leaf_side(b, d, false)),
                ));
            }
        }
    }
    let at = |side: &BlockSide, uv: [u32; 2]| {
        let axes: Vec<_> = (0..3).filter(|&i| i != side.axis).collect();
        let mut xyz = [0; 3];
        xyz[side.axis] = if side.upper { n - 1 } else { 0 };
        xyz[axes[0]] = uv[0];
        xyz[axes[1]] = uv[1];
        let l = cells[&(&side.block, xyz)];
        (uv, leaf_side(l, side.axis, side.upper))
    };
    for f in &inc.interfaces {
        for u in 0..n {
            for v in 0..n {
                let uv = [u, v];
                assert!(expected.insert(pair(
                    at(&f.first, uv),
                    at(&f.second, trans(uv, n, f.orientation))
                )));
            }
        }
    }
    let mut actual = BTreeSet::new();
    for p in plan.interfaces().unwrap() {
        assert!(p.first.leaf.level.abs_diff(p.second.leaf.level) <= 1);
        assert_eq!(p.master_rect, p.first_subrect);
        let [[x0, x1], [y0, y1]] = p.first_subrect.ranges;
        let other = p.second_subrect.ranges;
        for u in x0..x1 {
            for v in y0..y1 {
                let uv = [u, v];
                let xy = trans(uv, n, p.orientation);
                assert!((0..2).all(|d| other[d][0] <= xy[d] && xy[d] < other[d][1]));
                assert!(actual.insert(pair((uv, p.first.clone()), (xy, p.second.clone()))));
            }
        }
    }
    assert_eq!(expected.len(), actual.len());
    assert_eq!(expected.difference(&actual).next(), None);
    let mut boundary = BTreeSet::new();
    for side in &inc.exterior {
        for u in 0..n {
            for v in 0..n {
                boundary.insert(at(side, [u, v]));
            }
        }
    }
    for p in plan.exterior_patches() {
        let [[x0, x1], [y0, y1]] = p.root_rect.ranges;
        for u in x0..x1 {
            for v in y0..y1 {
                assert!(boundary.remove(&([u, v], p.side.clone())));
            }
        }
    }
    assert!(boundary.is_empty());
}
#[test]
fn each_d4_transports_a_refined_slave_and_propagates_to_its_master() {
    let roots = ["left".to_string(), "right".to_string()];
    for op in 0..8 {
        let first = BlockSide {
            block: roots[0].clone(),
            axis: 0,
            upper: true,
        };
        let second = BlockSide {
            block: roots[1].clone(),
            axis: 2,
            upper: false,
        };
        let mut exterior = Vec::new();
        for root in &roots {
            for axis in 0..3 {
                for upper in [false, true] {
                    let s = BlockSide {
                        block: root.clone(),
                        axis,
                        upper,
                    };
                    if s != first && s != second {
                        exterior.push(s);
                    }
                }
            }
        }
        let incidence = Incidence {
            interfaces: vec![Interface {
                first,
                second,
                orientation: Orientation {
                    swap: op & 1 != 0,
                    flip_first: op & 2 != 0,
                    flip_second: op & 4 != 0,
                },
            }],
            exterior,
        };
        let mut p = AdaptivePartition::new(&roots, &incidence, 3, 500).unwrap();
        verify(&p, &incidence);
        p.refine(&[LeafKey::root("right")]).unwrap();
        verify(&p, &incidence);
        let r = p
            .refine(&[LeafKey {
                root: "right".into(),
                level: 1,
                index: [0, 1, 0],
            }])
            .unwrap();
        assert!(r.induced.contains(&LeafKey::root("left")));
        verify(&p, &incidence);
        let before = p.snapshot().unwrap();
        assert!(p.refine(&[LeafKey::root("right")]).is_err());
        assert_eq!(p.snapshot().unwrap(), before);
    }
}
#[test]
fn local_bracket_refinement_preserves_all_internal_and_boundary_microfaces() {
    let blocks = bracket_blocks().unwrap();
    let inc = audit_bracket(&blocks).unwrap().incidence;
    let names = blocks.iter().map(|b| b.id.clone()).collect::<Vec<_>>();
    let target = names
        .iter()
        .find(|n| n.starts_with("hole_x_sector_east"))
        .unwrap()
        .clone();
    let mut plan = AdaptivePartition::new(&names, &inc, 3, 1600).unwrap();
    plan.refine(&[LeafKey::root(&target)]).unwrap();
    let second = LeafKey {
        root: target.clone(),
        level: 1,
        index: [1, 0, 0],
    };
    plan.refine(&[second]).unwrap();
    verify(&plan, &inc);
    plan.refine(&[LeafKey {
        root: target,
        level: 2,
        index: [3, 0, 0],
    }])
    .unwrap();
    verify(&plan, &inc);
}
