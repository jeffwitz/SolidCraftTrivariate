use solidcraft_trivariate::affine_overlap::{
    audit_declared_contact, classify_pair, AffineParallelepiped, ContactDecision, ContactFeature,
    PairRelation,
};

fn box_block(id: &str, lo: [i64; 3], hi: [i64; 3], id_base: u64) -> AffineParallelepiped {
    let edges = [
        [hi[0] - lo[0], 0, 0],
        [0, hi[1] - lo[1], 0],
        [0, 0, hi[2] - lo[2]],
    ];
    let vertices_um = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            lo[j]
                + ((i & 1) as i64) * edges[0][j]
                + (((i >> 1) & 1) as i64) * edges[1][j]
                + (((i >> 2) & 1) as i64) * edges[2][j]
        })
    });
    AffineParallelepiped {
        id: id.into(),
        vertices_um,
        vertex_ids: std::array::from_fn(|i| id_base + i as u64),
    }
}

fn box_relation(a: ([i64; 3], [i64; 3]), b: ([i64; 3], [i64; 3])) -> PairRelation {
    let mut touching = 0;
    let mut overlap_dimensions = 0;
    for axis in 0..3 {
        let lo = a.0[axis].max(b.0[axis]);
        let hi = a.1[axis].min(b.1[axis]);
        if lo > hi {
            return PairRelation::StrictlySeparated;
        }
        if lo == hi {
            touching += 1;
        } else {
            overlap_dimensions += 1;
        }
    }
    if touching > 0 {
        return PairRelation::BoundaryContact {
            dimension: overlap_dimensions,
        };
    }
    let a_contains_b = (0..3).all(|i| a.0[i] <= b.0[i] && a.1[i] >= b.1[i]);
    let b_contains_a = (0..3).all(|i| b.0[i] <= a.0[i] && b.1[i] >= a.1[i]);
    match (a_contains_b, b_contains_a) {
        (true, true) => PairRelation::Coincident,
        (true, false) => PairRelation::FirstContainsSecond,
        (false, true) => PairRelation::SecondContainsFirst,
        (false, false) => PairRelation::InteriorOverlap,
    }
}

#[test]
fn exhaustive_small_axis_boxes_match_independent_integer_interval_oracle() {
    let mut boxes = Vec::new();
    for x in -1..=0 {
        for y in -1..=0 {
            for z in -1..=0 {
                for dx in 1..=2 {
                    for dy in 1..=2 {
                        for dz in 1..=2 {
                            boxes.push(([x, y, z], [x + dx, y + dy, z + dz]));
                        }
                    }
                }
            }
        }
    }
    assert_eq!(boxes.len(), 64);
    for (i, a) in boxes.iter().enumerate() {
        for (j, b) in boxes.iter().enumerate() {
            let expected = box_relation(*a, *b);
            let aa = box_block("a", a.0, a.1, (i * 16) as u64);
            let bb = box_block("b", b.0, b.1, (j * 16) as u64 + 10_000);
            assert_eq!(
                classify_pair(&aa, &bb).unwrap(),
                expected,
                "pair {i},{j}: {a:?} / {b:?}"
            );
        }
    }
}

fn transform_common(block: &AffineParallelepiped) -> AffineParallelepiped {
    // Integer unimodular shear (x,y,z) -> (x+y,y,z), followed by translation.
    let vertices_um = block
        .vertices_um
        .map(|p| [p[0] + p[1] + 17, p[1] - 9, p[2] + 5]);
    AffineParallelepiped {
        id: format!("{}-sheared", block.id),
        vertices_um,
        vertex_ids: block.vertex_ids,
    }
}

#[test]
fn common_unimodular_shear_preserves_pair_relation() {
    let a = box_block("a", [0, 0, 0], [1000, 1000, 1000], 0);
    let b = box_block("b", [1000, 0, 0], [2000, 1000, 1000], 20);
    let before = classify_pair(&a, &b).unwrap();
    assert_eq!(before, PairRelation::BoundaryContact { dimension: 2 });
    let after = classify_pair(&transform_common(&a), &transform_common(&b)).unwrap();
    assert_eq!(after, before);
}

fn cross(a: [i128; 3], b: [i128; 3]) -> [i128; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [i128; 3], b: [i128; 3]) -> i128 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn projections(block: &AffineParallelepiped, axis: [i128; 3]) -> (i128, i128) {
    let values = block
        .vertices_um
        .iter()
        .map(|p| dot(p.map(i128::from), axis))
        .collect::<Vec<_>>();
    (*values.iter().min().unwrap(), *values.iter().max().unwrap())
}

#[test]
fn sole_edge_cross_axis_separates_oblique_blocks() {
    // Fixed independently enumerated integer fixture. Its 6 face axes all
    // overlap; exactly one of the 9 edge-pair axes separates the solids.
    let a = AffineParallelepiped {
        id: "a".into(),
        vertices_um: [
            [0, 0, 0],
            [2, 0, 0],
            [0, 2, 0],
            [2, 2, 0],
            [0, 0, 2],
            [2, 0, 2],
            [0, 2, 2],
            [2, 2, 2],
        ],
        vertex_ids: std::array::from_fn(|i| i as u64),
    };
    let b_edges = [[-2, -2, -2], [-2, -2, -1], [-2, -1, -2]];
    let b = AffineParallelepiped {
        id: "b".into(),
        vertices_um: std::array::from_fn(|i| {
            std::array::from_fn(|k| {
                [1, 0, 3][k]
                    + ((i & 1) as i64) * b_edges[0][k]
                    + (((i >> 1) & 1) as i64) * b_edges[1][k]
                    + (((i >> 2) & 1) as i64) * b_edges[2][k]
            })
        }),
        vertex_ids: std::array::from_fn(|i| 20 + i as u64),
    };
    let ea = [[2_i128, 0, 0], [0, 2, 0], [0, 0, 2]];
    let eb = b_edges.map(|e| e.map(i128::from));
    let face_axes = [
        cross(ea[1], ea[2]),
        cross(ea[2], ea[0]),
        cross(ea[0], ea[1]),
        cross(eb[1], eb[2]),
        cross(eb[2], eb[0]),
        cross(eb[0], eb[1]),
    ];
    let edge_axes = std::array::from_fn::<_, 9, _>(|i| cross(ea[i / 3], eb[i % 3]));
    assert!(face_axes.iter().all(|axis| {
        let pa = projections(&a, *axis);
        let pb = projections(&b, *axis);
        pa.1 >= pb.0 && pb.1 >= pa.0
    }));
    let separating_edge_axes = edge_axes
        .iter()
        .filter(|axis| {
            let pa = projections(&a, **axis);
            let pb = projections(&b, **axis);
            pa.1 < pb.0 || pb.1 < pa.0
        })
        .count();
    assert_eq!(separating_edge_axes, 1);
    assert_eq!(
        classify_pair(&a, &b).unwrap(),
        PairRelation::StrictlySeparated
    );
}

#[test]
fn omitted_and_partial_contacts_are_not_accepted_as_declared_features() {
    let a = box_block("a", [0, 0, 0], [1000, 1000, 1000], 0);
    let mut b = box_block("b", [1000, 0, 0], [2000, 1000, 1000], 100);
    assert_eq!(
        classify_pair(&a, &b).unwrap(),
        PairRelation::BoundaryContact { dimension: 2 }
    );
    assert_eq!(
        audit_declared_contact(&a, &b, &ContactFeature::Face([1, 3, 5, 7]))
            .unwrap()
            .decision,
        ContactDecision::Rejected,
        "a geometric contact without shared vertex identities is omitted"
    );

    // Partial inclusion and overlap are classified as volumetric, not as an
    // accepted face interface.
    b = box_block("b", [500, 0, 0], [1500, 1000, 1000], 100);
    assert_eq!(
        classify_pair(&a, &b).unwrap(),
        PairRelation::InteriorOverlap
    );
    let inner = box_block("inner", [100, 100, 100], [900, 900, 900], 200);
    assert_eq!(
        classify_pair(&a, &inner).unwrap(),
        PairRelation::FirstContainsSecond
    );
}
