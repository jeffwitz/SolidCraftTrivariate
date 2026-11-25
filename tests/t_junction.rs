//! Three-dimensional Q1 fixtures with separate local numbering in each brick.
use solidcraft_trivariate::constraint_system::{resolve, Lift, Ratio};
use solidcraft_trivariate::face_constraints::{FaceNet, Row, Span};
use solidcraft_trivariate::face_partition::{
    partition_master, FacePatch, FaceUse, MasterSide, PartitionUse, TopologicalFaceId,
};
use solidcraft_trivariate::incidence::{BlockSide, Orientation};

#[derive(Clone, Copy)]
struct Brick {
    id: &'static str,
    start: usize,
    x: [i64; 2],
    y: [i64; 2],
}

impl Brick {
    fn point(self, local: usize) -> [i64; 3] {
        [
            self.x[local / 4],
            self.y[(local / 2) % 2],
            (local % 2) as i64,
        ]
    }

    fn face(self, axis: usize, upper: bool) -> ([usize; 4], FaceNet, BlockSide) {
        let local: Vec<_> = (0..8)
            .filter(|i| ((i >> (2 - axis)) & 1) == usize::from(upper))
            .collect();
        let controls = local
            .iter()
            .map(|&i| {
                let [x, y, z] = self.point(i);
                [x, y, z, 1]
            })
            .collect();
        (
            std::array::from_fn(|i| self.start + local[i]),
            FaceNet {
                degrees: [1, 1],
                controls,
            },
            BlockSide {
                block: self.id.into(),
                axis,
                upper,
            },
        )
    }
}

fn attach(
    equations: &mut Vec<(usize, Row)>,
    master: Brick,
    axis: usize,
    neighbors: &[(Brick, Span)],
) {
    let (ids, net, side) = master.face(axis, true);
    let master_side = MasterSide::validate(side, net).unwrap();
    let mut patches = vec![];
    for (child, span) in neighbors {
        let (_, net, slave) = child.face(axis, false);
        patches.push(FacePatch {
            id: TopologicalFaceId::new(child.id).unwrap(),
            spans: [*span, Span::Full],
            usage: FaceUse::Interface {
                slave,
                net,
                orientation: Orientation {
                    swap: false,
                    flip_first: false,
                    flip_second: false,
                },
            },
        });
    }
    let partition = partition_master(&master_side, &patches).unwrap();
    for c in partition.cells {
        let PartitionUse::Interface {
            slave, constraints, ..
        } = c.usage
        else {
            panic!()
        };
        let neighbor = neighbors
            .iter()
            .find(|(b, _)| b.id == slave.block)
            .unwrap()
            .0;
        let (slave_ids, _, _) = neighbor.face(axis, false);
        for (dof, row) in slave_ids.into_iter().zip(constraints) {
            equations.push((
                dof,
                Row {
                    terms: row.terms.into_iter().map(|(i, w)| (ids[i], w)).collect(),
                    denominator: row.denominator,
                },
            ));
        }
    }
}

fn affine([x, y, z]: [i64; 3]) -> i64 {
    1 + 3 * x - 5 * y + 7 * z
}

fn verify_field(lift: &Lift, bricks: &[Brick]) {
    let prescribed: Vec<_> = bricks
        .iter()
        .flat_map(|b| (0..8).map(|i| affine(b.point(i)) as f64))
        .collect();
    let value = |r: Ratio| r.numerator() as f64 / r.denominator() as f64;
    for line in &lift.lines {
        let reconstructed = line.terms.iter().fold(value(line.offset), |s, (i, r)| {
            s + value(*r) * prescribed[lift.free_roots[*i]]
        });
        assert_eq!(reconstructed, prescribed[line.dof]);
    }
}

#[test]
fn three_bricks_resolve_duplicate_edge_constraints_and_reproduce_an_affine_field() {
    let master = Brick {
        id: "west",
        start: 0,
        x: [0, 1],
        y: [0, 2],
    };
    let low = Brick {
        id: "southeast",
        start: 8,
        x: [1, 2],
        y: [0, 1],
    };
    let high = Brick {
        id: "northeast",
        start: 16,
        x: [1, 2],
        y: [1, 2],
    };
    let mut eq = vec![];
    attach(
        &mut eq,
        master,
        0,
        &[(low, Span::Lower), (high, Span::Upper)],
    );
    attach(&mut eq, low, 1, &[(high, Span::Full)]);
    assert_eq!(eq.len(), 12);
    let lift = resolve(24, &eq, &[]).unwrap();
    // Two equations on the common hanging edge are redundant.
    assert_eq!(lift.free_roots.len(), 14);
    verify_field(&lift, &[master, low, high]);
    eq.reverse();
    assert_eq!(lift, resolve(24, &eq, &[]).unwrap());
}

#[test]
fn splitting_master_analysis_support_leaves_lower_interface_free_under_partial_dirichlet() {
    let low = Brick {
        id: "southwest",
        start: 0,
        x: [0, 1],
        y: [0, 1],
    };
    let high = Brick {
        id: "northwest",
        start: 8,
        x: [0, 1],
        y: [1, 2],
    };
    let east = Brick {
        id: "southeast",
        start: 16,
        x: [1, 2],
        y: [0, 1],
    };
    let mut eq = vec![];
    attach(&mut eq, low, 1, &[(high, Span::Full)]);
    attach(&mut eq, low, 0, &[(east, Span::Full)]);
    let unrestricted = resolve(24, &eq, &[]).unwrap();
    assert_eq!(unrestricted.free_roots.len(), 16);
    let boundary = high.face(0, true).0;
    // Prescribe independent representatives of the four upper-face values.
    let mut fixed = vec![];
    for i in boundary {
        let line = &unrestricted.lines[i];
        let [(index, weight)] = line.terms[..] else {
            panic!()
        };
        assert_eq!(weight, Ratio::one());
        fixed.push((
            unrestricted.free_roots[index],
            Ratio::new(affine(high.point(i - high.start)).into(), 1).unwrap(),
        ));
    }
    let lift = resolve(24, &eq, &fixed).unwrap();
    assert_eq!(lift.free_roots.len(), 12);
    verify_field(&lift, &[low, high, east]);
    for i in boundary {
        assert!(lift.lines[i].terms.is_empty());
    }
    for i in &low.face(0, true).0[..2] {
        assert_eq!(lift.lines[*i].terms.len(), 1);
    }
    fixed.push((boundary[0], Ratio::new(1234, 1).unwrap()));
    assert!(resolve(24, &eq, &fixed).is_err());
}

#[test]
fn long_constraint_chains_use_a_bounded_call_stack() {
    let size = 20_000;
    let eq: Vec<_> = (0..size - 1)
        .map(|i| {
            (
                i,
                Row {
                    terms: vec![(i + 1, 1)],
                    denominator: 1,
                },
            )
        })
        .collect();
    let reduced = resolve(size, &eq, &[]).unwrap();
    assert_eq!(reduced.free_roots, vec![size - 1]);
    assert_eq!(reduced.lines[0].terms, vec![(0, Ratio::one())]);
}
