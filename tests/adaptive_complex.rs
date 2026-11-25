use solidcraft_trivariate::{
    adaptive_complex::build_oriented_complex,
    adaptive_partition::{AdaptivePartition, LeafKey},
    adaptive_space::build_adaptive_space,
    analysis_space::{AnalysisBlock, AnalysisSpace},
    bracket::bracket_blocks,
    bracket_audit::audit_bracket,
    incidence::{BlockSide, Incidence, Interface, Orientation},
};
use std::collections::BTreeMap;

fn bracket_case(refine_index: usize, depth: u8) {
    let blocks = bracket_blocks().unwrap();
    let certificate = audit_bracket(&blocks).unwrap();
    let parent = solidcraft_trivariate::analysis_space::build_analysis_space(&blocks).unwrap();
    let roots = parent
        .blocks
        .iter()
        .map(|block| block.id.clone())
        .collect::<Vec<_>>();
    let mut partition =
        AdaptivePartition::new(&roots, &certificate.incidence, depth, 5000).unwrap();
    if depth > 0 {
        let leaf = LeafKey {
            root: roots[refine_index].clone(),
            level: 0,
            index: [0; 3],
        };
        partition.refine(&[leaf]).unwrap();
    }
    if depth > 1 {
        let root = roots[refine_index].clone();
        let leaves = partition
            .leaves()
            .iter()
            .filter(|leaf| leaf.root == root)
            .cloned()
            .collect::<Vec<_>>();
        partition.refine(&leaves).unwrap();
    }
    let adaptive = build_adaptive_space(&parent, &partition).unwrap();
    let complex = build_oriented_complex(&adaptive, &partition).unwrap();
    let euler = complex.vertices.len() as isize - complex.edges.len() as isize
        + complex.faces.len() as isize
        - complex.cells.len() as isize;
    assert_eq!(
        euler,
        -1,
        "refinement case {refine_index}, depth {depth}: V={}, E={}, F={}, C={}",
        complex.vertices.len(),
        complex.edges.len(),
        complex.faces.len(),
        complex.cells.len()
    );
    assert_eq!(
        complex.boundary_faces.len(),
        complex
            .faces
            .iter()
            .filter(|face| face.incident_cells.len() == 1)
            .count()
    );
    assert!(complex
        .faces
        .iter()
        .all(|face| face.incident_cells.len() == 1 || face.incident_cells.len() == 2));
}

#[test]
fn six_bracket_refinement_cases_have_verified_oriented_boundary() {
    for (index, depth) in [(0, 0), (2, 1), (8, 1), (14, 2), (22, 1), (30, 2)] {
        bracket_case(index, depth);
    }
}

fn two_cubes() -> (AnalysisSpace, Incidence) {
    let mut blocks = Vec::new();
    for (id, second) in [("a", false), ("b", true)] {
        let mut controls = Vec::with_capacity(27);
        let mut dofs = Vec::with_capacity(27);
        for i in 0..3_i128 {
            for j in 0..3_i128 {
                for k in 0..3_i128 {
                    controls.push(if second {
                        [2 + i, k, 2 - j, 2]
                    } else {
                        [i, j, k, 2]
                    });
                    dofs.push(if second { 27 + dofs.len() } else { dofs.len() });
                }
            }
        }
        blocks.push(AnalysisBlock {
            id: id.to_owned(),
            controls,
            dofs,
        });
    }
    let mut exterior = Vec::new();
    for block in ["a", "b"] {
        for axis in 0..3 {
            for upper in [false, true] {
                if block == "a" && axis == 0 && upper || block == "b" && axis == 0 && !upper {
                    continue;
                }
                exterior.push(BlockSide {
                    block: block.to_owned(),
                    axis,
                    upper,
                });
            }
        }
    }
    (
        AnalysisSpace {
            ndof: 54,
            blocks,
            boundary_dofs: BTreeMap::new(),
        },
        Incidence {
            interfaces: vec![Interface {
                first: BlockSide {
                    block: "a".into(),
                    axis: 0,
                    upper: true,
                },
                second: BlockSide {
                    block: "b".into(),
                    axis: 0,
                    upper: false,
                },
                orientation: Orientation {
                    swap: true,
                    flip_first: false,
                    flip_second: true,
                },
            }],
            exterior,
        },
    )
}

#[test]
fn d4_two_cube_junction_keeps_hanging_boundary_subdivision() {
    let (parent, incidence) = two_cubes();
    let mut partition =
        AdaptivePartition::new(&["a".into(), "b".into()], &incidence, 2, 100).unwrap();
    let root = LeafKey {
        root: "a".into(),
        level: 0,
        index: [0; 3],
    };
    partition.refine(&[root]).unwrap();
    let adaptive = build_adaptive_space(&parent, &partition).unwrap();
    let complex = build_oriented_complex(&adaptive, &partition).unwrap();
    assert_eq!(complex.cells.len(), 9);
    assert_eq!(
        complex.vertices.len() as isize - complex.edges.len() as isize
            + complex.faces.len() as isize
            - complex.cells.len() as isize,
        1
    );
    assert!(complex.boundary_faces.iter().any(|&face_id| {
        let face = &complex.faces[face_id];
        matches!(&face.source, solidcraft_trivariate::adaptive_complex::FaceSource::Exterior { side, .. }
            if side.root_side.axis == 0 && side.root_side.block == "b")
    }));
}

#[test]
fn optionally_write_d4_complex_for_independent_oracle() {
    let Ok(path) = std::env::var("ADAPTIVE_COMPLEX_D4_OUTPUT") else {
        return;
    };
    let (parent, incidence) = two_cubes();
    let mut partition =
        AdaptivePartition::new(&["a".into(), "b".into()], &incidence, 2, 100).unwrap();
    partition
        .refine(&[LeafKey {
            root: "a".into(),
            level: 0,
            index: [0; 3],
        }])
        .unwrap();
    let adaptive = build_adaptive_space(&parent, &partition).unwrap();
    let complex = build_oriented_complex(&adaptive, &partition).unwrap();
    std::fs::write(path, complex.to_json().to_string()).unwrap();
}
