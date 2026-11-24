use solidcraft_trivariate::{
    bracket::{bracket_blocks, BracketBlock},
    incidence::{enumerate_faces, BlockSide, Orientation},
};

#[derive(Debug)]
struct ReferenceBlock {
    id: String,
    degrees: [usize; 3],
    controls: Vec<[i128; 4]>,
}

fn reference_blocks() -> Vec<ReferenceBlock> {
    // Frozen output of the independent Python E7 prototype. Controls are
    // blockwise integerized; comparisons below permit one common scale only.
    include_str!("references/bracket_controls.txt")
        .lines()
        .map(|line| {
            let mut fields = line.split('|');
            let id = fields.next().expect("block id").to_owned();
            let degrees: Vec<usize> = fields
                .next()
                .expect("degrees")
                .split(',')
                .map(|value| value.parse().expect("integer degree"))
                .collect();
            let degrees: [usize; 3] = degrees.try_into().expect("three degrees");
            let controls = fields
                .next()
                .expect("control net")
                .split(';')
                .map(|point| {
                    let values: Vec<i128> = point
                        .split(',')
                        .map(|value| value.parse().expect("integer coefficient"))
                        .collect();
                    values.try_into().expect("homogeneous quadruple")
                })
                .collect();
            assert!(fields.next().is_none(), "unexpected fixture field");
            ReferenceBlock { id, degrees, controls }
        })
        .collect()
}

fn positively_scaled_equal(actual: &[[i64; 4]], expected: &[[i128; 4]]) -> bool {
    if actual.len() != expected.len() {
        return false;
    }
    let base = actual
        .iter()
        .zip(expected)
        .flat_map(|(a, e)| a.iter().zip(e))
        .find(|(_, e)| **e != 0);
    let Some((actual_base, expected_base)) = base else {
        return false;
    };
    let actual_base = i128::from(*actual_base);
    let expected_base = *expected_base;
    if actual_base * expected_base <= 0 {
        return false;
    }
    actual.iter().zip(expected).all(|(a, e)| {
        a.iter()
            .zip(e)
            .all(|(av, ev)| i128::from(*av) * expected_base == *ev * actual_base)
    })
}

#[derive(Debug)]
struct FaceTrace {
    degree: [usize; 2],
    controls: Vec<[i128; 4]>,
}

fn face_trace(block: &BracketBlock, side: &BlockSide) -> FaceTrace {
    let [p, q, r] = block.degrees;
    let mut points = Vec::new();
    for i in 0..=p {
        for j in 0..=q {
            for k in 0..=r {
                if [i, j, k][side.axis] == if side.upper { block.degrees[side.axis] } else { 0 } {
                    let index = (i * (q + 1) + j) * (r + 1) + k;
                    points.push(block.controls[index].map(i128::from));
                }
            }
        }
    }
    let degree = match side.axis {
        0 => [q, r],
        1 => [p, r],
        _ => [p, q],
    };
    FaceTrace { degree, controls: points }
}

fn elevate_trace(mut trace: FaceTrace) -> FaceTrace {
    for axis in 0..2 {
        if trace.degree[axis] == 2 {
            continue;
        }
        let old_degree = trace.degree;
        let mut new_degree = old_degree;
        new_degree[axis] = 2;
        let point = |i: usize, j: usize| trace.controls[i * (old_degree[1] + 1) + j];
        let mut controls = Vec::new();
        for i in 0..=new_degree[0] {
            for j in 0..=new_degree[1] {
                let (index, first, second) = if axis == 0 {
                    (i, point(0, j), point(1, j))
                } else {
                    (j, point(i, 0), point(i, 1))
                };
                controls.push(std::array::from_fn(|component| match index {
                    0 => 2 * first[component],
                    1 => first[component] + second[component],
                    _ => 2 * second[component],
                }));
            }
        }
        trace = FaceTrace { degree: new_degree, controls };
    }
    trace
}

fn orient_trace(trace: &FaceTrace, orientation: Orientation) -> FaceTrace {
    let [p, q] = trace.degree;
    let degree = if orientation.swap { [q, p] } else { [p, q] };
    let mut controls = Vec::new();
    for i in 0..=degree[0] {
        for j in 0..=degree[1] {
            let a = if orientation.flip_first { degree[0] - i } else { i };
            let b = if orientation.flip_second { degree[1] - j } else { j };
            let (source_i, source_j) = if orientation.swap { (b, a) } else { (a, b) };
            controls.push(trace.controls[source_i * (q + 1) + source_j]);
        }
    }
    FaceTrace { degree, controls }
}

fn matches_scaled(actual: &[[i128; 4]], expected: &[[i128; 4]]) -> bool {
    if actual.len() != expected.len() {
        return false;
    }
    let base = actual.iter().zip(expected).flat_map(|(a, e)| a.iter().zip(e))
        .find(|(_, e)| **e != 0);
    let Some((actual_base, expected_base)) = base else {
        return false;
    };
    if actual_base * expected_base <= 0 {
        return false;
    }
    actual.iter().zip(expected).all(|(a, e)| {
        a.iter().zip(e).all(|(av, ev)| av * expected_base == ev * actual_base)
    })
}

#[test]
fn rust_blocks_match_frozen_python_prototype_by_one_scale_per_block() {
    let actual = bracket_blocks().expect("fixed bracket construction");
    let expected = reference_blocks();
    assert_eq!(actual.len(), 36);
    assert_eq!(expected.len(), 36);
    for reference in &expected {
        let block = actual.iter().find(|block| block.id == reference.id)
            .unwrap_or_else(|| panic!("missing Rust block {}", reference.id));
        assert_eq!(block.degrees, reference.degrees, "{}", block.id);
        assert!(positively_scaled_equal(&block.controls, &reference.controls),
                "{} differs beyond a common positive homogeneous scale", block.id);
    }
}

#[test]
fn rust_bracket_has_expected_whole_face_incidence_and_orientations() {
    let blocks = bracket_blocks().unwrap();
    let incidence = enumerate_faces(&blocks).expect("manifold whole-face partition");
    assert_eq!(incidence.interfaces.len(), 55);
    assert_eq!(incidence.exterior.len(), 106);
    for interface in &incidence.interfaces {
        let first_block = blocks.iter().find(|block| block.id == interface.first.block).unwrap();
        let second_block = blocks.iter().find(|block| block.id == interface.second.block).unwrap();
        let first = elevate_trace(face_trace(first_block, &interface.first));
        let second = elevate_trace(face_trace(second_block, &interface.second));
        let oriented_second = orient_trace(&second, interface.orientation);
        assert_eq!(first.degree, oriented_second.degree);
        assert!(matches_scaled(&first.controls, &oriented_second.controls),
                "trace mismatch: {:?} / {:?}", interface.first, interface.second);
    }
}
