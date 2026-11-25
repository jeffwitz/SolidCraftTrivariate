use solidcraft_trivariate::{
    constraint_system::Ratio,
    sketch_arrangement::{arrange, FaceClass, Role, Segment},
};

fn edge(id: &str, a: [i32; 2], b: [i32; 2], role: Role) -> Segment {
    Segment {
        id: id.into(),
        a,
        b,
        role,
    }
}

fn square(prefix: &str, lo: i32, hi: i32) -> Vec<Segment> {
    let p = [[lo, lo], [hi, lo], [hi, hi], [lo, hi]];
    (0..4)
        .map(|i| {
            edge(
                &format!("{prefix}{i}"),
                p[i],
                p[(i + 1) % 4],
                Role::Boundary,
            )
        })
        .collect()
}

#[test]
fn construction_chords_split_a_closed_material_loop() {
    let mut segments = square("outer", 0, 10);
    segments.push(edge("divider", [0, 4], [10, 6], Role::Construction));
    let graph = arrange(&segments).unwrap();
    assert_eq!(graph.faces.len(), 2);
    assert!(graph.faces.iter().all(|f| f.class == FaceClass::Material));
    assert!(graph.faces.iter().all(|f| f.vertices.len() == 4));
    assert!(graph.faces.iter().all(|face| {
        let area = face
            .vertices
            .iter()
            .enumerate()
            .map(|(i, &a)| {
                let b = face.vertices[(i + 1) % face.vertices.len()];
                let p = graph.vertices[a];
                let q = graph.vertices[b];
                p[0].numerator() as f64 / p[0].denominator() as f64
                    * (q[1].numerator() as f64 / q[1].denominator() as f64)
                    - (p[1].numerator() as f64 / p[1].denominator() as f64)
                        * (q[0].numerator() as f64 / q[0].denominator() as f64)
            })
            .sum::<f64>();
        area > 0.
    }));
}

#[test]
fn nested_boundaries_classify_material_void_and_exterior_cycles() {
    let mut segments = square("outer", 0, 20);
    segments.extend(square("inner", 5, 15));
    let graph = arrange(&segments).unwrap();
    assert_eq!(graph.faces.len(), 2);
    assert!(graph.faces.iter().any(|f| f.class == FaceClass::Material));
    assert!(graph.faces.iter().any(|f| f.class == FaceClass::Void));
}

#[test]
fn crossing_construction_lines_create_exact_rational_nodes_and_faces() {
    let mut segments = square("outer", 0, 6);
    segments.push(edge("diag_a", [0, 1], [6, 4], Role::Construction));
    segments.push(edge("diag_b", [1, 0], [5, 6], Role::Construction));
    let graph = arrange(&segments).unwrap();
    assert_eq!(graph.faces.len(), 4);
    assert!(graph.faces.iter().all(|f| f.class == FaceClass::Material));
    assert!(graph
        .vertices
        .contains(&[Ratio::new(5, 2).unwrap(), Ratio::new(9, 4).unwrap(),]));
}

#[test]
fn coincident_partial_edges_are_split_once_with_all_provenance() {
    let mut segments = square("outer", 0, 10);
    segments.push(edge(
        "construction_overlap",
        [2, 0],
        [8, 0],
        Role::Construction,
    ));
    let graph = arrange(&segments).unwrap();
    let shared = graph
        .edges
        .iter()
        .find(|e| {
            let a = graph.vertices[e.vertices[0]];
            let b = graph.vertices[e.vertices[1]];
            a[1] == Ratio::zero()
                && b[1] == Ratio::zero()
                && a[0] >= Ratio::new(2, 1).unwrap()
                && b[0] <= Ratio::new(8, 1).unwrap()
        })
        .unwrap();
    assert_eq!(shared.provenance.len(), 2);
    assert!(shared.provenance.iter().any(|p| p.role == Role::Boundary));
    assert!(shared
        .provenance
        .iter()
        .any(|p| p.role == Role::Construction));
}

#[test]
fn invalid_inputs_are_rejected() {
    assert!(arrange(&[]).is_err());
    assert!(arrange(&[edge("x", [0, 0], [0, 0], Role::Boundary)]).is_err());
    assert!(arrange(&[edge("x", [0, 0], [1, 0], Role::Boundary)]).is_err());
    let mut duplicate = square("e", 0, 1);
    duplicate[1].id = duplicate[0].id.clone();
    assert!(arrange(&duplicate).is_err());
}
