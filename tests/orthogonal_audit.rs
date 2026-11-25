use solidcraft_trivariate::{
    orthogonal_audit::audit_orthogonal_program,
    orthogonal_program::{build_orthogonal_program, BoxFeature, BoxOp, OrthogonalProgram},
};

fn feature(id: &str, operation: BoxOp, bounds_mm: [[i32; 2]; 3]) -> BoxFeature {
    BoxFeature {
        id: id.into(),
        operation,
        bounds_mm,
    }
}
fn program(base_mm: [[i32; 2]; 3], features: Vec<BoxFeature>) -> OrthogonalProgram {
    OrthogonalProgram { base_mm, features }
}
fn build_and_audit(
    p: &OrthogonalProgram,
) -> solidcraft_trivariate::orthogonal_audit::OrthogonalCertificate {
    let blocks = build_orthogonal_program(p).unwrap();
    audit_orthogonal_program(p, &blocks).unwrap()
}

#[test]
fn orthogonal_program_corpus_has_exact_volume_topology_and_provenance() {
    let cube = program([[0, 4]; 3], vec![]);
    let a = build_and_audit(&cube);
    assert_eq!((a.volume_mm3, a.cavities, a.minimum_x_mm), (64, 0, 0));
    assert_eq!(a.certificate.incidence.exterior.len(), 6);
    assert_eq!(a.certificate.design_faces.len(), 6);
    assert!(a
        .certificate
        .design_faces
        .iter()
        .all(|(_, source)| source == "base"));
    assert_eq!(
        a.block_provenance,
        vec![("box_0_4_0_4_0_4".into(), "base".into())]
    );

    let cavity = program([[0, 4]; 3], vec![feature("bore", BoxOp::Cut, [[1, 3]; 3])]);
    let c = build_and_audit(&cavity);
    assert_eq!((c.volume_mm3, c.cavities), (56, 1));

    let tunnel = program(
        [[0, 4]; 3],
        vec![feature("tunnel", BoxOp::Cut, [[0, 4], [1, 3], [1, 3]])],
    );
    let t = build_and_audit(&tunnel);
    assert_eq!((t.volume_mm3, t.cavities), (48, 0));

    let crossing = program(
        [[0, 4]; 3],
        vec![
            feature("tunnel_x", BoxOp::Cut, [[0, 4], [1, 3], [1, 3]]),
            feature("tunnel_y", BoxOp::Cut, [[1, 3], [0, 4], [1, 3]]),
        ],
    );
    let x = build_and_audit(&crossing);
    assert_eq!((x.volume_mm3, x.cavities), (40, 0));

    let two_cavities = program(
        [[0, 6]; 3],
        vec![
            feature("cavity_a", BoxOp::Cut, [[1, 2]; 3]),
            feature("cavity_b", BoxOp::Cut, [[4, 5]; 3]),
        ],
    );
    let two = build_and_audit(&two_cavities);
    assert_eq!((two.volume_mm3, two.cavities), (214, 2));

    let face_add = program(
        [[0, 4]; 3],
        vec![feature("face_pad", BoxOp::Add, [[4, 6], [1, 3], [1, 3]])],
    );
    let add = build_and_audit(&face_add);
    assert_eq!((add.volume_mm3, add.cavities, add.minimum_x_mm), (72, 0, 0));

    let refill = program(
        [[0, 6]; 3],
        vec![
            feature("rough", BoxOp::Cut, [[1, 5]; 3]),
            feature("refill", BoxOp::Add, [[1, 3], [2, 4], [2, 4]]),
        ],
    );
    let fill = build_and_audit(&refill);
    assert_eq!((fill.volume_mm3, fill.cavities), (160, 1));
    assert!(fill
        .block_provenance
        .iter()
        .any(|(_, source)| source == "refill"));
    let face_sources = fill
        .certificate
        .design_faces
        .iter()
        .map(|(_, source)| source.as_str())
        .collect::<std::collections::BTreeSet<_>>();
    assert!(face_sources.contains("base"));
    assert!(face_sources.contains("rough"));
    assert!(face_sources.contains("refill"));

    let redundant = program(
        [[0, 4]; 3],
        vec![feature("redundant", BoxOp::Add, [[1, 2]; 3])],
    );
    let red = build_and_audit(&redundant);
    assert_eq!(red.volume_mm3, 64);
    assert!(red
        .block_provenance
        .iter()
        .all(|(_, source)| source == "base"));
}

#[test]
fn local_link_and_material_refusals_are_independently_reproduced() {
    let base = [[0, 4]; 3];
    for (program, reason) in [
        (
            program(
                base,
                vec![feature("edge", BoxOp::Add, [[4, 6], [4, 6], [1, 3]])],
            ),
            "orthogonal_nonmanifold_edge",
        ),
        (
            program(base, vec![feature("vertex", BoxOp::Add, [[4, 6]; 3])]),
            "orthogonal_nonmanifold_vertex",
        ),
        (
            program(base, vec![feature("floating", BoxOp::Add, [[5, 6]; 3])]),
            "orthogonal_disconnected_material",
        ),
        (
            program(base, vec![feature("empty", BoxOp::Cut, base)]),
            "orthogonal_empty_material",
        ),
    ] {
        assert_eq!(audit_orthogonal_program(&program, &[]), Err(reason));
    }
}

#[test]
fn malformed_program_and_altered_or_duplicated_blocks_are_refused() {
    let dup = program(
        [[0, 4]; 3],
        vec![
            feature("repeat", BoxOp::Add, [[1, 2]; 3]),
            feature("repeat", BoxOp::Cut, [[2, 3]; 3]),
        ],
    );
    assert_eq!(
        audit_orthogonal_program(&dup, &[]),
        Err("orthogonal_feature_id")
    );
    let flat = program([[0, 0], [0, 4], [0, 4]], vec![]);
    assert_eq!(
        audit_orthogonal_program(&flat, &[]),
        Err("orthogonal_program_range")
    );

    let p = program([[0, 4]; 3], vec![]);
    let mut blocks = build_orthogonal_program(&p).unwrap();
    blocks[0].controls[0][0] = 1;
    assert_eq!(
        audit_orthogonal_program(&p, &blocks),
        Err("orthogonal_altered_block")
    );
    let mut blocks = build_orthogonal_program(&p).unwrap();
    blocks.push(blocks[0].clone());
    assert_eq!(
        audit_orthogonal_program(&p, &blocks),
        Err("orthogonal_coverage_mismatch")
    );
}
