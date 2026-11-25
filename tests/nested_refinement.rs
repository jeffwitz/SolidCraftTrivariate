#[path = "../examples/support/t_mesh.rs"]
mod t_mesh;

use solidcraft_trivariate::analysis_space::AnalysisSpace;

type RefinementScenario = (
    &'static str,
    &'static [&'static str],
    usize,
    usize,
    usize,
    usize,
);
use solidcraft_trivariate::constraint_system::{resolve, Lift};

fn field([x, y, z]: [f64; 3], component: usize) -> f64 {
    let rows = [[1.0, 2.0, -3.0], [-2.0, 1.0, 4.0], [3.0, -5.0, 2.0]];
    let offsets = [2.0, -1.0, 3.0];
    offsets[component] + rows[component][0] * x + rows[component][1] * y + rows[component][2] * z
}

fn coordinate(space: &AnalysisSpace, dof: usize) -> [f64; 3] {
    let p = space
        .blocks
        .iter()
        .find_map(|block| {
            block
                .dofs
                .iter()
                .position(|&value| value == dof)
                .map(|local| block.controls[local])
        })
        .expect("every global DOF has a control");
    std::array::from_fn(|axis| p[axis] as f64 / p[3] as f64)
}

fn assert_field_reproduction(space: &AnalysisSpace, lift: &Lift, component: usize) {
    let roots: Vec<_> = lift
        .free_roots
        .iter()
        .map(|&dof| field(coordinate(space, dof), component))
        .collect();
    for line in &lift.lines {
        let reconstructed = line.terms.iter().fold(
            line.offset.numerator() as f64 / line.offset.denominator() as f64,
            |value, (root, coefficient)| {
                value
                    + coefficient.numerator() as f64 / coefficient.denominator() as f64
                        * roots[*root]
            },
        );
        let expected = field(coordinate(space, line.dof), component);
        assert!(
            (reconstructed - expected).abs() <= 2e-12,
            "component={component} dof={} reconstructed={reconstructed} expected={expected}",
            line.dof
        );
    }
}

#[test]
fn staged_master_slave_and_neighbor_refinements_rebuild_all_traces() {
    let base = t_mesh::base_three_cell().unwrap();
    assert_eq!(base.space.blocks.len(), 3);
    assert_eq!(base.relations.len(), 3);
    let scenarios: [RefinementScenario; 4] = [
        ("master", &["master"], 10, 17, 4, 145),
        (
            "low_slave_with_master",
            &["master", "neighbor_low"],
            17,
            34,
            10,
            217,
        ),
        (
            "high_slave_with_master",
            &["master", "neighbor_high"],
            17,
            34,
            10,
            217,
        ),
        (
            "two_neighbors_with_master",
            &["master", "neighbor_low", "neighbor_high"],
            24,
            48,
            8,
            305,
        ),
    ];
    for (name, selected, expected_blocks, expected_relations, expected_partial, expected_roots) in
        scenarios
    {
        let refined = t_mesh::replace_refined(&base, selected).unwrap();
        assert_eq!(refined.space.blocks.len(), expected_blocks, "{name}");
        assert_eq!(
            refined.equations.len(),
            refined.relations.len() * 9,
            "{name}"
        );
        assert_eq!(refined.relations.len(), expected_relations, "{name}");
        assert!(refined.relations.len() > base.relations.len(), "{name}");
        assert!(refined
            .relations
            .iter()
            .all(|relation| relation.master_side.0 == relation.slave_side.0));
        let partial_spans: Vec<_> = refined
            .relations
            .iter()
            .filter(|relation| {
                relation
                    .spans
                    .iter()
                    .any(|span| *span != solidcraft_trivariate::face_constraints::Span::Full)
            })
            .map(|relation| (&relation.master_id, &relation.slave_id, relation.spans))
            .collect();
        println!(
            "{name}: relations={}, partial={}",
            refined.relations.len(),
            partial_spans.len()
        );
        assert_eq!(partial_spans.len(), expected_partial, "{name}");
        for id in selected {
            assert_eq!(
                refined
                    .space
                    .blocks
                    .iter()
                    .filter(|block| block.id.starts_with(&format!("{id}/")))
                    .count(),
                8,
                "{name}:{id}"
            );
        }
        for block in &base.space.blocks {
            if !selected.contains(&block.id.as_str()) {
                assert!(
                    refined
                        .space
                        .blocks
                        .iter()
                        .any(|candidate| candidate.id == block.id),
                    "unaffected block {} disappeared",
                    block.id
                );
            }
        }
        let lift = resolve(refined.space.ndof, &refined.equations, &[]).unwrap();
        assert_eq!(lift.free_roots.len(), expected_roots, "{name}");
        assert!(
            lift.free_roots.len() < refined.space.ndof,
            "{name} did not constrain interfaces"
        );
        for component in 0..3 {
            assert_field_reproduction(&refined.space, &lift, component);
        }
        println!(
            "{name}: blocks={} ndof={} relations={} rows={} free_roots={}",
            refined.space.blocks.len(),
            refined.space.ndof,
            refined.relations.len(),
            refined.equations.len(),
            lift.free_roots.len()
        );
    }
}

#[test]
fn successive_uniform_refinements_reconstruct_every_tjunction_trace() {
    let base = t_mesh::base_three_cell().unwrap();
    let base_lift = resolve(base.space.ndof, &base.equations, &[]).unwrap();
    assert_eq!(base_lift.free_roots.len(), 57);
    let first_space = solidcraft_trivariate::h_refine::refine_uniform(&base.space).unwrap();
    assert_eq!(first_space.blocks.len(), 24);
    let first = t_mesh::rebuild(first_space).unwrap();
    let first_lift = resolve(first.space.ndof, &first.equations, &[]).unwrap();
    assert_eq!(first.relations.len(), 48);
    assert_eq!(first_lift.free_roots.len(), 305);
    for component in 0..3 {
        assert_field_reproduction(&first.space, &first_lift, component);
    }
    let second_space = solidcraft_trivariate::h_refine::refine_uniform(&first.space).unwrap();
    assert_eq!(second_space.blocks.len(), 192);
    let second = t_mesh::rebuild(second_space).unwrap();
    let second_lift = resolve(second.space.ndof, &second.equations, &[]).unwrap();
    assert_eq!(second.relations.len(), 480);
    assert_eq!(second_lift.free_roots.len(), 1953);
    for component in 0..3 {
        assert_field_reproduction(&second.space, &second_lift, component);
    }
    println!(
        "successive: base={} first={} second={} relations={}/{}/{} free_roots={}/{}/{}",
        base.space.blocks.len(),
        first.space.blocks.len(),
        second.space.blocks.len(),
        base.relations.len(),
        first.relations.len(),
        second.relations.len(),
        resolve(base.space.ndof, &base.equations, &[])
            .unwrap()
            .free_roots
            .len(),
        first_lift.free_roots.len(),
        second_lift.free_roots.len()
    );
}
