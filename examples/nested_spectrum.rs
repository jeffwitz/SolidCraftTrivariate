//! Independent spectrum export for the three-cell Q2 T-junction.

mod support {
    pub mod t_mesh;
}

use solidcraft_trivariate::constraint_system::{resolve, Lift, Ratio};
use solidcraft_trivariate::elasticity::assemble;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn value(ratio: Ratio) -> f64 {
    ratio.numerator() as f64 / ratio.denominator() as f64
}

fn coordinate(mesh: &support::t_mesh::TMesh, dof: usize) -> [f64; 3] {
    let point = mesh
        .space
        .blocks
        .iter()
        .find_map(|block| {
            block
                .dofs
                .iter()
                .position(|&entry| entry == dof)
                .map(|local| block.controls[local])
        })
        .expect("each scalar node comes from a block control");
    std::array::from_fn(|axis| point[axis] as f64 / point[3] as f64)
}

fn rigid_velocity(point: [f64; 3], mode: usize) -> [f64; 3] {
    match mode {
        0 => [1.0, 0.0, 0.0],
        1 => [0.0, 1.0, 0.0],
        2 => [0.0, 0.0, 1.0],
        3 => [0.0, -point[2], point[1]],
        4 => [point[2], 0.0, -point[0]],
        _ => [-point[1], point[0], 0.0],
    }
}

fn verify_rigid_lifts(mesh: &support::t_mesh::TMesh, lifts: [&Lift; 3]) {
    for mode in 0..6 {
        for (component, lift) in lifts.iter().enumerate() {
            let roots: Vec<_> = lift
                .free_roots
                .iter()
                .map(|&dof| rigid_velocity(coordinate(mesh, dof), mode)[component])
                .collect();
            for line in &lift.lines {
                let reconstructed = line
                    .terms
                    .iter()
                    .fold(value(line.offset), |sum, (root, weight)| {
                        sum + value(*weight) * roots[*root]
                    });
                let expected = rigid_velocity(coordinate(mesh, line.dof), mode)[component];
                assert!(
                    (reconstructed - expected).abs() <= 2e-13,
                    "rigid mode={mode} component={component} dof={} error={}",
                    line.dof,
                    (reconstructed - expected).abs()
                );
            }
        }
    }
}

fn reduced_matrix(
    k: &solidcraft_trivariate::elasticity::Stiffness,
    lifts: [&Lift; 3],
) -> Vec<Vec<f64>> {
    let mut component_offset = [0; 3];
    for c in 1..3 {
        component_offset[c] = component_offset[c - 1] + lifts[c - 1].free_roots.len();
    }
    let size = component_offset[2] + lifts[2].free_roots.len();
    let mut matrix = vec![vec![0.0; size]; size];
    for i in 0..k.offsets.len() - 1 {
        for position in k.offsets[i]..k.offsets[i + 1] {
            let j = k.columns[position];
            let block = k.values[position];
            for a in 0..3 {
                let row: Vec<_> = lifts[a].lines[i]
                    .terms
                    .iter()
                    .map(|(root, weight)| (component_offset[a] + root, value(*weight)))
                    .collect();
                for b in 0..3 {
                    let column: Vec<_> = lifts[b].lines[j]
                        .terms
                        .iter()
                        .map(|(root, weight)| (component_offset[b] + root, value(*weight)))
                        .collect();
                    for &(r, wr) in &row {
                        for &(c, wc) in &column {
                            matrix[r][c] += wr * block[a][b] * wc;
                        }
                    }
                }
            }
        }
    }
    matrix
}

fn write_matrix(path: &Path, matrix: &[Vec<f64>]) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    let n = matrix.len();
    let entries = (0..n)
        .flat_map(|i| {
            (i..n).filter_map(move |j| (matrix[i][j] != 0.0).then_some((i, j, matrix[i][j])))
        })
        .collect::<Vec<_>>();
    writeln!(file, "%%MatrixMarket matrix coordinate real symmetric")?;
    writeln!(file, "% reduced C^T K C, upper triangle")?;
    writeln!(file, "{n} {n} {}", entries.len())?;
    for (i, j, value) in entries {
        writeln!(file, "{} {} {:.17e}", i + 1, j + 1, value)?;
    }
    Ok(())
}

fn write_metadata(path: &Path, rows: &[(usize, usize, usize, usize)]) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "{{\"family\":\"three_cell_q2_t_junction_v1\",\"expected_free_rigid_modes\":6,\"threshold_ratio\":1e-10,\"orders\":[")?;
    for (index, (order, free_size, pinned_size, seconds_millis)) in rows.iter().enumerate() {
        writeln!(file, "{{\"gauss_order\":{order},\"free_size\":{free_size},\"pinned_size\":{pinned_size},\"elapsed_ms\":{seconds_millis}}}{}", if index + 1 == rows.len() { "" } else { "," })?;
    }
    writeln!(file, "]}}")?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output = PathBuf::from(
        std::env::args()
            .nth(1)
            .unwrap_or_else(|| "build/nested_spectrum".into()),
    );
    fs::create_dir_all(&output)?;
    let mesh = support::t_mesh::base_three_cell()?;
    assert_eq!(mesh.space.blocks.len(), 3);
    assert_eq!(mesh.relations.len(), 3);
    assert_eq!(mesh.space.ndof, 81);
    assert!(mesh.relations.iter().all(|relation| {
        relation.master_side.0 == relation.slave_side.0
            && relation.spans.iter().all(|span| {
                matches!(
                    span,
                    solidcraft_trivariate::face_constraints::Span::Full
                        | solidcraft_trivariate::face_constraints::Span::Lower
                        | solidcraft_trivariate::face_constraints::Span::Upper
                )
            })
    }));
    let free_lift = resolve(mesh.space.ndof, &mesh.equations, &[])?;
    assert_eq!(free_lift.free_roots.len(), 57);
    verify_rigid_lifts(&mesh, [&free_lift, &free_lift, &free_lift]);

    let master = mesh
        .space
        .blocks
        .iter()
        .find(|block| block.id == "master")
        .unwrap();
    let pins = [master.dofs[0], master.dofs[18], master.dofs[3]];
    let mut prescribed: [Vec<(usize, Ratio)>; 3] = std::array::from_fn(|_| Vec::new());
    for component_prescriptions in &mut prescribed {
        component_prescriptions.push((pins[0], Ratio::zero()));
    }
    prescribed[1].push((pins[1], Ratio::zero()));
    prescribed[2].push((pins[1], Ratio::zero()));
    prescribed[2].push((pins[2], Ratio::zero()));
    let fixed_lifts = [
        resolve(mesh.space.ndof, &mesh.equations, &prescribed[0])?,
        resolve(mesh.space.ndof, &mesh.equations, &prescribed[1])?,
        resolve(mesh.space.ndof, &mesh.equations, &prescribed[2])?,
    ];
    assert_eq!(
        [
            fixed_lifts[0].free_roots.len(),
            fixed_lifts[1].free_roots.len(),
            fixed_lifts[2].free_roots.len(),
        ],
        [56, 55, 54]
    );

    let mut metadata = Vec::new();
    for order in [3, 5] {
        let started = Instant::now();
        let stiffness = assemble(&mesh.space, 10.0, 0.25, order)?;
        let free = reduced_matrix(&stiffness, [&free_lift, &free_lift, &free_lift]);
        let fixed = reduced_matrix(
            &stiffness,
            [&fixed_lifts[0], &fixed_lifts[1], &fixed_lifts[2]],
        );
        let directory = output.join(format!("q{order}"));
        fs::create_dir_all(&directory)?;
        write_matrix(&directory.join("free.mtx"), &free)?;
        write_matrix(&directory.join("pinned.mtx"), &fixed)?;
        let elapsed_ms = started.elapsed().as_millis() as usize;
        metadata.push((order, free.len(), fixed.len(), elapsed_ms));
        println!("q={order} dofs={} reduced_free={} reduced_pinned={} interfaces={} elapsed_ms={elapsed_ms}",
            mesh.space.ndof, free.len(), fixed.len(), mesh.relations.len());
    }
    write_metadata(&output.join("metadata.json"), &metadata)?;
    Ok(())
}
