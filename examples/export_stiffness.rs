use solidcraft_trivariate::{
    analysis_space::build_analysis_space, bracket::bracket_blocks, elasticity::assemble,
};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

fn run() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let order: usize = args
        .next()
        .ok_or("usage: export_stiffness <quadrature_order> <output_directory>")?
        .parse()?;
    let directory = PathBuf::from(
        args.next()
            .ok_or("usage: export_stiffness <quadrature_order> <output_directory>")?,
    );
    if args.next().is_some() {
        return Err("usage: export_stiffness <quadrature_order> <output_directory>".into());
    }

    let blocks = bracket_blocks().map_err(|reason| format!("bracket build failed: {reason}"))?;
    let space = build_analysis_space(&blocks)
        .map_err(|reason| format!("analysis-space build failed: {reason}"))?;
    let started = Instant::now();
    let stiffness = assemble(&space, 210_000.0, 0.3, order)
        .map_err(|reason| format!("stiffness assembly failed: {reason}"))?;
    let assembly_seconds = started.elapsed().as_secs_f64();

    fs::create_dir_all(&directory)?;
    let vector_dofs = 3 * space.ndof;
    let nnz: usize = stiffness
        .values
        .iter()
        .map(|block| {
            block
                .iter()
                .flatten()
                .filter(|&&value| value != 0.0)
                .count()
        })
        .sum();
    let matrix_path = directory.join("stiffness.mtx");
    let mut matrix = BufWriter::new(File::create(&matrix_path)?);
    writeln!(matrix, "%%MatrixMarket matrix coordinate real general")?;
    writeln!(matrix, "% SolidCraftTrivariate fixed bracket Q2 elasticity")?;
    writeln!(
        matrix,
        "% E=210000 MPa, nu=0.3; row/column order is (dof,x/y/z)"
    )?;
    writeln!(matrix, "{vector_dofs} {vector_dofs} {nnz}")?;
    for scalar_row in 0..space.ndof {
        for block_entry in stiffness.offsets[scalar_row]..stiffness.offsets[scalar_row + 1] {
            let scalar_column = stiffness.columns[block_entry];
            let block = &stiffness.values[block_entry];
            for (row_component, row) in block.iter().enumerate() {
                for (column_component, &value) in row.iter().enumerate() {
                    if value != 0.0 {
                        let row_id = 3 * scalar_row + row_component + 1;
                        let column_id = 3 * scalar_column + column_component + 1;
                        writeln!(matrix, "{row_id} {column_id} {value:.17e}")?;
                    }
                }
            }
        }
    }
    matrix.flush()?;

    let y50 = space
        .boundary_dofs
        .get("y50")
        .ok_or("audited analysis space has no y50 boundary")?;
    let y50_path = directory.join("y50_dofs.txt");
    let mut dofs = BufWriter::new(File::create(&y50_path)?);
    writeln!(
        dofs,
        "# zero-based scalar H1 dof ids on the fixed design face y50"
    )?;
    for dof in y50 {
        writeln!(dofs, "{dof}")?;
    }
    dofs.flush()?;

    println!(
        "{{\"quadrature_order\":{order},\"scalar_dofs\":{},\"vector_dofs\":{vector_dofs},\"matrix_nnz\":{nnz},\"fixed_y50_scalar_dofs\":{},\"integrated_volume_mm3\":{:.17e},\"assembly_seconds\":{:.9},\"matrix\":{},\"fixed_dofs\":{}}}",
        space.ndof,
        y50.len(),
        stiffness.integrated_volume,
        assembly_seconds,
        matrix_path.display(),
        y50_path.display(),
    );
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("export_stiffness: {error}");
        std::process::exit(2);
    }
}
