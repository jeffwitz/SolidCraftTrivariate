//! Write a lossless binary64 checkpoint of the rational Q2 displacement field.
//! Little endian: magic SCQ2FLD1, u8 refinement level, x and y radius ticks,
//! Gauss order, u32 block count. Each block carries a u16 byte length, UTF-8
//! ID, then 27 rows of (weight, ux, uy, uz) as f64 in w-fast order.
use solidcraft_trivariate::analysis_space::AnalysisSpace;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn write(
    path: &Path,
    level: usize,
    radii: [u8; 2],
    q: usize,
    space: &AnalysisSpace,
    displacement: &[f64],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = BufWriter::new(std::fs::File::create(path)?);
    file.write_all(b"SCQ2FLD1")?;
    file.write_all(&[level.try_into()?, radii[0], radii[1], q.try_into()?])?;
    file.write_all(&u32::try_from(space.blocks.len())?.to_le_bytes())?;
    for block in &space.blocks {
        file.write_all(&u16::try_from(block.id.len())?.to_le_bytes())?;
        file.write_all(block.id.as_bytes())?;
        for (&d, h) in block.dofs.iter().zip(&block.controls) {
            file.write_all(&(h[3] as f64).to_le_bytes())?;
            for value in &displacement[3 * d..3 * d + 3] {
                file.write_all(&value.to_le_bytes())?;
            }
        }
    }
    file.flush()?;
    Ok(())
}
