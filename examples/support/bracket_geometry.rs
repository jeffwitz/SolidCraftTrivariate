//! Exact homogeneous geometry of an adaptive field checkpoint, in its leaf order.
use solidcraft_trivariate::analysis_space::AnalysisSpace;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn write(path: &Path, space: &AnalysisSpace) -> Result<(), Box<dyn std::error::Error>> {
    let mut output = BufWriter::new(std::fs::File::create(path)?);
    output.write_all(b"SCQ2GEO1")?;
    output.write_all(&u32::try_from(space.blocks.len())?.to_le_bytes())?;
    for block in &space.blocks {
        output.write_all(&u16::try_from(block.id.len())?.to_le_bytes())?;
        output.write_all(block.id.as_bytes())?;
        for point in &block.controls {
            for value in point {
                output.write_all(&value.to_le_bytes())?;
            }
        }
    }
    output.flush()?;
    Ok(())
}
