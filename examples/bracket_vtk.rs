//! Export a diagnostic hexahedral sampling of the certified bracket maps.

use solidcraft_trivariate::{
    bracket::bracket_blocks, bracket_audit::audit_bracket, RationalBezierMap,
};
use std::{
    error::Error,
    fs::File,
    io::{BufWriter, Write},
};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .ok_or("usage: bracket_vtk OUTPUT.vtk SUBDIVISIONS")?;
    let n: usize = args.next().ok_or("missing subdivisions")?.parse()?;
    if n == 0 || args.next().is_some() {
        return Err("expected positive subdivision count".into());
    }
    let mut blocks = bracket_blocks()?;
    audit_bracket(&blocks)?;
    blocks.sort_by(|a, b| a.id.cmp(&b.id));
    let m = n.checked_add(1).ok_or("size overflow")?;
    let points = m
        .checked_pow(3)
        .and_then(|x| x.checked_mul(blocks.len()))
        .ok_or("size overflow")?;
    let cells = n
        .checked_pow(3)
        .and_then(|x| x.checked_mul(blocks.len()))
        .ok_or("size overflow")?;
    let entries = cells.checked_mul(9).ok_or("size overflow")?;
    let maps = blocks
        .iter()
        .map(|b| {
            RationalBezierMap::new(
                b.degrees,
                b.controls.iter().map(|p| p.map(|x| x as f64)).collect(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out,"# vtk DataFile Version 3.0\nCertified bracket: diagnostic map samples in mm\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS {points} double")?;
    for map in &maps {
        for i in 0..=n {
            for j in 0..=n {
                for k in 0..=n {
                    let [x, y, z] = map
                        .evaluate([i, j, k].map(|i| i as f64 / n as f64))?
                        .position;
                    writeln!(out, "{x} {y} {z}")?;
                }
            }
        }
    }
    writeln!(out, "CELLS {cells} {entries}")?;
    for block in 0..blocks.len() {
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let index = |di, dj, dk| ((block * m + i + di) * m + j + dj) * m + k + dk;
                    let vs: [[usize; 3]; 8] = [
                        [0, 0, 0],
                        [1, 0, 0],
                        [1, 1, 0],
                        [0, 1, 0],
                        [0, 0, 1],
                        [1, 0, 1],
                        [1, 1, 1],
                        [0, 1, 1],
                    ];
                    write!(out, "8")?;
                    for [di, dj, dk] in vs {
                        write!(out, " {}", index(di, dj, dk))?;
                    }
                    writeln!(out)?;
                }
            }
        }
    }
    writeln!(out, "CELL_TYPES {cells}")?;
    for _ in 0..cells {
        writeln!(out, "12")?;
    }
    writeln!(
        out,
        "CELL_DATA {cells}\nSCALARS block_index int 1\nLOOKUP_TABLE default"
    )?;
    for id in 0..blocks.len() {
        for _ in 0..n.pow(3) {
            writeln!(out, "{id}")?;
        }
    }
    for normalized in [false, true] {
        let name = if normalized {
            "sampled_scaled_jacobian"
        } else {
            "sampled_det_mm3"
        };
        writeln!(out, "SCALARS {name} double 1\nLOOKUP_TABLE default")?;
        for map in &maps {
            for i in 0..n {
                for j in 0..n {
                    for k in 0..n {
                        let [a, b, c] = map
                            .evaluate([i, j, k].map(|i| (i as f64 + 0.5) / n as f64))?
                            .jacobian_columns;
                        let det = a[0] * (b[1] * c[2] - b[2] * c[1])
                            + a[1] * (b[2] * c[0] - b[0] * c[2])
                            + a[2] * (b[0] * c[1] - b[1] * c[0]);
                        let norm = |v: [f64; 3]| v.into_iter().map(|x| x * x).sum::<f64>().sqrt();
                        let value = if normalized {
                            det / (norm(a) * norm(b) * norm(c))
                        } else {
                            det
                        };
                        writeln!(out, "{value}")?;
                    }
                }
            }
        }
    }
    out.flush()?;
    eprintln!(
        "{points} points, {cells} sampled hexahedra across {} blocks",
        blocks.len()
    );
    Ok(())
}
