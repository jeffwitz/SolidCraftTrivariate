//! Independent audit of ordered coaxial cut/add programs.
//!
//! This module reconstructs the occupancy history from the recipe and compares
//! the saved integer nets with independently derived circular and square
//! charts. It never calls the geometry constructor.
//!
//! The accepted template has a square core of half-width `r_min/4`, concentric
//! rational circular bands with strictly increasing integer radii, and an
//! outer square at half-width 10 mm. Since `sqrt(2) r_min/4 < r_min <= 8 < 10`,
//! these boundaries are strictly nested. Each angular chart is one of the
//! eight monotone half-quadrant arcs; radial charts join consecutive nested
//! boundaries. Exact whole-face matching checks their seams, and the planar
//! Bernstein certificate checks a strictly positive Jacobian on every chart.
//! Axial slabs are disjoint open intervals whose ends form the full sorted
//! partition. Together these facts establish coverage and injectivity for this
//! fixed coaxial template, rather than for arbitrary trivariate nets. The
//! volume formula sums exact disk/annulus/square areas over the audited cells;
//! `f64` is used only to report multiples of pi.

use crate::bracket::BracketBlock;
use crate::bracket_audit::BracketCertificate;
use crate::coaxial_program::{CoaxialProgram, CylinderOp};
use crate::incidence::{enumerate_faces, BlockSide, Incidence};
use crate::planar::certify_planar;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

#[derive(Debug, Clone, PartialEq)]
pub struct CoaxialCertificate {
    pub certificate: BracketCertificate,
    pub ideal_volume_mm3: f64,
    pub cavities: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Cell {
    solid: bool,
    changed_by: usize,
}

#[derive(Clone)]
struct Grid {
    z: Vec<i32>,
    r: Vec<i32>,
    cells: Vec<Vec<Cell>>,
}

fn history(program: &CoaxialProgram) -> Result<Grid, &'static str> {
    if !(1..=16).contains(&program.base_height_mm) || program.features.len() > 8 {
        return Err("coaxial_program_range");
    }
    let mut ids = BTreeSet::new();
    let mut zs = BTreeSet::from([0, program.base_height_mm]);
    let mut rs = BTreeSet::new();
    for f in &program.features {
        let [z0, z1] = f.z_mm;
        if f.id.is_empty()
            || f.id.len() > 48
            || f.id == "base"
            || !f
                .id
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-')
            || !ids.insert(f.id.as_str())
        {
            return Err("coaxial_feature_id");
        }
        if !(1..=8).contains(&f.radius_mm) || z0 < 0 || z1 > 24 || z0 >= z1 {
            return Err("coaxial_feature_range");
        }
        zs.extend([z0, z1]);
        rs.insert(f.radius_mm);
    }
    let z: Vec<_> = zs.into_iter().collect();
    let r: Vec<_> = rs.into_iter().collect();
    let mut cells = Vec::new();
    for slab in z.windows(2) {
        let mut row = Vec::with_capacity(r.len() + 1);
        for band in 0..=r.len() {
            let mut cell = Cell {
                solid: slab[1] <= program.base_height_mm,
                changed_by: 0,
            };
            for (index, f) in program.features.iter().enumerate() {
                let active_z = slab[0] >= f.z_mm[0] && slab[1] <= f.z_mm[1];
                let active_r = r.get(band).is_some_and(|radius| *radius <= f.radius_mm);
                if active_z && active_r {
                    let next_solid = f.op == CylinderOp::Add;
                    if cell.solid != next_solid {
                        cell = Cell {
                            solid: next_solid,
                            changed_by: index + 1,
                        };
                    }
                }
            }
            row.push(cell);
        }
        cells.push(row);
    }
    let grid = Grid { z, r, cells };
    validate_connectivity(&grid)?;
    Ok(grid)
}

fn validate_connectivity(grid: &Grid) -> Result<(), &'static str> {
    let nz = grid.cells.len();
    let nr = grid.r.len() + 1;
    for z in 0..nz.saturating_sub(1) {
        for r in 0..nr.saturating_sub(1) {
            let a = grid.cells[z][r].solid;
            let b = grid.cells[z][r + 1].solid;
            if a != b && a == grid.cells[z + 1][r + 1].solid && b == grid.cells[z + 1][r].solid {
                return Err("coaxial_nonmanifold_edge_contact");
            }
        }
    }
    let total = grid.cells.iter().flatten().filter(|c| c.solid).count();
    let start = grid
        .cells
        .iter()
        .enumerate()
        .find_map(|(z, row)| row.iter().position(|cell| cell.solid).map(|r| (z, r)))
        .ok_or("coaxial_empty_material")?;
    let mut seen = vec![vec![false; nr]; nz];
    seen[start.0][start.1] = true;
    let mut queue = VecDeque::from([start]);
    let mut count = 0;
    while let Some((z, r)) = queue.pop_front() {
        count += 1;
        for (zz, rr) in [
            (z.wrapping_sub(1), r),
            (z + 1, r),
            (z, r.wrapping_sub(1)),
            (z, r + 1),
        ] {
            if zz < nz && rr < nr && grid.cells[zz][rr].solid && !seen[zz][rr] {
                seen[zz][rr] = true;
                queue.push_back((zz, rr));
            }
        }
    }
    if count != total {
        return Err("coaxial_disconnected_material");
    }
    Ok(())
}

type Curve = [[i64; 3]; 3];

fn circle(radius: i32, sector: usize) -> Curve {
    let quadrant = sector / 2;
    let mut half = sector % 2;
    if quadrant % 2 == 1 {
        half = 1 - half;
    }
    let mut net = if half == 0 {
        [[4, 0, 4], [4, 2, 4], [3, 4, 5]]
    } else {
        [[3, 4, 5], [2, 6, 6], [0, 8, 8]]
    };
    if quadrant % 2 == 1 {
        net.reverse();
    }
    for p in &mut net {
        let [x, y, w] = *p;
        let rotated = match quadrant {
            0 => (x, y),
            1 => (-x, y),
            2 => (-x, -y),
            _ => (x, -y),
        };
        *p = [
            4 * i64::from(radius) * rotated.0,
            4 * i64::from(radius) * rotated.1,
            4 * w,
        ];
    }
    net
}

fn square(quarter_mm: i64, sector: usize) -> Curve {
    const P: [[i64; 2]; 8] = [
        [1, 0],
        [1, 1],
        [0, 1],
        [-1, 1],
        [-1, 0],
        [-1, -1],
        [0, -1],
        [1, -1],
    ];
    let [x, y] = P[sector];
    let [xx, yy] = P[(sector + 1) % 8];
    [
        [2 * quarter_mm * x, 2 * quarter_mm * y, 8],
        [quarter_mm * (x + xx), quarter_mm * (y + yy), 8],
        [2 * quarter_mm * xx, 2 * quarter_mm * yy, 8],
    ]
}

fn rectangle(id: String, x: [i64; 2], y: [i64; 2], z: [i32; 2]) -> BracketBlock {
    let mut controls = Vec::with_capacity(8);
    for xx in x {
        for yy in y {
            for zz in z {
                controls.push([xx, yy, 4 * i64::from(zz), 4]);
            }
        }
    }
    BracketBlock {
        id,
        degrees: [1, 1, 1],
        controls,
    }
}

fn sector(id: String, outer: Curve, inner: Curve, z: [i32; 2]) -> BracketBlock {
    let mut controls = Vec::with_capacity(12);
    for i in 0..3 {
        for [x, y, w] in [outer[i], inner[i]] {
            for zz in z {
                controls.push([x, y, i64::from(zz) * w, w]);
            }
        }
    }
    BracketBlock {
        id,
        degrees: [2, 1, 1],
        controls,
    }
}

fn expected_blocks(grid: &Grid) -> Vec<BracketBlock> {
    let mut blocks = Vec::new();
    for z in grid.z.windows(2) {
        let zz = [z[0], z[1]];
        if grid.r.is_empty() {
            blocks.push(rectangle(
                format!("plate_{}_{}", z[0], z[1]),
                [-40, 40],
                [-40, 40],
                zz,
            ));
            continue;
        }
        let core = i64::from(grid.r[0]);
        for (band, cell) in grid.cells[grid.z.windows(2).position(|w| w == z).unwrap()]
            .iter()
            .enumerate()
        {
            if !cell.solid {
                continue;
            }
            if band == 0 {
                for (i, x) in [[-core, 0], [0, core]].into_iter().enumerate() {
                    for (j, y) in [[-core, 0], [0, core]].into_iter().enumerate() {
                        blocks.push(rectangle(
                            format!("core_{}_{}_{}_{}", z[0], z[1], i, j),
                            x,
                            y,
                            zz,
                        ));
                    }
                }
            }
            for s in 0..8 {
                let inner = if band == 0 {
                    square(core, s)
                } else {
                    circle(grid.r[band - 1], s)
                };
                let outer = if band == grid.r.len() {
                    square(40, s)
                } else {
                    circle(grid.r[band], s)
                };
                blocks.push(sector(
                    format!("band_{}_{}_{}_{}", z[0], z[1], band, s),
                    outer,
                    inner,
                    zz,
                ));
            }
        }
    }
    blocks
}

fn cavity_count(grid: &Grid) -> usize {
    let nz = grid.cells.len();
    let nr = grid.r.len() + 1;
    let mut seen = vec![vec![false; nr]; nz];
    let mut cavities = 0;
    for z0 in 0..nz {
        for r0 in 0..nr {
            if grid.cells[z0][r0].solid || seen[z0][r0] {
                continue;
            }
            let mut open = false;
            let mut queue = VecDeque::from([(z0, r0)]);
            seen[z0][r0] = true;
            while let Some((z, r)) = queue.pop_front() {
                if z == 0 || z + 1 == nz || r + 1 == nr {
                    open = true;
                }
                for (zz, rr) in [
                    (z.wrapping_sub(1), r),
                    (z + 1, r),
                    (z, r.wrapping_sub(1)),
                    (z, r + 1),
                ] {
                    if zz < nz && rr < nr && !grid.cells[zz][rr].solid && !seen[zz][rr] {
                        seen[zz][rr] = true;
                        queue.push_back((zz, rr));
                    }
                }
            }
            if !open {
                cavities += 1;
            }
        }
    }
    cavities
}

fn volume(grid: &Grid) -> f64 {
    let mut areas = Vec::with_capacity(grid.r.len() + 1);
    if let Some(&first) = grid.r.first() {
        areas.push(std::f64::consts::PI * f64::from(first).powi(2));
    }
    for pair in grid.r.windows(2) {
        areas
            .push(std::f64::consts::PI * (f64::from(pair[1]).powi(2) - f64::from(pair[0]).powi(2)));
    }
    let last_area = grid.r.last().map_or(400.0, |r| {
        400.0 - std::f64::consts::PI * f64::from(*r).powi(2)
    });
    areas.push(last_area);
    grid.cells
        .iter()
        .enumerate()
        .map(|(z, row)| {
            let dz = f64::from(grid.z[z + 1] - grid.z[z]);
            row.iter()
                .zip(&areas)
                .filter(|(cell, _)| cell.solid)
                .map(|(_, area)| dz * area)
                .sum::<f64>()
        })
        .sum()
}

fn normalize_block(block: &BracketBlock) -> Result<BracketBlock, &'static str> {
    let mut gcd = 0_i64;
    for point in &block.controls {
        for value in point {
            gcd = gcd_i64(gcd, *value).ok_or("coaxial_scale_overflow")?;
        }
    }
    if gcd <= 0 {
        return Err("coaxial_invalid_scale");
    }
    Ok(BracketBlock {
        id: block.id.clone(),
        degrees: block.degrees,
        controls: block.controls.iter().map(|p| p.map(|x| x / gcd)).collect(),
    })
}

fn gcd_i64(mut a: i64, mut b: i64) -> Option<i64> {
    a = a.checked_abs()?;
    b = b.checked_abs()?;
    while b != 0 {
        (a, b) = (b, a % b);
    }
    Some(a)
}

fn provenance(
    program: &CoaxialProgram,
    grid: &Grid,
    side: &BlockSide,
) -> Result<String, &'static str> {
    let block = side.block.as_str();
    let fields: Vec<_> = block.split('_').collect();
    if fields.first() == Some(&"plate") {
        return Ok("base".into());
    }
    let (z0, z1, band) = if fields.first() == Some(&"core") {
        (
            fields
                .get(1)
                .and_then(|s| s.parse::<i32>().ok())
                .ok_or("bad_id")?,
            fields
                .get(2)
                .and_then(|s| s.parse::<i32>().ok())
                .ok_or("bad_id")?,
            0,
        )
    } else if fields.first() == Some(&"band") {
        (
            fields
                .get(1)
                .and_then(|s| s.parse::<i32>().ok())
                .ok_or("bad_id")?,
            fields
                .get(2)
                .and_then(|s| s.parse::<i32>().ok())
                .ok_or("bad_id")?,
            fields
                .get(3)
                .and_then(|s| s.parse::<usize>().ok())
                .ok_or("bad_id")?,
        )
    } else {
        return Err("bad_id");
    };
    let zi = grid
        .z
        .windows(2)
        .position(|w| w == [z0, z1])
        .ok_or("bad_id")?;
    let cell = grid.cells[zi].get(band).ok_or("bad_band")?;
    let mut source = cell.changed_by;
    let neighbor = match side.axis {
        0 => None,
        1 => {
            let other = if side.upper {
                band.checked_sub(1)
            } else {
                Some(band + 1)
            };
            other.and_then(|r| grid.cells.get(zi).and_then(|row| row.get(r)))
        }
        2 => {
            let other = if side.upper {
                Some(zi + 1)
            } else {
                zi.checked_sub(1)
            };
            other.and_then(|z| grid.cells.get(z).and_then(|row| row.get(band)))
        }
        _ => return Err("bad_side"),
    };
    source = source.max(neighbor.map_or(0, |c| c.changed_by));
    if source > 0 {
        Ok(program.features[source - 1].id.clone())
    } else {
        Ok("base".into())
    }
}

/// Re-audit a serialized coaxial model without invoking its constructor.
pub fn audit_coaxial_program(
    program: &CoaxialProgram,
    blocks: &[BracketBlock],
) -> Result<CoaxialCertificate, &'static str> {
    let grid = history(program)?;
    let expected = expected_blocks(&grid);
    let expected_map: BTreeMap<_, _> = expected
        .iter()
        .map(normalize_block)
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .map(|b| (b.id.clone(), b))
        .collect();
    let actual_map: BTreeMap<_, _> = blocks
        .iter()
        .map(normalize_block)
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .map(|b| (b.id.clone(), b))
        .collect();
    if actual_map.len() != blocks.len() || actual_map != expected_map {
        return Err("coaxial_geometry_mismatch");
    }
    let incidence = enumerate_faces(blocks)?;
    let mut det = f64::INFINITY;
    // Geometry arithmetic runs on one canonical representative per block.
    // This avoids overflow for an otherwise valid large projective gauge.
    for block in actual_map.values() {
        if block.degrees == [2, 1, 1] {
            let mut xyw = Vec::with_capacity(6);
            for i in 0..3 {
                for v in 0..2 {
                    let index = (i * 2 + v) * 2;
                    let p = block.controls[index];
                    let q = block.controls[index + 1];
                    if q[0] != p[0] || q[1] != p[1] || q[3] != p[3] {
                        return Err("coaxial_not_extruded");
                    }
                    xyw.push([p[0], p[1], p[3]]);
                }
            }
            let bottom = block.controls[0];
            let top = block.controls[1];
            for (index, point) in block.controls.iter().enumerate() {
                let (reference, error) = if index % 2 == 0 {
                    (bottom, "coaxial_lower_cap_not_planar")
                } else {
                    (top, "coaxial_upper_cap_not_planar")
                };
                if i128::from(point[2]) * i128::from(reference[3])
                    != i128::from(reference[2]) * i128::from(point[3])
                {
                    return Err(error);
                }
            }
            let lower = certify_planar([2, 1], &xyw).map_err(|_| "coaxial_planar_jacobian")?;
            let z0 = bottom[2] as f64 / bottom[3] as f64;
            let z1 = top[2] as f64 / top[3] as f64;
            let dz = z1 - z0;
            det = det.min(lower * dz);
        } else {
            let mut lo = [i64::MAX; 3];
            let mut hi = [i64::MIN; 3];
            for p in &block.controls {
                if p[3] != block.controls[0][3] || p[3] <= 0 {
                    return Err("coaxial_rectangle_gauge");
                }
                for a in 0..3 {
                    let x = p[a];
                    lo[a] = lo[a].min(x);
                    hi[a] = hi[a].max(x);
                }
            }
            let gauge = block.controls[0][3] as f64;
            let v = ((hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2])) as f64 / gauge.powi(3);
            det = det.min(v);
        }
    }
    if !det.is_finite() || det <= 0.0 {
        return Err("coaxial_jacobian");
    }
    let design_faces = incidence
        .exterior
        .iter()
        .map(|side| Ok((side.clone(), provenance(program, &grid, side)?)))
        .collect::<Result<Vec<_>, &'static str>>()?;
    let certificate = BracketCertificate {
        incidence: Incidence {
            interfaces: incidence.interfaces,
            exterior: incidence.exterior,
        },
        determinant_lower_bound_mm3: det,
        radial_error_upper_mm: 0.0,
        design_faces,
    };
    Ok(CoaxialCertificate {
        certificate,
        ideal_volume_mm3: volume(&grid),
        cavities: cavity_count(&grid),
    })
}
