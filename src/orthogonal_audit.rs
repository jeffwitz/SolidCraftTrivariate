//! Independent exact audit for ordered orthogonal box programs.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    incidence::{enumerate_faces, Incidence},
    orthogonal_program::{Bounds, BoxOp, OrthogonalProgram},
};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

#[derive(Debug, Clone, PartialEq)]
pub struct OrthogonalCertificate {
    pub certificate: BracketCertificate,
    pub volume_mm3: i64,
    pub cavities: usize,
    pub minimum_x_mm: i32,
    /// Stable block ID and the last operation that changed its material state.
    pub block_provenance: Vec<(String, String)>,
}

fn valid_bounds(b: &Bounds) -> bool {
    b.iter()
        .all(|&[lo, hi]| (-64..hi).contains(&lo) && hi <= 64)
}
fn encloses(a: &Bounds, b: &Bounds) -> bool {
    (0..3).all(|d| a[d][0] <= b[d][0] && b[d][1] <= a[d][1])
}
fn axes(program: &OrthogonalProgram) -> [Vec<i32>; 3] {
    std::array::from_fn(|d| {
        let mut set = BTreeSet::from(program.base_mm[d]);
        for f in &program.features {
            set.extend(f.bounds_mm[d]);
        }
        set.into_iter().collect()
    })
}
fn bounds_at(axes: &[Vec<i32>; 3], p: [usize; 3]) -> Bounds {
    std::array::from_fn(|d| [axes[d][p[d]], axes[d][p[d] + 1]])
}
fn index(p: [usize; 3], dims: [usize; 3]) -> usize {
    (p[0] * dims[1] + p[1]) * dims[2] + p[2]
}
fn coords(i: usize, dims: [usize; 3]) -> [usize; 3] {
    [i / (dims[1] * dims[2]), i / dims[2] % dims[1], i % dims[2]]
}
fn occupied(solid: &[bool], dims: [usize; 3], p: [isize; 3]) -> bool {
    if (0..3).any(|d| p[d] < 0 || p[d] >= dims[d] as isize) {
        false
    } else {
        solid[index(p.map(|x| x as usize), dims)]
    }
}
fn validate_program(program: &OrthogonalProgram) -> Result<(), &'static str> {
    if !valid_bounds(&program.base_mm) || program.features.len() > 8 {
        return Err("orthogonal_program_range");
    }
    let mut ids = BTreeSet::new();
    for f in &program.features {
        if f.id.is_empty()
            || f.id.len() > 48
            || f.id == "base"
            || !f.id.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_')
            || !ids.insert(&f.id)
        {
            return Err("orthogonal_feature_id");
        }
        if !valid_bounds(&f.bounds_mm) {
            return Err("orthogonal_feature_bounds");
        }
    }
    Ok(())
}
fn check_local_links(solid: &[bool], dims: [usize; 3]) -> Result<(), &'static str> {
    for i in 0..=dims[0] {
        for j in 0..=dims[1] {
            for k in 0..=dims[2] {
                let v = [i, j, k];
                let mask = (0..8).fold(0u8, |m, o| {
                    let p = std::array::from_fn(|a| v[a] as isize - 1 + ((o >> a) & 1) as isize);
                    m | (u8::from(occupied(solid, dims, p)) << o)
                });
                if mask == 0 || mask == 255 {
                    continue;
                }
                for a in 0..3 {
                    let b = (a + 1) % 3;
                    let c = (a + 2) % 3;
                    for s in 0..2 {
                        let base = s << a;
                        let q = [
                            base,
                            base | (1 << b),
                            base | (1 << b) | (1 << c),
                            base | (1 << c),
                        ]
                        .map(|x| (mask >> x) & 1);
                        if q == [1, 0, 1, 0] || q == [0, 1, 0, 1] {
                            return Err("orthogonal_nonmanifold_edge");
                        }
                    }
                }
                if !octants_connected(mask) || !octants_connected(!mask) {
                    return Err("orthogonal_nonmanifold_vertex");
                }
            }
        }
    }
    Ok(())
}
fn octants_connected(mask: u8) -> bool {
    if mask == 0 {
        return true;
    }
    let mut seen = 1u8 << mask.trailing_zeros();
    loop {
        let mut next = seen;
        for i in 0..8 {
            if seen & (1 << i) != 0 {
                for a in 0..3 {
                    next |= (1 << (i ^ (1 << a))) & mask
                }
            }
        }
        if next == seen {
            return seen == mask;
        }
        seen = next
    }
}
fn check_components(solid: &[bool], dims: [usize; 3]) -> Result<usize, &'static str> {
    let seed = solid
        .iter()
        .position(|x| *x)
        .ok_or("orthogonal_empty_material")?;
    let mut seen = vec![false; solid.len()];
    seen[seed] = true;
    let mut q = VecDeque::from([seed]);
    let mut material = 0;
    while let Some(i) = q.pop_front() {
        material += 1;
        let p = coords(i, dims);
        for a in 0..3 {
            for d in [-1, 1] {
                let mut n = p.map(|x| x as isize);
                n[a] += d;
                if occupied(solid, dims, n) {
                    let j = index(n.map(|x| x as usize), dims);
                    if !seen[j] {
                        seen[j] = true;
                        q.push_back(j)
                    }
                }
            }
        }
    }
    if material != solid.iter().filter(|x| **x).count() {
        return Err("orthogonal_disconnected_material");
    }
    let mut air_seen = vec![false; solid.len()];
    let mut cavities = 0;
    for start in 0..solid.len() {
        if solid[start] || air_seen[start] {
            continue;
        }
        let mut touches = false;
        air_seen[start] = true;
        let mut q = VecDeque::from([start]);
        while let Some(i) = q.pop_front() {
            let p = coords(i, dims);
            touches |= (0..3).any(|a| p[a] == 0 || p[a] + 1 == dims[a]);
            for a in 0..3 {
                for d in [-1, 1] {
                    let mut n = p.map(|x| x as isize);
                    n[a] += d;
                    if (0..3).all(|axis| n[axis] >= 0 && n[axis] < dims[axis] as isize)
                        && !occupied(solid, dims, n)
                    {
                        let j = index(n.map(|x| x as usize), dims);
                        if !air_seen[j] {
                            air_seen[j] = true;
                            q.push_back(j)
                        }
                    }
                }
            }
        }
        if !touches {
            cavities += 1
        }
    }
    Ok(cavities)
}
fn actual_bounds(block: &BracketBlock) -> Result<Bounds, &'static str> {
    if block.degrees != [1, 1, 1] || block.controls.len() != 8 {
        return Err("orthogonal_altered_block");
    }
    let mut lo = [i64::MAX; 3];
    let mut hi = [i64::MIN; 3];
    for p in &block.controls {
        if p[3] != 1 {
            return Err("orthogonal_altered_block");
        }
        for d in 0..3 {
            lo[d] = lo[d].min(p[d]);
            hi[d] = hi[d].max(p[d]);
        }
    }
    if (0..3).any(|d| lo[d] < -64 || hi[d] > 64 || lo[d] >= hi[d]) {
        return Err("orthogonal_altered_block");
    }
    let b = std::array::from_fn(|d| [lo[d] as i32, hi[d] as i32]);
    let expected = (0..2)
        .flat_map(|a| {
            (0..2).flat_map(move |b| {
                (0..2).map(move |c| {
                    [
                        if a == 0 { lo[0] } else { hi[0] },
                        if b == 0 { lo[1] } else { hi[1] },
                        if c == 0 { lo[2] } else { hi[2] },
                        1,
                    ]
                })
            })
        })
        .collect::<Vec<_>>();
    if block.controls != expected {
        return Err("orthogonal_altered_block");
    }
    Ok(b)
}

/// Audit supplied boxes against an independently reconstructed exact occupancy grid.
pub fn audit_orthogonal_program(
    program: &OrthogonalProgram,
    blocks: &[BracketBlock],
) -> Result<OrthogonalCertificate, &'static str> {
    validate_program(program)?;
    let planes = axes(program);
    let dims = planes.each_ref().map(|a| a.len() - 1);
    let count = dims.iter().product::<usize>();
    let mut solid = vec![false; count];
    let mut last = vec!["".to_owned(); count];
    let mut last_order = vec![None::<usize>; count];
    for n in 0..count {
        let b = bounds_at(&planes, coords(n, dims));
        let mut inside = encloses(&program.base_mm, &b);
        if inside {
            last[n] = "base".into();
            last_order[n] = Some(0);
        }
        for (operation, f) in program.features.iter().enumerate() {
            if encloses(&f.bounds_mm, &b) {
                let next = f.operation == BoxOp::Add;
                if next != inside {
                    last[n] = f.id.clone();
                    last_order[n] = Some(operation + 1);
                }
                inside = next
            }
        }
        solid[n] = inside;
    }
    check_local_links(&solid, dims)?;
    let cavities = check_components(&solid, dims)?;
    let mut expected = BTreeMap::<String, (Bounds, String)>::new();
    let mut volume = 0_i64;
    let mut minimum_x = i32::MAX;
    let mut min_cell_volume = i64::MAX;
    for (i, &is_solid) in solid.iter().enumerate() {
        if !is_solid {
            continue;
        }
        let b = bounds_at(&planes, coords(i, dims));
        let id = format!(
            "box_{}_{}_{}_{}_{}_{}",
            b[0][0], b[0][1], b[1][0], b[1][1], b[2][0], b[2][1]
        );
        let v = (0..3).try_fold(1_i64, |acc, d| {
            acc.checked_mul(i64::from(b[d][1] - b[d][0]))
                .ok_or("orthogonal_volume_overflow")
        })?;
        volume = volume.checked_add(v).ok_or("orthogonal_volume_overflow")?;
        min_cell_volume = min_cell_volume.min(v);
        minimum_x = minimum_x.min(b[0][0]);
        expected.insert(id, (b, last[i].clone()));
    }
    if blocks.len() != expected.len() {
        return Err("orthogonal_coverage_mismatch");
    }
    let mut seen = BTreeSet::new();
    let mut provenance = Vec::new();
    for block in blocks {
        if !seen.insert(block.id.clone()) {
            return Err("orthogonal_duplicate_block");
        }
        let b = actual_bounds(block)?;
        let key = format!(
            "box_{}_{}_{}_{}_{}_{}",
            b[0][0], b[0][1], b[1][0], b[1][1], b[2][0], b[2][1]
        );
        let Some((want, source)) = expected.get(&key) else {
            return Err("orthogonal_coverage_mismatch");
        };
        if want != &b || block.id != key {
            return Err("orthogonal_altered_block");
        }
        provenance.push((block.id.clone(), source.clone()));
    }
    let incidence: Incidence = enumerate_faces(blocks)?;
    let mut design_faces = Vec::with_capacity(incidence.exterior.len());
    let bounds_by_id = expected
        .iter()
        .map(|(id, (bounds, _))| (id.as_str(), *bounds))
        .collect::<BTreeMap<_, _>>();
    for side in &incidence.exterior {
        let b = bounds_by_id
            .get(side.block.as_str())
            .ok_or("orthogonal_incidence_unknown_block")?;
        let mut p = [0_usize; 3];
        for d in 0..3 {
            p[d] = planes[d]
                .binary_search(&b[d][0])
                .map_err(|_| "orthogonal_incidence_plane")?;
        }
        if side.upper {
            p[side.axis] += 1;
        }
        let neighbor = if side.upper {
            p[side.axis] as isize
        } else {
            p[side.axis] as isize - 1
        };
        p[side.axis] = if side.upper {
            p[side.axis] - 1
        } else {
            p[side.axis]
        };
        let inside_id = index(p, dims);
        if !solid[inside_id] {
            return Err("orthogonal_incidence_material_mismatch");
        }
        let outside_id = if neighbor >= 0 && neighbor < dims[side.axis] as isize {
            let mut q = p;
            q[side.axis] = neighbor as usize;
            let id = index(q, dims);
            if solid[id] {
                return Err("orthogonal_incidence_material_mismatch");
            }
            Some(id)
        } else {
            None
        };
        let inner_order = last_order[inside_id].ok_or("orthogonal_missing_face_provenance")?;
        let source = if let Some(outside_id) = outside_id {
            match last_order[outside_id] {
                Some(order) if order > inner_order => last[outside_id].as_str(),
                _ => last[inside_id].as_str(),
            }
        } else {
            last[inside_id].as_str()
        };
        if source.is_empty() {
            return Err("orthogonal_missing_face_provenance");
        }
        design_faces.push((side.clone(), source.to_owned()));
    }
    design_faces.sort_by(|a, b| a.0.cmp(&b.0));
    let cert = BracketCertificate {
        incidence,
        determinant_lower_bound_mm3: min_cell_volume as f64,
        radial_error_upper_mm: 0.,
        design_faces,
    };
    provenance.sort();
    Ok(OrthogonalCertificate {
        certificate: cert,
        volume_mm3: volume,
        cavities,
        minimum_x_mm: minimum_x,
        block_provenance: provenance,
    })
}
