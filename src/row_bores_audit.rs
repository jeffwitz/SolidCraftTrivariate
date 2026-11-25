//! Continuous audit of a row of separately contained circular holes.
use crate::bracket::BracketBlock;
use crate::bracket_audit::{BracketCertificate, HOLE_PROFILE};
use crate::incidence::enumerate_faces;
use crate::planar::certify_planar_refined;
use std::collections::{BTreeMap, BTreeSet};

const NAMES: [&str; 4] = ["east", "north", "west", "south"];

fn turn([x, y]: [i128; 2], q: usize) -> [i128; 2] {
    match q {
        0 => [x, y],
        1 => [-y, x],
        2 => [-x, -y],
        _ => [y, -x],
    }
}

fn certify_block(
    block: &BracketBlock,
    half_y: i128,
    x: [i128; 2],
    center: [i128; 2],
    r: i128,
    q: usize,
) -> Result<f64, &'static str> {
    if block.degrees != [2, 1, 1] || block.controls.len() != 12 {
        return Err("row_bore_layout");
    }
    let corners = [
        [x[1], -half_y],
        [x[1], half_y],
        [x[0], half_y],
        [x[0], -half_y],
    ];
    let (a, b) = (corners[q], corners[(q + 1) % 4]);
    let mut profile = [[0; 3]; 6];
    let common = i128::from(block.controls[0][3]);
    for i in 0..3 {
        for j in 0..2 {
            let p = block.controls[4 * i + 2 * j];
            let upper = block.controls[4 * i + 2 * j + 1];
            if p[3] <= 0
                || p[2] != 0
                || upper[..2] != p[..2]
                || upper[3] != p[3]
                || i128::from(upper[2]) != 8 * i128::from(p[3])
            {
                return Err("row_bore_extrusion");
            }
            profile[2 * i + j] = [p[0], p[1], p[3]];
            let weight = i128::from(p[3]);
            let h = HOLE_PROFILE[2 * i + j].map(i128::from);
            if weight * i128::from(HOLE_PROFILE[0][2]) != common * h[2] {
                return Err("row_bore_weight_net");
            }
            let radial = turn([h[0], h[1]], q);
            for k in 0..2 {
                if j == 0 {
                    // Endpoints are in half-mm ticks and the middle is their mean.
                    let target = (2 - i as i128) * a[k] + i as i128 * b[k];
                    if 4 * i128::from(p[k]) != target * weight {
                        return Err("row_bore_outer_edge");
                    }
                } else if 8 * i128::from(p[k]) * h[2]
                    != (4 * center[k] * h[2] + r * radial[k]) * weight
                {
                    return Err("row_bore_arc_identity");
                }
            }
        }
    }
    certify_planar_refined([2, 1], &profile)
        .map(|v| 8. * v)
        .map_err(|e| match e {
            crate::planar::PlanarError::Invalid => "row_bore_nonpositive_jacobian",
            crate::planar::PlanarError::Indeterminate => "row_bore_indeterminate_bernstein",
        })
}

pub fn audit_row(
    half_extent_ticks: [i32; 2],
    centers_ticks: &[[i32; 2]],
    radii_ticks: &[i32],
    blocks: &[BracketBlock],
) -> Result<BracketCertificate, &'static str> {
    let [a, b] = half_extent_ticks.map(i128::from);
    let n = centers_ticks.len();
    if !(1..=16).contains(&n)
        || n != radii_ticks.len()
        || blocks.len() != 4 * n
        || !(1..=160).contains(&a)
        || !(1..=20).contains(&b)
        || (2 * a) % n as i128 != 0
    {
        return Err("row_bore_recipe_or_inventory");
    }
    let width = 2 * a / n as i128;
    let mut expected = BTreeMap::new();
    for (i, (&c, &r)) in centers_ticks.iter().zip(radii_ticks).enumerate() {
        let lo = -a + i as i128 * width;
        let center = c.map(i128::from);
        let r = i128::from(r);
        if r <= 0 || center[0] - r <= lo || center[0] + r >= lo + width || center[1].abs() + r >= b
        {
            return Err("row_bore_disk_not_strictly_inside_strip");
        }
        for (q, name) in NAMES.iter().enumerate() {
            expected.insert(
                format!("row_bore_cell_{i:02}_sector_{name}_z0_8"),
                (i, q, [lo, lo + width], center, r),
            );
        }
    }
    let mut seen = BTreeSet::new();
    let mut det = f64::INFINITY;
    for block in blocks {
        let &(_, q, x, c, r) = expected.get(&block.id).ok_or("row_bore_identity")?;
        if !seen.insert(&block.id) {
            return Err("row_bore_duplicate_identity");
        }
        det = det.min(certify_block(block, b, x, c, r, q)?);
    }
    // Within each strip the parameter annulus has one simple outer rectangle
    // and one simple internal circle with E6's inherited error bound. Positive
    // maps have positive degree and their shared interfaces cancel. Different
    // strips have disjoint interiors by construction.
    let incidence = enumerate_faces(blocks)?;
    if incidence.interfaces.len() != 5 * n - 1 || incidence.exterior.len() != 14 * n + 2 {
        return Err("row_bore_incidence");
    }
    let mut design = vec![];
    for s in &incidence.exterior {
        let &(i, q, _, _, _) = expected.get(&s.block).ok_or("row_bore_identity")?;
        let label = match (s.axis, s.upper) {
            (1, true) => format!("bore_{i:02}"),
            (2, false) => "bottom".into(),
            (2, true) => "top".into(),
            (1, false) => match q {
                0 if i + 1 == n => "right",
                1 => "back",
                2 if i == 0 => "left",
                3 => "front",
                _ => return Err("missing strip interface"),
            }
            .into(),
            _ => return Err("unexpected row boundary"),
        };
        design.push((s.clone(), label));
    }
    design.sort_by(|a, b| a.0.cmp(&b.0));
    let rmax = f64::from(*radii_ticks.iter().max().unwrap());
    Ok(BracketCertificate {
        incidence,
        determinant_lower_bound_mm3: det,
        radial_error_upper_mm: (4.20e-15 * rmax / 8.).next_up(),
        design_faces: design,
    })
}
