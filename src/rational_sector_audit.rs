//! An exact certificate for the stored homogeneous nets of a rational angular sector.
//!
//! At three angular nodes and both radial and axial ends, check all four
//! homogeneous polynomials against the analytic polar chart, up to one positive
//! factor per chart. This tensor grid is unisolvent for degree (2,1,1).
//! The weight is 1+s², the polar angle is 2 atan(s), and radius decreases
//! linearly. These identities on the closed cube establish an inverse and
//! show that distinct chart interiors have disjoint angular intervals.

use crate::{
    bracket::BracketBlock, bracket_audit::BracketCertificate, incidence::enumerate_faces,
    rational_sector::RationalSector,
};

pub fn audit_rational_sector(
    recipe: &RationalSector,
    blocks: &[BracketBlock],
) -> Result<BracketCertificate, &'static str> {
    let [n, d] = recipe.half_angle_tangent.map(i64::from);
    let [ri, ro] = recipe.radii_mm.map(i64::from);
    let h = i64::from(recipe.height_mm);
    if !(1..=256).contains(&n)
        || !(1..=256).contains(&d)
        || 4 * n < d
        || n > 4 * d
        || !(1..ro).contains(&ri)
        || ro > 64
        || !(1..=64).contains(&h)
    {
        return Err("sector_parameters_outside_certificate");
    }
    // Compare against 1+sqrt(2) without rounding the breakpoint.
    let split = n > d && (n - d) * (n - d) > 2 * d * d;
    let charts: &[(&str, i64, i64)] = if split {
        &[("first", 0, 1), ("second", 2, 1)]
    } else {
        &[("whole", 0, 2)]
    };
    if blocks.len() != charts.len() {
        return Err("sector_chart_inventory");
    }
    let q = i128::from(4 * d);
    for &(id, first, stride) in charts {
        let block = blocks
            .iter()
            .find(|b| b.id == id)
            .ok_or("sector_chart_identity")?;
        if block.degrees != [2, 1, 1]
            || block.controls.len() != 12
            || block.controls.iter().any(|p| p[3] <= 0)
        {
            return Err("sector_control_net_dimensions_or_weight");
        }
        let w0 = i128::from(block.controls[0][3]);
        let s0 = i128::from(first * n);
        let w_ref = q * q + s0 * s0;
        for (j, blend) in [[4_i128, 0, 0], [1, 2, 1], [0, 0, 4]].iter().enumerate() {
            let s = i128::from((first + j as i64 * stride) * n);
            let weight = q * q + s * s;
            for (v, radius) in [ro, ri].into_iter().enumerate() {
                let r = i128::from(radius);
                for (w, z) in [0, h].into_iter().enumerate() {
                    let wanted = [
                        r * (q * q - s * s),
                        2 * r * q * s,
                        i128::from(z) * weight,
                        weight,
                    ];
                    for (c, ideal) in wanted.into_iter().enumerate() {
                        let value: i128 = blend
                            .iter()
                            .enumerate()
                            .map(|(i, &a)| a * i128::from(block.controls[4 * i + 2 * v + w][c]))
                            .sum();
                        if value * w_ref != 4 * w0 * ideal {
                            return Err("sector_homogeneous_polynomial_identity");
                        }
                    }
                }
            }
        }
    }
    let incidence = enumerate_faces(blocks)?;
    if incidence.interfaces.len() + 1 != charts.len()
        || incidence.exterior.len() != 4 * charts.len() + 2
    {
        return Err("sector_chart_incidence");
    }
    let mut design_faces = Vec::new();
    for side in &incidence.exterior {
        let tag = match (side.axis, side.upper) {
            (0, false) if side.block != "second" => "angle_start",
            (0, true) if side.block != "first" => "angle_end",
            (1, false) => "outer",
            (1, true) => "inner",
            (2, false) => "start",
            (2, true) => "end",
            _ => return Err("sector_external_side"),
        };
        design_faces.push((side.clone(), tag.to_owned()));
    }
    // Integer products in this bounded domain are exactly representable in binary64.
    let numerator = 2 * n * d * ri * (ro - ri) * h;
    let denominator = charts.len() as i64 * (d * d + n * n);
    Ok(BracketCertificate {
        incidence,
        design_faces,
        determinant_lower_bound_mm3: (numerator as f64 / denominator as f64).next_down(),
        radial_error_upper_mm: 0.,
    })
}
