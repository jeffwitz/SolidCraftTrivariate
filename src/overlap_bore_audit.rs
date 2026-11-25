//! Continuous certificate for a square prism minus the union of the two fixed bores.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    incidence::enumerate_faces,
    planar::{certify_planar_refined, PlanarError},
};
use std::collections::{BTreeMap, BTreeSet};

// A common vertical scale maps the ideal union onto the represented union.
// P/Q approximates sqrt(3); the same affine map applies to both conics.
const P: i128 = 1_732_050_807_568_877;
const Q: i128 = 1_000_000_000_000_000;
const KNOT_NUMERATORS: [i128; 5] = [-3, -1, 0, 1, 3];
const OUTER: [[i128; 2]; 5] = [[0, -10], [10, -10], [10, 0], [10, 10], [0, 10]];

fn certify_chart(block: &BracketBlock, sign: i128, interval: usize) -> Result<f64, &'static str> {
    if block.degrees != [2, 1, 1] || block.controls.len() != 12 {
        return Err("overlap_bore_layout");
    }
    let mut profile = [[0; 3]; 6];
    for i in 0..3 {
        for j in 0..2 {
            let p = block.controls[4 * i + 2 * j];
            let cap = block.controls[4 * i + 2 * j + 1];
            if p[2] != 0
                || p[3] <= 0
                || p[..2] != cap[..2]
                || p[3] != cap[3]
                || i128::from(cap[2]) != 8 * i128::from(p[3])
            {
                return Err("overlap_bore_extrusion");
            }
            profile[2 * i + j] = [p[0], p[1], p[3]];
        }
        let [x, y, w] = profile[2 * i].map(i128::from);
        if w != i128::from(profile[2 * i + 1][2]) {
            return Err("overlap_bore_radial_weights");
        }
        for (k, h) in [x, y].into_iter().enumerate() {
            let start = OUTER[interval][k];
            let end = OUTER[interval + 1][k];
            if 2 * h != sign * w * ((2 - i as i128) * start + i as i128 * end) {
                return Err("overlap_bore_outer_face");
            }
        }
    }
    // Clearing the denominators gives two degree-four polynomials in u.
    // Exact equality at five distinct parameters determines each polynomial.
    let d = 12i128;
    for n in 0..=4i128 {
        let t = (4 - n) * KNOT_NUMERATORS[interval] + n * KNOT_NUMERATORS[interval + 1];
        let denominator = d * d + 3 * t * t;
        let weights = [(4 - n).pow(2), 2 * n * (4 - n), n * n];
        let [x, y, w]: [i128; 3] = std::array::from_fn(|c| {
            (0..3)
                .map(|i| weights[i] * i128::from(profile[2 * i + 1][c]))
                .sum()
        });
        if x * denominator != sign * w * (denominator + 2 * (d * d - 3 * t * t)) {
            return Err("overlap_bore_horizontal_conic");
        }
        let y_left = (y * denominator)
            .checked_mul(Q)
            .ok_or("overlap_bore_integer_range")?;
        let y_right = (sign * w * 4 * t * d)
            .checked_mul(P)
            .ok_or("overlap_bore_integer_range")?;
        if y_left != y_right {
            return Err("overlap_bore_vertical_conic");
        }
    }
    certify_planar_refined([2, 1], &profile)
        .map(|x| 8. * x)
        .map_err(|e| match e {
            PlanarError::Invalid => "overlap_bore_nonpositive_jacobian",
            PlanarError::Indeterminate => "overlap_bore_indeterminate_bernstein",
        })
}

/// Audit the full partition and bound the radial deviation of the materialized void.
pub fn audit_overlap(blocks: &[BracketBlock]) -> Result<BracketCertificate, &'static str> {
    if blocks.len() != 8 {
        return Err("overlap_bore_inventory");
    }
    let mut identities = BTreeMap::new();
    for (half, sign) in [("right", 1), ("left", -1)] {
        for i in 0..4 {
            identities.insert(format!("overlap_bore_{half}_sector_{i:02}_z0_8"), (sign, i));
        }
    }
    let mut seen = BTreeSet::new();
    let mut lower = f64::INFINITY;
    for block in blocks {
        if !seen.insert(&block.id) {
            return Err("overlap_bore_duplicate");
        }
        let &(sign, i) = identities.get(&block.id).ok_or("overlap_bore_identity")?;
        lower = lower.min(certify_chart(block, sign, i)?);
    }
    let incidence = enumerate_faces(blocks)?;
    if incidence.interfaces.len() != 8
        || incidence.exterior.len() != 32
        || incidence
            .interfaces
            .iter()
            .any(|f| f.first.axis != 0 || f.second.axis != 0)
    {
        return Err("overlap_bore_incidence");
    }
    // The right boundary visits (0,-s),(2,-s),(3,0),(2,s),(0,s)
    // in order. Its relative conic angle increases from -120 to +120 degrees.
    // It has x>=0 so this arc is outside the left void; rotation gives the left arc.
    // Seams at x=0 and y=0 point outwards along these axes. The two other right
    // seams are in opposite open quadrants, directed away from the right ellipse
    // (both components have the sign of the outward conic normal). They meet
    // neither each other nor any seam in x<0. The outer square and inner union
    // thus have simple boundaries. Cancellation on shared traces gives degree
    // one on the material and zero elsewhere. Positive local Jacobians count
    // all interior preimages with the same sign. Extrusion is direct on [0,8].
    let mut design_faces = Vec::new();
    for side in &incidence.exterior {
        let (sign, i) = identities[&side.block];
        let tag = match (side.axis, side.upper) {
            (2, false) => "z0",
            (2, true) => "z8",
            (1, true) => {
                if sign == 1 {
                    "bore_right"
                } else {
                    "bore_left"
                }
            }
            (1, false) => {
                if sign == 1 {
                    ["y-10", "x10", "x10", "y10"][i]
                } else {
                    ["y10", "x-10", "x-10", "y-10"][i]
                }
            }
            _ => return Err("overlap_bore_unexpected_exterior"),
        };
        design_faces.push((side.clone(), tag.into()));
    }
    // Since sqrt(3)*(s+sqrt(3)) >= 3, the radial error is bounded by
    // 2*|P²-3Q²|/(3Q²). Integer arithmetic verifies 7e-16 mm strictly.
    const {
        assert!(2 * (P * P - 3 * Q * Q).abs() * 10_000_000_000_000_000 < 7 * 3 * Q * Q);
    }
    Ok(BracketCertificate {
        incidence,
        design_faces,
        determinant_lower_bound_mm3: lower,
        radial_error_upper_mm: 7e-16,
    })
}
