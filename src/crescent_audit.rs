//! Independent certificate for the fixed intersecting through-bore and offset pocket.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    incidence::enumerate_faces,
    overlap_bore_audit::audit_overlap,
    planar::{certify_planar_refined, PlanarError},
};
use std::collections::{BTreeMap, BTreeSet};

const P: i128 = 1_732_050_807_568_877;
const Q: i128 = 1_000_000_000_000_000;

/// Restore the planar section to its reference extrusion after checking all
/// stored lower and upper controls against the physical half-mm slice planes.
fn reference_extrusion(
    block: &BracketBlock,
    bounds: [i32; 2],
) -> Result<BracketBlock, &'static str> {
    if !matches!(block.degrees, [2, 1, 1] | [2, 2, 1])
        || block.controls.len() != 6 * (block.degrees[1] + 1)
    {
        return Err("crescent_control_layout");
    }
    let mut restored = block.clone();
    for pair in restored.controls.as_chunks_mut::<2>().0 {
        if pair[0][..2] != pair[1][..2] || pair[0][3] != pair[1][3] || pair[0][3] <= 0 {
            return Err("crescent_extrusion");
        }
        for (k, p) in pair.iter_mut().enumerate() {
            if 2 * i128::from(p[2]) != i128::from(bounds[k]) * i128::from(p[3]) {
                return Err("crescent_layer_height");
            }
            p[2] = (8 * k as i64)
                .checked_mul(p[3])
                .ok_or("crescent_height_range")?;
        }
    }
    Ok(restored)
}

fn conic_edge(
    nodes: [[i128; 3]; 3],
    start: i128,
    end: i128,
    reflected_x: bool,
) -> Result<(), &'static str> {
    // t varies from start/6 to end/6. Five exact nodes prove the complete
    // degree-four identity after clearing the two quadratic denominators.
    for n in 0..=4i128 {
        let t = (4 - n) * start + n * end;
        let d = 24;
        let w = d * d + 3 * t * t;
        let sigma = [(4 - n).pow(2), 2 * n * (4 - n), n * n];
        let [x, y, weight]: [i128; 3] =
            std::array::from_fn(|k| (0..3).map(|i| sigma[i] * nodes[i][k]).sum());
        // Left branch is the pi-rotation of the reference right arc.
        // The hidden right branch is the reflection in x=1 of its visible branch.
        let right = w + 2 * (d * d - 3 * t * t);
        let target_x = if reflected_x { 2 * w - right } else { -right };
        let y_sign = if reflected_x { 1 } else { -1 };
        if x * w != weight * target_x || ratio(y, weight) != ratio(y_sign * 4 * t * d * P, w * Q) {
            return Err("crescent_boundary_conic");
        }
    }
    Ok(())
}

fn ratio(n: i128, d: i128) -> (i128, i128) {
    let (mut a, mut b) = (n.abs(), d);
    while b != 0 {
        (a, b) = (b, a % b);
    }
    (n / a, d / a)
}

fn linear(edge: [[i128; 3]; 3]) -> Result<(), &'static str> {
    if (0..3).any(|k| 2 * edge[1][k] != edge[0][k] + edge[2][k]) {
        return Err("crescent_noncanonical_straight_seam");
    }
    Ok(())
}

fn certify_crescent(block: &BracketBlock, lower: bool, patch: &str) -> Result<f64, &'static str> {
    if block.degrees != [2, 2, 1] {
        return Err("crescent_biquadratic_layout");
    }
    let profile: Vec<_> = block
        .controls
        .as_chunks::<2>()
        .0
        .iter()
        .map(|p| [p[0][0], p[0][1], p[0][3]])
        .collect();
    let h = |u: usize, v: usize| {
        let [x, y, w] = profile[3 * (if lower { 2 - u } else { u }) + v].map(i128::from);
        [x, if lower { -y } else { y }, w]
    };
    let left = std::array::from_fn(|v| h(0, v));
    let right = std::array::from_fn(|v| h(2, v));
    let bottom = std::array::from_fn(|u| h(u, 0));
    let top = std::array::from_fn(|u| h(u, 2));
    match patch {
        "a" => {
            conic_edge(left, 0, -2, false)?;
            linear(right)?;
            linear(top)?;
        }
        "b" => {
            conic_edge(right, 0, 1, true)?;
            linear(left)?;
            linear(top)?;
        }
        "c" => {
            conic_edge(right, 1, 2, true)?;
            conic_edge(top, -2, -6, false)?;
            linear(left)?;
            linear(bottom)?;
        }
        _ => return Err("crescent_patch_identity"),
    }
    if patch != "c" {
        let e = if patch == "a" { bottom[2] } else { bottom[0] };
        if e[0] != -2 * e[2] {
            return Err("crescent_baseline_midpoint");
        }
        if bottom.iter().any(|p| p[1] != 0)
            || bottom
                .windows(2)
                .any(|p| p[0][0] * p[1][2] >= p[1][0] * p[0][2])
            || bottom[0][0] < -3 * bottom[0][2]
            || bottom[2][0] > -bottom[2][2]
        {
            return Err("crescent_symmetry_edge");
        }
    }
    // D,E,F meet the physical boundary, while O lies strictly inside the
    // triangle of the three cuts. The three straight incident seams have
    // different directions and meet only at O. With positive chart
    // orientation, cancellation of matching seams gives winding number
    // one in this half-crescent and zero outside.
    let interior = match patch {
        "a" => h(2, 2),
        "b" => h(0, 2),
        _ => h(0, 0),
    };
    let [x, y, w] = interior;
    // Three exact oriented half-plane tests for D=(-2,s), E=(-2,0),
    // F=(-9/13,8s/13) with s=P/Q.
    if x <= -2 * w || 17 * Q * y <= 8 * P * (x + 2 * w) || 17 * Q * y >= P * (7 * w - 5 * x) {
        return Err("crescent_center_outside_cut_triangle");
    }
    certify_planar_refined([2, 2], &profile).map_err(|e| match e {
        PlanarError::Invalid => "crescent_nonpositive_jacobian",
        PlanarError::Indeterminate => "crescent_indeterminate_bernstein",
    })
}

pub fn audit_crescent_pocket(
    depth_half_mm: i32,
    blocks: &[BracketBlock],
) -> Result<BracketCertificate, &'static str> {
    if !(0..=16).contains(&depth_half_mm) {
        return Err("crescent_depth_range");
    }
    let cut = 16 - depth_half_mm;
    let mut lower_ring = Vec::new();
    let mut upper_ring = Vec::new();
    let mut crescent_count = 0;
    let mut seen = BTreeSet::new();
    let mut minimum = f64::INFINITY;
    let mut features = BTreeMap::new();
    for block in blocks {
        if !seen.insert(&block.id) {
            return Err("crescent_duplicate_identity");
        }
        if let Some(name) = block.id.strip_prefix("g10_crescent_") {
            let (half, patch) = name.split_once('_').ok_or("crescent_identity")?;
            if !matches!(half, "upper" | "lower") || cut == 0 {
                return Err("crescent_unexpected_half");
            }
            let lower = half == "lower";
            let source = reference_extrusion(block, [0, cut])?;
            minimum =
                minimum.min(certify_crescent(&source, lower, patch)? * (f64::from(cut) * 0.5));
            if patch != "a" {
                features.insert((block.id.clone(), 0, !lower), "bore_right".to_owned());
            }
            features.insert((block.id.clone(), 2, false), "z0".into());
            features.insert(
                (block.id.clone(), 2, true),
                if depth_half_mm == 0 {
                    "z8"
                } else {
                    "pocket_floor"
                }
                .into(),
            );
            crescent_count += 1;
            continue;
        }
        let (base, upper) = if let Some(s) = block.id.strip_prefix("g10_ring_upper_") {
            (s.strip_suffix("_split_8").ok_or("crescent_layer_id")?, true)
        } else if let Some(s) = block.id.strip_prefix("g10_ring_lower_") {
            (
                s.strip_suffix("_z0_split").ok_or("crescent_layer_id")?,
                false,
            )
        } else {
            return Err("crescent_unknown_block");
        };
        let z = if upper { [cut, 16] } else { [0, cut] };
        if z[0] == z[1] {
            return Err("crescent_empty_layer");
        }
        let mut source = reference_extrusion(block, z)?;
        source.id = format!("overlap_bore_{base}_z0_8");
        if upper {
            upper_ring.push(source);
        } else {
            lower_ring.push(source);
        }
    }
    if crescent_count != if cut > 0 { 6 } else { 0 } {
        return Err("crescent_inventory");
    }
    for (ring, upper, bounds) in [
        (&upper_ring, true, [cut, 16]),
        (&lower_ring, false, [0, cut]),
    ] {
        if bounds[0] == bounds[1] {
            continue;
        }
        let cert = audit_overlap(ring)?;
        minimum = minimum
            .min(cert.determinant_lower_bound_mm3 * (f64::from(bounds[1] - bounds[0]) / 16.));
        for (side, tag) in cert.design_faces {
            if side.axis == 2
                && ((side.upper && bounds[1] != 16) || (!side.upper && bounds[0] != 0))
            {
                continue;
            }
            if tag == "bore_left" && !upper {
                continue;
            }
            let base = side
                .block
                .strip_prefix("overlap_bore_")
                .and_then(|s| s.strip_suffix("_z0_8"))
                .ok_or("crescent_source_id")?;
            let layer = if upper { "upper" } else { "lower" };
            let feature = if tag == "bore_left" {
                "pocket_wall".to_owned()
            } else {
                tag
            };
            features.insert(
                (
                    format!(
                        "g10_ring_{layer}_{base}_{}",
                        if upper { "split_8" } else { "z0_split" }
                    ),
                    side.axis,
                    side.upper,
                ),
                feature,
            );
        }
    }
    let incidence = enumerate_faces(blocks)?;
    let expected = match depth_half_mm {
        0 => 20,
        16 => 8,
        _ => 36,
    };
    if incidence.interfaces.len() != expected || features.len() != incidence.exterior.len() {
        return Err("crescent_incidence");
    }
    let mut design_faces = Vec::new();
    for side in &incidence.exterior {
        let tag = features
            .remove(&(side.block.clone(), side.axis, side.upper))
            .ok_or("crescent_unexpected_side")?;
        design_faces.push((side.clone(), tag));
    }
    Ok(BracketCertificate {
        incidence,
        design_faces,
        determinant_lower_bound_mm3: minimum.next_down(),
        radial_error_upper_mm: 7e-16,
    })
}
