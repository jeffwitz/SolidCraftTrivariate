//! Exact conic identity and continuous regularity of the eight-chart bore.
use crate::bracket::BracketBlock;
use crate::bracket_audit::BracketCertificate;
use crate::incidence::enumerate_faces;
use crate::planar::{certify_planar, PlanarError};
use std::collections::{BTreeMap, BTreeSet};

const AXES: [[i128; 2]; 4] = [[1, 0], [0, 1], [-1, 0], [0, -1]];

fn rotate([x, y]: [i128; 2], quarter: usize) -> [i128; 2] {
    match quarter {
        0 => [x, y],
        1 => [-y, x],
        2 => [-x, -y],
        _ => [y, -x],
    }
}

fn wall_point(center: [i128; 2], normal: [i128; 2]) -> [i128; 2] {
    std::array::from_fn(|i| {
        if normal[i] == 0 {
            center[i]
        } else {
            10000 * normal[i]
        }
    })
}

fn check_chart(
    block: &BracketBlock,
    center: [i128; 2],
    radius: i128,
    quarter: usize,
    half: usize,
) -> Result<f64, &'static str> {
    if block.degrees != [2, 1, 1] || block.controls.len() != 12 {
        return Err("exact_bore_control_layout");
    }
    let dir = AXES[quarter];
    let next = AXES[(quarter + 1) % 4];
    let corner = std::array::from_fn(|i| 10000 * (dir[i] + next[i]));
    let vertices = [wall_point(center, dir), corner, wall_point(center, next)];
    let (a, b) = (vertices[half], vertices[half + 1]);
    let mut profile = [[0; 3]; 6];
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
                return Err("exact_bore_extrusion");
            }
            profile[2 * i + j] = [p[0], p[1], p[3]];
        }
        let [x, y, w] = profile[2 * i].map(i128::from);
        if w != i128::from(profile[2 * i + 1][2]) {
            return Err("exact_bore_radial_weight_mismatch");
        }
        for (axis, coordinate) in [x, y].into_iter().enumerate() {
            let expected = (2 - i as i128) * a[axis] + i as i128 * b[axis];
            if 2000 * coordinate != w * expected {
                return Err("exact_bore_outer_segment");
            }
        }
    }
    // After clearing W and 1+q², both coordinate identities have degree
    // at most four in u. Five distinct exact rational nodes establish both
    // identities everywhere, including the declared angular parameterization.
    for n in 0..=4i128 {
        let basis = [(4 - n).pow(2), 2 * n * (4 - n), n * n];
        let h: [i128; 3] = std::array::from_fn(|c| {
            (0..3)
                .map(|i| basis[i] * i128::from(profile[2 * i + 1][c]))
                .sum()
        });
        let q = 4 * half as i128 + n; // q/8 is the global quarter-circle parameter.
        let arc = rotate([64 - q * q, 16 * q], quarter);
        for axis in 0..2 {
            if (1000 * h[axis] - center[axis] * h[2]) * (64 + q * q) != radius * arc[axis] * h[2] {
                return Err("exact_bore_conic_identity");
            }
        }
    }
    certify_planar([2, 1], &profile)
        .map(|bound| 8. * bound)
        .map_err(|e| match e {
            PlanarError::Invalid => "exact_bore_nonpositive_jacobian",
            PlanarError::Indeterminate => "exact_bore_indeterminate_bernstein",
        })
}

/// Audit the physical partition of a square prism with one circular void.
/// Coordinates in the recipe are integer micrometres; controls evaluate in mm.
pub fn audit_exact_bore(
    center_um: [i32; 2],
    radius_um: i32,
    blocks: &[BracketBlock],
) -> Result<BracketCertificate, &'static str> {
    let center = center_um.map(i128::from);
    let r = i128::from(radius_um);
    if r <= 0 || center.iter().any(|c| c.abs() + r >= 10000) {
        return Err("exact_bore_outside_admissible_domain");
    }
    if blocks.len() != 8 {
        return Err("exact_bore_inventory");
    }
    let mut expected = BTreeMap::new();
    for (quarter, name) in ["east", "north", "west", "south"].into_iter().enumerate() {
        for half in 0..2 {
            expected.insert(format!("exact_bore_sector_{name}_half_{half}_z0_8000"), (quarter, half));
        }
    }
    let mut seen = BTreeSet::new();
    let mut lower = f64::INFINITY;
    for b in blocks {
        let &(quarter, half) = expected.get(&b.id).ok_or("exact_bore_identity")?;
        if !seen.insert(&b.id) {
            return Err("exact_bore_duplicate_identity");
        }
        lower = lower.min(check_chart(b, center, r, quarter, half)?);
    }
    let incidence = enumerate_faces(blocks)?;
    if incidence.interfaces.len() != 8
        || incidence.exterior.len() != 32
        || incidence
            .interfaces
            .iter()
            .any(|f| f.first.axis != 0 || f.second.axis != 0)
    {
        return Err("exact_bore_incidence");
    }
    // The oriented parameter complex is an annulus. The outer boundary maps
    // once around the simple square and the inner boundary once (with reverse
    // orientation) around the exact circle. Their winding difference is one in
    // the material and zero in the void/exterior. Each interior preimage counts
    // positively by the closed Jacobian bounds; paired traces are identical.
    // The eight seam interiors are straight connectors. Check them independently
    // below to include the lower-dimensional seams in the injectivity claim.
    let mut segments = vec![];
    for q in 0..4 {
        let a = AXES[q];
        let b = AXES[(q + 1) % 4];
        let outside = wall_point(center, a).map(|x| 5 * x);
        let inside = std::array::from_fn(|i| 5 * (center[i] + r * a[i]));
        segments.push((inside, outside));
        let direction = rotate([3, 4], q);
        segments.push((
            std::array::from_fn(|i| 5 * center[i] + r * direction[i]),
            std::array::from_fn(|i| 50000 * (a[i] + b[i])),
        ));
    }
    for (i, &segment) in segments.iter().enumerate() {
        let (start, end) = segment;
        let tangent_normal: [i128; 2] = std::array::from_fn(|a| start[a] - 5 * center[a]);
        let outward: i128 = (0..2)
            .map(|a| (end[a] - start[a]) * tangent_normal[a])
            .sum();
        if outward <= 0
            || segments[i + 1..]
                .iter()
                .any(|&other| intersects(segment, other))
        {
            return Err("exact_bore_crossing_connector");
        }
    }
    let mut design_faces = vec![];
    for side in &incidence.exterior {
        let &(q, half) = expected.get(&side.block).ok_or("exact_bore_identity")?;
        let tag = match (side.axis, side.upper) {
            (1, true) => "bore",
            (1, false) => ["x10", "y10", "x-10", "y-10"][(q + half) % 4],
            (2, false) => "z0",
            (2, true) => "z8",
            _ => return Err("exact_bore_unexpected_exterior"),
        };
        design_faces.push((side.clone(), tag.to_owned()));
    }
    design_faces.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(BracketCertificate {
        incidence,
        determinant_lower_bound_mm3: lower,
        radial_error_upper_mm: 0.,
        design_faces,
    })
}

type Segment = ([i128; 2], [i128; 2]);
fn intersects((a, b): Segment, (c, d): Segment) -> bool {
    fn orient(a: [i128; 2], b: [i128; 2], p: [i128; 2]) -> i128 {
        (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
    }
    if (0..2).any(|i| a[i].max(b[i]) < c[i].min(d[i]) || c[i].max(d[i]) < a[i].min(b[i])) {
        return false;
    }
    orient(a, b, c).signum() * orient(a, b, d).signum() <= 0
        && orient(c, d, a).signum() * orient(c, d, b).signum() <= 0
}
