//! Certify the stored charts of the equal-radius orthogonal through-bore junction.
//!
//! The exact nonnegative cylinder polynomials confine all charts to target material.
//! The oriented volume boundaries cancel pairwise on interior interfaces. Physical
//! exterior traces lie on the target stock planes or the actual cylinder surfaces.
//! Thus the sum of local degrees is constant on the connected interior of the target.
//! A unique affine brick occupies a neighborhood of each stock corner, fixing that
//! degree at one. Strict positive Jacobians then exclude missing or multiply covered
//! open volumes. No quadrature or point cloud is used for this conclusion.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    incidence::{enumerate_faces, BlockSide},
    orthogonal_sketch::{OrthogonalSketch, SweepPartition},
    spatial_certificate::{certify_spatial, CertificateBudget, CertificateStatus},
};
use std::collections::BTreeMap;

fn binomial(d: usize, k: usize) -> i128 {
    if d == 2 && k == 1 {
        2
    } else {
        1
    }
}
fn axis_index(i: usize, d: [usize; 3]) -> [usize; 3] {
    [
        i / ((d[1] + 1) * (d[2] + 1)),
        i / (d[2] + 1) % (d[1] + 1),
        i % (d[2] + 1),
    ]
}
fn cylinder_exclusion(
    b: &BracketBlock,
    axis: usize,
    center: [i32; 3],
    r: i32,
) -> Result<(), &'static str> {
    let d = b.degrees;
    let dims = d.map(|v| 2 * v + 1);
    let n = dims.iter().product();
    let mut coeff = vec![0i128; n];
    let a = (axis + 1) % 3;
    let c = (axis + 2) % 3;
    let mut terms = Vec::new();
    for (i, h) in b.controls.iter().enumerate() {
        let q = axis_index(i, d);
        let weight = (0..3).map(|k| binomial(d[k], q[k])).product::<i128>();
        let w = i128::from(h[3]);
        // The validated design window and normalized constructor nets keep these products in i128.
        let p = [
            i128::from(h[a]) - i128::from(center[a]) * w,
            i128::from(h[c]) - i128::from(center[c]) * w,
            i128::from(r) * w,
        ];
        terms.push((q, p, weight));
    }
    for &(i, p, v) in &terms {
        for &(j, q, w) in &terms {
            let k = ((i[0] + j[0]) * dims[1] + i[1] + j[1]) * dims[2] + i[2] + j[2];
            let prod = |k: usize| p[k].checked_mul(q[k]).ok_or("steinmetz_exclusion_overflow");
            let t = prod(0)?
                .checked_add(prod(1)?)
                .and_then(|s| s.checked_sub(prod(2).ok()?))
                .and_then(|s| s.checked_mul(v * w))
                .ok_or("steinmetz_exclusion_overflow")?;
            coeff[k] = coeff[k]
                .checked_add(t)
                .ok_or("steinmetz_exclusion_overflow")?;
        }
    }
    if coeff.iter().any(|v| *v < 0) {
        return Err("steinmetz_cylinder_exclusion_indeterminate");
    }
    Ok(())
}
fn corner_degree(blocks: &[BracketBlock], bounds: [[i32; 2]; 3]) -> Result<(), &'static str> {
    for octant in 0..8 {
        let p: [i32; 3] = std::array::from_fn(|k| bounds[k][(octant >> k) & 1]);
        let mut candidates = blocks.iter().filter(|b| {
            (0..3).all(|k| {
                b.controls
                    .iter()
                    .any(|h| i128::from(h[k]) == i128::from(p[k]) * i128::from(h[3]))
            })
        });
        let Some(b) = candidates.next() else {
            return Err("steinmetz_missing_corner");
        };
        if candidates.next().is_some()
            || b.degrees != [1; 3]
            || b.controls.iter().any(|h| h[3] != 1)
        {
            return Err("steinmetz_ambiguous_corner_degree");
        }
        let mut used = [false; 3];
        for k in 0..3 {
            let p0 = b.controls[0][k];
            let mut direction = None;
            for (axis, step) in [4, 2, 1].into_iter().enumerate() {
                if b.controls[step][k] != p0 && direction.replace((axis, step)).is_some() {
                    return Err("steinmetz_corner_not_affine");
                }
            }
            let (axis, step) = direction.ok_or("steinmetz_degenerate_corner")?;
            if used[axis]
                || b.controls.iter().enumerate().any(|(i, h)| {
                    h[k] != if i & step != 0 {
                        b.controls[step][k]
                    } else {
                        p0
                    }
                })
            {
                return Err("steinmetz_corner_not_affine");
            }
            used[axis] = true;
        }
    }
    Ok(())
}
fn face_on_plane(b: &BracketBlock, s: &BlockSide, k: usize, plane: i32) -> bool {
    b.controls
        .iter()
        .enumerate()
        .filter(|(i, _)| {
            axis_index(*i, b.degrees)[s.axis] == usize::from(s.upper) * b.degrees[s.axis]
        })
        .all(|(_, h)| i128::from(h[k]) == i128::from(plane) * i128::from(h[3]))
}

pub(crate) fn try_build(j: &OrthogonalSketch) -> Result<Option<SweepPartition>, &'static str> {
    let [a, b] = j.bores.as_slice() else {
        return Ok(None);
    };
    if a.axis == b.axis
        || a.axis >= 3
        || b.axis >= 3
        || !j.fillets.is_empty()
        || a.radius != b.radius
    {
        return Ok(None);
    }
    let normal = 3 - a.axis - b.axis;
    if a.center[normal] != b.center[normal] {
        return Ok(None);
    }
    if a.radius <= 0
        || a.radius > 1024
        || a.center
            .into_iter()
            .chain(b.center)
            .any(|x| !(-4096..=4096).contains(&x))
    {
        return Err("invalid_steinmetz_bore");
    }
    let bounds = std::array::from_fn(|k| {
        if k == 2 {
            [0, j.height]
        } else {
            let pts = j.outer.vertices.iter().map(|v| v[k]);
            [pts.clone().min().unwrap(), pts.max().unwrap()]
        }
    });
    let size = bounds.map(|v| i64::from(v[1] - v[0]));
    if crate::sketch_layers::area2(&j.outer) != 2 * size[0] * size[1] {
        return Ok(None);
    }
    let mut center = a.center;
    center[a.axis] = b.center[a.axis];
    for p in [a.center, b.center] {
        if (0..3).any(|k| p[k] <= bounds[k][0] || p[k] >= bounds[k][1]) {
            return Err("steinmetz_axis_seed_outside_stock");
        }
    }
    let guard = a.radius + 1;
    let axes = [a.axis, b.axis];
    let mut sorted_axes = axes;
    sorted_axes.sort();
    let mut blocks = crate::steinmetz_frame::intersecting_cylinders(
        bounds,
        center,
        a.radius,
        guard,
        sorted_axes,
        "junction",
    )?;
    blocks.sort_by(|u, v| u.id.cmp(&v.id));
    if blocks.len() != 40 {
        return Err("steinmetz_reported_partition_size");
    }
    let mut det = f64::INFINITY;
    for block in &blocks {
        if block.controls.iter().any(|h| {
            h[3] <= 0
                || (0..3).any(|k| {
                    let x = i128::from(h[k]);
                    let w = i128::from(h[3]);
                    x < i128::from(bounds[k][0]) * w || x > i128::from(bounds[k][1]) * w
                })
        }) {
            return Err("steinmetz_stock_convex_hull_violation");
        }
        for ax in axes {
            cylinder_exclusion(block, ax, center, a.radius)?;
        }
        let cert = certify_spatial(block.degrees, &block.controls, CertificateBudget::default());
        if cert.status != CertificateStatus::Valid {
            return Err("steinmetz_spatial_certificate_incomplete");
        }
        det = det.min(
            cert.determinant_lower_bound_mm3
                .ok_or("steinmetz_positive_bound_missing")?,
        );
    }
    corner_degree(&blocks, bounds)?;
    let incidence = enumerate_faces(&blocks)?;
    for link in &incidence.interfaces {
        let side = |s: &BlockSide| s.upper ^ (s.axis == 1);
        let o = link.orientation;
        if !(side(&link.first) ^ side(&link.second) ^ o.swap ^ o.flip_first ^ o.flip_second) {
            return Err("steinmetz_oriented_interfaces_do_not_cancel");
        }
    }
    let index: BTreeMap<_, _> = blocks.iter().map(|c| (&c.id, c)).collect();
    let mut design_faces = Vec::new();
    for s in &incidence.exterior {
        let block = index[&s.block];
        let mut label = None;
        for (k, ends) in bounds.iter().enumerate() {
            for (i, &h) in ends.iter().enumerate() {
                if face_on_plane(block, s, k, h) {
                    if label.is_some() {
                        return Err("steinmetz_ambiguous_plane");
                    }
                    label = Some(if k == 2 {
                        ["start", "end"][i].to_owned()
                    } else {
                        format!("outer/{}", j.outer.id)
                    });
                }
            }
        }
        if label.is_none() {
            for feature in [a, b] {
                if crate::steinmetz_audit::certify_boundary(
                    block,
                    s,
                    feature.axis,
                    center,
                    feature.radius,
                )? {
                    if label.is_some() {
                        return Err("steinmetz_ambiguous_cylinder");
                    }
                    label = Some(format!("bore/{}/wall", feature.id));
                }
            }
        }
        design_faces.push((
            s.clone(),
            label.ok_or("steinmetz_unclassified_exterior_face")?,
        ));
    }
    let box_volume = size.iter().product::<i64>();
    let r = i64::from(a.radius);
    let planes = std::array::from_fn(|k| vec![bounds[k][0], center[k], bounds[k][1]]);
    Ok(Some(SweepPartition {
        compact: false,
        blocks,
        certificate: BracketCertificate {
            incidence,
            design_faces,
            determinant_lower_bound_mm3: det,
            radial_error_upper_mm: 0.,
        },
        twice_constant: 2 * box_volume,
        four_pi: -4 * r * r * (size[a.axis] + size[b.axis]),
        planes,
        guards: Vec::new(),
        junction: Some((center, sorted_axes, a.radius, guard)),
    }))
}
