//! Independent exact audit for polygonal prisms with open through-holes.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    incidence::{enumerate_faces, Incidence},
    perforated_prism::{LoopProfile, PerforatedPrism},
    polygon_prism::Vertex,
    polygon_prism_audit::{
        check_quad_jacobian, edge_key, expected_controls, fields_for_index, next_down,
        parse_triangle_id, validate_triangle_coverage, DirectedBoundary, PolygonCertificate, Tri,
    },
};
use std::collections::{BTreeMap, BTreeSet};

type ValidatedProfiles = (Vec<LoopProfile>, Vec<Vertex>, i64, DirectedBoundary);

fn signed_area2(vertices: &[Vertex]) -> i64 {
    (0..vertices.len())
        .map(|i| {
            let a = vertices[i].xy_mm;
            let b = vertices[(i + 1) % vertices.len()].xy_mm;
            i64::from(a[0]) * i64::from(b[1]) - i64::from(a[1]) * i64::from(b[0])
        })
        .sum()
}

fn strict_inside(point2: [i64; 2], vertices: &[Vertex]) -> bool {
    let mut inside = false;
    for i in 0..vertices.len() {
        let a = vertices[i].xy_mm.map(i64::from);
        let b = vertices[(i + 1) % vertices.len()].xy_mm.map(i64::from);
        let cross = (b[0] - a[0]) * (point2[1] - 2 * a[1]) - (b[1] - a[1]) * (point2[0] - 2 * a[0]);
        if cross == 0
            && (2 * a[0]).min(2 * b[0]) <= point2[0]
            && point2[0] <= (2 * a[0]).max(2 * b[0])
            && (2 * a[1]).min(2 * b[1]) <= point2[1]
            && point2[1] <= (2 * a[1]).max(2 * b[1])
        {
            return false;
        }
        if (2 * a[1] > point2[1]) != (2 * b[1] > point2[1]) {
            let lhs = (point2[0] - 2 * a[0]) * (b[1] - a[1]);
            let rhs = (point2[1] - 2 * a[1]) * (b[0] - a[0]);
            if (b[1] > a[1]) == (lhs < rhs) {
                inside = !inside;
            }
        }
    }
    inside
}

fn validate_profiles(recipe: &PerforatedPrism) -> Result<ValidatedProfiles, &'static str> {
    if recipe.holes.len() > 4
        || !(1..=64).contains(&recipe.height_mm)
        || recipe.top_shift_mm.iter().any(|x| !(-32..=32).contains(x))
    {
        return Err("perforated_range");
    }
    let mut loops = vec![recipe.outer.clone()];
    loops.extend(recipe.holes.clone());
    let mut loop_ids = BTreeSet::new();
    let mut vertex_ids = BTreeMap::new();
    let mut total = 0;
    for (li, lp) in loops.iter().enumerate() {
        if lp.id.is_empty()
            || !lp
                .id
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || c == b'_')
            || !loop_ids.insert(lp.id.as_str())
            || lp.vertices.len() < 3
        {
            return Err("perforated_loop");
        }
        total += lp.vertices.len();
        for (vi, v) in lp.vertices.iter().enumerate() {
            if v.id.is_empty()
                || !v.id.bytes().all(|c| c.is_ascii_alphanumeric() || c == b'_')
                || v.xy_mm.iter().any(|x| !(-64..=64).contains(x))
                || vertex_ids
                    .insert(
                        v.id.as_str(),
                        vi + loops[..li].iter().map(|x| x.vertices.len()).sum::<usize>(),
                    )
                    .is_some()
            {
                return Err("perforated_vertex");
            }
            if crate::polygon_prism_audit::cross(
                lp.vertices[vi].xy_mm,
                lp.vertices[(vi + 1) % lp.vertices.len()].xy_mm,
                lp.vertices[(vi + 2) % lp.vertices.len()].xy_mm,
            ) == 0
            {
                return Err("perforated_collinear_adjacent_edges");
            }
            for j in vi + 1..lp.vertices.len() {
                if j == vi + 1 || (vi == 0 && j + 1 == lp.vertices.len()) {
                    continue;
                }
                if crate::polygon_prism_audit::segments_contact(
                    v.xy_mm,
                    lp.vertices[(vi + 1) % lp.vertices.len()].xy_mm,
                    lp.vertices[j].xy_mm,
                    lp.vertices[(j + 1) % lp.vertices.len()].xy_mm,
                ) {
                    return Err("perforated_loop_contact");
                }
            }
        }
        if signed_area2(&lp.vertices) == 0 {
            return Err("perforated_zero_area");
        }
    }
    if !(3..=32).contains(&total) {
        return Err("perforated_vertex_count");
    }
    for i in 0..loops.len() {
        for j in i + 1..loops.len() {
            for a in 0..loops[i].vertices.len() {
                for b in 0..loops[j].vertices.len() {
                    if crate::polygon_prism_audit::segments_contact(
                        loops[i].vertices[a].xy_mm,
                        loops[i].vertices[(a + 1) % loops[i].vertices.len()].xy_mm,
                        loops[j].vertices[b].xy_mm,
                        loops[j].vertices[(b + 1) % loops[j].vertices.len()].xy_mm,
                    ) {
                        return Err("perforated_boundary_contact");
                    }
                }
            }
        }
    }
    let mut normalized = loops;
    for (i, lp) in normalized.iter_mut().enumerate() {
        let area = signed_area2(&lp.vertices);
        if (area > 0) != (i == 0) {
            lp.vertices.reverse();
        }
    }
    for h in 1..normalized.len() {
        for v in &normalized[h].vertices {
            if !strict_inside(
                [2 * i64::from(v.xy_mm[0]), 2 * i64::from(v.xy_mm[1])],
                &normalized[0].vertices,
            ) {
                return Err("perforated_hole_outside");
            }
        }
        for other in 1..normalized.len() {
            if other != h
                && strict_inside(
                    [
                        2 * i64::from(normalized[h].vertices[0].xy_mm[0]),
                        2 * i64::from(normalized[h].vertices[0].xy_mm[1]),
                    ],
                    &normalized[other].vertices,
                )
            {
                return Err("perforated_nested_holes");
            }
        }
    }
    let area2 = normalized
        .iter()
        .map(|lp| signed_area2(&lp.vertices))
        .sum::<i64>();
    if area2 <= 0 {
        return Err("perforated_empty_material");
    }
    let vertices = normalized
        .iter()
        .flat_map(|lp| lp.vertices.iter().cloned())
        .collect::<Vec<_>>();
    let mut boundary = DirectedBoundary::new();
    let mut offset = 0;
    for lp in &normalized {
        for i in 0..lp.vertices.len() {
            let a = offset + i;
            let b = offset + (i + 1) % lp.vertices.len();
            let (lo, hi) = if lp.vertices[i].id < lp.vertices[(i + 1) % lp.vertices.len()].id {
                (
                    &lp.vertices[i].id,
                    &lp.vertices[(i + 1) % lp.vertices.len()].id,
                )
            } else {
                (
                    &lp.vertices[(i + 1) % lp.vertices.len()].id,
                    &lp.vertices[i].id,
                )
            };
            let tag = format!("edge/{}/{lo}/{hi}", lp.id);
            if boundary.insert(edge_key(a, b), ((a, b), tag)).is_some() {
                return Err("perforated_duplicate_boundary_edge");
            }
        }
        offset += lp.vertices.len();
    }
    Ok((normalized, vertices, area2, boundary))
}

/// Audit a perforated extrusion without relying on the constrained triangulator.
pub fn audit_perforated_prism(
    recipe: &PerforatedPrism,
    blocks: &[BracketBlock],
) -> Result<PolygonCertificate, &'static str> {
    let (_loops, vertices, area2, boundary) = validate_profiles(recipe)?;
    let idmap = vertices
        .iter()
        .enumerate()
        .map(|(i, v)| (v.id.as_str(), i))
        .collect::<BTreeMap<_, _>>();
    let triangles_expected = vertices.len() + 2 * recipe.holes.len() - 2;
    if blocks.len() != 3 * triangles_expected {
        return Err("perforated_block_count");
    }
    let mut unique = BTreeSet::new();
    let mut groups = BTreeMap::<[usize; 3], BTreeMap<usize, &BracketBlock>>::new();
    let mut min_det = i128::MAX;
    for block in blocks {
        if !unique.insert(&block.id) {
            return Err("perforated_duplicate_block");
        }
        if block.degrees != [1, 1, 1] || block.controls.len() != 8 {
            return Err("perforated_altered_block");
        }
        let (tri, at) = parse_triangle_id(&block.id, &idmap)?;
        if tri
            .iter()
            .map(|i| fields_for_index(&idmap, *i))
            .min()
            .unwrap()
            != fields_for_index(&idmap, tri[0])
            || block.id
                != format!(
                    "tri/{}/{}/{}/{}",
                    fields_for_index(&idmap, tri[0]),
                    fields_for_index(&idmap, tri[1]),
                    fields_for_index(&idmap, tri[2]),
                    fields_for_index(&idmap, at)
                )
        {
            return Err("perforated_block_id");
        }
        let signed = crate::polygon_prism_audit::cross(
            vertices[tri[0]].xy_mm,
            vertices[tri[1]].xy_mm,
            vertices[tri[2]].xy_mm,
        );
        if signed <= 0 {
            return Err("perforated_triangle_orientation");
        }
        if block.controls
            != expected_controls(tri, at, &vertices, recipe.height_mm, recipe.top_shift_mm)
                .ok_or("perforated_block_id")?
        {
            return Err("perforated_altered_block");
        }
        min_det = min_det.min(check_quad_jacobian(block)?);
        if groups.entry(tri).or_default().insert(at, block).is_some() {
            return Err("perforated_duplicate_patch");
        }
    }
    let mut triangles = Vec::new();
    for (ids, patches) in groups {
        if patches.len() != 3 || ids.iter().any(|i| !patches.contains_key(i)) {
            return Err("perforated_missing_patch");
        }
        triangles.push(Tri { ids });
    }
    validate_triangle_coverage(&triangles, &vertices, triangles_expected, area2, &boundary)?;
    let incidence: Incidence = enumerate_faces(blocks)?;
    let mut design_faces = Vec::new();
    for side in &incidence.exterior {
        let parts = side.block.split('/').collect::<Vec<_>>();
        if parts.len() != 5 || parts[0] != "tri" {
            return Err("perforated_boundary_provenance");
        }
        let tri = [
            *idmap
                .get(parts[1])
                .ok_or("perforated_boundary_provenance")?,
            *idmap
                .get(parts[2])
                .ok_or("perforated_boundary_provenance")?,
            *idmap
                .get(parts[3])
                .ok_or("perforated_boundary_provenance")?,
        ];
        let at = *idmap
            .get(parts[4])
            .ok_or("perforated_boundary_provenance")?;
        let place = tri
            .iter()
            .position(|i| *i == at)
            .ok_or("perforated_boundary_provenance")?;
        let prev = tri[(place + 2) % 3];
        let next = tri[(place + 1) % 3];
        let feature = match (side.axis, side.upper) {
            (2, false) => "start".to_owned(),
            (2, true) => "end".to_owned(),
            (0, false) => boundary
                .get(&edge_key(prev, at))
                .map(|x| x.1.clone())
                .ok_or("perforated_unexpected_exterior")?,
            (1, false) => boundary
                .get(&edge_key(at, next))
                .map(|x| x.1.clone())
                .ok_or("perforated_unexpected_exterior")?,
            _ => return Err("perforated_unexpected_exterior"),
        };
        design_faces.push((side.clone(), feature));
    }
    design_faces.sort_by(|a, b| a.0.cmp(&b.0));
    let volume = area2
        .checked_mul(3)
        .and_then(|x| x.checked_mul(i64::from(recipe.height_mm)))
        .ok_or("perforated_volume_overflow")?;
    let det_num = min_det
        .checked_mul(i128::from(recipe.height_mm))
        .ok_or("perforated_jacobian_overflow")?;
    Ok(PolygonCertificate {
        certificate: BracketCertificate {
            incidence,
            determinant_lower_bound_mm3: next_down(det_num as f64 / 36.),
            radial_error_upper_mm: 0.,
            design_faces,
        },
        volume_sixths_mm3: volume,
        projected_area_twice_mm2: area2,
    })
}
