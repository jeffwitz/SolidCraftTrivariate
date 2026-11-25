//! Independent continuous audit for polygonal prisms with exact circular bores.
//!
//! The exterior material is audited as a G18 polygon with disjoint guard
//! diamonds removed. Each guard diamond is then filled by eight rational
//! sectors. Exact coefficient reconstruction proves the circle identity and
//! the half-interval de Casteljau restrictions. Outward-rounded Bernstein
//! bounds prove a positive planar Jacobian on every closed sector.
//!
//! The guard diamond strictly contains its circle because `guard >= 2*r` and
//! `max(|x|+|y|)` on the radius-r disk is `sqrt(2)*r`. Both curves are convex.
//! Their cardinal and half-edge endpoints occur in the same cyclic order;
//! the straight connector segments therefore cannot cross, and lie outside
//! the disk except at their circle endpoints. Thus each sector boundary is a
//! Jordan curve. Positive Jacobians, exact whole-face gluing, the validated
//! G18 exterior partition, and the one-to-one outer and bore boundary cycles
//! give degree one and global coverage/injectivity of this fixed template.
//! This proof does not cover arbitrary partial intersections or other charts.

use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    circular_polygon::{CircularPolygonBore, CircularPolygonPrism},
    incidence::{enumerate_faces, BlockSide, Incidence},
    perforated_prism::{LoopProfile, PerforatedPrism},
    perforated_prism_audit::audit_perforated_prism,
    planar::{certify_planar_refined, PlanarError},
    polygon_prism::Vertex,
    polygon_prism_audit::next_down,
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq)]
pub struct CircularPolygonCertificate {
    pub certificate: BracketCertificate,
    /// Twice the signed area of the input outer polygon, in mm^2.
    pub outer_area_twice_mm2: i64,
    /// Exact sum of squared bore radii, in mm^2.
    pub sum_squared_radii_mm2: i64,
}

const QUADRANTS: [&str; 4] = ["east", "north", "west", "south"];

fn add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("circular_polygon_integer_overflow")
}

fn mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("circular_polygon_integer_overflow")
}

fn checked_i64(value: i128) -> Result<i64, &'static str> {
    i64::try_from(value).map_err(|_| "circular_polygon_integer_overflow")
}

fn valid_id(id: &str) -> bool {
    !id.is_empty() && id.len() <= 18 && id.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_')
}

fn guard_profiles(recipe: &CircularPolygonPrism) -> Result<Vec<LoopProfile>, &'static str> {
    if !(1..=4).contains(&recipe.bores.len())
        || recipe.outer.vertices.len() + 4 * recipe.bores.len() > 32
        || !(1..=64).contains(&recipe.height_mm)
        || recipe.top_shift_mm.iter().any(|x| !(-32..=32).contains(x))
    {
        return Err("circular_polygon_recipe_range");
    }
    let mut bore_ids = BTreeSet::new();
    let mut identities = BTreeSet::new();
    if recipe.outer.id.is_empty() || !identities.insert(recipe.outer.id.clone()) {
        return Err("circular_polygon_identity_collision");
    }
    for vertex in &recipe.outer.vertices {
        if vertex.id.is_empty() || !identities.insert(vertex.id.clone()) {
            return Err("circular_polygon_identity_collision");
        }
    }
    let mut loops = Vec::with_capacity(recipe.bores.len());
    for bore in &recipe.bores {
        if !valid_id(&bore.id) || !bore_ids.insert(bore.id.clone()) {
            return Err("circular_polygon_bore_identity");
        }
        if !identities.insert(bore.id.clone()) {
            return Err("circular_polygon_identity_collision");
        }
        let r = i64::from(bore.radius_mm);
        let g = i64::from(bore.guard_mm);
        if r < 1 || g < 2 * r {
            return Err("circular_polygon_guard_radius");
        }
        let loop_id = format!("guard_{}", bore.id);
        if !identities.insert(loop_id.clone()) {
            return Err("circular_polygon_identity_collision");
        }
        let [cx, cy] = bore.center_mm.map(i64::from);
        let points = [[cx + g, cy], [cx, cy + g], [cx - g, cy], [cx, cy - g]];
        let mut vertices = Vec::with_capacity(4);
        for (point, suffix) in points.into_iter().zip(["e", "n", "w", "s"]) {
            if point.iter().any(|x| !(-64..=64).contains(x)) {
                return Err("circular_polygon_guard_outside_domain");
            }
            let id = format!("g_{}_{}", bore.id, suffix);
            if !identities.insert(id.clone()) {
                return Err("circular_polygon_identity_collision");
            }
            vertices.push(Vertex {
                id,
                xy_mm: [point[0] as i32, point[1] as i32],
            });
        }
        loops.push(LoopProfile {
            id: loop_id,
            vertices,
        });
    }
    Ok(loops)
}

/// Canonical profiles in a counterclockwise cycle E->N->W->S->E.
/// Values are homogeneous controls at common gauge two (the weight-two
/// endpoints of the standard rational circle are retained exactly).
fn profiles(bore: &CircularPolygonBore, quadrant: usize) -> ([[i128; 3]; 3], [[i128; 3]; 3]) {
    let g = i128::from(bore.guard_mm);
    let r = i128::from(bore.radius_mm);
    let outer = [[2 * g, 0, 2], [g, g, 2], [0, 2 * g, 2]];
    let circle = [[2 * r, 0, 2], [2 * r, 2 * r, 2], [0, 4 * r, 4]];
    match quadrant {
        0 => (outer, circle),
        1 => (
            [[0, 2 * g, 2], [-g, g, 2], [-2 * g, 0, 2]],
            [[0, 4 * r, 4], [-2 * r, 2 * r, 2], [-2 * r, 0, 2]],
        ),
        2 => (
            [[-2 * g, 0, 2], [-g, -g, 2], [0, -2 * g, 2]],
            [[-2 * r, 0, 2], [-2 * r, -2 * r, 2], [0, -4 * r, 4]],
        ),
        _ => (
            [[0, -2 * g, 2], [g, -g, 2], [2 * g, 0, 2]],
            [[0, -4 * r, 4], [2 * r, -2 * r, 2], [2 * r, 0, 2]],
        ),
    }
}

fn expected_ring_block(
    bore: &CircularPolygonBore,
    height: i32,
    shift: [i32; 2],
    quadrant: usize,
    half: usize,
) -> Result<BracketBlock, &'static str> {
    let (outer, circle) = profiles(bore, quadrant);
    let mut net = [[[[0_i128; 4]; 2]; 2]; 3]; // u, v, w
    for i in 0..3 {
        for (v, raw) in [(0, outer[i]), (1, circle[i])] {
            let [rx, ry, weight] = raw;
            let cxw = mul(i128::from(bore.center_mm[0]), weight)?;
            let cyw = mul(i128::from(bore.center_mm[1]), weight)?;
            let x = add(cxw, rx)?;
            let y = add(cyw, ry)?;
            for (w, z) in [0_i128, i128::from(height)].into_iter().enumerate() {
                net[i][v][w] = [
                    add(x, mul(mul(i128::from(shift[0]), weight)?, w as i128)?)?,
                    add(y, mul(mul(i128::from(shift[1]), weight)?, w as i128)?)?,
                    mul(z, weight)?,
                    weight,
                ];
            }
        }
    }
    let mut controls = Vec::with_capacity(12);
    for i in 0..3 {
        for (v, v_layer) in net[0].iter().enumerate() {
            for (w, _) in v_layer.iter().enumerate() {
                let a = net[0][v][w];
                let b = net[1][v][w];
                let c = net[2][v][w];
                // Exact restriction of a quadratic homogeneous polynomial to
                // [0,1/2] or [1/2,1], written at common gauge eight.
                let h = match (half, i) {
                    (0, 0) => a
                        .map(|x| mul(4, x))
                        .into_iter()
                        .collect::<Result<Vec<_>, _>>()?,
                    (0, 1) => std::array::from_fn::<_, 4, _>(|k| 2 * (a[k] + b[k])).to_vec(),
                    (0, 2) | (1, 0) => {
                        std::array::from_fn::<_, 4, _>(|k| a[k] + 2 * b[k] + c[k]).to_vec()
                    }
                    (1, 1) => std::array::from_fn::<_, 4, _>(|k| 2 * (b[k] + c[k])).to_vec(),
                    (1, 2) => c
                        .map(|x| mul(4, x))
                        .into_iter()
                        .collect::<Result<Vec<_>, _>>()?,
                    _ => return Err("circular_polygon_half_range"),
                };
                controls.push([
                    checked_i64(h[0])?,
                    checked_i64(h[1])?,
                    checked_i64(h[2])?,
                    checked_i64(h[3])?,
                ]);
            }
        }
    }
    Ok(BracketBlock {
        id: format!(
            "bore_{}_sector_{}_half_{}_z0_{}",
            bore.id, QUADRANTS[quadrant], half, height
        ),
        degrees: [2, 1, 1],
        controls,
    })
}

fn power_quadratic(values: [i128; 3]) -> Result<[i128; 3], &'static str> {
    Ok([
        values[0],
        mul(
            2,
            values[1]
                .checked_sub(values[0])
                .ok_or("circular_polygon_integer_overflow")?,
        )?,
        add(values[0], values[2])?
            .checked_sub(mul(2, values[1])?)
            .ok_or("circular_polygon_integer_overflow")?,
    ])
}

fn product(a: &[i128], b: &[i128]) -> Result<Vec<i128>, &'static str> {
    let mut out = vec![0_i128; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            out[i + j] = add(out[i + j], mul(x, y)?)?;
        }
    }
    Ok(out)
}

fn prove_arc(block: &BracketBlock, center: [i32; 2], radius: i32) -> Result<(), &'static str> {
    let mut x_controls = [0_i128; 3];
    let mut y_controls = [0_i128; 3];
    let mut w_controls = [0_i128; 3];
    for i in 0..3 {
        let p = block.controls[4 * i + 2]; // v=1, w=0
        if p[3] <= 0 {
            return Err("circular_polygon_nonpositive_weight");
        }
        x_controls[i] = i128::from(p[0]) - i128::from(center[0]) * i128::from(p[3]);
        y_controls[i] = i128::from(p[1]) - i128::from(center[1]) * i128::from(p[3]);
        w_controls[i] = i128::from(p[3]);
    }
    let x = power_quadratic(x_controls)?;
    let y = power_quadratic(y_controls)?;
    let w = power_quadratic(w_controls)?;
    let xx = product(&x, &x)?;
    let yy = product(&y, &y)?;
    let ww = product(&w, &w)?;
    let rr = mul(i128::from(radius), i128::from(radius))?;
    for k in 0..5 {
        if add(xx[k], yy[k])? != mul(rr, ww[k])? {
            return Err("circular_polygon_not_exact_circle");
        }
    }
    // The angular derivative numerator X*Y' - Y*X' is cubic in power
    // coefficients. Its degree-three Bernstein controls are checked exactly.
    let mut angle = [0_i128; 4];
    for i in 0..3 {
        for j in 1..3 {
            let k = i + j - 1;
            angle[k] = add(angle[k], mul(x[i], mul(j as i128, y[j])?)?)?;
            angle[k] = angle[k]
                .checked_sub(mul(y[i], mul(j as i128, x[j])?)?)
                .ok_or("circular_polygon_integer_overflow")?;
        }
    }
    let bernstein3_times3 = [
        mul(3, angle[0])?,
        add(mul(3, angle[0])?, angle[1])?,
        add(add(mul(3, angle[0])?, mul(2, angle[1])?)?, angle[2])?,
        mul(3, add(add(angle[0], angle[1])?, add(angle[2], angle[3])?)?)?,
    ];
    if bernstein3_times3.iter().any(|&v| v <= 0) {
        return Err("circular_polygon_nonmonotone_arc");
    }
    Ok(())
}

fn side(block: &BracketBlock, axis: usize, upper: bool) -> BlockSide {
    BlockSide {
        block: block.id.clone(),
        axis,
        upper,
    }
}

/// Independently audit the exact circular-bore template, G18 guard partition,
/// continuous local regularity, whole-face incidence, and resulting boundary.
pub fn audit_circular_polygon_prism(
    recipe: &CircularPolygonPrism,
    blocks: &[BracketBlock],
) -> Result<CircularPolygonCertificate, &'static str> {
    let guards = guard_profiles(recipe)?;
    let outer_recipe = PerforatedPrism {
        outer: recipe.outer.clone(),
        holes: guards,
        height_mm: recipe.height_mm,
        top_shift_mm: recipe.top_shift_mm,
    };
    let total_vertices = recipe.outer.vertices.len() + 4 * recipe.bores.len();
    let triangle_count = total_vertices + 2 * recipe.bores.len() - 2;
    let exterior_count = 3 * triangle_count;
    if blocks.len() != exterior_count + 8 * recipe.bores.len() {
        return Err("circular_polygon_block_count");
    }
    let mut ids = BTreeSet::new();
    for block in blocks {
        if !ids.insert(block.id.as_str()) {
            return Err("circular_polygon_duplicate_block");
        }
    }
    let outer_blocks = blocks
        .iter()
        .filter(|b| b.id.starts_with("tri/"))
        .cloned()
        .collect::<Vec<_>>();
    if outer_blocks.len() != exterior_count {
        return Err("circular_polygon_exterior_block_count");
    }
    let outer_certificate = audit_perforated_prism(&outer_recipe, &outer_blocks)?;
    let mut expected = BTreeMap::<String, BracketBlock>::new();
    let mut bore_by_block = BTreeMap::<String, String>::new();
    let mut minimum_determinant = outer_certificate.certificate.determinant_lower_bound_mm3;
    for bore in &recipe.bores {
        for quadrant in 0..4 {
            for half in 0..2 {
                let block = expected_ring_block(
                    bore,
                    recipe.height_mm,
                    recipe.top_shift_mm,
                    quadrant,
                    half,
                )?;
                bore_by_block.insert(block.id.clone(), bore.id.clone());
                if expected.insert(block.id.clone(), block).is_some() {
                    return Err("circular_polygon_duplicate_expected_block");
                }
            }
        }
    }
    let mut actual_ring = BTreeMap::<String, &BracketBlock>::new();
    for block in blocks.iter().filter(|b| b.id.starts_with("bore_")) {
        if actual_ring.insert(block.id.clone(), block).is_some() {
            return Err("circular_polygon_duplicate_bore_block");
        }
    }
    if actual_ring.len() != expected.len() {
        return Err("circular_polygon_bore_block_count");
    }
    for (id, target) in &expected {
        let block = actual_ring
            .get(id)
            .ok_or("circular_polygon_missing_bore_block")?;
        if *block != target {
            return Err("circular_polygon_altered_bore_block");
        }
        let bore_id = bore_by_block
            .get(id)
            .ok_or("circular_polygon_bore_identity")?;
        let bore = recipe
            .bores
            .iter()
            .find(|b| b.id == *bore_id)
            .ok_or("circular_polygon_bore_identity")?;
        prove_arc(block, bore.center_mm, bore.radius_mm)?;
        let mut profile = Vec::with_capacity(6);
        for i in 0..3 {
            for v in 0..2 {
                let p = block.controls[4 * i + 2 * v];
                profile.push([p[0], p[1], p[3]]);
            }
        }
        let planar = certify_planar_refined([2, 1], &profile).map_err(|error| match error {
            PlanarError::Invalid => "circular_polygon_nonpositive_jacobian",
            PlanarError::Indeterminate => "circular_polygon_indeterminate_jacobian",
        })?;
        let determinant = next_down(planar * f64::from(recipe.height_mm));
        if !determinant.is_finite() || determinant <= 0.0 {
            return Err("circular_polygon_indeterminate_jacobian");
        }
        minimum_determinant = minimum_determinant.min(determinant);
    }
    if actual_ring.keys().any(|id| !expected.contains_key(id)) {
        return Err("circular_polygon_unexpected_bore_block");
    }

    let incidence: Incidence = enumerate_faces(blocks)?;
    let guard_tags = recipe
        .bores
        .iter()
        .map(|b| format!("edge/guard_{}/", b.id))
        .collect::<Vec<_>>();
    let expected_guard_sides = outer_certificate
        .certificate
        .design_faces
        .iter()
        .filter(|(_, tag)| guard_tags.iter().any(|prefix| tag.starts_with(prefix)))
        .map(|(s, _)| s.clone())
        .collect::<BTreeSet<_>>();
    let mut ring_guard_interfaces = BTreeMap::<String, usize>::new();
    let mut ring_sector_interfaces = BTreeMap::<String, usize>::new();
    for interface in &incidence.interfaces {
        let a = interface.first.block.starts_with("bore_");
        let b = interface.second.block.starts_with("bore_");
        if a ^ b {
            let ring_side = if a {
                &interface.first
            } else {
                &interface.second
            };
            let exterior_side = if a {
                &interface.second
            } else {
                &interface.first
            };
            let id = bore_by_block
                .get(&ring_side.block)
                .ok_or("circular_polygon_bore_identity")?;
            if !expected_guard_sides.contains(exterior_side)
                || ring_side.axis != 1
                || ring_side.upper
            {
                return Err("circular_polygon_guard_interface");
            }
            *ring_guard_interfaces.entry(id.clone()).or_default() += 1;
        } else if a && b {
            let id_a = bore_by_block
                .get(&interface.first.block)
                .ok_or("circular_polygon_bore_identity")?;
            let id_b = bore_by_block
                .get(&interface.second.block)
                .ok_or("circular_polygon_bore_identity")?;
            if id_a != id_b {
                return Err("circular_polygon_cross_bore_interface");
            }
            *ring_sector_interfaces.entry(id_a.clone()).or_default() += 1;
        }
    }
    for bore in &recipe.bores {
        if ring_guard_interfaces.get(&bore.id) != Some(&8)
            || ring_sector_interfaces.get(&bore.id) != Some(&8)
        {
            return Err("circular_polygon_interface_inventory");
        }
    }

    let mut design_faces = outer_certificate
        .certificate
        .design_faces
        .iter()
        .filter(|(_, tag)| !guard_tags.iter().any(|prefix| tag.starts_with(prefix)))
        .cloned()
        .collect::<Vec<_>>();
    for block in actual_ring.values() {
        design_faces.push((
            side(block, 1, true),
            format!(
                "bore/{}",
                bore_by_block
                    .get(&block.id)
                    .ok_or("circular_polygon_bore_identity")?
            ),
        ));
        design_faces.push((side(block, 2, false), "start".to_owned()));
        design_faces.push((side(block, 2, true), "end".to_owned()));
    }
    design_faces.sort_by(|a, b| a.0.cmp(&b.0));
    let actual_exterior = incidence.exterior.iter().cloned().collect::<BTreeSet<_>>();
    let expected_exterior = design_faces
        .iter()
        .map(|(s, _)| s.clone())
        .collect::<BTreeSet<_>>();
    if actual_exterior != expected_exterior || expected_exterior.len() != design_faces.len() {
        return Err("circular_polygon_boundary_inventory");
    }

    let mut outer_area_twice = 0_i64;
    for i in 0..recipe.outer.vertices.len() {
        let a = recipe.outer.vertices[i].xy_mm;
        let b = recipe.outer.vertices[(i + 1) % recipe.outer.vertices.len()].xy_mm;
        outer_area_twice = outer_area_twice
            .checked_add(i64::from(a[0]) * i64::from(b[1]) - i64::from(a[1]) * i64::from(b[0]))
            .ok_or("circular_polygon_area_overflow")?;
    }
    outer_area_twice = outer_area_twice
        .checked_abs()
        .ok_or("circular_polygon_area_overflow")?;
    let mut sum_squared_radii = 0_i64;
    for bore in &recipe.bores {
        let r = i64::from(bore.radius_mm);
        sum_squared_radii = sum_squared_radii
            .checked_add(r.checked_mul(r).ok_or("circular_polygon_area_overflow")?)
            .ok_or("circular_polygon_area_overflow")?;
    }
    if !minimum_determinant.is_finite() || minimum_determinant <= 0.0 {
        return Err("circular_polygon_indeterminate_jacobian");
    }
    design_faces.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(CircularPolygonCertificate {
        certificate: BracketCertificate {
            incidence,
            determinant_lower_bound_mm3: minimum_determinant,
            radial_error_upper_mm: 0.0,
            design_faces,
        },
        outer_area_twice_mm2: outer_area_twice,
        sum_squared_radii_mm2: sum_squared_radii,
    })
}
