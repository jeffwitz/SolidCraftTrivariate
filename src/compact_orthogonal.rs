//! Compact ten-block candidate for the orthogonal filleted bracket.
//!
//! The candidate is deliberately restricted to a six-vertex L section with one
//! equal-clearance re-entrant corner, a longitudinal round on that corner, and
//! one through-bore in each straight wing. Each wing and the elbow is a complete
//! untrimmed partition of its portion of the design domain.

use crate::{
    bracket::BracketBlock,
    orthogonal_sketch::{self, OrthogonalSketch, SweepPartition},
    sketch_layers::{area2, edges},
};
use std::collections::{BTreeMap, BTreeSet};

type Bounds = [[i32; 2]; 3];
type CutClosureResult = (
    Vec<BracketBlock>,
    BTreeMap<(String, usize, bool), String>,
    Vec<[usize; 3]>,
    Vec<orthogonal_sketch::CutSource>,
);

fn ordered(a: i32, b: i32) -> [i32; 2] {
    [a.min(b), a.max(b)]
}

fn standard_l(j: &OrthogonalSketch, c: &orthogonal_sketch::Corner<'_>) -> Option<[i32; 2]> {
    if j.outer.vertices.len() != 6 {
        return None;
    }
    let f = c.f.vertex;
    let x_ext = if c.sign[0] > 0 {
        j.outer.vertices.iter().map(|p| p[0]).max()
    } else {
        j.outer.vertices.iter().map(|p| p[0]).min()
    };
    let y_ext = if c.sign[1] > 0 {
        j.outer.vertices.iter().map(|p| p[1]).max()
    } else {
        j.outer.vertices.iter().map(|p| p[1]).min()
    };
    let (Some(x_ext), Some(y_ext)) = (x_ext, y_ext) else {
        return None;
    };
    let back = c.origin;
    let expected = [
        back,
        [x_ext, back[1]],
        [x_ext, f[1]],
        f,
        [f[0], y_ext],
        [back[0], y_ext],
    ];
    let actual: BTreeSet<_> = j.outer.vertices.iter().copied().collect();
    (actual.len() == 6 && expected.into_iter().collect::<BTreeSet<_>>() == actual)
        .then_some([x_ext, y_ext])
}

fn make_arm_bounds(
    j: &OrthogonalSketch,
    c: &orthogonal_sketch::Corner<'_>,
    ends: [i32; 2],
) -> Result<[Bounds; 2], &'static str> {
    let f = c.f.vertex;
    let [x_far, y_far] = std::array::from_fn(|k| c.origin[k] + c.sign[k] * c.size_xy[k]);
    let [x_end, y_end] = ends;
    if (x_end - x_far) * c.sign[0] <= 0 || (y_end - y_far) * c.sign[1] <= 0 {
        return Err("compact_elbow_nonpositive_wing_length");
    }
    let [x_band, y_band] = c.origin;
    let east = [ordered(x_far, x_end), ordered(y_band, f[1]), [0, j.height]];
    let north = [ordered(x_band, f[0]), ordered(y_far, y_end), [0, j.height]];
    if east[0][0] == east[0][1]
        || east[1][0] == east[1][1]
        || north[0][0] == north[0][1]
        || north[1][0] == north[1][1]
    {
        return Err("compact_elbow_empty_wing");
    }
    Ok([east, north])
}

fn supported_bore(b: &orthogonal_sketch::Bore, arm: &Bounds, expected_axis: usize) -> bool {
    if b.axis != expected_axis || !(1..=1024).contains(&b.radius) {
        return false;
    }
    if !(0..3).all(|k| arm[k][0] < b.center[k] && b.center[k] < arm[k][1]) {
        return false;
    }
    // This template represents a complete transverse through-hole. Its axis
    // spans the entire wing interval, and its circular section stays strictly
    // inside the other two sides of the rectangular arm.
    (0..3)
        .filter(|&k| k != b.axis)
        .all(|k| b.center[k] - b.radius > arm[k][0] && b.center[k] + b.radius < arm[k][1])
}

fn elbow_faces(
    j: &OrthogonalSketch,
    c: &orthogonal_sketch::Corner<'_>,
    physical: &mut BTreeMap<(String, usize, bool), String>,
    blocks: &mut Vec<BracketBlock>,
) -> Result<(), &'static str> {
    let z = [0, j.height];
    let id = "elbow";
    for mut block in
        crate::sketch_fillet_chart::fillet_charts_rect([0, 0], c.size_xy, c.f.radius, z, id)?
    {
        orthogonal_sketch::transform_fillet(&mut block, c)?;
        let center = [
            c.f.vertex[0] + c.sign[0] * c.f.radius,
            c.f.vertex[1] + c.sign[1] * c.f.radius,
            0,
        ];
        orthogonal_sketch::certify_circular_side(&block, 2, center, c.f.radius)?;
        physical.insert(
            (block.id.clone(), 1, true),
            format!("fillet/{}/wall", c.f.id),
        );
        blocks.push(block);
    }
    Ok(())
}

/// Try the compact elbow candidate. `Ok(None)` means the geometry is outside
/// this template's declared class and the caller may use its established
/// partitioner. A construction or certification error discards this candidate.
pub(crate) fn try_build(j: &OrthogonalSketch) -> Result<Option<SweepPartition>, &'static str> {
    if j.fillets.len() != 1 || !(2..=8).contains(&j.bores.len()) || j.height <= 0 {
        return Ok(None);
    }
    if edges(&j.outer).any(|(a, b)| a[0] != b[0] && a[1] != b[1]) {
        return Ok(None);
    }
    let c = match orthogonal_sketch::corner(j, &j.fillets[0]) {
        Ok(c) => c,
        Err(_) => return Ok(None),
    };
    let Some(ends) = standard_l(j, &c) else {
        return Ok(None);
    };
    let arms = make_arm_bounds(j, &c, ends)?;
    if j.bores.len() > 2 && j.bores.iter().any(|b| b.axis != j.bores[0].axis) {
        return try_multi_arm_hosts(j, &c, arms);
    }
    if j.bores.iter().all(|b| b.axis == j.bores[0].axis) {
        return try_multi_arm_hosts(j, &c, arms);
    }
    let mut east_bore = None;
    let mut north_bore = None;
    for bore in &j.bores {
        match bore.axis {
            1 if east_bore.replace(bore).is_none() => {}
            0 if north_bore.replace(bore).is_none() => {}
            _ => return Ok(None),
        }
    }
    let (Some(east_bore), Some(north_bore)) = (east_bore, north_bore) else {
        return Ok(None);
    };
    if [east_bore, north_bore].iter().any(|b| {
        !(1..=1024).contains(&b.radius) || b.center.iter().any(|v| !(-4096..=4096).contains(v))
    }) {
        return Ok(None);
    }
    if !supported_bore(east_bore, &arms[0], 1) || !supported_bore(north_bore, &arms[1], 0) {
        return Ok(None);
    }
    let mut blocks = Vec::with_capacity(10);
    let mut physical = BTreeMap::new();
    elbow_faces(j, &c, &mut physical, &mut blocks)?;
    for (b, bounds, wing) in [(east_bore, arms[0], "east"), (north_bore, arms[1], "north")] {
        let charts = crate::compact_wing::build(
            b.axis,
            bounds,
            b.center,
            b.radius,
            &format!("wing_{wing}"),
        )?;
        if charts.len() != 4 {
            return Err("compact_wing_chart_count");
        }
        for block in charts {
            orthogonal_sketch::certify_circular_side(&block, b.axis, b.center, b.radius)?;
            // Compact wing charts expose their bore wall on the radial v=1
            // side. The circle identity is checked by the wing builder.
            physical.insert((block.id.clone(), 1, true), format!("bore/{}/wall", b.id));
            blocks.push(block);
        }
    }
    if blocks.len() != 10 {
        return Err("compact_orthogonal_block_count");
    }
    let mut planes: [BTreeSet<i32>; 3] = std::array::from_fn(|_| BTreeSet::new());
    for bounds in [c.bounds, arms[0], arms[1]] {
        for axis in 0..3 {
            planes[axis].extend(bounds[axis]);
        }
    }
    let planes = planes.map(|axis| axis.into_iter().collect());
    let guards = vec![
        (c.f.id.clone(), c.bounds),
        (east_bore.id.clone(), arms[0]),
        (north_bore.id.clone(), arms[1]),
    ];
    let h = i128::from(j.height);
    let r_f = i128::from(c.f.radius);
    let mut twice_constant = i128::from(area2(&j.outer)) * h;
    let mut four_pi = -r_f * r_f * h;
    twice_constant = twice_constant
        .checked_add(2 * r_f * r_f * h)
        .ok_or("compact_orthogonal_volume_overflow")?;
    for b in [east_bore, north_bore] {
        let r = i128::from(b.radius);
        let wing = if b.axis == east_bore.axis {
            arms[0]
        } else {
            arms[1]
        };
        let length = i128::from(wing[b.axis][1] - wing[b.axis][0]);
        four_pi = four_pi
            .checked_sub(4 * r * r * length)
            .ok_or("compact_orthogonal_volume_overflow")?;
    }
    // Retain the exact area identity as a precondition check: the corner
    // decomposition adds r²(1-π/4), and each through-bore removes πr²L.
    if twice_constant < 0 || four_pi >= 0 {
        return Err("compact_orthogonal_volume_sign");
    }
    let partition = orthogonal_sketch::finish_partition(
        j,
        blocks,
        physical,
        twice_constant,
        four_pi,
        planes,
        guards,
    )?;
    // The per-block Bernstein bounds returned by finish_partition certify
    // regularity; exact shared traces and the prescribed physical labels
    // certify conformity and provenance for this assembled candidate.
    Ok(Some(partition))
}

/// Candidate using the mapped-cut closure policy. The caller compares it with
/// the frozen compact plan before selection; this entry point keeps failures
/// separate from the historical builder.
pub(crate) fn try_propagated(j: &OrthogonalSketch) -> Result<Option<SweepPartition>, &'static str> {
    if j.fillets.len() != 1
        || j.bores.len() < 2
        || j.height <= 0
        || j.bores
            .iter()
            .any(|b| b.axis > 2 || !(1..=1024).contains(&b.radius))
        || edges(&j.outer).any(|(a, b)| a[0] != b[0] && a[1] != b[1])
    {
        return Ok(None);
    }
    let c = match orthogonal_sketch::corner(j, &j.fillets[0]) {
        Ok(c) => c,
        Err(_) => return Ok(None),
    };
    let Some(ends) = standard_l(j, &c) else {
        return Ok(None);
    };
    let arms = make_arm_bounds(j, &c, ends)?;
    let mut by_arm: [Vec<&orthogonal_sketch::Bore>; 2] = std::array::from_fn(|_| Vec::new());
    for bore in &j.bores {
        let arm_index = match bore.axis {
            1 => Some(0),
            0 => Some(1),
            2 => {
                let fits = (0..2)
                    .filter(|&i| {
                        bore.center[2] > arms[i][2][0]
                            && bore.center[2] < arms[i][2][1]
                            && bore.center[0] - bore.radius > arms[i][0][0]
                            && bore.center[0] + bore.radius < arms[i][0][1]
                            && bore.center[1] - bore.radius > arms[i][1][0]
                            && bore.center[1] + bore.radius < arms[i][1][1]
                    })
                    .collect::<Vec<_>>();
                if fits.len() != 1 {
                    return Ok(None);
                }
                fits.first().copied()
            }
            _ => None,
        };
        let Some(index) = arm_index else {
            return Ok(None);
        };
        by_arm[index].push(bore);
    }
    // A host is one tensor extrusion direction. Mixed bore axes in the same
    // arm require a non-isoparametric T-junction and are deliberately refused.
    let mut blocks = Vec::new();
    let mut physical = BTreeMap::new();
    elbow_faces(j, &c, &mut physical, &mut blocks)?;
    let mut guards = vec![(c.f.id.clone(), c.bounds)];
    for arm_index in 0..2 {
        let arm_axis = if arm_index == 0 { 1 } else { 0 };
        let group = &by_arm[arm_index];
        if group.is_empty() {
            blocks.push(crate::compact_hosts::plain_box(
                format!("propagated_wing_free_{arm_index}"),
                arms[arm_index],
            ));
            continue;
        }
        let bore_axes = group.iter().map(|b| b.axis).collect::<BTreeSet<_>>();
        let hosted = if bore_axes.len() > 1 {
            crate::mixed_axis_hosts::build(
                arms[arm_index],
                group,
                &format!("propagated_wing_{arm_index}"),
            )?
        } else {
            let bore_axis = group[0].axis;
            crate::compact_hosts::holed_box_propagated(
                bore_axis,
                arms[arm_index],
                group,
                &format!("propagated_wing_{arm_index}"),
                bore_axis == arm_axis,
            )?
        };
        let Some(crate::compact_hosts::HoledBox {
            blocks: mut hosted,
            physical: mut labels,
            guards: mut host_guards,
        }) = hosted
        else {
            return Ok(None);
        };
        blocks.append(&mut hosted);
        physical.append(&mut labels);
        guards.append(&mut host_guards);
    }

    let mut planes: [BTreeSet<i32>; 3] = std::array::from_fn(|_| BTreeSet::new());
    for bounds in [c.bounds, arms[0], arms[1]] {
        for axis in 0..3 {
            planes[axis].extend(bounds[axis]);
        }
    }
    for (_, bounds) in guards.iter().skip(1) {
        for axis in 0..3 {
            planes[axis].extend(bounds[axis]);
        }
    }
    let h = i128::from(j.height);
    let rf = i128::from(c.f.radius);
    let mut twice_constant = i128::from(area2(&j.outer))
        .checked_mul(h)
        .ok_or("mapped_cut_volume_overflow")?;
    twice_constant = twice_constant
        .checked_add(2 * rf * rf * h)
        .ok_or("mapped_cut_volume_overflow")?;
    let mut four_pi = -rf * rf * h;
    for bore in &j.bores {
        let length = if bore.axis == 2 {
            j.height
        } else {
            arms[if bore.axis == 1 { 0 } else { 1 }][bore.axis][1]
                - arms[if bore.axis == 1 { 0 } else { 1 }][bore.axis][0]
        };
        let r = i128::from(bore.radius);
        four_pi = four_pi
            .checked_sub(4 * r * r * i128::from(length))
            .ok_or("mapped_cut_volume_overflow")?;
    }
    let (blocks, physical, cut_rounds, cut_sources) = close_mapped_cuts(blocks, physical, 1024)?;
    let mut result = orthogonal_sketch::finish_partition(
        j,
        blocks,
        physical,
        twice_constant,
        four_pi,
        planes.map(|p| p.into_iter().collect()),
        guards,
    )?;
    result.cut_rounds = cut_rounds;
    result.cut_sources = cut_sources;
    Ok(Some(result))
}

fn close_mapped_cuts(
    mut blocks: Vec<BracketBlock>,
    mut physical: BTreeMap<(String, usize, bool), String>,
    budget: usize,
) -> Result<CutClosureResult, &'static str> {
    let mut rounds = Vec::new();
    let mut sources = Vec::new();
    for round in 0..8 {
        let (faces, requests) = crate::mapped_seams::discover(&blocks)?;
        if requests.is_empty() {
            return Ok((blocks, physical, rounds, sources));
        }
        let previous_count = blocks.len();
        let plan = crate::mapped_cut_propagation::propagate_exact_cuts(
            &blocks, &faces, &requests, budget,
        )?;
        if plan.child_count <= previous_count || plan.child_count > budget {
            return Err("mapped_cut_no_progress_or_budget");
        }
        rounds.push([
            plan.requested_cut_count,
            plan.propagated_cut_count,
            plan.child_count,
        ]);
        for cut in &plan.cuts {
            sources.push(orthogonal_sketch::CutSource {
                round,
                block: blocks[cut.block].id.clone(),
                axis: cut.axis,
                numerator: cut.parameter.numerator,
                denominator: cut.parameter.denominator,
                sources: cut.provenance.clone(),
            });
        }
        let mut transferred = BTreeMap::new();
        for ((id, axis, upper), label) in &physical {
            let source_index = blocks
                .iter()
                .position(|b| b.id == *id)
                .ok_or("mapped_cut_lost_physical_source")?;
            let mut count = 0;
            for child in plan
                .boundary_children
                .iter()
                .filter(|c| c.source_block == source_index && c.axis == *axis && c.upper == *upper)
            {
                let child_id = plan.blocks[child.child_block].id.clone();
                transferred.insert((child_id, *axis, *upper), label.clone());
                count += 1;
            }
            if count == 0 {
                return Err("mapped_cut_lost_physical_boundary");
            }
        }
        blocks = plan.blocks;
        physical = transferred;
    }
    Err("mapped_cut_round_budget")
}

fn try_multi_arm_hosts(
    j: &OrthogonalSketch,
    c: &orthogonal_sketch::Corner<'_>,
    arms: [Bounds; 2],
) -> Result<Option<SweepPartition>, &'static str> {
    if j.bores.iter().any(|b| {
        b.axis > 1
            || !(1..=1024).contains(&b.radius)
            || b.center.iter().any(|v| !(-4096..=4096).contains(v))
    }) {
        return Ok(None);
    }
    let mut blocks = Vec::new();
    let mut physical = BTreeMap::new();
    elbow_faces(j, c, &mut physical, &mut blocks)?;
    let mut guards = vec![(c.f.id.clone(), c.bounds)];
    for (arm_index, arm_bounds) in arms.iter().copied().enumerate() {
        let axis = if arm_index == 0 { 1 } else { 0 };
        let bores = j
            .bores
            .iter()
            .filter(|b| b.axis == axis)
            .collect::<Vec<_>>();
        if bores.is_empty() {
            blocks.push(crate::compact_hosts::plain_box(
                format!("wing_free_{arm_index}"),
                arm_bounds,
            ));
        } else {
            let Some(crate::compact_hosts::HoledBox {
                blocks: mut hosted,
                physical: mut labels,
                guards: mut host_guards,
            }) = crate::compact_hosts::holed_box(
                axis,
                arm_bounds,
                &bores,
                &format!("wing_{arm_index}"),
                true,
            )?
            else {
                return Ok(None);
            };
            blocks.append(&mut hosted);
            physical.append(&mut labels);
            guards.append(&mut host_guards);
        }
    }
    let mut planes: [BTreeSet<i32>; 3] = std::array::from_fn(|_| BTreeSet::new());
    for b in [c.bounds, arms[0], arms[1]] {
        for k in 0..3 {
            planes[k].extend(b[k]);
        }
    }
    for (_, host) in guards.iter().skip(1) {
        for k in 0..3 {
            planes[k].extend(host[k]);
        }
    }
    let mut twice_constant = i128::from(area2(&j.outer)) * i128::from(j.height);
    let r_f = i128::from(c.f.radius);
    let h = i128::from(j.height);
    twice_constant = twice_constant
        .checked_add(2 * r_f * r_f * h)
        .ok_or("compact_orthogonal_volume_overflow")?;
    let mut four_pi = -r_f * r_f * h;
    for bore in &j.bores {
        let arm_index = if bore.axis == 1 { 0 } else { 1 };
        let r = i128::from(bore.radius);
        four_pi = four_pi
            .checked_sub(
                4 * r
                    * r
                    * i128::from(arms[arm_index][bore.axis][1] - arms[arm_index][bore.axis][0]),
            )
            .ok_or("compact_orthogonal_volume_overflow")?;
    }
    if twice_constant < 0 || four_pi >= 0 {
        return Err("compact_orthogonal_volume_sign");
    }
    let partition = orthogonal_sketch::finish_partition(
        j,
        blocks,
        physical,
        twice_constant,
        four_pi,
        planes.map(|p| p.into_iter().collect()),
        guards,
    )?;
    Ok(Some(partition))
}
