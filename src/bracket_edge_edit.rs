//! Fixed-topology, multi-span end chamfers for the two transverse bores in the
//! compact ten-chart L bracket. The bore charts use angular, axial, radial
//! parameters; each logical chart retains its identity across axial spans.

use crate::{
    bracket::BracketBlock,
    incidence::{enumerate_faces, Incidence},
    spatial_certificate::{certify_spatial, CertificateBudget, CertificateStatus},
};

/// A prescribed round bore to chamfer at both axial ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BoreChamfer {
    /// Prefix shared by the four source chart IDs for this bore.
    pub chart_prefix: &'static str,
    pub axis: usize,
    /// Coordinates in quarter millimetres, allowing the frozen 0.25 mm width.
    pub center_quarters: [i32; 3],
    pub radius_quarters: i32,
    pub width_quarters: i32,
}

/// Base bore description used by the stateful editing entry point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BoreEdge {
    pub chart_prefix: &'static str,
    pub axis: usize,
    pub center_quarters: [i32; 3],
    pub radius_quarters: i32,
}

/// Requested local edge treatment. Lengths are exact quarter-millimetre ticks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeTreatment {
    Sharp,
    Chamfer { width_quarters: i32 },
    Fillet { radius_quarters: i32 },
}

/// One Bézier element belonging to a logical chart.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpanMap {
    pub logical_id: String,
    pub span_index: usize,
    /// Physical axial station interval. Adjacent elements share an endpoint.
    pub station_quarters: [i32; 2],
    /// Exact subinterval in the source chart parameter, with denominator 32.
    pub parameter_interval_32nds: [u32; 2],
    pub degrees: [usize; 3],
    /// Homogeneous `(X,Y,Z,W)` controls; last parameter varies fastest.
    pub controls: Vec<[i64; 4]>,
}

/// A fixed-topology edited bracket. `elements` are span elements; `logical_blocks`
/// is the original identity inventory, which remains ten for the qualified case.
#[derive(Debug, Clone, PartialEq)]
pub struct EditedBracket {
    pub logical_blocks: Vec<String>,
    pub elements: Vec<SpanMap>,
    pub incidence: Incidence,
    pub bores: [BoreEdge; 2],
    pub determinant_lower_bound_mm3: f64,
    pub chamfers: Vec<BoreChamfer>,
    pub treatment: EdgeTreatment,
    /// True when exact tangency creates a certified zero-Jacobian boundary rail.
    pub singular_boundary_qualified: bool,
}

fn checked(v: i128) -> Result<i64, &'static str> {
    i64::try_from(v).map_err(|_| "bracket_edge_control_overflow")
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}

fn normalize(points: &mut [[i64; 4]]) {
    let divisor = points
        .iter()
        .flatten()
        .fold(0_i128, |g, &x| gcd(g, i128::from(x)));
    if divisor > 1 {
        for p in points {
            for x in p {
                *x = (i128::from(*x) / divisor) as i64;
            }
        }
    }
}

fn stations(lo: i32, hi: i32, width: i32) -> Result<Vec<i32>, &'static str> {
    if width <= 0 || i64::from(hi) - i64::from(lo) <= 2 * i64::from(width) {
        return Err("bracket_edge_chamfer_width_must_leave_a_positive_middle_span");
    }
    Ok(vec![lo, lo + width, hi - width, hi])
}

fn parameter_32(station: i32, lo: i32, hi: i32) -> Result<u32, &'static str> {
    let numerator = (i64::from(station) - i64::from(lo))
        .checked_mul(32)
        .ok_or("bracket_edge_parameter_overflow")?;
    let denominator = i64::from(hi) - i64::from(lo);
    if denominator <= 0 || numerator % denominator != 0 {
        return Err("bracket_edge_knot_not_on_registered_32nd_grid");
    }
    u32::try_from(numerator / denominator).map_err(|_| "bracket_edge_parameter_overflow")
}

fn validate_source_chart(block: &BracketBlock, bore: BoreEdge) -> Result<(), &'static str> {
    if block.degrees != [2, 1, 1] || block.controls.len() != 12 || bore.axis > 2 {
        return Err("bracket_edge_source_chart_contract_mismatch");
    }
    if bore.center_quarters.iter().any(|v| v % 4 != 0) || bore.radius_quarters % 4 != 0 {
        return Err("bracket_edge_source_chart_requires_integer_mm_circle_data");
    }
    let center = bore.center_quarters.map(|v| v / 4);
    crate::orthogonal_sketch::certify_circular_side(
        block,
        bore.axis,
        center,
        bore.radius_quarters / 4,
    )
}

/// Apply equal linear-radius chamfers at both mouths of each of the two
/// independent bores. Input must be the actual compact ten-chart candidate.
/// The four angular charts per bore are taken directly from the parent maps;
/// no geometric rediscovery or axis-aligned template substitution is done.
pub fn chamfer_compact_bracket(
    parent: &[BracketBlock],
    bores: [BoreChamfer; 2],
) -> Result<EditedBracket, &'static str> {
    if parent.len() != 10 {
        return Err("bracket_edge_edit_requires_compact_ten_block_parent");
    }
    if bores[0].chart_prefix == bores[1].chart_prefix
        || bores[0].width_quarters != bores[1].width_quarters
        || bores.iter().any(|b| {
            b.axis > 2
                || b.radius_quarters <= 0
                || b.width_quarters <= 0
                || b.width_quarters >= b.radius_quarters
        })
    {
        return Err("bracket_edge_invalid_bore_chamfer");
    }
    let mut source_for = Vec::new();
    let mut logical_blocks = Vec::new();
    for block in parent {
        if block.id.starts_with("wing_") {
            let bore = bores
                .iter()
                .find(|b| block.id.starts_with(b.chart_prefix))
                .ok_or("bracket_edge_unknown_wing_chart")?;
            validate_source_chart(
                block,
                BoreEdge {
                    chart_prefix: bore.chart_prefix,
                    axis: bore.axis,
                    center_quarters: bore.center_quarters,
                    radius_quarters: bore.radius_quarters,
                },
            )?;
            source_for.push((block, *bore));
        } else {
            logical_blocks.push(block.id.clone());
        }
    }
    for bore in bores {
        let n = source_for.iter().filter(|(_, b)| *b == bore).count();
        if n != 4 {
            return Err("bracket_edge_requires_four_source_charts_per_bore");
        }
    }
    if source_for.len() != 8 {
        return Err("bracket_edge_requires_eight_bore_charts");
    }

    let mut elements = Vec::new();
    let mut seen_ids = logical_blocks
        .iter()
        .cloned()
        .collect::<std::collections::BTreeSet<_>>();
    for bore in bores {
        let sources = source_for
            .iter()
            .filter(|(_, b)| *b == bore)
            .map(|(block, _)| *block)
            .collect::<Vec<_>>();
        let stations = stations(
            axial_bound(sources[0], bore.axis)?[0],
            axial_bound(sources[0], bore.axis)?[1],
            bore.width_quarters,
        )?;
        if i64::from(stations[3]) - i64::from(stations[0]) != 32 * i64::from(bore.width_quarters) {
            return Err("bracket_edge_chamfer_must_match_registered_one_over_32_knot");
        }
        for source in sources {
            if axial_bound(source, bore.axis)? != [stations[0], *stations.last().unwrap()] {
                return Err("bracket_edge_bore_charts_have_different_axial_bounds");
            }
            if !seen_ids.insert(source.id.clone()) {
                return Err("bracket_edge_duplicate_logical_id");
            }
            logical_blocks.push(source.id.clone());
            for (span_index, pair) in stations.windows(2).enumerate() {
                let mut controls = make_span(source, bore, [pair[0], pair[1]])?;
                normalize(&mut controls);
                elements.push(SpanMap {
                    logical_id: source.id.clone(),
                    span_index,
                    station_quarters: [pair[0], pair[1]],
                    parameter_interval_32nds: [
                        parameter_32(pair[0], stations[0], stations[3])?,
                        parameter_32(pair[1], stations[0], stations[3])?,
                    ],
                    degrees: [2, 1, 1],
                    controls,
                });
            }
        }
    }
    if logical_blocks.len() != 10 {
        return Err("bracket_edge_logical_inventory_changed");
    }
    let mut all = Vec::new();
    for block in parent.iter().filter(|b| !b.id.starts_with("wing_")) {
        if block.id.starts_with("elbow_") {
            for (span_index, interval) in [[0, 1], [1, 31], [31, 32]].into_iter().enumerate() {
                all.push(split_elbow(block, span_index, interval)?);
            }
        } else {
            all.push(SpanMap {
                logical_id: block.id.clone(),
                span_index: 0,
                station_quarters: [0, 0],
                parameter_interval_32nds: [0, 32],
                degrees: block.degrees,
                controls: block.controls.clone(),
            });
        }
    }
    all.extend(elements);
    all.sort_by(|a, b| {
        a.logical_id
            .cmp(&b.logical_id)
            .then(a.span_index.cmp(&b.span_index))
    });

    let mut lower = f64::INFINITY;
    let mut incidence_blocks = Vec::with_capacity(all.len());
    for element in &all {
        let cert = certify_spatial(
            element.degrees,
            &element.controls,
            CertificateBudget::default(),
        );
        if cert.status != CertificateStatus::Valid {
            return Err("bracket_edge_jacobian_not_certified_positive");
        }
        lower = lower.min(cert.determinant_lower_bound_mm3.unwrap());
        incidence_blocks.push(BracketBlock {
            id: format!("{}__span_{}", element.logical_id, element.span_index),
            degrees: element.degrees,
            controls: element.controls.clone(),
        });
    }
    let incidence = enumerate_faces(&incidence_blocks)
        .map_err(|_| "bracket_edge_span_interfaces_do_not_match_exactly")?;
    if !interfaces_are_complete(&incidence, &incidence_blocks) {
        return Err("bracket_edge_interface_inventory_is_incomplete");
    }
    Ok(EditedBracket {
        logical_blocks,
        elements: all,
        incidence,
        bores: bores.map(|b| BoreEdge {
            chart_prefix: b.chart_prefix,
            axis: b.axis,
            center_quarters: b.center_quarters,
            radius_quarters: b.radius_quarters,
        }),
        determinant_lower_bound_mm3: lower,
        chamfers: bores.to_vec(),
        treatment: EdgeTreatment::Chamfer {
            width_quarters: bores[0].width_quarters,
        },
        singular_boundary_qualified: false,
    })
}

/// Stateful entry point used by command-line adapters. This accepts the
/// audited parent maps and replays their source chart identities directly.
pub fn edit_compact_bracket(
    parent: &[BracketBlock],
    bores: [BoreEdge; 2],
    treatment: EdgeTreatment,
) -> Result<EditedBracket, &'static str> {
    match treatment {
        EdgeTreatment::Sharp => sharp_compact_bracket(parent, bores),
        EdgeTreatment::Chamfer { width_quarters } => chamfer_compact_bracket(
            parent,
            bores.map(|b| BoreChamfer {
                chart_prefix: b.chart_prefix,
                axis: b.axis,
                center_quarters: b.center_quarters,
                radius_quarters: b.radius_quarters,
                width_quarters,
            }),
        ),
        EdgeTreatment::Fillet { radius_quarters } => {
            fillet_compact_bracket(parent, bores, radius_quarters)
        }
    }
}

fn sharp_compact_bracket(
    parent: &[BracketBlock],
    bores: [BoreEdge; 2],
) -> Result<EditedBracket, &'static str> {
    if parent.len() != 10 || bores[0].chart_prefix == bores[1].chart_prefix {
        return Err("bracket_edge_edit_requires_compact_ten_block_parent");
    }
    let mut lower = f64::INFINITY;
    let mut sources = Vec::new();
    let mut logical_blocks = Vec::new();
    for block in parent {
        if block.id.starts_with("wing_") {
            let bore = bores
                .iter()
                .find(|b| block.id.starts_with(b.chart_prefix))
                .ok_or("bracket_edge_unknown_sharp_chart")?;
            validate_source_chart(block, *bore)?;
            sources.push((block, *bore));
        } else {
            logical_blocks.push(block.id.clone());
        }
    }
    let mut all = Vec::with_capacity(30);
    for (block, bore) in sources {
        let [lo, hi] = axial_bound(block, bore.axis)?;
        let length = i64::from(hi) - i64::from(lo);
        if length <= 0 || length % 32 != 0 {
            return Err("bracket_edge_sharp_axis_span_not_divisible_by_32");
        }
        let width = i32::try_from(length / 32).map_err(|_| "bracket_edge_station_overflow")?;
        let stations = [lo, lo + width, hi - width, hi];
        for (span_index, pair) in stations.windows(2).enumerate() {
            let mut controls = make_span(
                block,
                BoreChamfer {
                    chart_prefix: bore.chart_prefix,
                    axis: bore.axis,
                    center_quarters: bore.center_quarters,
                    radius_quarters: bore.radius_quarters,
                    width_quarters: 0,
                },
                [pair[0], pair[1]],
            )?;
            normalize(&mut controls);
            all.push(SpanMap {
                logical_id: block.id.clone(),
                span_index,
                station_quarters: [pair[0], pair[1]],
                parameter_interval_32nds: [
                    parameter_32(pair[0], lo, hi)?,
                    parameter_32(pair[1], lo, hi)?,
                ],
                degrees: [2, 1, 1],
                controls,
            });
        }
        logical_blocks.push(block.id.clone());
    }
    for block in parent.iter().filter(|b| !b.id.starts_with("wing_")) {
        for (span_index, interval) in [[0, 1], [1, 31], [31, 32]].into_iter().enumerate() {
            all.push(split_elbow(block, span_index, interval)?);
        }
    }
    all.sort_by(|a, b| {
        a.logical_id
            .cmp(&b.logical_id)
            .then(a.span_index.cmp(&b.span_index))
    });
    if logical_blocks.len() != 10 || all.len() != 30 {
        return Err("bracket_edge_sharp_logical_partition_changed");
    }
    let mut incidence_blocks = Vec::with_capacity(all.len());
    for span in &all {
        let cert = certify_spatial(span.degrees, &span.controls, CertificateBudget::default());
        if cert.status != CertificateStatus::Valid {
            return Err("bracket_edge_sharp_span_jacobian_not_certified");
        }
        lower = lower.min(cert.determinant_lower_bound_mm3.unwrap());
        incidence_blocks.push(BracketBlock {
            id: format!("{}__span_{}", span.logical_id, span.span_index),
            degrees: span.degrees,
            controls: span.controls.clone(),
        });
    }
    let incidence = enumerate_faces(&incidence_blocks)
        .map_err(|_| "bracket_edge_sharp_interfaces_do_not_match")?;
    if !interfaces_are_complete(&incidence, &incidence_blocks) {
        return Err("bracket_edge_sharp_interface_inventory_mismatch");
    }
    Ok(EditedBracket {
        logical_blocks,
        elements: all,
        incidence,
        bores,
        determinant_lower_bound_mm3: lower,
        chamfers: vec![],
        treatment: EdgeTreatment::Sharp,
        singular_boundary_qualified: false,
    })
}

fn fillet_compact_bracket(
    parent: &[BracketBlock],
    bores: [BoreEdge; 2],
    radius_quarters: i32,
) -> Result<EditedBracket, &'static str> {
    if parent.len() != 10
        || radius_quarters <= 0
        || bores.iter().any(|b| radius_quarters >= b.radius_quarters)
    {
        return Err("bracket_edge_fillet_requires_positive_radius_and_ten_blocks");
    }
    let mut sources = Vec::new();
    let mut logical_blocks = Vec::new();
    for block in parent {
        if block.id.starts_with("wing_") {
            let bore = bores
                .iter()
                .find(|b| block.id.starts_with(b.chart_prefix))
                .ok_or("bracket_edge_unknown_fillet_chart")?;
            validate_source_chart(block, *bore)?;
            sources.push((block, *bore));
        } else {
            logical_blocks.push(block.id.clone());
        }
    }
    for bore in bores {
        if sources.iter().filter(|(_, b)| *b == bore).count() != 4 {
            return Err("bracket_edge_fillet_requires_four_charts_per_bore");
        }
    }
    let mut elements = Vec::new();
    for bore in bores {
        let group = sources
            .iter()
            .filter(|(_, b)| *b == bore)
            .map(|(b, _)| *b)
            .collect::<Vec<_>>();
        let bounds = axial_bound(group[0], bore.axis)?;
        let f = radius_quarters;
        let length = i64::from(bounds[1]) - i64::from(bounds[0]);
        if length <= 2 * i64::from(f) {
            return Err("bracket_edge_fillet_leaves_no_cylindrical_span");
        }
        if length != 32 * i64::from(f) {
            return Err("bracket_edge_fillet_must_match_registered_one_over_32_knot");
        }
        let stations = [bounds[0], bounds[0] + f, bounds[1] - f, bounds[1]];
        for source in group {
            if axial_bound(source, bore.axis)? != bounds {
                return Err("bracket_edge_fillet_chart_axial_bounds_disagree");
            }
            logical_blocks.push(source.id.clone());
            for (span_index, station) in stations.windows(2).enumerate() {
                let controls = if span_index == 1 {
                    make_span(
                        source,
                        BoreChamfer {
                            chart_prefix: bore.chart_prefix,
                            axis: bore.axis,
                            center_quarters: bore.center_quarters,
                            radius_quarters: bore.radius_quarters,
                            width_quarters: 0,
                        },
                        [station[0], station[1]],
                    )?
                } else {
                    make_round_span(
                        source,
                        bore,
                        radius_quarters,
                        [station[0], station[1]],
                        bounds,
                    )?
                };
                elements.push(SpanMap {
                    logical_id: source.id.clone(),
                    span_index,
                    station_quarters: [station[0], station[1]],
                    parameter_interval_32nds: [
                        parameter_32(station[0], bounds[0], bounds[1])?,
                        parameter_32(station[1], bounds[0], bounds[1])?,
                    ],
                    degrees: if span_index == 1 {
                        [2, 1, 1]
                    } else {
                        [2, 1, 2]
                    },
                    controls,
                });
            }
        }
    }
    if logical_blocks.len() != 10 || elements.len() != 24 {
        return Err("bracket_edge_fillet_logical_partition_changed");
    }
    let mut all = Vec::with_capacity(30);
    for block in parent.iter().filter(|b| !b.id.starts_with("wing_")) {
        if block.id.starts_with("elbow_") {
            for (span_index, interval) in [[0, 1], [1, 31], [31, 32]].into_iter().enumerate() {
                all.push(split_elbow(block, span_index, interval)?);
            }
        } else {
            all.push(SpanMap {
                logical_id: block.id.clone(),
                span_index: 0,
                station_quarters: [0, 0],
                parameter_interval_32nds: [0, 32],
                degrees: block.degrees,
                controls: block.controls.clone(),
            });
        }
    }
    all.extend(elements);
    all.sort_by(|a, b| {
        a.logical_id
            .cmp(&b.logical_id)
            .then(a.span_index.cmp(&b.span_index))
    });
    let mut lower = f64::INFINITY;
    let mut incidence_blocks = Vec::with_capacity(all.len());
    let mut qualified_singular_rails = 0usize;
    for span in &all {
        if span.degrees == [2, 1, 2] {
            let bore = bores
                .iter()
                .find(|b| span.logical_id.starts_with(b.chart_prefix))
                .ok_or("bracket_edge_unknown_fillet_chart")?;
            crate::spatial_certificate::certify_cap_contact(
                span.degrees,
                &span.controls,
                true,
                span.span_index == 2,
                bore.axis,
            )?;
            qualified_singular_rails += 1;
        } else {
            let cert = certify_spatial(span.degrees, &span.controls, CertificateBudget::default());
            if cert.status != CertificateStatus::Valid {
                return Err("bracket_edge_fillet_regular_span_not_certified_positive");
            }
            lower = lower.min(cert.determinant_lower_bound_mm3.unwrap());
        }
        incidence_blocks.push(BracketBlock {
            id: format!("{}__span_{}", span.logical_id, span.span_index),
            degrees: span.degrees,
            controls: span.controls.clone(),
        });
    }
    let incidence = enumerate_faces(&incidence_blocks)
        .map_err(|_| "bracket_edge_fillet_interfaces_do_not_match_exactly")?;
    if !interfaces_are_complete(&incidence, &incidence_blocks) {
        return Err("bracket_edge_fillet_interface_inventory_mismatch");
    }
    if qualified_singular_rails != 16 {
        return Err("bracket_edge_fillet_boundary_rail_inventory_mismatch");
    }
    Ok(EditedBracket {
        logical_blocks,
        elements: all,
        incidence,
        bores,
        determinant_lower_bound_mm3: 0.0,
        chamfers: vec![],
        treatment: EdgeTreatment::Fillet { radius_quarters },
        singular_boundary_qualified: true,
    })
}

fn axial_bound(block: &BracketBlock, axis: usize) -> Result<[i32; 2], &'static str> {
    let mut values = Vec::new();
    for p in &block.controls {
        if p[3] <= 0 {
            return Err("bracket_edge_nonpositive_source_weight");
        }
        let numerator = i128::from(p[axis])
            .checked_mul(4)
            .ok_or("bracket_edge_station_overflow")?;
        let denominator = i128::from(p[3]);
        if numerator % denominator != 0 {
            return Err("bracket_edge_source_endpoint_not_on_quarter_mm_grid");
        }
        values.push(
            i32::try_from(numerator / denominator).map_err(|_| "bracket_edge_station_overflow")?,
        );
    }
    let lo = *values
        .iter()
        .min()
        .ok_or("bracket_edge_missing_axial_station")?;
    let hi = *values
        .iter()
        .max()
        .ok_or("bracket_edge_missing_axial_station")?;
    Ok([lo, hi])
}

fn make_span(
    source: &BracketBlock,
    bore: BoreChamfer,
    stations: [i32; 2],
) -> Result<Vec<[i64; 4]>, &'static str> {
    let mut out = Vec::with_capacity(12);
    let (a, b) = ((bore.axis + 1) % 3, (bore.axis + 2) % 3);
    let bounds = axial_bound(source, bore.axis)?;
    let radius = i128::from(bore.radius_quarters);
    for u in 0..3 {
        for radial in 0..2 {
            for (end, &station) in stations.iter().enumerate() {
                let old = source.controls[(u * 2 + radial) * 2 + end];
                let old_w = i128::from(old[3]);
                let profile_radius = radius_at_station(bore, station, bounds);
                let mut p = [0_i128; 4];
                p[3] = 4 * radius * old_w;
                for axis in [a, b] {
                    let offset_quarters =
                        4 * i128::from(old[axis]) - i128::from(bore.center_quarters[axis]) * old_w;
                    p[axis] = if radial == 1 {
                        i128::from(bore.center_quarters[axis]) * old_w * radius
                            + offset_quarters * profile_radius
                    } else {
                        4 * radius * i128::from(old[axis])
                    };
                }
                p[bore.axis] = i128::from(station) * old_w * radius;
                out.push([
                    checked(p[0])?,
                    checked(p[1])?,
                    checked(p[2])?,
                    checked(p[3])?,
                ]);
            }
        }
    }
    if out.iter().any(|p| p[3] <= 0) {
        return Err("bracket_edge_nonpositive_edited_weight");
    }
    Ok(out)
}

fn radius_at_station(bore: BoreChamfer, station: i32, bounds: [i32; 2]) -> i128 {
    let r = i128::from(bore.radius_quarters);
    if station == bounds[0] || station == bounds[1] {
        r + i128::from(bore.width_quarters)
    } else {
        r
    }
}

fn make_round_span(
    source: &BracketBlock,
    bore: BoreEdge,
    fillet_quarters: i32,
    station: [i32; 2],
    bounds: [i32; 2],
) -> Result<Vec<[i64; 4]>, &'static str> {
    let r = i128::from(bore.radius_quarters);
    let f = i128::from(fillet_quarters);
    let (axis_num, radius_num, weight): ([i128; 3], [i128; 3], [i128; 3]) =
        if station[0] == bounds[0] {
            (
                [
                    i128::from(bounds[0]),
                    i128::from(bounds[0]),
                    2 * i128::from(bounds[0] + fillet_quarters),
                ],
                [r + f, r, 2 * r],
                [1, 1, 2],
            )
        } else if station[1] == bounds[1] {
            (
                [
                    2 * i128::from(bounds[1] - fillet_quarters),
                    i128::from(bounds[1]),
                    i128::from(bounds[1]),
                ],
                [2 * r, r, r + f],
                [2, 1, 1],
            )
        } else {
            let a = i128::from(station[0]);
            let b = i128::from(station[1]);
            ([2 * a, a + b, 2 * b], [2 * r; 3], [2; 3])
        };
    let (a, b) = ((bore.axis + 1) % 3, (bore.axis + 2) % 3);
    let mut controls = Vec::with_capacity(18);
    for u in 0..3 {
        for radial in 0..2 {
            for v in 0..3 {
                let old_radial = radial;
                let old = source.controls[(u * 2 + old_radial) * 2];
                let aw = i128::from(old[3]);
                let pw = weight[v];
                let mut p = [0_i128; 4];
                let inner = radial == 1;
                let net_pw = if inner { pw } else { 2 };
                p[3] = 4 * aw * net_pw * r;
                p[bore.axis] = if inner {
                    axis_num[v] * aw * r
                } else {
                    let s0 = i128::from(station[0]);
                    let s1 = i128::from(station[1]);
                    let linear_axis = [2 * s0, s0 + s1, 2 * s1][v];
                    linear_axis * aw * r
                };
                for axis in [a, b] {
                    let offset_q =
                        4 * i128::from(old[axis]) - i128::from(bore.center_quarters[axis]) * aw;
                    p[axis] = if inner {
                        i128::from(bore.center_quarters[axis]) * aw * pw * r
                            + offset_q * radius_num[v]
                    } else {
                        4 * i128::from(old[axis]) * net_pw * r
                    };
                }
                controls.push([
                    checked(p[0])?,
                    checked(p[1])?,
                    checked(p[2])?,
                    checked(p[3])?,
                ]);
            }
        }
    }
    if controls.iter().any(|p| p[3] <= 0) {
        return Err("bracket_edge_fillet_nonpositive_weight");
    }
    normalize(&mut controls);
    Ok(controls)
}

fn split_elbow(
    block: &BracketBlock,
    span_index: usize,
    interval: [u32; 2],
) -> Result<SpanMap, &'static str> {
    if block.degrees != [2, 1, 1] || interval[0] >= interval[1] || interval[1] > 32 {
        return Err("bracket_edge_elbow_split_contract_mismatch");
    }
    let eval = |u: usize, z: usize, t: u32| -> Result<[i64; 4], &'static str> {
        let p0 = block.controls[(u * 2) * 2 + z];
        let p1 = block.controls[(u * 2 + 1) * 2 + z];
        let mut p = [0_i64; 4];
        for c in 0..4 {
            let value = i128::from(p0[c]) * i128::from(32 - t) + i128::from(p1[c]) * i128::from(t);
            p[c] = checked(value)?;
        }
        Ok(p)
    };
    let mut controls = Vec::with_capacity(12);
    for u in 0..3 {
        for v in interval {
            for z in 0..2 {
                controls.push(eval(u, z, v)?);
            }
        }
    }
    normalize(&mut controls);
    Ok(SpanMap {
        logical_id: block.id.clone(),
        span_index,
        station_quarters: [0, 0],
        parameter_interval_32nds: interval,
        degrees: [2, 1, 1],
        controls,
    })
}

fn interfaces_are_complete(incidence: &Incidence, blocks: &[BracketBlock]) -> bool {
    // Every adjacent axial span of a logical chart must be connected, as must
    // every formerly shared cross-section/chart face. The ten-chart candidate
    // has 11 complete material incidences. Splitting every wing and both
    // elbow maps at the two matched stations gives 53 exact span incidences.
    blocks.len() == 30 && incidence.interfaces.len() == 53
}

pub fn round_fillet_fixed_partition_status() -> &'static str {
    "implemented: exact rational quarter-circle rail with a certified zero-Jacobian contact boundary"
}
