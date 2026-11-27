//! Native P322 mechanics on the ten-chart crossed-tool S2 union.
//!
//! The four wing-to-core traces are imposed by preimage-aware numerical
//! collocation, matching the published independent study. This is a
//! deliberately empirical coupling route, not an exact mortar certificate.

use crate::extracted_space::{
    assemble_sampled, AxisExtraction, ElementBasisSample, ExtractedElement,
};
use crate::json::{object, Json};
use crate::{constraint_system::Ratio, elasticity::Stiffness};
use std::collections::{BTreeMap, BTreeSet};

const PU: usize = 3;
const PV: usize = 2;
const PW: usize = 2;
const NLOCAL: usize = (PU + 1) * (PV + 1) * (PW + 1);
const E_MPA: f64 = 1000.0;
const NU: f64 = 0.3;
const TRACE_RCOND: f64 = 1e-12;
const MAX_MODEL_BYTES: usize = 256 * 1024;
const MAX_RECIPE_BYTES: usize = 64 * 1024;
const MAX_SCALAR_DOFS: usize = 1600;

pub(crate) type Point = [f64; 4];
type CellKey = (usize, usize, usize, usize);
type DenseTransform = Vec<Vec<f64>>;
type ReducedSpace = (usize, usize, DenseTransform);

#[derive(Clone)]
struct Cell {
    key: CellKey,
    element: ExtractedElement,
}

struct Built {
    cells: Vec<Cell>,
    ids: BTreeMap<CellKey, usize>,
    scalar_dofs: usize,
    occurrence_compact: Vec<usize>,
}

fn local_index(i: usize, j: usize, k: usize) -> usize {
    (i * (PV + 1) + j) * (PW + 1) + k
}

fn array<const N: usize>(value: &Json) -> Result<[f64; N], &'static str> {
    let items = value.items()?;
    if items.len() != N {
        return Err("unexpected array length in crossed mechanics input");
    }
    let mut out: [f64; N] = [0.; N];
    for (target, source) in out.iter_mut().zip(items) {
        *target = source.number()?;
        if !(*target).is_finite() {
            return Err("nonfinite crossed mechanics control");
        }
    }
    Ok(out)
}

pub(crate) fn source_net(
    model: &Json,
) -> Result<BTreeMap<(usize, usize), Vec<Point>>, &'static str> {
    if model.to_string().len() > MAX_MODEL_BYTES {
        return Err("crossed bending model exceeds the input size budget");
    }
    if model.field("schema")?.text()? != "crossed_tool_charts_v1"
        || model.field("operation")?.text()? != "union"
    {
        return Err("crossed bending requires a native S2 union model");
    }
    let source_recipe = model.field("source_recipe")?;
    if source_recipe.to_string().len() > MAX_RECIPE_BYTES
        || source_recipe.field("family")?.text()? != "inclined_crossed_tools_v1"
    {
        return Err("crossed bending source recipe exceeds its supported form");
    }
    let supported_recipe = crate::json::parse(include_str!(
        "../tests/jobs/crossed_partial_rise3.recipe.json"
    ))
    .map_err(|_| "invalid built-in S2 source recipe")?;
    if source_recipe.field("sketch_pair")? != supported_recipe.field("sketch_pair")?
        || ![-3, 0, 3].contains(&source_recipe.field("transverse_rise_mm")?.number::<i64>()?)
    {
        return Err("crossed bending accepts only the recorded S2 sketch and tested rises");
    }
    let replayed = crate::g2_cross_plate::boolean::report(source_recipe, "union")
        .map_err(|_| "crossed bending could not replay its native source recipe")?;
    if &replayed != model {
        return Err("crossed bending source replay differs from saved S2 geometry");
    }
    let blocks = model.field("blocks")?.items()?;
    if blocks.len() != 10 {
        return Err("crossed bending requires exactly ten native S2 charts");
    }
    let mut nets = BTreeMap::new();
    for block in blocks {
        let id = block.field("id")?.text()?;
        let parts: Vec<_> = id.split('_').collect();
        if parts.len() != 4 || !parts[0].starts_with("band") || !parts[1].starts_with("col") {
            return Err("unsupported crossed chart identifier");
        }
        let band = parts[0][4..].parse().map_err(|_| "invalid S2 band id")?;
        let col = parts[1][3..].parse().map_err(|_| "invalid S2 column id")?;
        if array::<3>(block.field("degrees")?)? != [3., 1., 1.] {
            return Err("crossed bending requires source degree (3,1,1)");
        }
        let controls = block.field("controls")?.items()?;
        if controls.len() != 16 {
            return Err("crossed S2 chart requires sixteen homogeneous controls");
        }
        let net: Vec<Point> = controls
            .iter()
            .map(|point| array::<4>(point))
            .collect::<Result<_, _>>()?;
        if net.iter().any(|p| p[3] <= 0.) {
            return Err("crossed S2 chart has a nonpositive weight");
        }
        if nets.insert((band, col), net).is_some() {
            return Err("duplicate crossed S2 chart key");
        }
    }
    if nets.len() != 10 {
        return Err("incomplete crossed S2 chart set");
    }
    Ok(nets)
}

fn source_point(net: &[Point], i: usize, v: f64, w: f64) -> Point {
    let mut p = [0.; 4];
    for a in 0..2 {
        for b in 0..2 {
            let weight = (if a == 0 { 1. - v } else { v }) * (if b == 0 { 1. - w } else { w });
            for c in 0..4 {
                p[c] += weight * net[(i * 2 + a) * 2 + b][c];
            }
        }
    }
    p
}

fn elevated_child(net: &[Point], iv: usize, iw: usize, nv: usize, nw: usize) -> Vec<Point> {
    let v0 = iv as f64 / nv as f64;
    let v1 = (iv + 1) as f64 / nv as f64;
    let w0 = iw as f64 / nw as f64;
    let w1 = (iw + 1) as f64 / nw as f64;
    let mut child = vec![[0.; 4]; NLOCAL];
    for i in 0..=PU {
        let corners = [
            source_point(net, i, v0, w0),
            source_point(net, i, v1, w0),
            source_point(net, i, v0, w1),
            source_point(net, i, v1, w1),
        ];
        for j in 0..=PV {
            let s = j as f64 / PV as f64;
            let left: Point = std::array::from_fn(|c| (1. - s) * corners[0][c] + s * corners[1][c]);
            let right: Point =
                std::array::from_fn(|c| (1. - s) * corners[2][c] + s * corners[3][c]);
            for k in 0..=PW {
                let t = k as f64 / PW as f64;
                child[local_index(i, j, k)] =
                    std::array::from_fn(|c| (1. - t) * left[c] + t * right[c]);
            }
        }
    }
    child
}

fn identity_axis(degree: usize) -> AxisExtraction {
    AxisExtraction {
        degree,
        basis_count: degree + 1,
        row_offsets: (0..=degree + 1).collect(),
        columns: (0..=degree).collect(),
        coefficients: (0..=degree).map(|_| Ratio::one()).collect(),
    }
}

fn build(model: &Json, nv: usize, nw: usize) -> Result<Built, &'static str> {
    if !(1..=4).contains(&nv) || !(1..=2).contains(&nw) {
        return Err("crossed bending refinements must satisfy 1<=v<=4 and 1<=w<=2");
    }
    let parent_nets = source_net(model)?;
    let ncell = parent_nets.len() * nv * nw;
    let mut cells = Vec::with_capacity(ncell);
    let mut ids = BTreeMap::new();
    for (&(band, col), net) in &parent_nets {
        for iv in 0..nv {
            for iw in 0..nw {
                let key = (band, col, iv, iw);
                let geometry_controls = elevated_child(net, iv, iw, nv, nw);
                let base = cells.len() * NLOCAL;
                let element = ExtractedElement {
                    axes: [identity_axis(PU), identity_axis(PV), identity_axis(PW)],
                    dofs: (base..base + NLOCAL).collect(),
                    field_weights: geometry_controls.iter().map(|p| p[3]).collect(),
                    geometry_controls,
                };
                element.validate()?;
                ids.insert(key, cells.len());
                cells.push(Cell { key, element });
            }
        }
    }

    // Identify exact whole-face and h-face coefficients. Element evaluation
    // keeps unique local occurrence IDs; assembly merges shared pole IDs.
    let occurrences = cells.len() * NLOCAL;
    let mut uf = crate::analysis_space::UnionFind::new(occurrences);
    let node =
        |element: usize, i: usize, j: usize, k: usize| element * NLOCAL + local_index(i, j, k);
    for cell in &cells {
        let (band, col, iv, iw) = cell.key;
        let e = ids[&cell.key];
        if iv + 1 < nv {
            let other = ids[&(band, col, iv + 1, iw)];
            for i in 0..=PU {
                for k in 0..=PW {
                    uf.join(node(e, i, PV, k), node(other, i, 0, k));
                }
            }
        }
        if iw + 1 < nw {
            let other = ids[&(band, col, iv, iw + 1)];
            for i in 0..=PU {
                for j in 0..=PV {
                    uf.join(node(e, i, j, PW), node(other, i, j, 0));
                }
            }
        }
        if let Some(&other) = ids.get(&(band + 1, col, iv, iw)) {
            for j in 0..=PV {
                for k in 0..=PW {
                    uf.join(node(e, PU, j, k), node(other, 0, j, k));
                }
            }
        }
    }
    for (band, axis_i) in [(0, 0), (5, PU)] {
        for global_w in 0..=nw * PW {
            let iw = (global_w / PW).min(nw - 1);
            let k = global_w - iw * PW;
            let Some(&first) = ids.get(&(band, 1, 0, iw)) else {
                continue;
            };
            let base = node(first, axis_i, 0, k);
            for iv in 0..nv {
                let e = ids[&(band, 1, iv, iw)];
                for j in 0..=PV {
                    uf.join(base, node(e, axis_i, j, k));
                }
            }
        }
    }
    let roots: Vec<_> = (0..occurrences).map(|i| uf.find(i)).collect();
    let unique: BTreeSet<_> = roots.iter().copied().collect();
    let compact: BTreeMap<_, _> = unique
        .into_iter()
        .enumerate()
        .map(|(index, root)| (root, index))
        .collect();
    let occurrence_compact = roots.iter().map(|root| compact[root]).collect();
    if compact.len() > MAX_SCALAR_DOFS {
        return Err("crossed bending model exceeds the dense trace-space budget");
    }
    Ok(Built {
        cells,
        ids,
        scalar_dofs: compact.len(),
        occurrence_compact,
    })
}

fn partial_trace_rows(
    built: &Built,
    nv: usize,
    nw: usize,
) -> Result<(Vec<Vec<f64>>, f64, usize), &'static str> {
    let mut rows = Vec::new();
    let mut max_gap: f64 = 0.;
    let mut interface_count = 0;
    for band in [2usize, 3] {
        for col in [0usize, 2] {
            interface_count += 1;
            let side = if col == 0 { 0.0 } else { 1.0 };
            let slave_iv = if col == 0 { 0 } else { nv - 1 };
            let master_iv = if col == 0 { nv - 1 } else { 0 };
            let mut us = Vec::new();
            for iw in 0..=nw {
                let w = iw as f64 / nw as f64;
                let samples: Vec<_> = (0..=64).map(|i| i as f64 / 64.).collect();
                let z: Vec<f64> = samples
                    .iter()
                    .map(|&u| chart_position(built, band, col, slave_iv, nv, nw, u, side, w))
                    .collect::<Result<Vec<_>, _>>()?
                    .into_iter()
                    .map(|point| point[2])
                    .collect();
                for level in 1..nw {
                    let target = 12. * level as f64 / nw as f64;
                    if z.iter()
                        .map(|value| (value - target).abs())
                        .fold(0.0, f64::max)
                        < 1e-12
                    {
                        continue;
                    }
                    for pair in 0..samples.len() - 1 {
                        let f0 = z[pair] - target;
                        let f1 = z[pair + 1] - target;
                        if f0.abs() < 1e-12 {
                            us.push(samples[pair]);
                        } else if f0 * f1 < 0. {
                            us.push(find_u_root(
                                built,
                                band,
                                col,
                                slave_iv,
                                nv,
                                nw,
                                side,
                                w,
                                target,
                                samples[pair],
                                samples[pair + 1],
                            )?);
                        }
                    }
                }
            }
            us.extend((0..=20).map(|i| i as f64 / 20.));
            for u in &mut us {
                *u = (*u * 1e14).round() / 1e14;
            }
            us.sort_by(f64::total_cmp);
            us.dedup_by(|a, b| *a == *b);

            for u in us {
                let z0 = chart_position(built, band, col, slave_iv, nv, nw, u, side, 0.)?[2];
                let z1 = chart_position(built, band, col, slave_iv, nv, nw, u, side, 1.)?[2];
                let mut cuts: Vec<_> = (0..=nw).map(|i| i as f64 / nw as f64).collect();
                if (z1 - z0).abs() > 1e-13 {
                    cuts.extend((1..nw).filter_map(|i| {
                        let cut = (12. * i as f64 / nw as f64 - z0) / (z1 - z0);
                        (-1e-12..=1. + 1e-12)
                            .contains(&cut)
                            .then_some(cut.clamp(0., 1.))
                    }));
                }
                cuts.sort_by(f64::total_cmp);
                cuts.dedup_by(|a, b| (*a - *b).abs() <= 1e-11);
                for pair in cuts.windows(2) {
                    if pair[1] - pair[0] < 1e-12 {
                        continue;
                    }
                    for sample in 0..=3 {
                        let global_w = pair[0] + (pair[1] - pair[0]) * sample as f64 / 3.;
                        let (gs, bs) =
                            chart_sample(built, band, col, slave_iv, nv, nw, u, side, global_w)?;
                        let global_mw = (gs.position[2] / 12.).clamp(0., 1.);
                        let miw = ((global_mw * nw as f64).floor() as usize).min(nw - 1);
                        let mw = global_mw * nw as f64 - miw as f64;
                        let master_side = if col == 0 { 1.0 } else { 0.0 };
                        let (gm, bm) = built.cells[built.ids[&(band, 1, master_iv, miw)]]
                            .element
                            .sample([u, master_side, mw])?;
                        max_gap = max_gap.max(
                            (0..3)
                                .map(|c| (gs.position[c] - gm.position[c]).powi(2))
                                .sum::<f64>()
                                .sqrt(),
                        );
                        let mut row = vec![0.; built.scalar_dofs];
                        for (a, &dof) in bs.dofs.iter().enumerate() {
                            row[built.occurrence_compact[dof]] += bs.values[a];
                        }
                        for (a, &dof) in bm.dofs.iter().enumerate() {
                            row[built.occurrence_compact[dof]] -= bm.values[a];
                        }
                        rows.push(row);
                    }
                }
            }
        }
    }
    Ok((rows, max_gap, interface_count))
}

#[allow(clippy::too_many_arguments)] // Global chart and parameter coordinates stay explicit at this boundary.
fn chart_sample(
    built: &Built,
    band: usize,
    col: usize,
    v_owner: usize,
    nv: usize,
    nw: usize,
    u: f64,
    global_v: f64,
    global_w: f64,
) -> Result<
    (
        crate::extracted_space::ElementGeometrySample,
        ElementBasisSample,
    ),
    &'static str,
> {
    let iv = ((global_v * nv as f64).floor() as usize).min(nv - 1);
    let iw = ((global_w * nw as f64).floor() as usize).min(nw - 1);
    let local_v = global_v * nv as f64 - iv as f64;
    let local_w = global_w * nw as f64 - iw as f64;
    let iv = if global_v == 1. { v_owner } else { iv };
    let element = &built.cells[built.ids[&(band, col, iv, iw)]].element;
    element.sample([u, local_v, local_w])
}

#[allow(clippy::too_many_arguments)] // This is the scalar-coordinate convenience wrapper for chart_sample.
fn chart_position(
    built: &Built,
    band: usize,
    col: usize,
    v_owner: usize,
    nv: usize,
    nw: usize,
    u: f64,
    global_v: f64,
    global_w: f64,
) -> Result<[f64; 3], &'static str> {
    Ok(
        chart_sample(built, band, col, v_owner, nv, nw, u, global_v, global_w)?
            .0
            .position,
    )
}

#[allow(clippy::too_many_arguments)]
fn find_u_root(
    built: &Built,
    band: usize,
    col: usize,
    v_owner: usize,
    nv: usize,
    nw: usize,
    global_v: f64,
    global_w: f64,
    target_z: f64,
    mut lo: f64,
    mut hi: f64,
) -> Result<f64, &'static str> {
    let mut f_lo =
        chart_position(built, band, col, v_owner, nv, nw, lo, global_v, global_w)?[2] - target_z;
    for _ in 0..64 {
        let mid = 0.5 * (lo + hi);
        let f_mid = chart_position(built, band, col, v_owner, nv, nw, mid, global_v, global_w)?[2]
            - target_z;
        if f_mid.abs() <= 1e-14 || hi - lo <= 2e-14 {
            return Ok(mid);
        }
        if f_lo * f_mid <= 0. {
            hi = mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

#[allow(clippy::needless_range_loop)] // Householder updates align vector offsets with matrix rows.
fn nullspace(
    rows: &[Vec<f64>],
    columns: usize,
    threshold: f64,
) -> Result<(usize, Vec<Vec<f64>>), &'static str> {
    if rows
        .iter()
        .any(|r| r.len() != columns || r.iter().any(|v| !v.is_finite()))
    {
        return Err("invalid crossed trace matrix");
    }
    let sigma_max = estimate_largest_singular_value(rows, columns)?;
    let cutoff = sigma_max * threshold;
    let nrows = rows.len();
    let limit = columns.min(nrows);

    // Pivoted Householder QR of A^T. The column pivots are the trace rows;
    // the retained reflectors give an orthonormal basis for range(A^T).
    let mut qr = vec![vec![0.; nrows]; columns];
    for (row_id, row) in rows.iter().enumerate() {
        for (column, &value) in row.iter().enumerate() {
            qr[column][row_id] = value;
        }
    }
    let mut reflectors = Vec::<(usize, Vec<f64>, f64)>::new();
    for k in 0..limit {
        let pivot = (k..nrows)
            .max_by(|&a, &b| tail_norm_squared(&qr, k, a).total_cmp(&tail_norm_squared(&qr, k, b)))
            .ok_or("empty crossed QR pivot")?;
        let pivot_norm = tail_norm_squared(&qr, k, pivot).sqrt();
        if !pivot_norm.is_finite() || pivot_norm <= cutoff {
            break;
        }
        for row in &mut qr {
            row.swap(k, pivot);
        }
        let x0 = qr[k][k];
        let alpha = if x0 >= 0. { -pivot_norm } else { pivot_norm };
        let mut vector: Vec<_> = (k..columns).map(|i| qr[i][k]).collect();
        vector[0] -= alpha;
        let norm_squared = vector.iter().map(|v| v * v).sum::<f64>();
        if !norm_squared.is_finite() || norm_squared <= f64::MIN_POSITIVE {
            return Err("crossed QR reflector lost numerical rank");
        }
        let beta = 2. / norm_squared;
        for j in k..nrows {
            let projection = vector
                .iter()
                .enumerate()
                .map(|(offset, v)| v * qr[k + offset][j])
                .sum::<f64>()
                * beta;
            for (offset, v) in vector.iter().enumerate() {
                qr[k + offset][j] -= projection * v;
            }
        }
        qr[k][k] = alpha;
        for row in qr.iter_mut().skip(k + 1) {
            row[k] = 0.;
        }
        reflectors.push((k, vector, beta));
    }
    let rank = reflectors.len();
    let free = columns - rank;
    let mut z = vec![vec![0.; free]; columns];
    for free_column in 0..free {
        let mut vector = vec![0.; columns];
        vector[rank + free_column] = 1.;
        for (k, reflector, beta) in reflectors.iter().rev() {
            let projection = reflector
                .iter()
                .enumerate()
                .map(|(offset, value)| value * vector[*k + offset])
                .sum::<f64>()
                * beta;
            for (offset, value) in reflector.iter().enumerate() {
                vector[*k + offset] -= projection * value;
            }
        }
        for (row, value) in vector.into_iter().enumerate() {
            z[row][free_column] = value;
        }
    }

    // Check every original equation against the computed complement. This
    // catches a rank cutoff that would otherwise silently release a trace mode.
    let worst_residual = rows
        .iter()
        .map(|row| {
            let norm = row.iter().map(|value| value * value).sum::<f64>().sqrt();
            (0..free)
                .map(|column| {
                    row.iter()
                        .enumerate()
                        .map(|(i, value)| value * z[i][column])
                        .sum::<f64>()
                        .abs()
                        / norm.max(f64::MIN_POSITIVE)
                })
                .fold(0.0_f64, f64::max)
        })
        .fold(0.0_f64, f64::max);
    if !worst_residual.is_finite() || worst_residual > threshold * 1.25 {
        return Err("crossed trace nullspace failed its original-equation residual check");
    }
    Ok((rank, z))
}

fn tail_norm_squared(matrix: &[Vec<f64>], first_row: usize, column: usize) -> f64 {
    matrix
        .iter()
        .skip(first_row)
        .map(|row| row[column] * row[column])
        .sum()
}

fn estimate_largest_singular_value(rows: &[Vec<f64>], columns: usize) -> Result<f64, &'static str> {
    let mut vector: Vec<_> = (0..columns)
        .map(|index| ((index + 1) as f64 * 0.731).sin())
        .collect();
    let initial_norm = vector.iter().map(|value| value * value).sum::<f64>().sqrt();
    if initial_norm <= f64::MIN_POSITIVE {
        return Err("empty crossed singular-value estimate");
    }
    for value in &mut vector {
        *value /= initial_norm;
    }
    for _ in 0..80 {
        let mut next = vec![0.; columns];
        for row in rows {
            let projection = row.iter().zip(&vector).map(|(a, b)| a * b).sum::<f64>();
            for (target, &coefficient) in next.iter_mut().zip(row) {
                *target += coefficient * projection;
            }
        }
        let norm = next.iter().map(|value| value * value).sum::<f64>().sqrt();
        if !norm.is_finite() || norm <= f64::MIN_POSITIVE {
            return Err("crossed trace matrix has no resolvable scale");
        }
        for (target, value) in vector.iter_mut().zip(next) {
            *target = value / norm;
        }
    }
    let sigma = rows
        .iter()
        .map(|row| {
            row.iter()
                .zip(&vector)
                .map(|(a, b)| a * b)
                .sum::<f64>()
                .powi(2)
        })
        .sum::<f64>()
        .sqrt();
    if !sigma.is_finite() || sigma <= f64::MIN_POSITIVE {
        return Err("crossed trace matrix has no nonzero singular value");
    }
    Ok(sigma)
}

fn reduced_rows(
    built: &Built,
    partial: &[Vec<f64>],
    nv: usize,
    rank_cutoff: f64,
) -> Result<ReducedSpace, &'static str> {
    let (trace_rank, ztrace) = nullspace(partial, built.scalar_dofs, rank_cutoff)?;
    let mut rows = partial.to_vec();
    let mut clamp = BTreeSet::new();
    for cell in &built.cells {
        let (_, col, iv, _) = cell.key;
        if col == 0 && iv + 1 == nv {
            for i in 0..=PU {
                for k in 0..=PW {
                    let occurrence = cell.element.dofs[local_index(i, PV, k)];
                    clamp.insert(built.occurrence_compact[occurrence]);
                }
            }
        }
    }
    for dof in &clamp {
        let mut row = vec![0.; built.scalar_dofs];
        row[*dof] = 1.;
        rows.push(row);
    }
    let (final_rank, z) = nullspace(&rows, built.scalar_dofs, rank_cutoff)?;
    if final_rank < trace_rank || ztrace.is_empty() {
        return Err("invalid crossed trace/clamp rank");
    }
    let clamped_rank = final_rank - trace_rank;
    Ok((trace_rank, clamped_rank, z))
}

fn multiply_stiffness(
    stiffness: &Stiffness,
    t: &[Vec<f64>],
    x: &[f64],
) -> Result<Vec<f64>, &'static str> {
    let n = t.len();
    let nf = t.first().ok_or("empty crossed transform")?.len();
    if x.len() != 3 * nf {
        return Err("wrong crossed reduced vector length");
    }
    let mut full = vec![0.; 3 * n];
    for i in 0..n {
        for k in 0..nf {
            for c in 0..3 {
                full[3 * i + c] += t[i][k] * x[3 * k + c];
            }
        }
    }
    let applied = stiffness.apply(&full)?;
    let mut reduced = vec![0.; 3 * nf];
    for i in 0..n {
        for k in 0..nf {
            for c in 0..3 {
                reduced[3 * k + c] += t[i][k] * applied[3 * i + c];
            }
        }
    }
    Ok(reduced)
}

fn project(t: &[Vec<f64>], full: &[f64]) -> Vec<f64> {
    let nf = t[0].len();
    let mut reduced = vec![0.; 3 * nf];
    for i in 0..t.len() {
        for k in 0..nf {
            for c in 0..3 {
                reduced[3 * k + c] += t[i][k] * full[3 * i + c];
            }
        }
    }
    reduced
}

fn expand(t: &[Vec<f64>], reduced: &[f64]) -> Vec<f64> {
    let nf = t[0].len();
    let mut full = vec![0.; 3 * t.len()];
    for i in 0..t.len() {
        for k in 0..nf {
            for c in 0..3 {
                full[3 * i + c] += t[i][k] * reduced[3 * k + c];
            }
        }
    }
    full
}

#[allow(clippy::needless_range_loop)] // Elimination updates a pivot row against the same columns.
fn solve_dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Result<Vec<f64>, &'static str> {
    let n = b.len();
    if a.len() != n || a.iter().any(|row| row.len() != n) {
        return Err("invalid dense crossed system");
    }
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
            .ok_or("empty crossed pivot")?;
        if a[pivot][k].abs() <= 1e-14 {
            return Err("singular crossed reduced stiffness");
        }
        a.swap(k, pivot);
        b.swap(k, pivot);
        let d = a[k][k];
        for j in k..n {
            a[k][j] /= d;
        }
        b[k] /= d;
        for i in 0..n {
            if i != k {
                let f = a[i][k];
                if f != 0. {
                    for j in k..n {
                        a[i][j] -= f * a[k][j];
                    }
                    b[i] -= f * b[k];
                }
            }
        }
    }
    Ok(b)
}

fn vector_dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn reduced_diagonal(
    stiffness: &Stiffness,
    transform: &[Vec<f64>],
) -> Result<Vec<f64>, &'static str> {
    let n = transform.len();
    let nf = transform[0].len();
    let mut diagonal = vec![0.; 3 * nf];
    for (scalar, transform_row) in transform.iter().enumerate().take(n) {
        let mut local = [0.; 3];
        for entry in stiffness.offsets[scalar]..stiffness.offsets[scalar + 1] {
            if stiffness.columns[entry] == scalar {
                for (component, value) in local.iter_mut().enumerate() {
                    *value = stiffness.values[entry][component][component];
                }
                break;
            }
        }
        for (master, &coefficient) in transform_row.iter().enumerate().take(nf) {
            let square = coefficient.powi(2);
            for component in 0..3 {
                diagonal[3 * master + component] += square * local[component];
            }
        }
    }
    if diagonal.iter().any(|x| !x.is_finite() || *x <= 0.) {
        return Err("crossed reduced Jacobi diagonal is not positive");
    }
    Ok(diagonal)
}

fn solve_projected_cg(
    stiffness: &Stiffness,
    transform: &[Vec<f64>],
    rhs: &[f64],
    tolerance: f64,
    max_iterations: usize,
) -> Result<(Vec<f64>, usize), &'static str> {
    let diagonal = reduced_diagonal(stiffness, transform)?;
    let mut solution = vec![0.; rhs.len()];
    let mut residual = rhs.to_vec();
    let mut preconditioned: Vec<_> = residual.iter().zip(&diagonal).map(|(r, d)| r / d).collect();
    let mut direction = preconditioned.clone();
    let mut rz = vector_dot(&residual, &preconditioned);
    let rhs_norm = vector_dot(rhs, rhs).sqrt().max(f64::MIN_POSITIVE);
    for iteration in 0..max_iterations {
        if vector_dot(&residual, &residual).sqrt() / rhs_norm <= tolerance {
            return Ok((solution, iteration));
        }
        let applied = multiply_stiffness(stiffness, transform, &direction)?;
        let curvature = vector_dot(&direction, &applied);
        if !curvature.is_finite() || curvature <= 0. {
            return Err("crossed reduced CG encountered nonpositive curvature");
        }
        let alpha = rz / curvature;
        for i in 0..solution.len() {
            solution[i] += alpha * direction[i];
            residual[i] -= alpha * applied[i];
        }
        if vector_dot(&residual, &residual).sqrt() / rhs_norm <= tolerance {
            return Ok((solution, iteration + 1));
        }
        for i in 0..preconditioned.len() {
            preconditioned[i] = residual[i] / diagonal[i];
        }
        let next_rz = vector_dot(&residual, &preconditioned);
        let beta = next_rz / rz;
        for i in 0..direction.len() {
            direction[i] = preconditioned[i] + beta * direction[i];
        }
        rz = next_rz;
    }
    Err("crossed reduced CG iteration budget exhausted")
}

fn reduced_symmetry_probe(
    stiffness: &Stiffness,
    transform: &[Vec<f64>],
) -> Result<f64, &'static str> {
    let n = 3 * transform[0].len();
    let mut worst: f64 = 0.;
    for p in 0..6 {
        let x: Vec<_> = (0..n)
            .map(|i| ((i + 1) as f64 * (p + 1) as f64 * 0.013).sin())
            .collect();
        let y: Vec<_> = (0..n)
            .map(|i| ((i + 3) as f64 * (p + 2) as f64 * 0.017).cos())
            .collect();
        let ax = multiply_stiffness(stiffness, transform, &x)?;
        let ay = multiply_stiffness(stiffness, transform, &y)?;
        let defect = (vector_dot(&x, &ay) - vector_dot(&y, &ax)).abs();
        let scale = (vector_dot(&x, &x) * vector_dot(&ay, &ay)).sqrt()
            + (vector_dot(&y, &y) * vector_dot(&ax, &ax)).sqrt();
        worst = worst.max(defect / scale.max(f64::MIN_POSITIVE));
    }
    Ok(worst)
}

fn traction_force(built: &Built, order: usize, force: &mut [f64]) -> Result<(), &'static str> {
    let q = crate::elasticity::gauss(order)?;
    for cell in &built.cells {
        let (band, col, iv, _) = cell.key;
        if col != 2 || iv != 0 {
            continue;
        }
        let _ = band;
        for &(u, wu) in &q {
            for &(w, ww) in &q {
                let (g, b) = cell.element.sample([u, 0., w])?;
                let ju = g.jacobian_columns[0];
                let jw = g.jacobian_columns[2];
                let area = [
                    -(jw[1] * ju[2] - jw[2] * ju[1]),
                    -(jw[2] * ju[0] - jw[0] * ju[2]),
                    -(jw[0] * ju[1] - jw[1] * ju[0]),
                ];
                let weight = wu * ww;
                for (a, &occurrence) in b.dofs.iter().enumerate() {
                    let dof = built.occurrence_compact[occurrence];
                    force[3 * dof] += weight * b.values[a] * area[0];
                }
            }
        }
    }
    Ok(())
}

fn number(v: impl ToString) -> Json {
    Json::Number(v.to_string())
}

/// Analyze one already-built native S2 union at a specified Q322 refinement.
/// The partial traces use the explicitly recorded collocation contract.
#[allow(clippy::needless_range_loop)] // Matrix symmetrization updates mirrored entries by index.
pub fn analyze(model: &Json, nv: usize, nw: usize, order: usize) -> Result<Json, &'static str> {
    if !(4..=14).contains(&order) {
        return Err("crossed bending quadrature order must be in 4..=14");
    }
    let built = build(model, nv, nw)?;
    let (partial_rows, geometry_gap, partial_interfaces) = partial_trace_rows(&built, nv, nw)?;
    let (trace_rank, clamp_rank, transform) = reduced_rows(&built, &partial_rows, nv, TRACE_RCOND)?;
    if transform.is_empty() {
        return Err("crossed bending has no free coefficients");
    }
    let lambda = E_MPA * NU / ((1. + NU) * (1. - 2. * NU));
    let mu = E_MPA / (2. * (1. + NU));
    let assembled = assemble_sampled(
        built.scalar_dofs,
        built.cells.len(),
        lambda,
        mu,
        order,
        |cell_index, uvw| {
            let (geometry, local_basis) = built.cells[cell_index].element.sample(uvw)?;
            let mut merged = BTreeMap::<usize, (f64, [f64; 3])>::new();
            for (local, &occurrence) in local_basis.dofs.iter().enumerate() {
                let dof = built.occurrence_compact[occurrence];
                let entry = merged.entry(dof).or_default();
                entry.0 += local_basis.values[local];
                for axis in 0..3 {
                    entry.1[axis] += local_basis.physical_gradients[local][axis];
                }
            }
            let mut dofs = Vec::with_capacity(merged.len());
            let mut values = Vec::with_capacity(merged.len());
            let mut physical_gradients = Vec::with_capacity(merged.len());
            for (dof, (value, gradient)) in merged {
                dofs.push(dof);
                values.push(value);
                physical_gradients.push(gradient);
            }
            let basis = ElementBasisSample {
                dofs,
                values,
                physical_gradients,
            };
            Ok((geometry, basis))
        },
        |_| [0.; 3],
    )?;
    let mut force = assembled.force.clone();
    traction_force(&built, order, &mut force)?;
    let nfree = transform[0].len();
    let nd = 3 * nfree;
    if nd > 3000 || built.cells.len() > 80 {
        return Err("crossed bending request exceeds the native analysis budget");
    }
    let rhs = project(&transform, &force);
    let (reduced, iterations, symmetry, symmetry_method) = if nd <= 700 {
        let mut matrix = vec![vec![0.; nd]; nd];
        for j in 0..nd {
            let mut unit = vec![0.; nd];
            unit[j] = 1.;
            let column = multiply_stiffness(&assembled.stiffness, &transform, &unit)?;
            for i in 0..nd {
                matrix[i][j] = column[i];
            }
        }
        let matrix_scale = matrix
            .iter()
            .flatten()
            .fold(0.0_f64, |largest, value| largest.max(value.abs()));
        let skew = (0..nd)
            .flat_map(|i| (0..i).map(move |j| (i, j)))
            .fold(0.0_f64, |largest, (i, j)| {
                largest.max((matrix[i][j] - matrix[j][i]).abs())
            });
        let symmetry = skew / matrix_scale.max(f64::MIN_POSITIVE);
        if !symmetry.is_finite() || symmetry > 1e-9 {
            return Err("crossed reduced stiffness failed its pre-symmetrization check");
        }
        // Only remove roundoff skew after recording and checking the raw matrix.
        for i in 0..nd {
            for j in 0..i {
                let average = 0.5 * (matrix[i][j] + matrix[j][i]);
                matrix[i][j] = average;
                matrix[j][i] = average;
            }
        }
        (
            solve_dense(matrix, rhs.clone())?,
            0,
            symmetry,
            "full_reduced_matrix",
        )
    } else {
        let symmetry = reduced_symmetry_probe(&assembled.stiffness, &transform)?;
        let maximum = (8 * nd).max(2000);
        let (solution, iterations) =
            solve_projected_cg(&assembled.stiffness, &transform, &rhs, 1e-10, maximum)?;
        (
            solution,
            iterations,
            symmetry,
            "six_matrix_free_bilinear_probes",
        )
    };
    let displacement = expand(&transform, &reduced);
    let internal = assembled.stiffness.apply(&displacement)?;
    let residual = project(
        &transform,
        &internal
            .iter()
            .zip(&force)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>(),
    );
    let residual_norm = residual.iter().map(|x| x * x).sum::<f64>().sqrt()
        / rhs
            .iter()
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt()
            .max(f64::MIN_POSITIVE);
    let applied_force = force.iter().step_by(3).sum::<f64>();
    let tip = force
        .iter()
        .zip(&displacement)
        .map(|(f, u)| f * u)
        .sum::<f64>()
        / applied_force;
    let energy = 0.5
        * displacement
            .iter()
            .zip(&internal)
            .map(|(u, f)| u * f)
            .sum::<f64>();
    // Check affine patch energy on the unreduced assembly using homogeneous
    // geometry weights as rational field weights.
    let mut affine_global = vec![0.; 3 * built.scalar_dofs];
    for cell in &built.cells {
        for (a, &occurrence) in cell.element.dofs.iter().enumerate() {
            let dof = built.occurrence_compact[occurrence];
            let p = cell.element.geometry_controls[a];
            affine_global[3 * dof] = p[0] / p[3] / E_MPA;
            affine_global[3 * dof + 1] = -NU * p[1] / p[3] / E_MPA;
            affine_global[3 * dof + 2] = -NU * p[2] / p[3] / E_MPA;
        }
    }
    let affine_internal = assembled.stiffness.apply(&affine_global)?;
    let affine_energy = affine_global
        .iter()
        .zip(&affine_internal)
        .map(|(u, f)| u * f)
        .sum::<f64>();
    let affine_error = (affine_energy / (assembled.stiffness.integrated_volume / E_MPA) - 1.).abs();
    let scalar_dofs_after_full = built.scalar_dofs;
    Ok(object([
        (
            "schema",
            Json::String("native_crossed_bending_result_v1".into()),
        ),
        (
            "degree",
            Json::Array(vec![number(PU), number(PV), number(PW)]),
        ),
        ("refine_v", number(nv)),
        ("refine_w", number(nw)),
        ("quadrature_order", number(order)),
        ("charts", number(10)),
        ("elements", number(built.cells.len())),
        (
            "scalar_dofs_after_full_identifications",
            number(scalar_dofs_after_full),
        ),
        ("partial_trace_constraint_rows", number(partial_rows.len())),
        ("partial_trace_rank", number(trace_rank)),
        (
            "scalar_dofs_after_partial_trace",
            number(scalar_dofs_after_full - trace_rank),
        ),
        ("clamped_scalar_dofs", number(clamp_rank)),
        ("free_scalar_dofs", number(nfree)),
        ("free_vector_dofs", number(3 * nfree)),
        ("partial_interface_count", number(partial_interfaces)),
        (
            "max_partial_face_geometry_mismatch_mm",
            number(geometry_gap),
        ),
        ("affine_energy_error", number(affine_error)),
        ("reduced_equilibrium_residual", number(residual_norm)),
        ("volume_mm3", number(assembled.stiffness.integrated_volume)),
        ("force_N", number(applied_force)),
        ("loaded_end_mean_ux_mm", number(tip)),
        ("strain_energy_Nmm", number(energy)),
        ("reduced_stiffness_dimension", number(nd)),
        ("reduced_symmetry_relative", number(symmetry)),
        (
            "symmetry_check_method",
            Json::String(symmetry_method.into()),
        ),
        ("cg_iterations", number(iterations)),
    ]))
}
