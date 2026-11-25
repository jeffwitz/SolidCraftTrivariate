//! Budgeted adaptive quadrature for surface tractions on certified block faces.
//!
//! The panel estimator compares a parent Gauss panel with the sum over its four
//! dyadic children. It checks area, the complete local nodal load vector,
//! resultant, and moment on each selected patch. This is an a-posteriori
//! convergence diagnostic; it is not a mathematical error bound.

use crate::{
    analysis_space::{AnalysisBlock, AnalysisSpace},
    boundary_traction::{LoadPatch, SurfaceLoad},
    elasticity::{basis_geometry, gauss, Vector},
    incidence::BlockSide,
};
use std::collections::{BTreeMap, BTreeSet};

/// Explicit work and convergence limits for surface quadrature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceQuadrature {
    pub relative_tolerance: f64,
    pub max_depth: u8,
    /// Counts every evaluated panel, including coarse panels used by the estimator.
    pub max_panels: usize,
}

impl SurfaceQuadrature {
    /// Check the limits independently of a mesh, including for kinematic-only jobs.
    pub fn validate(&self) -> Result<(), &'static str> {
        if !self.relative_tolerance.is_finite()
            || self.relative_tolerance <= 0.0
            || self.relative_tolerance > 1.0
            || self.max_depth == 0
            || self.max_depth > 30
            || self.max_panels == 0
        {
            return Err("invalid adaptive surface quadrature configuration");
        }
        Ok(())
    }
}

/// Adaptively integrated traction and convergence diagnostics.
pub struct AdaptiveSurfaceLoad {
    pub load: SurfaceLoad,
    /// Number of panel integrals actually evaluated, including estimator panels.
    pub panels: usize,
    /// Accepted dyadic depths relative to each requested patch in input order.
    pub patch_depths: Vec<[u8; 2]>,
    /// Maximum dyadic subdivision depth of accepted patches.
    pub deepest: u8,
    /// Largest final parent-vs-children estimate across selected patches/metrics.
    pub maximum_estimated_relative_error: f64,
}

type Rect = [[f64; 2]; 2];

#[derive(Clone, Default)]
struct PanelIntegral {
    nodal: [[f64; 3]; 27],
    area: f64,
    resultant: Vector,
    moment: Vector,
    absolute_force_scale: f64,
    absolute_moment_scale: f64,
}

fn cross(a: Vector, b: Vector) -> Vector {
    std::array::from_fn(|k| a[(k + 1) % 3] * b[(k + 2) % 3] - a[(k + 2) % 3] * b[(k + 1) % 3])
}

fn norm(vector: Vector) -> f64 {
    vector[0].hypot(vector[1]).hypot(vector[2])
}

fn patch_rect(patch: &LoadPatch) -> Result<Rect, &'static str> {
    if patch.depth > 30 {
        return Err("surface patch depth exceeds the dyadic coordinate range");
    }
    let scale = 1_u32 << patch.depth;
    if patch
        .intervals
        .iter()
        .any(|&[lo, hi]| lo >= hi || hi > scale)
    {
        return Err("empty or invalid surface patch rectangle");
    }
    Ok(patch
        .intervals
        .map(|span| span.map(|x| f64::from(x) / f64::from(scale))))
}

fn split(rect: Rect) -> [Rect; 4] {
    let mid = rect.map(|[lo, hi]| (lo + hi) * 0.5);
    [
        [[rect[0][0], mid[0]], [rect[1][0], mid[1]]],
        [[mid[0], rect[0][1]], [rect[1][0], mid[1]]],
        [[rect[0][0], mid[0]], [mid[1], rect[1][1]]],
        [[mid[0], rect[0][1]], [mid[1], rect[1][1]]],
    ]
}

fn integrate_panel<F: Fn(Vector, Vector) -> Vector>(
    block: &AnalysisBlock,
    face: &BlockSide,
    rect: Rect,
    quadrature: &[(f64, f64)],
    traction: &F,
) -> Result<PanelIntegral, &'static str> {
    let mut result = PanelIntegral {
        nodal: [[0.0; 3]; 27],
        area: 0.0,
        resultant: [0.0; 3],
        moment: [0.0; 3],
        absolute_force_scale: 0.0,
        absolute_moment_scale: 0.0,
    };
    let axes = match face.axis {
        0 => [1, 2],
        1 => [0, 2],
        2 => [0, 1],
        _ => return Err("invalid traction face axis"),
    };
    let sign = if face.upper ^ (face.axis == 1) {
        1.0
    } else {
        -1.0
    };
    let width = rect.map(|[lo, hi]| hi - lo);
    for &(u, wu) in quadrature {
        for &(v, wv) in quadrature {
            let mut parameters = [0.0; 3];
            parameters[face.axis] = f64::from(u8::from(face.upper));
            parameters[axes[0]] = rect[0][0] + width[0] * u;
            parameters[axes[1]] = rect[1][0] + width[1] * v;
            let sample = basis_geometry(block, parameters)?;
            let area_vector = cross(
                sample.jacobian_columns[axes[0]],
                sample.jacobian_columns[axes[1]],
            );
            let jacobian = norm(area_vector);
            if !jacobian.is_finite() || jacobian <= 0.0 {
                return Err("invalid surface Jacobian");
            }
            let measure = jacobian * wu * wv * width[0] * width[1];
            let normal = area_vector.map(|x| sign * x / jacobian);
            let density = traction(sample.position, normal);
            if density.iter().any(|x| !x.is_finite()) {
                return Err("nonfinite integrated traction");
            }
            let force = density.map(|x| x * measure);
            let moment = cross(sample.position, force);
            result.area += measure;
            result.absolute_force_scale += norm(density) * measure;
            result.absolute_moment_scale += norm(cross(sample.position, density)) * measure;
            for (c, (&force_component, &moment_component)) in force.iter().zip(&moment).enumerate()
            {
                result.resultant[c] += force_component;
                result.moment[c] += moment_component;
            }
            for n in 0..27 {
                for (c, &component) in force.iter().enumerate() {
                    result.nodal[n][c] += sample.values[n] * component;
                }
            }
        }
    }
    if !result.area.is_finite()
        || !result.absolute_force_scale.is_finite()
        || !result.absolute_moment_scale.is_finite()
        || result
            .resultant
            .iter()
            .chain(&result.moment)
            .any(|x| !x.is_finite())
        || result.nodal.iter().flatten().any(|x| !x.is_finite())
    {
        return Err("nonfinite adaptive surface panel integral");
    }
    Ok(result)
}

fn add(a: &mut PanelIntegral, b: &PanelIntegral) {
    for n in 0..27 {
        for c in 0..3 {
            a.nodal[n][c] += b.nodal[n][c];
        }
    }
    a.area += b.area;
    for c in 0..3 {
        a.resultant[c] += b.resultant[c];
        a.moment[c] += b.moment[c];
    }
    a.absolute_force_scale += b.absolute_force_scale;
    a.absolute_moment_scale += b.absolute_moment_scale;
}

fn sum_panels(panels: &[(Rect, PanelIntegral)]) -> PanelIntegral {
    let mut total = PanelIntegral {
        nodal: [[0.0; 3]; 27],
        area: 0.0,
        resultant: [0.0; 3],
        moment: [0.0; 3],
        absolute_force_scale: 0.0,
        absolute_moment_scale: 0.0,
    };
    for (_, panel) in panels {
        add(&mut total, panel);
    }
    total
}

fn relative_difference(difference: f64, magnitude: f64, physical_scale: f64) -> f64 {
    let denominator = if physical_scale == 0.0 || magnitude > physical_scale * 1e-12 {
        magnitude
    } else {
        physical_scale
    };
    if denominator == 0.0 {
        if difference == 0.0 {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        difference / denominator
    }
}

fn estimate_error(coarse: &PanelIntegral, fine: &PanelIntegral) -> Result<f64, &'static str> {
    let area_error = relative_difference(
        (fine.area - coarse.area).abs(),
        fine.area.abs(),
        fine.area.abs(),
    );
    let mut nodal_diff_squared = 0.0;
    let mut nodal_norm_squared = 0.0;
    for n in 0..27 {
        for c in 0..3 {
            nodal_diff_squared += (fine.nodal[n][c] - coarse.nodal[n][c]).powi(2);
            nodal_norm_squared += fine.nodal[n][c].powi(2);
        }
    }
    let nodal_error = relative_difference(
        nodal_diff_squared.sqrt(),
        nodal_norm_squared.sqrt(),
        fine.absolute_force_scale,
    );
    let force_error = relative_difference(
        norm(std::array::from_fn(|c| {
            fine.resultant[c] - coarse.resultant[c]
        })),
        norm(fine.resultant),
        fine.absolute_force_scale,
    );
    let moment_error = relative_difference(
        norm(std::array::from_fn(|c| fine.moment[c] - coarse.moment[c])),
        norm(fine.moment),
        fine.absolute_moment_scale,
    );
    let maximum = [area_error, nodal_error, force_error, moment_error]
        .into_iter()
        .fold(0.0_f64, f64::max);
    if maximum.is_finite() {
        Ok(maximum)
    } else {
        Err("adaptive surface estimate has zero local scale")
    }
}

/// Choose which face parameters can be refined independently.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceSplit {
    Isotropic,
    Directional,
}

fn composite_grid(
    rect: Rect,
    depths: [u8; 2],
    limit: usize,
    evaluated: &mut usize,
    sample: impl Fn(Rect) -> Result<PanelIntegral, &'static str>,
) -> Result<PanelIntegral, &'static str> {
    let [nu, nv] = depths.map(|depth| 1_usize << depth);
    let count = nu.checked_mul(nv).ok_or("surface panel count overflow")?;
    if count > limit.saturating_sub(*evaluated) {
        return Err("adaptive surface quadrature panel budget exceeded");
    }
    let span = |axis: usize, i: usize, n: usize| {
        let [lo, hi] = rect[axis];
        [i, i + 1].map(|j| lo + (hi - lo) * (j as f64 / n as f64))
    };
    let mut sum = PanelIntegral::default();
    for i in 0..nu {
        let u = span(0, i, nu);
        for j in 0..nv {
            add(&mut sum, &sample([u, span(1, j, nv)])?);
        }
    }
    *evaluated += count;
    Ok(sum)
}

/// Reuse the accepted grid from the immediately preceding primal integration.
/// Job loads select whole certified faces; patch ordering is shared with the primal.
pub(crate) fn replay_traction_rate(
    space: &AnalysisSpace,
    patches: &[LoadPatch],
    depths: &[[u8; 2]],
    order: usize,
    budget: usize,
    traction: impl Fn(Vector, Vector) -> Vector,
) -> Result<(Vec<f64>, usize), &'static str> {
    if patches.len() != depths.len() {
        return Err("frozen surface plan does not match load patches");
    }
    let blocks: BTreeMap<_, _> = space.blocks.iter().map(|b| (b.id.as_str(), b)).collect();
    let quadrature = gauss(order)?;
    let mut count = 0;
    let mut nodal = vec![
        0.;
        space
            .ndof
            .checked_mul(3)
            .ok_or("force vector dimension overflow")?
    ];
    for (patch, &depth) in patches.iter().zip(depths) {
        if depth.into_iter().any(|level| level > 30) {
            return Err("invalid frozen surface plan depth");
        }
        let block = blocks
            .get(patch.side.block.as_str())
            .ok_or("unknown traction block")?;
        let integral = composite_grid(patch_rect(patch)?, depth, budget, &mut count, |r| {
            integrate_panel(block, &patch.side, r, &quadrature, &traction)
        })?;
        for (&id, force) in block.dofs.iter().zip(&integral.nodal) {
            for (component, value) in force.iter().enumerate() {
                nodal[3 * id + component] += value;
            }
        }
    }
    if nodal.iter().any(|value| !value.is_finite()) {
        return Err("nonfinite frozen surface load rate");
    }
    Ok((nodal, count))
}

fn directional_patch(
    rect: Rect,
    first: PanelIntegral,
    config: SurfaceQuadrature,
    evaluated: &mut usize,
    sample: impl Fn(Rect) -> Result<PanelIntegral, &'static str>,
) -> Result<(PanelIntegral, [u8; 2], f64), &'static str> {
    let mut depths = [0_u8; 2];
    let mut coarse = first;
    let mut composite =
        |depths| composite_grid(rect, depths, config.max_panels, evaluated, &sample);
    loop {
        let [du, dv] = depths;
        if du >= config.max_depth || dv >= config.max_depth {
            return Err("adaptive surface quadrature tolerance not met within depth budget");
        }
        let u = composite([du + 1, dv])?;
        let v = composite([du, dv + 1])?;
        let uv = composite([du + 1, dv + 1])?;
        let eu = estimate_error(&coarse, &u)?.max(estimate_error(&v, &uv)?);
        let ev = estimate_error(&coarse, &v)?.max(estimate_error(&u, &uv)?);
        let error = eu.max(ev).max(estimate_error(&coarse, &uv)?);
        if error <= config.relative_tolerance {
            return Ok((uv, [du + 1, dv + 1], error));
        }
        if eu > 2.0 * ev {
            depths[0] += 1;
            coarse = u;
        } else if ev > 2.0 * eu {
            depths[1] += 1;
            coarse = v;
        } else {
            depths = [du + 1, dv + 1];
            coarse = uv;
        }
    }
}

/// Integrate surface tractions by face-wise dyadic composite Gauss refinement.
///
/// The returned discrepancy is an estimator based on adjacent composite levels,
/// not a certified error bound. Pressure/resultant cancellation is normalized
/// against the integrated absolute traction magnitude when the resultant is
/// smaller than `1e-12` of that physical scale.
pub fn integrate_adaptive_traction<F: Fn(Vector, Vector) -> Vector>(
    space: &AnalysisSpace,
    exterior: &BTreeSet<BlockSide>,
    patches: &[LoadPatch],
    order: usize,
    config: SurfaceQuadrature,
    traction: F,
) -> Result<AdaptiveSurfaceLoad, &'static str> {
    integrate_surface_traction(
        space,
        exterior,
        patches,
        order,
        (config, SurfaceSplit::Isotropic),
        traction,
    )
}

/// Integrate using the selected empirical refinement strategy.
pub fn integrate_surface_traction<F: Fn(Vector, Vector) -> Vector>(
    space: &AnalysisSpace,
    exterior: &BTreeSet<BlockSide>,
    patches: &[LoadPatch],
    order: usize,
    options: (SurfaceQuadrature, SurfaceSplit),
    traction: F,
) -> Result<AdaptiveSurfaceLoad, &'static str> {
    let (config, strategy) = options;
    config.validate()?;
    let quadrature = gauss(order)?;
    if patches.is_empty() {
        return Err("empty traction selection");
    }
    let mut known = BTreeMap::<&str, &AnalysisBlock>::new();
    for block in &space.blocks {
        if known.insert(&block.id, block).is_some() {
            return Err("duplicate analysis block");
        }
    }
    let force_len = space
        .ndof
        .checked_mul(3)
        .ok_or("force vector dimension overflow")?;
    let mut result = SurfaceLoad {
        nodal_force_n: vec![0.0; force_len],
        area_mm2: 0.0,
        resultant_n: [0.0; 3],
        moment_n_mm: [0.0; 3],
    };
    let mut used = BTreeMap::<BlockSide, Vec<Rect>>::new();
    let mut validated_patches = Vec::with_capacity(patches.len());
    for patch in patches {
        let face = &patch.side;
        if face.axis >= 3 || !exterior.contains(face) {
            return Err("traction requires a certified exterior side");
        }
        let block = known
            .get(face.block.as_str())
            .ok_or("unknown traction block")?;
        if block.dofs.len() != 27 || block.dofs.iter().any(|&d| d >= space.ndof) {
            return Err("invalid traction block DOFs");
        }
        let rect = patch_rect(patch)?;
        let old_rects = used.entry(face.clone()).or_default();
        if old_rects
            .iter()
            .any(|old| (0..2).all(|i| rect[i][0] < old[i][1] && old[i][0] < rect[i][1]))
        {
            return Err("overlapping traction patches");
        }
        old_rects.push(rect);
        validated_patches.push((face.clone(), *block, rect));
    }
    let mut panel_evaluations = 0usize;
    let mut deepest = 0u8;
    let mut patch_depths = Vec::with_capacity(patches.len());
    let mut maximum_error = 0.0_f64;

    for (face, block, rect) in validated_patches {
        if panel_evaluations >= config.max_panels {
            return Err("adaptive surface quadrature panel budget exceeded");
        }
        let first = integrate_panel(block, &face, rect, &quadrature, &traction)?;
        panel_evaluations = panel_evaluations
            .checked_add(1)
            .ok_or("surface panel count overflow")?;
        let (accepted, depth, error) = if strategy == SurfaceSplit::Directional {
            directional_patch(rect, first, config, &mut panel_evaluations, |r| {
                integrate_panel(block, &face, r, &quadrature, &traction)
            })?
        } else {
            let mut current = vec![(rect, first)];
            let mut current_depth = 0u8;
            loop {
                let next_len = current
                    .len()
                    .checked_mul(4)
                    .ok_or("surface panel count overflow")?;
                if next_len > config.max_panels.saturating_sub(panel_evaluations) {
                    return Err("adaptive surface quadrature panel budget exceeded");
                }
                let mut fine_panels = Vec::with_capacity(next_len);
                for (parent_rect, _) in &current {
                    for child in split(*parent_rect) {
                        let integral =
                            integrate_panel(block, &face, child, &quadrature, &traction)?;
                        fine_panels.push((child, integral));
                    }
                }
                panel_evaluations = panel_evaluations
                    .checked_add(next_len)
                    .ok_or("surface panel count overflow")?;
                let coarse_sum = sum_panels(&current);
                let fine_sum = sum_panels(&fine_panels);
                let error = estimate_error(&coarse_sum, &fine_sum)?;
                if error <= config.relative_tolerance {
                    break (fine_sum, [current_depth + 1; 2], error);
                }
                if current_depth + 1 >= config.max_depth {
                    return Err(
                        "adaptive surface quadrature tolerance not met within depth budget",
                    );
                }
                current = fine_panels;
                current_depth += 1;
            }
        };
        deepest = deepest.max(depth[0]).max(depth[1]);
        patch_depths.push(depth);
        maximum_error = maximum_error.max(error);
        result.area_mm2 += accepted.area;
        for c in 0..3 {
            result.resultant_n[c] += accepted.resultant[c];
            result.moment_n_mm[c] += accepted.moment[c];
        }
        for (&dof, force) in block.dofs.iter().zip(&accepted.nodal) {
            for (c, &component) in force.iter().enumerate() {
                result.nodal_force_n[3 * dof + c] += component;
            }
        }
    }
    if !result.area_mm2.is_finite()
        || result
            .nodal_force_n
            .iter()
            .chain(&result.resultant_n)
            .chain(&result.moment_n_mm)
            .any(|x| !x.is_finite())
    {
        return Err("nonfinite assembled adaptive surface load");
    }
    Ok(AdaptiveSurfaceLoad {
        load: result,
        panels: panel_evaluations,
        patch_depths,
        deepest,
        maximum_estimated_relative_error: maximum_error,
    })
}
