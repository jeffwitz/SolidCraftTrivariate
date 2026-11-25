//! Transactional mechanical edits of the audited parameterized bracket.
//!
//! Every candidate is regenerated and audited before its adaptive space is
//! assembled. Element stiffnesses are cached by stable leaf identity, exact
//! homogeneous controls, material, and quadrature. Global numbering is rebuilt
//! for every candidate before local CSR entries are merged.

use crate::{
    adaptive_partition::AdaptivePartition,
    adaptive_space::{build_adaptive_space, AdaptiveSpace},
    analysis_space::{from_certified_faces, AnalysisSpace},
    bracket::BracketBlock,
    bracket_parameters::build_bracket_parameters,
    bracket_parameters_audit::{audit_bracket_parameters, ParameterizedBracketCertificate},
    bracket_solve::{bracket_load, bracket_tip_displacement},
    constrained_elasticity::{solve_constrained, solve_constrained_from, ConstrainedSolution},
    constraint_system::{resolve, Lift, Ratio},
    elasticity::{assemble, gauss, Matrix, Stiffness},
};
use std::{collections::BTreeMap, mem::size_of, sync::Arc, time::Instant};

const CACHE_CAP_BYTES: usize = 384 * 1024 * 1024;
const MAX_LEVEL: u8 = 3;
const YOUNG_MPA: f64 = 210_000.0;
const POISSON: f64 = 0.3;
const QUADRATURE_ORDER: usize = 8;
const SOLVE_TOLERANCE: f64 = 1e-10;
const MAX_SOLVE_ITERATIONS: usize = 20_000;

/// Configuration fixed for one edit session.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BracketMechanicsConfig {
    pub uniform_level: u8,
    pub young_mpa: f64,
    pub poisson: f64,
    pub quadrature_order: usize,
    pub cache_cap_bytes: usize,
    pub max_cells: usize,
}

impl BracketMechanicsConfig {
    pub fn campaign(uniform_level: u8) -> Self {
        Self {
            uniform_level,
            young_mpa: YOUNG_MPA,
            poisson: POISSON,
            quadrature_order: QUADRATURE_ORDER,
            cache_cap_bytes: CACHE_CAP_BYTES,
            max_cells: 36usize.saturating_mul(8usize.saturating_pow(uniform_level as u32)),
        }
    }
}

/// Conservative preflight for all per-leaf local stiffness cache entries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CachePreflight {
    pub level: u8,
    pub cells: usize,
    pub bytes_per_entry_upper_bound: usize,
    pub total_bytes_upper_bound: usize,
    pub cap_bytes: usize,
    pub fits: bool,
}

/// Timings and cache accounting for one committed candidate.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct BracketMechanicsTiming {
    pub generation_s: f64,
    pub audit_and_space_s: f64,
    pub contribution_lookup_s: f64,
    pub local_integration_s: f64,
    pub merge_s: f64,
    pub comparison_s: f64,
    pub cold_assembly_s: f64,
    pub load_s: f64,
    pub cold_solve_s: f64,
    pub warm_solve_s: f64,
    pub reused_contributions: usize,
    pub integrated_contributions: usize,
    pub cache_entries: usize,
    pub cache_bytes_estimated: usize,
    pub candidate_cells: usize,
    pub cold_solve_iterations: usize,
    pub warm_solve_iterations: Option<usize>,
    pub warm_start_used: bool,
    pub cold_relative_residual: f64,
    pub warm_relative_residual: Option<f64>,
    pub stiffness_action_relative_error: f64,
    pub warm_cold_displacement_relative_error: Option<f64>,
    pub warm_cold_energy_relative_error: Option<f64>,
}

/// Immutable, fully audited and solved state published by a transaction.
pub struct BracketMechanicsState {
    pub radii_sixteenths_mm: [u8; 2],
    pub blocks: Vec<BracketBlock>,
    pub certificate: ParameterizedBracketCertificate,
    pub adaptive: AdaptiveSpace,
    pub stiffness: Arc<Stiffness>,
    pub force_n: Vec<f64>,
    pub load_area_mm2: f64,
    pub displacement_mm: Vec<f64>,
    pub strain_energy_n_mm: f64,
    pub tip_displacement_mm: [f64; 3],
    pub timing: BracketMechanicsTiming,
    unconstrained_lift: Lift,
    clamped_lift: Lift,
}

impl std::fmt::Debug for BracketMechanicsState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BracketMechanicsState")
            .field("radii_sixteenths_mm", &self.radii_sixteenths_mm)
            .field("ndof", &self.adaptive.space.ndof)
            .field("cells", &self.adaptive.space.blocks.len())
            .field("strain_energy_n_mm", &self.strain_energy_n_mm)
            .field("timing", &self.timing)
            .finish()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ElementKey {
    id: String,
    controls: Vec<[i128; 4]>,
    young_bits: u64,
    poisson_bits: u64,
    quadrature_order: usize,
}

struct CachedElement {
    id: String,
    controls: Vec<[i128; 4]>,
    stiffness: Stiffness,
    estimated_bytes: usize,
}

#[derive(Default)]
pub struct BracketMechanicsSession {
    config: Option<BracketMechanicsConfig>,
    state: Option<Arc<BracketMechanicsState>>,
    cache: BTreeMap<ElementKey, Arc<CachedElement>>,
    cache_bytes: usize,
}

fn key_for(id: &str, controls: &[[i128; 4]], config: BracketMechanicsConfig) -> ElementKey {
    ElementKey {
        id: id.to_owned(),
        controls: controls.to_vec(),
        young_bits: config.young_mpa.to_bits(),
        poisson_bits: config.poisson.to_bits(),
        quadrature_order: config.quadrature_order,
    }
}

fn bytes_per_entry_upper_bound() -> usize {
    // An affine Q2 element has at most 27^2 block entries. Account for CSR
    // buffers, controls, stable ID, map node and allocator slack.
    // Vec growth can round CSR capacities up (e.g. 729 entries to 1024).
    let matrix_bytes = 2 * 27 * 27 * (size_of::<Matrix>() + size_of::<usize>());
    let fixed_bytes = 28 * size_of::<usize>() + 27 * size_of::<[i128; 4]>() + 512;
    matrix_bytes.saturating_add(fixed_bytes).max(128 * 1024)
}

pub fn cache_preflight(level: u8, cap_bytes: usize) -> Result<CachePreflight, &'static str> {
    if level > MAX_LEVEL {
        return Err("bracket mechanics level exceeds supported preflight range");
    }
    let cells = 36usize
        .checked_mul(
            8usize
                .checked_pow(level as u32)
                .ok_or("cell estimate overflow")?,
        )
        .ok_or("cell estimate overflow")?;
    let bytes_per_entry_upper_bound = bytes_per_entry_upper_bound();
    let total_bytes_upper_bound = cells
        .checked_mul(bytes_per_entry_upper_bound)
        .ok_or("cache estimate overflow")?;
    Ok(CachePreflight {
        level,
        cells,
        bytes_per_entry_upper_bound,
        total_bytes_upper_bound,
        cap_bytes,
        fits: total_bytes_upper_bound <= cap_bytes,
    })
}

fn make_partition(
    blocks: &[BracketBlock],
    incidence: &crate::incidence::Incidence,
    config: BracketMechanicsConfig,
) -> Result<AdaptivePartition, &'static str> {
    let roots = blocks
        .iter()
        .map(|block| block.id.clone())
        .collect::<Vec<_>>();
    let mut partition =
        AdaptivePartition::new(&roots, incidence, config.uniform_level, config.max_cells)?;
    for _ in 0..config.uniform_level {
        let leaves = partition.snapshot()?.leaves;
        partition.refine(&leaves)?;
    }
    Ok(partition)
}

fn mounted_lifts(adaptive: &AdaptiveSpace) -> Result<(Lift, Lift), &'static str> {
    let unconstrained = resolve(adaptive.space.ndof, &adaptive.constraints, &[])?;
    let fixed = adaptive
        .space
        .boundary_dofs
        .get("y50")
        .ok_or("adaptive bracket has no y50 Dirichlet feature")?
        .iter()
        .map(|&dof| Ok((dof, Ratio::zero())))
        .collect::<Result<Vec<_>, &'static str>>()?;
    let clamped = resolve(adaptive.space.ndof, &adaptive.constraints, &fixed)?;
    Ok((unconstrained, clamped))
}

fn transfer_guess(
    previous: &BracketMechanicsState,
    adaptive: &AdaptiveSpace,
    unconstrained: &Lift,
    clamped: &Lift,
) -> Option<Vec<f64>> {
    // Geometry can change while leaf identities and the equation/essential
    // maps remain identical. Only then are prior full coefficients valid
    // initial values in the new reduced basis.
    if previous.adaptive.space.ndof != adaptive.space.ndof
        || previous.adaptive.lift != adaptive.lift
        || previous.unconstrained_lift != *unconstrained
        || previous.clamped_lift != *clamped
        || previous
            .adaptive
            .space
            .blocks
            .iter()
            .map(|block| (&block.id, &block.dofs))
            .ne(adaptive
                .space
                .blocks
                .iter()
                .map(|block| (&block.id, &block.dofs)))
    {
        return None;
    }
    let n_free = clamped.free_roots.len();
    let mut guess = vec![0.0; 3 * n_free];
    for component in 0..3 {
        for (index, &dof) in clamped.free_roots.iter().enumerate() {
            guess[component * n_free + index] = previous.displacement_mm[3 * dof + component];
        }
    }
    (guess.len() == 3 * n_free && guess.iter().all(|value| value.is_finite())).then_some(guess)
}

fn relative_error(a: &[f64], b: &[f64]) -> Result<f64, &'static str> {
    if a.len() != b.len() || a.iter().chain(b).any(|value| !value.is_finite()) {
        return Err("invalid vectors for relative comparison");
    }
    let numerator = a
        .iter()
        .zip(b)
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt();
    let denominator = b.iter().map(|value| value * value).sum::<f64>().sqrt();
    Ok(numerator / denominator.max(f64::MIN_POSITIVE))
}

fn deterministic_vector(size: usize) -> Vec<f64> {
    (0..size)
        .map(|index| ((index.wrapping_mul(37) % 101) as f64 - 50.0) / 31.0)
        .collect()
}

fn action_error(a: &Stiffness, b: &Stiffness) -> Result<f64, &'static str> {
    let vector = deterministic_vector(a.size());
    let av = a.apply(&vector)?;
    let bv = b.apply(&vector)?;
    relative_error(&av, &bv)
}

fn merge_contributions(
    adaptive: &AdaptiveSpace,
    contributions: &BTreeMap<String, Arc<CachedElement>>,
) -> Result<Stiffness, &'static str> {
    let mut rows = vec![BTreeMap::<usize, Matrix>::new(); adaptive.space.ndof];
    let mut integrated_volume = 0.0;
    for block in &adaptive.space.blocks {
        let element = contributions
            .get(&block.id)
            .ok_or("missing cached bracket contribution")?;
        if element.id != block.id || element.controls != block.controls {
            return Err("bracket contribution key collision");
        }
        let local = &element.stiffness;
        integrated_volume += local.integrated_volume;
        for (local_row, &global_row) in block.dofs.iter().enumerate() {
            for entry in local.offsets[local_row]..local.offsets[local_row + 1] {
                let global_column = block.dofs[local.columns[entry]];
                let target = rows[global_row].entry(global_column).or_default();
                for (target_row, source_row) in target.iter_mut().zip(&local.values[entry]) {
                    for (target_value, source_value) in target_row.iter_mut().zip(source_row) {
                        *target_value += source_value;
                    }
                }
            }
        }
    }
    let mut merged = Stiffness {
        offsets: vec![0],
        columns: Vec::new(),
        values: Vec::new(),
        integrated_volume,
    };
    for row in rows {
        for (column, value) in row {
            merged.columns.push(column);
            merged.values.push(value);
        }
        merged.offsets.push(merged.columns.len());
    }
    merged.validate()?;
    if !integrated_volume.is_finite() || integrated_volume <= 0.0 {
        return Err("invalid cached bracket volume");
    }
    Ok(merged)
}

fn reduced_residual(
    stiffness: &Stiffness,
    force: &[f64],
    solution: &ConstrainedSolution,
    lifts: [&Lift; 3],
) -> Result<f64, &'static str> {
    let residual = stiffness.equilibrium_residual(&solution.displacement_mm, force)?;
    let mut projected = Vec::new();
    for component in 0..3 {
        for root in 0..lifts[component].free_roots.len() {
            let mut value = 0.0;
            for line in &lifts[component].lines {
                for &(master, ratio) in &line.terms {
                    if master == root {
                        value += ratio.numerator() as f64 / ratio.denominator() as f64
                            * residual[3 * line.dof + component];
                    }
                }
            }
            projected.push(value);
        }
    }
    let rhs_norm = force.iter().map(|value| value * value).sum::<f64>().sqrt();
    let residual_norm = projected
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    if rhs_norm <= 0.0 || !rhs_norm.is_finite() || !residual_norm.is_finite() {
        return Err("invalid true bracket residual");
    }
    Ok(residual_norm / rhs_norm)
}

impl BracketMechanicsSession {
    pub fn new(config: BracketMechanicsConfig) -> Result<Self, &'static str> {
        validate_config(config)?;
        let preflight = cache_preflight(config.uniform_level, config.cache_cap_bytes)?;
        if !preflight.fits {
            return Err("cache budget preflight exceeded");
        }
        Ok(Self {
            config: Some(config),
            state: None,
            cache: BTreeMap::new(),
            cache_bytes: 0,
        })
    }

    pub fn state(&self) -> Option<&Arc<BracketMechanicsState>> {
        self.state.as_ref()
    }

    pub fn cache_bytes_estimated(&self) -> usize {
        self.cache_bytes
    }

    pub fn cache_entries(&self) -> usize {
        self.cache.len()
    }

    /// Audit, assemble, and solve a candidate; publish state/cache atomically.
    pub fn update(
        &mut self,
        radii_sixteenths_mm: [u8; 2],
    ) -> Result<&Arc<BracketMechanicsState>, &'static str> {
        let config = self.config.ok_or("uninitialized bracket session")?;
        let preflight = cache_preflight(config.uniform_level, config.cache_cap_bytes)?;
        if !preflight.fits {
            return Err("cache budget preflight exceeded");
        }
        let mut timing = BracketMechanicsTiming::default();

        let started = Instant::now();
        let parameterized =
            build_bracket_parameters(radii_sixteenths_mm[0], radii_sixteenths_mm[1])?;
        let blocks = parameterized.blocks;
        timing.generation_s = started.elapsed().as_secs_f64();

        let started = Instant::now();
        let certificate =
            audit_bracket_parameters(radii_sixteenths_mm[0], radii_sixteenths_mm[1], &blocks)?;
        let root = from_certified_faces(
            &blocks,
            &certificate.certificate.incidence,
            &certificate.certificate.design_faces,
        )?;
        let partition = make_partition(&blocks, &certificate.certificate.incidence, config)?;
        let adaptive = build_adaptive_space(&root, &partition)?;
        let (unconstrained_lift, clamped_lift) = mounted_lifts(&adaptive)?;
        timing.audit_and_space_s = started.elapsed().as_secs_f64();

        let missing = adaptive
            .space
            .blocks
            .iter()
            .filter(|block| {
                let key = key_for(&block.id, &block.controls, config);
                !self.cache.get(&key).is_some_and(|cached| {
                    cached.id == block.id && cached.controls == block.controls
                })
            })
            .collect::<Vec<_>>();
        let prospective_bytes = missing
            .len()
            .checked_mul(bytes_per_entry_upper_bound())
            .and_then(|bytes| self.cache_bytes.checked_add(bytes))
            .ok_or("cache estimate overflow")?;
        if prospective_bytes > config.cache_cap_bytes {
            return Err("cache budget preflight exceeded");
        }

        let started = Instant::now();
        let mut candidate = BTreeMap::<String, Arc<CachedElement>>::new();
        for block in &adaptive.space.blocks {
            let key = key_for(&block.id, &block.controls, config);
            if let Some(cached) = self
                .cache
                .get(&key)
                .filter(|cached| cached.id == block.id && cached.controls == block.controls)
            {
                timing.reused_contributions += 1;
                candidate.insert(block.id.clone(), Arc::clone(cached));
            }
        }
        timing.contribution_lookup_s = started.elapsed().as_secs_f64();

        let started = Instant::now();
        let mut staged = Vec::<Arc<CachedElement>>::with_capacity(missing.len());
        for block in missing {
            let local = AnalysisSpace {
                ndof: 27,
                blocks: vec![crate::analysis_space::AnalysisBlock {
                    id: block.id.clone(),
                    controls: block.controls.clone(),
                    dofs: (0..27).collect(),
                }],
                boundary_dofs: BTreeMap::new(),
            };
            let stiffness = assemble(
                &local,
                config.young_mpa,
                config.poisson,
                config.quadrature_order,
            )?;
            let estimated_bytes = stiffness.offsets.capacity() * size_of::<usize>()
                + stiffness.columns.capacity() * size_of::<usize>()
                + stiffness.values.capacity() * size_of::<Matrix>()
                + block.controls.capacity() * size_of::<[i128; 4]>()
                + block.id.capacity()
                + block.controls.len() * size_of::<[i128; 4]>()
                + block.id.len()
                + 512;
            if estimated_bytes > bytes_per_entry_upper_bound() {
                return Err("local stiffness cache entry exceeds preflight bound");
            }
            let element = Arc::new(CachedElement {
                id: block.id.clone(),
                controls: block.controls.clone(),
                stiffness,
                estimated_bytes,
            });
            staged.push(Arc::clone(&element));
            candidate.insert(block.id.clone(), element);
        }
        timing.integrated_contributions = staged.len();
        timing.local_integration_s = started.elapsed().as_secs_f64();

        let started = Instant::now();
        let merged = merge_contributions(&adaptive, &candidate)?;
        timing.merge_s = started.elapsed().as_secs_f64();

        let started = Instant::now();
        let cold = assemble(
            &adaptive.space,
            config.young_mpa,
            config.poisson,
            config.quadrature_order,
        )?;
        timing.cold_assembly_s = started.elapsed().as_secs_f64();
        let started = Instant::now();
        timing.stiffness_action_relative_error = action_error(&merged, &cold)?;
        if timing.stiffness_action_relative_error > 1e-11 {
            return Err("incremental stiffness differs from cold assembly");
        }
        timing.comparison_s = started.elapsed().as_secs_f64();

        let started = Instant::now();
        let (force_n, load_area_mm2) = bracket_load(&adaptive.space, config.quadrature_order)?;
        if (load_area_mm2 - 320.0).abs() > 1e-8 {
            return Err("audited bracket load area differs from 320 mm^2");
        }
        timing.load_s = started.elapsed().as_secs_f64();

        // `y50` is a vector Dirichlet condition: all three displacement
        // components are prescribed on the selected face.
        let lifts = [&clamped_lift, &clamped_lift, &clamped_lift];
        let started = Instant::now();
        let cold_solution = solve_constrained(
            &cold,
            lifts,
            &force_n,
            SOLVE_TOLERANCE,
            MAX_SOLVE_ITERATIONS,
        )
        .map_err(|_| "cold bracket solve failed")?;
        timing.cold_solve_s = started.elapsed().as_secs_f64();
        timing.cold_solve_iterations = cold_solution.iterations;
        timing.cold_relative_residual = reduced_residual(&cold, &force_n, &cold_solution, lifts)?;
        if timing.cold_relative_residual > 1e-10 {
            return Err("cold bracket solve true residual exceeds tolerance");
        }

        let previous = self.state.as_deref();
        let warm_guess = previous
            .and_then(|state| transfer_guess(state, &adaptive, &unconstrained_lift, &clamped_lift));
        timing.warm_start_used = warm_guess.is_some();
        let mut warm_solution = None;
        if let Some(initial) = warm_guess {
            let started = Instant::now();
            let solved = solve_constrained_from(
                &cold,
                lifts,
                &force_n,
                SOLVE_TOLERANCE,
                MAX_SOLVE_ITERATIONS,
                &initial,
            )
            .map_err(|_| "warm bracket solve failed")?;
            timing.warm_solve_s = started.elapsed().as_secs_f64();
            timing.warm_solve_iterations = Some(solved.iterations);
            timing.warm_relative_residual =
                Some(reduced_residual(&merged, &force_n, &solved, lifts)?);
            if timing.warm_relative_residual.unwrap() > 1e-10 {
                return Err("warm bracket solve true residual exceeds tolerance");
            }
            let displacement_error =
                relative_error(&solved.displacement_mm, &cold_solution.displacement_mm)?;
            let energy_error = (solved.strain_energy_n_mm - cold_solution.strain_energy_n_mm).abs()
                / cold_solution
                    .strain_energy_n_mm
                    .abs()
                    .max(f64::MIN_POSITIVE);
            timing.warm_cold_displacement_relative_error = Some(displacement_error);
            timing.warm_cold_energy_relative_error = Some(energy_error);
            if displacement_error > 1e-8 {
                return Err("warm and cold bracket displacements disagree");
            }
            if energy_error > 1e-8 {
                return Err("warm and cold bracket energies disagree");
            }
            warm_solution = Some(solved);
        }
        let accepted_solution = warm_solution.as_ref().unwrap_or(&cold_solution);
        let tip_displacement_mm =
            bracket_tip_displacement(&adaptive.space, &accepted_solution.displacement_mm)?;

        timing.candidate_cells = adaptive.space.blocks.len();
        timing.cache_entries = self.cache.len() + staged.len();
        timing.cache_bytes_estimated = self.cache_bytes
            + staged
                .iter()
                .map(|entry| entry.estimated_bytes)
                .sum::<usize>();
        if timing.cache_bytes_estimated > config.cache_cap_bytes {
            return Err("measured contribution cache exceeds configured cap");
        }

        // No mutation to session state or cache has occurred before this point.
        for element in staged {
            let key = key_for(&element.id, &element.controls, config);
            self.cache_bytes = self.cache_bytes.saturating_add(element.estimated_bytes);
            self.cache.insert(key, element);
        }
        let stiffness = Arc::new(merged);
        let state = Arc::new(BracketMechanicsState {
            radii_sixteenths_mm,
            blocks,
            certificate,
            adaptive,
            stiffness,
            force_n,
            load_area_mm2,
            displacement_mm: accepted_solution.displacement_mm.clone(),
            strain_energy_n_mm: accepted_solution.strain_energy_n_mm,
            tip_displacement_mm,
            timing,
            unconstrained_lift,
            clamped_lift,
        });
        self.state = Some(state);
        Ok(self.state.as_ref().expect("state committed"))
    }
}

fn validate_config(config: BracketMechanicsConfig) -> Result<(), &'static str> {
    if config.uniform_level > MAX_LEVEL {
        return Err("bracket mechanics level exceeds supported range");
    }
    if !config.young_mpa.is_finite()
        || config.young_mpa <= 0.0
        || !config.poisson.is_finite()
        || !(-1.0 < config.poisson && config.poisson < 0.5)
    {
        return Err("invalid bracket mechanics material");
    }
    gauss(config.quadrature_order)?;
    let min_cells = 36usize
        .checked_mul(
            8usize
                .checked_pow(config.uniform_level as u32)
                .ok_or("cell estimate overflow")?,
        )
        .ok_or("cell estimate overflow")?;
    if config.max_cells < min_cells {
        return Err("adaptive cell budget is below uniform mesh size");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cache_preflight_enforces_cap_before_session_creation() {
        let h2 = cache_preflight(2, CACHE_CAP_BYTES).unwrap();
        assert!(h2.fits);
        let h3 = cache_preflight(3, CACHE_CAP_BYTES).unwrap();
        assert_eq!(h3.cells, 18_432);
        assert!(!h3.fits);
        assert!(BracketMechanicsSession::new(BracketMechanicsConfig {
            uniform_level: 3,
            ..BracketMechanicsConfig::campaign(3)
        })
        .is_err());
    }

    #[test]
    fn parameter_rejection_preserves_published_state_and_cache() {
        let mut session = BracketMechanicsSession::new(BracketMechanicsConfig {
            uniform_level: 0,
            quadrature_order: 8,
            ..BracketMechanicsConfig::campaign(0)
        })
        .unwrap();
        session.update([50, 74]).unwrap();
        let state = Arc::clone(session.state().unwrap());
        let bytes = session.cache_bytes_estimated();
        assert!(session.update([81, 74]).is_err());
        assert!(Arc::ptr_eq(session.state().unwrap(), &state));
        assert_eq!(session.cache_bytes_estimated(), bytes);
    }

    #[test]
    fn changed_geometry_reuses_unchanged_elements_and_matches_cold_k() {
        let mut session = BracketMechanicsSession::new(BracketMechanicsConfig {
            uniform_level: 0,
            quadrature_order: 8,
            ..BracketMechanicsConfig::campaign(0)
        })
        .unwrap();
        session.update([50, 74]).unwrap();
        let changed = session.update([51, 74]).unwrap();
        assert!(changed.timing.reused_contributions > 0);
        assert!(changed.timing.integrated_contributions > 0);
        assert!(changed.timing.stiffness_action_relative_error < 1e-11);
        assert!(
            changed
                .timing
                .warm_cold_displacement_relative_error
                .unwrap()
                < 1e-8
        );
        assert!(changed.timing.warm_cold_energy_relative_error.unwrap() < 1e-8);
    }
}
