//! Shape derivative of a rational sector at fixed chart topology.
//! Whole-face matching shares field coefficients directly, so its unrefined
//! assembly has a constant lift even though the field basis weights vary.
use crate::{
    analysis_space::{AnalysisBlock, AnalysisSpace},
    boundary_traction::LoadPatch,
    constraint_system::Lift,
    elasticity::{gauss, Vector},
    elasticity_sensitivity::evaluate_geometry_jet,
    mechanical_job::Traction,
    mechanical_sensitivity::LiftRate,
    rational_sector::RationalSector,
};

pub(crate) struct SectorMotion(pub RationalSector);

impl SectorMotion {
    fn chart_and_leaf_bounds(&self, id: &str) -> Result<[(f64, f64); 3], &'static str> {
        let mut parts = id.split('/');
        let chart = parts.next().ok_or("invalid sector leaf identity")?;
        let [n, d] = self.0.half_angle_tangent;
        let geometric_split = n > d && i128::from(n - d).pow(2) > 2 * i128::from(d).pow(2);
        let mut bounds = match (geometric_split, chart) {
            (false, "whole") => [(0., 1.), (0., 1.), (0., 1.)],
            (_, "first") => [(0., 0.5), (0., 1.), (0., 1.)],
            (_, "second") => [(0.5, 1.), (0., 1.), (0., 1.)],
            _ => return Err("sector derivative chart does not match the fixed regime"),
        };
        for path in parts {
            let bytes = path.as_bytes();
            if bytes.len() != 3 || bytes.iter().any(|b| !matches!(b, b'0' | b'1')) {
                return Err("invalid dyadic sector leaf identity");
            }
            for axis in 0..3 {
                let (lo, hi) = bounds[axis];
                let middle = (lo + hi) * 0.5;
                bounds[axis] = if bytes[axis] == b'0' {
                    (lo, middle)
                } else {
                    (middle, hi)
                };
            }
        }
        Ok(bounds)
    }

    pub(crate) fn jets(&self, block: &AnalysisBlock) -> Result<Vec<[f64; 4]>, &'static str> {
        let [n, d] = self.0.half_angle_tangent;
        let t = f64::from(n) / f64::from(d);
        let bounds = self.chart_and_leaf_bounds(&block.id)?;
        let (a, b) = bounds[0];
        let [ri, ro] = self.0.radii_mm.map(f64::from);
        if block.controls.len() != 27 {
            return Err("sector derivative requires Q2 space");
        }
        let mut jets = Vec::with_capacity(27);
        for (idx, p) in block.controls.iter().enumerate() {
            let angular = idx / 9;
            let [l, r] = match angular {
                0 => [a, a],
                1 => [a, b],
                _ => [b, b],
            };
            let v = bounds[1].0 + (bounds[1].1 - bounds[1].0) * (idx / 3 % 3) as f64 / 2.;
            let w = bounds[2].0 + (bounds[2].1 - bounds[2].0) * (idx % 3) as f64 / 2.;
            let radius = ro + (ri - ro) * v;
            let z = f64::from(self.0.height_mm) * w;
            if p[3] <= 0 {
                return Err("sector derivative requires positive homogeneous weights");
            }
            let scale = p[3] as f64 / (1. + t * t * l * r);
            let dw = 2. * t * l * r;
            jets.push([-radius * dw, radius * (l + r), z * dw, dw].map(|v| v * scale));
        }
        Ok(jets)
    }

    /// Differentiate the exact rational field lift on an unbalanced Q2 space.
    /// Local homogeneous gauges cancel in each logarithmic weight rate.
    pub(crate) fn lift_rates(
        &self,
        space: &AnalysisSpace,
        lift: &Lift,
    ) -> Result<Vec<LiftRate>, &'static str> {
        if space.ndof == 0 || lift.lines.len() != space.ndof {
            return Err("sector lift-rate dimension mismatch");
        }
        let mut log_weights: Vec<Option<f64>> = vec![None; space.ndof];
        for block in &space.blocks {
            if block.controls.len() != 27 || block.dofs.len() != 27 {
                return Err("sector lift-rate requires numbered Q2 blocks");
            }
            let jets = self.jets(block)?;
            for ((control, jet), &dof) in block.controls.iter().zip(&jets).zip(&block.dofs) {
                if dof >= space.ndof || control[3] <= 0 {
                    return Err("invalid sector homogeneous DOF representation");
                }
                let rate = jet[3] / control[3] as f64;
                if !rate.is_finite() {
                    return Err("nonfinite sector logarithmic weight rate");
                }
                if let Some(previous) = log_weights[dof] {
                    let tolerance = 128. * f64::EPSILON * previous.abs().max(rate.abs()).max(1.);
                    if (previous - rate).abs() > tolerance {
                        return Err("inconsistent sector weight rate at shared DOF");
                    }
                } else {
                    log_weights[dof] = Some(rate);
                }
            }
        }
        let log_weights = log_weights
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .ok_or("sector analysis DOF lacks a homogeneous representation")?;
        if lift.free_roots.iter().any(|&dof| dof >= space.ndof)
            || lift
                .lines
                .iter()
                .enumerate()
                .any(|(dof, line)| line.dof != dof)
        {
            return Err("sector lift rates require a canonical lift");
        }
        let mut result = Vec::with_capacity(space.ndof);
        for line in &lift.lines {
            let local_rate = log_weights[line.dof];
            let mut terms = Vec::with_capacity(line.terms.len());
            for &(root_index, coefficient) in &line.terms {
                let root_dof = *lift
                    .free_roots
                    .get(root_index)
                    .ok_or("sector lift references an invalid free root")?;
                let value = coefficient.numerator() as f64 / coefficient.denominator() as f64;
                let derivative = value * (log_weights[root_dof] - local_rate);
                if !derivative.is_finite() {
                    return Err("nonfinite sector lift coefficient rate");
                }
                terms.push((root_index, derivative));
            }
            result.push(LiftRate { terms, offset: 0. });
        }
        Ok(result)
    }

    /// Differentiate the lift of a spatially affine prescribed field on moving
    /// sector geometry, including its geometry-dependent affine offset.
    pub(crate) fn lift_rates_with_motion(
        &self,
        space: &AnalysisSpace,
        lift: &Lift,
        component: usize,
        motion: crate::prescribed_motion::Motion,
    ) -> Result<Vec<LiftRate>, &'static str> {
        if component >= 3 || lift.lines.len() != space.ndof {
            return Err("invalid moving sector support component or lift size");
        }
        let rates = self.lift_rates(space, lift)?;
        let affine = motion.as_f64()[component];
        let mut field_values: Vec<Option<crate::constraint_system::Ratio>> = vec![None; space.ndof];
        let mut field_rates: Vec<Option<f64>> = vec![None; space.ndof];
        let mut positions: Vec<Option<[i128; 4]>> = vec![None; space.ndof];
        for block in &space.blocks {
            if block.controls.len() != 27 || block.dofs.len() != 27 {
                return Err("moving sector support requires numbered Q2 blocks");
            }
            let jets = self.jets(block)?;
            for ((p, dp), &dof) in block.controls.iter().zip(&jets).zip(&block.dofs) {
                if dof >= space.ndof || p[3] <= 0 {
                    return Err("invalid moving sector support control");
                }
                if let Some(position) = positions[dof] {
                    for axis in 0..3 {
                        if crate::constraint_system::Ratio::new(position[axis], position[3])?
                            != crate::constraint_system::Ratio::new(p[axis], p[3])?
                        {
                            return Err("inconsistent physical point at shared sector DOF");
                        }
                    }
                } else {
                    positions[dof] = Some(*p);
                }
                let prescribed = motion.value(component, *p)?;
                if field_values[dof].is_some_and(|previous| previous != prescribed) {
                    return Err("inconsistent affine support value at shared sector DOF");
                }
                field_values[dof] = Some(prescribed);
                let weight = p[3] as f64;
                let x = std::array::from_fn::<_, 3, _>(|axis| p[axis] as f64 / weight);
                let dx = std::array::from_fn::<_, 3, _>(|axis| {
                    dp[axis] / weight - x[axis] * dp[3] / weight
                });
                let rate = (0..3).map(|axis| affine[axis + 1] * dx[axis]).sum::<f64>();
                if !rate.is_finite() {
                    return Err("nonfinite material support value rate");
                }
                if let Some(previous) = field_rates[dof] {
                    let tolerance = 128. * f64::EPSILON * previous.abs().max(rate.abs()).max(1.);
                    if (previous - rate).abs() > tolerance {
                        return Err("inconsistent material velocity at shared sector DOF");
                    }
                } else {
                    field_rates[dof] = Some(rate);
                }
            }
        }
        let field_values = field_values
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .ok_or("sector DOF has no affine support value")?;
        let field_rates = field_rates
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .ok_or("sector DOF has no material support rate")?;
        let mut output = rates;
        for line in &lift.lines {
            let mut expected_offset = field_values[line.dof];
            let mut derivative_offset = field_rates[line.dof];
            for &(root_index, coefficient) in &line.terms {
                let root = *lift
                    .free_roots
                    .get(root_index)
                    .ok_or("sector lift references invalid free root")?;
                let value = coefficient.checked_mul(field_values[root])?;
                expected_offset = expected_offset.checked_sub(value)?;
                let c = coefficient.numerator() as f64 / coefficient.denominator() as f64;
                let c_rate = output[line.dof]
                    .terms
                    .iter()
                    .find(|(j, _)| *j == root_index)
                    .map_or(0., |(_, v)| *v);
                derivative_offset -= c_rate * field_values[root].numerator() as f64
                    / field_values[root].denominator() as f64
                    + c * field_rates[root];
            }
            if expected_offset != line.offset {
                return Err("prescribed affine lift offset is inconsistent with the motion field");
            }
            if !derivative_offset.is_finite() {
                return Err("nonfinite moving sector lift offset rate");
            }
            output[line.dof].offset = derivative_offset;
        }
        Ok(output)
    }

    pub(crate) fn surface_force_rate(
        &self,
        space: &AnalysisSpace,
        patches: &[LoadPatch],
        order: usize,
        law: &Traction,
    ) -> Result<Vec<f64>, &'static str> {
        let rule = gauss(order)?;
        let mut rate = vec![0.; 3 * space.ndof];
        for patch in patches {
            if patch.depth != 0 || patch.intervals != [[0, 1]; 2] {
                return Err("sector derivative requires whole faces");
            }
            let side = &patch.side;
            let block = space
                .blocks
                .iter()
                .find(|b| b.id == side.block)
                .ok_or("unknown sector load block")?;
            let jets = self.jets(block)?;
            let [a, b] = match side.axis {
                0 => [1, 2],
                1 => [0, 2],
                2 => [0, 1],
                _ => return Err("invalid sector load side"),
            };
            let sign = if side.upper ^ (side.axis == 1) {
                1.
            } else {
                -1.
            };
            for &(u, wu) in &rule {
                for &(v, wv) in &rule {
                    let mut uvw = [0.; 3];
                    uvw[side.axis] = f64::from(side.upper);
                    uvw[a] = u;
                    uvw[b] = v;
                    let g = evaluate_geometry_jet(block, &jets, uvw)?;
                    let j = g.jacobian_columns;
                    let jp = g.jacobian_column_rates;
                    let area = cross(j[a], j[b]).map(|x| sign * x);
                    let da = cross(jp[a], j[b]);
                    let db = cross(j[a], jp[b]);
                    let area_rate: Vector = std::array::from_fn(|k| sign * (da[k] + db[k]));
                    let measure = dot(area, area).sqrt();
                    if !measure.is_finite() || measure <= 0. {
                        return Err("invalid sector face area");
                    }
                    let normal = area.map(|x| x / measure);
                    let density = law.value(normal).map(|x| x * measure);
                    let density_rate = match law {
                        Traction::Vector(v) => v.map(|x| x * dot(normal, area_rate)),
                        _ => law.value(area_rate),
                    };
                    for (index, &dof) in block.dofs.iter().enumerate() {
                        for k in 0..3 {
                            rate[3 * dof + k] += wu
                                * wv
                                * (g.basis_derivative[index] * density[k]
                                    + g.basis[index] * density_rate[k]);
                        }
                    }
                }
            }
        }
        Ok(rate)
    }
}

fn cross(a: Vector, b: Vector) -> Vector {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn dot(a: Vector, b: Vector) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        adaptive_partition::{AdaptivePartition, LeafKey},
        adaptive_space::build_adaptive_space,
        analysis_space::from_certified_faces,
        rational_sector::{build_rational_sector, RationalSector},
        rational_sector_audit::audit_rational_sector,
    };

    fn data(
        t: [i32; 2],
    ) -> (
        crate::analysis_space::AnalysisSpace,
        crate::incidence::Incidence,
    ) {
        let recipe = RationalSector {
            half_angle_tangent: t,
            radii_mm: [2, 6],
            height_mm: 8,
        };
        let blocks = build_rational_sector(&recipe).unwrap();
        let certificate = audit_rational_sector(&recipe, &blocks).unwrap();
        let space =
            from_certified_faces(&blocks, &certificate.incidence, &certificate.design_faces)
                .unwrap();
        (space, certificate.incidence)
    }

    fn coefficient(line: &crate::constraint_system::LiftLine, index: usize) -> f64 {
        line.terms
            .iter()
            .find(|(j, _)| *j == index)
            .map_or(0., |(_, r)| r.numerator() as f64 / r.denominator() as f64)
    }

    #[test]
    fn dyadic_leaf_weight_jets_match_the_full_rebuilt_lift() {
        let (parent, incidence) = data([1, 2]);
        let mut partition = AdaptivePartition::new(
            &parent
                .blocks
                .iter()
                .map(|b| b.id.clone())
                .collect::<Vec<_>>(),
            &incidence,
            2,
            128,
        )
        .unwrap();
        partition.refine(&[LeafKey::root("whole")]).unwrap();
        partition
            .refine(&[LeafKey {
                root: "whole".into(),
                level: 1,
                index: [0, 0, 0],
            }])
            .unwrap();
        let adaptive = build_adaptive_space(&parent, &partition).unwrap();
        let rates = SectorMotion(RationalSector {
            half_angle_tangent: [1, 2],
            radii_mm: [2, 6],
            height_mm: 8,
        })
        .lift_rates(&adaptive.space, &adaptive.lift)
        .unwrap();
        assert!(
            rates
                .iter()
                .flat_map(|r| &r.terms)
                .any(|(_, x)| x.abs() > 1e-8),
            "local refinement must create a nonzero C' coefficient"
        );

        let step = 1. / 510.;
        let (plus_parent, _) = data([128, 255]);
        let (minus_parent, _) = data([127, 255]);
        let plus = build_adaptive_space(&plus_parent, &partition).unwrap();
        let minus = build_adaptive_space(&minus_parent, &partition).unwrap();
        assert_eq!(adaptive.lift.free_roots, plus.lift.free_roots);
        assert_eq!(adaptive.lift.free_roots, minus.lift.free_roots);
        let motion = SectorMotion(RationalSector {
            half_angle_tangent: [1, 2],
            radii_mm: [2, 6],
            height_mm: 8,
        });
        let plus_blocks = plus
            .space
            .blocks
            .iter()
            .map(|b| (b.id.as_str(), b))
            .collect::<std::collections::BTreeMap<_, _>>();
        let minus_blocks = minus
            .space
            .blocks
            .iter()
            .map(|b| (b.id.as_str(), b))
            .collect::<std::collections::BTreeMap<_, _>>();
        for block in &adaptive.space.blocks {
            let jets = motion.jets(block).unwrap();
            let bp = plus_blocks[block.id.as_str()];
            let bm = minus_blocks[block.id.as_str()];
            for (i, jet) in jets.iter().enumerate() {
                for c in 0..3 {
                    let p = block.controls[i];
                    let point = p[c] as f64 / p[3] as f64;
                    let analytic = jet[c] / p[3] as f64 - point * jet[3] / p[3] as f64;
                    let fd = (bp.controls[i][c] as f64 / bp.controls[i][3] as f64
                        - bm.controls[i][c] as f64 / bm.controls[i][3] as f64)
                        / (2. * step);
                    assert!(
                        (analytic - fd).abs() < 1e-4,
                        "{} control {i} coordinate {c}: jet {analytic}, fd {fd}",
                        block.id
                    );
                }
            }
        }
        for (i, rate) in rates.iter().enumerate() {
            for (j, value) in &rate.terms {
                let finite_difference = (coefficient(&plus.lift.lines[i], *j)
                    - coefficient(&minus.lift.lines[i], *j))
                    / (2. * step);
                assert!(
                    (finite_difference - value).abs() < 5e-5,
                    "row {i}, master {j}: analytic {value}, finite difference {finite_difference}"
                );
            }
        }
    }

    #[test]
    fn unrefined_whole_face_lift_has_zero_constraint_rate() {
        let (parent, incidence) = data([1, 2]);
        let partition = AdaptivePartition::new(
            &parent
                .blocks
                .iter()
                .map(|b| b.id.clone())
                .collect::<Vec<_>>(),
            &incidence,
            1,
            8,
        )
        .unwrap();
        let adaptive = build_adaptive_space(&parent, &partition).unwrap();
        let rates = SectorMotion(RationalSector {
            half_angle_tangent: [1, 2],
            radii_mm: [2, 6],
            height_mm: 8,
        })
        .lift_rates(&adaptive.space, &adaptive.lift)
        .unwrap();
        assert!(rates
            .iter()
            .flat_map(|r| &r.terms)
            .all(|(_, x)| x.abs() < 1e-14));
    }

    #[test]
    fn fixed_analysis_halves_are_valid_below_the_geometric_split_threshold() {
        let recipe = RationalSector {
            half_angle_tangent: [1, 2],
            radii_mm: [2, 6],
            height_mm: 8,
        };
        let [whole] = build_rational_sector(&recipe).unwrap().try_into().unwrap();
        let charts = crate::rational_sector::split_for_analysis(&whole).unwrap();
        let incidence = crate::incidence::enumerate_faces(&charts).unwrap();
        let space = from_certified_faces(&charts, &incidence, &[]).unwrap();
        let motion = SectorMotion(recipe);
        assert_eq!(
            space
                .blocks
                .iter()
                .map(|b| b.id.as_str())
                .collect::<Vec<_>>(),
            ["first", "second"]
        );
        for block in &space.blocks {
            let jets = motion.jets(block).unwrap();
            assert_eq!(jets.len(), 27);
            assert!(jets.iter().flatten().all(|v| v.is_finite()));
        }
    }

    fn moving_support_lift(
        recipe: &RationalSector,
        partition: &AdaptivePartition,
        motion: crate::prescribed_motion::Motion,
    ) -> (
        crate::adaptive_space::AdaptiveSpace,
        crate::constraint_system::Lift,
    ) {
        let (parent, _) = data(recipe.half_angle_tangent);
        let adaptive = build_adaptive_space(&parent, partition).unwrap();
        let mut prescribed = std::collections::BTreeMap::new();
        for side in &adaptive.boundary_patches["angle_start"] {
            let id = &adaptive.leaf_blocks[&side.leaf];
            let block = adaptive.space.blocks.iter().find(|b| &b.id == id).unwrap();
            for v in 0..=2 {
                for w in 0..=2 {
                    let index = v * 3 + w;
                    let dof = block.dofs[index];
                    let value = motion.value(0, block.controls[index]).unwrap();
                    if let Some(previous) = prescribed.insert(dof, value) {
                        assert_eq!(previous, value);
                    }
                }
            }
        }
        let prescribed = prescribed.into_iter().collect::<Vec<_>>();
        let lift = crate::constraint_system::impose_on_lift(&adaptive.lift, &prescribed).unwrap();
        (adaptive, lift)
    }

    #[test]
    fn moving_affine_support_lift_rates_match_rebuilt_lifts() {
        use crate::{
            json::{object, Json},
            prescribed_motion::Motion,
        };

        let recipe = RationalSector {
            half_angle_tangent: [1, 2],
            radii_mm: [2, 6],
            height_mm: 8,
        };
        let (parent, incidence) = data([1, 2]);
        let mut partition = AdaptivePartition::new(
            &parent
                .blocks
                .iter()
                .map(|b| b.id.clone())
                .collect::<Vec<_>>(),
            &incidence,
            2,
            128,
        )
        .unwrap();
        partition.refine(&[LeafKey::root("whole")]).unwrap();
        partition
            .refine(&[LeafKey {
                root: "whole".into(),
                level: 1,
                index: [0, 0, 0],
            }])
            .unwrap();
        let field = Motion::parse(&object([
            ("denominator", Json::Number("100".into())),
            (
                "numerators",
                Json::Array(vec![
                    Json::Array(vec![
                        Json::Number("10".into()),
                        Json::Number("15".into()),
                        Json::Number("-7".into()),
                        Json::Number("2".into()),
                    ]),
                    Json::Array(vec![
                        Json::Number("0".into()),
                        Json::Number("0".into()),
                        Json::Number("0".into()),
                        Json::Number("0".into()),
                    ]),
                    Json::Array(vec![
                        Json::Number("0".into()),
                        Json::Number("0".into()),
                        Json::Number("0".into()),
                        Json::Number("0".into()),
                    ]),
                ]),
            ),
        ]))
        .unwrap();
        let (base_space, base_lift) = moving_support_lift(&recipe, &partition, field);
        let motion = SectorMotion(recipe);
        let rates = motion
            .lift_rates_with_motion(&base_space.space, &base_lift, 0, field)
            .unwrap();
        let step = 1. / 510.;
        let plus_recipe = RationalSector {
            half_angle_tangent: [128, 255],
            ..recipe
        };
        let minus_recipe = RationalSector {
            half_angle_tangent: [127, 255],
            ..recipe
        };
        let (_, plus_lift) = moving_support_lift(&plus_recipe, &partition, field);
        let (_, minus_lift) = moving_support_lift(&minus_recipe, &partition, field);
        assert_eq!(base_lift.free_roots, plus_lift.free_roots);
        assert_eq!(base_lift.free_roots, minus_lift.free_roots);
        for (i, rate) in rates.iter().enumerate() {
            let plus = plus_lift.lines[i]
                .terms
                .iter()
                .map(|(j, c)| (*j, c.numerator() as f64 / c.denominator() as f64))
                .collect::<std::collections::BTreeMap<_, _>>();
            let minus = minus_lift.lines[i]
                .terms
                .iter()
                .map(|(j, c)| (*j, c.numerator() as f64 / c.denominator() as f64))
                .collect::<std::collections::BTreeMap<_, _>>();
            for j in 0..base_lift.free_roots.len() {
                let analytic = rate
                    .terms
                    .iter()
                    .find(|(k, _)| *k == j)
                    .map_or(0., |(_, x)| *x);
                let finite = (plus.get(&j).copied().unwrap_or(0.)
                    - minus.get(&j).copied().unwrap_or(0.))
                    / (2. * step);
                assert!(
                    (analytic - finite).abs() < 5e-5,
                    "lift row {i}, root {j}: analytic {analytic}, FD {finite}"
                );
            }
            let finite_offset = (plus_lift.lines[i].offset.numerator() as f64
                / plus_lift.lines[i].offset.denominator() as f64
                - minus_lift.lines[i].offset.numerator() as f64
                    / minus_lift.lines[i].offset.denominator() as f64)
                / (2. * step);
            assert!(
                (rate.offset - finite_offset).abs() < 5e-5,
                "lift row {i} offset: analytic {}, FD {finite_offset}",
                rate.offset
            );
        }
        let mut invalid_lift = base_lift.clone();
        invalid_lift.lines[0].offset = crate::constraint_system::Ratio::new(100, 1).unwrap();
        assert!(motion
            .lift_rates_with_motion(&base_space.space, &invalid_lift, 0, field)
            .is_err());
    }
}
