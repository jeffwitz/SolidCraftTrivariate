//! Local energy norm of the difference between nested Q2 fields.

use solidcraft_trivariate::{
    adaptive_partition::{AdaptivePartition, LeafKey},
    adaptive_space::AdaptiveSpace,
    elasticity::{basis_geometry, gauss},
};
use std::collections::BTreeMap;

pub fn difference_energy_by_parent(
    target: &AdaptiveSpace,
    target_partition: &AdaptivePartition,
    source_partition: &AdaptivePartition,
    delta: &[f64],
    young: f64,
    poisson: f64,
    order: usize,
) -> Result<BTreeMap<LeafKey, f64>, &'static str> {
    if delta.len() != 3 * target.space.ndof
        || delta.iter().any(|value| !value.is_finite())
        || target_partition.leaves().len() != target.space.blocks.len()
        || !young.is_finite()
        || young <= 0.0
        || !(-1.0 < poisson && poisson < 0.5)
    {
        return Err("invalid enrichment energy inputs");
    }
    let quadrature = gauss(order)?;
    let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let mu = young / (2.0 * (1.0 + poisson));
    let leaf_by_block = target
        .leaf_blocks
        .iter()
        .map(|(leaf, block)| (block.as_str(), leaf))
        .collect::<BTreeMap<_, _>>();
    let mut result = source_partition
        .leaves()
        .iter()
        .cloned()
        .map(|leaf| (leaf, 0.0))
        .collect::<BTreeMap<_, _>>();
    for block in &target.space.blocks {
        let child = leaf_by_block
            .get(block.id.as_str())
            .ok_or("missing target leaf-to-block map")?;
        if child.level == 0 {
            return Err("enriched target contains a root leaf");
        }
        let source = LeafKey {
            root: child.root.clone(),
            level: child.level - 1,
            index: child.index.map(|coordinate| coordinate / 2),
        };
        if !source_partition.leaves().contains(&source) {
            return Err("target leaf is not one uniform level below its source leaf");
        }
        let mut energy = 0.0;
        for &(u, wu) in &quadrature {
            for &(v, wv) in &quadrature {
                for &(w, ww) in &quadrature {
                    let basis = basis_geometry(block, [u, v, w])?;
                    if basis.determinant <= 0.0 {
                        return Err("nonpositive enrichment geometry determinant");
                    }
                    let mut gradient = [[0.0; 3]; 3];
                    for (local, &dof) in block.dofs.iter().enumerate() {
                        for component in 0..3 {
                            for (axis, value) in gradient[component].iter_mut().enumerate() {
                                *value += delta[3 * dof + component] * basis.gradients[local][axis];
                            }
                        }
                    }
                    let strain = std::array::from_fn::<_, 3, _>(|i| {
                        std::array::from_fn::<_, 3, _>(|j| 0.5 * (gradient[i][j] + gradient[j][i]))
                    });
                    let trace = strain[0][0] + strain[1][1] + strain[2][2];
                    let norm_squared = strain
                        .iter()
                        .flatten()
                        .map(|value| value * value)
                        .sum::<f64>();
                    let density = lambda * trace * trace + 2.0 * mu * norm_squared;
                    energy += wu * wv * ww * basis.determinant * density;
                }
            }
        }
        if !energy.is_finite() || energy < 0.0 {
            return Err("invalid local enrichment energy");
        }
        *result
            .get_mut(&source)
            .ok_or("missing source leaf energy accumulator")? += energy;
    }
    Ok(result)
}

pub fn bulk_mark(
    contribution_by_leaf: &BTreeMap<LeafKey, f64>,
    theta: f64,
) -> Result<Vec<LeafKey>, &'static str> {
    if !theta.is_finite()
        || !(0.0 < theta && theta <= 1.0)
        || contribution_by_leaf
            .values()
            .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err("invalid enrichment Dörfler contributions");
    }
    let total = contribution_by_leaf.values().sum::<f64>();
    if total <= 0.0 {
        return Ok(Vec::new());
    }
    let mut ordered = contribution_by_leaf
        .iter()
        .map(|(leaf, value)| (leaf.clone(), *value))
        .collect::<Vec<_>>();
    ordered.sort_by(|left, right| {
        right
            .1
            .total_cmp(&left.1)
            .then_with(|| left.0.cmp(&right.0))
    });
    let mut selected = Vec::new();
    let mut accumulated = 0.0;
    for (leaf, value) in ordered {
        selected.push(leaf);
        accumulated += value;
        if accumulated >= theta * total {
            break;
        }
    }
    selected.sort();
    Ok(selected)
}

pub fn galerkin_identity_relative_error(
    difference_energy: f64,
    current_energy: f64,
    enriched_energy: f64,
) -> Result<f64, &'static str> {
    let expected = 2.0 * (enriched_energy - current_energy);
    if [difference_energy, current_energy, enriched_energy]
        .iter()
        .any(|value| !value.is_finite())
        || difference_energy < 0.0
        || expected <= 0.0
    {
        return Err("invalid nested Galerkin energy identity data");
    }
    Ok((difference_energy - expected).abs() / expected.abs())
}
