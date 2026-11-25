//! Exact unisolvent audit for a stored cylindrical boundary face.

use crate::{bracket::BracketBlock, constraint_system::Ratio, incidence::BlockSide};
use std::collections::BTreeMap;

/// Summary of exact cylindrical boundary identities checked in an audited model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundaryAudit {
    pub total_faces: usize,
    pub faces_by_bore: BTreeMap<String, usize>,
}

fn add(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_add(b)
}

fn sub(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_sub(b)
}

fn mul(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_mul(b)
}

fn basis(degree: usize, parameter: Ratio) -> Result<Vec<Ratio>, &'static str> {
    let one_minus = sub(Ratio::one(), parameter)?;
    match degree {
        1 => Ok(vec![one_minus, parameter]),
        2 => Ok(vec![
            mul(one_minus, one_minus)?,
            mul(Ratio::new(2, 1)?, mul(parameter, one_minus)?)?,
            mul(parameter, parameter)?,
        ]),
        _ => Err("unsupported Steinmetz face degree"),
    }
}

fn face_axes(axis: usize) -> [usize; 2] {
    match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}

fn flat(index: [usize; 3], degrees: [usize; 3]) -> usize {
    (index[0] * (degrees[1] + 1) + index[1]) * (degrees[2] + 1) + index[2]
}

fn evaluate_homogeneous_face(
    block: &BracketBlock,
    side: &BlockSide,
    a: Ratio,
    b: Ratio,
) -> Result<[Ratio; 4], &'static str> {
    let free = face_axes(side.axis);
    let basis_a = basis(block.degrees[free[0]], a)?;
    let basis_b = basis(block.degrees[free[1]], b)?;
    let fixed = if side.upper {
        block.degrees[side.axis]
    } else {
        0
    };
    let mut value = [Ratio::zero(); 4];
    for (i, &wa) in basis_a.iter().enumerate() {
        for (j, &wb) in basis_b.iter().enumerate() {
            let mut index = [0; 3];
            index[side.axis] = fixed;
            index[free[0]] = i;
            index[free[1]] = j;
            let weight = mul(wa, wb)?;
            let control = block.controls[flat(index, block.degrees)];
            for component in 0..4 {
                value[component] = add(
                    value[component],
                    mul(weight, Ratio::new(i128::from(control[component]), 1)?)?,
                )?;
            }
        }
    }
    Ok(value)
}

/// Return true exactly when the complete stored face lies on the requested
/// physical circular cylinder. The identity is tested at a tensor unisolvent
/// grid for degree `(2p, 2q)`, where `(p,q)` are the two face degrees.
///
/// `axis` is the physical cylinder axis. The cylinder center uses all three
/// physical coordinates; its two components normal to `axis` define the circle.
/// `side` selects the parametric block face that is being audited.
pub fn certify_boundary(
    block: &BracketBlock,
    side: &BlockSide,
    axis: usize,
    center: [i32; 3],
    radius: i32,
) -> Result<bool, &'static str> {
    if axis > 2
        || side.axis > 2
        || side.block != block.id
        || radius <= 0
        || block.degrees.iter().any(|&d| !(1..=2).contains(&d))
        || block.controls.len() != block.degrees.iter().map(|&d| d + 1).product::<usize>()
        || block.controls.iter().any(|p| p[3] <= 0)
    {
        return Err("invalid Steinmetz boundary input");
    }
    let free = face_axes(side.axis);
    let degrees = [block.degrees[free[0]], block.degrees[free[1]]];
    let samples = [2 * degrees[0] + 1, 2 * degrees[1] + 1];
    let transverse: Vec<_> = (0..3).filter(|&coordinate| coordinate != axis).collect();
    let radius_squared = i128::from(radius)
        .checked_mul(i128::from(radius))
        .ok_or("Steinmetz radius overflow")?;
    for i in 0..samples[0] {
        let a = Ratio::new(i as i128, (samples[0] - 1) as i128)?;
        for j in 0..samples[1] {
            let b = Ratio::new(j as i128, (samples[1] - 1) as i128)?;
            let homogeneous = evaluate_homogeneous_face(block, side, a, b)?;
            if homogeneous[3].numerator() <= 0 {
                return Err("nonpositive Steinmetz face weight");
            }
            let mut radius_identity = Ratio::zero();
            for &coordinate in &transverse {
                let centered = sub(
                    homogeneous[coordinate],
                    mul(
                        Ratio::new(i128::from(center[coordinate]), 1)?,
                        homogeneous[3],
                    )?,
                )?;
                radius_identity = add(radius_identity, mul(centered, centered)?)?;
            }
            let target = mul(
                Ratio::new(radius_squared, 1)?,
                mul(homogeneous[3], homogeneous[3])?,
            )?;
            if radius_identity != target {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Re-audit the model, then verify every selected `bore/<id>/wall` face by
/// exact unisolvent evaluation of its stored homogeneous trace.
pub fn certify_wall_features(model: &crate::json::Json) -> Result<BoundaryAudit, &'static str> {
    let schema = model.field("schema")?.text()?;
    if schema != "orthogonal_sketch_features_v1" {
        return Err("Steinmetz wall audit requires an orthogonal sketch model");
    }
    let audit = crate::interchange::audit(model)?;
    if audit.field("status")?.text()? != "valid" {
        return Err("Steinmetz source model is not certified valid");
    }
    let blocks = crate::bracket_json::decode_blocks(model)?;
    let block_map: BTreeMap<_, _> = blocks
        .iter()
        .map(|block| (block.id.as_str(), block))
        .collect();
    let mut bore_map = BTreeMap::<String, (usize, [i32; 3], i32)>::new();
    for bore in model.field("bores")?.items()? {
        let id = bore.field("id")?.text()?.to_owned();
        let axis = match bore.field("axis")?.text()? {
            "x" => 0,
            "y" => 1,
            "z" => 2,
            _ => return Err("unsupported Steinmetz cylinder axis"),
        };
        let center: [i32; 3] = bore
            .field("center_mm")?
            .items()?
            .iter()
            .map(|value| value.number())
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .map_err(|_| "Steinmetz cylinder center must have three coordinates")?;
        let radius: i32 = bore.field("radius_mm")?.number()?;
        if radius <= 0 || bore_map.insert(id, (axis, center, radius)).is_some() {
            return Err("invalid or duplicate Steinmetz cylinder");
        }
    }
    let mut faces_by_bore = BTreeMap::<String, usize>::new();
    for face in model.field("design_faces")?.items()? {
        let feature = face.field("feature")?.text()?;
        let Some(rest) = feature.strip_prefix("bore/") else {
            continue;
        };
        let Some(id) = rest.strip_suffix("/wall") else {
            continue;
        };
        let (axis, center, radius) = *bore_map.get(id).ok_or("wall selection has no cylinder")?;
        let side_name = face.field("side")?.text()?;
        let (side_axis, upper) = match side_name {
            "u0" => (0, false),
            "u1" => (0, true),
            "v0" => (1, false),
            "v1" => (1, true),
            "w0" => (2, false),
            "w1" => (2, true),
            _ => return Err("unknown Steinmetz design face side"),
        };
        let block_id = face.field("block")?.text()?;
        let block = block_map
            .get(block_id)
            .ok_or("wall face references unknown block")?;
        let side = BlockSide {
            block: block_id.to_owned(),
            axis: side_axis,
            upper,
        };
        if !certify_boundary(block, &side, axis, center, radius)? {
            return Err("stored bore wall is not on its declared circular cylinder");
        }
        *faces_by_bore.entry(id.to_owned()).or_default() += 1;
    }
    if bore_map
        .keys()
        .any(|id| faces_by_bore.get(id).copied().unwrap_or(0) == 0)
    {
        return Err("a Steinmetz cylinder has no selected wall face");
    }
    let total_faces = faces_by_bore.values().sum();
    Ok(BoundaryAudit {
        total_faces,
        faces_by_bore,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quarter_cylinder() -> BracketBlock {
        // u follows the rational quarter circle, v is the audited face axis,
        // and w is the cylinder axis. The other radial layer has radius 3.
        let mut controls = Vec::new();
        for angle_control in 0..3 {
            let weight = [4, 4, 8][angle_control];
            for radius in [3, 2] {
                let (x, z) = match angle_control {
                    0 => (4 * radius, 0),
                    1 => (4 * radius, 4 * radius),
                    _ => (0, 8 * radius),
                };
                for axial in [0, 3] {
                    controls.push([x, axial * weight, z, weight]);
                }
            }
        }
        BracketBlock {
            id: "quarter".into(),
            degrees: [2, 1, 1],
            controls,
        }
    }

    #[test]
    fn rational_unisolvent_grid_certifies_exact_circle_and_detects_mutation() {
        let block = quarter_cylinder();
        let side = BlockSide {
            block: block.id.clone(),
            axis: 1,
            upper: true,
        };
        assert!(certify_boundary(&block, &side, 1, [0, 0, 0], 2).unwrap());

        let mut changed = block.clone();
        changed.controls[2][0] += 1;
        assert!(!certify_boundary(&changed, &side, 1, [0, 0, 0], 2).unwrap());
    }

    #[test]
    fn malformed_axes_side_identity_and_weights_are_rejected() {
        let block = quarter_cylinder();
        let side = BlockSide {
            block: block.id.clone(),
            axis: 1,
            upper: true,
        };
        assert!(certify_boundary(&block, &side, 3, [0, 0, 0], 2).is_err());
        let wrong_side = BlockSide {
            block: "other".into(),
            ..side.clone()
        };
        assert!(certify_boundary(&block, &wrong_side, 1, [0, 0, 0], 2).is_err());
        let mut bad = block;
        bad.controls[0][3] = 0;
        assert!(certify_boundary(&bad, &side, 1, [0, 0, 0], 2).is_err());
    }
}
