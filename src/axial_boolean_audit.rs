//! Re-audit the geometry, topology and boundary ownership of a stored axial boolean.
use crate::{
    axial_boolean::{
        AxialAxis, AxialBooleanModel, AxialBooleanOperation as Op, AxialOperand,
        AxialSectionFamily as Family, OperandLabel,
    },
    bracket::BracketBlock,
    incidence::{enumerate_faces, BlockSide},
    section_audit::audit_section,
};
use std::collections::{BTreeMap, BTreeSet};

fn occupies(op: Op, a: bool, b: bool) -> bool {
    match op {
        Op::Union => a || b,
        Op::Intersection => a && b,
        Op::Difference => a && !b,
    }
}

fn memberships(z: i32, a: AxialOperand, b: AxialOperand) -> Vec<OperandLabel> {
    [(a, OperandLabel::A), (b, OperandLabel::B)]
        .into_iter()
        .filter(|(o, _)| (o.interval_mm[0]..o.interval_mm[1]).contains(&z))
        .map(|(_, label)| label)
        .collect()
}

/// Compare a candidate with an independently classified interval partition and the G2
/// section certificate. All coordinates of all output cells participate in the check.
pub fn audit_axial_boolean(
    op: Op,
    a: AxialOperand,
    b: AxialOperand,
    model: &AxialBooleanModel,
) -> Result<(), &'static str> {
    let (family, area) = match a.section {
        Family::ReferenceRectangleV1 => ("reference_rectangle_v1", 400.),
        Family::PiercedSquareV1 => ("pierced_square_v1", 400. - 16. * std::f64::consts::PI),
        _ => return Err("axial_audit_unsupported_section"),
    };
    if a.section != b.section
        || [a.axis, b.axis, model.axis] != [AxialAxis::Z; 3]
        || a.section != model.section
        || op != model.operation
    {
        return Err("axial_audit_recipe");
    }
    for [lo, hi] in [a.interval_mm, b.interval_mm] {
        if lo >= hi || lo < -1000 || hi > 1000 {
            return Err("axial_audit_interval");
        }
    }
    let cuts: Vec<_> = a
        .interval_mm
        .into_iter()
        .chain(b.interval_mm)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut expected = BTreeMap::new();
    let mut connected = 0;
    let mut prior = false;
    let mut length = 0;
    for (index, pair) in cuts.windows(2).enumerate() {
        let owners = memberships(pair[0], a, b);
        let keep = occupies(
            op,
            owners.contains(&OperandLabel::A),
            owners.contains(&OperandLabel::B),
        );
        if keep {
            connected += usize::from(!prior);
            length += pair[1] - pair[0];
            expected.insert([pair[0], pair[1]], (index, owners));
        }
        prior = keep;
    }
    if model.component_count != connected
        || !model.volume_mm3.is_finite()
        || (model.volume_mm3 - f64::from(length) * area).abs()
            > 1e-12 * (f64::from(length) * area).max(1.)
        || model.selected_intervals.len() != expected.len()
    {
        return Err("axial_audit_measure");
    }
    for (interval, (z, (_, owners))) in model.selected_intervals.iter().zip(&expected) {
        if interval.bounds_mm != *z || &interval.operands != owners {
            return Err("axial_audit_membership");
        }
    }
    let mut layers = BTreeMap::<[i32; 2], Vec<BracketBlock>>::new();
    let mut block_info = BTreeMap::new();
    for block in &model.blocks {
        if !matches!(block.degrees, [1, 1, 1] | [2, 1, 1])
            || block.controls.len() != (block.degrees[0] + 1) * 4
        {
            return Err("axial_audit_layout");
        }
        let mut bounds = [0; 2];
        for (k, value) in bounds.iter_mut().enumerate() {
            let [_, _, z, w] = block.controls[k];
            if w <= 0 || z % w != 0 {
                return Err("axial_audit_noninteger_cap");
            }
            *value = i32::try_from(z / w).map_err(|_| "axial_audit_cap_range")?;
        }
        let (index, _) = expected
            .get(&bounds)
            .ok_or("axial_audit_unselected_layer")?;
        let source_id = block
            .id
            .strip_prefix(&format!("axial_boolean_cell_{index:02}_"))
            .ok_or("axial_audit_block_identity")?;
        if block_info
            .insert(block.id.clone(), (source_id.to_owned(), bounds))
            .is_some()
        {
            return Err("axial_audit_duplicate_block");
        }
        let mut restored = block.clone();
        restored.id = source_id.to_owned();
        for (j, p) in restored.controls.iter_mut().enumerate() {
            if p[3] <= 0 || i128::from(p[2]) != i128::from(bounds[j % 2]) * i128::from(p[3]) {
                return Err("axial_audit_nonplanar_cap");
            }
            p[2] = (if j % 2 == 1 { 8i64 } else { 0 })
                .checked_mul(p[3])
                .ok_or("axial_audit_weight_range")?;
        }
        layers.entry(bounds).or_default().push(restored);
    }
    if layers.len() != expected.len() {
        return Err("axial_audit_missing_layer");
    }
    let mut source_labels = BTreeMap::new();
    for (z, blocks) in layers {
        let certificate = audit_section(family, &blocks)?;
        for (side, feature) in certificate.design_faces {
            source_labels.insert((z, side), feature);
        }
    }
    let incidence = enumerate_faces(&model.blocks)?;
    if incidence != model.incidence || model.faces.len() != incidence.exterior.len() {
        return Err("axial_audit_incidence");
    }
    let exterior: BTreeSet<_> = incidence.exterior.iter().collect();
    let mut seen = BTreeSet::new();
    for face in &model.faces {
        let side = &face.side;
        if !exterior.contains(side) || !seen.insert(side) {
            return Err("axial_audit_labelled_side");
        }
        let (source, z) = block_info
            .get(&side.block)
            .ok_or("axial_audit_unknown_block")?;
        let mut truth = BTreeSet::new();
        if side.axis == 2 {
            let cut = z[usize::from(side.upper)];
            for (operand, label) in [(a, OperandLabel::A), (b, OperandLabel::B)] {
                let reverse = label == OperandLabel::B && op == Op::Difference;
                let source_upper = side.upper ^ reverse;
                if cut == operand.interval_mm[usize::from(source_upper)] {
                    truth.insert((label, if source_upper { "z8" } else { "z0" }, reverse));
                }
            }
        } else {
            let original = BlockSide {
                block: source.clone(),
                axis: side.axis,
                upper: side.upper,
            };
            let feature = source_labels
                .get(&(*z, original))
                .ok_or("axial_audit_source_side")?;
            for owner in memberships(z[0], a, b) {
                truth.insert((owner, feature, false));
            }
        }
        let observed: BTreeSet<_> = face
            .provenance
            .iter()
            .map(|p| (p.operand, p.source_face.as_str(), p.orientation_reversed))
            .collect();
        if observed != truth || observed.len() != face.provenance.len() || truth.is_empty() {
            return Err("axial_audit_provenance");
        }
    }
    Ok(())
}
