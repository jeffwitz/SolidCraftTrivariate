//! Compose two native axial booleans while retaining the first result as provenance.
use crate::axial_boolean::{
    build_axial_boolean, AxialAxis, AxialBooleanModel, AxialBooleanOperation as Op, AxialOperand,
    AxialSectionFamily as Family, OperandLabel,
};
use crate::axial_boolean_audit::audit_axial_boolean;
use crate::bracket::BracketBlock;

#[derive(Debug)]
pub struct ComposedResult {
    pub first: AxialBooleanModel,
    pub second: AxialBooleanModel,
    pub composite_operand: AxialOperand,
    /// Maps copied from the first model and exactly restricted to surviving z ranges.
    pub carried_maps: Vec<CarriedMap>,
    /// Original recipe operands represented on each final exterior face.
    pub lineage: Vec<(String, Vec<&'static str>)>,
}

#[derive(Debug)]
pub struct CarriedMap {
    pub block: BracketBlock,
    pub source_block_id: String,
    pub source_interval_mm: [i32; 2],
    pub result_interval_mm: [i32; 2],
    pub original_operands: Vec<&'static str>,
}

fn block_z_bounds(model: &AxialBooleanModel, id: &str) -> Option<[i32; 2]> {
    let block = model.blocks.iter().find(|block| block.id == id)?;
    let mut values = block.controls.iter().map(|p| p[2] / p[3]);
    let first = values.next()?;
    let (lo, hi) = values.fold((first, first), |(lo, hi), z| (lo.min(z), hi.max(z)));
    Some([i32::try_from(lo).ok()?, i32::try_from(hi).ok()?])
}

fn restrict_block(
    block: &BracketBlock,
    source: [i32; 2],
    target: [i32; 2],
    id: String,
) -> Result<BracketBlock, &'static str> {
    if source[0] > target[0] || target[1] > source[1] || target[0] >= target[1] {
        return Err("restriction lies outside source block");
    }
    let mut output = block.clone();
    output.id = id;
    for point in &mut output.controls {
        let weight = point[3];
        let z = point[2];
        point[2] = if z == i64::from(source[0]) * weight {
            i64::from(target[0]) * weight
        } else if z == i64::from(source[1]) * weight {
            i64::from(target[1]) * weight
        } else {
            return Err("source map is not an affine axial product");
        };
    }
    Ok(output)
}

fn carry_first_maps(
    first: &AxialBooleanModel,
    second: &AxialBooleanModel,
) -> Result<Vec<CarriedMap>, &'static str> {
    let mut carried = Vec::new();
    for source_cell in &first.selected_intervals {
        for result_cell in &second.selected_intervals {
            let range = [
                source_cell.bounds_mm[0].max(result_cell.bounds_mm[0]),
                source_cell.bounds_mm[1].min(result_cell.bounds_mm[1]),
            ];
            if range[0] >= range[1] {
                continue;
            }
            for block in &first.blocks {
                if block_z_bounds(first, &block.id) != Some(source_cell.bounds_mm) {
                    continue;
                }
                let restricted = restrict_block(
                    block,
                    source_cell.bounds_mm,
                    range,
                    format!("{}@z{}_{}", block.id, range[0], range[1]),
                )?;
                let original_operands = source_cell
                    .operands
                    .iter()
                    .map(|label| match label {
                        OperandLabel::A => "A",
                        OperandLabel::B => "B",
                    })
                    .collect();
                carried.push(CarriedMap {
                    block: restricted,
                    source_block_id: block.id.clone(),
                    source_interval_mm: source_cell.bounds_mm,
                    result_interval_mm: range,
                    original_operands,
                });
            }
        }
    }
    if carried.is_empty() {
        return Err("no first-stage maps survived the second operation");
    }
    Ok(carried)
}

fn assert_matches_cold_restrictions(
    carried: &[CarriedMap],
    cold: &AxialBooleanModel,
) -> Result<(), &'static str> {
    for reused in carried {
        let cold_source = cold
            .blocks
            .iter()
            .find(|block| {
                block_z_bounds(cold, &block.id).is_some_and(|range| {
                    range[0] <= reused.result_interval_mm[0]
                        && reused.result_interval_mm[1] <= range[1]
                })
            })
            .ok_or("cold builder has no matching result cell")?;
        let reference = restrict_block(
            cold_source,
            block_z_bounds(cold, &cold_source.id).ok_or("bad cold block bounds")?,
            reused.result_interval_mm,
            "cold-restricted-reference".into(),
        )?;
        if reused.block.degrees != reference.degrees
            || reused.block.controls.len() != reference.controls.len()
            || reused
                .block
                .controls
                .iter()
                .zip(&reference.controls)
                .any(|(actual, expected)| actual != expected)
        {
            return Err("reused map differs from the cold-builder restriction");
        }
    }
    Ok(())
}

fn original_owners(
    first: &AxialBooleanModel,
    second: &AxialBooleanModel,
    carried: &[CarriedMap],
    face_index: usize,
) -> Result<Vec<&'static str>, &'static str> {
    let face = second.faces.get(face_index).ok_or("missing final face")?;
    let block_range = block_z_bounds(second, &face.side.block).ok_or("bad final block bounds")?;
    let mut result = Vec::new();
    for provenance in &face.provenance {
        if provenance.operand == OperandLabel::B {
            result.push("C");
            continue;
        }
        for map in carried {
            for source in first.faces.iter().filter(|source| {
                source.side.block == map.source_block_id
                    && source.side.axis == face.side.axis
                    && source.side.upper == face.side.upper
            }) {
                let compatible = if face.side.axis == 2 {
                    let z = if face.side.upper {
                        block_range[1]
                    } else {
                        block_range[0]
                    };
                    (map.result_interval_mm[0] == z || map.result_interval_mm[1] == z)
                        && (map.source_interval_mm[0] == z || map.source_interval_mm[1] == z)
                } else {
                    map.result_interval_mm[0] < block_range[1]
                        && block_range[0] < map.result_interval_mm[1]
                };
                if compatible {
                    result.extend(
                        source
                            .provenance
                            .iter()
                            .filter(|p| p.source_face == provenance.source_face)
                            .map(|p| match p.operand {
                                OperandLabel::A => "A",
                                OperandLabel::B => "B",
                            }),
                    );
                }
            }
        }
    }
    result.sort_unstable();
    result.dedup();
    if result.is_empty() {
        return Err("final face has no original provenance");
    }
    Ok(result)
}

/// Construct `(A union B) minus C`; the second stage consumes the first result's
/// certified occupied interval, rather than rebuilding it from the recipe inputs.
pub fn compose(
    a_interval: [i32; 2],
    b_interval: [i32; 2],
    c_interval: [i32; 2],
) -> Result<ComposedResult, &'static str> {
    let operand = |interval_mm| AxialOperand {
        section: Family::ReferenceRectangleV1,
        axis: AxialAxis::Z,
        interval_mm,
    };
    let a = operand(a_interval);
    let b = operand(b_interval);
    let c = operand(c_interval);
    let first = build_axial_boolean(Op::Union, a, b).map_err(|_| "first-stage union failed")?;
    audit_axial_boolean(Op::Union, a, b, &first).map_err(|_| "first-stage audit failed")?;

    let first_intervals = &first.selected_intervals;
    let bounds = [
        first_intervals
            .first()
            .ok_or("empty first-stage result")?
            .bounds_mm[0],
        first_intervals
            .last()
            .ok_or("empty first-stage result")?
            .bounds_mm[1],
    ];
    if first_intervals
        .windows(2)
        .any(|pair| pair[0].bounds_mm[1] != pair[1].bounds_mm[0])
    {
        return Err(
            "first-stage result is disconnected; a single axial operand cannot represent it",
        );
    }
    let composite_operand = operand(bounds);
    let second = build_axial_boolean(Op::Difference, composite_operand, c)
        .map_err(|_| "second-stage difference failed")?;
    audit_axial_boolean(Op::Difference, composite_operand, c, &second)
        .map_err(|_| "second-stage audit failed")?;
    let carried_maps = carry_first_maps(&first, &second)?;
    assert_matches_cold_restrictions(&carried_maps, &second)?;

    let lineage = (0..second.faces.len())
        .map(|index| {
            let owners = original_owners(&first, &second, &carried_maps, index)?;
            Ok((format!("face:{index}"), owners))
        })
        .collect::<Result<_, &'static str>>()?;
    Ok(ComposedResult {
        first,
        second,
        composite_operand,
        carried_maps,
        lineage,
    })
}

/// Evaluate a two-stage axial recipe and serialize surviving source restrictions
/// and original face owners. Coordinates are integer millimetres along z.
pub fn build_document(recipe: &crate::json::Json) -> Result<crate::json::Json, &'static str> {
    use crate::json::{object, Json};
    fn nums<T: ToString>(values: impl IntoIterator<Item = T>) -> Json {
        Json::Array(
            values
                .into_iter()
                .map(|v| Json::Number(v.to_string()))
                .collect(),
        )
    }
    fn labels(values: &[&str]) -> Json {
        Json::Array(values.iter().map(|s| Json::String((*s).into())).collect())
    }
    recipe.keys(&["family", "a_mm", "b_mm", "c_mm"])?;
    if recipe.field("family")?.text()? != "composed_axial_job_v1" {
        return Err("expected composed_axial_job_v1");
    }
    let interval = |key| -> Result<[i32; 2], &'static str> {
        let [lo, hi] = recipe.field(key)?.items()? else {
            return Err("expected two axial endpoints");
        };
        Ok([lo.number()?, hi.number()?])
    };
    let result = compose(interval("a_mm")?, interval("b_mm")?, interval("c_mm")?)?;
    let carried = result.carried_maps.iter().map(|m| {
        object([
            ("id", Json::String(m.block.id.clone())),
            ("degrees", nums(m.block.degrees)),
            (
                "controls",
                Json::Array(m.block.controls.iter().map(nums).collect()),
            ),
            ("source_block_id", Json::String(m.source_block_id.clone())),
            ("source_interval_mm", nums(m.source_interval_mm)),
            ("result_interval_mm", nums(m.result_interval_mm)),
            ("original_operands", labels(&m.original_operands)),
        ])
    });
    Ok(object([
        ("schema", Json::String("composed_axial_result_v1".into())),
        ("recipe", recipe.clone()),
        (
            "volume_mm3",
            Json::Number(result.second.volume_mm3.to_string()),
        ),
        ("maps", Json::Array(carried.collect())),
        (
            "face_lineage",
            Json::Array(
                result
                    .lineage
                    .iter()
                    .map(|(id, owners)| {
                        object([
                            ("id", Json::String(id.clone())),
                            ("original_operands", labels(owners)),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]))
}
