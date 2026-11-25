//! Canonical interchange for certified axial-product operations on the fixed z axis.
use crate::axial_boolean::{
    build_axial_boolean, AxialAxis, AxialBooleanOperation as Op, AxialOperand,
    AxialSectionFamily as Family, OperandLabel,
};
use crate::axial_boolean_audit::audit_axial_boolean;
use crate::bracket_json::{encode_blocks, side_json, topology};
use crate::json::{object, Json};

fn string(s: &str) -> Json {
    Json::String(s.to_owned())
}
fn number(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn label(l: OperandLabel) -> Json {
    string(match l {
        OperandLabel::A => "A",
        OperandLabel::B => "B",
    })
}

fn canonical(recipe: &Json) -> Result<(Json, Json), &'static str> {
    recipe.keys(&["family", "section", "operation", "a_z_mm", "b_z_mm"])?;
    if recipe.field("family")?.text()? != "axial_boolean_v1" {
        return Err("axial_boolean_recipe");
    }
    let section = match recipe.field("section")?.text()? {
        "reference_rectangle_v1" => Family::ReferenceRectangleV1,
        "pierced_square_v1" => Family::PiercedSquareV1,
        _ => return Err("axial_boolean_unsupported_section"),
    };
    let op = match recipe.field("operation")?.text()? {
        "union" => Op::Union,
        "intersection" => Op::Intersection,
        "difference" => Op::Difference,
        _ => return Err("axial_boolean_operation"),
    };
    let operand = |field| -> Result<AxialOperand, &'static str> {
        let [a, b] = recipe.field(field)?.items()? else {
            return Err("axial_boolean_interval");
        };
        Ok(AxialOperand {
            section,
            axis: AxialAxis::Z,
            interval_mm: [a.number()?, b.number()?],
        })
    };
    let (a, b) = (operand("a_z_mm")?, operand("b_z_mm")?);
    let m = build_axial_boolean(op, a, b).map_err(|e| match e {
        crate::axial_boolean::AxialBooleanError::NullInterval => "axial_boolean_zero_height",
        crate::axial_boolean::AxialBooleanError::InvertedInterval => {
            "axial_boolean_inverted_interval"
        }
        crate::axial_boolean::AxialBooleanError::IntervalOutOfRange => {
            "axial_boolean_interval_out_of_range"
        }
        _ => "axial_boolean_failed",
    })?;
    audit_axial_boolean(op, a, b, &m)?;
    let report = object([
        ("status", string("valid")),
        ("method", string("certified_axial_product_partition")),
        ("components", number(m.component_count)),
        ("blocks", number(m.blocks.len())),
        ("shared_faces", number(m.incidence.interfaces.len())),
        ("exterior_sides", number(m.incidence.exterior.len())),
        ("ideal_volume_mm3", number(m.volume_mm3)),
    ]);
    let intervals = m
        .selected_intervals
        .iter()
        .map(|s| {
            object([
                ("z_mm", Json::Array(s.bounds_mm.map(number).into())),
                (
                    "owners",
                    Json::Array(s.operands.iter().copied().map(label).collect()),
                ),
            ])
        })
        .collect();
    let faces = m
        .faces
        .iter()
        .map(|f| {
            object([
                ("side", side_json(&f.side)),
                (
                    "sources",
                    Json::Array(
                        f.provenance
                            .iter()
                            .map(|s| {
                                object([
                                    ("operand", label(s.operand)),
                                    ("feature", string(&s.source_face)),
                                    ("orientation_reversed", Json::Bool(s.orientation_reversed)),
                                ])
                            })
                            .collect(),
                    ),
                ),
            ])
        })
        .collect();
    let model = object([
        ("schema", string("axial_boolean_v1")),
        ("recipe", recipe.clone()),
        ("intervals", Json::Array(intervals)),
        ("topology", topology(&m.incidence)),
        ("design_faces", Json::Array(faces)),
        ("blocks", encode_blocks(m.blocks)),
    ]);
    Ok((model, report))
}

pub(super) fn build(recipe: &Json) -> Result<Json, &'static str> {
    Ok(canonical(recipe)?.0)
}

pub(super) fn audit(model: &Json) -> Result<Json, &'static str> {
    let (expected, report) = canonical(model.field("recipe")?)?;
    if model != &expected {
        return Err("axial_boolean_canonical_model_mismatch");
    }
    Ok(report)
}
