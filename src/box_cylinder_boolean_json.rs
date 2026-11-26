//! Versioned JSON exchange for the two-operand box-minus-cylinder constructor.
use crate::{
    box_cylinder_boolean::{subtract_through_cylinder, BoxOperand, ThroughCylinderOperand},
    bracket_json::{design, encode_blocks, topology},
    json::{object, Json},
};

pub(crate) const FAMILY: &str = "box_cylinder_boolean_v1";
pub(crate) const METHOD: &str = "certified_box_cylinder_boolean";

fn number(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}
fn string(value: &str) -> Json {
    Json::String(value.to_owned())
}
fn int_array<const N: usize>(value: &Json) -> Result<[i32; N], &'static str> {
    value
        .items()?
        .iter()
        .map(Json::number::<i32>)
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| "incorrect operand coordinate count")
}

#[derive(Clone, Copy)]
struct Inputs {
    stock: BoxOperand,
    tool: ThroughCylinderOperand,
}

fn inputs(value: &Json, is_model: bool) -> Result<Inputs, &'static str> {
    value.keys(if is_model {
        &[
            "schema",
            "operation",
            "stock",
            "tool",
            "blocks",
            "topology",
            "design_faces",
            "partition_plan",
            "volume_mm3",
        ]
    } else {
        &["family", "operation", "stock", "tool"]
    })?;
    if value
        .field(if is_model { "schema" } else { "family" })?
        .text()?
        != FAMILY
    {
        return Err("unexpected_box_cylinder_schema");
    }
    let stock = value.field("stock")?;
    stock.keys(&["kind", "min_mm", "max_mm"])?;
    if stock.field("kind")?.text()? != "box" {
        return Err("stock_operand_must_be_box");
    }
    let tool = value.field("tool")?;
    tool.keys(&["kind", "axis", "center_mm", "radius_mm", "interval_mm"])?;
    if tool.field("kind")?.text()? != "cylinder" {
        return Err("tool_operand_must_be_cylinder");
    }
    if tool.field("axis")?.text()? != "z" {
        return Err("cylinder_axis_must_be_z");
    }
    if value.field("operation")?.text()? != "difference" {
        return Err("only_box_minus_cylinder_is_supported");
    }
    let [cx, cy] = int_array(tool.field("center_mm")?)?;
    Ok(Inputs {
        stock: BoxOperand {
            min_mm: int_array(stock.field("min_mm")?)?,
            max_mm: int_array(stock.field("max_mm")?)?,
        },
        tool: ThroughCylinderOperand {
            center_mm: [cx, cy],
            radius_mm: tool.field("radius_mm")?.number()?,
            interval_mm: int_array(tool.field("interval_mm")?)?,
        },
    })
}

fn build_result(input: Inputs) -> Result<Json, &'static str> {
    let built = subtract_through_cylinder(input.stock, input.tool)?;
    let mut out = object([
        ("schema", string(FAMILY)),
        ("operation", string("difference")),
        (
            "stock",
            object([
                ("kind", string("box")),
                (
                    "min_mm",
                    Json::Array(input.stock.min_mm.into_iter().map(number).collect()),
                ),
                (
                    "max_mm",
                    Json::Array(input.stock.max_mm.into_iter().map(number).collect()),
                ),
            ]),
        ),
        (
            "tool",
            object([
                ("kind", string("cylinder")),
                ("axis", string("z")),
                (
                    "center_mm",
                    Json::Array(input.tool.center_mm.into_iter().map(number).collect()),
                ),
                ("radius_mm", number(input.tool.radius_mm)),
                (
                    "interval_mm",
                    Json::Array(input.tool.interval_mm.into_iter().map(number).collect()),
                ),
            ]),
        ),
        ("blocks", encode_blocks(built.blocks)),
        ("topology", topology(&built.certificate.incidence)),
        ("design_faces", design(&built.certificate)),
        (
            "partition_plan",
            object([
                ("method", string(built.planner_method)),
                ("boolean", string("stock_minus_through_cylinder")),
                ("tool_sector_count", number(4)),
                (
                    "determinant_lower_bound_mm3",
                    number(built.certificate.determinant_lower_bound_mm3),
                ),
            ]),
        ),
    ]);
    if let Json::Object(fields) = &mut out {
        fields.insert(
            "volume_mm3".into(),
            object([
                ("constant", number(built.volume_constant_mm3)),
                ("pi_coefficient", number(built.volume_pi_coefficient_mm3)),
            ]),
        );
    }
    Ok(out)
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    build_result(inputs(recipe, false)?)
}

pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    let input = inputs(model, true)?;
    let expected = build_result(input)?;
    for key in [
        "operation",
        "stock",
        "tool",
        "blocks",
        "topology",
        "design_faces",
        "partition_plan",
        "volume_mm3",
    ] {
        if model.field(key)? != expected.field(key)? {
            return Err("box-cylinder model differs from exact operand rebuild");
        }
    }
    let interfaces = model.field("topology")?.field("interfaces")?.items()?.len();
    Ok(object([
        ("status", string("valid")),
        ("method", string(METHOD)),
        ("blocks", number(model.field("blocks")?.items()?.len())),
        ("interfaces", number(interfaces)),
        (
            "determinant_lower_bound_mm3",
            expected
                .field("partition_plan")?
                .field("determinant_lower_bound_mm3")?
                .clone(),
        ),
        ("volume_mm3", expected.field("volume_mm3")?.clone()),
    ]))
}
