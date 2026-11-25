//! Recipe interface to native sketch primitives and their construction-time plans.
use crate::{
    constraint_system::Ratio,
    json::{object, Json},
    primitive_template::{Operation, Primitive, PrimitivePlan},
    sketch_junction::Continuity,
};
type Error = &'static str;
fn number(v: impl ToString) -> Json {
    Json::Number(v.to_string())
}
fn coordinate(v: &Json) -> Result<Ratio, Error> {
    let [n, d] = v.items()? else {
        return Err("expected exact coordinate [numerator,denominator]");
    };
    Ratio::new(n.number()?, d.number()?)
}
fn point(v: &Json) -> Result<[Ratio; 2], Error> {
    let [x, y] = v.items()? else {
        return Err("expected two sketch coordinates");
    };
    Ok([coordinate(x)?, coordinate(y)?])
}
pub(crate) fn primitive(v: &Json) -> Result<Primitive, Error> {
    match v.field("kind")?.text()? {
        "rectangle" => {
            v.keys(&["kind", "min", "max"])?;
            Ok(Primitive::Rectangle {
                min: point(v.field("min")?)?,
                max: point(v.field("max")?)?,
            })
        }
        "circle" => {
            v.keys(&["kind", "center", "radius"])?;
            Ok(Primitive::Circle {
                center: point(v.field("center")?)?,
                radius: coordinate(v.field("radius")?)?,
            })
        }
        "ellipse" => {
            v.keys(&["kind", "center", "radii"])?;
            Ok(Primitive::Ellipse {
                center: point(v.field("center")?)?,
                radii: point(v.field("radii")?)?,
            })
        }
        _ => Err("unknown sketch primitive"),
    }
}
pub(crate) fn contour(plan: &PrimitivePlan) -> Json {
    let c = &plan.contour;
    object([
        ("degree", number(c.degree)),
        ("knots", Json::Array(c.knots.iter().map(number).collect())),
        (
            "controls",
            Json::Array(
                c.controls
                    .iter()
                    .map(|p| Json::Array(p.map(number).into()))
                    .collect(),
            ),
        ),
    ])
}
pub(crate) fn continuity(c: Continuity) -> &'static str {
    match c {
        Continuity::Corner => "G0",
        Continuity::Tangent => "G1",
        Continuity::Curvature => "G2",
        Continuity::Disconnected => "disconnected",
        Continuity::Singular => "singular",
    }
}
fn template(recipe: &Json) -> Result<PrimitivePlan, Error> {
    recipe.keys(&["family", "primitive", "operation"])?;
    if recipe.field("family")?.text()? != "primitive_sweep_v1" {
        return Err("unknown primitive planning recipe");
    }
    let op = recipe.field("operation")?;
    let operation = match op.field("kind")?.text()? {
        "extrude" => {
            op.keys(&["kind", "length_mm"])?;
            Operation::Extrude {
                length: op.field("length_mm")?.number()?,
            }
        }
        "revolve" => {
            op.keys(&["kind"])?;
            Operation::Revolve
        }
        _ => return Err("unknown primitive sweep operation"),
    };
    primitive(recipe.field("primitive")?)?.plan(operation)
}

pub fn plan(recipe: &Json) -> Result<Json, Error> {
    let p = template(recipe)?;
    let op = recipe.field("operation")?;
    let patches = p
        .section_patches
        .iter()
        .map(|s| {
            object([
                ("degrees", Json::Array(s.degrees.map(number).into())),
                (
                    "controls",
                    Json::Array(
                        s.controls
                            .iter()
                            .map(|c| Json::Array(c.map(number).into()))
                            .collect(),
                    ),
                ),
            ])
        })
        .collect();
    Ok(object([
        ("source_primitive", recipe.field("primitive")?.clone()),
        ("operation", op.clone()),
        ("contour", contour(&p)),
        (
            "junctions",
            Json::Array(
                p.junctions
                    .iter()
                    .map(|c| Json::String(continuity(*c).into()))
                    .collect(),
            ),
        ),
        ("section_patches", Json::Array(patches)),
        ("block_count", number(p.block_count())),
    ]))
}

type Geometry = (
    Vec<crate::bracket::BracketBlock>,
    Vec<(crate::incidence::BlockSide, String)>,
);
pub(crate) fn construct(recipe: &Json) -> Result<Geometry, Error> {
    let template = template(recipe)?;
    let blocks = template.build()?;
    let exterior = crate::incidence::enumerate_faces(&blocks)?.exterior;
    let mut faces = Vec::with_capacity(exterior.len());
    for side in exterior {
        let feature = if side.axis == 2 {
            if recipe.field("operation")?.field("kind")?.text()? != "extrude" {
                return Err("unpaired revolution section face");
            }
            if side.upper {
                "end".to_string()
            } else {
                "start".to_string()
            }
        } else {
            let section = side.block.split("_turn_").next().unwrap_or(&side.block);
            format!(
                "{section}:{}{}",
                ["u", "v"][side.axis],
                usize::from(side.upper)
            )
        };
        faces.push((side, feature));
    }
    Ok((blocks, faces))
}
