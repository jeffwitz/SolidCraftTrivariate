//! Strict interchange for a certified ordered cartesian box history.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    json::{object, Json},
    orthogonal_audit::audit_orthogonal_program,
    orthogonal_program::{build_orthogonal_program, Bounds, BoxFeature, BoxOp, OrthogonalProgram},
};
pub(crate) const FAMILY: &str = "orthogonal_feature_program_v1";
fn n(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn bounds(x: &Json) -> Result<Bounds, &'static str> {
    let mut result = [[0; 2]; 3];
    let axes = x.items()?;
    if axes.len() != 3 {
        return Err("orthogonal bounds require three intervals");
    }
    for (out, input) in result.iter_mut().zip(axes) {
        let [a, b] = input.items()? else {
            return Err("orthogonal interval needs two endpoints");
        };
        *out = [a.number()?, b.number()?];
    }
    Ok(result)
}
fn parse(j: &Json) -> Result<OrthogonalProgram, &'static str> {
    let input = j.field("features")?.items()?;
    if input.len() > 8 {
        return Err("orthogonal_program_range");
    }
    let mut features = Vec::new();
    for f in input {
        f.keys(&["id", "operation", "bounds_mm"])?;
        let operation = match f.field("operation")?.text()? {
            "add" => BoxOp::Add,
            "cut" => BoxOp::Cut,
            _ => return Err("unknown orthogonal feature operation"),
        };
        features.push(BoxFeature {
            id: f.field("id")?.text()?.into(),
            operation,
            bounds_mm: bounds(f.field("bounds_mm")?)?,
        });
    }
    Ok(OrthogonalProgram {
        base_mm: bounds(j.field("base_mm")?)?,
        features,
    })
}
pub(crate) fn build(j: &Json) -> Result<Json, &'static str> {
    j.keys(&["family", "base_mm", "features"])?;
    let program = parse(j)?;
    let blocks = build_orthogonal_program(&program)?;
    let cert = audit_orthogonal_program(&program, &blocks)?;
    let mut result = encode(FAMILY, blocks, &cert.certificate);
    if let Json::Object(f) = &mut result {
        for key in ["base_mm", "features"] {
            f.insert(key.into(), j.field(key)?.clone());
        }
    }
    Ok(result)
}
pub(crate) fn audit(j: &Json) -> Result<Json, &'static str> {
    j.keys(&[
        "schema",
        "base_mm",
        "features",
        "blocks",
        "topology",
        "design_faces",
    ])?;
    let program = parse(j)?;
    let blocks = decode_blocks(j)?;
    let result = audit_orthogonal_program(&program, &blocks)?;
    let c = &result.certificate;
    if j.field("topology")? != &topology(&c.incidence) {
        return Err("orthogonal incidence mismatch");
    }
    if j.field("design_faces")? != &design(c) {
        return Err("orthogonal boundary provenance mismatch");
    }
    Ok(object([
        ("status", Json::String("valid".into())),
        (
            "method",
            Json::String("certified_orthogonal_feature_program".into()),
        ),
        ("blocks", n(blocks.len())),
        ("shared_faces", n(c.incidence.interfaces.len())),
        ("exterior_sides", n(c.incidence.exterior.len())),
        (
            "determinant_lower_bound_mm3",
            n(c.determinant_lower_bound_mm3),
        ),
        ("volume_mm3", n(result.volume_mm3)),
        ("cavities", n(result.cavities)),
        ("minimum_x_mm", n(result.minimum_x_mm)),
    ]))
}
