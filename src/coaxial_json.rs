//! Versioned recipes and checked persistence of ordered coaxial feature histories.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    coaxial_audit::audit_coaxial_program,
    coaxial_program::{build_coaxial_program, CoaxialProgram, CylinderFeature, CylinderOp},
    json::{object, Json},
};

pub(crate) const FAMILY: &str = "coaxial_feature_program_v1";
fn n(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}

fn program(data: &Json) -> Result<CoaxialProgram, &'static str> {
    let items = data.field("features")?.items()?;
    if items.len() > 8 {
        return Err("coaxial_program_range");
    }
    let mut features = Vec::with_capacity(items.len());
    for item in items {
        item.keys(&["id", "operation", "radius_mm", "z_mm"])?;
        let op = match item.field("operation")?.text()? {
            "add" => CylinderOp::Add,
            "cut" => CylinderOp::Cut,
            _ => return Err("coaxial_unknown_operation"),
        };
        let [lo, hi] = item.field("z_mm")?.items()? else {
            return Err("coaxial_axial_interval");
        };
        features.push(CylinderFeature {
            id: item.field("id")?.text()?.into(),
            op,
            radius_mm: item.field("radius_mm")?.number()?,
            z_mm: [lo.number()?, hi.number()?],
        });
    }
    Ok(CoaxialProgram {
        base_height_mm: data.field("base_height_mm")?.number()?,
        features,
    })
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&["family", "base_height_mm", "features"])?;
    let source = program(recipe)?;
    let blocks = build_coaxial_program(&source)?;
    let report = audit_coaxial_program(&source, &blocks)?;
    let mut model = encode(FAMILY, blocks, &report.certificate);
    if let Json::Object(fields) = &mut model {
        for key in ["base_height_mm", "features"] {
            fields.insert(key.into(), recipe.field(key)?.clone());
        }
    }
    Ok(model)
}

pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    model.keys(&[
        "schema",
        "base_height_mm",
        "features",
        "blocks",
        "topology",
        "design_faces",
    ])?;
    let source = program(model)?;
    let blocks = decode_blocks(model)?;
    let report = audit_coaxial_program(&source, &blocks)?;
    let c = &report.certificate;
    if model.field("topology")? != &topology(&c.incidence) {
        return Err("coaxial_topology_mismatch");
    }
    if model.field("design_faces")? != &design(c) {
        return Err("coaxial_design_face_mismatch");
    }
    Ok(object([
        ("status", Json::String("valid".into())),
        (
            "method",
            Json::String("certified_coaxial_feature_program".into()),
        ),
        ("blocks", n(blocks.len())),
        ("shared_faces", n(c.incidence.interfaces.len())),
        ("exterior_sides", n(c.incidence.exterior.len())),
        (
            "determinant_lower_bound_mm3",
            n(c.determinant_lower_bound_mm3),
        ),
        ("radial_error_upper_bound_mm", n(c.radial_error_upper_mm)),
        ("ideal_volume_mm3", n(report.ideal_volume_mm3)),
        ("cavities", n(report.cavities)),
    ]))
}
