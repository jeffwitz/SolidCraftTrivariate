//! Exchange and exact replay of the sketch arrangement and its propagated layers.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    constraint_system::Ratio,
    json::{object, Json},
    polygon_json::pair,
    sketch_arrangement::{Role, Segment},
    sketch_layers::{partition, Partition, SketchLayers, SketchLoop},
};
pub(crate) const FAMILY: &str = "sketch_layer_partition_v1";
pub(crate) const METHOD: &str = "certified_sketch_layer_partition";
fn n<T: ToString>(x: T) -> Json {
    Json::Number(x.to_string())
}
fn fraction(v: Ratio) -> Json {
    Json::Array(vec![n(v.numerator()), n(v.denominator())])
}
fn curve(j: &Json) -> Result<SketchLoop, &'static str> {
    Ok(SketchLoop {
        id: j.field("id")?.text()?.into(),
        vertices: j
            .field("vertices")?
            .items()?
            .iter()
            .map(pair)
            .collect::<Result<_, _>>()?,
    })
}
fn parse(j: &Json, model: bool) -> Result<SketchLayers, &'static str> {
    let mut keys = vec!["outer", "cuts", "construction", "height_mm"];
    if j.field("top_shift_mm").is_ok() {
        keys.push("top_shift_mm");
    }
    let force_triangulation = if let Ok(v) = j.field("partition_strategy") {
        keys.push("partition_strategy");
        match v.text()? {
            "automatic" => false,
            "triangulate" => true,
            _ => return Err("unknown_sketch_partition_strategy"),
        }
    } else {
        false
    };
    if model {
        keys.extend([
            "schema",
            "blocks",
            "design_faces",
            "topology",
            "partition_plan",
        ]);
    } else {
        keys.push("family");
    }
    j.keys(&keys)?;
    let outer = j.field("outer")?;
    outer.keys(&["id", "vertices"])?;
    let mut cuts = Vec::new();
    for c in j.field("cuts")?.items()? {
        c.keys(&["id", "vertices", "z_interval_mm"])?;
        cuts.push((curve(c)?, pair(c.field("z_interval_mm")?)?));
    }
    let mut guides = Vec::new();
    for g in j.field("construction")?.items()? {
        g.keys(&["id", "a", "b"])?;
        guides.push(Segment {
            id: g.field("id")?.text()?.into(),
            a: pair(g.field("a")?)?,
            b: pair(g.field("b")?)?,
            role: Role::Construction,
        });
    }
    Ok(SketchLayers {
        outer: curve(outer)?,
        force_triangulation,
        cuts,
        guides,
        height: j.field("height_mm")?.number()?,
        shift: j
            .field("top_shift_mm")
            .ok()
            .map(pair)
            .transpose()?
            .unwrap_or([0, 0]),
    })
}
fn plan(p: &Partition) -> Json {
    let g = &p.graph;
    let vertices = g
        .vertices
        .iter()
        .map(|pt| {
            Json::Array(
                pt.iter()
                    .copied()
                    .map(|v| {
                        fraction(
                            v.checked_mul(Ratio::new(1, 2).unwrap())
                                .expect("bounded sketch coordinate"),
                        )
                    })
                    .collect(),
            )
        })
        .collect();
    let edges = g
        .edges
        .iter()
        .map(|e| {
            object([
                ("vertices", Json::Array(e.vertices.iter().map(n).collect())),
                (
                    "sources",
                    Json::Array(
                        e.provenance
                            .iter()
                            .map(|s| {
                                object([
                                    ("id", Json::String(s.segment_id.clone())),
                                    ("forward", Json::Bool(s.forward)),
                                    (
                                        "role",
                                        Json::String(
                                            match s.role {
                                                Role::Boundary => "boundary",
                                                Role::Construction => "construction",
                                            }
                                            .into(),
                                        ),
                                    ),
                                ])
                            })
                            .collect(),
                    ),
                ),
            ])
        })
        .collect();
    let faces = g
        .faces
        .iter()
        .map(|f| Json::Array(f.vertices.iter().map(n).collect()))
        .collect();
    let layers = p
        .layer_counts
        .iter()
        .map(|(z, b)| {
            object([
                ("interval_mm", Json::Array(z.iter().map(n).collect())),
                ("material_blocks", n(b)),
            ])
        })
        .collect();
    object([
        (
            "method",
            Json::String(
                if p.rectangular_grid {
                    "orthogonal_feature_footprint"
                } else {
                    "triangulated_feature_footprint"
                }
                .into(),
            ),
        ),
        ("vertices_mm", Json::Array(vertices)),
        ("edges", Json::Array(edges)),
        ("bounded_faces", Json::Array(faces)),
        (
            "generated_horizontal_levels_twice_mm",
            Json::Array(p.generated_levels.iter().map(n).collect()),
        ),
        ("layers", Json::Array(layers)),
    ])
}
pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    let p = partition(&parse(recipe, false)?)?;
    let ledger = plan(&p);
    let mut model = encode(FAMILY, p.blocks, &p.certificate);
    if let Json::Object(fields) = &mut model {
        for k in [
            "outer",
            "cuts",
            "construction",
            "height_mm",
            "top_shift_mm",
            "partition_strategy",
        ] {
            if let Ok(v) = recipe.field(k) {
                fields.insert(k.into(), v.clone());
            }
        }
        fields.insert("partition_plan".into(), ledger);
    }
    Ok(model)
}
pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    let job = parse(model, true)?;
    let cavities = job
        .cuts
        .iter()
        .filter(|(_, z)| z[0] > 0 && z[1] < job.height)
        .count();
    let p = partition(&job)?;
    if decode_blocks(model)? != p.blocks
        || model.field("topology")? != &topology(&p.certificate.incidence)
        || model.field("design_faces")? != &design(&p.certificate)
        || model.field("partition_plan")? != &plan(&p)
    {
        return Err("sketch_partition_replay_mismatch");
    }
    Ok(object([
        ("status", Json::String("valid".into())),
        ("method", Json::String(METHOD.into())),
        ("blocks", n(p.blocks.len())),
        ("shared_faces", n(p.certificate.incidence.interfaces.len())),
        ("exterior_sides", n(p.certificate.incidence.exterior.len())),
        (
            "determinant_lower_bound_mm3",
            n(p.certificate.determinant_lower_bound_mm3),
        ),
        ("volume_exact_mm3", fraction(p.volume)),
        ("bounded_sketch_faces", n(p.planar_faces)),
        ("cavities", n(cavities)),
        ("hausdorff_error_upper_bound_mm", n(0)),
    ]))
}
