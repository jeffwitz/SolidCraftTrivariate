//! Exchange and exact replay of the sketch arrangement and its propagated layers.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    constraint_system::Ratio,
    json::{object, Json},
    polygon_json::pair,
    sketch_arrangement::{Role, Segment},
    sketch_layers::{partition, Partition, PlannerStrategy, SketchLayers, SketchLoop},
};
pub(crate) const FAMILY: &str = "sketch_layer_partition_v1";
pub(crate) const BOOLEAN_FAMILY: &str = "sketch_boolean_layers_v1";
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
    let boolean = j.field(if model { "schema" } else { "family" })?.text()? == BOOLEAN_FAMILY;
    let features = if boolean { "operations" } else { "cuts" };
    let mut operations = boolean.then(Vec::new);
    let mut keys = vec!["outer", features, "height_mm"];
    if j.field("construction").is_ok() {
        keys.push("construction");
    }
    if j.field("top_shift_mm").is_ok() {
        keys.push("top_shift_mm");
    }
    let strategy = if let Ok(v) = j.field("partition_strategy") {
        keys.push("partition_strategy");
        match v.text()? {
            "automatic" => PlannerStrategy::Automatic,
            "triangulate" => PlannerStrategy::Triangulate,
            "ring_planar_rerank" => PlannerStrategy::RingPlanarRerank,
            "ring_planar_lookahead" => PlannerStrategy::RingPlanarLookahead,
            _ => return Err("unknown_sketch_partition_strategy"),
        }
    } else {
        PlannerStrategy::Automatic
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
    for c in j.field(features)?.items()? {
        if let Some(ops) = &mut operations {
            c.keys(&["id", "vertices", "z_interval_mm", "operation"])?;
            ops.push(match c.field("operation")?.text()? {
                "add" => true,
                "remove" => false,
                _ => return Err("invalid_sketch_material_operation"),
            });
        } else {
            c.keys(&["id", "vertices", "z_interval_mm"])?;
        }
        cuts.push((curve(c)?, pair(c.field("z_interval_mm")?)?));
    }
    let mut guides = Vec::new();
    for g in j
        .field("construction")
        .ok()
        .map(Json::items)
        .transpose()?
        .unwrap_or(&[])
    {
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
        force_triangulation: strategy == PlannerStrategy::Triangulate,
        operations,
        cuts,
        guides,
        height: j.field("height_mm")?.number()?,
        shift: j
            .field("top_shift_mm")
            .ok()
            .map(pair)
            .transpose()?
            .unwrap_or([0, 0]),
        strategy,
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
    let mut result = object([
        ("method", Json::String(p.method.into())),
        ("vertices_mm", Json::Array(vertices)),
        ("edges", Json::Array(edges)),
        ("bounded_faces", Json::Array(faces)),
        (
            "generated_horizontal_levels_twice_mm",
            Json::Array(p.generated_levels.iter().map(n).collect()),
        ),
        ("layers", Json::Array(layers)),
    ]);
    if let Some((baseline_valid, baseline_blocks, guided_valid, guided_blocks)) =
        p.candidate_comparison
    {
        if let Json::Object(fields) = &mut result {
            fields.insert(
                "candidate_comparison".into(),
                object([
                    (
                        "baseline",
                        object([
                            ("valid_certificate", Json::Bool(baseline_valid)),
                            ("blocks", n(baseline_blocks)),
                        ]),
                    ),
                    (
                        "guided",
                        object([
                            ("valid_certificate", Json::Bool(guided_valid)),
                            ("blocks", n(guided_blocks)),
                        ]),
                    ),
                ]),
            );
        }
    }
    result
}
pub(crate) fn plan_sketch(recipe: &Json) -> Result<Json, &'static str> {
    Ok(plan(&partition(&parse(recipe, false)?)?))
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    let p = partition(&parse(recipe, false)?)?;
    let ledger = plan(&p);
    let mut model = encode(recipe.field("family")?.text()?, p.blocks, &p.certificate);
    if let Json::Object(fields) = &mut model {
        for k in [
            "outer",
            "cuts",
            "operations",
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
    let blocks = decode_blocks(model)?;
    let replay = |p: &Partition| -> Result<bool, &'static str> {
        Ok(blocks == p.blocks
            && model.field("topology")? == &topology(&p.certificate.incidence)
            && model.field("design_faces")? == &design(&p.certificate)
            && model.field("partition_plan")? == &plan(p))
    };
    // Older files can store a fully valid choice displaced by a newer proposal.
    // Replay all coefficients, design labels, incidence and ledger together.
    let planners = [
        partition,
        crate::sketch_layers::partition_before_rings,
        crate::sketch_layers::partition_rectangular_hosts,
        crate::sketch_layers::partition_classic,
        crate::sketch_layers::partition_legacy,
    ];
    let mut found = None;
    for planner in planners {
        if let Ok(part) = planner(&job) {
            if replay(&part)? {
                found = Some(part);
                break;
            }
        }
    }
    let p = found.ok_or("sketch_partition_replay_mismatch")?;
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
        ("cavities", n(p.cavities)),
        ("hausdorff_error_upper_bound_mm", n(0)),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn audit_replays_a_model_from_the_pre_guided_partition() {
        let recipe = crate::json::parse(
            r#"{"family":"sketch_layer_partition_v1","height_mm":20,"outer":{"id":"body","vertices":[[0,0],[60,0],[60,40],[0,40]]},"cuts":[{"id":"pocket","vertices":[[20,15],[40,15],[40,25],[20,25]],"z_interval_mm":[14,20]}],"construction":[]}"#,
        )
        .unwrap();
        let job = parse(&recipe, false).unwrap();
        let old = crate::sketch_layers::partition_legacy(&job).unwrap();
        let mut model = encode(FAMILY, old.blocks.clone(), &old.certificate);
        if let Json::Object(fields) = &mut model {
            for key in ["outer", "cuts", "construction", "height_mm"] {
                fields.insert(key.into(), recipe.field(key).unwrap().clone());
            }
            fields.insert("partition_plan".into(), plan(&old));
        }
        assert_eq!(old.blocks.len(), 17);
        let report = audit(&model).unwrap();
        assert_eq!(report.field("status").unwrap().text().unwrap(), "valid");
        assert_eq!(
            report.field("blocks").unwrap().number::<usize>().unwrap(),
            17
        );
    }
}
