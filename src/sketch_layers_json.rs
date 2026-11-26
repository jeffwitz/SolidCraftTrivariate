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
    if nurbs_mode(j)?.is_some() {
        keys.push("nurbs_partition");
        if model {
            keys.push("nurbs_patches");
        }
    }
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
        force_triangulation,
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
#[derive(Clone, Copy, PartialEq)]
enum NurbsMode {
    Assemble,
    Optimize,
}

fn nurbs_mode(j: &Json) -> Result<Option<NurbsMode>, &'static str> {
    let Ok(mode) = j.field("nurbs_partition") else {
        return Ok(None);
    };
    match mode.text()? {
        "assemble" => Ok(Some(NurbsMode::Assemble)),
        "optimize" => Ok(Some(NurbsMode::Optimize)),
        _ => Err("unknown_nurbs_partition_mode"),
    }
}

fn select(job: &SketchLayers, mode: Option<NurbsMode>) -> Result<Partition, &'static str> {
    if mode == Some(NurbsMode::Optimize) {
        crate::sketch_layers::partition_nurbs_aware(job)
    } else {
        partition(job)
    }
}

/// Recover the orientation and verify the entire extracted homogeneous net up to
/// one positive common factor. Four rotations preserve the extrusion axis and sign.
fn source_rotation(span: &[[[i64; 4]; 4]], source: &[[i64; 4]; 8]) -> Result<usize, &'static str> {
    for rotation in 0..4 {
        let source_at = |s: usize, r: usize, z: usize| {
            let (u, v) = match rotation {
                0 => (s, r),
                1 => (1 - r, s),
                2 => (1 - s, 1 - r),
                _ => (r, 1 - s),
            };
            source[4 * u + 2 * v + z]
        };
        let w_src = source_at(0, 0, 0)[3];
        let w_span = span[0][0][3];
        if w_src <= 0 || w_span <= 0 {
            return Err("nonpositive_sketch_spline_weight");
        }
        let matches = span.iter().enumerate().all(|(s, section)| {
            section.iter().enumerate().all(|(j, point)| {
                let other = source_at(s, j / 2, j % 2);
                point.iter().zip(other).all(|(a, b)| {
                    i128::from(*a) * i128::from(w_src) == i128::from(b) * i128::from(w_span)
                })
            })
        });
        if matches {
            return Ok(rotation);
        }
    }
    Err("sketch_spline_extraction_mismatch")
}

fn patches(p: &Partition) -> Result<Json, &'static str> {
    let cells = p
        .blocks
        .iter()
        .map(|b| {
            if b.degrees != [1, 1, 1] {
                return Err("unsupported_sketch_spline_degree");
            }
            b.controls
                .as_slice()
                .try_into()
                .map_err(|_| "invalid_sketch_spline_control_count")
        })
        .collect::<Result<Vec<[[i64; 4]; 8]>, _>>()?;
    let assembled = crate::sketch_nurbs::assemble(&cells);
    let mut covered = vec![false; cells.len()];
    let list = assembled
        .iter()
        .map(|patch| {
            let spans = patch.cells.len();
            if spans == 0 || patch.sections.len() != spans + 1 {
                return Err("invalid_sketch_spline_shape");
            }
            let mut rotations = Vec::with_capacity(spans);
            for (span, &index) in patch.sections.windows(2).zip(&patch.cells) {
                let used = covered.get_mut(index).ok_or("unknown_sketch_spline_cell")?;
                if *used {
                    return Err("repeated_sketch_spline_cell");
                }
                *used = true;
                rotations.push(n(source_rotation(span, &cells[index])?));
            }
            let knots = |count| {
                Json::Array(
                    std::iter::once(0)
                        .chain(0..=count)
                        .chain(std::iter::once(count))
                        .map(n)
                        .collect(),
                )
            };
            let sections = patch
                .sections
                .iter()
                .map(|row| {
                    Json::Array(
                        row.iter()
                            .map(|point| Json::Array(point.iter().map(n).collect()))
                            .collect(),
                    )
                })
                .collect();
            Ok(object([
                ("source_xy_rotations", Json::Array(rotations)),
                (
                    "degrees",
                    Json::Array([1, 1, 1].into_iter().map(n).collect()),
                ),
                ("knots", Json::Array(vec![knots(spans), knots(1), knots(1)])),
                ("sections", Json::Array(sections)),
                (
                    "cells",
                    Json::Array(
                        patch
                            .cells
                            .iter()
                            .map(|&i| Json::String(p.blocks[i].id.clone()))
                            .collect(),
                    ),
                ),
            ]))
        })
        .collect::<Result<Vec<_>, _>>()?;
    if covered.contains(&false) {
        return Err("missing_sketch_spline_cell");
    }
    Ok(Json::Array(list))
}

fn ledger(p: &Partition, patch_list: Option<&Json>) -> Result<Json, &'static str> {
    let mut result = plan(p);
    if let (Some(list), Json::Object(fields)) = (patch_list, &mut result) {
        fields.insert("nurbs_patch_count".into(), n(list.items()?.len()));
        fields.insert("bezier_span_count".into(), n(p.blocks.len()));
    }
    Ok(result)
}

fn representation(p: &Partition, mode: Option<NurbsMode>) -> Result<Option<Json>, &'static str> {
    mode.map(|_| patches(p)).transpose()
}

pub(crate) fn plan_sketch(recipe: &Json) -> Result<Json, &'static str> {
    let mode = nurbs_mode(recipe)?;
    let p = select(&parse(recipe, false)?, mode)?;
    ledger(&p, representation(&p, mode)?.as_ref())
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    let mode = nurbs_mode(recipe)?;
    let p = select(&parse(recipe, false)?, mode)?;
    let patch_list = representation(&p, mode)?;
    let ledger = ledger(&p, patch_list.as_ref())?;
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
            "nurbs_partition",
        ] {
            if let Ok(v) = recipe.field(k) {
                fields.insert(k.into(), v.clone());
            }
        }
        fields.insert("partition_plan".into(), ledger);
        if let Some(list) = patch_list {
            fields.insert("nurbs_patches".into(), list);
        }
    }
    Ok(model)
}
pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    let job = parse(model, true)?;
    let mode = nurbs_mode(model)?;
    let blocks = decode_blocks(model)?;
    let replay = |p: &Partition| -> Result<bool, &'static str> {
        let patch_list = representation(p, mode)?;
        Ok(model.field("nurbs_patches").ok() == patch_list.as_ref()
            && blocks == p.blocks
            && model.field("topology")? == &topology(&p.certificate.incidence)
            && model.field("design_faces")? == &design(&p.certificate)
            && model.field("partition_plan")? == &ledger(p, patch_list.as_ref())?)
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
    if mode.is_some() {
        let p = select(&job, mode)?;
        if replay(&p)? {
            found = Some(p);
        }
    } else {
        for planner in planners {
            if let Ok(part) = planner(&job) {
                if replay(&part)? {
                    found = Some(part);
                    break;
                }
            }
        }
    }
    let p = found.ok_or("sketch_partition_replay_mismatch")?;
    let mut report = object([
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
    ]);
    if let (Some(_), Json::Object(fields)) = (mode, &mut report) {
        fields.insert(
            "nurbs_patch_count".into(),
            n(model.field("nurbs_patches")?.items()?.len()),
        );
        fields.insert("bezier_span_count".into(), n(p.blocks.len()));
    }
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn octagon_replays_four_patches_and_rejects_modified_knots() {
        let mut recipe = crate::json::parse(include_str!("../docs/site/data/bench-poly-through_octagon_high_aspect/bench_poly_through_octagon_high_aspect.recipe.json")).unwrap();
        let Json::Object(fields) = &mut recipe else {
            unreachable!()
        };
        fields.insert("nurbs_partition".into(), Json::String("assemble".into()));
        let mut model = build(&recipe).unwrap();
        let report = plan_sketch(&recipe).unwrap();
        assert_eq!(report.field("nurbs_patch_count").unwrap(), &n(4));
        assert_eq!(report.field("bezier_span_count").unwrap(), &n(8));
        assert_eq!(model.field("partition_plan").unwrap(), &report);
        audit(&model).unwrap();
        let Json::Object(fields) = &mut model else {
            unreachable!()
        };
        let Some(Json::Array(list)) = fields.get_mut("nurbs_patches") else {
            unreachable!()
        };
        let Json::Object(first) = &mut list[0] else {
            unreachable!()
        };
        first.insert("knots".into(), Json::Array(vec![]));
        assert_eq!(audit(&model), Err("sketch_partition_replay_mismatch"));
    }

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
