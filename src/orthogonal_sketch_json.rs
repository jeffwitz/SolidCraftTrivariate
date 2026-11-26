//! Versioned design-feature input and inspectable multi-axis sweep plan.
use crate::{
    bracket_json::{decode_blocks, design, encode, topology},
    json::{object, Json},
    orthogonal_sketch::{self, Bore, Fillet, OrthogonalSketch, SweepPartition},
    polygon_json::pair,
    sketch_layers::SketchLoop,
};
pub(crate) const FAMILY: &str = "orthogonal_sketch_features_v1";
pub(crate) const METHOD: &str = "certified_orthogonal_sketch_features";
pub(crate) const COMPACT_METHOD: &str = "four_sector_local_sweeps";
pub(crate) const PROPAGATED_METHOD: &str = "propagated_isoparametric_cuts";
const PARAMETERS: [&str; 4] = ["outer", "height_mm", "fillets", "bores"];
fn n(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn parse(j: &Json) -> Result<OrthogonalSketch, &'static str> {
    let p = j.field("outer")?;
    p.keys(&["id", "vertices"])?;
    let outer = SketchLoop {
        id: p.field("id")?.text()?.to_owned(),
        vertices: p
            .field("vertices")?
            .items()?
            .iter()
            .map(pair)
            .collect::<Result<_, _>>()?,
    };
    let mut fillets = Vec::new();
    for v in j.field("fillets")?.items()? {
        v.keys(&["id", "vertex", "radius_mm"])?;
        fillets.push(Fillet {
            id: v.field("id")?.text()?.into(),
            vertex: pair(v.field("vertex")?)?,
            radius: v.field("radius_mm")?.number()?,
        });
    }
    let mut bores = Vec::new();
    for v in j.field("bores")?.items()? {
        v.keys(&["id", "axis", "center_mm", "radius_mm"])?;
        let axis = match v.field("axis")?.text()? {
            "x" => 0,
            "y" => 1,
            "z" => 2,
            _ => return Err("unsupported_bore_axis"),
        };
        let [x, y, z] = v.field("center_mm")?.items()? else {
            return Err("expected_bore_seed_xyz");
        };
        bores.push(Bore {
            id: v.field("id")?.text()?.into(),
            axis,
            center: [x.number()?, y.number()?, z.number()?],
            radius: v.field("radius_mm")?.number()?,
        });
    }
    Ok(OrthogonalSketch {
        outer,
        height: j.field("height_mm")?.number()?,
        fillets,
        bores,
    })
}
fn plan(p: &SweepPartition) -> Json {
    let mut out = object([
        (
            "method",
            Json::String(
                if p.compact {
                    COMPACT_METHOD
                } else {
                    "axis_planes_and_disjoint_local_sweeps"
                }
                .into(),
            ),
        ),
        (
            "axis_planes_mm",
            Json::Array(
                p.planes
                    .iter()
                    .map(|axis| Json::Array(axis.iter().map(n).collect()))
                    .collect(),
            ),
        ),
        (
            "feature_guards",
            Json::Array(
                p.guards
                    .iter()
                    .map(|(name, b)| {
                        object([
                            ("feature", Json::String(name.clone())),
                            (
                                "bounds_mm",
                                Json::Array(
                                    b.iter()
                                        .map(|v| Json::Array(v.iter().map(n).collect()))
                                        .collect(),
                                ),
                            ),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    if let (Some(comparison), Json::Object(fields)) = (&p.cut_comparison, &mut out) {
        let score = |count: Result<usize, &'static str>| match count {
            Ok(count) => object([("blocks", n(count))]),
            Err(reason) => object([("refusal", Json::String(reason.into()))]),
        };
        fields.insert(
            "cut_comparison".into(),
            object([
                ("previous", score(comparison.prior)),
                ("propagated", score(comparison.propagated)),
            ]),
        );
        if !p.cut_rounds.is_empty() {
            fields.insert(
                "cut_rounds".into(),
                Json::Array(
                    p.cut_rounds
                        .iter()
                        .map(|r| {
                            object([
                                ("requests", n(r[0])),
                                ("closed_cuts", n(r[1])),
                                ("blocks", n(r[2])),
                            ])
                        })
                        .collect(),
                ),
            );
            fields.insert(
                "cut_sources".into(),
                Json::Array(
                    p.cut_sources
                        .iter()
                        .map(|c| {
                            object([
                                ("round", n(c.round)),
                                ("block", Json::String(c.block.clone())),
                                ("axis", n(c.axis)),
                                (
                                    "parameter",
                                    Json::Array(vec![n(c.numerator), n(c.denominator)]),
                                ),
                                (
                                    "sources",
                                    Json::Array(
                                        c.sources.iter().cloned().map(Json::String).collect(),
                                    ),
                                ),
                            ])
                        })
                        .collect(),
                ),
            );
        }
        if comparison.selected {
            let method = if p.cut_rounds.is_empty() {
                "local_bore_hosts"
            } else {
                PROPAGATED_METHOD
            };
            fields.insert("method".into(), Json::String(method.into()));
        }
    }
    if let (Some(reason), Json::Object(fields)) = (p.fallback_reason.as_ref(), &mut out) {
        fields.insert("fallback_reason".into(), Json::String(reason.clone()));
    }
    if let (Some((center, axes, r, g)), Json::Object(f)) = (p.junction, &mut out) {
        f.insert(
            "method".into(),
            Json::String("equal_radius_intersecting_axis_sweeps".into()),
        );
        f.insert(
            "axis_intersection_mm".into(),
            Json::Array(center.iter().map(n).collect()),
        );
        f.insert(
            "bore_axes".into(),
            Json::Array(axes.iter().map(n).collect()),
        );
        f.insert("radius_mm".into(), n(r));
        f.insert("guard_mm".into(), n(g));
    }
    out
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&[&["family"][..], &PARAMETERS].concat())?;
    let p = orthogonal_sketch::build(&parse(recipe)?)?;
    let ledger = plan(&p);
    let mut out = encode(FAMILY, p.blocks, &p.certificate);
    if let Json::Object(f) = &mut out {
        for k in PARAMETERS {
            f.insert(k.into(), recipe.field(k)?.clone());
        }
        f.insert("partition_plan".into(), ledger);
    }
    Ok(out)
}
pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    model.keys(
        &[
            &[
                "schema",
                "blocks",
                "topology",
                "design_faces",
                "partition_plan",
            ][..],
            &PARAMETERS,
        ]
        .concat(),
    )?;
    let j = parse(model)?;
    // Models written before fallback reasons were recorded retain their historical replay path.
    let compact = model.field("partition_plan")?.field("method")?.text()? == COMPACT_METHOD;
    let has_fallback_reason = model
        .field("partition_plan")?
        .field("fallback_reason")
        .is_ok();
    let compared = model
        .field("partition_plan")?
        .field("cut_comparison")
        .is_ok();
    let p = orthogonal_sketch::build_with_strategy(
        &j,
        compact || has_fallback_reason || compared,
        compared,
    )?;
    let mut expected_plan = plan(&p);
    if !has_fallback_reason {
        if let Json::Object(fields) = &mut expected_plan {
            fields.remove("fallback_reason");
        }
    }
    if decode_blocks(model)? != p.blocks
        || model.field("partition_plan")? != &expected_plan
        || model.field("topology")? != &topology(&p.certificate.incidence)
        || model.field("design_faces")? != &design(&p.certificate)
    {
        return Err("orthogonal_feature_partition_replay_mismatch");
    }
    let mut volume = object([
        ("twice_rational_part_mm3", n(p.twice_constant)),
        ("pi_coefficient_mm3", Json::Array(vec![n(p.four_pi), n(4)])),
    ]);
    if let (Some((_, _, r, _)), Json::Object(v)) = (p.junction, &mut volume) {
        let thirds = 3 * (p.twice_constant / 2) + 16 * i64::from(r).pow(3);
        v.insert(
            "rational_part_mm3".into(),
            Json::Array(vec![n(thirds), n(3)]),
        );
        if thirds % 3 == 0 {
            v.insert("twice_rational_part_mm3".into(), n(2 * (thirds / 3)));
        } else {
            v.remove("twice_rational_part_mm3");
        }
    }
    Ok(object([
        ("status", Json::String("valid".into())),
        ("method", Json::String(METHOD.into())),
        ("blocks", n(p.blocks.len())),
        ("shared_faces", n(p.certificate.incidence.interfaces.len())),
        ("exterior_sides", n(p.certificate.incidence.exterior.len())),
        ("cavities", n(0)),
        ("through_bores", n(j.bores.len())),
        (
            "determinant_lower_bound_mm3",
            n(p.certificate.determinant_lower_bound_mm3),
        ),
        ("hausdorff_error_upper_bound_mm", n(0)),
        ("volume_exact", volume),
    ]))
}

#[cfg(test)]
mod replay_tests {
    use super::*;

    #[test]
    fn archived_fallback_keeps_its_grid_and_rejects_a_forged_reason() {
        let recipe = crate::json::parse(
            r#"{
            "family":"orthogonal_sketch_features_v1",
            "outer":{"id":"body","vertices":[[0,0],[60,0],[60,8],[8,8],[8,50],[0,50]]},
            "height_mm":40,"fillets":[{"id":"f","vertex":[8,8],"radius_mm":6}],
            "bores":[{"id":"low","axis":"y","center_mm":[38,4,10],"radius_mm":3},
                     {"id":"high","axis":"y","center_mm":[38,4,30],"radius_mm":3},
                     {"id":"upright","axis":"x","center_mm":[4,32,20],"radius_mm":4}]}"#,
        )
        .unwrap();
        let job = parse(&recipe).unwrap();
        let old = orthogonal_sketch::build_with_strategy(&job, true, false).unwrap();
        assert_eq!(old.blocks.len(), 112);
        let ledger = plan(&old);
        assert!(ledger.field("fallback_reason").is_ok());
        let mut model = encode(FAMILY, old.blocks, &old.certificate);
        if let Json::Object(f) = &mut model {
            for k in PARAMETERS {
                f.insert(k.into(), recipe.field(k).unwrap().clone());
            }
            f.insert("partition_plan".into(), ledger);
        }
        audit(&model).unwrap();
        if let Json::Object(f) = &mut model {
            if let Some(Json::Object(p)) = f.get_mut("partition_plan") {
                p.insert("fallback_reason".into(), Json::String("made_up".into()));
            }
        }
        assert_eq!(
            audit(&model),
            Err("orthogonal_feature_partition_replay_mismatch")
        );
    }
}
