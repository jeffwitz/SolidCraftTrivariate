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
            Json::String("axis_planes_and_disjoint_local_sweeps".into()),
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
    let p = orthogonal_sketch::build(&j)?;
    if decode_blocks(model)? != p.blocks
        || model.field("partition_plan")? != &plan(&p)
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
