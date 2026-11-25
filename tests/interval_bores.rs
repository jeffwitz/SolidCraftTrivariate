use solidcraft_trivariate::{
    circular_polygon::CircularPolygonBore,
    interval_bores::{build_interval_bore_prism, IntervalBore, IntervalBorePrism},
    interval_bores_audit::audit_interval_bore_prism,
    json::{self, Json},
    perforated_prism::LoopProfile,
    polygon_prism::Vertex,
};
use std::collections::BTreeSet;

fn integer(value: &Json) -> Result<i32, &'static str> {
    value.number()
}

fn pair(value: &Json) -> Result<[i32; 2], &'static str> {
    let [a, b] = value.items().unwrap() else {
        return Err("expected integer pair");
    };
    Ok([integer(a)?, integer(b)?])
}

fn recipe(value: &Json) -> Result<IntervalBorePrism, &'static str> {
    let outer = value.field("outer").unwrap();
    let vertices = outer
        .field("vertices")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|vertex| -> Result<_, &'static str> {
            Ok(Vertex {
                id: vertex.field("id")?.text()?.into(),
                xy_mm: pair(vertex.field("xy_mm")?)?,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let bores = value
        .field("bores")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|bore| -> Result<_, &'static str> {
            Ok(IntervalBore {
                bore: CircularPolygonBore {
                    id: bore.field("id")?.text()?.into(),
                    center_mm: pair(bore.field("center_mm")?)?,
                    radius_mm: integer(bore.field("radius_mm")?)?,
                    guard_mm: integer(bore.field("guard_mm")?)?,
                },
                z_interval_mm: pair(bore.field("z_interval_mm")?)?,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(IntervalBorePrism {
        outer: LoopProfile {
            id: outer.field("id")?.text()?.into(),
            vertices,
        },
        bores,
        height_mm: integer(value.field("height_mm")?)?,
        top_shift_mm: pair(value.field("top_shift_mm")?)?,
    })
}

#[test]
fn frozen_interval_bore_corpus_decisions_and_partitions() {
    let corpus = json::parse(include_str!("references/interval_bores_corpus.json")).unwrap();
    let mut digests = BTreeSet::new();
    for case in corpus.items().unwrap() {
        let id = case.field("id").unwrap().text().unwrap();
        let expected = case.field("expected").unwrap() == &Json::Bool(true);
        let Ok(model) = recipe(case.field("recipe").unwrap()) else {
            assert!(
                !expected,
                "nonintegral/malformed case {id} was expected to be admissible"
            );
            continue;
        };
        let result = build_interval_bore_prism(&model);
        assert_eq!(result.is_ok(), expected, "case {id}: {result:?}");
        let Ok(blocks) = result else { continue };
        let certificate = audit_interval_bore_prism(&model, &blocks)
            .unwrap_or_else(|reason| panic!("audit rejected admissible case {id}: {reason}"));
        assert!(
            certificate.certificate.determinant_lower_bound_mm3 > 0.0,
            "case {id}"
        );
        let expected_removed = model
            .bores
            .iter()
            .map(|bore| {
                i64::from(bore.bore.radius_mm).pow(2)
                    * i64::from(bore.z_interval_mm[1] - bore.z_interval_mm[0])
            })
            .sum::<i64>();
        assert_eq!(
            certificate.removed_pi_volume_coefficient_mm3, expected_removed,
            "case {id}"
        );
        assert_eq!(
            certificate.through_bores,
            model
                .bores
                .iter()
                .filter(|b| b.z_interval_mm == [0, model.height_mm])
                .count(),
            "case {id}"
        );
        assert_eq!(
            certificate.enclosed_voids,
            model
                .bores
                .iter()
                .filter(|b| b.z_interval_mm[0] > 0 && b.z_interval_mm[1] < model.height_mm)
                .count(),
            "case {id}"
        );
        assert!(!blocks.is_empty(), "case {id} has no cells");
        assert!(blocks.windows(2).all(|pair| pair[0].id < pair[1].id));
        assert!(blocks.iter().all(|block| {
            let expected_count = block.degrees.into_iter().map(|d| d + 1).product::<usize>();
            block.controls.len() == expected_count && block.controls.iter().all(|p| p[3] > 0)
        }));
        let digest = blocks
            .iter()
            .map(|block| format!("{}:{:?}:{:?}", block.id, block.degrees, block.controls))
            .collect::<Vec<_>>()
            .join("\n");
        if id.ends_with("_reorder") {
            let base = id.trim_end_matches("_reorder");
            assert!(digests.contains(base), "missing base partition for {id}");
        } else {
            digests.insert(id.to_owned());
        }
        if id == "through" {
            assert!(blocks
                .iter()
                .all(|b| !b.id.contains("fill_sector") && !b.id.contains("core_sector")));
        }
        if id == "blind_top" {
            assert!(blocks
                .iter()
                .any(|b| b.id.contains("fill_sector") && b.id.ends_with("_z0_4")));
            assert!(!blocks
                .iter()
                .any(|b| b.id.contains("fill_sector") && b.id.ends_with("_z4_12")));
        }
        if id == "blind_bottom" {
            assert!(!blocks
                .iter()
                .any(|b| b.id.contains("fill_sector") && b.id.ends_with("_z0_7")));
            assert!(blocks
                .iter()
                .any(|b| b.id.contains("fill_sector") && b.id.ends_with("_z7_12")));
        }
        if id.ends_with("_reorder") {
            let base = id.trim_end_matches("_reorder");
            let base_case = corpus
                .items()
                .unwrap()
                .iter()
                .find(|candidate| candidate.field("id").unwrap().text().unwrap() == base)
                .unwrap();
            let base_recipe = recipe(base_case.field("recipe").unwrap()).unwrap();
            let base_blocks = build_interval_bore_prism(&base_recipe).unwrap();
            assert_eq!(blocks, base_blocks, "input bore ordering changed {id}");
            assert_eq!(
                digest,
                base_blocks
                    .iter()
                    .map(|block| format!("{}:{:?}:{:?}", block.id, block.degrees, block.controls))
                    .collect::<Vec<_>>()
                    .join("\n")
            );
        }
    }
}

#[test]
fn core_orientation_regression_has_positive_certified_jacobian() {
    let regression = json::parse(include_str!(
        "regressions/interval_bore_core_orientation.json"
    ))
    .unwrap();
    let model = recipe(regression.field("recipe").unwrap()).unwrap();
    let blocks = build_interval_bore_prism(&model).unwrap();
    let certificate = audit_interval_bore_prism(&model, &blocks).unwrap();
    assert!(
        certificate.certificate.determinant_lower_bound_mm3 > 0.0,
        "old constructor failed this assertion with a nonpositive/uncertified core Jacobian"
    );
    assert_eq!(certificate.removed_pi_volume_coefficient_mm3, 144);
    assert_eq!(certificate.through_bores, 0);
    assert_eq!(certificate.enclosed_voids, 1);
    assert_eq!(
        blocks
            .iter()
            .filter(|b| b.id.contains("core_sector"))
            .count(),
        8
    );
}

#[test]
fn ring_orientation_regression_keeps_exact_radial_seams() {
    let regression = json::parse(include_str!(
        "regressions/interval_bore_ring_orientation.json"
    ))
    .unwrap();
    let model = recipe(regression.field("recipe").unwrap()).unwrap();
    let blocks = build_interval_bore_prism(&model).unwrap();
    assert!(blocks.iter().all(|b| !b.id.contains("core_sector")));
    let certificate = audit_interval_bore_prism(&model, &blocks).unwrap();
    assert!(certificate.certificate.determinant_lower_bound_mm3 > 0.0);
    assert_eq!(certificate.through_bores, 1);
}
