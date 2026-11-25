use solidcraft_trivariate::{
    g2_cross_plate::construct,
    incidence::enumerate_faces,
    json::{parse, Json},
    spatial_certificate::{certify_spatial, CertificateStatus},
};

fn recipe() -> Json {
    parse(include_str!("references/g2_cross_sketches.recipe.json")).unwrap()
}

fn set_bounds(recipe: &mut Json, x: [i64; 2], y: [i64; 2]) {
    let Json::Object(fields) = recipe else {
        unreachable!()
    };
    fields.insert(
        "bounds_xyz_mm".into(),
        Json::Array(vec![
            Json::Array(x.into_iter().map(|v| Json::Number(v.to_string())).collect()),
            Json::Array(y.into_iter().map(|v| Json::Number(v.to_string())).collect()),
            Json::Array(vec![Json::Number("0".into()), Json::Number("12".into())]),
        ]),
    );
    let Json::Array(sketches) = fields.get_mut("sketches").unwrap() else {
        unreachable!()
    };
    let Json::Object(cross) = &mut sketches[1] else {
        unreachable!()
    };
    cross.insert("length_mm".into(), Json::Number(x[1].to_string()));
}

#[test]
fn constructs_the_crossed_through_bores_with_certified_blocks_and_provenance() {
    let (blocks, faces) = construct(&recipe()).unwrap();
    let incidence = enumerate_faces(&blocks).unwrap();
    assert_eq!(blocks.len(), 32);
    assert_eq!(incidence.interfaces.len(), 46);
    assert_eq!(faces.len(), incidence.exterior.len());
    assert!(blocks.iter().all(|b| {
        certify_spatial(b.degrees, &b.controls, Default::default()).status
            == CertificateStatus::Valid
    }));
    for label in [
        "nurbs_z_hole/wall",
        "rect_x/ymin",
        "rect_x/ymax",
        "rect_x/zmin",
        "rect_x/zmax",
        "stock/xmin",
        "stock/xmax",
        "stock/ymin",
        "stock/ymax",
        "stock/zmin",
        "stock/zmax",
    ] {
        assert!(
            faces.iter().any(|(_, feature)| feature == label),
            "missing {label}"
        );
    }
}

#[test]
fn stock_extent_is_derived_from_the_recipe() {
    let mut input = recipe();
    set_bounds(&mut input, [0, 120], [0, 55]);
    let (blocks, faces) = construct(&input).unwrap();
    assert_eq!(blocks.len(), 32);
    assert!(faces.iter().any(|(_, f)| f == "stock/xmax"));
    assert!(faces.iter().any(|(_, f)| f == "stock/ymax"));
    assert!(blocks.iter().all(|b| {
        certify_spatial(b.degrees, &b.controls, Default::default()).status
            == CertificateStatus::Valid
    }));
}

fn replace_controls(recipe: &mut Json, controls: Vec<Json>) {
    let Json::Object(fields) = recipe else {
        unreachable!()
    };
    let Json::Array(sketches) = fields.get_mut("sketches").unwrap() else {
        unreachable!()
    };
    let Json::Object(sketch) = &mut sketches[0] else {
        unreachable!()
    };
    let Json::Array(loops) = sketch.get_mut("loops").unwrap() else {
        unreachable!()
    };
    let Json::Object(loop_fields) = &mut loops[0] else {
        unreachable!()
    };
    let Json::Object(profile) = loop_fields.get_mut("curve").unwrap() else {
        unreachable!()
    };
    profile.insert("controls".into(), Json::Array(controls));
}

fn replace_knots(recipe: &mut Json, knots: Vec<i64>) {
    let Json::Object(fields) = recipe else {
        unreachable!()
    };
    let Json::Array(sketches) = fields.get_mut("sketches").unwrap() else {
        unreachable!()
    };
    let Json::Object(sketch) = &mut sketches[0] else {
        unreachable!()
    };
    let Json::Array(loops) = sketch.get_mut("loops").unwrap() else {
        unreachable!()
    };
    let Json::Object(loop_fields) = &mut loops[0] else {
        unreachable!()
    };
    let Json::Object(profile) = loop_fields.get_mut("curve").unwrap() else {
        unreachable!()
    };
    profile.insert(
        "knots".into(),
        Json::Array(
            knots
                .into_iter()
                .map(|v| Json::Number(v.to_string()))
                .collect(),
        ),
    );
}

fn controls(recipe: &Json) -> Vec<Vec<i64>> {
    recipe.field("sketches").unwrap().items().unwrap()[0]
        .field("loops")
        .unwrap()
        .items()
        .unwrap()[0]
        .field("curve")
        .unwrap()
        .field("controls")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|p| {
            p.items()
                .unwrap()
                .iter()
                .map(|v| v.number().unwrap())
                .collect()
        })
        .collect()
}

fn json_controls(values: &[Vec<i64>]) -> Vec<Json> {
    values
        .iter()
        .map(|p| Json::Array(p.iter().map(|v| Json::Number(v.to_string())).collect()))
        .collect()
}

fn unlabeled(
    mut blocks: Vec<solidcraft_trivariate::bracket::BracketBlock>,
) -> Vec<solidcraft_trivariate::bracket::BracketBlock> {
    for block in &mut blocks {
        block.id.clear();
    }
    blocks.sort_by(|a, b| a.controls.cmp(&b.controls));
    blocks
}

#[test]
fn cyclic_reindexing_and_orientation_reversal_preserve_the_partition() {
    let original = recipe();
    let expected = unlabeled(construct(&original).unwrap().0);
    for reverse in [false, true] {
        let mut changed = original.clone();
        let mut poles = controls(&original);
        poles.truncate(8);
        if reverse {
            poles[1..].reverse();
        } else {
            poles.rotate_left(1);
        }
        let periodic: Vec<_> = (0..3).map(|i| poles[i].clone()).collect();
        poles.extend(periodic);
        replace_controls(&mut changed, json_controls(&poles));
        assert_eq!(expected, unlabeled(construct(&changed).unwrap().0));
    }
}

#[test]
fn exact_internal_knot_insertion_preserves_the_partition() {
    let original = recipe();
    let expected = unlabeled(construct(&original).unwrap().0);
    let p = controls(&original);
    // Insert u=3 into the degree-3 periodic source. The resulting rational
    // blends have denominators 3, so scale the complete homogeneous net by 3.
    let mut q = vec![vec![0_i64; 3]; p.len() + 1];
    for i in 0..=3 {
        q[i] = p[i].iter().map(|v| 3 * v).collect();
    }
    let mut knots: Vec<i64> = (-3..=11).collect();
    for i in 5..=10 {
        q[i + 1] = p[i].iter().map(|v| 3 * v).collect();
    }
    for (i, alpha_num, alpha_den) in [(4, 2_i64, 3_i64), (5, 1, 3)] {
        for d in 0..3 {
            q[i][d] = alpha_num * p[i][d] + (alpha_den - alpha_num) * p[i - 1][d];
        }
    }
    knots.insert(7, 3);
    let mut changed = original;
    replace_controls(&mut changed, json_controls(&q));
    replace_knots(&mut changed, knots);
    assert_eq!(expected, unlabeled(construct(&changed).unwrap().0));
}
