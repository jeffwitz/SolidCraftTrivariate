use solidcraft_trivariate::{json::parse, nurbs_edge_loops::normalize, nurbs_multiloop::construct};

fn fixture() -> solidcraft_trivariate::json::Json {
    parse(include_str!("references/nurbs_edge_loops.recipe.json")).unwrap()
}

#[test]
fn discovers_two_components_from_unordered_seeded_edges() {
    let input = fixture();
    let normalized = normalize(&input).unwrap();
    assert_eq!(
        normalized.field("family").unwrap().text().unwrap(),
        "nurbs_multiloop_v1"
    );
    let loops = normalized.field("loops").unwrap().items().unwrap();
    assert_eq!(loops.len(), 2);
    assert_eq!(
        loops[0].field("id").unwrap().text().unwrap(),
        "north-free-profile"
    );
    assert_eq!(
        loops[1].field("id").unwrap().text().unwrap(),
        "south-curved-pocket"
    );
    let lower = loops[0].field("lower").unwrap();
    assert_eq!(
        lower.field("knots").unwrap(),
        &parse("[0,0,0,0,1,2,4,4,4,4]").unwrap()
    );
    assert_eq!(
        lower.field("controls").unwrap(),
        &parse("[[10,12,1],[26,24,2],[17,13,1],[66,39,3],[52,24,2],[30,12,1]]").unwrap()
    );
    let (blocks, faces) = construct(&normalized).unwrap();
    assert!(!blocks.is_empty());
    assert!(faces
        .iter()
        .any(|(_, name)| name == "south-curved-pocket:floor"));
}

#[test]
fn edge_and_operation_permutations_seeds_names_and_homogeneous_gauges_are_invariant() {
    let original = fixture();
    let expected = normalize(&original).unwrap();
    let mut changed = match original {
        solidcraft_trivariate::json::Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    let mut edges = changed.remove("edges").unwrap().items().unwrap().to_vec();
    edges.reverse();
    for edge in &mut edges {
        if let solidcraft_trivariate::json::Json::Object(fields) = edge {
            if let solidcraft_trivariate::json::Json::Object(curve) =
                fields.get_mut("curve").unwrap()
            {
                if let solidcraft_trivariate::json::Json::Array(controls) =
                    curve.get_mut("controls").unwrap()
                {
                    for point in controls {
                        if let solidcraft_trivariate::json::Json::Array(values) = point {
                            for value in values {
                                let n: i64 = value.number().unwrap();
                                *value =
                                    solidcraft_trivariate::json::Json::Number((n * 9).to_string());
                            }
                        }
                    }
                }
            }
        }
    }
    changed.insert(
        "edges".into(),
        solidcraft_trivariate::json::Json::Array(edges),
    );
    let mut cuts = changed.remove("cuts").unwrap().items().unwrap().to_vec();
    cuts.reverse();
    for cut in &mut cuts {
        if let solidcraft_trivariate::json::Json::Object(fields) = cut {
            let seed = fields.get("seed_edge").unwrap().text().unwrap();
            let replacement = match seed {
                "s_left" => "s_upper",
                "n_upper_reversed" => "n_lower",
                other => other,
            };
            fields.insert(
                "seed_edge".into(),
                solidcraft_trivariate::json::Json::String(replacement.into()),
            );
        }
    }
    changed.insert(
        "cuts".into(),
        solidcraft_trivariate::json::Json::Array(cuts),
    );
    assert_eq!(
        normalize(&solidcraft_trivariate::json::Json::Object(changed)).unwrap(),
        expected
    );
}

#[test]
fn rejects_open_branched_oblique_and_duplicate_component_operations() {
    let original = fixture();

    let mutate =
        |f: &dyn Fn(&mut std::collections::BTreeMap<String, solidcraft_trivariate::json::Json>)| {
            let mut fields = match original.clone() {
                solidcraft_trivariate::json::Json::Object(fields) => fields,
                _ => unreachable!(),
            };
            f(&mut fields);
            solidcraft_trivariate::json::Json::Object(fields)
        };

    let open = mutate(&|fields| {
        let mut edges = fields.get("edges").unwrap().items().unwrap().to_vec();
        edges.retain(|edge| edge.field("id").unwrap().text().unwrap() != "n_left");
        fields.insert(
            "edges".into(),
            solidcraft_trivariate::json::Json::Array(edges),
        );
    });
    assert!(normalize(&open).is_err());

    let branched = mutate(&|fields| {
        let mut edges = fields.get("edges").unwrap().items().unwrap().to_vec();
        let mut branch = edges[0].clone();
        if let solidcraft_trivariate::json::Json::Object(fields) = &mut branch {
            fields.insert(
                "id".into(),
                solidcraft_trivariate::json::Json::String("branch".into()),
            );
        }
        edges.push(branch);
        fields.insert(
            "edges".into(),
            solidcraft_trivariate::json::Json::Array(edges),
        );
    });
    assert!(normalize(&branched).is_err());

    let oblique = mutate(&|fields| {
        let mut edges = fields.get("edges").unwrap().items().unwrap().to_vec();
        for edge in &mut edges {
            if edge.field("id").unwrap().text().unwrap() == "s_right" {
                if let solidcraft_trivariate::json::Json::Object(edge_fields) = edge {
                    if let solidcraft_trivariate::json::Json::Object(curve) =
                        edge_fields.get_mut("curve").unwrap()
                    {
                        if let solidcraft_trivariate::json::Json::Array(controls) =
                            curve.get_mut("controls").unwrap()
                        {
                            if let solidcraft_trivariate::json::Json::Array(last) = &mut controls[1]
                            {
                                last[0] = solidcraft_trivariate::json::Json::Number("81".into());
                            }
                        }
                    }
                }
            }
        }
        fields.insert(
            "edges".into(),
            solidcraft_trivariate::json::Json::Array(edges),
        );
    });
    assert!(normalize(&oblique).is_err());

    let duplicate = mutate(&|fields| {
        let mut cuts = fields.get("cuts").unwrap().items().unwrap().to_vec();
        if let solidcraft_trivariate::json::Json::Object(first) = &mut cuts[0] {
            first.insert(
                "seed_edge".into(),
                solidcraft_trivariate::json::Json::String("n_lower".into()),
            );
        }
        fields.insert(
            "cuts".into(),
            solidcraft_trivariate::json::Json::Array(cuts),
        );
    });
    assert!(normalize(&duplicate).is_err());
}
