use solidcraft_trivariate::{
    circular_polygon::CircularPolygonBore,
    interval_bores::{build_interval_bore_prism, IntervalBore, IntervalBorePrism},
    interval_bores_session::IntervalBoreSession,
    perforated_prism::LoopProfile,
    polygon_prism::Vertex,
};

#[test]
fn failed_cap_binding_preserves_the_published_interval_and_cached_blocks() {
    let recipe = IntervalBorePrism {
        outer: LoopProfile {
            id: "profile".into(),
            vertices: [[-12, -12], [12, -12], [12, 12], [-12, 12]]
                .into_iter()
                .enumerate()
                .map(|(k, xy_mm)| Vertex {
                    id: format!("p{k}"),
                    xy_mm,
                })
                .collect(),
        },
        bores: vec![IntervalBore {
            bore: CircularPolygonBore {
                id: "pilot".into(),
                center_mm: [1, 1],
                radius_mm: 2,
                guard_mm: 4,
            },
            z_interval_mm: [2, 6],
        }],
        height_mm: 8,
        top_shift_mm: [3, -1],
    };
    let mut session = IntervalBoreSession::new(recipe.clone()).unwrap();
    let actual = session.current().blocks.as_ptr();
    assert!(session
        .apply("pilot", [0, 6], &["bore/pilot/lower".into()])
        .is_err());
    assert_eq!(session.current().blocks.as_ptr(), actual);
    assert_eq!(session.current().recipe, recipe);
    let no_change = session.apply("pilot", [2, 6], &[]).unwrap();
    assert_eq!((no_change.built_layers, no_change.reused_layers), (0, 3));
    assert_eq!(no_change.reused_blocks, session.current().blocks.len());
    let moved = session
        .apply("pilot", [0, 6], &["bore/pilot/wall".into()])
        .unwrap();
    assert_eq!(moved.old_class, "enclosed");
    assert_eq!(moved.new_class, "open_start");
    assert_eq!(moved.lost_features, ["bore/pilot/lower"]);
    assert_eq!(moved.volume_change_pi_mm3, -8);
    assert_eq!(
        session.current().blocks,
        build_interval_bore_prism(&session.current().recipe).unwrap()
    );
}
