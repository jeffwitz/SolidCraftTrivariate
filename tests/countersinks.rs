use solidcraft_trivariate::{
    counterbores::{Counterbore, CounterboreOpening, CounterborePrism},
    countersinks::build_countersinks,
    perforated_prism::LoopProfile,
    polygon_prism::Vertex,
};

fn recipe(opening: CounterboreOpening) -> CounterborePrism {
    CounterborePrism {
        outer: LoopProfile {
            id: "outer".into(),
            vertices: [[-20, -20], [20, -20], [20, 20], [-20, 20]]
                .into_iter()
                .enumerate()
                .map(|(i, xy_mm)| Vertex {
                    id: format!("v{i}"),
                    xy_mm,
                })
                .collect(),
        },
        bores: vec![Counterbore {
            id: "pilot".into(),
            center_mm: [0, 0],
            radius_mm: 2,
            guard_mm: 8,
            counterbore_radius_mm: 4,
            depth_mm: 3,
            opening,
        }],
        height_mm: 8,
        top_shift_mm: [0, 0],
    }
}

fn check_radius_at_layer(opening: CounterboreOpening, layer: &str, z: i32, expected: f64) {
    let blocks = build_countersinks(&recipe(opening)).unwrap();
    let ring = blocks
        .iter()
        .find(|b| b.id == format!("bore_pilot_sector_east_half_0_{layer}"))
        .unwrap();
    // The v=1, u=0 endpoint is a cardinal point on the circle. w is the
    // layer endpoint selected by the suffix; its homogeneous radius is exact.
    let w = if ring.id.ends_with(&format!("_{z}")) {
        1
    } else {
        0
    };
    let p = ring.controls[2 + w];
    let radius = (p[0] as f64 / p[3] as f64).hypot(p[1] as f64 / p[3] as f64);
    assert!(
        (radius - expected).abs() < 1e-12,
        "{} at z={z}: {radius} != {expected}",
        ring.id
    );
}

#[test]
fn start_and_end_tapers_retain_g21_guards_and_exact_layer_radii() {
    let start = build_countersinks(&recipe(CounterboreOpening::Start)).unwrap();
    let end = build_countersinks(&recipe(CounterboreOpening::End)).unwrap();
    assert_eq!(start.len(), end.len());
    assert!(start
        .iter()
        .any(|b| b.id == "bore_pilot_sector_east_half_0_z0_3"));
    assert!(start
        .iter()
        .any(|b| b.id == "bore_pilot_sector_east_half_0_z3_8"));
    assert!(start.iter().all(|b| b.controls.iter().all(|p| p[3] > 0)));
    check_radius_at_layer(CounterboreOpening::Start, "z0_3", 0, 4.);
    check_radius_at_layer(CounterboreOpening::Start, "z0_3", 3, 2.);
    check_radius_at_layer(CounterboreOpening::End, "z5_8", 5, 2.);
    check_radius_at_layer(CounterboreOpening::End, "z5_8", 8, 4.);
}

#[test]
fn invalid_taper_radii_and_depths_are_rejected() {
    let mut bad = recipe(CounterboreOpening::Start);
    bad.bores[0].counterbore_radius_mm = 2;
    assert!(build_countersinks(&bad).is_err());
    let mut bad = recipe(CounterboreOpening::End);
    bad.bores[0].depth_mm = 8;
    assert!(build_countersinks(&bad).is_err());
}
