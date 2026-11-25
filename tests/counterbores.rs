use solidcraft_trivariate::{
    counterbores::{build_counterbores, Counterbore, CounterboreOpening, CounterborePrism},
    perforated_prism::LoopProfile,
    polygon_prism::Vertex,
};
use std::collections::BTreeSet;

fn outer() -> LoopProfile {
    LoopProfile {
        id: "outer".into(),
        vertices: [[-24, -16], [24, -16], [24, 16], [-24, 16]]
            .into_iter()
            .enumerate()
            .map(|(i, xy_mm)| Vertex {
                id: format!("v{i}"),
                xy_mm,
            })
            .collect(),
    }
}

fn bore(id: &str, center_mm: [i32; 2], opening: CounterboreOpening, depth_mm: i32) -> Counterbore {
    Counterbore {
        id: id.into(),
        center_mm,
        radius_mm: 3,
        guard_mm: 12,
        counterbore_radius_mm: 6,
        depth_mm,
        opening,
    }
}

fn small_bore(
    id: &str,
    center_mm: [i32; 2],
    opening: CounterboreOpening,
    depth_mm: i32,
) -> Counterbore {
    Counterbore {
        id: id.into(),
        center_mm,
        radius_mm: 2,
        guard_mm: 6,
        counterbore_radius_mm: 3,
        depth_mm,
        opening,
    }
}

fn recipe(bores: Vec<Counterbore>) -> CounterborePrism {
    CounterborePrism {
        outer: outer(),
        bores,
        height_mm: 12,
        top_shift_mm: [3, -2],
    }
}

fn pilot_layers(blocks: &[solidcraft_trivariate::bracket::BracketBlock]) -> Vec<String> {
    blocks
        .iter()
        .filter(|block| block.id.contains("_pilot_sector_"))
        .map(|block| block.id.clone())
        .collect()
}

#[test]
fn start_and_end_counterbores_create_the_correct_pilot_only_layers() {
    let start = build_counterbores(&recipe(vec![bore(
        "start",
        [0, 0],
        CounterboreOpening::Start,
        4,
    )]))
    .unwrap();
    let start_pilot = pilot_layers(&start);
    assert_eq!(start_pilot.len(), 8);
    assert!(start_pilot.iter().all(|id| id.ends_with("_z4_12")));
    assert!(
        start
            .iter()
            .filter(|block| block.id.contains("_counter_sector_"))
            .count()
            == 16
    );

    let end = build_counterbores(&recipe(vec![bore(
        "end",
        [0, 0],
        CounterboreOpening::End,
        5,
    )]))
    .unwrap();
    let end_pilot = pilot_layers(&end);
    assert_eq!(end_pilot.len(), 8);
    assert!(end_pilot.iter().all(|id| id.ends_with("_z0_7")));
    assert!(end
        .iter()
        .all(|block| block.controls.iter().all(|control| control[3] > 0)));
}

#[test]
fn global_shoulder_events_split_all_g21_shell_cells() {
    let blocks = build_counterbores(&recipe(vec![
        Counterbore {
            id: "a".into(),
            center_mm: [-8, 0],
            radius_mm: 2,
            guard_mm: 6,
            counterbore_radius_mm: 3,
            depth_mm: 3,
            opening: CounterboreOpening::Start,
        },
        Counterbore {
            id: "b".into(),
            center_mm: [8, 0],
            radius_mm: 2,
            guard_mm: 6,
            counterbore_radius_mm: 3,
            depth_mm: 5,
            opening: CounterboreOpening::End,
        },
    ]))
    .unwrap();
    assert!(blocks.windows(2).all(|pair| pair[0].id < pair[1].id));
    let shell_layers = blocks
        .iter()
        .filter(|block| !block.id.starts_with("bore_"))
        .filter_map(|block| block.id.rsplit_once("_z"))
        .map(|(_, interval)| interval.to_owned())
        .collect::<BTreeSet<_>>();
    assert!(shell_layers.contains("0_3"));
    assert!(shell_layers.contains("3_7"));
    assert!(shell_layers.contains("7_12"));
    assert_eq!(pilot_layers(&blocks).len(), 32);
}

#[test]
fn bore_order_is_stable_and_invalid_counterbores_are_rejected() {
    let a = small_bore("a", [-8, 0], CounterboreOpening::Start, 4);
    let b = small_bore("b", [8, 0], CounterboreOpening::End, 5);
    assert_eq!(
        build_counterbores(&recipe(vec![a.clone(), b.clone()])).unwrap(),
        build_counterbores(&recipe(vec![b, a.clone()])).unwrap()
    );

    for invalid in [
        Counterbore {
            counterbore_radius_mm: 2,
            ..a.clone()
        },
        Counterbore {
            guard_mm: 5,
            ..a.clone()
        },
        Counterbore {
            depth_mm: 0,
            ..a.clone()
        },
        Counterbore {
            depth_mm: 12,
            ..a.clone()
        },
    ] {
        assert!(build_counterbores(&recipe(vec![invalid])).is_err());
    }
}
