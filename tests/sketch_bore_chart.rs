use solidcraft_trivariate::{
    bracket::BracketBlock, incidence::enumerate_faces, planar::certify_planar_refined,
    sketch_bore_chart::bore_charts,
};

fn point(block: &solidcraft_trivariate::bracket::BracketBlock, index: usize) -> [f64; 3] {
    let p = block.controls[index];
    std::array::from_fn(|i| p[i] as f64 / p[3] as f64)
}

fn square_side_box(
    axis: usize,
    center: [i32; 3],
    interval: [i32; 2],
    half_width: i32,
    quadrant: usize,
    half: usize,
) -> BracketBlock {
    let g = i128::from(half_width);
    let (start, corner, end) = match quadrant {
        0 => ([g, 0], [g, g], [0, g]),
        1 => ([0, g], [-g, g], [-g, 0]),
        2 => ([-g, 0], [-g, -g], [0, -g]),
        _ => ([0, -g], [g, -g], [g, 0]),
    };
    let (a, b) = if half == 0 {
        (start, corner)
    } else {
        (corner, end)
    };
    let outward = [b[1] - a[1], a[0] - b[0]];
    let along = [
        [2 * b[0], 2 * b[1]],
        [a[0] + b[0], a[1] + b[1]],
        [2 * a[0], 2 * a[1]],
    ];
    let axes = [(axis + 1) % 3, (axis + 2) % 3];
    let mut controls = Vec::new();
    for u in along {
        for v in 0..=1 {
            for z in interval {
                let mut p = [0_i128; 4];
                p[axes[0]] =
                    i128::from(center[axes[0]]) * 2 + u[0] + i128::from(v) * 2 * outward[0];
                p[axes[1]] =
                    i128::from(center[axes[1]]) * 2 + u[1] + i128::from(v) * 2 * outward[1];
                p[axis] = i128::from(z) * 2;
                p[3] = 2;
                controls.push(p.map(|value| i64::try_from(value).unwrap()));
            }
        }
    }
    BracketBlock {
        id: format!("square_{axis}_{quadrant}_{half}"),
        degrees: [2, 1, 1],
        controls,
    }
}

#[test]
fn all_axes_have_eight_matching_exact_half_charts() {
    for axis in 0..3 {
        let center = [3, -4, 11];
        let interval = [-2, 7];
        let radius = 3;
        let blocks = bore_charts(axis, center, interval, 9, radius, "probe").unwrap();
        assert_eq!(blocks.len(), 8);
        assert!(blocks.iter().all(|b| b.degrees == [2, 1, 1]));
        assert!(blocks.iter().all(|b| b.controls.iter().all(|p| p[3] > 0)));
        let incidence = enumerate_faces(&blocks).unwrap();
        assert_eq!(incidence.interfaces.len(), 8);

        let transverse = [(axis + 1) % 3, (axis + 2) % 3];
        for block in &blocks {
            let planar = block
                .controls
                .iter()
                .step_by(2)
                .map(|p| [p[transverse[0]], p[transverse[1]], p[3]])
                .collect::<Vec<_>>();
            assert!(
                certify_planar_refined([2, 1], &planar).is_ok(),
                "{}",
                block.id
            );
            for u in [0, 2] {
                for (w, &z) in interval.iter().enumerate() {
                    let p = point(block, (u * 2 + 1) * 2 + w);
                    let dx = p[transverse[0]] - f64::from(center[transverse[0]]);
                    let dy = p[transverse[1]] - f64::from(center[transverse[1]]);
                    assert!((dx.hypot(dy) - f64::from(radius)).abs() < 1e-12);
                    assert_eq!(p[axis], f64::from(z));
                }
            }
        }
    }
}

#[test]
fn outer_half_sides_match_linear_box_faces_exactly() {
    let axis = 2;
    let center = [3, -4, 0];
    let interval = [1, 9];
    let mut blocks = bore_charts(axis, center, interval, 8, 3, "join").unwrap();
    for index in 0..8 {
        blocks.push(square_side_box(
            axis,
            center,
            interval,
            8,
            index / 2,
            index % 2,
        ));
    }
    let incidence = enumerate_faces(&blocks).unwrap();
    assert_eq!(incidence.interfaces.len(), 20);
    let chart_box_interfaces = incidence
        .interfaces
        .iter()
        .filter(|interface| {
            let first_chart = interface.first.block.starts_with("bore_join_sector_");
            let second_chart = interface.second.block.starts_with("bore_join_sector_");
            first_chart != second_chart
                && (interface.first.block.starts_with("square_")
                    || interface.second.block.starts_with("square_"))
        })
        .count();
    assert_eq!(chart_box_interfaces, 8);
}

#[test]
fn invalid_domains_and_ids_are_refused_and_i32_extremes_remain_representable() {
    assert!(bore_charts(3, [0; 3], [0, 1], 5, 2, "b").is_err());
    assert!(bore_charts(0, [0; 3], [1, 1], 5, 2, "b").is_err());
    assert!(bore_charts(0, [0; 3], [0, 1], 2, 2, "b").is_err());
    assert!(bore_charts(0, [0; 3], [0, 1], 5, 0, "b").is_err());
    assert!(bore_charts(0, [0; 3], [0, 1], 5, 2, "bad/id").is_err());
    let extreme = bore_charts(
        2,
        [i32::MAX, i32::MAX, i32::MAX],
        [0, 1],
        i32::MAX,
        i32::MAX - 1,
        "b",
    )
    .unwrap();
    assert_eq!(extreme.len(), 8);
}
