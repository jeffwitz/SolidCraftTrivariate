use solidcraft_trivariate::{
    bracket::BracketBlock,
    bracket_edge_edit::{
        chamfer_compact_bracket, edit_compact_bracket, round_fillet_fixed_partition_status,
        BoreChamfer, BoreEdge, EdgeTreatment,
    },
    interchange,
    json::{parse, Json},
};

fn blocks(model: &Json) -> Vec<BracketBlock> {
    model
        .field("blocks")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|b| BracketBlock {
            id: b.field("id").unwrap().text().unwrap().to_owned(),
            degrees: b
                .field("degrees")
                .unwrap()
                .items()
                .unwrap()
                .iter()
                .map(|v| v.number::<usize>().unwrap())
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),
            controls: b
                .field("controls")
                .unwrap()
                .items()
                .unwrap()
                .iter()
                .map(|p| {
                    p.items()
                        .unwrap()
                        .iter()
                        .map(|v| v.number::<i64>().unwrap())
                        .collect::<Vec<_>>()
                        .try_into()
                        .unwrap()
                })
                .collect(),
        })
        .collect()
}

fn de_casteljau(values: &[f64], t: f64) -> f64 {
    let mut work = values.to_vec();
    for n in (1..work.len()).rev() {
        for i in 0..n {
            work[i] = (1.0 - t) * work[i] + t * work[i + 1];
        }
    }
    work[0]
}

fn eval(block: &BracketBlock, q: [f64; 3]) -> [f64; 3] {
    let [du, dv, dw] = block.degrees;
    let mut h = [0.0; 4];
    for i in 0..=du {
        let mut along_v = vec![[0.0; 4]; dw + 1];
        for j in 0..=dv {
            for (k, accumulated) in along_v.iter_mut().enumerate() {
                let p = block.controls[(i * (dv + 1) + j) * (dw + 1) + k];
                for c in 0..4 {
                    accumulated[c] += p[c] as f64 * bernstein(dv, j, q[1]);
                }
            }
        }
        for c in 0..4 {
            let row = along_v.iter().map(|p| p[c]).collect::<Vec<_>>();
            let value = de_casteljau(&row, q[2]);
            h[c] += value * bernstein(du, i, q[0]);
        }
    }
    std::array::from_fn(|a| h[a] / h[3])
}

fn bernstein(n: usize, i: usize, t: f64) -> f64 {
    let binomial = match (n, i) {
        (0, 0) | (1, 0) | (1, 1) => 1.0,
        (2, 0) | (2, 2) => 1.0,
        (2, 1) => 2.0,
        _ => unreachable!(),
    };
    binomial * t.powi(i as i32) * (1.0 - t).powi((n - i) as i32)
}

fn fixture() -> Json {
    parse(include_str!(
        "../docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json"
    ))
    .unwrap()
}

fn bore_inputs() -> [BoreChamfer; 2] {
    [
        BoreChamfer {
            chart_prefix: "wing_east_",
            axis: 1,
            center_quarters: [160, 16, 80],
            radius_quarters: 16,
            width_quarters: 1,
        },
        BoreChamfer {
            chart_prefix: "wing_north_",
            axis: 0,
            center_quarters: [16, 128, 80],
            radius_quarters: 16,
            width_quarters: 1,
        },
    ]
}

#[test]
fn fieldcraft_r6_r4_bracket_keeps_ten_logical_charts_and_certifies_every_chamfer_span() {
    let model = fixture();
    interchange::audit(&model).expect("published source model must pass its independent audit");
    let source = blocks(&model);
    let edited = chamfer_compact_bracket(&source, bore_inputs()).unwrap();
    assert_eq!(edited.logical_blocks.len(), 10);
    assert_eq!(edited.elements.len(), 30);
    assert_eq!(edited.incidence.interfaces.len(), 53);
    assert!(edited.determinant_lower_bound_mm3 > 0.0);

    for bore in bore_inputs() {
        let charts = edited
            .elements
            .iter()
            .filter(|span| span.logical_id.starts_with(bore.chart_prefix))
            .collect::<Vec<_>>();
        assert_eq!(charts.len(), 12);
        let [lo, hi] = charts.iter().fold([i32::MAX, i32::MIN], |r, s| {
            [
                r[0].min(s.station_quarters[0]),
                r[1].max(s.station_quarters[1]),
            ]
        });
        assert_eq!(charts.iter().filter(|s| s.span_index == 0).count(), 4);
        for chart in charts {
            let block = BracketBlock {
                id: chart.logical_id.clone(),
                degrees: chart.degrees,
                controls: chart.controls.clone(),
            };
            let at_edge = if chart.station_quarters[0] == lo {
                Some((0.0, 17.0 / 4.0))
            } else if chart.station_quarters[1] == hi {
                Some((1.0, 17.0 / 4.0))
            } else {
                None
            };
            if let Some((v, radius)) = at_edge {
                for u in [0.0, 0.5, 1.0] {
                    let p = eval(&block, [u, 1.0, v]); // inner bore wall
                    let [a, b] = [(bore.axis + 1) % 3, (bore.axis + 2) % 3];
                    let center = [
                        bore.center_quarters[a] as f64 / 4.0,
                        bore.center_quarters[b] as f64 / 4.0,
                    ];
                    let observed = ((p[a] - center[0]).powi(2) + (p[b] - center[1]).powi(2)).sqrt();
                    assert!(
                        (observed - radius).abs() < 2e-9,
                        "{} at {u}: {observed}",
                        chart.logical_id
                    );
                }
            }
        }
    }
}

#[test]
fn invalid_width_and_wrong_parent_are_rejected() {
    let source = blocks(&fixture());
    let mut bores = bore_inputs();
    bores[0].width_quarters = 16;
    assert!(chamfer_compact_bracket(&source, bores).is_err());
    assert!(chamfer_compact_bracket(&source[..9], bore_inputs()).is_err());
}

#[test]
fn sharp_chamfer_and_fillet_share_the_same_thirty_span_partition() {
    let source = blocks(&fixture());
    let bores = [
        BoreEdge {
            chart_prefix: "wing_east_",
            axis: 1,
            center_quarters: [160, 16, 80],
            radius_quarters: 16,
        },
        BoreEdge {
            chart_prefix: "wing_north_",
            axis: 0,
            center_quarters: [16, 128, 80],
            radius_quarters: 16,
        },
    ];
    let sharp = edit_compact_bracket(&source, bores, EdgeTreatment::Sharp).unwrap();
    let chamfer =
        edit_compact_bracket(&source, bores, EdgeTreatment::Chamfer { width_quarters: 1 }).unwrap();
    assert_eq!(sharp.elements.len(), 30);
    assert_eq!(chamfer.elements.len(), 30);
    assert_eq!(sharp.incidence.interfaces.len(), 53);
    assert_eq!(chamfer.incidence.interfaces.len(), 53);
    assert!(sharp.determinant_lower_bound_mm3 > 0.0);
    assert!(chamfer.determinant_lower_bound_mm3 > 0.0);
    assert!(!sharp.singular_boundary_qualified && !chamfer.singular_boundary_qualified);
}

#[test]
fn round_fillet_domain_is_explicitly_scoped_to_matching_seam_requirement() {
    let source = blocks(&fixture());
    let bores = [
        BoreEdge {
            chart_prefix: "wing_east_",
            axis: 1,
            center_quarters: [160, 16, 80],
            radius_quarters: 16,
        },
        BoreEdge {
            chart_prefix: "wing_north_",
            axis: 0,
            center_quarters: [16, 128, 80],
            radius_quarters: 16,
        },
    ];
    let fillet = edit_compact_bracket(&source, bores, EdgeTreatment::Fillet { radius_quarters: 1 });
    let fillet = fillet.unwrap();
    assert!(fillet.singular_boundary_qualified);
    assert_eq!(fillet.determinant_lower_bound_mm3, 0.0);
    assert_eq!(fillet.elements.len(), 30);
    assert_eq!(fillet.incidence.interfaces.len(), 53);
    assert!(round_fillet_fixed_partition_status().contains("zero-Jacobian contact boundary"));
}
