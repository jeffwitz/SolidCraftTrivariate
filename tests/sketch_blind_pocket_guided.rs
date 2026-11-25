use solidcraft_trivariate::json::{parse, Json};
use std::fs;
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT_DIR: AtomicUsize = AtomicUsize::new(0);

fn build_and_audit(recipe: &str) -> (Json, Json) {
    let dir = std::env::temp_dir().join(format!(
        "solidcraft-guided-pocket-{}-{}",
        std::process::id(),
        NEXT_DIR.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&dir).unwrap();
    let input = dir.join("recipe.json");
    let output = dir.join("model.json");
    fs::write(&input, recipe).unwrap();
    let built = Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
        .args(["build", input.to_str().unwrap(), output.to_str().unwrap()])
        .output()
        .expect("launch CLI build");
    assert!(
        built.status.success(),
        "{}",
        String::from_utf8_lossy(&built.stdout)
    );
    let model = parse(&fs::read_to_string(&output).unwrap()).unwrap();
    let audited = Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
        .args(["audit", output.to_str().unwrap()])
        .output()
        .expect("launch CLI audit");
    assert!(
        audited.status.success(),
        "{}",
        String::from_utf8_lossy(&audited.stdout)
    );
    let report = parse(std::str::from_utf8(&audited.stdout).unwrap()).unwrap();
    fs::remove_dir_all(dir).unwrap();
    (model, report)
}

fn recipe(outer: [[i32; 2]; 4], pocket: [[i32; 2]; 4], z0: i32) -> String {
    let points = |p: [[i32; 2]; 4]| {
        p.iter()
            .map(|v| format!("[{},{}]", v[0], v[1]))
            .collect::<Vec<_>>()
            .join(",")
    };
    format!(
        r#"{{"family":"sketch_layer_partition_v1","height_mm":20,"outer":{{"id":"body","vertices":[{}]}},"cuts":[{{"id":"pocket","vertices":[{}],"z_interval_mm":[{},20]}}],"construction":[]}}"#,
        points(outer),
        points(pocket),
        z0
    )
}

fn block_ids(model: &Json) -> Vec<String> {
    model
        .field("blocks")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|block| block.field("id").unwrap().text().unwrap().to_owned())
        .collect()
}

#[test]
fn blind_pocket_uses_nine_certified_blocks_and_records_the_baseline_comparison() {
    let canonical_recipe = recipe(
        [[0, 0], [60, 0], [60, 40], [0, 40]],
        [[20, 15], [40, 15], [40, 25], [20, 25]],
        14,
    );
    let (model, report) = build_and_audit(&canonical_recipe);
    assert_eq!(model.field("blocks").unwrap().items().unwrap().len(), 9);
    assert_eq!(report.field("status").unwrap().text().unwrap(), "valid");
    assert_eq!(
        report.field("blocks").unwrap().number::<usize>().unwrap(),
        9
    );
    let plan = model.field("partition_plan").unwrap();
    assert_eq!(
        plan.field("method").unwrap().text().unwrap(),
        "guided_rectangular_blind_pocket"
    );
    let comparison = plan.field("candidate_comparison").unwrap();
    let baseline = comparison.field("baseline").unwrap();
    let guided = comparison.field("guided").unwrap();
    assert_eq!(
        baseline.field("valid_certificate").unwrap(),
        &Json::Bool(true)
    );
    assert_eq!(
        baseline.field("blocks").unwrap().number::<usize>().unwrap(),
        17
    );
    assert_eq!(
        guided.field("valid_certificate").unwrap(),
        &Json::Bool(true)
    );
    assert_eq!(
        guided.field("blocks").unwrap().number::<usize>().unwrap(),
        9
    );
    let layers: Vec<_> = plan
        .field("layers")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|entry| {
            entry
                .field("material_blocks")
                .unwrap()
                .number::<u64>()
                .unwrap()
        })
        .collect();
    assert_eq!(layers, [5, 4]);
    let guide_edges = plan
        .field("edges")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .filter(|edge| {
            edge.field("sources")
                .unwrap()
                .items()
                .unwrap()
                .iter()
                .any(|source| {
                    source
                        .field("id")
                        .unwrap()
                        .text()
                        .is_ok_and(|id| id.contains("pocket-corner-"))
                })
        })
        .count();
    assert_eq!(guide_edges, 4);
    let volume = report.field("volume_exact_mm3").unwrap().items().unwrap();
    assert_eq!(volume[0].number::<i128>().unwrap(), 46800);
    assert_eq!(volume[1].number::<i128>().unwrap(), 1);

    let (rewritten, _) = build_and_audit(&recipe(
        [[60, 40], [0, 40], [0, 0], [60, 0]],
        [[40, 25], [20, 25], [20, 15], [40, 15]],
        14,
    ));
    assert_eq!(block_ids(&model), block_ids(&rewritten));
}

#[test]
fn shifted_and_resized_rectangular_blind_pockets_keep_the_nine_block_topology() {
    for recipe in [
        recipe(
            [[-10, 3], [70, 3], [70, 53], [-10, 53]],
            [[12, 17], [42, 17], [42, 29], [12, 29]],
            13,
        ),
        recipe(
            [[100, -30], [160, -30], [160, 10], [100, 10]],
            [[124, -7], [134, -7], [134, 1], [124, 1]],
            8,
        ),
    ] {
        let (model, report) = build_and_audit(&recipe);
        assert_eq!(model.field("blocks").unwrap().items().unwrap().len(), 9);
        assert_eq!(report.field("status").unwrap().text().unwrap(), "valid");
        assert_eq!(
            report.field("blocks").unwrap().number::<usize>().unwrap(),
            9
        );
        assert_eq!(
            model
                .field("partition_plan")
                .unwrap()
                .field("candidate_comparison")
                .unwrap()
                .field("baseline")
                .unwrap()
                .field("blocks")
                .unwrap()
                .number::<usize>()
                .unwrap(),
            17
        );
    }
}
