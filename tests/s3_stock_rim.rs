use solidcraft_trivariate::json::parse;
use solidcraft_trivariate::s3_stock_rim::{construct, RimState};

#[test]
fn all_stock_rim_states_keep_the_delivered_pocket_interfaces() {
    let source = parse(include_str!("../docs/results/s3-face-pocket.model.json")).unwrap();
    for state in [RimState::Sharp, RimState::Chamfer, RimState::Fillet] {
        let out = construct(&source, state).unwrap();
        assert_eq!(out.field("topology"), source.field("topology"));
        let vols = out.field("volumes").unwrap().items().unwrap();
        assert_eq!(vols.len(), 93);
        assert_eq!(
            vols.iter()
                .map(|v| v.field("spans").unwrap().items().unwrap().len())
                .sum::<usize>(),
            145
        );
    }
}

#[test]
fn unsupported_rim_request_preserves_the_destination() {
    let dest = std::env::temp_dir().join(format!(
        "solidcraft-rim-rejection-{}.json",
        std::process::id()
    ));
    std::fs::write(&dest, b"prior model").unwrap();
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
        .args([
            "edit-s3-stock-rim",
            "docs/results/s3-face-pocket.model.json",
            "unsupported",
        ])
        .arg(&dest)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(&dest).unwrap(), b"prior model");
    std::fs::remove_file(dest).unwrap();
}
