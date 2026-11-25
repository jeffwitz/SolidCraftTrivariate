use solidcraft_trivariate::{
    interchange::{audit, build},
    json::{parse, Json},
};

const POCKET: &str =
    include_str!("../docs/site/data/plate_blind_pocket/plate_blind_pocket.recipe.json");
const C_SECTION: &str =
    include_str!("../docs/site/data/c_section_two_holes/c_section_two_holes.recipe.json");

#[test]
fn curved_pocket_preserves_its_nine_blocks_and_rejects_altered_controls() {
    let mut model = build(&parse(POCKET).unwrap()).unwrap();
    let report = audit(&model).unwrap();
    assert_eq!(
        report.field("blocks").unwrap().number::<usize>().unwrap(),
        9
    );
    assert!(
        report
            .field("determinant_lower_bound_mm3")
            .unwrap()
            .number::<f64>()
            .unwrap()
            > 0.
    );
    let Json::Object(fields) = &mut model else {
        unreachable!()
    };
    let Json::Array(blocks) = fields.get_mut("blocks").unwrap() else {
        unreachable!()
    };
    let Json::Object(first) = &mut blocks[0] else {
        unreachable!()
    };
    let Json::Array(net) = first.get_mut("controls").unwrap() else {
        unreachable!()
    };
    let Json::Array(point) = &mut net[0] else {
        unreachable!()
    };
    point[0] = Json::Number((point[0].number::<i64>().unwrap() + 1).to_string());
    assert_eq!(audit(&model).unwrap_err(), "curved_sketch_replay_mismatch");
}

#[test]
fn swept_c_section_preserves_its_oblique_axis_and_uses_quadratic_wall_surfaces() {
    let model = build(&parse(C_SECTION).unwrap()).unwrap();
    let receipt = audit(&model).unwrap();
    assert_eq!(
        receipt
            .field("shared_faces")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        10
    );
    let blocks = model.field("blocks").unwrap().items().unwrap();
    for block in blocks {
        let net = block.field("controls").unwrap().items().unwrap();
        assert_eq!(net.len(), 18);
        for column in net.as_chunks::<2>().0 {
            let a = column[0].items().unwrap();
            let b = column[1].items().unwrap();
            let weight = a[3].number::<i64>().unwrap();
            for (i, shift) in [3, -2, 7].into_iter().enumerate() {
                assert_eq!(
                    b[i].number::<i64>().unwrap() - a[i].number::<i64>().unwrap(),
                    shift * weight
                );
            }
        }
    }
    let step = solidcraft_trivariate::step_boundary::export_step(&model).unwrap();
    assert!(step.contains("MANIFOLD_SOLID_BREP"));
}
