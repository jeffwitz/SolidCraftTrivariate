use solidcraft_trivariate::{
    bracket::BracketBlock,
    incidence::BlockSide,
    interchange,
    json::{parse, Json},
    steinmetz_audit::{certify_boundary, certify_wall_features},
};

fn array<T: std::str::FromStr>(json: &Json) -> Result<Vec<T>, &'static str> {
    json.items()?
        .iter()
        .map(Json::number)
        .collect::<Result<_, _>>()
}

fn decode_block(value: &Json) -> Result<BracketBlock, &'static str> {
    Ok(BracketBlock {
        id: value.field("id")?.text()?.to_owned(),
        degrees: array(value.field("degrees")?)?
            .try_into()
            .map_err(|_| "expected three degrees")?,
        controls: value
            .field("controls")?
            .items()?
            .iter()
            .map(|control| {
                array(control)?
                    .try_into()
                    .map_err(|_| "expected four homogeneous integers")
            })
            .collect::<Result<_, _>>()?,
    })
}

fn side_from_name(block: &str, name: &str) -> Result<BlockSide, &'static str> {
    let (axis, upper) = match name {
        "u0" => (0, false),
        "u1" => (0, true),
        "v0" => (1, false),
        "v1" => (1, true),
        "w0" => (2, false),
        "w1" => (2, true),
        _ => return Err("unknown face side"),
    };
    Ok(BlockSide {
        block: block.to_owned(),
        axis,
        upper,
    })
}

fn flat(index: [usize; 3], degrees: [usize; 3]) -> usize {
    (index[0] * (degrees[1] + 1) + index[1]) * (degrees[2] + 1) + index[2]
}

#[test]
fn independently_audits_both_exact_wall_families_and_mutation() {
    let recipe = parse(include_str!("references/steinmetz_box.recipe.json")).unwrap();
    let mut model = interchange::build(&recipe).unwrap();
    let report = certify_wall_features(&model).unwrap();
    assert_eq!(report.total_faces, 32);
    assert_eq!(report.faces_by_bore["a"], 16);
    assert_eq!(report.faces_by_bore["b"], 16);

    let face = model
        .field("design_faces")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .find(|face| face.field("feature").unwrap().text().unwrap() == "bore/a/wall")
        .unwrap();
    let block_id = face.field("block").unwrap().text().unwrap().to_owned();
    let side = side_from_name(&block_id, face.field("side").unwrap().text().unwrap()).unwrap();
    let block_json = model
        .field("blocks")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .find(|block| block.field("id").unwrap().text().unwrap() == block_id)
        .unwrap();
    let block = decode_block(block_json).unwrap();
    assert!(certify_boundary(&block, &side, 0, [30, 20, 15], 3).unwrap());

    let free = match side.axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    };
    let mut local = [0; 3];
    local[side.axis] = if side.upper {
        block.degrees[side.axis]
    } else {
        0
    };
    let control_index = flat(local, block.degrees);
    let transverse_coordinate = 2;
    assert!(free.contains(&transverse_coordinate));
    let Json::Object(fields) = &mut model else {
        panic!("model must be an object")
    };
    let Json::Array(blocks) = fields.get_mut("blocks").unwrap() else {
        panic!("blocks must be an array")
    };
    let target = blocks
        .iter_mut()
        .find(|candidate| candidate.field("id").unwrap().text().unwrap() == block_id)
        .unwrap();
    let Json::Object(block_fields) = target else {
        panic!("block must be an object")
    };
    let Json::Array(controls) = block_fields.get_mut("controls").unwrap() else {
        panic!("controls must be an array")
    };
    let Json::Array(control) = &mut controls[control_index] else {
        panic!("control must be a four-vector")
    };
    let value: i64 = control[transverse_coordinate].number().unwrap();
    control[transverse_coordinate] = Json::Number((value + 1).to_string());

    let changed = decode_block(
        blocks
            .iter()
            .find(|candidate| candidate.field("id").unwrap().text().unwrap() == block_id)
            .unwrap(),
    )
    .unwrap();
    assert!(!certify_boundary(&changed, &side, 0, [30, 20, 15], 3).unwrap());
}
