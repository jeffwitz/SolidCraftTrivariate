//! Stable byte presentation for archived native annulus example models.
use solidcraft_trivariate::json::Json;

fn int_array<const N: usize>(value: &Json) -> [i64; N] {
    value
        .items()
        .expect("integer array")
        .iter()
        .map(|item| item.number::<i64>().expect("integer entry"))
        .collect::<Vec<_>>()
        .try_into()
        .expect("expected fixed-size integer array")
}

fn controls(value: &Json) -> Vec<[i64; 4]> {
    value
        .items()
        .expect("control net")
        .iter()
        .map(int_array::<4>)
        .collect()
}

pub fn print_model(model: &Json) {
    if let Ok(blocks) = model.field("blocks") {
        println!(
            r#"{{"degrees":[2,1,1],"knots":[[0,0,0,1,1,1],[0,0,1,1],[0,0,0.25,7.75,8,8]],"fastest_index":"axial","blocks":["#
        );
        let blocks = blocks.items().expect("block list");
        for (index, block) in blocks.iter().enumerate() {
            let id = block
                .field("id")
                .expect("block id")
                .text()
                .expect("string id");
            let net = controls(block.field("xyzw").expect("block net"));
            let comma = if index + 1 == blocks.len() { "" } else { "," };
            println!(r#"{{"id":"{id}","xyzw":{net:?}}}{comma}"#);
        }
        println!("]}}");
        return;
    }

    let charts = model
        .field("charts")
        .expect("chart list")
        .items()
        .expect("charts");
    println!(r#"{{"unit":"mm","fastest_index":"transverse","charts":["#);
    for (index, chart) in charts.iter().enumerate() {
        let id = chart
            .field("id")
            .expect("chart id")
            .text()
            .expect("string id");
        let degrees = int_array::<3>(chart.field("degrees").expect("degrees"));
        let shape = int_array::<3>(chart.field("shape").expect("shape"));
        let knots = chart.field("knots").expect("knots");
        let knots: Vec<Vec<i64>> = knots
            .items()
            .expect("knot axes")
            .iter()
            .map(|axis| {
                axis.items()
                    .expect("knot vector")
                    .iter()
                    .map(|item| item.number::<i64>().expect("integer knot"))
                    .collect()
            })
            .collect();
        let net = controls(chart.field("xyzw").expect("control net"));
        let comma = if index + 1 == charts.len() { "" } else { "," };
        println!(
            r#"{{"id":"{id}","degrees":{degrees:?},"shape":{shape:?},"knots":{knots:?},"xyzw":{net:?}}}{comma}"#
        );
    }
    println!("]}}");
}
