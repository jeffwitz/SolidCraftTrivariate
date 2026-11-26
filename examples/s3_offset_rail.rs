use solidcraft_trivariate::s3_offset_rail;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let side = std::env::args().nth(1).unwrap_or_else(|| "both".to_owned());
    let sides: Vec<i8> = match side.as_str() {
        "left" => vec![-1],
        "right" => vec![1],
        "both" => vec![-1, 1],
        _ => return Err("expected left, right, or both".into()),
    };
    let mut json=String::from("{\"schema\":\"s3-offset-rail-v1\",\"certificate\":\"interval-fourth-derivative-cubic-hermite-bound\",\"rails\":[");
    for (ri, s) in sides.iter().enumerate() {
        if ri > 0 {
            json.push(',')
        }
        let r = s3_offset_rail::build(*s)?;
        json.push_str(&format!("{{\"source_block\":\"{}\",\"side\":{},\"radius_mm\":{},\"continuity\":\"{}\",\"max_error_bound_mm\":{:.17e},\"spans\":[",r.source_block,r.side,r.radius_mm,r.continuity,r.spans.iter().map(|x|x.error_bound_mm).fold(0.0,f64::max)));
        for (i, sp) in r.spans.iter().enumerate() {
            if i > 0 {
                json.push(',')
            }
            json.push_str(&format!("{{\"interval\":[{:.17e},{:.17e}],\"degree\":3,\"coefficient_error_bound_mm\":{:.17e},\"json_roundtrip_error_bound_mm\":0,\"controls\":[",sp.interval[0],sp.interval[1],sp.coefficient_error_bound_mm));
            for (j, p) in sp.controls.iter().enumerate() {
                if j > 0 {
                    json.push(',')
                }
                json.push_str(&format!("[{:.17e},{:.17e},1]", p[0], p[1]));
            }
            json.push_str(&format!(
                "],\"error_bound_mm\":{:.17e}}}",
                sp.error_bound_mm
            ));
        }
        json.push_str("]}");
    }
    json.push_str("]}\n");
    let out = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "build/s3-offset-rail.json".to_owned());
    if let Some(parent) = std::path::Path::new(&out).parent() {
        std::fs::create_dir_all(parent)?;
    }
    std::fs::write(out, json)?;
    Ok(())
}
