//! Compare local stretch diagnostics and continuous bounds down to a 1 µm ligament.
use solidcraft_trivariate::{
    exact_bore::build_exact_bore,
    exact_bore_audit::audit_exact_bore,
    json::{object, Json},
    quality::{at_jacobian, bound_block},
    RationalBezierMap,
};
fn main() -> Result<(), &'static str> {
    let mut rows = Vec::new();
    for clearance in [6000, 1000, 100, 10, 1] {
        let center = [6000 - clearance, 0];
        let model = build_exact_bore(center[0], 0, 4000).map_err(|_| "bore build failed")?;
        let cert = audit_exact_bore(center, 4000, &model.blocks)?;
        let mut upper = 0_f64;
        let mut worst = 0_f64;
        let mut min_scaled = 1_f64;
        for block in &model.blocks {
            upper =
                upper.max(bound_block(block, cert.determinant_lower_bound_mm3)?.condition_upper);
            let map = RationalBezierMap::new(
                block.degrees,
                block.controls.iter().map(|p| p.map(|x| x as f64)).collect(),
            )
            .map_err(|_| "invalid map")?;
            for i in 0..=16 {
                for j in 0..=8 {
                    let jac = map
                        .evaluate([i as f64 / 16., j as f64 / 8., 0.5])
                        .map_err(|_| "invalid sample")?;
                    let q = at_jacobian(jac.jacobian_columns)?;
                    worst = worst.max(q.condition_number);
                    min_scaled = min_scaled.min(q.scaled_determinant);
                }
            }
        }
        let n = |x: f64| Json::Number(x.to_string());
        rows.push(object([
            ("clearance_um", n(clearance as f64)),
            ("determinant_lower_mm3", n(cert.determinant_lower_bound_mm3)),
            ("sample_maximum_condition", n(worst)),
            ("full_domain_condition_upper", n(upper)),
            ("sample_minimum_scaled_determinant", n(min_scaled)),
        ]));
    }
    println!("{}", Json::Array(rows));
    Ok(())
}
