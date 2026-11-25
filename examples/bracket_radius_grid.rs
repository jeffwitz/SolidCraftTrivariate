//! Exhaust the discrete two-radius bracket domain and integrate its actual maps.
use solidcraft_trivariate::{
    bracket_parameters::build_bracket_parameters,
    bracket_parameters_audit::{audit_bracket_parameter_derivatives, audit_bracket_parameters},
    elasticity::gauss,
    json::{object, Json},
    RationalBezierMap,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    time::Instant,
};

fn n(v: f64) -> Json {
    Json::Number(v.to_string())
}
fn main() -> Result<(), &'static str> {
    let start = Instant::now();
    let q = gauss(8)?;
    let mut cache = BTreeMap::new();
    let mut max_rad = 0_f64;
    let mut max_vol = 0_f64;
    let mut min_det = f64::INFINITY;
    let mut vol_witness = [0; 2];
    let mut rad_witness = [0; 2];
    let mut det_witness = [0; 2];
    let mut count = 0_u32;
    for rx in 40..=80 {
        for ry in 40..=80 {
            let model = build_bracket_parameters(rx, ry)?;
            let cert = audit_bracket_parameters(rx, ry, &model.blocks)?;
            audit_bracket_parameter_derivatives(rx, ry, &model.blocks, &model.derivatives)?;
            let boundary = &cert.certificate.design_faces;
            let tags = boundary.iter().map(|(_, tag)| tag).collect::<BTreeSet<_>>();
            if model.blocks.len() != 36
                || cert.certificate.incidence.interfaces.len() != 55
                || boundary.len() != 106
                || tags.len() != 11
            {
                return Err("unexpected bracket topology or selections");
            }
            let det = cert.certificate.determinant_lower_bound_mm3;
            if det < min_det {
                min_det = det;
                det_witness = [rx, ry];
            }
            let mut sum = 0.;
            for b in &model.blocks {
                let key = (b.degrees, b.controls.clone());
                let map = RationalBezierMap::new(
                    b.degrees,
                    b.controls.iter().map(|v| v.map(|x| x as f64)).collect(),
                )
                .map_err(|_| "invalid map")?;
                let volume = if let Some(&v) = cache.get(&key) {
                    v
                } else {
                    let mut v = 0.;
                    for &(u, wu) in &q {
                        for &(t, wt) in &q {
                            for &(w, ww) in &q {
                                let [a, b, c] = map
                                    .evaluate([u, t, w])
                                    .map_err(|_| "invalid integration point")?
                                    .jacobian_columns;
                                let det: f64 = (0..3)
                                    .map(|i| {
                                        a[i] * (b[(i + 1) % 3] * c[(i + 2) % 3]
                                            - b[(i + 2) % 3] * c[(i + 1) % 3])
                                    })
                                    .sum();
                                if det <= 0. {
                                    return Err("sampled nonpositive Jacobian");
                                }
                                v += det * wu * wt * ww;
                            }
                        }
                    }
                    cache.insert(key, v);
                    v
                };
                sum += volume;
                for (side, tag) in boundary.iter().filter(|(s, tag)| {
                    s.block == b.id && matches!(tag.as_str(), "hole_x" | "hole_y")
                }) {
                    if side.axis != 1 || !side.upper {
                        return Err("unexpected bore trace parameter convention");
                    }
                    for j in 0..9 {
                        let [x, y, z] = map
                            .evaluate([f64::from(j) / 8., 1., 0.375])
                            .map_err(|_| "invalid bore point")?
                            .position;
                        let err = if tag == "hole_x" {
                            ((x - 40.).hypot(z - 20.) - f64::from(rx) / 16.).abs()
                        } else {
                            ((y - 32.).hypot(z - 20.) - f64::from(ry) / 16.).abs()
                        };
                        if err > max_rad {
                            max_rad = err;
                            rad_witness = [rx, ry];
                        }
                    }
                }
            }
            let ideal = 34080.
                - std::f64::consts::PI
                    * (360. + 8. * ((f64::from(rx) / 16.).powi(2) + (f64::from(ry) / 16.).powi(2)));
            let error = (sum / ideal - 1.).abs();
            if error > max_vol {
                max_vol = error;
                vol_witness = [rx, ry];
            }
            count += 1;
        }
    }
    let invalid = [[39, 64], [81, 64], [64, 39], [64, 81], [0, 0], [255, 255]];
    for [rx, ry] in invalid {
        if build_bracket_parameters(rx, ry).is_ok() {
            return Err("out-of-range pair accepted");
        }
    }
    if max_vol > 1e-8 || max_rad > 1e-12 || min_det <= 0. {
        return Err("exhaustive radius grid tolerance missed");
    }
    let pair = |ticks: [u8; 2]| Json::Array(ticks.into_iter().map(|x| n(f64::from(x))).collect());
    println!(
        "{}",
        object([
            ("accepted_pairs", n(f64::from(count))),
            ("refused_out_of_range_pairs", n(invalid.len() as f64)),
            (
                "distinct_complete_control_nets_integrated_q8",
                n(cache.len() as f64)
            ),
            ("bore_test_points_per_pair", n(72.)),
            ("max_bore_radius_error_mm", n(max_rad)),
            ("radius_witness_ticks", pair(rad_witness)),
            ("max_relative_volume_error", n(max_vol)),
            ("volume_witness_ticks", pair(vol_witness)),
            ("minimum_certified_determinant_mm3", n(min_det)),
            ("determinant_witness_ticks", pair(det_witness)),
            ("elapsed_s", n(start.elapsed().as_secs_f64())),
        ])
    );
    Ok(())
}
