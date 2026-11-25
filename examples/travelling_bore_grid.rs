//! Exhaust the integer-centre domain of the bounded travelling-bore family.
use solidcraft_trivariate::{
    bracket::BracketBlock,
    elasticity::gauss,
    json::{object, Json},
    travelling_bore_bounded::travelling_bore_bounded_blocks,
    travelling_bore_bounded_audit::audit_travelling_bore_bounded,
    travelling_bore_events::analyze_travelling_bore,
    RationalBezierMap,
};
use std::{collections::BTreeSet, time::Instant};

const TANGENCIES: [i32; 4] = [-4000, 0, 8000, 12000];
fn n(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}

fn volume(blocks: &[BracketBlock], q: &[(f64, f64)]) -> Result<f64, &'static str> {
    let mut volume = 0.;
    for b in blocks {
        let map = RationalBezierMap::new(
            b.degrees,
            b.controls.iter().map(|p| p.map(|v| v as f64)).collect(),
        )
        .map_err(|_| "invalid integration map")?;
        for &(u, wu) in q {
            for &(v, wv) in q {
                for &(w, ww) in q {
                    let [a, b, c] = map
                        .evaluate([u, v, w])
                        .map_err(|_| "invalid integration point")?
                        .jacobian_columns;
                    let det = a[0] * (b[1] * c[2] - b[2] * c[1])
                        - a[1] * (b[0] * c[2] - b[2] * c[0])
                        + a[2] * (b[0] * c[1] - b[1] * c[0]);
                    if !det.is_finite() || det <= 0. {
                        return Err("nonpositive evaluated determinant");
                    }
                    volume += det * wu * wv * ww;
                }
            }
        }
    }
    Ok(volume)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let mut rows = Vec::new();
    let mut failures = Vec::new();
    let mut refusals = Vec::new();
    let mut radial = (0, 0.);
    let mut hausdorff = (0, 0.);
    let mut det = (0, f64::INFINITY);
    let mut samples = BTreeSet::new();
    for center in -5000..=14000 {
        if center % 500 == 0 {
            eprintln!(
                "centre {center}, elapsed {:.1}s",
                start.elapsed().as_secs_f64()
            );
        }
        let now = Instant::now();
        let result = travelling_bore_bounded_blocks(center).and_then(|blocks| {
            let r = audit_travelling_bore_bounded(center, &blocks)?;
            Ok((blocks.len(), r))
        });
        match result {
            Ok((blocks, r)) => {
                if TANGENCIES.contains(&center) {
                    failures.push(object([
                        ("center_um", n(center)),
                        ("reason", Json::String("accepted tangent".into())),
                    ]));
                }
                if center % 100 == 0 || TANGENCIES.iter().any(|e| (center - e).abs() == 1) {
                    samples.insert(center);
                }
                let minimum = r.certificate.determinant_lower_bound_mm3;
                if radial.1 < r.radial_error_upper_mm {
                    radial = (center, r.radial_error_upper_mm);
                }
                if hausdorff.1 < r.hausdorff_error_upper_mm {
                    hausdorff = (center, r.hausdorff_error_upper_mm);
                }
                if minimum < det.1 {
                    det = (center, minimum);
                }
                rows.push(Json::Array(vec![
                    n(center),
                    n(blocks),
                    n(r.certificate.incidence.interfaces.len()),
                    n(r.certificate.incidence.exterior.len()),
                    n(minimum),
                    n(r.radial_error_upper_mm),
                    n(r.hausdorff_error_upper_mm),
                    n(now.elapsed().as_secs_f64()),
                ]));
            }
            Err(message) => {
                let record = object([
                    ("center_um", n(center)),
                    ("reason", Json::String(message.into())),
                ]);
                if TANGENCIES.contains(&center) && message.contains("tangent") {
                    refusals.push(record);
                } else {
                    failures.push(record);
                }
            }
        }
    }
    // Replace the rejected multiples of 100 by their admitted neighbours to
    // retain at least 191 volume controls distributed across the recipe range.
    samples.extend([radial.0, hausdorff.0, det.0]);
    let q = gauss(8)?;
    let mut integrals = Vec::new();
    let mut worst = (0, 0.);
    for center in samples {
        let blocks = travelling_bore_bounded_blocks(center)?;
        let measured = volume(&blocks, &q)?;
        let ideal = analyze_travelling_bore(center.into())?.volume_mm3;
        let error = (measured / ideal - 1.).abs();
        if error > worst.1 {
            worst = (center, error);
        }
        integrals.push(Json::Array(vec![
            n(center),
            n(measured),
            n(ideal),
            n(error),
        ]));
    }
    let witness = |(c, v): (i32, f64)| object([("center_um", n(c)), ("value", n(v))]);
    let qualified = failures.is_empty()
        && rows.len() == 18997
        && refusals.len() == 4
        && integrals.len() >= 191
        && worst.1 <= 1e-9;
    let report = object([
        (
            "schema",
            Json::String("travelling_bore_exhaustive_grid_v1".into()),
        ),
        ("domain_um", Json::Array(vec![n(-5000), n(14000)])),
        ("visits", n(19001)),
        ("admissions", n(rows.len())),
        ("expected_refusals", Json::Array(refusals)),
        ("unexpected_outcomes", Json::Array(failures)),
        ("max_radial_bound_mm", witness(radial)),
        ("max_hausdorff_bound_mm", witness(hausdorff)),
        ("minimum_determinant_bound_mm3", witness(det)),
        ("volume_samples", n(integrals.len())),
        ("maximum_relative_volume_error", witness(worst)),
        (
            "row_columns",
            Json::Array(
                [
                    "center_um",
                    "blocks",
                    "interfaces",
                    "exterior",
                    "determinant_lower_bound_mm3",
                    "radial_upper_mm",
                    "hausdorff_upper_mm",
                    "build_audit_s",
                ]
                .map(|s| Json::String(s.into()))
                .into(),
            ),
        ),
        (
            "volume_columns",
            Json::Array(
                ["center_um", "q8_mm3", "ideal_mm3", "relative_error"]
                    .map(|s| Json::String(s.into()))
                    .into(),
            ),
        ),
        ("rows", Json::Array(rows)),
        ("volume_checks", Json::Array(integrals)),
        ("wall_s", n(start.elapsed().as_secs_f64())),
        ("qualified", Json::Bool(qualified)),
    ]);
    println!("{report}");
    Ok(())
}
