//! Per-case adapter used only to compare independent P6/P7 references to the
//! public Rust APIs. It reads one JSON array from stdin and prints one array.
use solidcraft_trivariate::{
    json::{object, parse, Json},
    reparameterization::optimize_reparameterization,
    sketch_intersections::intersect_circles,
};
use std::io::{self, Read};

fn n(value: f64) -> Json {
    if value.is_finite() {
        Json::Number(value.to_string())
    } else {
        Json::Null
    }
}
fn u(value: usize) -> Json {
    Json::Number(value.to_string())
}
fn b(value: bool) -> Json {
    Json::Bool(value)
}
fn s(value: &str) -> Json {
    Json::String(value.to_owned())
}
fn array(values: impl IntoIterator<Item = Json>) -> Json {
    Json::Array(values.into_iter().collect())
}
fn f64_value(value: &Json) -> Result<f64, &'static str> {
    match value {
        Json::Number(v) => v.parse().map_err(|_| "invalid float"),
        Json::String(v) if v == "NaN" => Ok(f64::NAN),
        Json::String(v) if v == "Infinity" => Ok(f64::INFINITY),
        Json::String(v) if v == "-Infinity" => Ok(f64::NEG_INFINITY),
        _ => Err("expected number"),
    }
}
fn case_name(value: &Json) -> Result<&str, &'static str> {
    value.field("case")?.text()
}

fn p6(input: &Json) -> Result<Json, &'static str> {
    let theta = f64_value(input.field("theta")?)?;
    let weights = input.field("weights")?.items()?;
    if weights.len() != 2 {
        return Err("weights must contain two values");
    }
    let data_weights = [f64_value(&weights[0])?, f64_value(&weights[1])?];
    let regularization = f64_value(input.field("regularization")?)?;
    let result = optimize_reparameterization(theta, data_weights, regularization);
    let name = s(case_name(input)?);
    match result {
        Err(error) => Ok(object([
            ("case", name),
            ("accepted", b(false)),
            ("error", s(error)),
        ])),
        Ok(result) => {
            let coefficients = result
                .coefficients
                .map(|v| array(v.into_iter().map(n)))
                .unwrap_or(Json::Null);
            let derivative = result
                .coefficient_derivative_wrt_theta
                .map(|v| array(v.into_iter().map(n)))
                .unwrap_or(Json::Null);
            let speed = result.coefficients.map(|v| 1.0 - v[0].abs());
            let values = result
                .coefficients
                .map(|_| {
                    let vals = (0..2).flat_map(|brick| {
                        [0.125, 0.37, 0.5, 0.83]
                            .into_iter()
                            .map(move |u| n(result.map(brick, u).expect("valid probe point")))
                    });
                    array(vals)
                })
                .unwrap_or(Json::Null);
            let map_jets = result
                .coefficient_derivative_wrt_theta
                .map(|_| {
                    let vals = (0..2).flat_map(|brick| {
                        [0.125, 0.37, 0.5, 0.83].into_iter().map(move |u| {
                            n(result
                                .map_derivative_wrt_theta(brick, u)
                                .expect("smooth probe point"))
                        })
                    });
                    array(vals)
                })
                .unwrap_or(Json::Null);
            Ok(object([
                ("case", name),
                ("accepted", b(true)),
                ("unique", b(result.unique)),
                ("regime", s(&format!("{:?}", result.regime))),
                ("coefficients", coefficients),
                ("derivative", derivative),
                ("objective", n(result.objective)),
                ("kkt_residual", n(result.kkt.residual)),
                ("tangent_curvature", n(result.kkt.tangent_curvature)),
                ("equality_rank", u(result.kkt.equality_rank)),
                (
                    "activation_threshold",
                    result.activation_threshold.map(n).unwrap_or(Json::Null),
                ),
                ("minimum_map_speed", speed.map(n).unwrap_or(Json::Null)),
                ("map_values", values),
                ("map_jets", map_jets),
            ]))
        }
    }
}

fn p7(input: &Json) -> Result<Json, &'static str> {
    let read_circle = |key: &str| -> Result<[i64; 3], &'static str> {
        let values = input.field(key)?.items()?;
        if values.len() != 3 {
            return Err("circle must contain three integers");
        }
        Ok([
            values[0].number()?,
            values[1].number()?,
            values[2].number()?,
        ])
    };
    let name = s(case_name(input)?);
    let first = read_circle("first")?;
    let second = read_circle("second")?;
    match intersect_circles(first, second) {
        Err(error) => Ok(object([
            ("case", name),
            ("accepted", b(false)),
            ("error", s(error)),
        ])),
        Ok(result) => {
            let points = array(result.points.iter().map(|point| {
                object([
                    ("point_mm", array(point.point_mm.into_iter().map(n))),
                    (
                        "jet_mm_per_mm",
                        array(
                            point
                                .jet_mm_per_mm
                                .into_iter()
                                .map(|row| array(row.into_iter().map(n))),
                        ),
                    ),
                    ("condition", n(point.condition_estimate)),
                    ("constraint_residual_mm2", n(point.constraint_residual_mm2)),
                ])
            }));
            let tangent = result
                .tangent_point_mm
                .map(|p| array(p.into_iter().map(n)))
                .unwrap_or(Json::Null);
            Ok(object([
                ("case", name),
                ("accepted", b(true)),
                ("relation", s(&format!("{:?}", result.relation))),
                ("cardinality", s(&format!("{:?}", result.cardinality))),
                ("count", u(result.intersection_count)),
                ("rank", result.rank.map(u).unwrap_or(Json::Null)),
                ("points", points),
                ("tangent_point_mm", tangent),
                (
                    "derivative_unavailable",
                    result.derivative_unavailable.map(s).unwrap_or(Json::Null),
                ),
            ]))
        }
    }
}

fn run() -> Result<(), &'static str> {
    let mode = std::env::args().nth(1).ok_or("expected p6 or p7 mode")?;
    let mut source = String::new();
    io::stdin()
        .read_to_string(&mut source)
        .map_err(|_| "stdin read failed")?;
    let input = parse(&source)?;
    let output = array(
        input
            .items()?
            .iter()
            .map(|case| match mode.as_str() {
                "p6" => p6(case),
                "p7" => p7(case),
                _ => Err("expected p6 or p7 mode"),
            })
            .collect::<Result<Vec<_>, _>>()?,
    );
    println!("{output}");
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("p6b_oracle_probe: {error}");
        std::process::exit(2);
    }
}
