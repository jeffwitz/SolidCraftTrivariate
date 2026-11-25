//! Import the literal `solid.box` command subset into an audited G16 program.
//!
//! The importer handles one new base and up to eight named joins/cuts. Only
//! millimetre literals and explicit corner/center placements are accepted;
//! expressions, participants, ignored options, and multi-body operations are
//! rejected before any native model is produced.

use crate::{
    interchange,
    json::{object, Json},
    orthogonal_program::{build_orthogonal_program, Bounds, BoxFeature, BoxOp, OrthogonalProgram},
};
use std::collections::BTreeSet;

fn literal(value: &Json) -> Result<&str, &'static str> {
    match value {
        Json::Number(value) => Ok(value),
        _ => Err("solid.box values must be numeric literals"),
    }
}

/// Parse a strict decimal literal exactly as a reduced rational pair.
/// Leading plus signs and noncanonical leading zeroes are rejected, while
/// signed exponents remain valid JSON number syntax.
pub(crate) fn decimal_rational(token: &str) -> Result<(i128, i128), &'static str> {
    let (mantissa, exponent) = match token.split_once(['e', 'E']) {
        Some((mantissa, exponent)) => (
            mantissa,
            exponent
                .parse::<i32>()
                .map_err(|_| "solid.box numeric range")?,
        ),
        None => (token, 0),
    };
    let (negative, mantissa) = if let Some(value) = mantissa.strip_prefix('-') {
        (true, value)
    } else {
        (false, mantissa)
    };
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    if whole.is_empty()
        || !whole.bytes().all(|byte| byte.is_ascii_digit())
        || !fraction.bytes().all(|byte| byte.is_ascii_digit())
        || (whole.len() > 1 && whole.starts_with('0'))
        || (mantissa.ends_with('.') && fraction.is_empty())
    {
        return Err("invalid solid.box decimal literal");
    }
    if let Some((_, exponent)) = token.split_once(['e', 'E']) {
        let digits = exponent.trim_start_matches(['+', '-']);
        if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
            return Err("invalid solid.box decimal exponent");
        }
    }
    let mut digits = format!("{whole}{fraction}");
    let mut power = exponent
        .checked_sub(i32::try_from(fraction.len()).map_err(|_| "solid.box numeric range")?)
        .ok_or("solid.box numeric range")?;
    while digits.ends_with('0') {
        digits.pop();
        power = power.checked_add(1).ok_or("solid.box numeric range")?;
    }
    let numerator = if digits.is_empty() {
        0
    } else {
        digits
            .parse::<i128>()
            .map_err(|_| "solid.box numeric range")?
    };
    let (mut numerator, denominator) = if power >= 0 {
        let factor = 10_i128
            .checked_pow(power as u32)
            .ok_or("solid.box numeric range")?;
        (
            numerator
                .checked_mul(factor)
                .ok_or("solid.box numeric range")?,
            1,
        )
    } else {
        (
            numerator,
            10_i128
                .checked_pow(power.unsigned_abs())
                .ok_or("solid.box numeric range")?,
        )
    };
    if negative {
        numerator = numerator.checked_neg().ok_or("solid.box numeric range")?;
    }
    if numerator == 0 {
        return Ok((0, 1));
    }
    let gcd = gcd_u128(numerator.unsigned_abs(), denominator as u128);
    let divisor = i128::try_from(gcd).map_err(|_| "solid.box numeric range")?;
    Ok((numerator / divisor, denominator / divisor))
}

fn decimal(value: &Json) -> Result<(i128, i128), &'static str> {
    decimal_rational(literal(value)?)
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn exact_integer(value: &Json) -> Result<i32, &'static str> {
    let (numerator, denominator) = decimal(value)?;
    if numerator % denominator != 0 {
        return Err("solid.box dimension or corner must lie on an integer millimetre plane");
    }
    i32::try_from(numerator / denominator).map_err(|_| "solid.box coordinate range")
}

fn point3(value: &Json) -> Result<[Json; 3], &'static str> {
    match value {
        Json::Array(values) if values.len() == 3 => {
            Ok([values[0].clone(), values[1].clone(), values[2].clone()])
        }
        Json::Object(values)
            if values.len() == 3 && ["x", "y", "z"].iter().all(|key| values.contains_key(*key)) =>
        {
            Ok([
                values.get("x").unwrap().clone(),
                values.get("y").unwrap().clone(),
                values.get("z").unwrap().clone(),
            ])
        }
        _ => Err("solid.box placement must be a three-coordinate point"),
    }
}

fn integer_bounds(params: &Json) -> Result<Bounds, &'static str> {
    let dimensions = [
        exact_integer(params.field("length")?)?,
        exact_integer(params.field("width")?)?,
        exact_integer(params.field("height")?)?,
    ];
    if dimensions.iter().any(|&length| length <= 0) {
        return Err("solid.box dimensions must be positive");
    }
    let has_corner = params.field("corner").is_ok();
    let has_center = params.field("center").is_ok();
    if has_corner == has_center {
        return Err("solid.box requires exactly one explicit corner or center");
    }
    let placement = point3(if has_corner {
        params.field("corner")?
    } else {
        params.field("center")?
    })?;
    let mut bounds = [[0_i32; 2]; 3];
    for axis in 0..3 {
        let (numerator, denominator) = decimal(&placement[axis])?;
        if has_corner {
            if numerator % denominator != 0 {
                return Err("solid.box corner must lie on an integer millimetre plane");
            }
            let lower =
                i32::try_from(numerator / denominator).map_err(|_| "solid.box coordinate range")?;
            bounds[axis] = [
                lower,
                lower
                    .checked_add(dimensions[axis])
                    .ok_or("solid.box coordinate range")?,
            ];
        } else {
            // Keep centers rational until subtracting/adding half the integer size.
            let twice_center = numerator
                .checked_mul(2)
                .ok_or("solid.box coordinate range")?;
            if twice_center % denominator != 0 {
                return Err("solid.box center is outside the exact half-mm lattice");
            }
            let center_twice = twice_center / denominator;
            let lower_twice = center_twice
                .checked_sub(i128::from(dimensions[axis]))
                .ok_or("solid.box coordinate range")?;
            let upper_twice = center_twice
                .checked_add(i128::from(dimensions[axis]))
                .ok_or("solid.box coordinate range")?;
            if lower_twice % 2 != 0 || upper_twice % 2 != 0 {
                return Err("solid.box center and size do not define integer planes");
            }
            bounds[axis] = [
                i32::try_from(lower_twice / 2).map_err(|_| "solid.box coordinate range")?,
                i32::try_from(upper_twice / 2).map_err(|_| "solid.box coordinate range")?,
            ];
        }
    }
    Ok(bounds)
}

fn command(command: &Json) -> Result<&Json, &'static str> {
    command.keys(&["command", "params"])?;
    if command.field("command")?.text()? != "solid.box" {
        return Err("solidcraft box history contains an unsupported command");
    }
    command.field("params")
}

fn parse_feature(params: &Json, base: bool) -> Result<(String, BoxOp, Bounds), &'static str> {
    let mut required = vec!["height", "length", "name", "operation", "width"];
    let has_corner = params.field("corner").is_ok();
    let has_center = params.field("center").is_ok();
    if has_corner {
        required.push("corner");
    }
    if has_center {
        required.push("center");
    }
    params.keys(&required)?;
    let id = params.field("name")?.text()?.to_owned();
    if id.is_empty()
        || id.len() > 48
        || !id
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || byte == b'_')
    {
        return Err("solid.box feature name is not a stable native id");
    }
    let operation = match params.field("operation")?.text()? {
        "new" if base => BoxOp::Add,
        "join" if !base => BoxOp::Add,
        "cut" if !base => BoxOp::Cut,
        "new" if !base => return Err("only the initial solid.box may create a body"),
        _ => return Err("solid.box operation is outside the qualified new/join/cut subset"),
    };
    if base != (id == "base") {
        return Err("the first solid.box must be named base and later ids must be distinct");
    }
    Ok((id, operation, integer_bounds(params)?))
}

fn recipe(program: &OrthogonalProgram) -> Json {
    let bounds = |bounds: &Bounds| {
        Json::Array(
            bounds
                .iter()
                .map(|interval| {
                    Json::Array(
                        interval
                            .iter()
                            .map(|value| Json::Number(value.to_string()))
                            .collect(),
                    )
                })
                .collect(),
        )
    };
    object([
        (
            "family",
            Json::String("orthogonal_feature_program_v1".into()),
        ),
        ("base_mm", bounds(&program.base_mm)),
        (
            "features",
            Json::Array(
                program
                    .features
                    .iter()
                    .map(|feature| {
                        object([
                            ("id", Json::String(feature.id.clone())),
                            (
                                "operation",
                                Json::String(
                                    match feature.operation {
                                        BoxOp::Add => "add",
                                        BoxOp::Cut => "cut",
                                    }
                                    .into(),
                                ),
                            ),
                            ("bounds_mm", bounds(&feature.bounds_mm)),
                        ])
                    })
                    .collect(),
            ),
        ),
    ])
}

/// Convert a literal ordered `solid.box` command history into a native audited model.
pub fn import_boxes(script: &Json) -> Result<Json, &'static str> {
    script.keys(&["commands"])?;
    let commands = script.field("commands")?.items()?;
    if commands.is_empty() || commands.len() > 9 {
        return Err("solid.box program must contain one base and at most eight features");
    }
    let mut ids = BTreeSet::new();
    let mut program: Option<OrthogonalProgram> = None;
    for (index, command_json) in commands.iter().enumerate() {
        let is_base = index == 0;
        let params = command(command_json)?;
        let (id, operation, bounds_mm) = parse_feature(params, is_base)?;
        if !ids.insert(id.clone()) {
            return Err("solid.box feature names must be unique");
        }
        if is_base {
            let candidate = OrthogonalProgram {
                base_mm: bounds_mm,
                features: Vec::new(),
            };
            build_orthogonal_program(&candidate)?;
            program = Some(candidate);
        } else {
            let candidate = program.as_mut().ok_or("missing initial solid.box")?;
            candidate.features.push(BoxFeature {
                id,
                operation,
                bounds_mm,
            });
            // Each prefix must still denote one admissible native solid. This follows
            // the host's ordered feature semantics and prevents hidden multi-body steps.
            build_orthogonal_program(candidate)?;
        }
    }
    let model = interchange::build(&recipe(
        program.as_ref().ok_or("missing initial solid.box")?,
    ))?;
    if interchange::audit(&model)?.field("status")?.text()? != "valid" {
        return Err("solid.box history failed native model audit");
    }
    Ok(model)
}

#[cfg(test)]
mod decimal_tests {
    use super::decimal_rational;

    #[test]
    fn decimal_literals_have_shared_exact_exponent_and_zero_semantics() {
        assert_eq!(decimal_rational("1e+2"), Ok((100, 1)));
        assert_eq!(decimal_rational("12.50"), Ok((25, 2)));
        assert_eq!(decimal_rational("1e-1"), Ok((1, 10)));
        assert_eq!(decimal_rational("-0.000e+17"), Ok((0, 1)));
        assert!(decimal_rational("+1").is_err());
        assert!(decimal_rational("01").is_err());
        assert!(decimal_rational("1e+").is_err());
        assert!(decimal_rational("1+1").is_err());
    }
}
