//! Translate the supported SolidCraft command history into an audited bracket recipe.
//! Geometry is recognized from command semantics. Unrecognized geometric parameters
//! are rejected before any native model is adopted.
use crate::{
    interchange,
    json::{object, Json},
};

/// Convert a decimal JSON token to integral eighths of a millimeter, without a
/// floating-point comparison or rounding onto the supported parameter lattice.
fn eighths(x: &Json) -> Result<i64, &'static str> {
    let Json::Number(s) = x else {
        return Err("expected literal SolidCraft length");
    };
    let (mantissa, exp) = match s.split_once(['e', 'E']) {
        Some((m, e)) => (m, e.parse::<i32>().map_err(|_| "SolidCraft numeric range")?),
        None => (s.as_str(), 0),
    };
    let negative = mantissa.starts_with('-');
    let m = mantissa.trim_start_matches('-');
    let (whole, frac) = m.split_once('.').unwrap_or((m, ""));
    let mut digits = format!("{whole}{frac}");
    if digits.bytes().all(|b| b == b'0') {
        return Ok(0);
    }
    let mut power = exp
        .checked_sub(i32::try_from(frac.len()).map_err(|_| "SolidCraft numeric range")?)
        .ok_or("SolidCraft numeric range")?;
    while digits.ends_with('0') {
        digits.pop();
        power = power.checked_add(1).ok_or("SolidCraft numeric range")?;
    }
    let value = digits
        .parse::<i128>()
        .ok()
        .and_then(|x| x.checked_mul(8))
        .ok_or("SolidCraft numeric range")?;
    let scale = 10_i128
        .checked_pow(power.unsigned_abs())
        .ok_or("SolidCraft numeric range")?;
    let value = if power >= 0 {
        value.checked_mul(scale).ok_or("SolidCraft numeric range")?
    } else {
        if value % scale != 0 {
            return Err("SolidCraft length is outside the eighth-millimeter lattice");
        }
        value / scale
    };
    i64::try_from(if negative { -value } else { value }).map_err(|_| "SolidCraft numeric range")
}
fn vector<const N: usize>(x: &Json) -> Result<[i64; N], &'static str> {
    x.items()?
        .iter()
        .map(eighths)
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| "SolidCraft vector length")
}
fn command<'a>(x: &'a Json, expected: &str) -> Result<&'a Json, &'static str> {
    x.keys(&["command", "params"])?;
    if x.field("command")?.text()? != expected {
        return Err("unsupported SolidCraft operation sequence");
    }
    x.field("params")
}
fn named(params: &Json, required: &[&str], name: &str) -> Result<(), &'static str> {
    if let Ok(v) = params.field(name) {
        let _ = v.text()?;
        let mut keys = required.to_vec();
        keys.push(name);
        params.keys(&keys)
    } else {
        params.keys(required)
    }
}

/// Accept the seven-command, extruded L bracket with independent through holes.
/// The two holes commute because their material cuts are strictly disjoint.
pub fn bracket_recipe(script: &Json) -> Result<Json, &'static str> {
    script.keys(&["commands"])?;
    let [create, poly, finish, extrude, fillet, h1, h2] = script.field("commands")?.items()? else {
        return Err("unsupported SolidCraft history length");
    };
    let sketch = command(create, "sketch.create")?;
    named(sketch, &["plane"], "name")?;
    if sketch.field("plane")?.text()? != "XY" {
        return Err("unsupported SolidCraft sketch plane");
    }
    let poly = command(poly, "sketch.line")?;
    poly.keys(&["closed", "points"])?;
    if poly.field("closed")? != &Json::Bool(true) {
        return Err("SolidCraft bracket requires a closed contour");
    }
    let vertices = poly
        .field("points")?
        .items()?
        .iter()
        .map(vector::<2>)
        .collect::<Result<Vec<_>, _>>()?;
    let expected = [[0, 0], [60, 0], [60, 8], [8, 8], [8, 50], [0, 50]].map(|p| p.map(|x| 8 * x));
    if vertices.len() != 6 || !(0..6).any(|k| (0..6).all(|j| vertices[j] == expected[(j + k) % 6]))
    {
        return Err("unrecognized SolidCraft bracket contour");
    }
    command(finish, "sketch.finish")?.keys(&[])?;
    let extrude = command(extrude, "solid.extrude")?;
    named(extrude, &["distance"], "body_name")?;
    if eighths(extrude.field("distance")?)? != 320 {
        return Err("unsupported SolidCraft bracket extrusion");
    }
    let f = command(fillet, "solid.fillet")?;
    f.keys(&["edges", "radius"])?;
    let [edge] = f.field("edges")?.items()? else {
        return Err("unsupported SolidCraft fillet selection");
    };
    if eighths(f.field("radius")?)? != 48 || vector(edge)? != [64, 64, 160] {
        return Err("unsupported SolidCraft fillet");
    }
    let mut radii = [None; 2];
    for h in [h1, h2] {
        let h = command(h, "solid.hole")?;
        h.keys(&["position", "diameter"])?;
        let index = match vector(h.field("position")?)? {
            [320, 64, 160] => 0,
            [64, 256, 160] => 1,
            _ => return Err("unsupported SolidCraft through-hole support"),
        };
        // Radius in sixteenths equals diameter in eighths.
        let ticks = eighths(h.field("diameter")?)?;
        if !(40..=80).contains(&ticks) {
            return Err("SolidCraft through-hole radius outside bracket domain");
        }
        if radii[index].replace(ticks).is_some() {
            return Err("duplicate SolidCraft through-hole support");
        }
    }
    Ok(object([
        ("family", Json::String("parameterized_bracket_v1".into())),
        (
            "radii_sixteenths_mm",
            Json::Array(
                radii
                    .into_iter()
                    .map(|x| {
                        x.map(|n| Json::Number(n.to_string()))
                            .ok_or("missing SolidCraft through-hole")
                    })
                    .collect::<Result<_, _>>()?,
            ),
        ),
    ]))
}
/// Translate, construct, and independently audit before returning the native state.
pub fn import_script(script: &Json) -> Result<Json, &'static str> {
    if script.field("about").is_ok() {
        let first = script
            .field("commands")?
            .items()?
            .first()
            .ok_or("empty SolidCraft timeline")?;
        let first_type = first
            .field("params")?
            .field("feature")?
            .field("type")?
            .text()?;
        match first_type {
            "box" => return crate::solidcraft_box_timeline::import_box_timeline(script),
            "sketch" => {
                let commands = script.field("commands")?.items()?;
                let second_type = commands
                    .get(1)
                    .and_then(|step| step.field("params").ok())
                    .and_then(|params| params.field("feature").ok())
                    .and_then(|feature| feature.field("type").ok())
                    .and_then(|value| value.text().ok());
                if commands.len() == 2 && second_type == Some("extrude") {
                    return crate::solidcraft_polygon::import_polygon_timeline(script);
                }
            }
            _ => {}
        }
        return crate::solidcraft_timeline::translate_timeline(script);
    }
    let first = script
        .field("commands")?
        .items()?
        .first()
        .ok_or("empty SolidCraft script")?;
    if first.field("command")?.text()? == "solid.box" {
        return crate::solidcraft_boxes::import_boxes(script);
    }
    let candidate = interchange::build(&bracket_recipe(script)?)?;
    if interchange::audit(&candidate)?.field("status")?.text()? != "valid" {
        return Err("SolidCraft candidate failed native audit");
    }
    Ok(candidate)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_decimal_lattice_rejects_hidden_roundoff() {
        for (s, expected) in [
            ("8", 64),
            ("8.0000", 64),
            ("6.375", 51),
            ("6375e-3", 51),
            ("-0.00", 0),
            ("0.625E1", 50),
        ] {
            assert_eq!(eighths(&Json::Number(s.into())).unwrap(), expected);
        }
        for s in ["6.37500000000000000001", "8.1", "1e309", "1e-999"] {
            assert!(eighths(&Json::Number(s.into())).is_err());
        }
    }
}
