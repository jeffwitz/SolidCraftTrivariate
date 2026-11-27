//! Versioned construction and audit shared by the CLI and native callers.

use super::annular::{build_annular, AnnularBlock, AnnularPair};
use super::annular_audit::audit_annular;
use super::json::{object, Json};

/// Recipes using a fixed two-dimensional section extruded in one or two layers.
pub const SECTION_FAMILIES: [&str; 6] = [
    "reference_rectangle_v1",
    "core_disk_v1",
    "full_annulus_v1",
    "pierced_square_v1",
    "rounded_corner_v1",
    "layered_pocket_v1",
];

fn string(s: &str) -> Json {
    Json::String(s.to_owned())
}
fn num(n: impl ToString) -> Json {
    Json::Number(n.to_string())
}
fn array<T>(s: impl IntoIterator<Item = T>, f: impl FnMut(T) -> Json) -> Json {
    Json::Array(s.into_iter().map(f).collect())
}
fn dims<const N: usize, T: std::str::FromStr>(v: &Json) -> Result<[T; N], &'static str> {
    v.items()?
        .iter()
        .map(Json::number)
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| "incorrect array length")
}
fn id(s: &Json, values: &[&'static str]) -> Result<&'static str, &'static str> {
    let s = s.text()?;
    values
        .iter()
        .copied()
        .find(|v| s == *v)
        .ok_or("unknown identifier")
}

fn annular_json(m: AnnularPair) -> Json {
    object([
        ("schema", string("annular_pair_v1")),
        ("heights_mm", array(m.heights_mm, num)),
        (
            "blocks",
            array(m.blocks, |b| {
                object([
                    ("id", string(b.id)),
                    ("degrees", array(b.degrees, num)),
                    ("controls", array(b.controls, |p| array(p, num))),
                ])
            }),
        ),
        ("interface", array(m.interface, string)),
    ])
}

fn decode_annular(model: &Json) -> Result<AnnularPair, &'static str> {
    model.keys(&["schema", "heights_mm", "blocks", "interface"])?;
    let heights_mm = dims(model.field("heights_mm")?)?;
    let mut blocks = vec![];
    for b in model.field("blocks")?.items()? {
        b.keys(&["id", "degrees", "controls"])?;
        let controls = b
            .field("controls")?
            .items()?
            .iter()
            .map(dims::<4, f64>)
            .collect::<Result<_, _>>()?;
        blocks.push(AnnularBlock {
            id: id(b.field("id")?, &["lower", "upper"])?,
            degrees: dims(b.field("degrees")?)?,
            controls,
        });
    }
    let interface = model
        .field("interface")?
        .items()?
        .iter()
        .map(|v| id(v, &["lower.w1", "upper.w0", "identity_uv"]))
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| "incorrect incidence count")?;
    Ok(AnnularPair {
        heights_mm,
        blocks,
        interface,
    })
}

// The historical affine family admits arbitrary positive integer widths.
// Its only arithmetic operation is exact addition of two nonnegative integers.
fn natural(v: &Json) -> Result<&str, &'static str> {
    match v {
        Json::Number(n)
            if n.as_bytes()
                .first()
                .is_some_and(|b| matches!(b, b'1'..=b'9'))
                && n.bytes().all(|b| b.is_ascii_digit()) =>
        {
            Ok(n)
        }
        _ => Err("expected a positive integer width"),
    }
}
fn widths(v: &Json) -> Result<[&str; 2], &'static str> {
    let [a, b] = v.items()? else {
        return Err("expected two widths");
    };
    Ok([natural(a)?, natural(b)?])
}
fn add_decimal(a: &str, b: &str) -> String {
    let mut x = a.bytes().rev();
    let mut y = b.bytes().rev();
    let mut digits = vec![];
    let mut carry = 0;
    loop {
        let p = x.next();
        let q = y.next();
        if p.is_none() && q.is_none() && carry == 0 {
            break;
        }
        let v = p.map_or(0, |c| c - b'0') + q.map_or(0, |c| c - b'0') + carry;
        digits.push(b'0' + v % 10);
        carry = v / 10;
    }
    digits.reverse();
    String::from_utf8(digits).expect("decimal digits")
}
fn affine_model(w: [&str; 2]) -> Json {
    let total = add_decimal(w[0], w[1]);
    let blocks = [("left", "0", w[0]), ("right", w[0], total.as_str())];
    object([
        ("schema", string("affine_pair_v1")),
        ("widths_mm", array(w, num)),
        (
            "blocks",
            array(blocks, |(name, x0, x1)| {
                let mut points = vec![];
                for x in [x0, x1] {
                    for y in 0..2 {
                        for z in 0..2 {
                            points.push(Json::Array(vec![num(x), num(y), num(z), num(1)]));
                        }
                    }
                }
                object([("id", string(name)), ("controls", Json::Array(points))])
            }),
        ),
        (
            "interface",
            array(["left.u1", "right.u0", "identity_vw"], string),
        ),
    ])
}
fn audit_affine(m: &Json) -> Result<Json, &'static str> {
    m.keys(&["schema", "widths_mm", "blocks", "interface"])?;
    let w = widths(m.field("widths_mm")?)?;
    let sum = add_decimal(w[0], w[1]);
    let blocks = m.field("blocks")?.items()?;
    if blocks.len() != 2
        || *m.field("interface")? != array(["left.u1", "right.u0", "identity_vw"], string)
    {
        return Err("affine incidence");
    }
    for (index, b) in blocks.iter().enumerate() {
        b.keys(&["id", "controls"])?;
        if b.field("id")?.text()? != ["left", "right"][index] {
            return Err("affine block identity");
        }
        let pts = b.field("controls")?.items()?;
        if pts.len() != 8 {
            return Err("affine control count");
        }
        for (n, p) in pts.iter().enumerate() {
            let xs = ["0", w[0], sum.as_str()];
            if *p
                != Json::Array(vec![
                    num(xs[index + n / 4]),
                    num((n / 2) % 2),
                    num(n % 2),
                    num(1),
                ])
            {
                return Err("affine coefficient mismatch");
            }
        }
    }
    let min = w
        .into_iter()
        .min_by(|a, b| a.len().cmp(&b.len()).then(a.cmp(b)))
        .unwrap();
    Ok(object([
        ("status", string("valid")),
        ("method", string("exact_integer_affine_coefficients")),
        ("min_det_mm3", num(min)),
        ("volume_mm3", num(sum)),
        ("shared_faces", num(1)),
        ("exterior_sides", num(10)),
    ]))
}

pub fn build(recipe: &Json) -> Result<Json, &'static str> {
    match recipe.field("family")?.text()? {
        crate::curved_sketch_audit::FAMILY => crate::curved_sketch_audit::build(recipe),
        super::affine_image::FAMILY => super::affine_image::build(recipe),
        "ruled_section_sweep_v1" => super::ruled_json::build(recipe),
        super::coaxial_json::FAMILY => super::coaxial_json::build(recipe),
        super::circular_bend_json::FAMILY => super::circular_bend_json::build(recipe),
        super::orthogonal_json::FAMILY => super::orthogonal_json::build(recipe),
        super::rational_sector_json::FAMILY => super::rational_sector_json::build(recipe),
        super::nurbs_cross_json::FAMILY => super::nurbs_cross_json::build(recipe),
        super::nurbs_cross_json::BOUNDED_FAMILY => super::nurbs_cross_json::build(recipe),
        super::nurbs_cross_json::OVAL_FAMILY => super::nurbs_cross_json::build(recipe),
        super::counterbores_json::FAMILY => super::counterbores_json::build(recipe),
        super::countersinks_json::FAMILY => super::countersinks_json::build(recipe),
        super::sketch_layers_json::FAMILY | super::sketch_layers_json::BOOLEAN_FAMILY => {
            super::sketch_layers_json::build(recipe)
        }
        super::g2_cross_plate_json::FAMILY | super::g2_cross_plate_json::PRIMITIVE | super::g2_cross_plate_json::FACE_POCKET => super::g2_cross_plate_json::build(recipe),
        super::nurbs_multiloop_json::FAMILY | super::nurbs_multiloop_json::EDGE_FAMILY => {
            super::nurbs_multiloop_json::build(recipe)
        }
        super::nurbs_sketch_json::FAMILY => super::nurbs_sketch_json::build(recipe),
        super::orthogonal_sketch_json::FAMILY => super::orthogonal_sketch_json::build(recipe),
        super::oblique_pockets_json::FAMILY => super::oblique_pockets_json::build(recipe),
        super::box_cylinder_boolean_json::FAMILY => super::box_cylinder_boolean_json::build(recipe),
        super::box_cylinder_boolean::SUITE_FAMILY => super::box_cylinder_boolean_json::build_suite(recipe),
        super::interval_bores_json::FAMILY => super::interval_bores_json::build(recipe),
        super::circular_polygon_json::FAMILY => super::circular_polygon_json::build(recipe),
        super::polygon_json::FAMILY | super::polygon_json::PERFORATED => {
            super::polygon_json::build(recipe)
        }
        "axial_boolean_v1" => super::axial_json::build(recipe),
        family
            if family == "filleted_bracket_v1"
                || SECTION_FAMILIES.contains(&family)
                || super::bracket_json::PARAMETRIC_FAMILIES.contains(&family) =>
        {
            super::bracket_json::build(recipe)
        }
        "affine_pair_v1" => {
            recipe.keys(&["family", "widths_mm"])?;
            Ok(affine_model(widths(recipe.field("widths_mm")?)?))
        }
        "annular_pair_v1" => {
            recipe.keys(&["family", "heights_mm"])?;
            Ok(annular_json(build_annular(dims(
                recipe.field("heights_mm")?,
            )?)?))
        }
        _ => Err("unsupported family"),
    }
}

/// Report planar junctions and transverse cuts before creating any volume blocks.
pub fn plan_sketch(recipe: &Json) -> Result<Json, &'static str> {
    if matches!(
        recipe.field("family")?.text()?,
        crate::sketch_layers_json::FAMILY | crate::sketch_layers_json::BOOLEAN_FAMILY
    ) {
        return crate::sketch_layers_json::plan_sketch(recipe);
    }
    if recipe.field("family")?.text()? == "primitive_sweep_v1" {
        return crate::primitive_template_json::plan(recipe);
    }
    crate::sketch_partition_plan::plan(recipe).map(|(_, report)| report)
}

pub fn audit(model: &Json) -> Result<Json, &'static str> {
    match model.field("schema")?.text()? {
        crate::curved_sketch_audit::FAMILY => crate::curved_sketch_audit::audit(model),
        super::affine_image::FAMILY => super::affine_image::audit(model),
        "ruled_section_sweep_v1" => super::ruled_json::audit(model),
        super::coaxial_json::FAMILY => super::coaxial_json::audit(model),
        super::circular_bend_json::FAMILY => super::circular_bend_json::audit(model),
        super::orthogonal_json::FAMILY => super::orthogonal_json::audit(model),
        super::rational_sector_json::FAMILY => super::rational_sector_json::audit(model),
        super::nurbs_cross_json::FAMILY => super::nurbs_cross_json::audit(model),
        super::nurbs_cross_json::BOUNDED_FAMILY => super::nurbs_cross_json::audit(model),
        super::nurbs_cross_json::OVAL_FAMILY => super::nurbs_cross_json::audit(model),
        super::counterbores_json::FAMILY => super::counterbores_json::audit(model),
        super::countersinks_json::FAMILY => super::countersinks_json::audit(model),
        super::sketch_layers_json::FAMILY | super::sketch_layers_json::BOOLEAN_FAMILY => {
            super::sketch_layers_json::audit(model)
        }
        super::g2_cross_plate_json::FAMILY | super::g2_cross_plate_json::PRIMITIVE | super::g2_cross_plate_json::FACE_POCKET => super::g2_cross_plate_json::audit(model),
        super::nurbs_multiloop_json::FAMILY | super::nurbs_multiloop_json::EDGE_FAMILY => {
            super::nurbs_multiloop_json::audit(model)
        }
        super::nurbs_sketch_json::FAMILY => super::nurbs_sketch_json::audit(model),
        super::orthogonal_sketch_json::FAMILY => super::orthogonal_sketch_json::audit(model),
        super::oblique_pockets_json::FAMILY => super::oblique_pockets_json::audit(model),
        super::box_cylinder_boolean_json::FAMILY => super::box_cylinder_boolean_json::audit(model),
        super::box_cylinder_boolean::SUITE_FAMILY => super::box_cylinder_boolean_json::audit_suite(model),
        super::interval_bores_json::FAMILY => super::interval_bores_json::audit(model),
        super::circular_polygon_json::FAMILY => super::circular_polygon_json::audit(model),
        super::polygon_json::FAMILY | super::polygon_json::PERFORATED => {
            super::polygon_json::audit(model)
        }
        "axial_boolean_v1" => super::axial_json::audit(model),
        family
            if family == "filleted_bracket_v1"
                || SECTION_FAMILIES.contains(&family)
                || super::bracket_json::PARAMETRIC_FAMILIES.contains(&family) =>
        {
            super::bracket_json::audit(model)
        }
        "affine_pair_v1" => audit_affine(model),
        "annular_pair_v1" => {
            let cert = audit_annular(&decode_annular(model)?)?;
            Ok(object([
                ("status", string("valid")),
                (
                    "method",
                    string("exact_canonical_coefficients_and_closed_domain_bounds"),
                ),
                ("shared_faces", num(1)),
                ("exterior_sides", num(10)),
                ("min_weight", num(cert.min_weight)),
                (
                    "determinant_lower_bound_mm3",
                    num(cert.determinant_lower_bound_mm3),
                ),
                ("radial_error_upper_mm", num(cert.radial_error_upper_mm)),
                ("height_error_upper_mm", num(cert.height_error_upper_mm)),
            ]))
        }
        _ => Ok(object([
            ("status", string("indeterminate")),
            ("reason", string("unknown_schema")),
        ])),
    }
}
