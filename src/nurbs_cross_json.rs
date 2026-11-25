//! Versioned JSON exchange for a transversely perforated NURBS extrusion.

use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    bracket_json::{decode_blocks, design, encode, topology},
    constraint_system::Ratio,
    incidence::enumerate_faces,
    json::{object, Json},
    nurbs_cross::{build_cross, NurbsCurve2},
};

pub(crate) const FAMILY: &str = "nurbs_cross_bores_v1";
pub(crate) const BOUNDED_FAMILY: &str = "nurbs_cross_bores_bounded_v1";
pub(crate) const OVAL_FAMILY: &str = "nurbs_cross_oval_v1";
pub(crate) const METHOD: &str = "canonical_nurbs_cross_reconstruction_and_exact_incidence";
pub(crate) const BOUNDED_METHOD: &str = "certified_rational_nurbs_graph_approximation";
pub(crate) const OVAL_METHOD: &str = "exact_symmetric_nurbs_oval_reconstruction";
const PARAMETERS: [&str; 2] = ["bounds_mm", "profiles"];

fn number(n: impl ToString) -> Json {
    Json::Number(n.to_string())
}
fn text(s: &str) -> Json {
    Json::String(s.to_owned())
}
fn array<T>(xs: impl IntoIterator<Item = T>, f: impl FnMut(T) -> Json) -> Json {
    Json::Array(xs.into_iter().map(f).collect())
}
fn pair(v: &Json) -> Result<[i64; 2], &'static str> {
    v.items()?
        .iter()
        .map(Json::number)
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| "expected_two_bounds")
}
fn parse_bounds(v: &Json) -> Result<[[i64; 2]; 3], &'static str> {
    v.items()?
        .iter()
        .map(pair)
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| "expected_three_axis_bounds")
}
fn curve_json(c: &NurbsCurve2) -> Json {
    object([
        ("degree", number(c.degree)),
        ("knots", array(c.knots.iter().copied(), number)),
        (
            "controls",
            array(c.controls.iter(), |p| array(p.iter().copied(), number)),
        ),
    ])
}
fn parse_curve(v: &Json) -> Result<NurbsCurve2, &'static str> {
    v.keys(&["degree", "knots", "controls"])?;
    let curve = NurbsCurve2 {
        degree: v.field("degree")?.number()?,
        knots: v
            .field("knots")?
            .items()?
            .iter()
            .map(Json::number)
            .collect::<Result<_, _>>()?,
        controls: v
            .field("controls")?
            .items()?
            .iter()
            .map(|p| {
                p.items()?
                    .iter()
                    .map(Json::number)
                    .collect::<Result<Vec<i64>, _>>()?
                    .try_into()
                    .map_err(|_| "expected_homogeneous_2d_control")
            })
            .collect::<Result<_, _>>()?,
    };
    if !(1..=3).contains(&curve.degree) {
        return Err("unsupported_nurbs_profile_degree");
    }
    Ok(curve)
}
fn parse_profiles(v: &Json) -> Result<[NurbsCurve2; 4], &'static str> {
    v.keys(&["x_min", "x_max", "z_min", "z_max"])?;
    Ok([
        parse_curve(v.field("x_min")?)?,
        parse_curve(v.field("x_max")?)?,
        parse_curve(v.field("z_min")?)?,
        parse_curve(v.field("z_max")?)?,
    ])
}
fn parse_center(v: &Json) -> Result<[i64; 3], &'static str> {
    v.items()?
        .iter()
        .map(Json::number)
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| "expected_three_center_coordinates")
}
fn parse_meridians(v: &Json) -> Result<[NurbsCurve2; 2], &'static str> {
    v.keys(&["x", "z"])?;
    Ok([parse_curve(v.field("x")?)?, parse_curve(v.field("z")?)?])
}
fn parse_tolerance(v: &Json) -> Result<Ratio, &'static str> {
    let [numerator, denominator] = v.items()? else {
        return Err("expected_rational_lateral_tolerance");
    };
    let ratio = Ratio::new(numerator.number()?, denominator.number()?)?;
    if ratio.numerator() <= 0 {
        return Err("lateral_tolerance_must_be_positive");
    }
    Ok(ratio)
}
fn params_json(bounds: [[i64; 2]; 3], profiles: &[NurbsCurve2; 4]) -> (Json, Json) {
    let b = array(bounds, |axis| array(axis, number));
    let p = object([
        ("x_min", curve_json(&profiles[0])),
        ("x_max", curve_json(&profiles[1])),
        ("z_min", curve_json(&profiles[2])),
        ("z_max", curve_json(&profiles[3])),
    ]);
    (b, p)
}
fn build_parts(
    bounds: [[i64; 2]; 3],
    profiles: &[NurbsCurve2; 4],
) -> Result<Vec<BracketBlock>, &'static str> {
    for interval in bounds {
        if interval[0] >= interval[1] {
            return Err("cross_bounds_must_be_strictly_increasing");
        }
    }
    build_cross(bounds, std::array::from_fn(|i| &profiles[i]))
}
fn build_bounded_parts(
    bounds: [[i64; 2]; 3],
    profiles: &[NurbsCurve2; 4],
    tolerance: Ratio,
    max_depth: usize,
) -> Result<Vec<BracketBlock>, &'static str> {
    for interval in bounds {
        if interval[0] >= interval[1] {
            return Err("cross_bounds_must_be_strictly_increasing");
        }
    }
    if max_depth > 12 {
        return Err("max_depth_exceeds_12");
    }
    crate::nurbs_cross::build_cross_bounded(
        bounds,
        std::array::from_fn(|i| &profiles[i]),
        tolerance,
        max_depth,
    )
}
fn certificate(blocks: &[BracketBlock]) -> Result<BracketCertificate, &'static str> {
    let incidence = enumerate_faces(blocks)?;
    let design_faces = incidence
        .exterior
        .iter()
        .cloned()
        .map(|side| (side, "boundary".to_owned()))
        .collect();
    Ok(BracketCertificate {
        incidence,
        determinant_lower_bound_mm3: 0.0,
        radial_error_upper_mm: 0.0,
        design_faces,
    })
}

/// Construct the candidate exclusively through the native NURBS-cross builder.
pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    match recipe.field("family")?.text()? {
        FAMILY => build_exact(recipe),
        BOUNDED_FAMILY => build_bounded(recipe),
        OVAL_FAMILY => build_oval(recipe),
        _ => Err("unsupported_family"),
    }
}
fn build_exact(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&[&["family"][..], &PARAMETERS].concat())?;
    let bounds = parse_bounds(recipe.field("bounds_mm")?)?;
    let profiles = parse_profiles(recipe.field("profiles")?)?;
    let blocks = build_parts(bounds, &profiles)?;
    let cert = certificate(&blocks)?;
    let mut model = encode(FAMILY, blocks, &cert);
    if let Json::Object(fields) = &mut model {
        let (b, p) = params_json(bounds, &profiles);
        fields.insert("bounds_mm".into(), b);
        fields.insert("profiles".into(), p);
    }
    Ok(model)
}
fn build_bounded(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(
        &[
            &["family"][..],
            &PARAMETERS,
            &["lateral_tolerance_mm", "max_depth"][..],
        ]
        .concat(),
    )?;
    let bounds = parse_bounds(recipe.field("bounds_mm")?)?;
    let profiles = parse_profiles(recipe.field("profiles")?)?;
    let tolerance = parse_tolerance(recipe.field("lateral_tolerance_mm")?)?;
    let max_depth = recipe.field("max_depth")?.number()?;
    let blocks = build_bounded_parts(bounds, &profiles, tolerance, max_depth)?;
    let cert = certificate(&blocks)?;
    let mut model = encode(BOUNDED_FAMILY, blocks, &cert);
    if let Json::Object(fields) = &mut model {
        let (b, p) = params_json(bounds, &profiles);
        fields.insert("bounds_mm".into(), b);
        fields.insert("profiles".into(), p);
        fields.insert(
            "lateral_tolerance_mm".into(),
            array([tolerance.numerator(), tolerance.denominator()], number),
        );
        fields.insert("max_depth".into(), number(max_depth));
    }
    Ok(model)
}
fn build_oval(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&["family", "bounds_mm", "center_mm", "meridians"])?;
    let bounds = parse_bounds(recipe.field("bounds_mm")?)?;
    let center = parse_center(recipe.field("center_mm")?)?;
    let meridians = parse_meridians(recipe.field("meridians")?)?;
    for interval in bounds {
        if interval[0] >= interval[1] {
            return Err("cross_bounds_must_be_strictly_increasing");
        }
    }
    let blocks = crate::nurbs_cross::build_oval(bounds, center, &meridians[0], &meridians[1])?;
    let cert = certificate(&blocks)?;
    let mut model = encode(OVAL_FAMILY, blocks, &cert);
    if let Json::Object(fields) = &mut model {
        fields.insert("bounds_mm".into(), array(bounds, |v| array(v, number)));
        fields.insert("center_mm".into(), array(center, number));
        fields.insert(
            "meridians".into(),
            object([
                ("x", curve_json(&meridians[0])),
                ("z", curve_json(&meridians[1])),
            ]),
        );
    }
    Ok(model)
}

/// Rebuild from the stored recipe, compare the complete coefficient ledger, and
/// independently re-enumerate exact whole-face incidences. This is a structural
/// certificate; continuous geometric claims remain with the construction audit.
pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    match model.field("schema")?.text()? {
        FAMILY => audit_exact(model),
        BOUNDED_FAMILY => audit_bounded(model),
        OVAL_FAMILY => audit_oval(model),
        _ => Err("unsupported_schema"),
    }
}
fn audit_exact(model: &Json) -> Result<Json, &'static str> {
    model.keys(
        &[
            &["schema", "blocks", "topology", "design_faces"][..],
            &PARAMETERS,
        ]
        .concat(),
    )?;
    let bounds = parse_bounds(model.field("bounds_mm")?)?;
    let profiles = parse_profiles(model.field("profiles")?)?;
    let expected = build_parts(bounds, &profiles)?;
    let actual = decode_blocks(model)?;
    if actual != expected {
        return Err("nurbs_cross_block_reconstruction_mismatch");
    }
    let cert = certificate(&actual)?;
    if model.field("topology")? != &topology(&cert.incidence)
        || model.field("design_faces")? != &design(&cert)
    {
        return Err("nurbs_cross_incidence_or_provenance_mismatch");
    }
    Ok(object([
        ("status", text("valid")),
        ("method", text(METHOD)),
        ("fidelity", text("EXACT_REPRESENTATION")),
        ("blocks", number(actual.len())),
        ("shared_faces", number(cert.incidence.interfaces.len())),
        ("exterior_sides", number(cert.incidence.exterior.len())),
        ("weights_status", text("positive_on_closed_domain")),
        (
            "ordinates_status",
            text("strictly_increasing_on_closed_domain"),
        ),
        (
            "lateral_order_status",
            text("strict_order_certified_on_closed_domain"),
        ),
        (
            "structural_reconstruction",
            text("exact_integer_coefficients_and_whole_face_incidence"),
        ),
        (
            "continuous_geometry",
            text("exact_rational_bernstein_controls_checked_for_order_and_positive_weights"),
        ),
        (
            "certificate_scope",
            text("recipe_reconstruction_and_exact_whole_face_incidence_only"),
        ),
    ]))
}
fn audit_bounded(model: &Json) -> Result<Json, &'static str> {
    model.keys(
        &[
            &["schema", "blocks", "topology", "design_faces"][..],
            &PARAMETERS,
            &["lateral_tolerance_mm", "max_depth"][..],
        ]
        .concat(),
    )?;
    let bounds = parse_bounds(model.field("bounds_mm")?)?;
    let profiles = parse_profiles(model.field("profiles")?)?;
    let tolerance = parse_tolerance(model.field("lateral_tolerance_mm")?)?;
    let max_depth = model.field("max_depth")?.number()?;
    let expected = build_bounded_parts(bounds, &profiles, tolerance, max_depth)?;
    let actual = decode_blocks(model)?;
    if actual != expected {
        return Err("nurbs_cross_bounded_reconstruction_mismatch");
    }
    let cert = certificate(&actual)?;
    if model.field("topology")? != &topology(&cert.incidence)
        || model.field("design_faces")? != &design(&cert)
    {
        return Err("nurbs_cross_incidence_or_provenance_mismatch");
    }
    let upper = tolerance.checked_mul(Ratio::new(2, 1)?)?;
    Ok(object([
        ("status", text("valid")),
        ("method", text(BOUNDED_METHOD)),
        ("blocks", number(actual.len())),
        ("shared_faces", number(cert.incidence.interfaces.len())),
        ("exterior_sides", number(cert.incidence.exterior.len())),
        ("fidelity", text("APPROX_CERTIFIED")),
        ("weights_status", text("positive_on_closed_domain")),
        (
            "ordinates_status",
            text("strictly_increasing_on_closed_domain"),
        ),
        (
            "lateral_order_status",
            text("strict_order_certified_for_constructed_rational_envelope"),
        ),
        (
            "hausdorff_error_upper_bound_mm",
            array([upper.numerator(), upper.denominator()], number),
        ),
        (
            "structural_reconstruction",
            text("exact_integer_coefficients_and_whole_face_incidence"),
        ),
        (
            "certificate_scope",
            text("bounded_lateral_profile_fidelity_and_exact_whole_face_incidence"),
        ),
        ("max_depth", number(max_depth)),
    ]))
}
fn audit_oval(model: &Json) -> Result<Json, &'static str> {
    model.keys(&[
        "schema",
        "blocks",
        "topology",
        "design_faces",
        "bounds_mm",
        "center_mm",
        "meridians",
    ])?;
    let bounds = parse_bounds(model.field("bounds_mm")?)?;
    let center = parse_center(model.field("center_mm")?)?;
    let meridians = parse_meridians(model.field("meridians")?)?;
    let expected = crate::nurbs_cross::build_oval(bounds, center, &meridians[0], &meridians[1])?;
    let actual = decode_blocks(model)?;
    if actual != expected {
        return Err("nurbs_cross_oval_block_reconstruction_mismatch");
    }
    let cert = certificate(&actual)?;
    if model.field("topology")? != &topology(&cert.incidence)
        || model.field("design_faces")? != &design(&cert)
    {
        return Err("nurbs_cross_oval_incidence_or_provenance_mismatch");
    }
    Ok(object([
        ("status", text("valid")),
        ("method", text(OVAL_METHOD)),
        ("fidelity", text("EXACT_REPRESENTATION")),
        ("blocks", number(actual.len())),
        ("shared_faces", number(cert.incidence.interfaces.len())),
        ("exterior_sides", number(cert.incidence.exterior.len())),
        (
            "structural_reconstruction",
            text("exact_integer_coefficients_and_whole_face_incidence"),
        ),
        (
            "certificate_scope",
            text("exact_symmetric_nurbs_reconstruction_and_whole_face_incidence"),
        ),
    ]))
}
