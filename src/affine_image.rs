//! Exact rational affine images of independently certified whole-face models.
//! For source homogeneous coordinates (p,w), store (A p + t w, d w).
use crate::{
    bracket_json::{decode_blocks, encode_blocks},
    interchange,
    json::{object, Json},
};

pub(crate) const FAMILY: &str = "affine_image_v1";
pub(crate) const METHOD: &str = "certified_rational_affine_image";
const EXACT_BINARY64_LIMIT: i128 = 1 << 53;

fn number(n: impl ToString) -> Json {
    Json::Number(n.to_string())
}

fn triple(j: &Json, limit: i32) -> Result<[i128; 3], &'static str> {
    let [a, b, c] = j.items()? else {
        return Err("expected_three_affine_coefficients");
    };
    let checked = |x: &Json| {
        let n = x.number::<i32>()?;
        if (-limit..=limit).contains(&n) {
            Ok(i128::from(n))
        } else {
            Err("affine_coefficient_range")
        }
    };
    Ok([checked(a)?, checked(b)?, checked(c)?])
}

struct Transform {
    matrix: [[i128; 3]; 3],
    offset: [i128; 3],
    denominator: i128,
    exact_binary64: bool,
    adjugate: [[i128; 3]; 3],
    determinant: i128,
}

fn dot(a: [i128; 3], b: [i128; 3]) -> i128 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}

impl Transform {
    fn parse(j: &Json) -> Result<Self, &'static str> {
        let exact_binary64 = match j.field("homogeneous_range") {
            Ok(value) => {
                j.keys(&["matrix", "translation", "denominator", "homogeneous_range"])?;
                match value.text()? {
                    "i64" => false,
                    "binary64_exact" => true,
                    _ => return Err("unknown_affine_homogeneous_range"),
                }
            }
            Err(_) => {
                j.keys(&["matrix", "translation", "denominator"])?;
                true
            }
        };
        let [a, b, c] = j.field("matrix")?.items()? else {
            return Err("expected_three_affine_rows");
        };
        let matrix = [triple(a, 64)?, triple(b, 64)?, triple(c, 64)?];
        let offset = triple(j.field("translation")?, 4096)?;
        let denominator = i128::from(j.field("denominator")?.number::<i32>()?);
        if !(1..=64).contains(&denominator) {
            return Err("affine_denominator_range");
        }
        let adjugate = std::array::from_fn(|i| {
            std::array::from_fn(|k| {
                let (r, s) = ((k + 1) % 3, (k + 2) % 3);
                let (u, v) = ((i + 1) % 3, (i + 2) % 3);
                matrix[r][u] * matrix[s][v] - matrix[r][v] * matrix[s][u]
            })
        });
        let determinant = dot(matrix[0], adjugate.map(|row| row[0]));
        if determinant <= 0 {
            return Err("affine_image_requires_positive_determinant");
        }
        Ok(Self {
            matrix,
            offset,
            denominator,
            exact_binary64,
            adjugate,
            determinant,
        })
    }

    fn forward(&self, point: [i64; 4]) -> Result<[i64; 4], &'static str> {
        let [x, y, z, w] = point.map(i128::from);
        let xyz =
            std::array::from_fn::<_, 3, _>(|i| dot(self.matrix[i], [x, y, z]) + self.offset[i] * w);
        let image = [xyz[0], xyz[1], xyz[2], self.denominator * w];
        if self.exact_binary64 && image.iter().any(|v| v.abs() > EXACT_BINARY64_LIMIT) {
            return Err("affine_image_exact_coordinate_range");
        }
        let [a, b, c, w] = image.map(|v| i64::try_from(v).map_err(|_| "affine_image_i64_overflow"));
        Ok([a?, b?, c?, w?])
    }

    fn verifies_inverse(&self, source: [i64; 4], image: [i64; 4]) -> bool {
        let old = source.map(i128::from);
        let new = image.map(i128::from);
        let rhs = std::array::from_fn(|i| new[i] - self.offset[i] * old[3]);
        (!self.exact_binary64 || new.iter().all(|v| v.abs() <= EXACT_BINARY64_LIMIT))
            && new[3] == self.denominator * old[3]
            && (0..3).all(|i| dot(self.adjugate[i], rhs) == self.determinant * old[i])
    }
}

fn source_certificate(model: &Json) -> Result<Json, &'static str> {
    if model.field("schema")?.text()? == FAMILY {
        return Err("nested_affine_image");
    }
    let report = interchange::audit(model)?;
    if report.field("status")?.text()? != "valid"
        || !crate::boundary_shell::supports_source_method(report.field("method")?.text()?)
    {
        return Err("unsupported_affine_source");
    }
    Ok(report)
}

pub(crate) fn build(recipe: &Json) -> Result<Json, &'static str> {
    recipe.keys(&["family", "source", "transform"])?;
    let inner = recipe.field("source")?;
    if inner.field("family")?.text()? == FAMILY {
        return Err("nested_affine_image");
    }
    let transform = recipe.field("transform")?;
    let f = Transform::parse(transform)?;
    let source = interchange::build(inner)?;
    source_certificate(&source)?;
    let mut blocks = decode_blocks(&source)?;
    for block in &mut blocks {
        for p in &mut block.controls {
            *p = f.forward(*p)?;
        }
    }
    Ok(object([
        ("schema", Json::String(FAMILY.into())),
        ("blocks", encode_blocks(blocks)),
        ("topology", source.field("topology")?.clone()),
        ("design_faces", source.field("design_faces")?.clone()),
        ("transform", transform.clone()),
        ("source_model", source),
    ]))
}

pub(crate) fn audit(model: &Json) -> Result<Json, &'static str> {
    model.keys(&[
        "schema",
        "blocks",
        "topology",
        "design_faces",
        "transform",
        "source_model",
    ])?;
    let f = Transform::parse(model.field("transform")?)?;
    let source = model.field("source_model")?;
    let source_report = source_certificate(source)?;
    for key in ["topology", "design_faces"] {
        if model.field(key)? != source.field(key)? {
            return Err("affine_trace_or_provenance_mismatch");
        }
    }
    let base = decode_blocks(source)?;
    let image = decode_blocks(model)?;
    if base.len() != image.len() {
        return Err("affine_block_count_mismatch");
    }
    for (a, b) in base.iter().zip(&image) {
        if a.id != b.id
            || a.degrees != b.degrees
            || a.controls.len() != b.controls.len()
            || !a
                .controls
                .iter()
                .zip(&b.controls)
                .all(|(&s, &y)| f.verifies_inverse(s, y))
        {
            return Err("affine_inverse_coefficient_mismatch");
        }
    }
    // Inputs to these divisions are small exact integers. Round each bound outward.
    let scale = (f.determinant as f64 / f.denominator.pow(3) as f64).next_down();
    let lower = (scale
        * source_report
            .field("determinant_lower_bound_mm3")?
            .number::<f64>()?)
    .next_down();
    if !lower.is_finite() || lower <= 0.0 {
        return Err("affine_determinant_bound_indeterminate");
    }
    let max_row = f
        .matrix
        .iter()
        .map(|r| r.iter().map(|v| v.abs()).sum::<i128>())
        .max()
        .unwrap();
    let max_col = (0..3)
        .map(|i| f.matrix.iter().map(|r| r[i].abs()).sum::<i128>())
        .max()
        .unwrap();
    let norm = (((max_row * max_col) as f64).sqrt().next_up() / f.denominator as f64).next_up();
    let error = match source.field("schema")?.text()? {
        "orthogonal_feature_program_v1" | "polygon_prism_v1" | "polygonal_perforated_prism_v1" => {
            Some(0.0)
        }
        _ => source_report
            .field("hausdorff_error_upper_bound_mm")
            .ok()
            .map(Json::number::<f64>)
            .transpose()?,
    };
    Ok(object([
        ("status", Json::String("valid".into())),
        ("method", Json::String(METHOD.into())),
        ("source_certificate", source_report.clone()),
        ("blocks", number(image.len())),
        ("shared_faces", source_report.field("shared_faces")?.clone()),
        (
            "exterior_sides",
            source_report.field("exterior_sides")?.clone(),
        ),
        (
            "cavities",
            source_report
                .field("cavities")
                .cloned()
                .unwrap_or_else(|_| number(0)),
        ),
        ("determinant_numerator", number(f.determinant)),
        ("determinant_denominator", number(f.denominator.pow(3))),
        ("determinant_lower_bound_mm3", number(lower)),
        ("linear_operator_norm_upper", number(norm)),
        (
            "hausdorff_error_upper_bound_mm",
            error.map_or(Json::Null, |v| {
                number(if v == 0.0 { 0.0 } else { (v * norm).next_up() })
            }),
        ),
    ]))
}
