//! Exact generating support and source-span lineage for qualified spline families.
//!
//! This module recovers a tensor parent from constructive recipe data. It does
//! not infer a parent from unordered blocks and does not alter the source model.

use crate::{
    bracket_json,
    constraint_system::Ratio,
    json::{object, Json},
    nurbs_curve::extract_bezier,
};

type Result<T> = std::result::Result<T, &'static str>;
type HCurvePoint = [Ratio; 3];

const SOURCE_FAMILY: &str = "nurbs_sketch_strip_v1";
const SUPPORT_SCHEMA: &str = "tensor_support_v1";

fn number(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}

fn string(value: impl Into<String>) -> Json {
    Json::String(value.into())
}

fn arr<T>(values: impl IntoIterator<Item = T>, f: impl FnMut(T) -> Json) -> Json {
    Json::Array(values.into_iter().map(f).collect())
}

fn ratio(numerator: i128, denominator: i128) -> Result<Ratio> {
    Ratio::new(numerator, denominator)
}

fn checked_div(a: Ratio, b: Ratio) -> Result<Ratio> {
    let numerator = a
        .numerator()
        .checked_mul(b.denominator())
        .ok_or("spline support rational numerator overflow")?;
    let denominator = a
        .denominator()
        .checked_mul(b.numerator())
        .ok_or("spline support rational denominator overflow")?;
    Ratio::new(numerator, denominator)
}

fn parse_curve(curve: &Json) -> Result<CurveData> {
    curve.keys(&["degree", "knots", "controls"])?;
    let degree = curve.field("degree")?.number::<usize>()?;
    if !(1..=3).contains(&degree) {
        return Err("spline parent degree must be in 1..=3");
    }
    let knots = curve
        .field("knots")?
        .items()?
        .iter()
        .map(Json::number::<i64>)
        .collect::<Result<Vec<_>>>()?;
    let controls = curve
        .field("controls")?
        .items()?
        .iter()
        .map(|point| {
            let [x, y, w] = point.items()? else {
                return Err("spline parent expects homogeneous curve triples");
            };
            Ok([x.number::<i64>()?, y.number::<i64>()?, w.number::<i64>()?])
        })
        .collect::<Result<Vec<_>>>()?;
    if controls.len() < degree + 1 || knots.len() != controls.len() + degree + 1 {
        return Err("spline parent knot and control counts do not match");
    }
    if controls.iter().any(|point| point[2] <= 0) {
        return Err("spline parent weights must be strictly positive");
    }
    if knots.windows(2).any(|pair| pair[0] > pair[1]) {
        return Err("spline parent knots must be nondecreasing");
    }
    if knots[..=degree].iter().any(|&knot| knot != knots[degree])
        || knots[controls.len()..]
            .iter()
            .any(|&knot| knot != knots[controls.len()])
        || knots[degree] >= knots[controls.len()]
    {
        return Err("spline parent source must be clamped with a nonempty domain");
    }
    Ok(CurveData {
        degree,
        knots,
        controls,
    })
}

#[derive(Debug, Clone)]
struct CurveData {
    degree: usize,
    knots: Vec<i64>,
    controls: Vec<[i64; 3]>,
}

fn curve_json(curve: &Json) -> Result<CurveData> {
    parse_curve(curve)
}

fn hpoint(point: [i64; 3]) -> HCurvePoint {
    point.map(|value| Ratio::new(i128::from(value), 1).expect("integer ratio"))
}

fn blend(a: HCurvePoint, b: HCurvePoint, alpha: Ratio) -> Result<HCurvePoint> {
    let complement = Ratio::one().checked_sub(alpha)?;
    let mut result = [Ratio::zero(); 3];
    for component in 0..3 {
        result[component] = a[component]
            .checked_mul(alpha)?
            .checked_add(b[component].checked_mul(complement)?)?;
    }
    Ok(result)
}

fn insert_knot(
    degree: usize,
    knots: &[i64],
    controls: &[HCurvePoint],
    value: i64,
) -> Result<(Vec<i64>, Vec<HCurvePoint>)> {
    let last = controls
        .len()
        .checked_sub(1)
        .ok_or("empty spline parent curve")?;
    let multiplicity = knots.iter().filter(|&&knot| knot == value).count();
    if multiplicity >= degree {
        return Err("spline parent knot insertion reached full multiplicity");
    }
    let span = (degree..=last)
        .find(|&index| knots[index] <= value && value < knots[index + 1])
        .ok_or("spline parent insertion knot is outside the active domain")?;
    let mut next = vec![[Ratio::zero(); 3]; controls.len() + 1];
    let unchanged = span - degree + 1;
    next[..unchanged].copy_from_slice(&controls[..unchanged]);
    let blend_start = span - degree + 1;
    let blend_end = span - multiplicity;
    if blend_start <= blend_end {
        for index in blend_start..=blend_end {
            let denominator = i128::from(knots[index + degree]) - i128::from(knots[index]);
            let numerator = i128::from(value) - i128::from(knots[index]);
            if denominator <= 0 || numerator < 0 || numerator > denominator {
                return Err("invalid spline parent knot insertion interval");
            }
            next[index] = blend(
                controls[index],
                controls[index - 1],
                ratio(numerator, denominator)?,
            )?;
        }
    }
    let shifted = span - multiplicity + 1;
    next[shifted..last + 2].copy_from_slice(&controls[shifted - 1..=last]);
    let mut next_knots = knots.to_vec();
    next_knots.insert(span + 1, value);
    Ok((next_knots, next))
}

/// Extract Bézier controls while retaining the original homogeneous gauge.
/// The public curve extractor remains the canonical checked implementation;
/// this companion retains its exact rational scale for lineage comparison.
fn extract_parent_gauge(curve: &CurveData) -> Result<Vec<([i64; 2], Vec<HCurvePoint>)>> {
    // Reuse the native extractor as the input validity gate.
    extract_bezier(curve.degree, &curve.knots, &curve.controls)?;
    let mut knots = curve.knots.clone();
    let mut controls = curve
        .controls
        .iter()
        .copied()
        .map(hpoint)
        .collect::<Vec<_>>();
    let start = knots[curve.degree];
    let end = knots[curve.controls.len()];
    let internal = knots[curve.degree + 1..curve.controls.len()]
        .iter()
        .copied()
        .filter(|&knot| start < knot && knot < end)
        .collect::<std::collections::BTreeSet<_>>();
    for knot in internal {
        while knots.iter().filter(|&&entry| entry == knot).count() < curve.degree {
            (knots, controls) = insert_knot(curve.degree, &knots, &controls, knot)?;
        }
    }
    let last = controls.len() - 1;
    let mut spans = Vec::new();
    for index in curve.degree..=last {
        if knots[index] < knots[index + 1] {
            let first = index - curve.degree;
            spans.push((
                [knots[index], knots[index + 1]],
                controls[first..=index].to_vec(),
            ));
        }
    }
    Ok(spans)
}

/// The affine local-to-parent coordinate map for one source-cell restriction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParametricRestriction {
    pub parent_id: String,
    /// Increasing source intervals, retained in the source's actual units.
    pub intervals: [[i64; 2]; 3],
    /// Local coordinate axis `a` maps to source axis `axes[a]`.
    pub axes: [usize; 3],
    /// Reversal flags are applied to local coordinates before interval mapping.
    pub reversals: [bool; 3],
}

impl ParametricRestriction {
    /// Map local unit-cube coordinates into the parent's actual parameter domain.
    pub fn map_parameters(&self, local: [Ratio; 3]) -> Result<[Ratio; 3]> {
        let mut permutation = self.axes;
        permutation.sort();
        if permutation != [0, 1, 2] {
            return Err("spline parent restriction axes must be a permutation");
        }
        let mut source = [Ratio::zero(); 3];
        for (local_axis, local_parameter) in local.iter().enumerate() {
            let parent_axis = self.axes[local_axis];
            let [lo, hi] = self.intervals[parent_axis];
            if lo >= hi {
                return Err("spline parent restriction interval must be increasing");
            }
            let parameter = if self.reversals[local_axis] {
                Ratio::one().checked_sub(*local_parameter)?
            } else {
                *local_parameter
            };
            let length = i128::from(hi) - i128::from(lo);
            source[parent_axis] =
                ratio(i128::from(lo), 1)?.checked_add(parameter.checked_mul(ratio(length, 1)?)?)?;
        }
        Ok(source)
    }

    /// Return `d(source_axis)/d(local_axis)`, including orientation and scale.
    pub fn derivative_scale(&self, local_axis: usize) -> Result<Ratio> {
        if local_axis >= 3 {
            return Err("spline parent restriction axis is outside 0..3");
        }
        let mut permutation = self.axes;
        permutation.sort();
        if permutation != [0, 1, 2] {
            return Err("spline parent restriction axes must be a permutation");
        }
        let parent_axis = self.axes[local_axis];
        let [lo, hi] = self.intervals[parent_axis];
        if lo >= hi {
            return Err("spline parent restriction interval must be increasing");
        }
        let length = i128::from(hi) - i128::from(lo);
        ratio(
            if self.reversals[local_axis] {
                -length
            } else {
                length
            },
            1,
        )
    }
}

/// Exact evidence that a source interface inherits continuity from the parent.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContinuityEvidence {
    pub interface: String,
    pub coordinate_frame: String,
    pub direction: usize,
    pub knot: i64,
    pub multiplicity: usize,
    pub claimed_order: usize,
    pub method: String,
    pub assumptions: Vec<String>,
    pub status: String,
}

/// One exact source block restriction and its positive homogeneous gauge.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceSpan {
    pub source_block_id: String,
    pub restriction: ParametricRestriction,
    /// Positive reduced ratio: source-block gauge / recovered-parent gauge.
    pub gauge: [i128; 2],
    pub generating_operation: String,
    pub event_provenance: String,
}

/// Recovered tensor-product parent and its exact source-block lineage.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorSupport {
    pub id: String,
    pub source_revision: String,
    pub degrees: [usize; 3],
    pub knots: [Vec<i64>; 3],
    pub domain: [[i64; 2]; 3],
    pub control_shape: [usize; 3],
    /// Homogeneous `(wx, wy, wz, w)` controls, w-fast then v then u.
    pub homogeneous_controls: Vec<[i128; 4]>,
    pub source_spans: Vec<SourceSpan>,
    pub continuity_evidence: Vec<ContinuityEvidence>,
    pub generating_recipe: Json,
    pub certificate_references: Vec<String>,
}

fn recipe_from_model(model: &Json) -> Result<Json> {
    Ok(object([
        ("family", model.field("schema")?.clone()),
        ("lower", model.field("lower")?.clone()),
        ("upper", model.field("upper")?.clone()),
        ("height_mm", model.field("height_mm")?.clone()),
    ]))
}

fn fingerprint(recipe: &Json) -> String {
    // This is a deterministic lineage token, not a cryptographic integrity hash.
    let mut hash = 0xcbf29ce484222325_u64;
    for byte in recipe.to_string().bytes() {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x100000001b3);
    }
    format!("fnv1a64:{hash:016x}")
}

fn same_projective_vector(actual: &[i64; 4], expected: [Ratio; 4], gauge: Ratio) -> Result<bool> {
    for component in 0..4 {
        let scaled = expected[component].checked_mul(gauge)?;
        if scaled != Ratio::new(i128::from(actual[component]), 1)? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn source_restrictions(
    model: &Json,
    lower: &CurveData,
    upper: &CurveData,
    parent_id: &str,
    height: i32,
) -> Result<Vec<SourceSpan>> {
    let lower_spans = extract_parent_gauge(lower)?;
    let upper_spans = extract_parent_gauge(upper)?;
    if lower_spans.len() != upper_spans.len() {
        return Err("spline strip source span counts differ");
    }
    for ((lower_interval, _), (upper_interval, _)) in lower_spans.iter().zip(&upper_spans) {
        if lower_interval != upper_interval {
            return Err("spline strip source intervals differ");
        }
    }
    let intervals = model.field("span_intervals")?.items()?;
    let blocks = bracket_json::decode_blocks(model)?;
    if intervals.len() != lower_spans.len() || blocks.len() != lower_spans.len() {
        return Err("spline strip block and source span counts differ");
    }
    let mut result = Vec::with_capacity(blocks.len());
    for (index, (interval, block)) in intervals.iter().zip(&blocks).enumerate() {
        let [a, b] = interval.items()? else {
            return Err("spline strip interval must be a pair");
        };
        let a = a.number::<i64>()?;
        let b = b.number::<i64>()?;
        let reversed = a > b;
        let source_interval = if reversed { [b, a] } else { [a, b] };
        let lower_span = lower_spans
            .iter()
            .find(|span| span.0 == source_interval)
            .ok_or("spline strip source interval lineage mismatch")?;
        let upper_span = upper_spans
            .iter()
            .find(|span| span.0 == source_interval)
            .ok_or("spline strip source interval lineage mismatch")?;
        if lower_span.0 != source_interval || upper_span.0 != source_interval {
            return Err("spline strip source interval lineage mismatch");
        }
        let degree = lower.degree;
        if block.degrees != [degree, 1, 1] || block.controls.len() != (degree + 1) * 4 {
            return Err("spline strip source block degree or control count mismatch");
        }
        let lower_controls: Vec<_> = if reversed {
            lower_span.1.iter().rev().copied().collect()
        } else {
            lower_span.1.clone()
        };
        let upper_controls: Vec<_> = if reversed {
            upper_span.1.iter().rev().copied().collect()
        } else {
            upper_span.1.clone()
        };
        let gauge = checked_div(
            Ratio::new(i128::from(block.controls[0][3]), 1)?,
            lower_controls[0][2],
        )?;
        if gauge.numerator() <= 0 {
            return Err("spline strip source gauge must be positive");
        }
        for u in 0..=degree {
            for v in 0..=1 {
                let curve_control = if v == 0 {
                    lower_controls[u]
                } else {
                    upper_controls[u]
                };
                for w in 0..=1 {
                    let actual = &block.controls[(u * 2 + v) * 2 + w];
                    let expected = [
                        curve_control[0],
                        curve_control[1],
                        curve_control[2].checked_mul(Ratio::new(
                            i128::from(height)
                                .checked_mul(
                                    i128::try_from(w).map_err(|_| "invalid extrusion layer")?,
                                )
                                .ok_or("spline parent homogeneous z overflow")?,
                            1,
                        )?)?,
                        curve_control[2],
                    ];
                    if !same_projective_vector(actual, expected, gauge)? {
                        return Err("spline strip block is not the exact parent restriction");
                    }
                }
            }
        }
        let source_block_id = block.id.clone();
        result.push(SourceSpan {
            source_block_id: source_block_id.clone(),
            restriction: ParametricRestriction {
                parent_id: parent_id.to_owned(),
                intervals: [source_interval, [0, 1], [0, 1]],
                axes: [0, 1, 2],
                reversals: [reversed, false, false],
            },
            gauge: [gauge.numerator(), gauge.denominator()],
            generating_operation: "nurbs_strip_extrusion".into(),
            event_provenance: format!("source knot span {}..{}; model span {}", a, b, index),
        });
    }
    Ok(result)
}

fn continuity_evidence(curve: &CurveData) -> Vec<ContinuityEvidence> {
    let mut knots = std::collections::BTreeSet::new();
    for &knot in &curve.knots[curve.degree + 1..curve.controls.len()] {
        if curve.knots[curve.degree] < knot && knot < curve.knots[curve.controls.len()] {
            knots.insert(knot);
        }
    }
    knots
        .into_iter()
        .map(|knot| {
            let multiplicity = curve.knots.iter().filter(|&&entry| entry == knot).count();
            ContinuityEvidence {
                interface: format!("source-knot:{knot}"),
                coordinate_frame: "parent parameter coordinates in source knot units".into(),
                direction: 0,
                knot,
                multiplicity,
                claimed_order: curve.degree.saturating_sub(multiplicity),
                method: "exact clamped B-spline knot-multiplicity continuity".into(),
                assumptions: vec![
                    "source homogeneous controls are unchanged".into(),
                    "all source weights are strictly positive".into(),
                    "both strip boundaries share the same knot vector and weights".into(),
                ],
                status: "verified_from_source_support".into(),
            }
        })
        .collect()
}

/// Recover a complete tensor parent from an audited NURBS strip model.
///
/// The model is first rebuilt from its retained recipe. The recovery then
/// verifies source compatibility and every exact homogeneous block restriction.
pub fn recover_nurbs_strip_parent(model: &Json) -> Result<TensorSupport> {
    let schema = model.field("schema")?.text()?;
    if schema != SOURCE_FAMILY {
        return Err("unsupported source family for spline parent recovery");
    }
    let recipe = recipe_from_model(model)?;
    if crate::nurbs_sketch_json::build(&recipe)? != *model {
        return Err("spline parent source model does not replay exactly");
    }
    let lower = curve_json(recipe.field("lower")?)?;
    let upper = curve_json(recipe.field("upper")?)?;
    if lower.degree != upper.degree || lower.knots != upper.knots {
        return Err("spline parent boundary curves need matching degrees and knots");
    }
    if lower.controls.len() != upper.controls.len()
        || lower
            .controls
            .iter()
            .zip(&upper.controls)
            .any(|(a, b)| a[0] != b[0] || a[2] != b[2])
    {
        return Err("spline parent boundaries need a common x map and weights");
    }
    let height = recipe.field("height_mm")?.number::<i32>()?;
    if height <= 0 {
        return Err("spline parent extrusion height must be positive");
    }
    let u_count = lower.controls.len();
    let degrees = [lower.degree, 1, 1];
    let control_shape = [u_count, 2, 2];
    let domain = [
        [lower.knots[lower.degree], lower.knots[u_count]],
        [0, 1],
        [0, 1],
    ];
    let fingerprint = fingerprint(&recipe);
    let id = format!("nurbs-sketch-strip-parent-v1:{fingerprint}");
    let mut homogeneous_controls = Vec::with_capacity(u_count * 4);
    for u in 0..u_count {
        for v in 0..2 {
            let point = if v == 0 {
                lower.controls[u]
            } else {
                upper.controls[u]
            };
            for w in 0..2 {
                homogeneous_controls.push([
                    i128::from(point[0]),
                    i128::from(point[1]),
                    i128::from(point[2])
                        .checked_mul(i128::from(height) * i128::from(w))
                        .ok_or("spline parent homogeneous z overflow")?,
                    i128::from(point[2]),
                ]);
            }
        }
    }
    let source_spans = source_restrictions(model, &lower, &upper, &id, height)?;
    Ok(TensorSupport {
        id,
        source_revision: fingerprint,
        degrees,
        knots: [lower.knots.clone(), vec![0, 0, 1, 1], vec![0, 0, 1, 1]],
        domain,
        control_shape,
        homogeneous_controls,
        source_spans,
        continuity_evidence: continuity_evidence(&lower),
        generating_recipe: recipe,
        certificate_references: vec![
            "nurbs_sketch_strip_v1 exact strip separation and positive Jacobian audit".into(),
            "exact knot insertion and homogeneous restriction comparison".into(),
        ],
    })
}

/// A complete single affine cell is already a tensor parent on a unit cube.
/// Verify the affine parallelepiped on every control. Positive orientation
/// and equality with the complete material come from the native model audit.
pub fn recover_complete_parent(model: &Json) -> Result<TensorSupport> {
    if model.field("schema")?.text()? == SOURCE_FAMILY {
        return recover_nurbs_strip_parent(model);
    }
    let audit = crate::interchange::audit(model)?;
    if audit.field("status")?.text()? != "valid" {
        return Err("source geometry is not valid");
    }
    let cells = bracket_json::decode_blocks(model)?;
    let [cell] = cells.as_slice() else {
        return Err("unsupported: no complete parent witness for this material");
    };
    if cell.degrees != [1; 3] || cell.controls.len() != 8 {
        return Err("unsupported: single cell is outside affine parent recovery");
    }
    let net: Vec<_> = cell.controls.iter().map(|p| p.map(i128::from)).collect();
    let w = net[0][3];
    if w <= 0 || net.iter().any(|h| h[3] != w) {
        return Err("unsupported: affine parent requires constant positive weights");
    }
    for (i, h) in net.iter().enumerate() {
        for c in 0..3 {
            let affine = [4usize, 2, 1]
                .iter()
                .filter(|&&k| i & k != 0)
                .fold(net[0][c], |x, &k| x + net[k][c] - net[0][c]);
            if affine != h[c] {
                return Err("unsupported: trilinear control net lacks affine parent witness");
            }
        }
    }
    let edge = |i: usize, j: usize| Ratio::new(net[i][j] - net[0][j], w);
    let mut det = Ratio::zero();
    for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
        let minor = edge(2, j)?
            .checked_mul(edge(1, k)?)?
            .checked_sub(edge(2, k)?.checked_mul(edge(1, j)?)?)?;
        det = det.checked_add(edge(4, i)?.checked_mul(minor)?)?;
    }
    if det.numerator() <= 0 {
        return Err("unsupported: affine parent lacks a positive orientation");
    }
    let revision = fingerprint(model);
    let id = format!("affine-cell-parent-v1:{revision}");
    let intervals = [[0, 1]; 3];
    Ok(TensorSupport {
        id: id.clone(),
        source_revision: revision,
        degrees: [1; 3],
        knots: std::array::from_fn(|_| vec![0, 0, 1, 1]),
        domain: intervals,
        control_shape: [2; 3],
        homogeneous_controls: net,
        source_spans: vec![SourceSpan {
            source_block_id: cell.id.clone(),
            restriction: ParametricRestriction {
                parent_id: id,
                intervals,
                axes: [0, 1, 2],
                reversals: [false; 3],
            },
            gauge: [1, 1],
            generating_operation: "single audited affine source cell".into(),
            event_provenance: cell.id.clone(),
        }],
        continuity_evidence: Vec::new(),
        generating_recipe: model.clone(),
        certificate_references: vec![
            "source model audit and complete affine control-net identity".into(),
        ],
    })
}

pub(crate) fn initial_mesh_admission(
    support: &TensorSupport,
) -> std::result::Result<(), &'static str> {
    for (a, knots) in support.knots.iter().enumerate() {
        let [lo, hi] = support.domain[a];
        let mut breaks = vec![lo];
        for &k in knots.iter().filter(|&&k| lo < k && k < hi) {
            if breaks.last() == Some(&k) {
                return Err("repeated source interior knot");
            }
            breaks.push(k);
        }
        breaks.push(hi);
        if breaks.len() > 2 && support.degrees[a] < 3 {
            return Err("degree elevation has no proof of cubic simple-knot membership");
        }
        let widths: Vec<_> = breaks
            .windows(2)
            .map(|v| i128::from(v[1]) - i128::from(v[0]))
            .collect();
        if widths.iter().any(|&w| w != widths[0]) {
            return Err("source knot intervals are nonuniform under affine scaling");
        }
    }
    Ok(())
}
fn rational_json(value: i128) -> Json {
    object([
        ("numerator", string(value.to_string())),
        ("denominator", string("1")),
    ])
}

fn restriction_json(restriction: &ParametricRestriction) -> Json {
    object([
        ("parent_id", string(restriction.parent_id.clone())),
        (
            "intervals",
            arr(restriction.intervals, |interval| arr(interval, number)),
        ),
        ("axes", arr(restriction.axes, number)),
        ("reversals", arr(restriction.reversals, Json::Bool)),
        (
            "coordinate_map",
            string("source = interval_start + (reversed ? 1-local : local) * interval_length"),
        ),
    ])
}

/// Serialize exact support data without converting coefficients through binary64.
pub fn serialize_support(support: &TensorSupport) -> Json {
    object([
        ("schema", string(SUPPORT_SCHEMA)),
        ("id", string(support.id.clone())),
        ("source_revision", string(support.source_revision.clone())),
        ("degrees", arr(support.degrees, number)),
        ("knots", arr(support.knots.iter(), |axis| arr(axis, number))),
        (
            "domain",
            arr(support.domain, |interval| arr(interval, number)),
        ),
        ("control_shape", arr(support.control_shape, number)),
        (
            "homogeneous_controls",
            arr(support.homogeneous_controls.iter(), |control| {
                arr(*control, rational_json)
            }),
        ),
        (
            "source_spans",
            arr(&support.source_spans, |span| {
                object([
                    ("source_block_id", string(span.source_block_id.clone())),
                    ("restriction", restriction_json(&span.restriction)),
                    ("gauge", arr(span.gauge, |value| string(value.to_string()))),
                    (
                        "generating_operation",
                        string(span.generating_operation.clone()),
                    ),
                    ("event_provenance", string(span.event_provenance.clone())),
                ])
            }),
        ),
        (
            "continuity_evidence",
            arr(&support.continuity_evidence, |evidence| {
                object([
                    ("interface", string(evidence.interface.clone())),
                    (
                        "coordinate_frame",
                        string(evidence.coordinate_frame.clone()),
                    ),
                    ("direction", number(evidence.direction)),
                    ("knot", number(evidence.knot)),
                    ("multiplicity", number(evidence.multiplicity)),
                    ("claimed_order", number(evidence.claimed_order)),
                    ("method", string(evidence.method.clone())),
                    (
                        "assumptions",
                        arr(&evidence.assumptions, |item| string(item.clone())),
                    ),
                    ("status", string(evidence.status.clone())),
                ])
            }),
        ),
        ("generating_recipe", support.generating_recipe.clone()),
        (
            "certificate_references",
            arr(&support.certificate_references, |item| string(item.clone())),
        ),
    ])
}

/// Return a compact exact audit record for the recovered source parent.
pub fn audit_parent(support: &TensorSupport) -> Json {
    let source_control_occurrences = support.source_spans.len() * (support.degrees[0] + 1) * 4;
    let c0_unique_source_controls =
        source_control_occurrences.saturating_sub(support.source_spans.len().saturating_sub(1) * 4);
    object([
        ("status", string("valid")),
        ("method", string("exact_knot_insertion_parent_restriction")),
        ("parent_id", string(support.id.clone())),
        ("degrees", arr(support.degrees, number)),
        ("control_shape", arr(support.control_shape, number)),
        ("control_count", number(support.homogeneous_controls.len())),
        ("source_spans", number(support.source_spans.len())),
        (
            "source_control_occurrences",
            number(source_control_occurrences),
        ),
        (
            "c0_unique_source_controls",
            number(c0_unique_source_controls),
        ),
        (
            "inherited_continuity_orders",
            arr(&support.continuity_evidence, |item| {
                number(item.claimed_order)
            }),
        ),
        (
            "source_domain",
            arr(support.domain, |interval| arr(interval, number)),
        ),
    ])
}

fn exact_interface_record(interface: &Json) -> Result<Json> {
    let first = interface.field("first")?.items()?;
    let second = interface.field("second")?.items()?;
    if first.len() != 2 || second.len() != 2 {
        return Err("invalid shared-face interface record");
    }
    Ok(object([
        ("incidence", interface.clone()),
        ("role", string("retained_shared_face")),
        ("established_continuity", string("C0_exact_position_match")),
        ("higher_order_continuity", string("not_certified")),
        (
            "reason",
            string("source interface lacks an exact common-parent restriction certificate"),
        ),
    ]))
}

fn known_curve_sources(model: &Json, schema: &str) -> Result<Vec<Json>> {
    if schema != "nurbs_cross_bores_v1" {
        return Ok(Vec::new());
    }
    let profiles = model.field("profiles")?;
    ["x_min", "x_max", "z_min", "z_max"]
        .into_iter()
        .map(|name| {
            let curve = profiles.field(name)?;
            let degree = curve.field("degree")?.number::<usize>()?;
            let knots = curve
                .field("knots")?
                .items()?
                .iter()
                .map(Json::number::<i64>)
                .collect::<Result<Vec<_>>>()?;
            let controls = curve.field("controls")?.items()?.len();
            let start = knots[degree];
            let end = knots[controls];
            let spans = (degree..controls)
                .filter(|&index| knots[index] < knots[index + 1])
                .map(|index| [knots[index], knots[index + 1]])
                .collect::<Vec<_>>();
            Ok(object([
                ("source_id", string(name)),
                ("degree", number(degree)),
                ("knots", arr(knots, number)),
                ("domain", arr([start, end], number)),
                (
                    "source_span_intervals",
                    arr(spans, |pair| arr(pair, number)),
                ),
                (
                    "role",
                    string("generating sketch curve; not a complete 3D material parent"),
                ),
            ]))
        })
        .collect()
}

/// Produce a read-only conversion-eligibility report while preserving source validity.
///
/// A valid source model can be outside the first single-parent T-spline profile.
/// Every exact shared-face incidence is listed, and higher continuity is never
/// inferred from geometric coincidence alone.
pub fn eligibility_report(model: &Json) -> Result<Json> {
    let audit = crate::interchange::audit(model)?;
    let schema = model.field("schema")?.text()?;
    let single = model.field("blocks")?.items()?.len() == 1;
    let candidate = if schema == SOURCE_FAMILY {
        Some(recover_nurbs_strip_parent(model)?)
    } else if single {
        recover_complete_parent(model).ok()
    } else {
        None
    };
    if let Some(support) = candidate {
        let mesh = initial_mesh_admission(&support);
        let reason = mesh
            .err()
            .unwrap_or("one exact complete parent qualifies for the initial uniform tensor mesh");
        return Ok(object([
            ("schema", string("spline_eligibility_report_v1")),
            ("profile", string("structured_tricubic_m2_v1")),
            ("source_schema", string(schema)),
            ("source_validity", audit.field("status")?.clone()),
            (
                "conversion_status",
                string(if mesh.is_ok() {
                    "eligible_single_complete_parent"
                } else {
                    "unsupported_profile"
                }),
            ),
            ("parent_recovery", string("verified")),
            ("parent_support", serialize_support(&support)),
            (
                "retained_joins",
                arr(
                    support
                        .continuity_evidence
                        .iter()
                        .filter(|e| e.claimed_order == 0),
                    |e| {
                        object([
                            ("source_interface", string(&e.interface)),
                            ("basis", string("source_knot_multiplicity_retains_c0")),
                        ])
                    },
                ),
            ),
            ("reason", string(reason)),
        ]));
    }
    let interfaces = model
        .field("topology")?
        .field("interfaces")?
        .items()?
        .iter()
        .map(exact_interface_record)
        .collect::<Result<Vec<_>>>()?;
    let known_sources = known_curve_sources(model, schema)?;
    let common_parent_obstruction = if matches!(
        schema,
        "box_cylinder_boolean_v1" | "box_cylinder_boolean_suite_v1"
    ) {
        let stock = model.field("stock")?.field("max_mm")?.items()?;
        let stock_min = model.field("stock")?.field("min_mm")?.items()?;
        let tool = model.field("tool")?;
        let [z_min, z_max] = tool.field("interval_mm")?.items()? else {
            return Err("cylinder interval must have two endpoints");
        };
        let [_, _, stock_z_max] = stock else {
            return Err("box max point must have three coordinates");
        };
        let stock_z_min = stock_min[2].number::<i64>()?;
        let radius = tool.field("radius_mm")?.number::<i64>()?;
        let [center_x, center_y] = tool.field("center_mm")?.items()? else {
            return Err("cylinder center must have two coordinates");
        };
        let center_x = center_x.number::<i64>()?;
        let center_y = center_y.number::<i64>()?;
        let recognized_tool = tool.field("kind")?.text()? == "cylinder"
            && tool.field("axis")?.text()? == "z"
            && matches!(
                model.field("operation")?.text()?,
                "difference" | "a_minus_b"
            );
        let contained_in_xy = stock_min[0].number::<i64>()? < center_x - radius
            && center_x + radius < stock[0].number::<i64>()?
            && stock_min[1].number::<i64>()? < center_y - radius
            && center_y + radius < stock[1].number::<i64>()?;
        let through = z_min.number::<i64>()? < stock_z_min
            && z_max.number::<i64>()? > stock_z_max.number::<i64>()?
            && radius > 0
            && recognized_tool
            && contained_in_xy;
        if through {
            object([
                ("kind", string("exact_through_tunnel_topology")),
                (
                    "evidence",
                    string(
                        "positive-radius cylindrical tool strictly crosses both stock end planes",
                    ),
                ),
                (
                    "consequence",
                    string(
                        "material domain is not homeomorphic to one injective tensor-product cube",
                    ),
                ),
            ])
        } else {
            object([
                ("kind", string("not_demonstrated")),
                (
                    "evidence",
                    string("through-tunnel preconditions were not all established"),
                ),
                (
                    "consequence",
                    string("common-parent topology remains unclassified"),
                ),
            ])
        }
    } else if schema == "nurbs_cross_bores_v1" {
        object([
            ("kind", string("multiple_generating_supports")),
            (
                "evidence",
                string("four source profiles generate a multi-cell cross construction"),
            ),
            (
                "consequence",
                string("no exact complete 3D tensor parent restriction is recorded for every material block"),
            ),
        ])
    } else {
        object([
            ("kind", string("no_complete_parent_witness")),
            (
                "evidence",
                string("no complete source restriction proof in this recovery family"),
            ),
            (
                "consequence",
                string("global topology and mathematical representability remain unclassified"),
            ),
        ])
    };
    Ok(object([
        ("schema", string("spline_eligibility_report_v1")),
        ("profile", string("structured_tricubic_m2_v1")),
        ("source_schema", string(schema)),
        ("source_validity", audit.field("status")?.clone()),
        ("conversion_status", string("unsupported_profile")),
        ("source_model_preserved", Json::Bool(true)),
        ("known_generating_sources", Json::Array(known_sources)),
        (
            "known_source_restriction_status",
            string("recipe-level supports retained; complete per-block 3D parent restrictions unavailable"),
        ),
        ("retained_joins", Json::Array(interfaces)),
        ("regional_source_columns", recover_translation_columns(model, 2)?),
        ("common_parent_obstruction", common_parent_obstruction),
        (
            "reason",
            string("the first profile requires one regular injective parent covering the entire material body"),
        ),
    ]))
}

/// Recover exact translation columns for blocks extruded along any parametric axis.
///
/// For a degree-one local axis, equal endpoint weights and the identity
/// `H_1 = H_0 + t W` on every transverse control prove `X(q)=X(0)+q t` over
/// the complete block. The returned source net is the exact `q=0` restriction.
pub fn recover_translation_columns(model: &Json, local_axis: usize) -> Result<Json> {
    if local_axis >= 3 {
        return Err("translation-column axis is outside 0..3");
    }
    let blocks = bracket_json::decode_blocks(model)?;
    let mut reports = Vec::with_capacity(blocks.len());
    let mut groups: std::collections::BTreeMap<[Ratio; 3], Vec<String>> =
        std::collections::BTreeMap::new();
    for block in blocks {
        let shape = block.degrees.map(|degree| degree + 1);
        if block.degrees[local_axis] != 1 || shape.iter().product::<usize>() != block.controls.len()
        {
            reports.push(object([
                ("block_id", string(block.id)),
                ("status", string("unsupported_column_axis")),
                (
                    "reason",
                    string("selected parametric axis is not a valid degree-one extrusion"),
                ),
            ]));
            continue;
        }
        let mut transverse_shape = shape;
        transverse_shape[local_axis] = 1;
        let mut translation: Option<[Ratio; 3]> = None;
        let mut base_controls = Vec::new();
        let mut control_indices = Vec::new();
        let mut rejection = None;
        'transverse: for i in 0..transverse_shape[0] {
            for j in 0..transverse_shape[1] {
                for k in 0..transverse_shape[2] {
                    let low_index = [i, j, k];
                    let mut high_index = low_index;
                    high_index[local_axis] = 1;
                    let flat = |[a, b, c]: [usize; 3]| (a * shape[1] + b) * shape[2] + c;
                    let low_flat = flat(low_index);
                    let high_flat = flat(high_index);
                    let low = block.controls[low_flat];
                    let high = block.controls[high_flat];
                    if low[3] <= 0 || high[3] != low[3] {
                        rejection = Some("endpoint weights are not identical and positive");
                        break 'transverse;
                    }
                    let candidate: [Result<Ratio>; 3] = std::array::from_fn(|component| {
                        Ratio::new(
                            i128::from(high[component]) - i128::from(low[component]),
                            i128::from(low[3]),
                        )
                    });
                    let candidate = [candidate[0]?, candidate[1]?, candidate[2]?];
                    if translation.is_some_and(|previous| previous != candidate) {
                        rejection =
                            Some("homogeneous endpoint differences are not one translation vector");
                        break 'transverse;
                    }
                    translation = Some(candidate);
                    base_controls.push(low);
                    control_indices.push(low_flat);
                }
            }
        }
        if let Some(reason) = rejection {
            reports.push(object([
                ("block_id", string(block.id)),
                ("status", string("not_a_translation_column")),
                ("parametric_axis", number(local_axis)),
                ("reason", string(reason)),
            ]));
            continue;
        }
        let translation = translation.ok_or("empty translation-column control net")?;
        groups
            .entry(translation)
            .or_default()
            .push(block.id.clone());
        reports.push(object([
            ("block_id", string(block.id)),
            ("status", string("certified_exact_prismatic_column")),
            ("parametric_axis", number(local_axis)),
            ("degrees", arr(block.degrees, number)),
            ("source_column_shape", arr(transverse_shape, number)),
            ("source_control_indices", arr(control_indices, number)),
            (
                "source_homogeneous_controls",
                arr(base_controls, |control| arr(control, number)),
            ),
            (
                "translation_vector_mm",
                arr(translation, |component| {
                    object([
                        ("numerator", string(component.numerator().to_string())),
                        ("denominator", string(component.denominator().to_string())),
                    ])
                }),
            ),
            (
                "proof",
                string("all transverse controls satisfy equal positive weights and H1-H0=t*W"),
            ),
        ]));
    }
    let grouped = groups
        .into_iter()
        .map(|(translation, block_ids)| {
            object([
                (
                    "translation_vector_mm",
                    arr(translation, |component| {
                        object([
                            ("numerator", string(component.numerator().to_string())),
                            ("denominator", string(component.denominator().to_string())),
                        ])
                    }),
                ),
                ("block_ids", arr(block_ids, string)),
                (
                    "role",
                    string("parallel translation direction; distinct source columns remain distinct geometry"),
                ),
            ])
        })
        .collect();
    Ok(object([
        ("schema", string("exact_translation_columns_v1")),
        ("parametric_axis", number(local_axis)),
        ("blocks", Json::Array(reports)),
        ("parallel_column_groups", Json::Array(grouped)),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json::parse;

    fn fixture_model() -> Json {
        let recipe: Json =
            parse(&std::fs::read_to_string("tests/references/nurbs_strip.recipe.json").unwrap())
                .unwrap();
        crate::nurbs_sketch_json::build(&recipe).unwrap()
    }

    #[test]
    fn recovers_exact_twenty_control_parent_and_c2_source_interface() {
        let support = recover_nurbs_strip_parent(&fixture_model()).unwrap();
        assert_eq!(support.degrees, [3, 1, 1]);
        assert_eq!(support.control_shape, [5, 2, 2]);
        assert_eq!(support.homogeneous_controls.len(), 20);
        assert_eq!(support.domain, [[-2, 4], [0, 1], [0, 1]]);
        assert_eq!(support.source_spans.len(), 2);
        assert_eq!(support.continuity_evidence.len(), 1);
        assert_eq!(support.continuity_evidence[0].knot, 1);
        assert_eq!(support.continuity_evidence[0].multiplicity, 1);
        assert_eq!(support.continuity_evidence[0].claimed_order, 2);
        assert_eq!(support.source_spans[0].restriction.intervals[0], [-2, 1]);
        assert_eq!(support.source_spans[1].restriction.intervals[0], [1, 4]);
    }

    #[test]
    fn serializes_exact_integer_controls_as_canonical_rationals() {
        let support = recover_nurbs_strip_parent(&fixture_model()).unwrap();
        let encoded = serialize_support(&support);
        assert_eq!(
            encoded.field("schema").unwrap().text().unwrap(),
            SUPPORT_SCHEMA
        );
        let first = &encoded
            .field("homogeneous_controls")
            .unwrap()
            .items()
            .unwrap()[0];
        let first_component = &first.items().unwrap()[0];
        assert_eq!(
            first_component.field("numerator").unwrap().text().unwrap(),
            "0"
        );
        assert_eq!(
            first_component
                .field("denominator")
                .unwrap()
                .text()
                .unwrap(),
            "1"
        );
        assert_eq!(serialize_support(&support), encoded);
    }

    #[test]
    fn restriction_mapping_preserves_permutations_reversals_and_actual_metric() {
        let restriction = ParametricRestriction {
            parent_id: "parent".into(),
            intervals: [[-2, 4], [10, 14], [20, 22]],
            axes: [2, 0, 1],
            reversals: [true, false, true],
        };
        let mapped = restriction
            .map_parameters([
                Ratio::new(1, 4).unwrap(),
                Ratio::new(1, 2).unwrap(),
                Ratio::new(3, 4).unwrap(),
            ])
            .unwrap();
        assert_eq!(mapped[0], Ratio::new(1, 1).unwrap());
        assert_eq!(mapped[1], Ratio::new(11, 1).unwrap());
        assert_eq!(mapped[2], Ratio::new(43, 2).unwrap());
        assert_eq!(
            restriction.derivative_scale(0).unwrap(),
            Ratio::new(-2, 1).unwrap()
        );
        assert_eq!(
            restriction.derivative_scale(1).unwrap(),
            Ratio::new(6, 1).unwrap()
        );
        assert_eq!(
            restriction.derivative_scale(2).unwrap(),
            Ratio::new(-4, 1).unwrap()
        );
    }

    #[test]
    fn forged_block_coefficients_and_incompatible_curve_supports_are_rejected() {
        let mut model = fixture_model();
        let Json::Array(blocks) = model.field("blocks").unwrap().clone() else {
            panic!("blocks are an array");
        };
        let mut changed_blocks = blocks;
        let Json::Array(controls) = changed_blocks[0].field("controls").unwrap().clone() else {
            panic!("controls are an array");
        };
        let mut changed_controls = controls;
        changed_controls[0] = Json::Array(vec![number(1), number(0), number(0), number(1)]);
        let mut block_object = match changed_blocks[0].clone() {
            Json::Object(object) => object,
            _ => unreachable!(),
        };
        block_object.insert("controls".into(), Json::Array(changed_controls));
        changed_blocks[0] = Json::Object(block_object);
        let mut model_object = match model {
            Json::Object(object) => object,
            _ => unreachable!(),
        };
        model_object.insert("blocks".into(), Json::Array(changed_blocks));
        model = Json::Object(model_object);
        assert!(recover_nurbs_strip_parent(&model).is_err());
    }

    #[test]
    fn actual_nonuniform_knots_are_not_normalized_away() {
        let support = recover_nurbs_strip_parent(&fixture_model()).unwrap();
        assert_eq!(support.knots[0], [-2, -2, -2, -2, 1, 4, 4, 4, 4]);
        let first_span = &support.source_spans[0].restriction;
        let mapped = first_span
            .map_parameters([
                Ratio::new(1, 2).unwrap(),
                Ratio::new(1, 2).unwrap(),
                Ratio::new(1, 2).unwrap(),
            ])
            .unwrap();
        assert_eq!(mapped[0], Ratio::new(-1, 2).unwrap());
    }

    #[test]
    fn nonuniform_source_intervals_keep_signed_metric_scales() {
        let mut recipe: Json =
            parse(&std::fs::read_to_string("tests/references/nurbs_strip.recipe.json").unwrap())
                .unwrap();
        let mut recipe_object = match recipe.clone() {
            Json::Object(object) => object,
            _ => unreachable!(),
        };
        for key in ["lower", "upper"] {
            let mut curve = match recipe_object[key].clone() {
                Json::Object(object) => object,
                _ => unreachable!(),
            };
            curve.insert(
                "knots".into(),
                Json::Array(
                    [-5, -5, -5, -5, 0, 4, 4, 4, 4]
                        .into_iter()
                        .map(number)
                        .collect(),
                ),
            );
            recipe_object.insert(key.into(), Json::Object(curve));
        }
        recipe = Json::Object(recipe_object);
        let model = crate::nurbs_sketch_json::build(&recipe).unwrap();
        let support = recover_nurbs_strip_parent(&model).unwrap();
        assert_eq!(support.domain[0], [-5, 4]);
        assert_eq!(support.source_spans[0].restriction.intervals[0], [-5, 0]);
        assert_eq!(support.source_spans[1].restriction.intervals[0], [0, 4]);
        assert_eq!(
            support.source_spans[0]
                .restriction
                .derivative_scale(0)
                .unwrap(),
            Ratio::new(5, 1).unwrap()
        );
        assert_eq!(
            support.source_spans[1]
                .restriction
                .derivative_scale(0)
                .unwrap(),
            Ratio::new(4, 1).unwrap()
        );
    }

    #[test]
    fn reversed_source_curve_records_a_signed_span_restriction() {
        let mut recipe: Json =
            parse(&std::fs::read_to_string("tests/references/nurbs_strip.recipe.json").unwrap())
                .unwrap();
        let mut recipe_object = match recipe.clone() {
            Json::Object(object) => object,
            _ => unreachable!(),
        };
        for key in ["lower", "upper"] {
            let mut curve = match recipe_object[key].clone() {
                Json::Object(object) => object,
                _ => unreachable!(),
            };
            let controls = curve["controls"]
                .items()
                .unwrap()
                .iter()
                .cloned()
                .rev()
                .collect();
            curve.insert("controls".into(), Json::Array(controls));
            recipe_object.insert(key.into(), Json::Object(curve));
        }
        recipe = Json::Object(recipe_object);
        let model = crate::nurbs_sketch_json::build(&recipe).unwrap();
        let support = recover_nurbs_strip_parent(&model).unwrap();
        assert!(support
            .source_spans
            .iter()
            .all(|span| span.restriction.reversals[0]));
        assert_eq!(support.source_spans[0].restriction.intervals[0], [1, 4]);
        assert_eq!(
            support.source_spans[0]
                .restriction
                .derivative_scale(0)
                .unwrap(),
            Ratio::new(-3, 1).unwrap()
        );
    }

    #[test]
    fn support_revision_changes_when_recipe_geometry_changes() {
        let model = fixture_model();
        let first = recover_nurbs_strip_parent(&model).unwrap();
        let mut recipe: Json =
            parse(&std::fs::read_to_string("tests/references/nurbs_strip.recipe.json").unwrap())
                .unwrap();
        let mut recipe_object = match recipe.clone() {
            Json::Object(object) => object,
            _ => unreachable!(),
        };
        let mut lower = match recipe_object["lower"].clone() {
            Json::Object(object) => object,
            _ => unreachable!(),
        };
        let mut controls = lower["controls"].items().unwrap().to_vec();
        let Json::Array(mut first_control) = controls[1].clone() else {
            panic!("source control is an array");
        };
        first_control[1] = number(1);
        controls[1] = Json::Array(first_control);
        lower.insert("controls".into(), Json::Array(controls));
        recipe_object.insert("lower".into(), Json::Object(lower));
        recipe = Json::Object(recipe_object);
        let changed_model = crate::nurbs_sketch_json::build(&recipe).unwrap();
        let changed = recover_nurbs_strip_parent(&changed_model).unwrap();
        assert_ne!(first.id, changed.id);
        assert_ne!(first.source_revision, changed.source_revision);
    }

    fn recipe_model(path: &str) -> Json {
        let recipe: Json = parse(&std::fs::read_to_string(path).unwrap()).unwrap();
        crate::interchange::build(&recipe).unwrap()
    }

    #[test]
    fn eligibility_reports_preserve_box_cylinder_validity_and_all_four_joins() {
        let report = eligibility_report(&recipe_model(
            "tests/references/box_minus_cylinder.recipe.json",
        ))
        .unwrap();
        assert_eq!(
            report.field("conversion_status").unwrap().text().unwrap(),
            "unsupported_profile"
        );
        assert_eq!(
            report.field("source_validity").unwrap().text().unwrap(),
            "valid"
        );
        assert_eq!(
            report.field("source_model_preserved").unwrap(),
            &Json::Bool(true)
        );
        assert_eq!(
            report
                .field("retained_joins")
                .unwrap()
                .items()
                .unwrap()
                .len(),
            4
        );
        assert_eq!(
            report
                .field("common_parent_obstruction")
                .unwrap()
                .field("kind")
                .unwrap()
                .text()
                .unwrap(),
            "exact_through_tunnel_topology"
        );
    }

    #[test]
    fn eligibility_reports_all_nurbs_cross_joins_and_recipe_sources() {
        let report = eligibility_report(&recipe_model(
            "tests/references/nurbs_cross_bores.recipe.json",
        ))
        .unwrap();
        assert_eq!(
            report.field("conversion_status").unwrap().text().unwrap(),
            "unsupported_profile"
        );
        assert_eq!(
            report.field("source_validity").unwrap().text().unwrap(),
            "valid"
        );
        assert_eq!(
            report
                .field("known_generating_sources")
                .unwrap()
                .items()
                .unwrap()
                .len(),
            4
        );
        assert_eq!(
            report
                .field("retained_joins")
                .unwrap()
                .items()
                .unwrap()
                .len(),
            36
        );
    }

    #[test]
    fn source_columns_are_recovered_on_real_cylinder_and_blind_pocket_models() {
        let cylinder = recipe_model("tests/references/box_minus_cylinder.recipe.json");
        let cylinder_columns = recover_translation_columns(&cylinder, 2).unwrap();
        assert_eq!(
            cylinder_columns
                .field("blocks")
                .unwrap()
                .items()
                .unwrap()
                .len(),
            4
        );
        assert_eq!(
            cylinder_columns
                .field("parallel_column_groups")
                .unwrap()
                .items()
                .unwrap()
                .len(),
            1
        );
        let vector = cylinder_columns
            .field("parallel_column_groups")
            .unwrap()
            .items()
            .unwrap()[0]
            .field("translation_vector_mm")
            .unwrap();
        let values = vector.items().unwrap();
        assert_eq!(values[0].field("numerator").unwrap().text().unwrap(), "0");
        assert_eq!(values[1].field("numerator").unwrap().text().unwrap(), "0");
        assert_eq!(values[2].field("numerator").unwrap().text().unwrap(), "12");

        let pocket = recipe_model("tests/references/sketch_boolean_blind_pocket.recipe.json");
        let pocket_columns = recover_translation_columns(&pocket, 2).unwrap();
        assert_eq!(
            pocket_columns
                .field("blocks")
                .unwrap()
                .items()
                .unwrap()
                .len(),
            9
        );
        assert_eq!(
            pocket_columns
                .field("parallel_column_groups")
                .unwrap()
                .items()
                .unwrap()
                .len(),
            2
        );
        for item in pocket_columns.field("blocks").unwrap().items().unwrap() {
            assert_eq!(
                item.field("status").unwrap().text().unwrap(),
                "certified_exact_prismatic_column"
            );
        }
    }

    #[test]
    fn translation_column_proof_accepts_arbitrary_exact_world_vectors() {
        let mut controls = Vec::new();
        for i in 0..2_i64 {
            for j in 0..2_i64 {
                for k in 0..2_i64 {
                    controls.push(Json::Array(
                        [2 * i, 3 * j, 5 * k, 1].into_iter().map(number).collect(),
                    ));
                }
            }
        }
        for i in 0..2_i64 {
            for j in 0..2_i64 {
                let base = (i as usize * 2 + j as usize) * 2;
                let mut high = controls[base]
                    .items()
                    .unwrap()
                    .iter()
                    .map(Json::number::<i64>)
                    .collect::<Result<Vec<_>>>()
                    .unwrap();
                high[0] += 3;
                high[1] += 4;
                controls[base + 1] = Json::Array(high.into_iter().map(number).collect());
            }
        }
        let model = object([(
            "blocks",
            Json::Array(vec![object([
                ("id", string("diagonal")),
                ("degrees", arr([1, 1, 1], number)),
                ("controls", Json::Array(controls)),
            ])]),
        )]);
        let report = recover_translation_columns(&model, 2).unwrap();
        let vector = report
            .field("parallel_column_groups")
            .unwrap()
            .items()
            .unwrap()[0]
            .field("translation_vector_mm")
            .unwrap()
            .items()
            .unwrap();
        assert_eq!(vector[0].field("numerator").unwrap().text().unwrap(), "3");
        assert_eq!(vector[1].field("numerator").unwrap().text().unwrap(), "4");
        assert_eq!(vector[2].field("numerator").unwrap().text().unwrap(), "0");
    }
}
