//! Uniform section h-refinement on the polygonal blind-pocket benchmark.
//! C0 identifies oriented traces. C1 additionally matches physical gradients
//! on every open span, including both endpoints for vertex enclosure.
use crate::{
    bracket_json::{decode_blocks, decode_side_name},
    constraint_system::{resolve, Ratio as Q},
    extracted_space::{ElementBasisSample, ElementGeometrySample},
    json::{object, Json},
    section_manufactured::{solve_sampled, IntegrationDomain},
    spline_constraints::{add, conditioned_kernel, Result, Row},
    transverse_jets::bernstein,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::OnceLock,
};

type Point = [i64; 2];
type Quad = [Point; 4];
type Side = (usize, usize, usize);
type GeometrySignature = (Vec<String>, Vec<String>, Vec<String>);
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SectionJoins {
    None,
    Radial,
    Footprint,
    All,
}

/// One conservative or geometry-qualified choice for a single field solve.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AutoSelection {
    pub degree: usize,
    pub intervals: usize,
    pub joins: SectionJoins,
    pub qualified_geometry_reused: bool,
    pub reason: &'static str,
}

/// One exact radial seam in the qualified model's current identity namespace.
#[derive(Clone, Debug, PartialEq)]
pub struct QualifiedInterface {
    pub first_block: String,
    pub first_side: String,
    pub second_block: String,
    pub second_side: String,
    pub orientation: Json,
}

/// Select one solve from the measured curved-pocket family. Reuse is based on
/// exact block controls and incidence topology, not names or a sampled hash.
/// The interval count is caller-selected; no candidate field is solved here.
pub fn plan_auto(model: &Json, intervals: usize) -> Result<AutoSelection> {
    if ![1, 2, 4, 8].contains(&intervals) {
        return Err("unsupported automatic pocket interval count");
    }
    if is_qualified_geometry(model)? {
        Ok(AutoSelection {
            degree: 3,
            intervals,
            joins: SectionJoins::Radial,
            qualified_geometry_reused: true,
            reason: "exact published geometry and qualified cubic radial-join field family",
        })
    } else {
        Ok(AutoSelection {
            degree: 3,
            intervals,
            joins: SectionJoins::None,
            qualified_geometry_reused: false,
            reason:
                "geometry is outside the exact qualified model; use conservative C0 field family",
        })
    }
}

/// Exact identity predicate shared with other automatic analysis policies.
pub fn is_qualified_geometry(model: &Json) -> Result<bool> {
    if crate::interchange::audit(model)?.field("status")?.text()? != "valid" {
        return Err("source geometry is not certified");
    }
    matches_qualified_layout(model)
}

/// The caller must audit before using the selected field. The inherited
/// compiler performs that audit while constructing its space.
pub(crate) fn matches_qualified_layout(model: &Json) -> Result<bool> {
    if model.field("schema")?.text()? != "curved_sketch_layers_v1"
        || model.field("blocks")?.items()?.len() != 9
        || model.field("topology")?.field("interfaces")?.items()?.len() != 16
        || model.field("topology")?.field("exterior")?.items()?.len() != 22
    {
        return Ok(false);
    }
    let signature = match exact_geometry_signature(model) {
        Ok(signature) => signature,
        Err(_) => return Ok(false),
    };
    Ok(signature == *reference_geometry_signature()?)
}

fn reference_geometry_signature() -> Result<&'static GeometrySignature> {
    static SIGNATURE: OnceLock<Result<GeometrySignature>> = OnceLock::new();
    SIGNATURE
        .get_or_init(|| {
            let reference = crate::json::parse(include_str!(
                "../docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"
            ))?;
            exact_geometry_signature(&reference)
        })
        .as_ref()
        .map_err(|error| *error)
}

/// Return the straight lateral volume faces on the four qualified radial
/// section seams (eight block interfaces across the two axial layers). IDs are
/// resolved from the supplied model after exact geometry and topology match.
pub fn qualified_radial_interfaces(model: &Json) -> Result<Vec<QualifiedInterface>> {
    if !is_qualified_geometry(model)? {
        return Err("radial join reuse requires the exact qualified curved pocket geometry");
    }
    radial_interfaces_after_qualification(model)
}

/// Derive the radial faces after a caller has already qualified this model
/// with [`is_qualified_geometry`] or [`plan_auto`].
pub(crate) fn radial_interfaces_after_qualification(
    model: &Json,
) -> Result<Vec<QualifiedInterface>> {
    let blocks = decode_blocks(model)?;
    let mut charts = BTreeMap::new();
    for block in &blocks {
        charts.insert(block.id.clone(), chart_from_block(block)?);
    }
    let mut selected = Vec::new();
    for interface in model.field("topology")?.field("interfaces")?.items()? {
        let first = interface.field("first")?.items()?;
        let second = interface.field("second")?.items()?;
        if first.len() != 2 || second.len() != 2 {
            return Err("invalid interface side in qualified geometry");
        }
        let (first_axis, first_upper) = decode_side_name(first[1].text()?)?;
        let (second_axis, second_upper) = decode_side_name(second[1].text()?)?;
        if first_axis == 2 || second_axis == 2 {
            continue;
        }
        let first_id = first[0].text()?;
        let second_id = second[0].text()?;
        let a = charts
            .get(first_id)
            .ok_or("unknown first block in qualified interface")?;
        let b = charts
            .get(second_id)
            .ok_or("unknown second block in qualified interface")?;
        if a.edge_is_straight(first_axis, usize::from(first_upper))
            && b.edge_is_straight(second_axis, usize::from(second_upper))
        {
            selected.push(QualifiedInterface {
                first_block: first_id.to_owned(),
                first_side: first[1].text()?.to_owned(),
                second_block: second_id.to_owned(),
                second_side: second[1].text()?.to_owned(),
                orientation: interface.field("orientation")?.clone(),
            });
        }
    }
    if selected.len() != 8 {
        return Err("qualified radial interface inventory differs from eight lateral faces");
    }
    Ok(selected)
}

/// Run exactly one solve using [`plan_auto`], then attach the selection plan.
pub fn solve_model_auto(model: &Json, intervals: usize, job: &Json) -> Result<Json> {
    let plan = plan_auto(model, intervals)?;
    let mut result =
        solve_model_selected_audited(model, plan.degree, plan.intervals, plan.joins, job)?;
    if let Json::Object(fields) = &mut result {
        fields.insert(
            "selection_mode".into(),
            Json::String("automatic_single_solve".into()),
        );
        fields.insert("selection_reason".into(), Json::String(plan.reason.into()));
        fields.insert(
            "qualified_geometry_reused".into(),
            Json::Bool(plan.qualified_geometry_reused),
        );
        fields.insert("trial_solve_count".into(), Json::Number("0".into()));
        fields.insert("final_solve_count".into(), Json::Number("1".into()));
        fields.insert(
            "selection_model_block_count".into(),
            Json::Number(model.field("blocks")?.items()?.len().to_string()),
        );
        fields.insert(
            "selection_model_interface_count".into(),
            Json::Number(
                model
                    .field("topology")?
                    .field("interfaces")?
                    .items()?
                    .len()
                    .to_string(),
            ),
        );
        fields.insert(
            "selection_model_exterior_side_count".into(),
            Json::Number(
                model
                    .field("topology")?
                    .field("exterior")?
                    .items()?
                    .len()
                    .to_string(),
            ),
        );
        fields.insert(
            "selection_degree".into(),
            Json::Number(plan.degree.to_string()),
        );
        fields.insert(
            "selection_intervals".into(),
            Json::Number(plan.intervals.to_string()),
        );
        fields.insert(
            "selection_joins".into(),
            Json::String(format!("{:?}", plan.joins)),
        );
    }
    Ok(result)
}

/// Canonical exact geometry identity. Block IDs are replaced by their full
/// degree/control-net records before topology is compared. Duplicate maps are
/// rejected because their references would otherwise be ambiguous.
fn exact_geometry_signature(model: &Json) -> Result<GeometrySignature> {
    let mut blocks_by_id = Vec::new();
    for block in model.field("blocks")?.items()? {
        let id = block.field("id")?.text()?.to_owned();
        let signature = Json::Object(
            [
                ("degrees".into(), block.field("degrees")?.clone()),
                ("controls".into(), block.field("controls")?.clone()),
            ]
            .into_iter()
            .collect(),
        )
        .to_string();
        blocks_by_id.push((id, signature));
    }
    blocks_by_id.sort_by(|a, b| a.1.cmp(&b.1));
    if blocks_by_id.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Err("duplicate block identity in geometry signature");
    }
    let blocks: Vec<_> = blocks_by_id
        .iter()
        .map(|(_, signature)| signature.clone())
        .collect();
    if blocks.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err("duplicate block geometry in geometry signature");
    }
    let ids: BTreeMap<_, _> = blocks_by_id
        .into_iter()
        .enumerate()
        .map(|(index, (id, _))| (id, index))
        .collect();
    let side_signature = |side: &Json| -> Result<String> {
        let pair = side.items()?;
        if pair.len() != 2 {
            return Err("invalid side in geometry signature");
        }
        let key = ids
            .get(pair[0].text()?)
            .ok_or("unknown block in geometry signature")?;
        Ok(format!("{key}:{}", pair[1].text()?))
    };
    let topology = model.field("topology")?;
    let mut interfaces = Vec::new();
    for interface in topology.field("interfaces")?.items()? {
        let a = side_signature(interface.field("first")?)?;
        let b = side_signature(interface.field("second")?)?;
        let orientation = interface.field("orientation")?.to_string();
        interfaces.push(format!("{a}|{b}|{orientation}"));
    }
    interfaces.sort();
    let mut exterior = topology
        .field("exterior")?
        .items()?
        .iter()
        .map(side_signature)
        .collect::<Result<Vec<_>>>()?;
    exterior.sort();
    Ok((blocks, interfaces, exterior))
}
#[derive(Clone, PartialEq, Eq)]
struct Chart {
    corners: Quad,
    /// Tensor quadratic XY Bezier net, in i-major order. Bilinear source
    /// maps are exactly degree-elevated into this representation.
    xy: [[Q; 2]; 9],
}
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum AxisKind {
    Value,
    Slope,
    Middle,
}
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct AxisDof {
    kind: AxisKind,
    knot: usize,
}
fn integer(i: usize) -> Q {
    Q::new(i as i128, 1).unwrap()
}
fn number(q: Q) -> f64 {
    q.numerator() as f64 / q.denominator() as f64
}
fn endpoints(c: &Quad, a: usize, b: usize) -> [Point; 2] {
    if a == 0 {
        [c[2 * b], c[2 * b + 1]]
    } else {
        [c[b], c[2 + b]]
    }
}
impl Chart {
    fn edge_controls(&self, axis: usize, endpoint: usize) -> [[Q; 2]; 3] {
        std::array::from_fn(|k| {
            let (i, j) = if axis == 0 {
                (endpoint * 2, k)
            } else {
                (k, endpoint * 2)
            };
            self.xy[i * 3 + j]
        })
    }
    fn edge_is_straight(&self, axis: usize, endpoint: usize) -> bool {
        let points = self.edge_controls(axis, endpoint);
        let dx = points[2][0].checked_sub(points[0][0]);
        let dy = points[2][1].checked_sub(points[0][1]);
        let (Ok(dx), Ok(dy)) = (dx, dy) else {
            return false;
        };
        points[1][0]
            .checked_sub(points[0][0])
            .and_then(|x| x.checked_mul(dy))
            .and_then(|x| {
                points[1][1]
                    .checked_sub(points[0][1])
                    .and_then(|y| y.checked_mul(dx).and_then(|y| x.checked_sub(y)))
            })
            .is_ok_and(|x| x.numerator() == 0)
    }
    fn derivatives(&self, u: Q, v: Q) -> Result<[[Q; 2]; 2]> {
        let (bu, du) = bernstein(2, u)?;
        let (bv, dv) = bernstein(2, v)?;
        let mut out = [[Q::zero(); 2]; 2];
        for i in 0..3 {
            for j in 0..3 {
                let p = self.xy[i * 3 + j];
                for k in 0..2 {
                    out[0][k] =
                        out[0][k].checked_add(du[i].checked_mul(bv[j])?.checked_mul(p[k])?)?;
                    out[1][k] =
                        out[1][k].checked_add(bu[i].checked_mul(dv[j])?.checked_mul(p[k])?)?;
                }
            }
        }
        Ok(out)
    }
    fn evaluate(&self, u: f64, v: f64) -> ([f64; 2], [[f64; 2]; 2]) {
        let bu = [(1. - u).powi(2), 2. * u * (1. - u), u.powi(2)];
        let bv = [(1. - v).powi(2), 2. * v * (1. - v), v.powi(2)];
        let du = [-2. * (1. - u), 2. - 4. * u, 2. * u];
        let dv = [-2. * (1. - v), 2. - 4. * v, 2. * v];
        let mut position = [0.; 2];
        let mut jac = [[0.; 2]; 2];
        for i in 0..3 {
            for j in 0..3 {
                let p = self.xy[i * 3 + j].map(number);
                for k in 0..2 {
                    position[k] += bu[i] * bv[j] * p[k];
                    jac[0][k] += du[i] * bv[j] * p[k];
                    jac[1][k] += bu[i] * dv[j] * p[k];
                }
            }
        }
        (position, jac)
    }
}
fn chart_from_block(b: &crate::bracket::BracketBlock) -> Result<Chart> {
    if b.degrees != [1, 1, 1] && b.degrees != [2, 2, 1] {
        return Err("unsupported curved pocket chart degree");
    }
    let expected = if b.degrees == [1, 1, 1] { 8 } else { 18 };
    if b.controls.len() != expected {
        return Err("invalid curved pocket control count");
    }
    let weight = b.controls[0][3];
    if weight <= 0 || b.controls.iter().any(|h| h[3] != weight) {
        return Err("curved pocket chart requires one positive common weight");
    }
    let point = |h: &[i64; 4]| -> Result<[Q; 2]> {
        Ok([
            Q::new(i128::from(h[0]), i128::from(h[3]))?,
            Q::new(i128::from(h[1]), i128::from(h[3]))?,
        ])
    };
    let mut xy = [[Q::zero(); 2]; 9];
    if b.degrees == [2, 2, 1] {
        for i in 0..3 {
            for j in 0..3 {
                let lo = &b.controls[(i * 3 + j) * 2];
                let hi = &b.controls[(i * 3 + j) * 2 + 1];
                if lo[0] != hi[0] || lo[1] != hi[1] || lo[3] != hi[3] {
                    return Err("curved pocket XY chart changes along z");
                }
                xy[i * 3 + j] = point(lo)?;
            }
        }
    } else {
        let mut linear = [[Q::zero(); 2]; 4];
        for i in 0..2 {
            for j in 0..2 {
                let lo = &b.controls[(i * 2 + j) * 2];
                let hi = &b.controls[(i * 2 + j) * 2 + 1];
                if lo[0] != hi[0] || lo[1] != hi[1] || lo[3] != hi[3] {
                    return Err("bilinear pocket XY chart changes along z");
                }
                linear[i * 2 + j] = point(lo)?;
            }
        }
        let elevate = [
            [Q::one(), Q::zero()],
            [Q::new(1, 2)?, Q::new(1, 2)?],
            [Q::zero(), Q::one()],
        ];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..2 {
                    for a in 0..2 {
                        for d in 0..2 {
                            xy[i * 3 + j][k] = xy[i * 3 + j][k].checked_add(
                                elevate[i][a]
                                    .checked_mul(elevate[j][d])?
                                    .checked_mul(linear[a * 2 + d][k])?,
                            )?;
                        }
                    }
                }
            }
        }
    }
    let integer = |x: Q| -> Result<i64> {
        if x.denominator() != 1 {
            return Err("curved pocket corner is not an integer topology point");
        }
        i64::try_from(x.numerator()).map_err(|_| "curved pocket corner overflows")
    };
    let indices = [0, 2, 6, 8];
    let corners = [
        [integer(xy[indices[0]][0])?, integer(xy[indices[0]][1])?],
        [integer(xy[indices[1]][0])?, integer(xy[indices[1]][1])?],
        [integer(xy[indices[2]][0])?, integer(xy[indices[2]][1])?],
        [integer(xy[indices[3]][0])?, integer(xy[indices[3]][1])?],
    ];
    Ok(Chart { corners, xy })
}
fn axis_q(rows: &[Row], t: Q) -> Result<BTreeMap<usize, [Q; 2]>> {
    let (b, d) = bernstein(rows.len() - 1, t)?;
    let mut out = BTreeMap::new();
    for (k, row) in rows.iter().enumerate() {
        for (&i, &c) in row {
            let v = out.entry(i).or_insert([Q::zero(); 2]);
            v[0] = v[0].checked_add(b[k].checked_mul(c)?)?;
            v[1] = v[1].checked_add(d[k].checked_mul(c)?)?;
        }
    }
    Ok(out)
}
struct Basis {
    n: usize,
    dim: usize,
    charts: Vec<Chart>,
    ids: Vec<Vec<(usize, Q)>>,
    spans: Vec<Vec<Row>>,
}
impl Basis {
    fn value_trace(&self, side: Side, k: usize, t: Q, sign: Q, row: &mut Row) -> Result<()> {
        let (chart, axis, endpoint) = side;
        let mut cells = [k; 2];
        cells[axis] = endpoint * (self.n - 1);
        let mut pt = [t; 2];
        pt[axis] = integer(endpoint);
        let bx = axis_q(&self.spans[cells[0]], pt[0])?;
        let by = axis_q(&self.spans[cells[1]], pt[1])?;
        for (i, [x, _]) in bx {
            for (&j, &[y, _]) in &by {
                add(
                    row,
                    self.ids[chart][i * self.dim + j].0,
                    sign.checked_mul(x.checked_mul(y)?)?
                        .checked_mul(self.ids[chart][i * self.dim + j].1)?,
                )?;
            }
        }
        Ok(())
    }
    fn jet(
        &self,
        side: Side,
        k: usize,
        t: Q,
        x_component: bool,
        sign: Q,
        row: &mut Row,
    ) -> Result<()> {
        let (chart, a, b) = side;
        let mut cells = [k; 2];
        cells[a] = b * (self.n - 1);
        let mut pt = [t; 2];
        pt[a] = integer(b);
        let parent = |j: usize| {
            integer(cells[j])
                .checked_add(pt[j])?
                .checked_div(integer(self.n))
        };
        let source = self.charts[chart].derivatives(parent(0)?, parent(1)?)?;
        let u = [source[0][0], source[0][1]];
        let v = [source[1][0], source[1][1]];
        let det = u[0]
            .checked_mul(v[1])?
            .checked_sub(v[0].checked_mul(u[1])?)?;
        let scale = sign.checked_mul(integer(self.n))?.checked_div(det)?;
        let minus = |x: Q| Q::zero().checked_sub(x);
        let (du, dv) = if x_component {
            (v[1], minus(u[1])?)
        } else {
            (minus(v[0])?, u[0])
        };
        let bx = axis_q(&self.spans[cells[0]], pt[0])?;
        let by = axis_q(&self.spans[cells[1]], pt[1])?;
        for (i, [x, dx]) in bx {
            for (&j, &[y, dy]) in &by {
                let value = du
                    .checked_mul(dx)?
                    .checked_mul(y)?
                    .checked_add(dv.checked_mul(x)?.checked_mul(dy)?)?;
                add(
                    row,
                    self.ids[chart][i * self.dim + j].0,
                    scale
                        .checked_mul(value)?
                        .checked_mul(self.ids[chart][i * self.dim + j].1)?,
                )?;
            }
        }
        Ok(())
    }
}

fn lift_residual(rows: &[Row], lift: &[BTreeMap<usize, f64>]) -> f64 {
    rows.iter()
        .map(|row| {
            let mut sum = BTreeMap::<usize, f64>::new();
            let mut magnitude = BTreeMap::<usize, f64>::new();
            for (&id, &coefficient) in row {
                for (&free, &value) in &lift[id] {
                    let term = number(coefficient) * value;
                    *sum.entry(free).or_default() += term;
                    *magnitude.entry(free).or_default() += term.abs();
                }
            }
            let scale = magnitude.values().copied().fold(0., f64::max);
            if scale == 0. {
                0.
            } else {
                sum.values().map(|x| x.abs() / scale).fold(0., f64::max)
            }
        })
        .fold(0., f64::max)
}
#[derive(PartialEq, Eq, PartialOrd, Ord)]
enum Node {
    Vertex(Point),
    Edge([Point; 2], AxisDof),
    Interior(usize, usize, usize),
}
fn axis_dofs(p: usize, n: usize) -> Vec<AxisDof> {
    let mut out = Vec::new();
    for knot in 0..=n {
        out.push(AxisDof {
            kind: AxisKind::Value,
            knot,
        });
        out.push(AxisDof {
            kind: AxisKind::Slope,
            knot,
        });
        if p == 4 && knot < n {
            out.push(AxisDof {
                kind: AxisKind::Middle,
                knot,
            });
        }
    }
    out
}
/// Exact local Bernstein controls from endpoint Hermite data. Quartic spans
/// retain an independent middle control; adjacent spans share only knot jets.
fn hermite_spans(p: usize, n: usize) -> Vec<Vec<Row>> {
    let dofs = axis_dofs(p, n);
    let find = |kind, knot| {
        dofs.iter()
            .position(|d| d.kind == kind && d.knot == knot)
            .unwrap()
    };
    (0..n)
        .map(|k| {
            let mut rows = vec![Row::new(); p + 1];
            rows[0].insert(find(AxisKind::Value, k), Q::one());
            rows[1].insert(find(AxisKind::Value, k), Q::one());
            rows[1].insert(find(AxisKind::Slope, k), Q::new(1, p as i128).unwrap());
            if p == 3 {
                rows[2].insert(find(AxisKind::Value, k + 1), Q::one());
                rows[2].insert(find(AxisKind::Slope, k + 1), Q::new(-1, p as i128).unwrap());
            } else {
                rows[2].insert(find(AxisKind::Middle, k), Q::one());
                rows[3].insert(find(AxisKind::Value, k + 1), Q::one());
                rows[3].insert(find(AxisKind::Slope, k + 1), Q::new(-1, p as i128).unwrap());
            }
            rows[p].insert(find(AxisKind::Value, k + 1), Q::one());
            rows
        })
        .collect()
}
fn edge_node(c: &Quad, axis: usize, endpoint: usize, tangential: AxisDof, n: usize) -> (Node, Q) {
    let [start, end] = endpoints(c, axis, endpoint);
    if start <= end {
        (Node::Edge([start, end], tangential), Q::one())
    } else {
        // Canonical edge orientation mirrors knot and span indices; only an
        // oriented tangential slope changes sign under that reversal.
        let mapped = AxisDof {
            kind: tangential.kind,
            knot: match tangential.kind {
                AxisKind::Value | AxisKind::Slope => n - tangential.knot,
                AxisKind::Middle => n - 1 - tangential.knot,
            },
        };
        let sign = if tangential.kind == AxisKind::Slope {
            Q::new(-1, 1).unwrap()
        } else {
            Q::one()
        };
        (Node::Edge([end, start], mapped), sign)
    }
}
fn node(c: &Quad, chart: usize, i: usize, j: usize, n: usize, dofs: &[AxisDof]) -> (Node, Q) {
    let (u, v) = (dofs[i], dofs[j]);
    let ub = (u.kind == AxisKind::Value && (u.knot == 0 || u.knot == n)).then_some(u.knot / n);
    let vb = (v.kind == AxisKind::Value && (v.knot == 0 || v.knot == n)).then_some(v.knot / n);
    match (ub, vb) {
        (Some(a), Some(b)) => (Node::Vertex(c[2 * a + b]), Q::one()),
        (Some(b), None) => edge_node(c, 0, b, v, n),
        (None, Some(b)) => edge_node(c, 1, b, u, n),
        (None, None) => (Node::Interior(chart, i, j), Q::one()),
    }
}
struct Leaf {
    c: Chart,
    z: [f64; 2],
    ij: [usize; 2],
    dofs: Vec<usize>,
    rows: Vec<Vec<(usize, f64)>>,
}
impl Leaf {
    fn new(
        b: &Basis,
        col: usize,
        z: [f64; 2],
        ij: [usize; 2],
        lift: &[BTreeMap<usize, f64>],
    ) -> Result<Self> {
        let mut rows = Vec::new();
        let mut active = BTreeSet::new();
        for a in &b.spans[ij[0]] {
            for d in &b.spans[ij[1]] {
                let mut row = BTreeMap::new();
                for (&i, &ca) in a {
                    for (&j, &cb) in d {
                        let w = ca.checked_mul(cb)?;
                        let (local_id, orientation) = b.ids[col][i * b.dim + j];
                        for (&id, &c) in &lift[local_id] {
                            *row.entry(id).or_insert(0.) += number(w.checked_mul(orientation)?) * c;
                        }
                    }
                }
                row.retain(|_, v| *v != 0.);
                active.extend(row.keys().copied());
                rows.push(row);
            }
        }
        let dofs = Vec::from_iter(active);
        let rows = rows
            .into_iter()
            .map(|r| {
                r.into_iter()
                    .map(|(id, v)| (dofs.binary_search(&id).unwrap(), v))
                    .collect()
            })
            .collect();
        Ok(Self {
            c: b.charts[col].clone(),
            z,
            ij,
            dofs,
            rows,
        })
    }
    fn sample(
        &self,
        n: usize,
        p: usize,
        t: [f64; 3],
    ) -> Result<(ElementGeometrySample, ElementBasisSample)> {
        let [u, v] = std::array::from_fn(|k| (self.ij[k] as f64 + t[k]) / n as f64);
        let height = self.z[1] - self.z[0];
        let z = self.z[0] + height * t[2];
        let mut pos = [0., 0., z];
        let mut jac = [[0.; 3]; 3];
        jac[2][2] = height;
        let (xy, dxy) = self.c.evaluate(u, v);
        pos[..2].copy_from_slice(&xy);
        jac[0][..2].copy_from_slice(&dxy[0].map(|x| x / n as f64));
        jac[1][..2].copy_from_slice(&dxy[1].map(|x| x / n as f64));
        let det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        let axis = |x: f64| {
            // Recurrence at both endpoints, including quartic end functions.
            let mut b = [0.; 5];
            b[0] = 1.;
            let mut d = [0.; 5];
            for k in 1..=p {
                if k == p {
                    for i in 0..k {
                        d[i] -= k as f64 * b[i];
                        d[i + 1] += k as f64 * b[i];
                    }
                }
                for i in (1..=k).rev() {
                    b[i] = (1. - x) * b[i] + x * b[i - 1];
                }
                b[0] *= 1. - x;
            }
            (b, d)
        };
        let (bx, dx) = axis(t[0]);
        let (by, dy) = axis(t[1]);
        let mut values = vec![0.; self.dofs.len()];
        let mut gradients = vec![[0.; 3]; self.dofs.len()];
        for i in 0..=p {
            for j in 0..=p {
                let val = bx[i] * by[j] / 20.;
                let du = dx[i] * by[j];
                let dv = bx[i] * dy[j];
                let gx = z / 20. * (jac[1][1] * du - jac[0][1] * dv) / det;
                let gy = z / 20. * (jac[0][0] * dv - jac[1][0] * du) / det;
                for &(id, c) in &self.rows[i * (p + 1) + j] {
                    values[id] += c * val * z;
                    gradients[id][0] += c * gx;
                    gradients[id][1] += c * gy;
                    gradients[id][2] += c * val;
                }
            }
        }
        Ok((
            ElementGeometrySample {
                position: pos,
                jacobian_columns: jac,
                determinant: det * height,
            },
            ElementBasisSample {
                dofs: self.dofs.clone(),
                values,
                physical_gradients: gradients,
            },
        ))
    }
}
/// Keep the exact kernel in its conditioned Hermite interpolation coordinates.
fn jet_lift(count: usize, equations: &[Row]) -> Result<(Vec<Row>, Vec<usize>)> {
    conditioned_kernel(count, equations)
}

/// Preserve the original bilinear benchmark API and its source geometry.
pub fn solve(p: usize, n: usize, c1: bool, order: usize) -> Result<Json> {
    if ![8, 10].contains(&order) {
        return Err("unsupported pocket study parameters");
    }
    let recipe = crate::json::parse(include_str!(
        "../tests/references/blind_bilinear_pocket.recipe.json"
    ))?;
    let model = crate::interchange::build(&recipe)?;
    let job = object([
        ("schema", Json::String("section_quintic_mms_v1".into())),
        ("quadrature", Json::Number(order.to_string())),
        ("x_origin_mm", Json::Number("-10".into())),
        ("x_scale_mm", Json::Number("70".into())),
    ]);
    solve_model(&model, p, n, c1, &job)
}

/// Solve on the actual biquadratic model published with the blind-pocket case.
pub fn solve_published(p: usize, n: usize, c1: bool, job: &Json) -> Result<Json> {
    let model = crate::json::parse(include_str!(
        "../docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"
    ))?;
    solve_model(&model, p, n, c1, job)
}

/// Run the conforming manufactured-field study on a supplied trivariate model.
pub fn solve_model(model: &Json, p: usize, n: usize, c1: bool, job: &Json) -> Result<Json> {
    solve_model_selected(
        model,
        p,
        n,
        if c1 {
            SectionJoins::All
        } else {
            SectionJoins::None
        },
        job,
    )
}

/// Solve with an explicit subset of physical C1 joins for diagnostic controls.
pub fn solve_model_selected(
    model: &Json,
    p: usize,
    n: usize,
    joins: SectionJoins,
    job: &Json,
) -> Result<Json> {
    crate::interchange::audit(model)?;
    solve_model_selected_audited(model, p, n, joins, job)
}

fn solve_model_selected_audited(
    model: &Json,
    p: usize,
    n: usize,
    joins: SectionJoins,
    job: &Json,
) -> Result<Json> {
    if ![3, 4].contains(&p) || ![1, 2, 4, 8].contains(&n) {
        return Err("unsupported pocket study parameters");
    }
    let blocks = decode_blocks(model)?;
    let mut z_bounds = [f64::INFINITY, f64::NEG_INFINITY];
    for b in &blocks {
        for h in &b.controls {
            let z = h[2] as f64 / h[3] as f64;
            z_bounds[0] = z_bounds[0].min(z);
            z_bounds[1] = z_bounds[1].max(z);
        }
    }
    if (z_bounds[0] - 0.).abs() > 1e-12 || (z_bounds[1] - 20.).abs() > 1e-12 {
        return Err("curved pocket study requires z bounds [0,20] mm");
    }
    let mut charts = Vec::new();
    let mut columns = Vec::new();
    for b in &blocks {
        let c = chart_from_block(b)?;
        let col = charts.iter().position(|q| q == &c).unwrap_or_else(|| {
            charts.push(c.clone());
            charts.len() - 1
        });
        columns.push(col);
    }
    if charts.len() != 5 || blocks.len() != 9 {
        return Err("unexpected scaffold topology");
    }
    let dofs = axis_dofs(p, n);
    let dim = dofs.len();
    let spans = hermite_spans(p, n);
    let mut nodes = BTreeMap::new();
    let mut ids = Vec::new();
    let mut edges: BTreeMap<[Point; 2], Vec<Side>> = BTreeMap::new();
    for (k, c) in charts.iter().enumerate() {
        let mut local = Vec::new();
        for i in 0..dim {
            for j in 0..dim {
                let next = nodes.len();
                let (key, sign) = node(&c.corners, k, i, j, n, &dofs);
                let id = *nodes.entry(key).or_insert(next);
                local.push((id, sign));
            }
        }
        ids.push(local);
        for a in 0..2 {
            for b in 0..2 {
                let mut pair = endpoints(&c.corners, a, b);
                pair.sort();
                edges.entry(pair).or_default().push((k, a, b));
            }
        }
    }
    for sides in edges.values() {
        if sides.len() > 2 {
            return Err("nonmanifold chart edge");
        }
        if let [a, b] = sides[..] {
            let first = charts[a.0].edge_controls(a.1, a.2);
            let mut second = charts[b.0].edge_controls(b.1, b.2);
            if endpoints(&charts[a.0].corners, a.1, a.2)
                != endpoints(&charts[b.0].corners, b.1, b.2)
            {
                second.reverse();
            }
            if first != second {
                return Err("shared pocket edge has mismatched full Bezier trace");
            }
        }
    }
    let count = nodes.len();
    let basis = Basis {
        n,
        dim,
        charts,
        ids,
        spans,
    };
    let mut equations = Vec::new();
    let mut value_probe_rows = Vec::new();
    let mut gradient_probe_rows = Vec::new();
    if joins != SectionJoins::None {
        for ([start, end], sides) in edges {
            if sides.len() == 1 {
                continue;
            }
            let [a, b] = sides[..] else {
                return Err("nonmanifold chart edge");
            };
            let radial = basis.charts[a.0].edge_is_straight(a.1, a.2)
                && basis.charts[b.0].edge_is_straight(b.1, b.2);
            let selected = match joins {
                SectionJoins::None => false,
                SectionJoins::All => true,
                SectionJoins::Radial => radial,
                SectionJoins::Footprint => !radial,
            };
            if !selected {
                continue;
            }
            let reverse = endpoints(&basis.charts[a.0].corners, a.1, a.2)
                != endpoints(&basis.charts[b.0].corners, b.1, b.2);
            // A physical coordinate direction transverse to the common tangent.
            let x_component = start[1] != end[1];
            for k in 0..n {
                // Clearing the physical-map Jacobian denominators gives a
                // univariate numerator of degree at most p+5 on each span.
                // p+6 exact points therefore imply equality throughout it.
                // Bilinear pairs have the sharper bound p+1. Keeping their
                // original sample grid preserves the conditioned Hermite lift.
                let jet_degree = if blocks.iter().all(|block| block.degrees == [1, 1, 1]) {
                    p + 1
                } else {
                    p + 5
                };
                for r in 0..=jet_degree {
                    let t = Q::new(r as i128, jet_degree as i128)?;
                    let mut row = Row::new();
                    basis.jet(a, k, t, x_component, Q::one(), &mut row)?;
                    let (kb, tb) = if reverse {
                        (n - 1 - k, Q::one().checked_sub(t)?)
                    } else {
                        (k, t)
                    };
                    basis.jet(b, kb, tb, x_component, Q::new(-1, 1)?, &mut row)?;
                    equations.push(row);
                }
                for r in 1..=2 {
                    let t = Q::new(r as i128, (p + 7) as i128)?;
                    let (kb, tb) = if reverse {
                        (n - 1 - k, Q::one().checked_sub(t)?)
                    } else {
                        (k, t)
                    };
                    let mut value = Row::new();
                    basis.value_trace(a, k, t, Q::one(), &mut value)?;
                    basis.value_trace(b, kb, tb, Q::new(-1, 1)?, &mut value)?;
                    value_probe_rows.push(value);
                    for x_component in [false, true] {
                        let mut gradient = Row::new();
                        basis.jet(a, k, t, x_component, Q::one(), &mut gradient)?;
                        basis.jet(b, kb, tb, x_component, Q::new(-1, 1)?, &mut gradient)?;
                        gradient_probe_rows.push(gradient);
                    }
                }
            }
        }
    }
    let (exact_lift, free) = jet_lift(count, &equations)?;
    let (lift, interpolation_residual) =
        if !equations.is_empty() && blocks.iter().any(|b| b.degrees == [2, 2, 1]) {
            crate::spline_constraints::evaluation_kernel(count, &equations, free.len())?
        } else {
            (
                exact_lift
                    .iter()
                    .map(|r| r.iter().map(|(&i, &q)| (i, number(q))).collect())
                    .collect(),
                0.,
            )
        };
    let value_probe_residual = lift_residual(&value_probe_rows, &lift);
    let gradient_probe_residual = lift_residual(&gradient_probe_rows, &lift);
    if !equations.is_empty() && (value_probe_residual > 1e-8 || gradient_probe_residual > 1e-8) {
        return Err("indeterminate: independent physical continuity probe residual");
    }
    let mut leaves = Vec::new();
    let mut boundary = Vec::new();
    let mut has_bottom = false;
    let mut has_top = false;
    let mut has_floor = false;
    let exterior = model.field("topology")?.field("exterior")?.items()?;
    for (k, b) in blocks.iter().enumerate() {
        let mut sides = Vec::new();
        for entry in exterior {
            let e = entry.items()?;
            if e[0].text()? == b.id {
                sides.push(decode_side_name(e[1].text()?)?)
            }
        }
        for &(axis, upper) in &sides {
            if axis == 2 {
                let h = &b.controls[usize::from(upper)];
                let z = h[2] as f64 / h[3] as f64;
                has_bottom |= z.abs() < 1e-12;
                has_top |= (z - 20.).abs() < 1e-12;
                has_floor |= (z - 14.).abs() < 1e-12;
            }
        }
        for i in 0..n {
            for j in 0..n {
                let ij = [i, j];
                for &(a, hi) in &sides {
                    if a == 2 || ij[a] == if hi { n - 1 } else { 0 } {
                        boundary.push((leaves.len(), a, hi));
                    }
                }
                leaves.push(Leaf::new(
                    &basis,
                    columns[k],
                    [
                        b.controls[0][2] as f64 / b.controls[0][3] as f64,
                        b.controls[1][2] as f64 / b.controls[1][3] as f64,
                    ],
                    ij,
                    &lift,
                )?);
            }
        }
    }
    if !(has_bottom && has_top && has_floor) {
        return Err("curved pocket model lacks bottom, upper, or floor boundary");
    }
    let essential = resolve(free.len(), &[], &[])?;
    let mut result = solve_sampled(
        &IntegrationDomain {
            scalar_dofs: free.len(),
            cell_count: leaves.len(),
            z_bounds,
            boundary: &boundary,
            essential: &essential,
        },
        job,
        |i, t| leaves[i].sample(n, p, t),
    )?;
    if let Json::Object(fields) = &mut result {
        let exact_degree_bound = p + 5;
        fields.insert(
            "chart_count".into(),
            Json::Number(basis.charts.len().to_string()),
        );
        fields.insert("span_count_per_axis".into(), Json::Number(n.to_string()));
        fields.insert(
            "interpolation_jet_residual".into(),
            Json::Number(interpolation_residual.to_string()),
        );
        fields.insert(
            "independent_c0_probe_relative_residual".into(),
            Json::Number(value_probe_residual.to_string()),
        );
        fields.insert(
            "independent_c1_gradient_probe_relative_residual".into(),
            Json::Number(gradient_probe_residual.to_string()),
        );
        fields.insert(
            "c1_join_selection".into(),
            Json::String(format!("{joins:?}")),
        );
        fields.insert(
            "c1_constraint_rows".into(),
            Json::Number(equations.len().to_string()),
        );
        fields.insert(
            "c1_constraint_rank".into(),
            Json::Number(count.saturating_sub(free.len()).to_string()),
        );
        fields.insert(
            "c1_gradient_numerator_degree_bound".into(),
            Json::Number(exact_degree_bound.to_string()),
        );
        fields.insert(
            "c1_points_per_span".into(),
            Json::Number((exact_degree_bound + 1).to_string()),
        );
    }
    Ok(result)
}
