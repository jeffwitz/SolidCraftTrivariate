//! Uniform section h-refinement on the polygonal blind-pocket benchmark.
//! C0 identifies oriented traces. C1 additionally matches physical gradients
//! on every open span, including both endpoints for vertex enclosure.
use crate::{
    bracket_json::{decode_blocks, decode_side_name},
    constraint_system::{Ratio as Q, resolve},
    extracted_space::{ElementBasisSample, ElementGeometrySample},
    json::{Json, object},
    section_manufactured::{IntegrationDomain, solve_sampled},
    spline_constraints::{Result, Row, add, conditioned_kernel},
    transverse_jets::bernstein,
};
use std::collections::{BTreeMap, BTreeSet};

type Point = [i64; 2];
type Quad = [Point; 4];
type Side = (usize, usize, usize);
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
/// D(u,v) of the parent chart, with integer vertex differences.
fn tangent(c: &Quad, a: usize, t: Q) -> Result<[Q; 2]> {
    let lo = endpoints(c, 1 - a, 0);
    let hi = endpoints(c, 1 - a, 1);
    let mut v = [Q::zero(); 2];
    for j in 0..2 {
        let l = i128::from(lo[1][j]) - i128::from(lo[0][j]);
        let h = i128::from(hi[1][j]) - i128::from(hi[0][j]);
        v[j] = Q::new(l, 1)?.checked_add(t.checked_mul(Q::new(h - l, 1)?)?)?;
    }
    Ok(v)
}
struct Basis {
    n: usize,
    dim: usize,
    charts: Vec<Quad>,
    ids: Vec<Vec<(usize, Q)>>,
    spans: Vec<Vec<Row>>,
}
impl Basis {
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
        let u = tangent(&self.charts[chart], 0, parent(1)?)?;
        let v = tangent(&self.charts[chart], 1, parent(0)?)?;
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
    c: Quad,
    z: [f64; 2],
    ij: [usize; 2],
    dofs: Vec<usize>,
    rows: Vec<Vec<(usize, f64)>>,
}
impl Leaf {
    fn new(b: &Basis, col: usize, z: [f64; 2], ij: [usize; 2], lift: &[Row]) -> Result<Self> {
        let mut rows = Vec::new();
        let mut active = BTreeSet::new();
        for a in &b.spans[ij[0]] {
            for d in &b.spans[ij[1]] {
                let mut row = Row::new();
                for (&i, &ca) in a {
                    for (&j, &cb) in d {
                        let w = ca.checked_mul(cb)?;
                        let (local_id, orientation) = b.ids[col][i * b.dim + j];
                        for (&id, &c) in &lift[local_id] {
                            add(&mut row, id, w.checked_mul(c)?.checked_mul(orientation)?)?
                        }
                    }
                }
                active.extend(row.keys().copied());
                rows.push(row);
            }
        }
        let dofs = Vec::from_iter(active);
        let rows = rows
            .into_iter()
            .map(|r| {
                r.into_iter()
                    .map(|(id, v)| (dofs.binary_search(&id).unwrap(), number(v)))
                    .collect()
            })
            .collect();
        Ok(Self {
            c: b.charts[col],
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
        for k in 0..2 {
            let [a, b, c, d] = self.c.map(|xy| xy[k] as f64);
            pos[k] = (1. - u) * ((1. - v) * a + v * b) + u * ((1. - v) * c + v * d);
            jac[0][k] = ((1. - v) * (c - a) + v * (d - b)) / n as f64;
            jac[1][k] = ((1. - u) * (b - a) + u * (d - c)) / n as f64;
        }
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

/// Run the pre-registered conforming vector benchmark from its original sketch.
pub fn solve(p: usize, n: usize, c1: bool, order: usize) -> Result<Json> {
    if ![3, 4].contains(&p) || ![1, 2, 4, 8].contains(&n) || ![8, 10].contains(&order) {
        return Err("unsupported pocket study parameters");
    }
    let recipe = crate::json::parse(include_str!(
        "../tests/references/blind_bilinear_pocket.recipe.json"
    ))?;
    let model = crate::interchange::build(&recipe)?;
    let blocks = decode_blocks(&model)?;
    let mut charts = Vec::new();
    let mut columns = Vec::new();
    for b in &blocks {
        if b.degrees != [1; 3] || b.controls.iter().any(|h| h[3] != 1) {
            return Err("expected trilinear unit-weight scaffold");
        }
        let c = std::array::from_fn(|k| [b.controls[2 * k][0], b.controls[2 * k][1]]);
        let col = charts.iter().position(|q| q == &c).unwrap_or_else(|| {
            charts.push(c);
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
                let (key, sign) = node(c, k, i, j, n, &dofs);
                let id = *nodes.entry(key).or_insert(next);
                local.push((id, sign));
            }
        }
        ids.push(local);
        for a in 0..2 {
            for b in 0..2 {
                let mut pair = endpoints(c, a, b);
                pair.sort();
                edges.entry(pair).or_default().push((k, a, b));
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
    if c1 {
        for ([start, end], sides) in edges {
            if sides.len() == 1 {
                continue;
            }
            let [a, b] = sides[..] else {
                return Err("nonmanifold chart edge");
            };
            let reverse =
                endpoints(&basis.charts[a.0], a.1, a.2) != endpoints(&basis.charts[b.0], b.1, b.2);
            // A physical coordinate direction transverse to the common tangent.
            let x_component = start[1] != end[1];
            for k in 0..n {
                for r in 0..=p + 1 {
                    let t = Q::new(r as i128, (p + 1) as i128)?;
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
            }
        }
    }
    let (lift, free) = jet_lift(count, &equations)?;
    let mut leaves = Vec::new();
    let mut boundary = Vec::new();
    let exterior = model.field("topology")?.field("exterior")?.items()?;
    for (k, b) in blocks.iter().enumerate() {
        let mut sides = Vec::new();
        for entry in exterior {
            let e = entry.items()?;
            if e[0].text()? == b.id {
                sides.push(decode_side_name(e[1].text()?)?)
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
                    [b.controls[0][2] as f64, b.controls[1][2] as f64],
                    ij,
                    &lift,
                )?);
            }
        }
    }
    let essential = resolve(free.len(), &[], &[])?;
    let job = object([
        ("schema", Json::String("section_quintic_mms_v1".into())),
        ("quadrature", Json::Number(order.to_string())),
        ("x_origin_mm", Json::Number("-10".into())),
        ("x_scale_mm", Json::Number("70".into())),
    ]);
    solve_sampled(
        &IntegrationDomain {
            scalar_dofs: free.len(),
            cell_count: leaves.len(),
            z_bounds: [0., 20.],
            boundary: &boundary,
            essential: &essential,
        },
        &job,
        |i, t| leaves[i].sample(n, p, t),
    )
}
