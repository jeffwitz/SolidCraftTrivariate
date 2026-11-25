//! Conforming material charts for crossed extrusions with NURBS side walls.
//! All refinement, synchronization and acceptance predicates use checked rationals.
use crate::bracket::BracketBlock;
use crate::constraint_system::Ratio as R;
use crate::nurbs_curve::{extract_bezier, BezierSpan};
use std::collections::BTreeSet;

/// Planar homogeneous controls are (w * lateral coordinate, w * y, w).
#[derive(Clone, Debug)]
pub struct NurbsCurve2 {
    pub degree: usize,
    pub knots: Vec<i64>,
    pub controls: Vec<[i64; 3]>,
}

type Error = &'static str;
type Point = [R; 3];
fn int(n: i64) -> R {
    R::new(i128::from(n), 1).expect("integer")
}
fn div(a: R, b: R) -> Result<R, Error> {
    a.checked_mul(R::new(b.denominator(), b.numerator())?)
}
fn less(a: R, b: R) -> Result<bool, Error> {
    Ok(a.checked_sub(b)?.numerator() < 0)
}
fn ratio(a: i64, b: i64, x: i64) -> Result<R, Error> {
    R::new(i128::from(x) - i128::from(a), i128::from(b) - i128::from(a))
}
fn blend(a: Point, b: Point, t: R) -> Result<Point, Error> {
    let mut p = a;
    for c in 0..3 {
        p[c] = a[c].checked_add(t.checked_mul(b[c].checked_sub(a[c])?)?)?;
    }
    Ok(p)
}

fn split(net: &[Point], t: R) -> Result<(Vec<Point>, Vec<Point>), Error> {
    let mut row = net.to_vec();
    let mut left = Vec::new();
    let mut right = Vec::new();
    while !row.is_empty() {
        left.push(row[0]);
        right.push(row[row.len() - 1]);
        row = row
            .windows(2)
            .map(|pair| blend(pair[0], pair[1], t))
            .collect::<Result<_, _>>()?;
    }
    right.reverse();
    Ok((left, right))
}

fn restricted(span: &BezierSpan, interval: [i64; 2], degree: usize) -> Result<Vec<Point>, Error> {
    let [a, b] = span.interval;
    let [lo, hi] = interval;
    let mut net: Vec<_> = span.controls.iter().map(|p| p.map(int)).collect();
    if lo > a {
        net = split(&net, ratio(a, b, lo)?)?.1;
    }
    if hi < b {
        net = split(&net, ratio(lo, b, hi)?)?.0;
    }
    while net.len() <= degree {
        let n = net.len();
        let mut elevated = vec![net[0]];
        for i in 1..n {
            elevated.push(blend(net[i], net[i - 1], R::new(i as i128, n as i128)?)?);
        }
        elevated.push(net[n - 1]);
        net = elevated;
    }
    let scale = net[0][2];
    for p in &mut net {
        for value in p {
            *value = div(*value, scale)?;
        }
    }
    Ok(net)
}

fn endpoint(spans: &[BezierSpan], upper: bool) -> Result<[R; 2], Error> {
    let p = if upper {
        spans.last().and_then(|s| s.controls.last())
    } else {
        spans.first().and_then(|s| s.controls.first())
    }
    .ok_or("empty NURBS")?;
    Ok([
        R::new(p[0].into(), p[2].into())?,
        R::new(p[1].into(), p[2].into())?,
    ])
}

// Store one common integer scale per block; cancel common factors on conversion.
fn integer_net(net: &[[R; 4]]) -> Result<Vec<[i64; 4]>, Error> {
    fn gcd(mut a: i128, mut b: i128) -> i128 {
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a
    }
    let mut scale = 1i128;
    for c in net.iter().flatten() {
        scale = (scale / gcd(scale, c.denominator()))
            .checked_mul(c.denominator())
            .ok_or("NURBS integral scale exceeds i128")?;
    }
    let mut integers = Vec::with_capacity(net.len());
    let mut common = 0;
    for p in net {
        let mut q = [0i128; 4];
        for c in 0..4 {
            q[c] = p[c]
                .numerator()
                .checked_mul(scale / p[c].denominator())
                .ok_or("NURBS integer coefficient exceeds i128")?;
            common = gcd(
                common,
                q[c].checked_abs()
                    .ok_or("NURBS integer magnitude exceeds i128")?,
            );
        }
        integers.push(q);
    }
    integers
        .into_iter()
        .map(|p| {
            let mut q = [0; 4];
            for c in 0..4 {
                q[c] = i64::try_from(p[c] / common)
                    .map_err(|_| "NURBS controls exceed native i64 storage")?;
            }
            Ok(q)
        })
        .collect()
}

// One row at a Bernstein index contains lateral bounds X0..X3, Z0..Z3,
// the common ordinate and the common positive homogeneous weight.
struct Row {
    x: [R; 4],
    z: [R; 4],
    y: R,
    weight: R,
}

fn append_band(
    blocks: &mut Vec<BracketBlock>,
    band: usize,
    rows: &[Row],
    active: [bool; 2],
) -> Result<(), Error> {
    for (i, row) in rows.iter().enumerate() {
        if row.weight.numerator() <= 0 {
            return Err("nonpositive NURBS denominator");
        }
        for bounds in [&row.x, &row.z] {
            for pair in bounds.windows(2) {
                if !less(pair[0], pair[1])? {
                    return Err("unproved NURBS wall separation");
                }
            }
        }
        if i > 0
            && !less(
                div(rows[i - 1].y, rows[i - 1].weight)?,
                div(row.y, row.weight)?,
            )?
        {
            return Err("unproved strict ordinate monotonicity");
        }
    }
    for x in 0..3 {
        for z in 0..3 {
            if (active[0] && x == 1) || (active[1] && z == 1) {
                continue;
            }
            let mut net = Vec::with_capacity(4 * rows.len());
            for i in x..=x + 1 {
                for row in rows {
                    for k in z..=z + 1 {
                        net.push([row.x[i], row.y, row.z[k], row.weight]);
                    }
                }
            }
            blocks.push(BracketBlock {
                id: format!("span_{band}_x{x}_z{z}"),
                degrees: [1, rows.len() - 1, 1],
                controls: integer_net(&net)?,
            });
        }
    }
    Ok(())
}

fn make_row(bounds: [[i64; 2]; 3], h: [R; 4], y: R, w: R) -> Result<Row, Error> {
    Ok(Row {
        x: [
            int(bounds[0][0]).checked_mul(w)?,
            h[0],
            h[1],
            int(bounds[0][1]).checked_mul(w)?,
        ],
        z: [
            int(bounds[2][0]).checked_mul(w)?,
            h[2],
            h[3],
            int(bounds[2][1]).checked_mul(w)?,
        ],
        y,
        weight: w,
    })
}

/// Build from four clamped curves (x_min,x_max,z_min,z_max). In the active
/// spans the control ordinates and weights must coincide after exact refinement
/// and normalization. The lateral controls are independent.
pub fn build_cross(
    bounds: [[i64; 2]; 3],
    profiles: [&NurbsCurve2; 4],
) -> Result<Vec<BracketBlock>, Error> {
    if bounds.iter().any(|b| b[0] >= b[1]) {
        return Err("invalid NURBS stock bounds");
    }
    let mut curves = Vec::new();
    let mut breaks = BTreeSet::new();
    let mut ranges = Vec::new();
    for curve in profiles {
        let spans = extract_bezier(curve.degree, &curve.knots, &curve.controls)?;
        let first = spans.first().ok_or("empty NURBS support")?;
        let last = spans.last().unwrap();
        ranges.push([first.interval[0], last.interval[1]]);
        for span in &spans {
            breaks.extend(span.interval);
        }
        curves.push(spans);
    }
    if ranges[0] != ranges[1] || ranges[2] != ranges[3] {
        return Err("opposite sketch walls need common parameter endpoints");
    }
    let mut bands = Vec::new();
    let breaks: Vec<_> = breaks.into_iter().collect();
    for window in breaks.windows(2) {
        let [lo, hi] = [window[0], window[1]];
        let spans: Vec<_> = curves
            .iter()
            .map(|ss| {
                ss.iter()
                    .find(|s| s.interval[0] <= lo && hi <= s.interval[1])
            })
            .collect();
        let degree = spans
            .iter()
            .flatten()
            .map(|s| s.controls.len() - 1)
            .max()
            .unwrap_or(1);
        let nets: Vec<_> = spans
            .iter()
            .map(|s| s.map(|p| restricted(p, [lo, hi], degree)).transpose())
            .collect::<Result<_, _>>()?;
        let active = [nets[0].is_some(), nets[2].is_some()];
        let mut rows = Vec::new();
        if let Some(reference) = nets.iter().flatten().next() {
            for j in 0..=degree {
                let [_, y, w] = reference[j];
                let mut h = [R::zero(); 4];
                for k in 0..4 {
                    if let Some(net) = &nets[k] {
                        if net[j][1] != y || net[j][2] != w {
                            return Err(
                                "NURBS sketch ordinates require a different reparameterization",
                            );
                        }
                        h[k] = net[j][0];
                    } else {
                        h[k] = endpoint(&curves[k], lo >= ranges[k][1])?[0].checked_mul(w)?;
                    }
                }
                rows.push(make_row(bounds, h, y, w)?);
            }
        } else {
            let mut h = [R::zero(); 4];
            for k in 0..4 {
                h[k] = endpoint(&curves[k], lo >= ranges[k][1])?[0];
            }
            let left = ranges
                .iter()
                .position(|s| s[1] == lo)
                .ok_or("missing lower ordinate")?;
            let right = ranges
                .iter()
                .position(|s| s[0] == hi)
                .ok_or("missing upper ordinate")?;
            for y in [
                endpoint(&curves[left], true)?[1],
                endpoint(&curves[right], false)?[1],
            ] {
                rows.push(make_row(bounds, h, y, R::one())?);
            }
        }
        bands.push((rows, active));
    }
    // All event sections must meet exactly, including a start and an end
    // at the same parameter value belonging to otherwise disjoint sketches.
    fn section(row: &Row) -> Result<[R; 5], Error> {
        let mut plane = [R::zero(); 5];
        for (i, c) in [row.x[1], row.x[2], row.z[1], row.z[2], row.y]
            .into_iter()
            .enumerate()
        {
            plane[i] = div(c, row.weight)?;
        }
        Ok(plane)
    }
    for adjacent in bands.windows(2) {
        if section(adjacent[0].0.last().unwrap())? != section(&adjacent[1].0[0])? {
            return Err("incompatible NURBS event traces");
        }
    }
    let mut blocks = Vec::new();
    let mut band = 0;
    let cap = |row: &Row,
               stock: i64,
               upper: bool,
               band: usize,
               blocks: &mut Vec<BracketBlock>|
     -> Result<(), Error> {
        let [a, b, c, d, y] = section(row)?;
        let levels = if upper {
            [y, int(stock)]
        } else {
            [int(stock), y]
        };
        if !less(levels[0], levels[1])? {
            return Err("NURBS profile end must be inside stock");
        }
        let rows = [
            make_row(bounds, [a, b, c, d], levels[0], R::one())?,
            make_row(bounds, [a, b, c, d], levels[1], R::one())?,
        ];
        append_band(blocks, band, &rows, [false, false])
    };
    cap(&bands[0].0[0], bounds[1][0], false, band, &mut blocks)?;
    for (rows, active) in &bands {
        band += 1;
        append_band(&mut blocks, band, rows, *active)?;
    }
    cap(
        bands.last().unwrap().0.last().unwrap(),
        bounds[1][1],
        true,
        band + 1,
        &mut blocks,
    )?;
    Ok(blocks)
}

// Read the piecewise linear graph at a physical ordinate. A cut that is
// inactive here still leaves its end trace on the neighboring material.
fn graph_at(graph: &[[R; 2]], cursor: &mut usize, y: R) -> Result<R, Error> {
    let end = graph.len() - 1;
    if !less(graph[0][1], y)? {
        return Ok(graph[0][0]);
    }
    if !less(y, graph[end][1])? {
        return Ok(graph[end][0]);
    }
    while less(graph[*cursor + 1][1], y)? {
        *cursor += 1;
    }
    let [x0, y0] = graph[*cursor];
    let [x1, y1] = graph[*cursor + 1];
    x0.checked_add(
        x1.checked_sub(x0)?
            .checked_mul(div(y.checked_sub(y0)?, y1.checked_sub(y0)?)?)?,
    )
}

/// Certified approximation for unequal ordinate parameterizations. The wall
/// at each physical y differs laterally by at most tolerance. A simultaneous
/// monotone piecewise linear map in x and z bounds the 3D Hausdorff distance
/// in both directions by sqrt(2)*tolerance, hence by 2*tolerance.
pub fn build_cross_bounded(
    bounds: [[i64; 2]; 3],
    profiles: [&NurbsCurve2; 4],
    tolerance: R,
    max_depth: usize,
) -> Result<Vec<BracketBlock>, Error> {
    use crate::nurbs_graph_approx::certified_graph;
    if bounds.iter().any(|s| s[0] >= s[1]) || tolerance.numerator() <= 0 || max_depth > 12 {
        return Err("invalid bounded NURBS parameters");
    }
    let half = tolerance.checked_mul(R::new(1, 2)?)?;
    let two_tol = tolerance.checked_mul(int(2))?;
    let graphs = profiles
        .iter()
        .map(|c| certified_graph(c, half, max_depth))
        .collect::<Result<Vec<_>, _>>()?;
    let endpoints: Vec<_> = graphs
        .iter()
        .map(|g| [g[0][1], g[g.len() - 1][1]])
        .collect();
    if endpoints[0] != endpoints[1] || endpoints[2] != endpoints[3] {
        return Err("opposite sketch walls must meet their transverse caps");
    }
    let bottom = int(bounds[1][0]);
    let top = int(bounds[1][1]);
    if endpoints
        .iter()
        .any(|e| !matches!(less(bottom, e[0]), Ok(true)) || !matches!(less(e[1], top), Ok(true)))
    {
        return Err("NURBS caps require clearance inside stock");
    }
    // Merge four already ordered event lists; Ratio's derived Ord orders
    // numerator/denominator tuples, so use exact arithmetic comparisons.
    let mut cursors = [0; 4];
    let mut events = vec![bottom];
    loop {
        let mut next = None;
        for (k, g) in graphs.iter().enumerate() {
            if let Some(p) = g.get(cursors[k]) {
                if match next {
                    None => true,
                    Some(y) => less(p[1], y)?,
                } {
                    next = Some(p[1]);
                }
            }
        }
        let Some(y) = next else { break };
        if events.len() > 8192 {
            return Err("bounded NURBS event budget exceeded");
        }
        events.push(y);
        for (k, g) in graphs.iter().enumerate() {
            if g.get(cursors[k]).is_some_and(|p| p[1] == y) {
                cursors[k] += 1;
            }
        }
    }
    events.push(top);
    // Round all lateral values onto the same rational lattice to keep native
    // integer gauges small. The extra error is strictly below half the budget.
    let snap = |c: R| -> Result<R, Error> {
        let cells = div(c, half)?;
        R::new(cells.numerator().div_euclid(cells.denominator()), 1)?.checked_mul(half)
    };
    cursors = [0; 4];
    let mut rows = Vec::with_capacity(events.len());
    for &y in &events {
        let mut h = [R::zero(); 4];
        for k in 0..4 {
            h[k] = snap(graph_at(&graphs[k], &mut cursors[k], y)?)?;
        }
        let row = make_row(bounds, h, y, R::one())?;
        // The error bands of opposite walls and the stock must remain disjoint.
        for axis in [&row.x, &row.z] {
            for (i, pair) in axis.windows(2).enumerate() {
                let budget = if i == 1 { two_tol } else { tolerance };
                if !less(budget, pair[1].checked_sub(pair[0])?)? {
                    return Err("NURBS approximation cannot certify ligament clearance");
                }
            }
        }
        rows.push(row);
    }
    let mut blocks = Vec::new();
    for (band, pair) in rows.windows(2).enumerate() {
        let covers = |i: usize| -> Result<bool, Error> {
            Ok(!less(pair[0].y, endpoints[i][0])? && !less(endpoints[i][1], pair[1].y)?)
        };
        append_band(&mut blocks, band, pair, [covers(0)?, covers(2)?])?;
    }
    Ok(blocks)
}

#[path = "nurbs_oval.rs"]
mod oval;
pub use oval::build_oval;
