//! Reconstruct four-edge NURBS contour loops from unordered sketch edges.

use crate::{
    json::{object, Json},
    nurbs_curve::extract_bezier,
    nurbs_multiloop,
};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

type Error = &'static str;
type Point = [crate::constraint_system::Ratio; 2];

#[derive(Clone)]
struct Curve {
    degree: usize,
    knots: Vec<i64>,
    controls: Vec<[i64; 3]>,
}

#[derive(Clone)]
struct Edge {
    id: String,
    curve: Curve,
    start: Point,
    end: Point,
}

#[derive(Clone)]
struct Operation {
    id: String,
    seed: String,
    z: [i64; 2],
}

fn parse_curve(value: &Json) -> Result<Curve, Error> {
    value.keys(&["degree", "knots", "controls"])?;
    let degree: usize = value.field("degree")?.number()?;
    let knots = value
        .field("knots")?
        .items()?
        .iter()
        .map(Json::number)
        .collect::<Result<Vec<i64>, _>>()?;
    let mut controls = Vec::new();
    for point in value.field("controls")?.items()? {
        let [x, y, w] = point.items()? else {
            return Err("edge control must be a homogeneous triple");
        };
        controls.push([x.number()?, y.number()?, w.number()?]);
    }
    extract_bezier(degree, &knots, &controls)?;
    fn gcd(mut a: u128, mut b: u128) -> u128 {
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a
    }
    let common = controls
        .iter()
        .flatten()
        .fold(0_u128, |g, &v| gcd(g, i128::from(v).unsigned_abs()));
    if common == 0 {
        return Err("zero homogeneous edge curve");
    }
    let common = i128::try_from(common).map_err(|_| "edge homogeneous gauge exceeds i128")?;
    let controls = controls
        .into_iter()
        .map(|p| {
            Ok([
                i64::try_from(i128::from(p[0]) / common)
                    .map_err(|_| "normalized edge control exceeds i64")?,
                i64::try_from(i128::from(p[1]) / common)
                    .map_err(|_| "normalized edge control exceeds i64")?,
                i64::try_from(i128::from(p[2]) / common)
                    .map_err(|_| "normalized edge weight exceeds i64")?,
            ])
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let n = controls
        .len()
        .checked_sub(1)
        .ok_or("empty edge control net")?;
    if degree >= controls.len()
        || knots.len() != controls.len() + degree + 1
        || knots[..=degree].iter().any(|&k| k != knots[degree])
        || knots[n + 1..].iter().any(|&k| k != knots[n + 1])
    {
        return Err("edge curves must be clamped NURBS");
    }
    Ok(Curve {
        degree,
        knots,
        controls,
    })
}

fn endpoint(curve: &Curve, last: bool) -> Result<Point, Error> {
    let p = if last {
        *curve.controls.last().ok_or("empty edge curve")?
    } else {
        curve.controls[0]
    };
    if p[2] <= 0 {
        return Err("edge endpoint weight must be positive");
    }
    Ok([
        crate::constraint_system::Ratio::new(i128::from(p[0]), i128::from(p[2]))?,
        crate::constraint_system::Ratio::new(i128::from(p[1]), i128::from(p[2]))?,
    ])
}

fn reverse(curve: &Curve) -> Result<Curve, Error> {
    let n = curve.controls.len() - 1;
    let a = curve.knots[curve.degree];
    let b = curve.knots[n + 1];
    let knots = curve
        .knots
        .iter()
        .rev()
        .map(|&k| {
            let reflected = i128::from(a) + i128::from(b) - i128::from(k);
            i64::try_from(reflected).map_err(|_| "reversed knot exceeds i64")
        })
        .collect::<Result<Vec<_>, _>>()?;
    let controls = curve.controls.iter().copied().rev().collect();
    Ok(Curve {
        degree: curve.degree,
        knots,
        controls,
    })
}

fn curve_json(curve: &Curve) -> Json {
    object([
        ("degree", Json::Number(curve.degree.to_string())),
        (
            "knots",
            Json::Array(
                curve
                    .knots
                    .iter()
                    .map(|v| Json::Number(v.to_string()))
                    .collect(),
            ),
        ),
        (
            "controls",
            Json::Array(
                curve
                    .controls
                    .iter()
                    .map(|p| Json::Array(p.iter().map(|v| Json::Number(v.to_string())).collect()))
                    .collect(),
            ),
        ),
    ])
}

fn parse_edges(recipe: &Json) -> Result<Vec<Edge>, Error> {
    let input = recipe.field("edges")?.items()?;
    if input.is_empty() || input.len() > 128 {
        return Err("edge count outside supported budget");
    }
    let mut ids = BTreeSet::new();
    let mut edges = Vec::new();
    for item in input {
        item.keys(&["id", "curve"])?;
        let id = item.field("id")?.text()?.to_owned();
        if id.is_empty()
            || id.len() > 64
            || !id
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || matches!(c, b'_' | b'-'))
            || !ids.insert(id.clone())
        {
            return Err("invalid or duplicate sketch edge identifier");
        }
        let curve = parse_curve(item.field("curve")?)?;
        let start = endpoint(&curve, false)?;
        let end = endpoint(&curve, true)?;
        if start == end {
            return Err("collapsed sketch edge");
        }
        edges.push(Edge {
            id,
            curve,
            start,
            end,
        });
    }
    Ok(edges)
}

fn parse_operations(recipe: &Json, edge_ids: &BTreeSet<String>) -> Result<Vec<Operation>, Error> {
    let input = recipe.field("cuts")?.items()?;
    if input.is_empty() || input.len() > 32 {
        return Err("cut count outside supported budget");
    }
    let mut ids = BTreeSet::new();
    let mut seeds = BTreeSet::new();
    let mut operations = Vec::new();
    for item in input {
        item.keys(&["id", "seed_edge", "z_interval_mm"])?;
        let id = item.field("id")?.text()?.to_owned();
        let seed = item.field("seed_edge")?.text()?.to_owned();
        let [a, b] = item.field("z_interval_mm")?.items()? else {
            return Err("cut interval must have two endpoints");
        };
        let z = [a.number()?, b.number()?];
        if id.is_empty()
            || id.len() > 40
            || !id
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || matches!(c, b'_' | b'-'))
            || !ids.insert(id)
        {
            return Err("invalid or duplicate cut identifier");
        }
        if !edge_ids.contains(&seed) {
            return Err("cut seed does not identify an edge");
        }
        if !seeds.insert(seed.clone()) {
            return Err("two cuts cannot use the same seed edge");
        }
        operations.push(Operation {
            id: item.field("id")?.text()?.to_owned(),
            seed,
            z,
        });
    }
    Ok(operations)
}

fn components(edges: &[Edge]) -> Result<Vec<Vec<usize>>, Error> {
    let mut adjacency = BTreeMap::<Point, Vec<usize>>::new();
    for (i, edge) in edges.iter().enumerate() {
        adjacency.entry(edge.start).or_default().push(i);
        adjacency.entry(edge.end).or_default().push(i);
    }
    if adjacency.values().any(|incident| incident.len() != 2) {
        return Err("sketch edge graph must have degree two at every endpoint");
    }
    let mut unseen: BTreeSet<usize> = (0..edges.len()).collect();
    let mut result = Vec::new();
    while let Some(&seed) = unseen.iter().next() {
        let mut queue = VecDeque::from([seed]);
        let mut component = BTreeSet::new();
        while let Some(edge_index) = queue.pop_front() {
            if !component.insert(edge_index) {
                continue;
            }
            unseen.remove(&edge_index);
            let edge = &edges[edge_index];
            for point in [edge.start, edge.end] {
                for &neighbor in adjacency.get(&point).ok_or("missing graph endpoint")? {
                    if !component.contains(&neighbor) {
                        queue.push_back(neighbor);
                    }
                }
            }
        }
        if component.len() != 4 {
            return Err("each contour component must have exactly four edges");
        }
        let vertices: BTreeSet<_> = component
            .iter()
            .flat_map(|&i| [edges[i].start, edges[i].end])
            .collect();
        if vertices.len() != 4 {
            return Err("each contour must be a simple four-edge cycle");
        }
        result.push(component.into_iter().collect());
    }
    Ok(result)
}

fn vertical_line(curve: &Curve) -> Result<bool, Error> {
    if curve.degree != 1 || curve.controls.len() != 2 {
        return Ok(false);
    }
    let a = endpoint(curve, false)?;
    let b = endpoint(curve, true)?;
    if a[0] != b[0] || a[1] == b[1] {
        return Ok(false);
    }
    for control in &curve.controls {
        let x =
            crate::constraint_system::Ratio::new(i128::from(control[0]), i128::from(control[2]))?;
        if x != a[0] {
            return Ok(false);
        }
    }
    Ok(true)
}

fn canonical_flanks(edges: &[Edge], component: &[usize]) -> Result<(Curve, Curve), Error> {
    let mut flanks = Vec::new();
    let mut verticals = Vec::new();
    for &i in component {
        let edge = &edges[i];
        if edge.start[0] == edge.end[0] {
            if !vertical_line(&edge.curve)? {
                return Err("closure edges must be straight and vertical");
            }
            verticals.push(edge);
        } else {
            let curve = if edge.start[0] < edge.end[0] {
                edge.curve.clone()
            } else {
                reverse(&edge.curve)?
            };
            let start = endpoint(&curve, false)?;
            let end = endpoint(&curve, true)?;
            if start[0] >= end[0] {
                return Err("profile flank must be strictly oriented in increasing x");
            }
            flanks.push((start, end, curve));
        }
    }
    if flanks.len() != 2 || verticals.len() != 2 {
        return Err("contour must contain two profile flanks and two vertical closures");
    }
    let (a0, a1, ca) = flanks.remove(0);
    let (b0, b1, cb) = flanks.remove(0);
    if a0[0] != b0[0] || a1[0] != b1[0] {
        return Err("profile flanks must share the same x projection");
    }
    let lower_upper = if a0[1] < b0[1] && a1[1] < b1[1] {
        (ca, cb)
    } else if b0[1] < a0[1] && b1[1] < a1[1] {
        (cb, ca)
    } else {
        return Err("profile flanks must remain consistently ordered at both ends");
    };
    let expected = [a0, a1, b0, b1];
    for endpoint in expected {
        if !verticals
            .iter()
            .any(|edge| edge.start == endpoint || edge.end == endpoint)
        {
            return Err("vertical closures must connect both profile endpoints");
        }
    }
    // The point checks above are necessary; this exact set check excludes
    // miswired pairs that happen to share one endpoint.
    for edge in verticals {
        let actual = BTreeSet::from([edge.start, edge.end]);
        let valid = [BTreeSet::from([a0, b0]), BTreeSet::from([a1, b1])];
        if !valid.contains(&actual) {
            return Err("vertical closures must join lower and upper endpoints");
        }
    }
    Ok(lower_upper)
}

fn x_start(curve: &Curve) -> Result<crate::constraint_system::Ratio, Error> {
    Ok(endpoint(curve, false)?[0])
}

/// Convert an unordered, seed-selected four-edge contour graph into the
/// canonical multi-loop recipe consumed by the M1 constructor.
pub fn normalize(recipe: &Json) -> Result<Json, Error> {
    recipe.keys(&["family", "bounds_xy_mm", "height_mm", "edges", "cuts"])?;
    if recipe.field("family")?.text()? != "nurbs_edge_loops_v1" {
        return Err("unsupported NURBS edge-loop family");
    }
    let edges = parse_edges(recipe)?;
    let edge_ids: BTreeSet<_> = edges.iter().map(|edge| edge.id.clone()).collect();
    let operations = parse_operations(recipe, &edge_ids)?;
    let components = components(&edges)?;
    if components.len() != operations.len() {
        return Err("every contour component requires exactly one cut operation");
    }
    let mut seed_component = BTreeMap::new();
    for (component_index, component) in components.iter().enumerate() {
        for &edge_index in component {
            seed_component.insert(edges[edge_index].id.as_str(), component_index);
        }
    }
    let mut assigned = BTreeMap::<usize, &Operation>::new();
    for operation in &operations {
        let component = *seed_component
            .get(operation.seed.as_str())
            .ok_or("cut seed is not in a contour")?;
        if assigned.insert(component, operation).is_some() {
            return Err("two cuts cannot target the same contour component");
        }
    }
    if assigned.len() != components.len() {
        return Err("every contour component requires exactly one cut operation");
    }
    let mut loops = Vec::new();
    for (component_index, component) in components.iter().enumerate() {
        let operation = assigned[&component_index];
        let (lower, upper) = canonical_flanks(&edges, component)?;
        loops.push((
            x_start(&lower)?,
            object([
                ("id", Json::String(operation.id.clone())),
                ("lower", curve_json(&lower)),
                ("upper", curve_json(&upper)),
                (
                    "z_interval_mm",
                    Json::Array(operation.z.map(|v| Json::Number(v.to_string())).into()),
                ),
            ]),
        ));
    }
    loops.sort_by(|a, b| a.0.cmp(&b.0));
    let normalized = object([
        ("family", Json::String("nurbs_multiloop_v1".into())),
        ("bounds_xy_mm", recipe.field("bounds_xy_mm")?.clone()),
        ("height_mm", recipe.field("height_mm")?.clone()),
        (
            "loops",
            Json::Array(loops.into_iter().map(|(_, contour)| contour).collect()),
        ),
    ]);
    // Reuse the accepted M1 contract for all stock, monotonicity, separation,
    // z-interval, and complete-partition checks.
    nurbs_multiloop::construct(&normalized)?;
    Ok(normalized)
}
