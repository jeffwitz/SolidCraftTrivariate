//! Conforming 22-chart S2 elasticity controls with matching full-face traces.
//!
//! `C1` denotes knot removal only across the uniform interior v/w cell seams;
//! interfaces between the fitted charts remain coefficient-matched C0.

use crate::constrained_elasticity::solve_constrained;
use crate::constraint_system::{resolve, Ratio};
use crate::elasticity::gauss;
use crate::extracted_space::{assemble_sampled, AxisExtraction, ExtractedElement};
use crate::face_constraints::Row;
use crate::json::{object, Json};
use std::collections::{BTreeMap, BTreeSet};

type P = [f64; 4];
type Chart = (usize, usize, usize);
const N: usize = 36;
const YOUNG: f64 = 1000.;
const NU: f64 = 0.3;

#[derive(Clone)]
struct Cell {
    key: (Chart, usize, usize),
    dofs: Vec<usize>,
    element: ExtractedElement,
}
struct ChartNet {
    controls: Vec<P>,
    weights: Vec<Ratio>,
}
struct BuiltCells {
    cells: Vec<Cell>,
    labels: Vec<usize>,
    fixed: Vec<bool>,
    weights: Vec<Vec<Ratio>>,
}

fn axis(degree: usize) -> AxisExtraction {
    AxisExtraction {
        degree,
        basis_count: degree + 1,
        row_offsets: (0..=degree + 1).collect(),
        columns: (0..=degree).collect(),
        coefficients: (0..=degree).map(|_| Ratio::one()).collect(),
    }
}
fn ix(i: usize, j: usize, k: usize) -> usize {
    (i * 3 + j) * 3 + k
}
fn blend(a: P, b: P, t: f64) -> P {
    std::array::from_fn(|c| a[c] * (1. - t) + b[c] * t)
}
fn split(curve: &[P], t: f64) -> (Vec<P>, Vec<P>) {
    let p = curve.len() - 1;
    let mut levels = vec![curve.to_vec()];
    while levels.last().unwrap().len() > 1 {
        let last = levels.last().unwrap();
        levels.push(last.windows(2).map(|v| blend(v[0], v[1], t)).collect());
    }
    (
        (0..=p).map(|i| levels[i][0]).collect(),
        (0..=p).map(|i| levels[p - i][i]).collect(),
    )
}
fn restrict(curve: &[P], a: f64, b: f64) -> Vec<P> {
    if a == 0. && b == 1. {
        return curve.to_vec();
    }
    let tail = if a == 0. {
        curve.to_vec()
    } else {
        split(curve, a).1
    };
    if b == 1. {
        tail
    } else {
        split(&tail, (b - a) / (1. - a)).0
    }
}
fn eval_bilinear(net: &[P], i: usize, v: f64, w: f64) -> P {
    let a = blend(net[i * 4], net[i * 4 + 2], v);
    let b = blend(net[i * 4 + 1], net[i * 4 + 3], v);
    blend(a, b, w)
}
fn fitted(model: &Json) -> Result<BTreeMap<Chart, ChartNet>, &'static str> {
    let rise: i64 = model
        .field("source_recipe")?
        .field("transverse_rise_mm")?
        .number()?;
    if rise != 3 {
        return Err("the fitted bending control requires a 3 mm rise");
    }
    let source = crate::crossed_bending::source_net(model)?;
    let mut exact_weights = BTreeMap::new();
    for block in model.field("blocks")?.items()? {
        let id = block.field("id")?.text()?;
        let parts: Vec<_> = id.split('_').collect();
        let key = (
            parts[0][4..].parse().map_err(|_| "invalid S2 band")?,
            parts[1][3..].parse().map_err(|_| "invalid S2 column")?,
        );
        let mut weights = Vec::new();
        for control in block.field("controls")?.items()? {
            let w: i128 = control
                .items()?
                .get(3)
                .ok_or("missing homogeneous weight")?
                .number()?;
            weights.push(Ratio::new(w, 1)?);
        }
        exact_weights.insert(key, weights);
    }
    let mut out = BTreeMap::new();
    for (&(band, col), net) in &source {
        let layers: Vec<usize> = if col == 1 { vec![0, 1, 2] } else { vec![1] };
        for layer in layers {
            let mut fitted = net.clone();
            if col == 1 {
                for side in 0..2 {
                    let level = layer + side;
                    for i in 0..4 {
                        for j in 0..2 {
                            let p = &mut fitted[(i * 2 + j) * 2 + side];
                            p[2] = [0., 4., 8., 12.][level] * p[3];
                            if level == 1 || level == 2 {
                                p[2] += 0.03 * p[0];
                            }
                        }
                    }
                }
            }
            out.insert(
                (band, col, layer),
                ChartNet {
                    controls: fitted,
                    weights: exact_weights[&(band, col)].clone(),
                },
            );
        }
    }
    if out.len() != 22 {
        return Err("fitted S2 cover must contain 22 charts");
    }
    Ok(out)
}
fn average(a: Ratio, b: Ratio) -> Result<Ratio, &'static str> {
    a.checked_add(b)?.checked_div(Ratio::new(2, 1)?)
}
fn split_ratio(curve: &[Ratio], t: Ratio) -> Result<(Vec<Ratio>, Vec<Ratio>), &'static str> {
    let p = curve.len() - 1;
    let mut levels = vec![curve.to_vec()];
    while levels.last().unwrap().len() > 1 {
        let prev = levels.last().unwrap();
        let mut next = Vec::new();
        for pair in prev.windows(2) {
            next.push(
                pair[0]
                    .checked_mul(Ratio::one().checked_sub(t)?)?
                    .checked_add(pair[1].checked_mul(t)?)?,
            );
        }
        levels.push(next);
    }
    Ok((
        (0..=p).map(|i| levels[i][0]).collect(),
        (0..=p).map(|i| levels[p - i][i]).collect(),
    ))
}
fn restrict_ratio(curve: &[Ratio], a: Ratio, b: Ratio) -> Result<Vec<Ratio>, &'static str> {
    if a == Ratio::zero() && b == Ratio::one() {
        return Ok(curve.to_vec());
    }
    let tail = if a == Ratio::zero() {
        curve.to_vec()
    } else {
        split_ratio(curve, a)?.1
    };
    if b == Ratio::one() {
        Ok(tail)
    } else {
        Ok(split_ratio(
            &tail,
            b.checked_sub(a)?
                .checked_div(Ratio::one().checked_sub(a)?)?,
        )?
        .0)
    }
}
fn child(
    net: &ChartNet,
    iv: usize,
    iw: usize,
    nv: usize,
    nw: usize,
) -> Result<(Vec<P>, Vec<Ratio>), &'static str> {
    let mut elevated = vec![[0.; 4]; 4 * 3 * 3];
    for i in 0..4 {
        for j in 0..3 {
            for k in 0..3 {
                let v = j as f64 / 2.;
                let w = k as f64 / 2.;
                let a = eval_bilinear(&net.controls, i, v, w); // source has degree 1 in each transverse direction
                elevated[(i * 3 + j) * 3 + k] = a;
            }
        }
    }
    let mut sub = vec![[0.; 4]; 4 * 3 * 3];
    for i in 0..4 {
        for k in 0..3 {
            let c: Vec<P> = (0..3).map(|j| elevated[(i * 3 + j) * 3 + k]).collect();
            let r = restrict(&c, iv as f64 / nv as f64, (iv + 1) as f64 / nv as f64);
            for j in 0..3 {
                sub[(i * 3 + j) * 3 + k] = r[j];
            }
        }
        for j in 0..3 {
            let c: Vec<P> = (0..3).map(|k| sub[(i * 3 + j) * 3 + k]).collect();
            let r = restrict(&c, iw as f64 / nw as f64, (iw + 1) as f64 / nw as f64);
            for k in 0..3 {
                sub[(i * 3 + j) * 3 + k] = r[k];
            }
        }
    }
    let mut elevated_v = vec![Ratio::zero(); 4 * 3 * 2];
    for i in 0..4 {
        for j in 0..3 {
            for k in 0..2 {
                let a = net.weights[(i * 2) * 2 + k];
                let b = net.weights[(i * 2 + 1) * 2 + k];
                elevated_v[(i * 3 + j) * 2 + k] = if j == 0 {
                    a
                } else if j == 2 {
                    b
                } else {
                    average(a, b)?
                };
            }
        }
    }
    let mut elevated_w = vec![Ratio::zero(); 4 * 3 * 3];
    for i in 0..4 {
        for j in 0..3 {
            let a = elevated_v[(i * 3 + j) * 2];
            let b = elevated_v[(i * 3 + j) * 2 + 1];
            elevated_w[(i * 3 + j) * 3] = a;
            elevated_w[(i * 3 + j) * 3 + 1] = average(a, b)?;
            elevated_w[(i * 3 + j) * 3 + 2] = b;
        }
    }
    let mut sub_w = vec![Ratio::zero(); 4 * 3 * 3];
    for i in 0..4 {
        for k in 0..3 {
            let c: Vec<_> = (0..3).map(|j| elevated_w[(i * 3 + j) * 3 + k]).collect();
            let r = restrict_ratio(
                &c,
                Ratio::new(iv as i128, nv as i128)?,
                Ratio::new((iv + 1) as i128, nv as i128)?,
            )?;
            for j in 0..3 {
                sub_w[(i * 3 + j) * 3 + k] = r[j];
            }
        }
        for j in 0..3 {
            let c: Vec<_> = (0..3).map(|k| sub_w[(i * 3 + j) * 3 + k]).collect();
            let r = restrict_ratio(
                &c,
                Ratio::new(iw as i128, nw as i128)?,
                Ratio::new((iw + 1) as i128, nw as i128)?,
            )?;
            for k in 0..3 {
                sub_w[(i * 3 + j) * 3 + k] = r[k];
            }
        }
    }
    Ok((sub, sub_w))
}
fn compact_cells(
    charts: BTreeMap<Chart, ChartNet>,
    nv: usize,
    nw: usize,
) -> Result<BuiltCells, &'static str> {
    let mut cells = Vec::new();
    let mut lookup = BTreeMap::new();
    let mut exact_weights = Vec::with_capacity(charts.len() * nv * nw);
    for (&chart, net) in &charts {
        for iv in 0..nv {
            for iw in 0..nw {
                let (geom, weights) = child(net, iv, iw, nv, nw)?;
                let id = cells.len();
                exact_weights.push(weights);
                let dofs = (id * N..(id + 1) * N).collect::<Vec<_>>();
                let element = ExtractedElement {
                    axes: [axis(3), axis(2), axis(2)],
                    dofs: dofs.clone(),
                    field_weights: geom.iter().map(|p| p[3]).collect(),
                    geometry_controls: geom.clone(),
                };
                element.validate()?;
                lookup.insert((chart, iv, iw), id);
                cells.push(Cell {
                    key: (chart, iv, iw),
                    dofs,
                    element,
                });
            }
        }
    }
    let occ = cells.len() * N;
    let mut parent: Vec<usize> = (0..occ).collect();
    fn root(p: &mut [usize], mut a: usize) -> usize {
        while p[a] != a {
            p[a] = p[p[a]];
            a = p[a]
        }
        a
    }
    fn join(p: &mut [usize], a: usize, b: usize) {
        let ra = root(p, a);
        let rb = root(p, b);
        if ra != rb {
            p[ra] = rb
        }
    }
    for (e, c) in cells.iter().enumerate() {
        let (chart, iv, iw) = c.key;
        if iv + 1 < nv {
            let b = lookup[&(chart, iv + 1, iw)];
            for i in 0..4 {
                for k in 0..3 {
                    join(&mut parent, e * N + ix(i, 2, k), b * N + ix(i, 0, k));
                }
            }
        }
        if iw + 1 < nw {
            let b = lookup[&(chart, iv, iw + 1)];
            for i in 0..4 {
                for j in 0..3 {
                    join(&mut parent, e * N + ix(i, j, 2), b * N + ix(i, j, 0));
                }
            }
        }
        let (band, col, layer) = chart;
        for (other, axis) in [
            ((band + 1, col, layer), 0usize),
            ((band, col.checked_sub(1).unwrap_or(usize::MAX), layer), 1),
            ((band, col, layer + 1), 2),
        ] {
            let neighbor_key = match axis {
                1 if iv + 1 == nv => (other, 0, iw),
                2 if iw + 1 == nw => (other, iv, 0),
                1 | 2 => continue,
                _ => (other, iv, iw),
            };
            if let Some(&b) = lookup.get(&neighbor_key) {
                match axis {
                    0 => {
                        for j in 0..3 {
                            for k in 0..3 {
                                join(&mut parent, e * N + ix(3, j, k), b * N + ix(0, j, k));
                            }
                        }
                    }
                    1 => {
                        for i in 0..4 {
                            for k in 0..3 {
                                join(&mut parent, e * N + ix(i, 2, k), b * N + ix(i, 0, k));
                            }
                        }
                    }
                    _ => {
                        for i in 0..4 {
                            for j in 0..3 {
                                join(&mut parent, e * N + ix(i, j, 2), b * N + ix(i, j, 0));
                            }
                        }
                    }
                }
            }
        }
    }
    // Collapsed polar lines on the first and last cubic bands.
    for (band, end) in [(0, 0), (5, 3)] {
        for layer in 0..3 {
            let chart = (band, 1, layer);
            if !charts.contains_key(&chart) {
                continue;
            }
            for iw in 0..nw {
                for k in 0..3 {
                    let anchor = lookup[&(chart, 0, iw)] * N + ix(end, 0, k);
                    for iv in 0..nv {
                        let e = lookup[&(chart, iv, iw)];
                        for j in 0..3 {
                            join(&mut parent, anchor, e * N + ix(end, j, k));
                        }
                    }
                }
            }
        }
    }
    let roots: Vec<_> = (0..occ).map(|i| root(&mut parent, i)).collect();
    let uniq: BTreeSet<_> = roots.iter().copied().collect();
    let map: BTreeMap<_, _> = uniq.iter().enumerate().map(|(i, r)| (*r, i)).collect();
    let mut fixed = vec![false; occ];
    for (e, c) in cells.iter().enumerate() {
        let (_, col, iv) = (c.key.0 .0, c.key.0 .1, c.key.1);
        if col == 0 && iv + 1 == nv {
            for i in 0..4 {
                for k in 0..3 {
                    fixed[e * N + ix(i, 2, k)] = true;
                }
            }
        }
    }
    let labels = roots.iter().map(|r| map[r]).collect();
    Ok(BuiltCells {
        cells,
        labels,
        fixed,
        weights: exact_weights,
    })
}
fn solve(
    model: &Json,
    nv: usize,
    nw: usize,
    q: usize,
    continuity: &str,
) -> Result<Json, &'static str> {
    if !matches!(continuity, "C0" | "C1")
        || !(1..=6).contains(&nv)
        || !(1..=4).contains(&nw)
        || !(4..=16).contains(&q)
    {
        return Err("unsupported conforming bending parameters");
    }
    let built = compact_cells(fitted(model)?, nv, nw)?;
    let cells = built.cells;
    let labels = built.labels;
    let fixed = built.fixed;
    let exact_weights = built.weights;
    let nscalar = cells.len() * N;
    let mut representative = vec![usize::MAX; labels.iter().copied().max().unwrap_or(0) + 1];
    for (occ, &label) in labels.iter().enumerate() {
        if representative[label] == usize::MAX {
            representative[label] = occ;
        }
    }
    let mut fixed_label = vec![false; representative.len()];
    for (occ, &is_fixed) in fixed.iter().enumerate() {
        if is_fixed {
            fixed_label[labels[occ]] = true;
        }
    }
    let mut equations = Vec::new();
    for (occ, &label) in labels.iter().enumerate() {
        let root = representative[label];
        if occ != root {
            equations.push((
                occ,
                Row {
                    terms: vec![(root, 1)],
                    denominator: 1,
                },
            ));
        }
    }
    if continuity == "C1" {
        let lookup: BTreeMap<_, _> = cells.iter().enumerate().map(|(i, c)| (c.key, i)).collect();
        for (e, c) in cells.iter().enumerate() {
            let (chart, iv, iw) = c.key;
            for (axis, neighbor) in [
                (1usize, lookup.get(&(chart, iv + 1, iw)).copied()),
                (2, lookup.get(&(chart, iv, iw + 1)).copied()),
            ] {
                let Some(b) = neighbor else { continue };
                let pairs: Vec<(usize, usize, usize)> = if axis == 1 {
                    (0..4)
                        .flat_map(|i| (0..3).map(move |k| (ix(i, 1, k), ix(i, 2, k), ix(i, 1, k))))
                        .collect()
                } else {
                    (0..4)
                        .flat_map(|i| (0..3).map(move |j| (ix(i, j, 1), ix(i, j, 2), ix(i, j, 1))))
                        .collect()
                };
                for (left, center, right) in pairs {
                    let ga = cells[e].element.geometry_controls[left];
                    let gm = cells[e].element.geometry_controls[center];
                    let gb = cells[b].element.geometry_controls[right];
                    if (0..4).any(|c| {
                        ((ga[c] + gb[c]) * 0.5 - gm[c]).abs()
                            > 3e-13 * ga[c].abs().max(gb[c].abs()).max(gm[c].abs()) + 3e-12
                    }) {
                        return Err(
                            "C1 knot removal requires matching homogeneous geometry coefficients",
                        );
                    }
                    let a = exact_weights[e][left];
                    let m = exact_weights[e][center];
                    let bweight = exact_weights[b][right];
                    if a.checked_add(bweight)? != m.checked_mul(Ratio::new(2, 1)?)? {
                        return Err("C1 knot-removal denominator condition failed");
                    }
                    let ids = [
                        labels[c.dofs[left]],
                        labels[c.dofs[center]],
                        labels[cells[b].dofs[right]],
                    ];
                    let mapped = ids.map(|id| {
                        if fixed_label[id] {
                            usize::MAX
                        } else {
                            representative[id]
                        }
                    });
                    let [l, midx, r] = mapped;
                    if l == midx && midx == r {
                        continue;
                    }
                    if mapped.contains(&usize::MAX) {
                        return Err("C1 knot removal crosses the clamp");
                    }
                    if midx == l || midx == r {
                        return Err("unexpected degenerate C1 constraint");
                    }
                    let al = exact_weights[e][left]
                        .checked_div(exact_weights[e][center].checked_mul(Ratio::new(2, 1)?)?)?;
                    let ar = exact_weights[b][right]
                        .checked_div(exact_weights[e][center].checked_mul(Ratio::new(2, 1)?)?)?;
                    let common = al
                        .denominator()
                        .checked_mul(ar.denominator())
                        .ok_or("C1 denominator overflow")?;
                    let wa = al
                        .numerator()
                        .checked_mul(common / al.denominator())
                        .ok_or("C1 coefficient overflow")?;
                    let wb = ar
                        .numerator()
                        .checked_mul(common / ar.denominator())
                        .ok_or("C1 coefficient overflow")?;
                    let denominator = wa.checked_add(wb).ok_or("C1 denominator overflow")?;
                    equations.push((
                        midx,
                        Row {
                            terms: vec![(l, wa), (r, wb)],
                            denominator,
                        },
                    ));
                }
            }
        }
    }
    let prescribed: Vec<_> = fixed
        .iter()
        .enumerate()
        .filter_map(|(i, &is_fixed)| is_fixed.then_some((i, Ratio::zero())))
        .collect();
    let lift = resolve(nscalar, &equations, &prescribed)?;
    let cell_count = cells.len();
    let lambda = YOUNG * NU / ((1. + NU) * (1. - 2. * NU));
    let mu = YOUNG / (2. * (1. + NU));
    let assembled = assemble_sampled(
        nscalar,
        cell_count,
        lambda,
        mu,
        q,
        |i, p| cells[i].element.sample(p),
        |_| [0.; 3],
    )?;
    let mut force = vec![0.; 3 * nscalar];
    let gp = gauss(q)?;
    for c in &cells {
        let (chart, iv, _) = c.key;
        if chart.1 != 2 || iv != 0 {
            continue;
        }
        for &(u, wu) in &gp {
            for &(w, ww) in &gp {
                let (geom, basis) = c.element.sample([u, 0., w])?;
                let a = geom.jacobian_columns[0];
                let b = geom.jacobian_columns[2];
                let cross = [
                    a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0],
                ];
                let measure = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2])
                    .sqrt()
                    * wu
                    * ww;
                for (j, &dof) in basis.dofs.iter().enumerate() {
                    force[3 * dof] += basis.values[j] * measure;
                }
            }
        }
    }
    let solution = solve_constrained(
        &assembled.stiffness,
        [&lift, &lift, &lift],
        &force,
        1e-10,
        20_000,
    )?;
    let ftotal: f64 = force.iter().step_by(3).sum();
    let work: f64 = force
        .iter()
        .zip(&solution.displacement_mm)
        .map(|(a, b)| a * b)
        .sum();
    let dofs = 3 * lift.free_roots.len();
    Ok(object([
        (
            "schema",
            Json::String("conforming_crossed_bending_result_v1".into()),
        ),
        ("rise_mm", Json::Number("3".into())),
        (
            "degrees_uvw",
            Json::Array(
                [3, 2, 2]
                    .into_iter()
                    .map(|d| Json::Number(d.to_string()))
                    .collect(),
            ),
        ),
        (
            "subdivisions_vw",
            Json::Array(
                [nv, nw]
                    .into_iter()
                    .map(|d| Json::Number(d.to_string()))
                    .collect(),
            ),
        ),
        ("quadrature_order", Json::Number(q.to_string())),
        ("geometry_charts", Json::Number("22".into())),
        (
            "continuity_scope",
            Json::String(
                "C0 chart interfaces; selected continuity applies only to interior v/w cell seams"
                    .into(),
            ),
        ),
        ("continuity", Json::String(continuity.into())),
        ("free_vector_dofs", Json::Number(dofs.to_string())),
        ("mean_ux_mm", Json::Number((work / ftotal).to_string())),
        (
            "energy_N_mm",
            Json::Number(solution.strain_energy_n_mm.to_string()),
        ),
        ("force_N", Json::Number(ftotal.to_string())),
        (
            "volume_mm3",
            Json::Number(assembled.stiffness.integrated_volume.to_string()),
        ),
        (
            "relative_residual",
            Json::Number(solution.relative_residual.to_string()),
        ),
        ("cells", Json::Number(cell_count.to_string())),
    ]))
}

/// Solve the recorded fitted S2 conforming bending control.
pub fn analyze(
    model: &Json,
    nv: usize,
    nw: usize,
    order: usize,
    continuity: &str,
) -> Result<Json, &'static str> {
    solve(model, nv, nw, order, continuity)
}
