//! Complete a curved tool partition to the six planes of a rectangular stock.
use crate::bracket::BracketBlock;
use crate::face_profile_extrusion::{FaceProfileExtrusion, SectionPatch, lift_sections};

fn xy(p: [i64; 3]) -> [f64; 2] {
    [p[0] as f64 / p[2] as f64, p[1] as f64 / p[2] as f64]
}

fn direction(p: [f64; 2], center: [f64; 2]) -> [f64; 2] {
    let d = [p[0] - center[0], p[1] - center[1]];
    let length = d[0].hypot(d[1]);
    [d[0] / length, d[1] / length]
}

// Match the four corners to distinct curve vertices in cyclic order. Minimize
// angular displacement about the mean boundary vertex, then divide each side
// according to chord length. Final acceptance tests the maps themselves.
fn perimeter(ring: &[SectionPatch], bounds: [[f64; 2]; 2]) -> Result<Vec<[f64; 2]>, &'static str> {
    let n = ring.len();
    if !(4..=1024).contains(&n) {
        return Err("rectangular enclosure needs 4 to 1024 boundary spans");
    }
    let pts: Vec<_> = ring.iter().map(|p| xy(p.controls[0])).collect();
    for p in ring.iter().flat_map(|p| p.controls.iter().step_by(2)) {
        let v = xy(*p);
        if (0..2).any(|j| v[j] <= bounds[j][0] || v[j] >= bounds[j][1]) {
            return Err("profile control hull must be strictly inside the support rectangle");
        }
    }
    let [[x0, x1], [y0, y1]] = bounds;
    let corners = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
    let center = std::array::from_fn(|k| pts.iter().map(|p| p[k]).sum::<f64>() / n as f64);
    let dirs = pts
        .iter()
        .map(|&p| direction(p, center))
        .collect::<Vec<_>>();
    let cost = |i: usize, j: usize| {
        let a = direction(corners[j], center);
        let b = dirs[i % n];
        1. - (a[0] * b[0] + a[1] * b[1])
    };
    let mut best = f64::INFINITY;
    let mut winner = [0; 4];
    for start in 0..n {
        let mut states = [(f64::INFINITY, [0; 4]); 4];
        states[0] = (cost(start, 0), [start; 4]);
        for k in start + 1..start + n {
            for corner in (1..=3).rev() {
                let (prev, mut ids) = states[corner - 1];
                let score = prev + cost(k, corner);
                if score < states[corner].0 {
                    ids[corner] = k;
                    states[corner] = (score, ids);
                }
            }
        }
        if states[3].0 < best {
            (best, winner) = states[3];
        }
    }
    if !best.is_finite() {
        return Err("cannot find an ordered rectangular enclosure");
    }
    let mut out = vec![[0.; 2]; n];
    let length = |i: usize| {
        let a = pts[i % n];
        let b = pts[(i + 1) % n];
        (a[0] - b[0]).hypot(a[1] - b[1])
    };
    for j in 0..4 {
        let start = winner[j];
        let end = if j == 3 { winner[0] + n } else { winner[j + 1] };
        let total = (start..end).map(length).sum::<f64>();
        if total <= 0. {
            return Err("zero-length interval on support rectangle");
        }
        let mut dist = 0.;
        for i in start..end {
            let t = dist / total;
            out[i % n] = std::array::from_fn(|k| {
                corners[j][k] + t * (corners[(j + 1) % 4][k] - corners[j][k])
            });
            dist += length(i);
        }
    }
    Ok(out)
}
fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}

fn enclosing_patch(
    source: &SectionPatch,
    a: [f64; 2],
    b: [f64; 2],
) -> Result<SectionPatch, &'static str> {
    const SCALE: i128 = 1_000_000_000;
    let d = source.degree as i128;
    let curve = source.controls.iter().step_by(2).collect::<Vec<_>>();
    // Give both ends of a radial edge the same homogeneous scale as its
    // matching tool endpoint. Degree-elevate the straight outside edge exactly.
    let point = |v: [f64; 2], w: i64| {
        [
            ((v[0] * SCALE as f64).round() as i128) * i128::from(w),
            ((v[1] * SCALE as f64).round() as i128) * i128::from(w),
            SCALE * i128::from(w),
        ]
    };
    let l = point(a, curve[0][2]);
    let r = point(b, curve[curve.len() - 1][2]);
    let mut net = Vec::<[i128; 3]>::new();
    for (i, h) in curve.iter().enumerate() {
        net.push(std::array::from_fn(|k| {
            (d - i as i128) * l[k] + i as i128 * r[k]
        }));
        net.push(h.map(|v| i128::from(v) * SCALE * d));
    }
    let g = net.iter().flatten().copied().fold(0, gcd).max(1);
    let controls = net
        .into_iter()
        .map(|h| {
            let [x, y, w] = h.map(|v| i64::try_from(v / g));
            Ok([
                x.map_err(|_| "stock coordinate overflow")?,
                y.map_err(|_| "stock coordinate overflow")?,
                w.map_err(|_| "stock weight overflow")?,
            ])
        })
        .collect::<Result<Vec<_>, &'static str>>()?;
    Ok(SectionPatch {
        degree: source.degree,
        controls,
    })
}

pub(crate) fn layers(
    r: &FaceProfileExtrusion,
    core: &[SectionPatch],
    ring: &[SectionPatch],
) -> Result<Vec<BracketBlock>, &'static str> {
    let tool: Vec<_> = core.iter().chain(ring.iter()).cloned().collect();
    let points = perimeter(ring, r.bounds)?;
    let outside = ring
        .iter()
        .enumerate()
        .map(|(i, p)| enclosing_patch(p, points[i], points[(i + 1) % points.len()]))
        .collect::<Result<Vec<_>, _>>()?;
    let mut blocks = Vec::new();
    let mut row = |start: f64, end: f64, inner: bool, outer: bool| -> Result<(), &'static str> {
        if end <= start {
            return Ok(());
        }
        let mut f = r.frame;
        f[0] = std::array::from_fn(|k| f[0][k] + start * f[3][k]);
        for patches in [
            inner.then_some(tool.as_slice()),
            outer.then_some(outside.as_slice()),
        ]
        .into_iter()
        .flatten()
        {
            for mut block in lift_sections(patches, end - start, f, false)? {
                block.id = format!("face_stock_{}", blocks.len());
                blocks.push(block);
            }
        }
        Ok(())
    };
    match r.operation.as_str() {
        "add" => {
            row(-r.thickness, 0., true, true)?;
            row(0., r.depth, true, false)?;
        }
        "cut" => {
            let bottom = (-r.depth).max(-r.thickness);
            row(-r.thickness, bottom, true, true)?;
            row(bottom, 0., false, true)?;
        }
        _ => return Err("unsupported face stock operation"),
    }
    Ok(blocks)
}
