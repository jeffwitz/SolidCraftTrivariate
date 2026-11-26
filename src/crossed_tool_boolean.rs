//! Exact ruled charts for the Boolean combinations of two crossed sketch tools.
use super::*;
use crate::json::object;

#[derive(Clone, Copy)]
pub enum Operation {
    Intersection,
    ProfileMinusRectangle,
    RectangleMinusProfile,
    Union,
}

impl Operation {
    pub fn parse(name: &str) -> Result<Self, Error> {
        match name {
            "intersection" => Ok(Self::Intersection),
            "a-minus-b" => Ok(Self::ProfileMinusRectangle),
            "b-minus-a" => Ok(Self::RectangleMinusProfile),
            "union" => Ok(Self::Union),
            _ => Err("unknown crossed-tool Boolean operation"),
        }
    }

    fn keeps(self, profile: bool, rectangle: bool) -> bool {
        match self {
            Self::Intersection => profile && rectangle,
            Self::ProfileMinusRectangle => profile && !rectangle,
            Self::RectangleMinusProfile => rectangle && !profile,
            Self::Union => profile || rectangle,
        }
    }
}

struct Band {
    lower: R,
    upper: R,
    left: HNet,
    right: HNet,
    traces: [[R; 2]; 2],
}

fn paired_bands(curve: &Curve, spans: &[BezierSpan]) -> Result<Vec<Band>, Error> {
    let mut cuts = vec![Vec::new(); spans.len()];
    for level in [curve.y_low, curve.y_high] {
        for root in crossings(spans, level)? {
            if root.local != R::zero() && root.local != R::one() {
                cuts[root.span].push(root.local);
            }
        }
    }
    let mut groups = BTreeMap::<(R, R), Vec<(HNet, [R; 2])>>::new();
    for (span, mut cuts) in spans.iter().zip(cuts) {
        cuts.sort_by(|a, b| compare(*a, *b));
        for (mut net, a, b) in split_piece_at(ratio_net(span), &cuts)? {
            let [start, end] = span.interval.map(int);
            let knot = |t: R| start.checked_add(end.checked_sub(start)?.checked_mul(t)?);
            let mut trace = [knot(a)?, knot(b)?];
            let mut ys: Vec<_> = net
                .iter()
                .map(|p| point(*p).map(|v| v[1]))
                .collect::<Result<_, _>>()?;
            if compare(ys[0], ys[3]).is_gt() {
                ys.reverse();
                net.reverse();
                trace.reverse();
            }
            if ys[0] == ys[3] || ys.windows(2).any(|p| compare(p[0], p[1]).is_gt()) {
                return Err("profile span lacks a monotone ordinate certificate");
            }
            let scale = div(R::one(), net[0][2])?;
            scale_net(&mut net, scale)?;
            groups.entry((ys[0], ys[3])).or_default().push((net, trace));
        }
    }
    let mut bands = Vec::new();
    for ((lower, upper), rows) in groups {
        let [(mut left, mut left_trace), (mut right, mut right_trace)]: [(HNet, [R; 2]); 2] = rows
            .try_into()
            .map_err(|_| "profile arcs require a different ordinate pairing")?;
        if left.iter().zip(&right).any(|(l, r)| l[1..] != r[1..]) {
            return Err("paired arcs have different ordinate parametrizations");
        }
        let widths: Vec<_> = left
            .iter()
            .zip(&right)
            .map(|(l, r)| compare(l[0], r[0]))
            .collect();
        if widths.iter().all(|w| !w.is_lt()) {
            std::mem::swap(&mut left, &mut right);
            std::mem::swap(&mut left_trace, &mut right_trace);
        } else if widths.iter().any(|w| w.is_gt()) {
            return Err("profile band width lacks a Bernstein sign certificate");
        }
        if left == right {
            return Err("zero width profile band");
        }
        bands.push(Band {
            lower,
            upper,
            traces: [left_trace, right_trace],
            left,
            right,
        });
    }
    bands.sort_by(|a, b| compare(a.lower, b.lower));
    if bands.windows(2).any(|pair| {
        pair[0].upper != pair[1].lower
            || [
                (&pair[0].left, &pair[1].left),
                (&pair[0].right, &pair[1].right),
            ]
            .iter()
            .any(|(a, b)| point(a[3]) != point(b[0]))
    }) {
        return Err("profile bands do not join");
    }
    let first = bands.first().ok_or("empty profile")?;
    let last = bands.last().unwrap();
    if point(first.left[0])? != point(first.right[0])?
        || point(last.left[3])? != point(last.right[3])?
    {
        return Err("profile bands do not close");
    }
    Ok(bands)
}

fn at_x(net: &HNet, x: i64) -> Result<HNet, Error> {
    net.iter()
        .map(|p| Ok([int(x).checked_mul(p[2])?, p[1], p[2]]))
        .collect()
}

// Only the two faces inherited from the transverse tool follow its sweep.
// The height of each is z_at_origin + (rise / stock_width) * x.
fn incline(
    controls: &mut [[i64; 4]],
    levels: [i64; 2],
    curve: &Curve,
    rise: i64,
) -> Result<(), Error> {
    if rise == 0 {
        return Ok(());
    }
    let width = i128::from(curve.bounds[1]);
    let mut scaled = Vec::with_capacity(controls.len());
    let mut divisor = 0;
    for (i, p) in controls.iter().enumerate() {
        let mut h = p.map(|c| i128::from(c) * width);
        if curve.z_cut.contains(&levels[i % 2]) {
            h[2] = h[2]
                .checked_add(i128::from(rise) * i128::from(p[0]))
                .ok_or("inclined coordinate overflow")?;
        }
        for c in h {
            divisor = gcd(divisor, c.unsigned_abs());
        }
        scaled.push(h);
    }
    let divisor = i128::try_from(divisor).map_err(|_| "inclined gauge overflow")?;
    for (p, h) in controls.iter_mut().zip(scaled) {
        for (c, value) in p.iter_mut().zip(h) {
            *c = i64::try_from(value / divisor).map_err(|_| "inclined control exceeds i64")?;
        }
    }
    Ok(())
}

fn interval(values: [R; 2]) -> Json {
    Json::Array(
        values
            .into_iter()
            .map(|r| {
                Json::Array(
                    [r.numerator(), r.denominator()]
                        .map(|v| Json::Number(v.to_string()))
                        .into(),
                )
            })
            .collect(),
    )
}

/// Shared arrangement of the two source tools. Reuse it for any Boolean selection.
pub struct CrossedTools {
    curve: Curve,
    bands: Vec<Band>,
    rise: i64,
    sketch_plan: Json,
}

impl CrossedTools {
    pub fn from_recipe(recipe: &Json) -> Result<Self, Error> {
        let (pair, rise) = if recipe.field("family")?.text()? == "inclined_crossed_tools_v1" {
            recipe.keys(&["family", "sketch_pair", "transverse_rise_mm"])?;
            (
                recipe.field("sketch_pair")?,
                recipe.field("transverse_rise_mm")?.number::<i64>()?,
            )
        } else {
            (recipe, 0)
        };
        let (parameters, sketch_plan) = crate::sketch_partition_plan::plan(pair)?;
        let (curve, spans) = parse_curve(&parameters)?;
        if i128::from(curve.z_cut[0]) + i128::from(rise) <= 0
            || i128::from(curve.z_cut[1]) + i128::from(rise) >= i128::from(curve.height)
        {
            return Err("inclined transverse tool exits a stock cap");
        }
        let bands = paired_bands(&curve, &spans)?;
        Ok(Self {
            curve,
            bands,
            rise,
            sketch_plan,
        })
    }

    /// Emit positive-interior maps ordered by ordinate, decreasing abscissa and height.
    pub fn charts(&self, operation: Operation) -> Result<Vec<BracketBlock>, Error> {
        let curve = &self.curve;
        let mut blocks = Vec::new();
        for (band_id, band) in self.bands.iter().enumerate() {
            let in_strip = !compare(band.lower, curve.y_low).is_lt()
                && !compare(band.upper, curve.y_high).is_gt();
            let stock_right = at_x(&band.right, curve.bounds[1])?;
            let stock_left = at_x(&band.left, curve.bounds[0])?;
            let rows = [&stock_left, &band.left, &band.right, &stock_right];
            for col in 0..3 {
                // Merge consecutive retained z intervals, avoiding propagation
                // of a cut through an entire retained profile column.
                let levels = [0, curve.z_cut[0], curve.z_cut[1], curve.height];
                let mut start = None;
                for layer in 0..3 {
                    let keep = operation.keeps(col == 1, in_strip && layer == 1);
                    if keep && start.is_none() {
                        start = Some(levels[layer]);
                    }
                    if let Some(z0) = start {
                        if !keep || layer == 2 {
                            let z1 = levels[layer + usize::from(keep)];
                            let mut controls = block_controls(rows[col + 1], rows[col], z0, z1)?;
                            incline(&mut controls, [z0, z1], curve, self.rise)?;
                            blocks.push(BracketBlock {
                                id: format!("band{band_id}_col{col}_z{z0}_{z1}"),
                                degrees: [3, 1, 1],
                                controls,
                            });
                            start = None;
                        }
                    }
                }
            }
        }
        Ok(blocks)
    }
}

pub fn construct(recipe: &Json, operation: Operation) -> Result<Vec<BracketBlock>, Error> {
    CrossedTools::from_recipe(recipe)?.charts(operation)
}

/// Persist volume charts separately from the face-conforming solver schema.
pub fn report(recipe: &Json, operation: &str) -> Result<Json, Error> {
    let arrangement = CrossedTools::from_recipe(recipe)?;
    let blocks = arrangement.charts(Operation::parse(operation)?)?;
    let sources = arrangement
        .bands
        .iter()
        .enumerate()
        .map(|(band_id, band)| {
            object([
                ("band", Json::Number(band_id.to_string())),
                ("y_interval_mm", interval([band.lower, band.upper])),
                ("left_profile_parameters", interval(band.traces[0])),
                ("right_profile_parameters", interval(band.traces[1])),
            ])
        })
        .collect();
    Ok(object([
        ("schema", Json::String("crossed_tool_charts_v1".into())),
        ("operation", Json::String(operation.into())),
        ("source_recipe", recipe.clone()),
        ("base_sketch_plan", arrangement.sketch_plan),
        ("profile_bands", Json::Array(sources)),
        ("blocks", crate::bracket_json::encode_blocks(blocks)),
    ]))
}
