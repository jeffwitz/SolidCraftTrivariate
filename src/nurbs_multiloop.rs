//! Exact partition of a rectangular stock by x-monotone NURBS loop profiles.
//!
//! The admitted profiles have synchronized lower and upper rational Bezier
//! spans. Between disjoint loop projections, their endpoint ordinates are
//! joined by straight rational spans. Each x segment is partitioned into
//! three y bands and then extruded through the global z event layers.

use crate::{
    bracket::BracketBlock,
    incidence::{enumerate_faces, BlockSide},
    json::Json,
    nurbs_curve::BezierSpan,
    nurbs_extrusion::extrude_strips,
    nurbs_sketch_json::spans,
};
use std::collections::{BTreeSet, HashMap};

type R = crate::constraint_system::Ratio;
type Error = &'static str;
type Construction = (Vec<BracketBlock>, Vec<(BlockSide, String)>);

#[derive(Clone)]
struct Loop {
    id: String,
    lower: Vec<BezierSpan>,
    upper: Vec<BezierSpan>,
    z: [i64; 2],
    x: [R; 2],
}

#[derive(Clone)]
struct Segment {
    lower: Vec<BezierSpan>,
    upper: Vec<BezierSpan>,
    left_loop: Option<usize>,
    right_loop: Option<usize>,
    contour: Option<usize>,
    part: usize,
}

#[derive(Clone)]
struct Tag {
    segment: usize,
    band: usize,
    layer: usize,
}

fn less(a: R, b: R) -> Result<bool, Error> {
    Ok(a.checked_sub(b)?.numerator() < 0)
}

fn integer(value: i64) -> R {
    R::new(i128::from(value), 1).expect("integer ratio")
}

fn coordinate(point: [i64; 3], axis: usize) -> Result<R, Error> {
    if point[2] <= 0 {
        return Err("NURBS control weight must be positive");
    }
    R::new(i128::from(point[axis]), i128::from(point[2]))
}

fn endpoints(curve: &[BezierSpan]) -> Result<[R; 2], Error> {
    let first = curve.first().ok_or("empty NURBS profile")?;
    let last = curve.last().ok_or("empty NURBS profile")?;
    Ok([
        coordinate(first.controls[0], 0)?,
        coordinate(*last.controls.last().ok_or("empty NURBS profile")?, 0)?,
    ])
}

fn valid_id(id: &str) -> bool {
    !id.is_empty()
        && id.len() <= 40
        && id
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-'))
}

fn parse_pair(value: &Json) -> Result<[i64; 2], Error> {
    let [a, b] = value.items()? else {
        return Err("expected a two-value interval");
    };
    Ok([a.number()?, b.number()?])
}

fn constant_curve(curve: &[BezierSpan], y: i64) -> Result<Vec<BezierSpan>, Error> {
    curve
        .iter()
        .map(|span| {
            let controls = span
                .controls
                .iter()
                .map(|p| {
                    Ok([
                        p[0],
                        i64::try_from(
                            i128::from(y)
                                .checked_mul(i128::from(p[2]))
                                .ok_or("stock boundary control overflow")?,
                        )
                        .map_err(|_| "stock boundary control exceeds i64")?,
                        p[2],
                    ])
                })
                .collect::<Result<Vec<_>, Error>>()?;
            Ok(BezierSpan {
                interval: span.interval,
                controls,
            })
        })
        .collect()
}

fn rational_line(a: [R; 2], b: [R; 2]) -> Result<BezierSpan, Error> {
    fn gcd(mut a: i128, mut b: i128) -> i128 {
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a.abs()
    }
    let mut scale = 1_i128;
    for p in [a, b] {
        for v in p {
            let d = v.denominator();
            let g = gcd(scale, d);
            scale = (scale / g)
                .checked_mul(d)
                .ok_or("gap denominator overflow")?;
        }
    }
    let point = |p: [R; 2]| -> Result<[i64; 3], Error> {
        let x = p[0]
            .numerator()
            .checked_mul(scale / p[0].denominator())
            .ok_or("gap x control overflow")?;
        let y = p[1]
            .numerator()
            .checked_mul(scale / p[1].denominator())
            .ok_or("gap y control overflow")?;
        Ok([
            i64::try_from(x).map_err(|_| "gap x control exceeds i64")?,
            i64::try_from(y).map_err(|_| "gap y control exceeds i64")?,
            i64::try_from(scale).map_err(|_| "gap weight exceeds i64")?,
        ])
    };
    Ok(BezierSpan {
        interval: [0, 1],
        controls: vec![point(a)?, point(b)?],
    })
}

fn line_pair(
    x0: R,
    x1: R,
    y0: R,
    y1: R,
    v0: R,
    v1: R,
) -> Result<(Vec<BezierSpan>, Vec<BezierSpan>), Error> {
    if !less(x0, x1)? {
        return Err("loop x projections must be strictly separated");
    }
    let x0p = [x0, y0];
    let x1p = [x1, y1];
    Ok((
        vec![rational_line(x0p, x1p)?],
        vec![rational_line([x0, v0], [x1, v1])?],
    ))
}

fn shifted(block: &mut BracketBlock, z: i64) -> Result<(), Error> {
    for p in &mut block.controls {
        let offset = i128::from(z)
            .checked_mul(i128::from(p[3]))
            .ok_or("translated z control overflow")?;
        p[2] = i64::try_from(
            i128::from(p[2])
                .checked_add(offset)
                .ok_or("translated z control overflow")?,
        )
        .map_err(|_| "translated z control exceeds i64")?;
    }
    Ok(())
}

fn active(loop_: &Loop, z: [i64; 2]) -> bool {
    loop_.z[0] <= z[0] && z[1] <= loop_.z[1]
}

/// Construct the complete three-band partition and return one design label
/// for every exterior block side.
pub fn construct(recipe: &Json) -> Result<Construction, Error> {
    recipe.keys(&["family", "bounds_xy_mm", "height_mm", "loops"])?;
    if recipe.field("family")?.text()? != "nurbs_multiloop_v1" {
        return Err("unsupported NURBS multiloop family");
    }
    let bounds = recipe.field("bounds_xy_mm")?.items()?;
    let [xb, yb] = bounds else {
        return Err("stock bounds must be x and y intervals");
    };
    let [xmin, xmax] = parse_pair(xb)?;
    let [ymin, ymax] = parse_pair(yb)?;
    let height: i64 = recipe.field("height_mm")?.number()?;
    if xmin >= xmax || ymin >= ymax || height <= 0 || height > i64::from(i32::MAX) {
        return Err("invalid NURBS multiloop stock bounds");
    }
    let mut loops = Vec::new();
    let mut ids = BTreeSet::new();
    for item in recipe.field("loops")?.items()? {
        item.keys(&["id", "lower", "upper", "z_interval_mm"])?;
        let id = item.field("id")?.text()?.to_owned();
        if !valid_id(&id) || !ids.insert(id.clone()) {
            return Err("invalid or duplicate NURBS contour identifier");
        }
        let lower = spans(item.field("lower")?)?;
        let upper = spans(item.field("upper")?)?;
        if lower.is_empty() || lower.len() != upper.len() {
            return Err("NURBS contour boundaries need synchronized spans");
        }
        for (a, b) in lower.iter().zip(&upper) {
            if a.interval != b.interval {
                return Err("NURBS contour boundaries need matching knots");
            }
        }
        let x = endpoints(&lower)?;
        if !less(x[0], x[1])? {
            return Err("NURBS contour must increase strictly in x");
        }
        if !less(integer(xmin), x[0])? || !less(x[1], integer(xmax))? {
            return Err("interior contour must be strictly inside stock");
        }
        for p in lower.iter().chain(&upper).flat_map(|s| &s.controls) {
            if p[2] <= 0
                || i128::from(p[1]) <= i128::from(ymin) * i128::from(p[2])
                || i128::from(p[1]) >= i128::from(ymax) * i128::from(p[2])
            {
                return Err("NURBS contour must remain strictly inside stock in y");
            }
        }
        let z = parse_pair(item.field("z_interval_mm")?)?;
        if z[0] < 0 || z[0] >= z[1] || z[1] > height {
            return Err("z interval must have positive length inside stock");
        }
        loops.push(Loop {
            id,
            lower,
            upper,
            z,
            x,
        });
    }
    if loops.is_empty() {
        return Err("NURBS multiloop count outside supported range");
    }
    loops.sort_by(|a, b| a.x[0].cmp_value(&b.x[0]));
    for pair in loops.windows(2) {
        if !less(pair[0].x[1], pair[1].x[0])? {
            return Err("interior loop x-projections must be strictly disjoint");
        }
    }

    // Build a physical-x schedule. Original NURBS Bezier spans are retained;
    // gaps use exact linear spans interpolating both endpoint ordinates.
    let mut segments = Vec::new();
    let first = &loops[0];
    if less(integer(xmin), first.x[0])? {
        let l0 = coordinate(first.lower[0].controls[0], 1)?;
        let u0 = coordinate(first.upper[0].controls[0], 1)?;
        let (lower, upper) = line_pair(integer(xmin), first.x[0], l0, l0, u0, u0)?;
        segments.push(Segment {
            lower,
            upper,
            left_loop: None,
            right_loop: Some(0),
            contour: None,
            part: 0,
        });
    }
    for (index, loop_) in loops.iter().enumerate() {
        for (part, (lower_span, upper_span)) in loop_.lower.iter().zip(&loop_.upper).enumerate() {
            segments.push(Segment {
                lower: vec![lower_span.clone()],
                upper: vec![upper_span.clone()],
                left_loop: if part == 0 { Some(index) } else { None },
                right_loop: if part + 1 == loop_.lower.len() {
                    Some(index)
                } else {
                    None
                },
                contour: Some(index),
                part,
            });
        }
        if let Some(next) = loops.get(index + 1) {
            let l0 = coordinate(*loop_.lower.last().unwrap().controls.last().unwrap(), 1)?;
            let u0 = coordinate(*loop_.upper.last().unwrap().controls.last().unwrap(), 1)?;
            let l1 = coordinate(next.lower[0].controls[0], 1)?;
            let u1 = coordinate(next.upper[0].controls[0], 1)?;
            let (lower, upper) = line_pair(loop_.x[1], next.x[0], l0, l1, u0, u1)?;
            segments.push(Segment {
                lower,
                upper,
                left_loop: Some(index),
                right_loop: Some(index + 1),
                contour: None,
                part: 0,
            });
        }
    }
    let last_index = loops.len() - 1;
    let last = &loops[last_index];
    if less(last.x[1], integer(xmax))? {
        let l1 = coordinate(*last.lower.last().unwrap().controls.last().unwrap(), 1)?;
        let u1 = coordinate(*last.upper.last().unwrap().controls.last().unwrap(), 1)?;
        let (lower, upper) = line_pair(last.x[1], integer(xmax), l1, l1, u1, u1)?;
        segments.push(Segment {
            lower,
            upper,
            left_loop: Some(last_index),
            right_loop: None,
            contour: None,
            part: 0,
        });
    }

    // Reindex every scheduled span onto one common increasing parameter
    // interval and certify the complete chain in one exact extrusion call.
    let mut scheduled_lower = Vec::new();
    let mut scheduled_upper = Vec::new();
    let mut interval = 0_i64;
    for segment in &mut segments {
        if segment.lower.len() != segment.upper.len() {
            return Err("scheduled NURBS boundaries need paired spans");
        }
        for (lower, upper) in segment.lower.iter_mut().zip(&mut segment.upper) {
            let end = interval.checked_add(1).ok_or("span schedule overflow")?;
            lower.interval = [interval, end];
            upper.interval = [interval, end];
            scheduled_lower.push(lower.clone());
            scheduled_upper.push(upper.clone());
            interval = end;
        }
    }
    extrude_strips(&scheduled_lower, &scheduled_upper, 1)?;

    let mut events = BTreeSet::from([0_i64, height]);
    for l in &loops {
        events.insert(l.z[0]);
        events.insert(l.z[1]);
    }
    let events: Vec<_> = events.into_iter().collect();
    let mut blocks = Vec::new();
    let mut tags = HashMap::<String, Tag>::new();
    for (layer, zrange) in events.windows(2).enumerate() {
        let z = [zrange[0], zrange[1]];
        if z[0] == z[1] {
            continue;
        }
        let dz = i32::try_from(z[1] - z[0]).map_err(|_| "z layer exceeds i32")?;
        for (si, segment) in segments.iter().enumerate() {
            let lower_stock = constant_curve(&segment.lower, ymin)?;
            let upper_stock = constant_curve(&segment.upper, ymax)?;
            let bands = [
                extrude_strips(&lower_stock, &segment.lower, dz)?,
                extrude_strips(&segment.lower, &segment.upper, dz)?,
                extrude_strips(&segment.upper, &upper_stock, dz)?,
            ];
            for (band, band_blocks) in bands.into_iter().enumerate() {
                if band == 1 && segment.contour.is_some_and(|i| active(&loops[i], z)) {
                    continue;
                }
                for (sp, mut block) in band_blocks.into_iter().enumerate() {
                    let contour_id = segment
                        .contour
                        .map(|i| loops[i].id.as_str())
                        .unwrap_or("gap");
                    block.id = format!(
                        "ml_s{si}_p{}_y{band}_z{}_b{sp}_{contour_id}",
                        segment.part, layer
                    );
                    shifted(&mut block, z[0])?;
                    tags.insert(
                        block.id.clone(),
                        Tag {
                            segment: si,
                            band,
                            layer,
                        },
                    );
                    blocks.push(block);
                }
            }
        }
    }
    let incidence = enumerate_faces(&blocks)?;
    let mut faces = Vec::new();
    for side in incidence.exterior {
        let tag = tags
            .get(&side.block)
            .ok_or("missing multiloop block provenance")?;
        let segment = &segments[tag.segment];
        let z0 = events[tag.layer];
        let z1 = events[tag.layer + 1];
        let label = match side.axis {
            0 if side.upper && segment.right_loop.is_none() => "stock:xmax".to_owned(),
            0 if !side.upper && segment.left_loop.is_none() => "stock:xmin".to_owned(),
            0 if tag.band == 1 && side.upper => {
                let i = segment.right_loop.ok_or("unattributed right loop end")?;
                if !active(&loops[i], [z0, z1]) {
                    return Err("unexpected exposed x face away from an active loop");
                }
                format!("{}:end_lower", loops[i].id)
            }
            0 if tag.band == 1 && !side.upper => {
                let i = segment.left_loop.ok_or("unattributed left loop end")?;
                if !active(&loops[i], [z0, z1]) {
                    return Err("unexpected exposed x face away from an active loop");
                }
                format!("{}:end_upper", loops[i].id)
            }
            1 if !side.upper && tag.band == 0 => "stock:ymin".to_owned(),
            1 if side.upper && tag.band == 2 => "stock:ymax".to_owned(),
            1 if side.upper && tag.band == 0 => {
                let i = segment
                    .contour
                    .ok_or("unexpected exposed lower y face in a gap")?;
                if !active(&loops[i], [z0, z1]) {
                    return Err("unexpected exposed lower wall outside its active interval");
                }
                format!("{}:wall_lower", loops[i].id)
            }
            1 if !side.upper && tag.band == 2 => {
                let i = segment
                    .contour
                    .ok_or("unexpected exposed upper y face in a gap")?;
                if !active(&loops[i], [z0, z1]) {
                    return Err("unexpected exposed upper wall outside its active interval");
                }
                format!("{}:wall_upper", loops[i].id)
            }
            2 if !side.upper && z0 == 0 => "stock:zmin".to_owned(),
            2 if side.upper && z1 == height => "stock:zmax".to_owned(),
            2 => {
                if tag.band != 1 {
                    return Err("unexpected exposed multiloop layer face");
                }
                let i = segment.contour.ok_or("unattributed pocket floor")?;
                let l = &loops[i];
                if side.upper && z1 == l.z[0] {
                    format!("{}:floor", l.id)
                } else if !side.upper && z0 == l.z[1] {
                    format!("{}:ceiling", l.id)
                } else {
                    return Err("unattributed z-layer face");
                }
            }
            _ => return Err("unattributed multiloop exterior face"),
        };
        faces.push((side, label));
    }
    Ok((blocks, faces))
}
