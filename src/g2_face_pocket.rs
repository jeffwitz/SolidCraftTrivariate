//! Conforming S3 face pocket built from an audited S2 transverse-bore model.

use crate::{
    bracket::BracketBlock,
    constraint_system::Ratio,
    g2_cross_plate_json,
    incidence::{enumerate_faces, BlockSide},
    json::Json,
    s3_pocket_guide::radial_guide,
    s3_sector_plan::{certify_xmax_void, plan_sector, SectorPlan, SliceStatus},
    spatial_certificate::{certify_spatial, CertificateStatus},
};
use std::collections::BTreeMap;

type Error = &'static str;
type R = Ratio;
type Point2H = [R; 3];
type H4 = [R; 4];
type Output = (Vec<BracketBlock>, Vec<(BlockSide, String)>);

#[derive(Clone, Copy, PartialEq, Eq)]
enum RadialPart {
    Outer,
    Inner,
}

#[derive(Clone)]
struct BlockTag {
    parent_id: String,
    z_half: Option<usize>,
    z0: R,
    z1: R,
    u_interval: [R; 2],
    pocket_slice: bool,
    omit_outer: bool,
    wall_xmin: bool,
    radial: RadialPart,
}

fn int(value: i64) -> R {
    R::new(i128::from(value), 1).expect("integer ratio")
}

fn div(a: R, b: R) -> Result<R, Error> {
    if b.numerator() == 0 {
        return Err("division by zero in S3 exact geometry");
    }
    a.checked_mul(R::new(b.denominator(), b.numerator())?)
}

fn read_pair(value: &Json) -> Result<[i64; 2], Error> {
    let [a, b] = value.items()? else {
        return Err("expected a two-integer pair");
    };
    Ok([a.number()?, b.number()?])
}

fn read_ratio(value: &Json) -> Result<R, Error> {
    let [n, d] = read_pair(value)?;
    if d <= 0 {
        return Err("rational denominator must be positive");
    }
    R::new(i128::from(n), i128::from(d))
}

fn vector(value: &Json) -> Result<Vec<i64>, Error> {
    value.items()?.iter().map(Json::number).collect()
}

fn validate_parent(parent: &Json) -> Result<Output, Error> {
    let source = parent.field("source_recipe")?;
    if source.field("family")?.text()? != g2_cross_plate_json::FAMILY {
        return Err("S3 parent must be a delivered g2_cross_plate_v1 model");
    }
    g2_cross_plate_json::audit(parent)?;
    let blocks = crate::bracket_json::decode_blocks(parent)?;
    let faces = crate::bracket_json::decode_design_faces(parent)?;
    Ok((blocks, faces))
}

fn validate_operation(
    operation: &Json,
    support_z: R,
    stock_bounds: [[R; 2]; 2],
) -> Result<(R, R, R, R, R), Error> {
    operation.keys(&["support_feature", "frame", "direction", "depth_mm", "loop"])?;
    if operation.field("support_feature")?.text()? != "rect_x/zmin" {
        return Err("S3 pocket support must be rect_x/zmin");
    }
    let frame = operation.field("frame")?;
    frame.keys(&["origin_mm", "u", "v"])?;
    let origin = vector(frame.field("origin_mm")?)?;
    let support_z_i64 = i64::try_from(support_z.numerator() / support_z.denominator())
        .map_err(|_| "S3 support height exceeds integer frame range")?;
    if origin != [0, 0, support_z_i64]
        || vector(frame.field("u")?)? != [1, 0, 0]
        || vector(frame.field("v")?)? != [0, 1, 0]
        || vector(operation.field("direction")?)? != [0, 0, -1]
    {
        return Err("S3 pocket must use rect_x/zmin and run in the negative Z direction");
    }
    let depth = int(operation.field("depth_mm")?.number::<i64>()?);
    if !depth.cmp_value(&R::zero()).is_gt() || !depth.cmp_value(&support_z).is_lt() {
        return Err("S3 pocket depth must remain strictly inside the stock");
    }
    let loop_json = operation.field("loop")?;
    loop_json.keys(&["id", "primitive"])?;
    let primitive = loop_json.field("primitive")?;
    primitive.keys(&["kind", "min", "max"])?;
    if primitive.field("kind")?.text()? != "rectangle" {
        return Err("S3 pocket loop must be a native rectangle");
    }
    let [x_min_json, y_min_json] = primitive.field("min")?.items()? else {
        return Err("rectangle min must contain x and y");
    };
    let [x_max_json, y_max_json] = primitive.field("max")?.items()? else {
        return Err("rectangle max must contain x and y");
    };
    let x_min = read_ratio(x_min_json)?;
    let y_min = read_ratio(y_min_json)?;
    let x_max = read_ratio(x_max_json)?;
    let y_max = read_ratio(y_max_json)?;
    if !x_min.cmp_value(&stock_bounds[0][0]).is_ge()
        || !x_max.cmp_value(&stock_bounds[0][1]).is_le()
        || !y_min.cmp_value(&stock_bounds[1][0]).is_ge()
        || !y_max.cmp_value(&stock_bounds[1][1]).is_le()
        || !x_min.cmp_value(&x_max).is_lt()
        || !y_min.cmp_value(&y_max).is_lt()
    {
        return Err("S3 rectangle is outside the stock bounds");
    }
    Ok((x_min, x_max, y_min, y_max, depth))
}

fn metadata(block: &BracketBlock) -> Result<(usize, usize, usize, usize), Error> {
    // Stable S2 IDs encode sector, piece, Z layer, and source Bézier span.
    let fields: Vec<_> = block.id.split('_').collect();
    if fields.len() != 5 || fields[0] != "g2" || fields[1].len() < 2 || fields[2].len() < 2 {
        return Err("S3 requires canonical S2 block identifiers");
    }
    let sector = fields[1][1..].parse().map_err(|_| "invalid S2 sector ID")?;
    let piece = fields[2][1..].parse().map_err(|_| "invalid S2 piece ID")?;
    let layer = fields[3]
        .strip_prefix('z')
        .ok_or("invalid S2 layer ID")?
        .parse()
        .map_err(|_| "invalid S2 layer ID")?;
    let source = fields[4]
        .strip_prefix("src")
        .ok_or("invalid S2 span ID")?
        .parse()
        .map_err(|_| "invalid S2 span ID")?;
    Ok((sector, piece, layer, source))
}

fn flat_index(degree: [usize; 3], i: usize, j: usize, k: usize) -> usize {
    (i * (degree[1] + 1) + j) * (degree[2] + 1) + k
}

fn interpolate(a: H4, b: H4, t: R) -> Result<H4, Error> {
    let mut q = [R::zero(); 4];
    for d in 0..4 {
        q[d] = a[d].checked_add(b[d].checked_sub(a[d])?.checked_mul(t)?)?;
    }
    Ok(q)
}

fn split_axis(
    block: &BracketBlock,
    axis: usize,
    t: R,
) -> Result<(BracketBlock, BracketBlock), Error> {
    let degree = block.degrees[axis];
    if !(1..=3).contains(&degree) || t == R::zero() || t == R::one() {
        return Err("S3 subdivision must be interior to a linear or cubic axis");
    }
    let mut left = vec![[R::zero(); 4]; block.controls.len()];
    let mut right = vec![[R::zero(); 4]; block.controls.len()];
    let mut other = Vec::new();
    for a in 0..3 {
        if a != axis {
            other.push(a);
        }
    }
    let count0 = block.degrees[other[0]] + 1;
    let count1 = block.degrees[other[1]] + 1;
    for a in 0..count0 {
        for b in 0..count1 {
            let mut row: Vec<H4> = (0..=degree)
                .map(|s| {
                    let mut q = [0; 3];
                    q[axis] = s;
                    q[other[0]] = a;
                    q[other[1]] = b;
                    block.controls[flat_index(block.degrees, q[0], q[1], q[2])].map(int)
                })
                .collect();
            let mut l = vec![row[0]];
            let mut r = vec![*row.last().ok_or("empty S3 subdivision row")?];
            while row.len() > 1 {
                row = row
                    .windows(2)
                    .map(|pair| interpolate(pair[0], pair[1], t))
                    .collect::<Result<_, _>>()?;
                l.push(row[0]);
                r.push(*row.last().ok_or("empty S3 subdivision row")?);
            }
            r.reverse();
            for s in 0..=degree {
                for (target, value) in [(&mut left, l[s]), (&mut right, r[s])] {
                    let mut q = [0; 3];
                    q[axis] = s;
                    q[other[0]] = a;
                    q[other[1]] = b;
                    target[flat_index(block.degrees, q[0], q[1], q[2])] = value;
                }
            }
        }
    }
    let prefix = block.id.clone();
    Ok((
        BracketBlock {
            id: format!("{prefix}_cut{axis}a"),
            degrees: block.degrees,
            controls: integerize(left)?,
        },
        BracketBlock {
            id: format!("{prefix}_cut{axis}b"),
            degrees: block.degrees,
            controls: integerize(right)?,
        },
    ))
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn lcm(a: i128, b: i128) -> Result<i128, Error> {
    let g = i128::try_from(gcd(a.unsigned_abs(), b.unsigned_abs()))
        .map_err(|_| "S3 scale exceeds i128")?;
    (a / g).checked_mul(b).ok_or("S3 scale exceeds i128")
}

fn integerize(net: Vec<H4>) -> Result<Vec<[i64; 4]>, Error> {
    let mut scale = 1_i128;
    for p in &net {
        for c in p {
            scale = lcm(scale, c.denominator())?;
        }
    }
    let mut integers = Vec::with_capacity(net.len());
    let mut common = 0_u128;
    for p in net {
        let q = p.map(|v| {
            v.numerator()
                .checked_mul(scale / v.denominator())
                .ok_or("S3 control exceeds i128")
        });
        let q = q.into_iter().collect::<Result<Vec<_>, _>>()?;
        for v in &q {
            common = gcd(common, v.unsigned_abs());
        }
        integers.push(q);
    }
    let common = i128::try_from(common).map_err(|_| "S3 gauge exceeds i128")?;
    integers
        .into_iter()
        .map(|p| {
            p.into_iter()
                .map(|v| i64::try_from(v / common).map_err(|_| "S3 control exceeds i64"))
                .collect::<Result<Vec<_>, _>>()?
                .try_into()
                .map_err(|_| "invalid S3 homogeneous control")
        })
        .collect()
}

fn extract_curves(block: &BracketBlock) -> Result<(Vec<Point2H>, Vec<Point2H>, R, R), Error> {
    if block.degrees != [3, 1, 1] {
        return Err("S3 parent blocks must be cubic radial charts");
    }
    let mut outer = Vec::with_capacity(4);
    let mut inner = Vec::with_capacity(4);
    let mut z_bounds = None;
    for u in 0..4 {
        for v in 0..2 {
            let control = block.controls[flat_index(block.degrees, u, v, 0)];
            let h = [int(control[0]), int(control[1]), int(control[3])];
            if v == 0 {
                outer.push(h);
            } else {
                inner.push(h);
            }
            let lo = R::new(i128::from(control[2]), i128::from(control[3]))?;
            let high_control = block.controls[flat_index(block.degrees, u, v, 1)];
            let hi = R::new(i128::from(high_control[2]), i128::from(high_control[3]))?;
            if let Some((z0, z1)) = z_bounds {
                if z0 != lo || z1 != hi {
                    return Err("S3 source chart is not an axial extrusion");
                }
            } else {
                z_bounds = Some((lo, hi));
            }
        }
    }
    let (z0, z1) = z_bounds.ok_or("empty S3 parent chart")?;
    Ok((outer, inner, z0, z1))
}

fn patch_block(
    id: String,
    first: &[Point2H],
    second: &[Point2H],
    z0: R,
    z1: R,
) -> Result<BracketBlock, Error> {
    if first.len() != 4 || second.len() != 4 {
        return Err("S3 guide must be cubic");
    }
    let mut net = Vec::with_capacity(16);
    for i in 0..4 {
        for h in [first[i], second[i]] {
            for z in [z0, z1] {
                net.push([h[0], h[1], z.checked_mul(h[2])?, h[2]]);
            }
        }
    }
    Ok(BracketBlock {
        id,
        degrees: [3, 1, 1],
        controls: integerize(net)?,
    })
}

fn base_face_map(faces: &[(BlockSide, String)]) -> BTreeMap<(String, usize, bool), String> {
    faces
        .iter()
        .map(|(s, f)| ((s.block.clone(), s.axis, s.upper), f.clone()))
        .collect()
}

fn parent_label(
    tags: &BTreeMap<String, BlockTag>,
    old: &BTreeMap<(String, usize, bool), String>,
    id: &str,
    axis: usize,
    upper: bool,
) -> Option<String> {
    let tag = tags.get(id)?;
    old.get(&(tag.parent_id.clone(), axis, upper)).cloned()
}

fn homogeneous_face_y(block: &BracketBlock, side: &BlockSide) -> Result<Option<R>, Error> {
    if side.axis != 0 {
        return Ok(None);
    }
    let u = if side.upper { block.degrees[0] } else { 0 };
    let mut y = None;
    for v in 0..=block.degrees[1] {
        for w in 0..=block.degrees[2] {
            let p = block.controls[flat_index(block.degrees, u, v, w)];
            let value = R::new(i128::from(p[1]), i128::from(p[3]))?;
            if y.is_some_and(|old| old != value) {
                return Ok(None);
            }
            y = Some(value);
        }
    }
    Ok(y)
}

fn homogeneous_face_x_bounds(
    block: &BracketBlock,
    side: &BlockSide,
) -> Result<Option<[R; 2]>, Error> {
    if side.axis != 0 {
        return Ok(None);
    }
    let u = if side.upper { block.degrees[0] } else { 0 };
    let mut bounds: Option<[R; 2]> = None;
    for v in 0..=block.degrees[1] {
        for w in 0..=block.degrees[2] {
            let p = block.controls[flat_index(block.degrees, u, v, w)];
            let value = R::new(i128::from(p[0]), i128::from(p[3]))?;
            bounds = Some(bounds.map_or([value, value], |[lo, hi]| {
                [
                    if value.cmp_value(&lo).is_lt() {
                        value
                    } else {
                        lo
                    },
                    if value.cmp_value(&hi).is_gt() {
                        value
                    } else {
                        hi
                    },
                ]
            }));
        }
    }
    Ok(bounds)
}

/// Build the S3 model only from an audited, embedded S2 parent model.
pub fn construct(recipe: &Json) -> Result<Output, Error> {
    recipe.keys(&["family", "parent_model", "operation"])?;
    if recipe.field("family")?.text()? != "s3_face_pocket_v1" {
        return Err("unsupported S3 face-pocket family");
    }
    let parent = recipe.field("parent_model")?;
    // Reject non-S2 parents before invoking the comparatively expensive audit.
    if parent.field("source_recipe")?.field("family")?.text()? != g2_cross_plate_json::FAMILY {
        return Err("S3 parent must be a delivered g2_cross_plate_v1 model");
    }
    let (base_blocks, base_faces) = validate_parent(parent)?;
    if !base_faces
        .iter()
        .any(|(_, feature)| feature == "rect_x/zmin")
    {
        return Err("S3 parent has no rect_x/zmin support face");
    }

    let mut stock_bounds = [[R::zero(), R::zero()]; 2];
    let mut first_control = true;
    for block in &base_blocks {
        for p in &block.controls {
            let x = R::new(i128::from(p[0]), i128::from(p[3]))?;
            let y = R::new(i128::from(p[1]), i128::from(p[3]))?;
            if first_control {
                stock_bounds = [[x, x], [y, y]];
                first_control = false;
            } else {
                if x.cmp_value(&stock_bounds[0][0]).is_lt() {
                    stock_bounds[0][0] = x;
                }
                if x.cmp_value(&stock_bounds[0][1]).is_gt() {
                    stock_bounds[0][1] = x;
                }
                if y.cmp_value(&stock_bounds[1][0]).is_lt() {
                    stock_bounds[1][0] = y;
                }
                if y.cmp_value(&stock_bounds[1][1]).is_gt() {
                    stock_bounds[1][1] = y;
                }
            }
        }
    }
    let support = base_faces
        .iter()
        .find(|(_, feature)| feature == "rect_x/zmin")
        .ok_or("S3 parent has no rect_x/zmin support face")?;
    if support.0.axis != 2 || !support.0.upper {
        return Err("S3 support feature must be the upper-Z face of the recessed rectangle");
    }
    let support_block = base_blocks
        .iter()
        .find(|b| b.id == support.0.block)
        .ok_or("S3 support face references an unknown block")?;
    let (_, _, _, support_z) = extract_curves(support_block)?;
    let (x_min, x_max, y_min, y_max, depth) =
        validate_operation(recipe.field("operation")?, support_z, stock_bounds)?;
    let floor_z = support_z.checked_sub(depth)?;
    let support_blocks: std::collections::BTreeSet<_> = base_faces
        .iter()
        .filter(|(_, feature)| feature == "rect_x/zmin")
        .map(|(side, _)| side.block.clone())
        .collect();
    let mut right_outer = Vec::new();
    let mut right_inner = Vec::new();
    let mut target_plans = BTreeMap::<(usize, usize, usize), SectorPlan>::new();
    let mut auxiliary_inner_min: Option<R> = None;
    for block in &base_blocks {
        let (_, _, layer, _) = metadata(block)?;
        if layer != 0 || !support_blocks.contains(&block.id) {
            continue;
        }
        let (outer, inner, _, _) = extract_curves(block)?;
        for p in &inner {
            let x = div(p[0], p[2])?;
            auxiliary_inner_min = Some(auxiliary_inner_min.map_or(x, |old| {
                if x.cmp_value(&old).is_lt() {
                    x
                } else {
                    old
                }
            }));
        }
        let on_stock_xmax = outer
            .iter()
            .all(|p| div(p[0], p[2]).ok() == Some(stock_bounds[0][1]));
        if on_stock_xmax {
            right_outer.push(outer);
            right_inner.push(inner);
            continue;
        }
        if !outer
            .iter()
            .all(|p| div(p[0], p[2]).ok() == Some(stock_bounds[0][0]))
        {
            continue;
        }
        if !outer
            .iter()
            .zip(&inner)
            .all(|(a, b)| div(a[1], a[2]).ok() == div(b[1], b[2]).ok())
        {
            continue;
        }
        let (sector, piece, _, source_span) = metadata(block)?;
        let plan = plan_sector(
            &outer,
            &inner,
            stock_bounds[0][0],
            x_min,
            x_max,
            y_min,
            y_max,
        )?;
        if plan.remove_interval.is_some() {
            target_plans.insert((sector, piece, source_span), plan);
        }
    }
    certify_xmax_void(&right_outer, &right_inner, x_max, y_min, y_max)?;
    if target_plans.is_empty() {
        return Err("S3 rectangle does not meet a certified left-wing bore sector");
    }
    let auxiliary_guide_x = auxiliary_inner_min
        .ok_or("S3 parent has no inner profile")?
        .checked_add(stock_bounds[0][0])?
        .checked_mul(R::new(1, 2)?)?;
    let old_faces = base_face_map(&base_faces);
    let mut output = Vec::new();
    let mut tags = BTreeMap::<String, BlockTag>::new();
    for base in &base_blocks {
        let (sector, piece, layer, source_span) = metadata(base)?;
        let mut u_parts = vec![(
            base.clone(),
            None,
            [R::zero(), R::one()],
            false,
            false,
            false,
            false,
        )];
        if let Some(plan) = target_plans.get(&(sector, piece, source_span)) {
            let plan = plan.clone();
            let _interval = plan
                .remove_interval
                .ok_or("S3 pocket does not intersect its certified sector")?;
            let mut cuts = plan.u_cuts.clone();
            cuts.sort_by(|a, b| a.cmp_value(b));
            let mut remaining = base.clone();
            let mut previous = R::zero();
            let mut pieces = Vec::new();
            for (index, cut) in cuts
                .iter()
                .copied()
                .enumerate()
                .skip(1)
                .take(cuts.len().saturating_sub(2))
            {
                let local = div(cut.checked_sub(previous)?, R::one().checked_sub(previous)?)?;
                let (left, right) = split_axis(&remaining, 0, local)?;
                let inside = plan
                    .slices
                    .iter()
                    .any(|s| s.interval == [previous, cut] && s.status == SliceStatus::InsideBand);
                pieces.push((
                    left,
                    Some(index - 1),
                    [previous, cut],
                    inside,
                    plan.omit_outer,
                    plan.omit_inner,
                    plan.wall_xmin,
                ));
                remaining = right;
                previous = cut;
            }
            let last = *cuts.last().ok_or("S3 planner returned no cuts")?;
            let inside = plan
                .slices
                .iter()
                .any(|s| s.interval == [previous, last] && s.status == SliceStatus::InsideBand);
            pieces.push((
                remaining,
                Some(pieces.len()),
                [previous, last],
                inside,
                plan.omit_outer,
                plan.omit_inner,
                plan.wall_xmin,
            ));
            u_parts = pieces;
        }
        for (u_block, u_half, u_interval, pocket_slice, omit_outer, omit_inner, wall_xmin) in
            u_parts
        {
            let (outer, inner, z0, z1) = extract_curves(&u_block)?;
            let guide_x = if x_min == R::zero() {
                auxiliary_guide_x
            } else {
                x_min
            };
            let guide = radial_guide(&outer, &inner, stock_bounds[0][0], guide_x)?;
            if guide.len() != 4 {
                return Err("radial guide degree does not match its source chart");
            }
            let mut z_parts = vec![(z0, z1, None)];
            if layer == 0 {
                if z0 != R::zero() || z1 != support_z {
                    return Err("S3 expects the registered S2 lower axial layer");
                }
                z_parts = vec![(R::zero(), floor_z, Some(0)), (floor_z, support_z, Some(1))];
            }
            for (za, zb, z_half) in z_parts {
                for radial in [RadialPart::Outer, RadialPart::Inner] {
                    let omitted = pocket_slice
                        && z_half == Some(1)
                        && match radial {
                            RadialPart::Outer => omit_outer,
                            RadialPart::Inner => omit_inner,
                        };
                    if omitted {
                        continue;
                    }
                    let (a, b) = match radial {
                        RadialPart::Outer => (&outer, &guide),
                        RadialPart::Inner => (&guide, &inner),
                    };
                    let radial_name = if radial == RadialPart::Outer {
                        "o"
                    } else {
                        "i"
                    };
                    let half_name = u_half.map_or_else(|| "u".to_owned(), |i| format!("u{i}"));
                    let z_name = match z_half {
                        Some(0) => format!("z0_{}", floor_z.numerator() / floor_z.denominator()),
                        Some(1) => format!("z{}_4", floor_z.numerator() / floor_z.denominator()),
                        _ => "zbase".to_owned(),
                    };
                    let id = format!("s3_{radial_name}_{half_name}_{z_name}_{}", u_block.id);
                    let block = patch_block(id.clone(), a, b, za, zb)?;
                    let certificate =
                        certify_spatial(block.degrees, &block.controls, Default::default());
                    if certificate.status != CertificateStatus::Valid {
                        return Err("S3 chart lacks a positive Bernstein determinant certificate");
                    }
                    tags.insert(
                        id.clone(),
                        BlockTag {
                            parent_id: base.id.clone(),
                            z_half,
                            z0: za,
                            z1: zb,
                            u_interval,
                            pocket_slice,
                            omit_outer,
                            wall_xmin,
                            radial,
                        },
                    );
                    output.push(block);
                }
            }
        }
    }

    let incidence = enumerate_faces(&output)?;
    let by_id: BTreeMap<_, _> = output.iter().map(|b| (b.id.as_str(), b)).collect();
    let mut faces = Vec::new();
    for side in incidence.exterior {
        let tag = tags.get(&side.block).ok_or("missing S3 block provenance")?;
        let mut feature = None;
        // Preserve all original S2 exterior features through the exact chart split.
        match side.axis {
            0 => {
                let is_parent_end = (!side.upper && tag.u_interval[0] == R::zero())
                    || (side.upper && tag.u_interval[1] == R::one());
                if is_parent_end {
                    feature = parent_label(&tags, &old_faces, &side.block, 0, side.upper);
                }
                if feature.is_none()
                    && tag.z_half == Some(1)
                    && tag.z0 == floor_z
                    && tag.z1 == support_z
                {
                    let y = homogeneous_face_y(by_id[side.block.as_str()], &side)?;
                    let x_bounds = homogeneous_face_x_bounds(by_id[side.block.as_str()], &side)?;
                    let within_rectangle = x_bounds.is_some_and(|[lo, hi]| {
                        lo.cmp_value(&x_min).is_ge() && hi.cmp_value(&x_max).is_le()
                    });
                    if side.upper
                        && y == Some(y_max)
                        && within_rectangle
                        && (tag.radial == RadialPart::Inner || x_min == stock_bounds[0][0])
                    {
                        feature = Some("pocket/wall_ymax".to_owned());
                    } else if !side.upper
                        && y == Some(y_min)
                        && within_rectangle
                        && (tag.radial == RadialPart::Inner || x_min == stock_bounds[0][0])
                    {
                        feature = Some("pocket/wall_ymin".to_owned());
                    }
                }
            }
            1 => {
                let parent_side = match tag.radial {
                    RadialPart::Outer => !side.upper,
                    RadialPart::Inner => side.upper,
                };
                if parent_side {
                    feature = parent_label(&tags, &old_faces, &side.block, 1, side.upper);
                } else if tag.wall_xmin
                    && tag.pocket_slice
                    && tag.z_half == Some(1)
                    && tag.z0 == floor_z
                    && tag.z1 == support_z
                    && tag.radial == RadialPart::Outer
                    && side.upper
                {
                    feature = Some("pocket/wall_xmin".to_owned());
                }
            }
            2 => {
                if let Some(z_half) = tag.z_half {
                    if z_half == 0 && !side.upper {
                        feature = parent_label(&tags, &old_faces, &side.block, 2, false);
                    } else if z_half == 1 && side.upper {
                        feature = parent_label(&tags, &old_faces, &side.block, 2, true);
                    }
                    if feature.is_none()
                        && tag.pocket_slice
                        && z_half == 0
                        && side.upper
                        && (tag.radial == RadialPart::Inner || tag.omit_outer)
                    {
                        feature = Some("pocket/floor".to_owned());
                    }
                } else {
                    feature = parent_label(&tags, &old_faces, &side.block, 2, side.upper);
                }
            }
            _ => return Err("unexpected S3 exterior axis"),
        }
        faces.push((side, feature.ok_or("unclassified S3 exterior face")?));
    }
    Ok((output, faces))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registered_fixture_builds_a_conforming_certified_pocket() {
        let recipe = crate::json::parse(include_str!(
            "../tests/references/s3_face_pocket.recipe.json"
        ))
        .unwrap();
        let (blocks, faces) = construct(&recipe).unwrap();
        let incidence = enumerate_faces(&blocks).unwrap();
        assert_eq!(faces.len(), incidence.exterior.len());
        assert!(blocks.iter().all(
            |b| certify_spatial(b.degrees, &b.controls, Default::default()).status
                == CertificateStatus::Valid
        ));
        for feature in [
            "pocket/floor",
            "pocket/wall_xmin",
            "pocket/wall_ymin",
            "pocket/wall_ymax",
            "rect_x/zmin",
            "nurbs_z_hole/wall",
        ] {
            assert!(faces.iter().any(|(_, f)| f == feature), "missing {feature}");
        }
    }

    #[test]
    fn refuses_a_non_s2_parent_before_building() {
        let recipe = crate::json::parse(include_str!(
            "../tests/references/s3_face_pocket.recipe.json"
        ))
        .unwrap();
        let mut fields = match recipe {
            Json::Object(fields) => fields,
            _ => unreachable!(),
        };
        fields.insert("parent_model".into(), Json::Null);
        assert!(construct(&Json::Object(fields)).is_err());
    }

    #[test]
    fn open_side_pocket_removes_both_radial_rows_and_has_no_xmin_wall() {
        let recipe = crate::json::parse(include_str!(
            "../tests/references/s3_face_pocket_open_side.recipe.json"
        ))
        .unwrap();
        let (blocks, faces) = construct(&recipe).unwrap();
        assert_eq!(
            faces.len(),
            enumerate_faces(&blocks).unwrap().exterior.len()
        );
        assert!(faces.iter().any(|(_, f)| f == "pocket/floor"));
        assert!(faces.iter().any(|(_, f)| f == "pocket/wall_ymin"));
        assert!(faces.iter().any(|(_, f)| f == "pocket/wall_ymax"));
        assert!(!faces.iter().any(|(_, f)| f == "pocket/wall_xmin"));
    }
}
