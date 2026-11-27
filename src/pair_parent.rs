//! Exact recovery of one tensor Bézier parent from two adjacent restrictions.
//!
//! This is a pairwise algebraic test. It proves only that the two supplied
//! full-face restrictions are exact restrictions of one parent net; it does
//! not assign a historical feature identity or prove global material
//! coverage, regularity, or disjoint interiors.

use crate::{
    bracket::BracketBlock,
    constraint_system::Ratio,
    incidence::{BlockSide, Orientation},
    json::{object, Json},
};
use std::collections::BTreeMap;

type Result<T> = std::result::Result<T, &'static str>;
type HPoint = [Ratio; 4];
type Net = Vec<HPoint>;

fn string(value: impl Into<String>) -> Json {
    Json::String(value.into())
}

fn number(value: impl ToString) -> Json {
    Json::Number(value.to_string())
}

fn parse_side(value: &Json) -> Result<BlockSide> {
    let entries = value.items()?;
    let [block, side] = entries else {
        return Err("pair parent interface side requires block and face identifiers");
    };
    crate::mechanical_job::side(block.text()?, side.text()?)
}

fn parse_orientation(value: &Json) -> Result<Orientation> {
    value.keys(&["swap", "flip_first", "flip_second"])?;
    let get = |key| match value.field(key)? {
        Json::Bool(v) => Ok(*v),
        _ => Err("pair parent orientation fields must be booleans"),
    };
    Ok(Orientation {
        swap: get("swap")?,
        flip_first: get("flip_first")?,
        flip_second: get("flip_second")?,
    })
}

fn side_json(side: &BlockSide) -> Json {
    Json::Array(vec![
        string(side.block.clone()),
        string(format!(
            "{}{}",
            ["u", "v", "w"][side.axis],
            usize::from(side.upper)
        )),
    ])
}

fn orientation_json(orientation: Orientation) -> Json {
    object([
        ("swap", Json::Bool(orientation.swap)),
        ("flip_first", Json::Bool(orientation.flip_first)),
        ("flip_second", Json::Bool(orientation.flip_second)),
    ])
}

/// Audit all declared whole-face interfaces of an already audited model.
/// Recoverable pairs carry exact parent nets; every other pair remains in the
/// report as a refusal or an indeterminate arithmetic/ratio result.
pub fn recover_model_interfaces(model: &Json) -> Result<Json> {
    let mut blocks = BTreeMap::new();
    for block in crate::bracket_json::decode_blocks(model)? {
        validate(&block)?;
        if blocks.insert(block.id.clone(), block).is_some() {
            return Err("pair parent model has duplicate block identifiers");
        }
    }
    let interfaces = model.field("topology")?.field("interfaces")?.items()?;
    let mut rows = Vec::with_capacity(interfaces.len());
    let mut recovered = 0usize;
    let mut refused = 0usize;
    let mut indeterminate = 0usize;
    let mut axis_counts = [[0usize; 4]; 3];
    for (interface_index, entry) in interfaces.iter().enumerate() {
        let first_side = parse_side(entry.field("first")?)?;
        let second_side = parse_side(entry.field("second")?)?;
        let orientation = parse_orientation(entry.field("orientation")?)?;
        let first = blocks
            .get(&first_side.block)
            .ok_or("pair parent interface references an unknown first block")?;
        let second = blocks
            .get(&second_side.block)
            .ok_or("pair parent interface references an unknown second block")?;
        let attempt = recover_pair(
            first,
            first_side.clone(),
            second,
            second_side.clone(),
            orientation,
        );
        let (status, reason, recovery) = match attempt {
            Ok(parent) => {
                recovered += 1;
                axis_counts[first_side.axis][0] += 1;
                ("recovered", Json::Null, parent.document())
            }
            Err(reason)
                if reason.contains("indeterminate")
                    || reason.contains("overflow")
                    || reason.contains("exceeds i128")
                    || reason.contains("budget") =>
            {
                indeterminate += 1;
                axis_counts[first_side.axis][2] += 1;
                ("indeterminate", string(reason), Json::Null)
            }
            Err(reason) => {
                refused += 1;
                axis_counts[first_side.axis][1] += 1;
                ("refused", string(reason), Json::Null)
            }
        };
        rows.push(object([
            ("interface_index", number(interface_index)),
            ("first", side_json(&first_side)),
            ("second", side_json(&second_side)),
            ("orientation", orientation_json(orientation)),
            ("status", string(status)),
            ("reason", reason),
            ("recovery", recovery),
        ]));
    }
    let model_schema = model
        .field("schema")
        .and_then(Json::text)
        .ok()
        .map(string)
        .unwrap_or(Json::Null);
    let axis_summary = ["u", "v", "w"]
        .iter()
        .enumerate()
        .map(|(axis, name)| {
            object([
                ("axis", string(*name)),
                (
                    "interface_count",
                    number(axis_counts[axis].iter().sum::<usize>()),
                ),
                ("recovered_count", number(axis_counts[axis][0])),
                ("refused_count", number(axis_counts[axis][1])),
                ("indeterminate_count", number(axis_counts[axis][2])),
            ])
        })
        .collect();
    Ok(object([
        ("schema", string("pair_bezier_parent_report_v1")),
        ("model_schema", model_schema),
        ("claim_scope", string("pairwise exact mathematical support only; no historical feature provenance or global geometric certification")),
        ("interface_count", number(interfaces.len())),
        ("recovered_count", number(recovered)),
        ("refused_count", number(refused)),
        ("indeterminate_count", number(indeterminate)),
        ("counts_by_first_interface_axis", Json::Array(axis_summary)),
        ("interfaces", Json::Array(rows)),
    ]))
}

fn index(degrees: [usize; 3], p: [usize; 3]) -> usize {
    (p[0] * (degrees[1] + 1) + p[1]) * (degrees[2] + 1) + p[2]
}

fn checked_scale(point: HPoint, scale: Ratio) -> Result<HPoint> {
    let mut out = [Ratio::zero(); 4];
    for c in 0..4 {
        out[c] = point[c].checked_mul(scale)?;
    }
    Ok(out)
}

fn checked_add_scaled(a: HPoint, sa: Ratio, b: HPoint, sb: Ratio) -> Result<HPoint> {
    let mut out = [Ratio::zero(); 4];
    for c in 0..4 {
        out[c] = a[c].checked_mul(sa)?.checked_add(b[c].checked_mul(sb)?)?;
    }
    Ok(out)
}

fn validate(block: &BracketBlock) -> Result<()> {
    let expected = block
        .degrees
        .iter()
        .try_fold(1usize, |n, d| n.checked_mul(d.checked_add(1)?))
        .ok_or("pair parent control count overflow")?;
    if block.id.is_empty()
        || block.degrees.iter().any(|d| !(1..=3).contains(d))
        || block.controls.len() != expected
    {
        return Err("invalid pair parent input block");
    }
    if block.controls.iter().any(|p| p[3] <= 0) {
        return Err("pair parent requires positive homogeneous weights");
    }
    Ok(())
}

fn elevate_axis(mut net: Net, degrees: &mut [usize; 3], axis: usize, target: usize) -> Result<Net> {
    if axis >= 3 || target > 3 || target < degrees[axis] {
        return Err("unsupported pair parent degree elevation");
    }
    while degrees[axis] < target {
        let old_degrees = *degrees;
        let mut new_degrees = old_degrees;
        new_degrees[axis] += 1;
        let mut new_net = vec![[Ratio::zero(); 4]; (0..3).map(|a| new_degrees[a] + 1).product()];
        let other: Vec<_> = (0..3).filter(|&a| a != axis).collect();
        for a in 0..=old_degrees[other[0]] {
            for b in 0..=old_degrees[other[1]] {
                for k in 0..=new_degrees[axis] {
                    let mut at = [0; 3];
                    at[other[0]] = a;
                    at[other[1]] = b;
                    at[axis] = k;
                    let destination = index(new_degrees, at);
                    if k == 0 {
                        at[axis] = 0;
                        new_net[destination] = net[index(old_degrees, at)];
                    } else if k == new_degrees[axis] {
                        at[axis] = old_degrees[axis];
                        new_net[destination] = net[index(old_degrees, at)];
                    } else {
                        let alpha = Ratio::new(k as i128, new_degrees[axis] as i128)?;
                        let complement = Ratio::one().checked_sub(alpha)?;
                        at[axis] = k - 1;
                        let previous = net[index(old_degrees, at)];
                        at[axis] = k;
                        let current = net[index(old_degrees, at)];
                        new_net[destination] =
                            checked_add_scaled(previous, alpha, current, complement)?;
                    }
                }
            }
        }
        *degrees = new_degrees;
        net = new_net;
    }
    Ok(net)
}

fn reverse_axis(net: Net, degrees: [usize; 3], axis: usize) -> Net {
    let mut out = vec![[Ratio::zero(); 4]; net.len()];
    for i in 0..=degrees[0] {
        for j in 0..=degrees[1] {
            for k in 0..=degrees[2] {
                let output = [i, j, k];
                let mut source = output;
                source[axis] = degrees[axis] - output[axis];
                out[index(degrees, output)] = net[index(degrees, source)];
            }
        }
    }
    out
}

fn input_net(block: &BracketBlock) -> Net {
    block
        .controls
        .iter()
        .map(|p| p.map(|x| Ratio::new(i128::from(x), 1).expect("integer ratio")))
        .collect()
}

fn orient_second(
    block: &BracketBlock,
    first_axis: usize,
    second_side: &BlockSide,
    orientation: Orientation,
) -> Result<(Net, [usize; 3])> {
    if first_axis >= 3 || second_side.axis >= 3 {
        return Err("pair parent interface axis outside 0..3");
    }
    let mut source_axis = [usize::MAX; 3];
    source_axis[first_axis] = second_side.axis;
    let first_free: Vec<_> = (0..3).filter(|&a| a != first_axis).collect();
    let second_free: Vec<_> = (0..3).filter(|&a| a != second_side.axis).collect();
    let mut flip = [false; 3];
    if second_side.upper {
        flip[first_axis] = true;
    }
    for output_tangent in 0..2 {
        let source_tangent = if orientation.swap {
            1 - output_tangent
        } else {
            output_tangent
        };
        source_axis[first_free[output_tangent]] = second_free[source_tangent];
        flip[first_free[output_tangent]] = if orientation.swap {
            if output_tangent == 0 {
                orientation.flip_second
            } else {
                orientation.flip_first
            }
        } else if output_tangent == 0 {
            orientation.flip_first
        } else {
            orientation.flip_second
        };
    }
    if source_axis.iter().any(|&a| a >= 3) {
        return Err("invalid pair parent axis permutation");
    }
    let degrees = std::array::from_fn(|a| block.degrees[source_axis[a]]);
    let source = input_net(block);
    let mut out = vec![[Ratio::zero(); 4]; source.len()];
    for i in 0..=degrees[0] {
        for j in 0..=degrees[1] {
            for k in 0..=degrees[2] {
                let output = [i, j, k];
                let mut input = [0; 3];
                for axis in 0..3 {
                    let source_id = source_axis[axis];
                    input[source_id] = if flip[axis] {
                        block.degrees[source_id] - output[axis]
                    } else {
                        output[axis]
                    };
                }
                out[index(degrees, output)] = source[index(block.degrees, input)];
            }
        }
    }
    Ok((out, degrees))
}

fn face_points(net: &Net, degrees: [usize; 3], axis: usize, upper: bool) -> Vec<HPoint> {
    let free: Vec<_> = (0..3).filter(|&a| a != axis).collect();
    let fixed = if upper { degrees[axis] } else { 0 };
    let mut result = Vec::new();
    for i in 0..=degrees[free[0]] {
        for j in 0..=degrees[free[1]] {
            let mut p = [0; 3];
            p[axis] = fixed;
            p[free[0]] = i;
            p[free[1]] = j;
            result.push(net[index(degrees, p)]);
        }
    }
    result
}

fn common_gauge(a: &[HPoint], b: &[HPoint]) -> Result<Ratio> {
    if a.len() != b.len() || a.is_empty() {
        return Err("pair parent interface trace dimensions differ");
    }
    let gauge = a[0][3].checked_div(b[0][3])?;
    if gauge <= Ratio::zero() {
        return Err("pair parent interface gauge must be positive");
    }
    for (left, right) in a.iter().zip(b) {
        if *left != checked_scale(*right, gauge)? {
            return Err("pair parent traces differ beyond a constant homogeneous gauge");
        }
    }
    Ok(gauge)
}

fn infer_split(a: &Net, b: &Net, degrees: [usize; 3], axis: usize) -> Result<Ratio> {
    let p = degrees[axis];
    let free: Vec<_> = (0..3).filter(|&x| x != axis).collect();
    let mut t: Option<Ratio> = None;
    for i in 0..=degrees[free[0]] {
        for j in 0..=degrees[free[1]] {
            let mut seam = [0; 3];
            seam[free[0]] = i;
            seam[free[1]] = j;
            let mut left_inner = seam;
            left_inner[axis] = p - 1;
            seam[axis] = p;
            let left_delta = std::array::from_fn::<_, 4, _>(|c| {
                a[index(degrees, seam)][c].checked_sub(a[index(degrees, left_inner)][c])
            });
            let left_delta = left_delta.into_iter().collect::<Result<Vec<_>>>()?;
            let mut right_inner = seam;
            right_inner[axis] = 1;
            seam[axis] = 0;
            let right_delta = std::array::from_fn::<_, 4, _>(|c| {
                b[index(degrees, right_inner)][c].checked_sub(b[index(degrees, seam)][c])
            });
            let right_delta = right_delta.into_iter().collect::<Result<Vec<_>>>()?;
            for c in 0..4 {
                let sum = left_delta[c].checked_add(right_delta[c])?;
                if sum != Ratio::zero() {
                    let candidate = left_delta[c].checked_div(sum)?;
                    if t.is_some_and(|known| known != candidate) {
                        return Err("pair parent first jets imply inconsistent interval ratios");
                    }
                    t = Some(candidate);
                } else if left_delta[c] != Ratio::zero() || right_delta[c] != Ratio::zero() {
                    return Err("pair parent first jets imply no positive interval ratio");
                }
            }
        }
    }
    let t = t.ok_or("pair parent split ratio is indeterminate from zero first jets")?;
    if t.cmp_value(&Ratio::zero()).is_le() || t.cmp_value(&Ratio::one()).is_ge() {
        return Err("pair parent split ratio is outside the open unit interval");
    }
    for i in 0..=degrees[free[0]] {
        for j in 0..=degrees[free[1]] {
            let mut seam = [0; 3];
            seam[free[0]] = i;
            seam[free[1]] = j;
            let mut left_inner = seam;
            left_inner[axis] = p - 1;
            seam[axis] = p;
            let left = a[index(degrees, seam)];
            let left_previous = a[index(degrees, left_inner)];
            let mut right_inner = seam;
            right_inner[axis] = 1;
            seam[axis] = 0;
            let right = b[index(degrees, right_inner)];
            let right_start = b[index(degrees, seam)];
            for c in 0..4 {
                let da = left[c].checked_sub(left_previous[c])?;
                let db = right[c].checked_sub(right_start[c])?;
                if da != t.checked_mul(da.checked_add(db)?)? {
                    return Err("pair parent first jets are not one positive scaled trace");
                }
            }
        }
    }
    Ok(t)
}

fn binomial(n: usize, k: usize) -> i128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    (1..=k).fold(1_i128, |v, i| v * (n - k + i) as i128 / i as i128)
}

fn power(base: Ratio, exponent: usize) -> Result<Ratio> {
    let mut result = Ratio::one();
    for _ in 0..exponent {
        result = result.checked_mul(base)?;
    }
    Ok(result)
}

fn recover_parent(first: &Net, degrees: [usize; 3], axis: usize, t: Ratio) -> Result<Net> {
    let mut parent = vec![[Ratio::zero(); 4]; first.len()];
    let p = degrees[axis];
    let complement = Ratio::one().checked_sub(t)?;
    let free: Vec<_> = (0..3).filter(|&x| x != axis).collect();
    for i in 0..=degrees[free[0]] {
        for j in 0..=degrees[free[1]] {
            let mut at = [0; 3];
            at[free[0]] = i;
            at[free[1]] = j;
            let mut line = Vec::with_capacity(p + 1);
            for q in 0..=p {
                at[axis] = q;
                line.push(first[index(degrees, at)]);
            }
            let mut recovered = vec![[Ratio::zero(); 4]; p + 1];
            recovered[0] = line[0];
            for k in 1..=p {
                let mut known = [Ratio::zero(); 4];
                for (j, previous) in recovered.iter().take(k).enumerate() {
                    let coefficient = Ratio::new(binomial(k, j), 1)?
                        .checked_mul(power(t, j)?)?
                        .checked_mul(power(complement, k - j)?)?;
                    for c in 0..4 {
                        known[c] = known[c].checked_add(coefficient.checked_mul(previous[c])?)?;
                    }
                }
                let denominator = power(t, k)?;
                for c in 0..4 {
                    recovered[k][c] = line[k][c].checked_sub(known[c])?.checked_div(denominator)?;
                }
            }
            for (k, control) in recovered.iter().copied().enumerate() {
                at[axis] = k;
                parent[index(degrees, at)] = control;
            }
        }
    }
    Ok(parent)
}

fn split_line(line: &[HPoint], t: Ratio) -> Result<(Vec<HPoint>, Vec<HPoint>)> {
    let p = line.len() - 1;
    let complement = Ratio::one().checked_sub(t)?;
    let mut level = line.to_vec();
    let mut left = vec![line[0]];
    let mut right = vec![line[p]];
    while level.len() > 1 {
        let mut next = Vec::with_capacity(level.len() - 1);
        for pair in level.windows(2) {
            next.push(checked_add_scaled(pair[0], complement, pair[1], t)?);
        }
        left.push(next[0]);
        right.push(*next.last().expect("nonempty split level"));
        level = next;
    }
    right.reverse();
    Ok((left, right))
}

fn split_net(parent: &Net, degrees: [usize; 3], axis: usize, t: Ratio) -> Result<(Net, Net)> {
    let mut left = vec![[Ratio::zero(); 4]; parent.len()];
    let mut right = left.clone();
    let free: Vec<_> = (0..3).filter(|&a| a != axis).collect();
    for i in 0..=degrees[free[0]] {
        for j in 0..=degrees[free[1]] {
            let mut at = [0; 3];
            at[free[0]] = i;
            at[free[1]] = j;
            let mut line = Vec::with_capacity(degrees[axis] + 1);
            for k in 0..=degrees[axis] {
                at[axis] = k;
                line.push(parent[index(degrees, at)]);
            }
            let (l, r) = split_line(&line, t)?;
            for k in 0..=degrees[axis] {
                at[axis] = k;
                left[index(degrees, at)] = l[k];
                right[index(degrees, at)] = r[k];
            }
        }
    }
    Ok((left, right))
}

#[derive(Debug, Clone)]
pub struct PairParentRecovery {
    /// New mathematical support identity, not a historical feature identity.
    pub source_ids: [String; 2],
    pub degrees: [usize; 3],
    pub controls: Vec<[Ratio; 4]>,
    /// First interval's fraction of the recovered parent interval.
    pub split_parameter: Ratio,
    /// Exact first:second parameter-interval length ratio.
    pub first_to_second_interval_ratio: Ratio,
    /// Positive scale applied to the second input net to match the first gauge.
    pub second_to_first_gauge: Ratio,
}

impl PairParentRecovery {
    pub fn document(&self) -> Json {
        let exact = |x: Ratio| {
            Json::Array(vec![
                Json::String(x.numerator().to_string()),
                Json::String(x.denominator().to_string()),
            ])
        };
        object([
            ("schema",Json::String("pair_bezier_parent_v1".into())),
            ("identity_kind",Json::String("mathematical_parent_recovery_not_feature_provenance".into())),
            ("source_block_ids",Json::Array(self.source_ids.iter().cloned().map(Json::String).collect())),
            ("degrees",Json::Array(self.degrees.iter().map(|x|Json::Number(x.to_string())).collect())),
            ("controls_homogeneous",Json::Array(self.controls.iter().map(|p|Json::Array(p.iter().copied().map(exact).collect())).collect())),
            ("split_parameter",exact(self.split_parameter)),
            ("first_to_second_interval_ratio",exact(self.first_to_second_interval_ratio)),
            ("second_to_first_gauge",exact(self.second_to_first_gauge)),
            ("verification",Json::String("exact full-net de Casteljau restriction after orientation and positive constant gauge".into())),
        ])
    }
}

/// Recover one tensor Bézier map from two audited, adjacent full-face blocks.
///
/// The first block is oriented so its declared interface is the upper face;
/// the second is oriented so its interface is the lower face. `orientation`
/// maps the second tangential trace into the first trace convention. Tangential
/// and normal degrees are harmonized by exact Bernstein elevation (up to cubic).
/// The first-jet ratio is inferred from the complete homogeneous face net, then
/// every parent control line is recovered and both block nets are reconstructed
/// by exact de Casteljau restriction. Unsupported or inconsistent pairs refuse.
pub fn recover_pair(
    first: &BracketBlock,
    first_side: BlockSide,
    second: &BracketBlock,
    second_side: BlockSide,
    orientation: Orientation,
) -> Result<PairParentRecovery> {
    validate(first)?;
    validate(second)?;
    if first.id == second.id || first_side.block != first.id || second_side.block != second.id {
        return Err("pair parent side identities do not match input blocks");
    }
    if first_side.axis >= 3 || second_side.axis >= 3 {
        return Err("pair parent interface axis outside 0..3");
    }
    let axis = first_side.axis;
    let mut a = input_net(first);
    let mut degree_a = first.degrees;
    if !first_side.upper {
        a = reverse_axis(a, degree_a, axis);
    }
    let (mut b, mut degree_b) = orient_second(second, axis, &second_side, orientation)?;
    let mut target = degree_a;
    for ax in 0..3 {
        target[ax] = degree_a[ax].max(degree_b[ax]);
    }
    for (ax, &target_degree) in target.iter().enumerate() {
        a = elevate_axis(a, &mut degree_a, ax, target_degree)?;
        b = elevate_axis(b, &mut degree_b, ax, target_degree)?;
    }
    let face_a = face_points(&a, target, axis, true);
    let face_b = face_points(&b, target, axis, false);
    let gauge = common_gauge(&face_a, &face_b)?;
    for point in &mut b {
        *point = checked_scale(*point, gauge)?;
    }
    let t = infer_split(&a, &b, target, axis)?;
    let parent = recover_parent(&a, target, axis, t)?;
    if parent.iter().any(|p| p[3] <= Ratio::zero()) {
        return Err("recovered pair parent has a nonpositive homogeneous weight");
    }
    let (reconstructed_a, reconstructed_b) = split_net(&parent, target, axis, t)?;
    if reconstructed_a != a || reconstructed_b != b {
        return Err("full restrictions do not reproduce both supplied homogeneous nets");
    }
    let ratio = t.checked_div(Ratio::one().checked_sub(t)?)?;
    Ok(PairParentRecovery {
        source_ids: [first.id.clone(), second.id.clone()],
        degrees: target,
        controls: parent,
        split_parameter: t,
        first_to_second_interval_ratio: ratio,
        second_to_first_gauge: gauge,
    })
}
