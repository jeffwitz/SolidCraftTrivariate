//! Divide a rectangular wing into disjoint hosts for bores of different axes.
//! A separating slab boundary must clear every tool over its entire length.
use crate::compact_hosts::{holed_box, Bounds, HoledBox};
use crate::orthogonal_sketch::Bore;

pub(crate) fn build(
    bounds: Bounds,
    bores: &[&Bore],
    prefix: &str,
) -> Result<Option<HoledBox>, &'static str> {
    let Some(first) = bores.first() else {
        return Ok(None);
    };
    if bores.iter().all(|b| b.axis == first.axis) {
        return holed_box(first.axis, bounds, bores, prefix, false);
    }
    for axis in 0..3 {
        if bores.iter().any(|b| b.axis == axis) {
            continue;
        }
        let mut ordered = bores.to_vec();
        ordered.sort_by_key(|b| (b.center[axis], &b.id));
        for split in 1..ordered.len() {
            let (lower, upper) = ordered.split_at(split);
            let end = lower
                .iter()
                .map(|b| b.center[axis] + b.radius)
                .max()
                .unwrap();
            let start = upper
                .iter()
                .map(|b| b.center[axis] - b.radius)
                .min()
                .unwrap();
            let cut = end + (start - end) / 2;
            if !(end < cut && cut < start && bounds[axis][0] < cut && cut < bounds[axis][1]) {
                continue;
            }
            let mut a = bounds;
            let mut b = bounds;
            a[axis][1] = cut;
            b[axis][0] = cut;
            let Some(mut left) = build(a, lower, &format!("{prefix}_a{axis}_{cut}_lo"))? else {
                continue;
            };
            let Some(mut right) = build(b, upper, &format!("{prefix}_a{axis}_{cut}_hi"))? else {
                continue;
            };
            left.blocks.append(&mut right.blocks);
            left.physical.append(&mut right.physical);
            left.guards.append(&mut right.guards);
            return Ok(Some(left));
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn perpendicular_tools_require_a_separating_interval_in_the_common_cross_axis() {
        let z = Bore {
            id: "z_bore".into(),
            axis: 2,
            center: [25, 4, 20],
            radius: 2,
        };
        let mut y = Bore {
            id: "y_bore".into(),
            axis: 1,
            center: [40, 4, 20],
            radius: 4,
        };
        let bounds = [[14, 60], [0, 8], [0, 40]];
        // Allow clearance around the side-drilled bore within the host's XZ section.
        let candidate = build(bounds, &[&z, &y], "arm").unwrap().unwrap();
        assert_eq!(candidate.guards.len(), 2);
        y.center[0] = 25;
        assert!(build(bounds, &[&z, &y], "arm").unwrap().is_none());
    }
}
