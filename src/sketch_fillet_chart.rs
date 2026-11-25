//! Two exact rational Bézier charts for a square corner with a circular fillet.
//!
//! The circular arc uses the rational parametrization
//! `(-2t, -(1-t^2), 1+t^2)` on the unit quarter circle. Splitting at `t=1/2`
//! therefore keeps all homogeneous controls rational and places the shared
//! point at `(L-4r/5, L-3r/5)`. Both arc traces lie on `v=1`; the exterior
//! square paths lie on `v=0`.

use crate::{bracket::BracketBlock, planar::certify_planar_refined};

type HPoint = [i128; 3];

fn point(x: i128, y: i128, w: i128) -> HPoint {
    [x, y, w]
}

fn checked_add(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_add(b).ok_or("fillet chart arithmetic overflow")
}

fn checked_mul(a: i128, b: i128) -> Result<i128, &'static str> {
    a.checked_mul(b).ok_or("fillet chart arithmetic overflow")
}

fn translated(local: [i128; 2], weight: i128, origin: [i32; 2]) -> Result<HPoint, &'static str> {
    Ok([
        checked_add(local[0], checked_mul(i128::from(origin[0]), weight)?)?,
        checked_add(local[1], checked_mul(i128::from(origin[1]), weight)?)?,
        weight,
    ])
}

fn arc_controls(
    origin: [i32; 2],
    size: i32,
    radius: i32,
) -> Result<[[HPoint; 3]; 2], &'static str> {
    let l = i128::from(size);
    let r = i128::from(radius);
    // De Casteljau restriction at t=1/2, multiplied by four to clear the
    // denominators. The second chart reverses its half so both charts have
    // positive orientation while their shared arc endpoint remains exact.
    let left_forward = [
        point(4 * l, 4 * (l - r), 4),
        point(4 * l - 2 * r, 4 * l - 4 * r, 4),
        point(5 * l - 4 * r, 5 * l - 3 * r, 5),
    ];
    let right_forward = [
        point(5 * l - 4 * r, 5 * l - 3 * r, 5),
        point(6 * l - 6 * r, 6 * l - 2 * r, 6),
        point(8 * l - 8 * r, 8 * l, 8),
    ];
    let east = left_forward.into_iter().rev().collect::<Vec<_>>();
    let north = right_forward.into_iter().rev().collect::<Vec<_>>();
    let transform = |p: HPoint| translated([p[0], p[1]], p[2], origin);
    let east = [
        transform(east[0])?,
        transform(east[1])?,
        transform(east[2])?,
    ];
    let north = [
        transform(north[0])?,
        transform(north[1])?,
        transform(north[2])?,
    ];
    Ok([east, north])
}

fn outer_controls(size: i32, origin: [i32; 2], north: bool) -> Result<[HPoint; 3], &'static str> {
    let l = i128::from(size);
    let local = if north {
        // The north chart runs from (0,L) to (0,0). Collinear positive
        // controls define the same physical segment with a free parameter map.
        [point(0, 8 * l, 8), point(0, 4 * l, 8), point(0, 0, 5)]
    } else {
        // The east chart runs from (0,0) to (L,0).
        [point(0, 0, 5), point(2 * l, 0, 4), point(4 * l, 0, 4)]
    };
    Ok([
        translated([local[0][0], local[0][1]], local[0][2], origin)?,
        translated([local[1][0], local[1][1]], local[1][2], origin)?,
        translated([local[2][0], local[2][1]], local[2][2], origin)?,
    ])
}

fn make_block(
    id: String,
    outer: [HPoint; 3],
    inner: [HPoint; 3],
    z: [i32; 2],
) -> Result<BracketBlock, &'static str> {
    let mut planar = Vec::with_capacity(6);
    for u in 0..3 {
        planar.push([
            i64::try_from(outer[u][0]).map_err(|_| "fillet chart control exceeds i64")?,
            i64::try_from(outer[u][1]).map_err(|_| "fillet chart control exceeds i64")?,
            i64::try_from(outer[u][2]).map_err(|_| "fillet chart control exceeds i64")?,
        ]);
        planar.push([
            i64::try_from(inner[u][0]).map_err(|_| "fillet chart control exceeds i64")?,
            i64::try_from(inner[u][1]).map_err(|_| "fillet chart control exceeds i64")?,
            i64::try_from(inner[u][2]).map_err(|_| "fillet chart control exceeds i64")?,
        ]);
    }
    certify_planar_refined([2, 1], &planar).map_err(|error| match error {
        crate::planar::PlanarError::Invalid => "fillet chart invalid or nonpositive corner",
        crate::planar::PlanarError::Indeterminate => {
            "fillet chart Bernstein determinant indeterminate"
        }
    })?;

    let mut controls = Vec::with_capacity(12);
    for u in 0..3 {
        for v in 0..2 {
            let h = if v == 0 { outer[u] } else { inner[u] };
            for height in z {
                controls.push([
                    i64::try_from(h[0]).map_err(|_| "fillet chart control exceeds i64")?,
                    i64::try_from(h[1]).map_err(|_| "fillet chart control exceeds i64")?,
                    i64::try_from(checked_mul(i128::from(height), h[2])?)
                        .map_err(|_| "fillet chart control exceeds i64")?,
                    i64::try_from(h[2]).map_err(|_| "fillet chart control exceeds i64")?,
                ]);
            }
        }
    }
    Ok(BracketBlock {
        id,
        degrees: [2, 1, 1],
        controls,
    })
}

/// Build two positive-orientation rational charts for a square corner fillet.
///
/// `origin` and `size` define the closed square `[x,x+size] × [y,y+size]`;
/// the removed material is the quarter disk of radius `radius` at its upper
/// right corner. The output is extruded between `z[0]` and `z[1]` in mm.
/// Both circular traces are on `v=1`. The diagonal trace is exactly shared.
pub fn fillet_charts(
    origin: [i32; 2],
    size: i32,
    radius: i32,
    z: [i32; 2],
    id: &str,
) -> Result<Vec<BracketBlock>, &'static str> {
    if size <= 0 || radius <= 0 || radius >= size || z[0] >= z[1] {
        return Err("invalid fillet chart dimensions");
    }
    if id.is_empty()
        || id.len() > 48
        || !id
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, b'_' | b'-'))
    {
        return Err("invalid fillet chart id");
    }
    let [east_arc, north_arc] = arc_controls(origin, size, radius)?;
    let east_outer = outer_controls(size, origin, false)?;
    let north_outer = outer_controls(size, origin, true)?;
    Ok(vec![
        make_block(format!("{id}_east_fillet"), east_outer, east_arc, z)?,
        make_block(format!("{id}_north_fillet"), north_outer, north_arc, z)?,
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RationalBezierMap;

    fn map(block: &BracketBlock) -> RationalBezierMap {
        RationalBezierMap::new(
            block.degrees,
            block
                .controls
                .iter()
                .map(|point| point.map(|v| v as f64))
                .collect(),
        )
        .unwrap()
    }

    #[test]
    fn exact_arc_split_and_shared_diagonal_survive_translation_and_scale() {
        for (origin, size, radius, z) in [
            ([0, 0], 14, 6, [0, 15]),
            ([-3, 5], 12, 3, [-2, 7]),
            ([10, -4], 31, 17, [0, 8]),
            ([0, 0], 9, 8, [2, 19]),
        ] {
            let blocks = fillet_charts(origin, size, radius, z, "corner").unwrap();
            let [east, north] = blocks.as_slice() else {
                panic!("two charts expected")
            };
            assert_eq!(east.degrees, [2, 1, 1]);
            assert_eq!(north.degrees, [2, 1, 1]);
            assert_eq!(east.id, "corner_east_fillet");
            assert_eq!(north.id, "corner_north_fillet");
            for chart in [east, north] {
                let planar: Vec<_> = chart
                    .controls
                    .iter()
                    .step_by(2)
                    .map(|p| [p[0], p[1], p[3]])
                    .collect();
                assert!(certify_planar_refined([2, 1], &planar).is_ok());
            }
            let a = map(east);
            let b = map(north);
            let corner = [f64::from(origin[0] + size), f64::from(origin[1] + size)];
            for chart in [&a, &b] {
                for j in 0..=32 {
                    let u = f64::from(j) / 32.0;
                    let point = chart.evaluate([u, 1.0, 0.63]).unwrap().position;
                    let error =
                        (point[0] - corner[0]).hypot(point[1] - corner[1]) - f64::from(radius);
                    assert!(error.abs() <= 1e-12);
                    let expected_z = f64::from(z[0]) + 0.63 * f64::from(z[1] - z[0]);
                    assert!((point[2] - expected_z).abs() < 1e-12);
                }
            }
            for j in 0..=16 {
                let v = f64::from(j) / 16.0;
                let p = a.evaluate([0.0, v, 0.37]).unwrap().position;
                let q = b.evaluate([1.0, v, 0.37]).unwrap().position;
                for axis in 0..3 {
                    assert!((p[axis] - q[axis]).abs() < 1e-12);
                }
            }
            let mid = a.evaluate([0.0, 1.0, 0.0]).unwrap().position;
            assert!((mid[0] - (corner[0] - 4.0 * f64::from(radius) / 5.0)).abs() < 1e-12);
            assert!((mid[1] - (corner[1] - 3.0 * f64::from(radius) / 5.0)).abs() < 1e-12);
            for (chart, u, expected) in [(east, 1.0, 0), (north, 0.0, 1)] {
                let line_start = map(chart).evaluate([u, 0.0, 0.0]).unwrap().position;
                let line_end = map(chart).evaluate([u, 1.0, 0.0]).unwrap().position;
                for k in 0..=8 {
                    let t = f64::from(k) / 8.0;
                    let p = map(chart).evaluate([u, t, 0.0]).unwrap().position;
                    if expected == 0 {
                        assert!((p[0] - f64::from(origin[0] + size)).abs() < 1e-12);
                        let y = line_start[1] + t * (line_end[1] - line_start[1]);
                        assert!((p[1] - y).abs() < 1e-12);
                    } else {
                        assert!((p[1] - f64::from(origin[1] + size)).abs() < 1e-12);
                        let x = line_start[0] + t * (line_end[0] - line_start[0]);
                        assert!((p[0] - x).abs() < 1e-12);
                    }
                }
            }
        }
    }

    #[test]
    fn rejects_empty_or_out_of_range_dimensions_and_ids() {
        assert!(fillet_charts([0, 0], 4, 0, [0, 1], "x").is_err());
        assert!(fillet_charts([0, 0], 4, 4, [0, 1], "x").is_err());
        assert!(fillet_charts([0, 0], 4, 2, [1, 1], "x").is_err());
        assert!(fillet_charts([0, 0], 4, 2, [0, 1], "bad/id").is_err());
    }
}
