//! Exact moving circular bore in a two-layer square plate.
//!
//! The triangular reference patches are converted to three biquadratic
//! quadrilaterals by `triangle_quads`; this module only creates the exact
//! planar triangle nets and extrudes their returned controls.

use crate::{bracket::BracketBlock, triangle_quads::quadrangulate};
use std::convert::TryFrom;

const HALF_WIDTH_UM: i128 = 10_000;
const POCKET_LEFT_UM: i128 = -6_000;
const POCKET_RIGHT_UM: i128 = -2_000;
const RADIUS_UM: i128 = 2_000;
const HALF_HEIGHT_UM: i128 = 4_000;

#[derive(Clone, Copy)]
struct Rational {
    n: i128,
    d: i128,
}

impl Rational {
    fn new(mut n: i128, mut d: i128) -> Result<Self, &'static str> {
        if d == 0 {
            return Err("zero rational denominator");
        }
        if d < 0 {
            n = n.checked_neg().ok_or("rational overflow")?;
            d = d.checked_neg().ok_or("rational overflow")?;
        }
        let g = gcd(n, d);
        Ok(Self { n: n / g, d: d / g })
    }
    fn add(self, rhs: Self) -> Result<Self, &'static str> {
        Self::new(
            self.n
                .checked_mul(rhs.d)
                .and_then(|v| rhs.n.checked_mul(self.d).and_then(|w| v.checked_add(w)))
                .ok_or("rational overflow")?,
            self.d.checked_mul(rhs.d).ok_or("rational overflow")?,
        )
    }
    fn half(self) -> Result<Self, &'static str> {
        Self::new(self.n, self.d.checked_mul(2).ok_or("rational overflow")?)
    }
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.max(1)
}

fn lcm(a: i128, b: i128) -> Result<i128, &'static str> {
    (a / gcd(a, b))
        .checked_mul(b)
        .ok_or("homogeneous gauge overflow")
}

#[derive(Clone, Copy)]
struct Point {
    x_um: i128,
    y_um: i128,
    half_angle: Option<(i128, i128)>,
    bounded_y_fm: Option<i128>,
}

impl Point {
    fn plain(x_um: i128, y_um: i128) -> Self {
        Self {
            x_um,
            y_um,
            half_angle: None,
            bounded_y_fm: None,
        }
    }

    fn circle(x_um: i128, y_um: i128, cx_um: i128) -> Result<Self, &'static str> {
        let x = (x_um - cx_um).abs();
        if x > RADIUS_UM {
            return Err("circle point outside radius");
        }
        let radicand = RADIUS_UM
            .checked_mul(RADIUS_UM)
            .and_then(|v| x.checked_mul(x).and_then(|q| v.checked_sub(q)))
            .ok_or("circle radicand overflow")?;
        let y = isqrt(radicand);
        if y.checked_mul(y) != Some(radicand) || y != y_um.abs() {
            return Err("circle intersection is not an integer micrometre");
        }
        let den = RADIUS_UM + x;
        Ok(Self {
            x_um,
            y_um,
            half_angle: Some((y, den)),
            bounded_y_fm: None,
        })
    }
}

pub(crate) fn isqrt(n: i128) -> i128 {
    if n <= 0 {
        return 0;
    }
    let mut x = n;
    let mut y = (x + 1) / 2;
    while y < x {
        x = y;
        y = (x + n / x) / 2;
    }
    x
}

type RH = [Rational; 3];

fn vertex_h(p: Point, cx_um: i128) -> Result<RH, &'static str> {
    let (num, den) = p.half_angle.unwrap_or((0, 1));
    let q = den
        .checked_mul(den)
        .and_then(|v| num.checked_mul(num).and_then(|n| v.checked_add(n)))
        .ok_or("circle weight overflow")?;
    let d = den.checked_mul(den).ok_or("circle weight overflow")?;
    let _ = cx_um;
    Ok([
        Rational::new(p.x_um.checked_mul(q).ok_or("coordinate overflow")?, d)?,
        Rational::new(p.y_um.checked_mul(q).ok_or("coordinate overflow")?, d)?,
        Rational::new(1_000_i128.checked_mul(q).ok_or("weight overflow")?, d)?,
    ])
}

fn arc_mid_h(a: Point, b: Point, cx_um: i128) -> Result<RH, &'static str> {
    let (na, da) = a.half_angle.ok_or("arc endpoint has no angle")?;
    let (nb, db) = b.half_angle.ok_or("arc endpoint has no angle")?;
    let den = da.checked_mul(db).ok_or("arc denominator overflow")?;
    let product = na
        .checked_mul(nb)
        .and_then(|v| da.checked_mul(db).and_then(|d| v.checked_add(d)))
        .ok_or("arc weight overflow")?;
    let dx = da
        .checked_mul(db)
        .and_then(|v| na.checked_mul(nb).and_then(|n| v.checked_sub(n)))
        .ok_or("arc coordinate overflow")?;
    let dy = na
        .checked_mul(db)
        .and_then(|v| nb.checked_mul(da).and_then(|n| v.checked_add(n)))
        .ok_or("arc coordinate overflow")?;
    let sx: i128 = if a.x_um + b.x_um >= 2 * cx_um { 1 } else { -1 };
    let sy: i128 = if a.y_um + b.y_um >= 0 { 1 } else { -1 };
    let xnum = cx_um
        .checked_mul(product)
        .and_then(|v| {
            sx.checked_mul(RADIUS_UM)
                .and_then(|r| r.checked_mul(dx))
                .and_then(|r| v.checked_add(r))
        })
        .ok_or("arc coordinate overflow")?;
    Ok([
        Rational::new(xnum, den)?,
        Rational::new(
            sy.checked_mul(RADIUS_UM)
                .and_then(|v| v.checked_mul(dy))
                .ok_or("arc coordinate overflow")?,
            den,
        )?,
        Rational::new(
            1_000_i128
                .checked_mul(product)
                .ok_or("arc weight overflow")?,
            den,
        )?,
    ])
}

fn midpoint(a: RH, b: RH) -> Result<RH, &'static str> {
    Ok([
        a[0].add(b[0])?.half()?,
        a[1].add(b[1])?.half()?,
        a[2].add(b[2])?.half()?,
    ])
}

fn triangle_controls(points: [Point; 3], cx_um: i128) -> Result<[[i64; 3]; 6], &'static str> {
    let mut p = points;
    let cross = (p[1].x_um - p[0].x_um) * (p[2].y_um - p[0].y_um)
        - (p[1].y_um - p[0].y_um) * (p[2].x_um - p[0].x_um);
    if cross == 0 {
        return Err("degenerate macro-triangle");
    }
    if cross < 0 {
        p.swap(1, 2);
    }
    if p.iter().any(|v| v.bounded_y_fm.is_some()) {
        let points = p.map(|v| {
            [
                v.x_um * 1_000_000_000,
                v.bounded_y_fm.unwrap_or(v.y_um * 1_000_000_000),
            ]
        });
        return crate::travelling_bore_bounded::triangle_net(
            points,
            p.map(|v| v.bounded_y_fm.is_some()),
            cx_um,
        );
    }
    let v = [
        vertex_h(p[0], cx_um)?,
        vertex_h(p[1], cx_um)?,
        vertex_h(p[2], cx_um)?,
    ];
    let edge = |i: usize, j: usize| -> Result<RH, &'static str> {
        if p[i].half_angle.is_some() && p[j].half_angle.is_some() {
            arc_mid_h(p[i], p[j], cx_um)
        } else {
            midpoint(v[i], v[j])
        }
    };
    let h = [v[0], v[1], v[2], edge(0, 1)?, edge(1, 2)?, edge(2, 0)?];
    let mut gauge = 1_i128;
    for value in h.iter().flatten() {
        gauge = lcm(gauge, value.d)?;
    }
    let mut result = [[0_i64; 3]; 6];
    for i in 0..6 {
        for d in 0..3 {
            let n = h[i][d]
                .n
                .checked_mul(gauge / h[i][d].d)
                .ok_or("triangle common gauge overflow")?;
            result[i][d] = i64::try_from(n).map_err(|_| "triangle coefficient exceeds i64")?;
        }
    }
    Ok(result)
}

fn append_triangle(
    blocks: &mut Vec<BracketBlock>,
    points: [Point; 3],
    cx_um: i128,
    base_id: &str,
    z0_mm: i64,
    z1_mm: i64,
) -> Result<(), &'static str> {
    let tri = triangle_controls(points, cx_um)?;
    let patches = quadrangulate(tri)?;
    for (corner, planar) in patches.iter().enumerate() {
        crate::planar::certify_planar([2, 2], planar)
            .map_err(|_| "planar triangle restriction is not certified regular")?;
        let mut controls = Vec::with_capacity(18);
        for h in planar {
            for z in [z0_mm, z1_mm] {
                controls.push([
                    h[0],
                    h[1],
                    i64::try_from(
                        i128::from(z)
                            .checked_mul(i128::from(h[2]))
                            .ok_or("extrusion overflow")?,
                    )
                    .map_err(|_| "extruded coordinate exceeds i64")?,
                    h[2],
                ]);
            }
        }
        blocks.push(BracketBlock {
            id: format!("{base_id}_q{corner}"),
            degrees: [2, 2, 1],
            controls,
        });
    }
    Ok(())
}

fn grid_x(cx_um: i128) -> Vec<i128> {
    let mut xs = vec![
        -HALF_WIDTH_UM,
        HALF_WIDTH_UM,
        POCKET_LEFT_UM,
        POCKET_RIGHT_UM,
    ];
    if cx_um > -HALF_WIDTH_UM && cx_um < HALF_WIDTH_UM {
        xs.push(cx_um);
    }
    xs.sort_unstable();
    xs.dedup();
    xs
}

fn rect_triangles(x0: i128, x1: i128, y0: i128, y1: i128) -> [[Point; 3]; 2] {
    let a = Point::plain(x0, y0);
    let b = Point::plain(x1, y0);
    let c = Point::plain(x1, y1);
    let d = Point::plain(x0, y1);
    [[a, b, c], [a, c, d]]
}

fn append_rect(
    blocks: &mut Vec<BracketBlock>,
    cx_um: i128,
    bounds: [[i128; 2]; 3],
    id: &str,
) -> Result<(), &'static str> {
    for (triangle, points) in rect_triangles(bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1])
        .into_iter()
        .enumerate()
    {
        append_triangle(
            blocks,
            points,
            cx_um,
            &format!("{id}_t{triangle}"),
            i64::try_from(bounds[2][0]).map_err(|_| "invalid z bounds")?,
            i64::try_from(bounds[2][1]).map_err(|_| "invalid z bounds")?,
        )?;
    }
    Ok(())
}

fn sqrt_circle_y(x_distance: i128) -> Result<i128, &'static str> {
    let radicand = RADIUS_UM * RADIUS_UM - x_distance * x_distance;
    let y = isqrt(radicand);
    if y * y != radicand {
        return Err("circle intersection is not an integer micrometre");
    }
    Ok(y)
}

fn append_central_half(
    blocks: &mut Vec<BracketBlock>,
    cx_um: i128,
    x_bounds: [i128; 2],
    y_sign: i128,
    id: &str,
    z_bounds: [i64; 2],
    bounded: bool,
) -> Result<(), &'static str> {
    let [x0, x1] = x_bounds;
    let [z0, z1] = z_bounds;
    let twice_mid = x0 + x1;
    if twice_mid > 2 * POCKET_LEFT_UM && twice_mid < 2 * POCKET_RIGHT_UM && z0 == 4 {
        return Ok(());
    }
    let sx = if twice_mid >= 2 * cx_um { 1 } else { -1 };
    let a = (x0 - cx_um).abs().min((x1 - cx_um).abs());
    let b = (x0 - cx_um).abs().max((x1 - cx_um).abs());
    let y0 = 0_i128;
    let y4 = y_sign * HALF_HEIGHT_UM;
    if a >= RADIUS_UM {
        return append_rect(
            blocks,
            cx_um,
            [[x0, x1], [y0, y4], [i128::from(z0), i128::from(z1)]],
            id,
        );
    }

    let xa = cx_um + sx * a;
    let circle_point = |x_um: i128| {
        let distance = (x_um - cx_um).abs();
        if bounded {
            let fm = y_sign * crate::travelling_bore_bounded::circle_height_fm(distance);
            Ok(Point {
                x_um,
                y_um: fm / 1_000_000_000,
                half_angle: None,
                bounded_y_fm: Some(fm),
            })
        } else {
            Point::circle(x_um, y_sign * sqrt_circle_y(distance)?, cx_um)
        }
    };
    let p = circle_point(xa)?;
    let t_l = Point::plain(xa, y4);
    if b <= RADIUS_UM {
        let xb = cx_um + sx * b;
        let e = circle_point(xb)?;
        let t_r = Point::plain(xb, y4);
        append_triangle(blocks, [p, e, t_r], cx_um, &format!("{id}_t0"), z0, z1)?;
        append_triangle(blocks, [p, t_r, t_l], cx_um, &format!("{id}_t1"), z0, z1)?;
    } else {
        let e = circle_point(cx_um + sx * RADIUS_UM)?;
        let b_point = Point::plain(cx_um + sx * b, 0);
        let t_r = Point::plain(cx_um + sx * b, y4);
        append_triangle(blocks, [p, e, t_r], cx_um, &format!("{id}_t0"), z0, z1)?;
        append_triangle(
            blocks,
            [e, b_point, t_r],
            cx_um,
            &format!("{id}_t1"),
            z0,
            z1,
        )?;
        append_triangle(blocks, [p, t_r, t_l], cx_um, &format!("{id}_t2"), z0, z1)?;
    }
    Ok(())
}

/// Build the two-layer plate with a bore centered at `(cx_um, 0)`.
///
/// The outer plate is `[-10,10]² × [0,8]` mm. A radius-2 mm circular bore
/// crosses the lower layer; the upper layer also has the fixed blind pocket
/// `[-6,-2] × [-4,4]` mm. Every returned control is an exact integer multiple
/// of its block's positive homogeneous scale.
pub fn travelling_bore_blocks(cx_um: i32) -> Result<Vec<BracketBlock>, &'static str> {
    build_travel(cx_um, false)
}

pub(crate) fn build_travel(cx_um: i32, bounded: bool) -> Result<Vec<BracketBlock>, &'static str> {
    let cx_um = i128::from(cx_um);
    if !(-5_000..=14_000).contains(&cx_um) {
        return Err("bore center outside supported travel interval");
    }
    for wall in [POCKET_LEFT_UM, POCKET_RIGHT_UM, HALF_WIDTH_UM] {
        if (wall - cx_um).abs() == RADIUS_UM {
            return Err("bore tangent to a vertical pocket or boundary wall");
        }
    }
    let xs = grid_x(cx_um);
    let mut blocks = Vec::new();
    for (column, pair) in xs.windows(2).enumerate() {
        let (x0, x1) = (pair[0], pair[1]);
        for (row, (y0, y1)) in [(-10_000_i128, -4_000_i128), (4_000, 10_000)]
            .into_iter()
            .enumerate()
        {
            for layer in 0..2 {
                let (z0, z1) = if layer == 0 { (0, 4) } else { (4, 8) };
                let id = format!("travel_c{column:02}_r{row}_z{layer}");
                append_rect(
                    &mut blocks,
                    cx_um,
                    [[x0, x1], [y0, y1], [i128::from(z0), i128::from(z1)]],
                    &id,
                )?;
            }
        }
        for y_sign in [1_i128, -1] {
            for layer in 0..2 {
                let (z0, z1) = if layer == 0 { (0, 4) } else { (4, 8) };
                let id = format!(
                    "travel_c{column:02}_ym{}_z{layer}",
                    if y_sign > 0 { "p" } else { "n" }
                );
                append_central_half(&mut blocks, cx_um, [x0, x1], y_sign, &id, [z0, z1], bounded)?;
            }
        }
    }
    if blocks.is_empty()
        || blocks.iter().any(|b| {
            b.degrees != [2, 2, 1] || b.controls.len() != 18 || b.controls.iter().any(|p| p[3] <= 0)
        })
    {
        return Err("invalid travelling bore block net");
    }
    Ok(blocks)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn radius_intersections_and_tangent_centers_are_classified_exactly() {
        assert_eq!(sqrt_circle_y(0).unwrap(), 2_000);
        assert_eq!(sqrt_circle_y(2_000).unwrap(), 0);
        assert!(sqrt_circle_y(1_000).is_err());
        for cx in [-8_000, -4_000, 8_000] {
            assert!(travelling_bore_blocks(cx).is_err());
        }
        assert!(travelling_bore_blocks(-5_001).is_err());
        assert!(travelling_bore_blocks(14_001).is_err());
    }

    #[test]
    fn supported_center_builds_positive_exact_blocks() {
        let blocks = travelling_bore_blocks(3_000).unwrap();
        assert!(!blocks.is_empty());
        assert!(blocks
            .iter()
            .all(|b| b.degrees == [2, 2, 1] && b.controls.len() == 18));
        let mut ids = blocks.iter().map(|b| b.id.as_str()).collect::<Vec<_>>();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), blocks.len());
    }

    #[test]
    fn exact_travel_samples_are_certified_and_tangencies_are_refused() {
        for cx_mm in -5..=14 {
            let cx_um = cx_mm * 1_000;
            let tangent = [-8_000, -4_000, 0, 8_000, 12_000].contains(&cx_um);
            let result = travelling_bore_blocks(cx_um);
            if tangent {
                assert!(
                    result.is_err(),
                    "expected exact tangency refusal at {cx_mm} mm"
                );
                continue;
            }
            match result {
                Ok(blocks) => {
                    assert!(!blocks.is_empty());
                    let incidence = crate::incidence::enumerate_faces(&blocks)
                        .unwrap_or_else(|e| panic!("incidence at {cx_mm} mm: {e}"));
                    assert!(!incidence.interfaces.is_empty());
                    assert!(blocks.iter().all(|b| b.controls.iter().all(|h| h[3] > 0)));
                }
                Err(e) => {
                    // Some center positions put a circle endpoint on the fixed
                    // plate edge at a nonintegral micrometre ordinate.
                    assert_eq!(e, "circle intersection is not an integer micrometre");
                }
            }
        }
        for cx_mm in [10_000, 11_920, 13_000] {
            let blocks = travelling_bore_blocks(cx_mm).unwrap();
            assert!(!blocks.is_empty());
        }
    }
}
