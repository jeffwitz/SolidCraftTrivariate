//! Exact charts for one bounded pair of intersecting, oblique blind pockets.
//!
//! The fixture is deliberately narrow: a ten-millimetre cube carries two
//! interior footprints on its top face. The first has a quadratic wall with a
//! C1 span join; the second is a linear rectangle with G0 corners. Their
//! separate sweep vectors are `(0, 1/5, -1)` and `(0, -1/10, -1)` mm per mm
//! of depth. This is not a common affine shear of a pre-existing Boolean.
//!
//! Crossed wall pairs are partitioned by their exact rational intersection
//! graph. The two side wings sweep the complete section arrangement, including
//! cells occupied by the cutters, so the pocket side walls have whole-face
//! traces. The source-loop corners remain G0; the quadratic profile also has
//! an internal C1 join at x=5.

use crate::{bracket::BracketBlock, constraint_system::Ratio as Q};

type Error = &'static str;

const DEGREES: [usize; 3] = [2, 1, 1];

fn q(n: i128, d: i128) -> Result<Q, Error> {
    Q::new(n, d)
}
fn add(a: Q, b: Q) -> Result<Q, Error> {
    a.checked_add(b)
}
fn sub(a: Q, b: Q) -> Result<Q, Error> {
    a.checked_sub(b)
}
fn mul(a: Q, b: Q) -> Result<Q, Error> {
    a.checked_mul(b)
}
fn div(a: Q, b: Q) -> Result<Q, Error> {
    a.checked_div(b)
}
fn int(n: i64) -> Q {
    Q::new(i128::from(n), 1).expect("integer ratio")
}

fn profile(x: Q) -> Result<Q, Error> {
    let d = sub(x, int(5))?;
    let divisor = if x.cmp_value(&int(5)).is_gt() {
        500
    } else {
        1000
    };
    add(int(2), div(mul(d, d)?, int(divisor))?)
}

fn crossing_depth(x: Q) -> Result<Q, Error> {
    div(sub(int(3), profile(x)?)?, q(3, 10)?)
}

fn a_lower(x: Q, t: Q) -> Result<Q, Error> {
    add(profile(x)?, div(t, int(5))?)
}
fn a_upper(x: Q, t: Q) -> Result<Q, Error> {
    add(a_lower(x, t)?, int(2))
}
fn b_lower(t: Q) -> Result<Q, Error> {
    sub(int(3), div(t, int(10))?)
}
fn b_upper(t: Q) -> Result<Q, Error> {
    add(b_lower(t)?, int(2))
}

#[derive(Clone, Copy)]
enum SweepRegion {
    /// Both pockets are active above the crossing graph.
    BeforeCross,
    /// Both pockets are active below the crossing graph.
    AfterCross,
    /// The first pocket remains after the second pocket reaches its floor.
    FirstPocketFloor,
    /// Material below the first pocket floor; retain its separator grid.
    FullMaterial,
}

#[derive(Clone, Copy)]
enum TBound {
    Zero,
    Crossing,
    Five,
    Seven,
    Ten,
}

#[derive(Clone, Copy)]
enum Boundary {
    Zero,
    Ten,
    EnvelopeLow,
    EnvelopeHigh,
    MiddleOne,
    MiddleTwo,
    BLow,
    ALow,
    BHigh,
    AHigh,
}

#[derive(Clone, Copy)]
struct ChartSpec {
    physical_x: [i64; 2],
    profile_x: [i64; 2],
    region: SweepRegion,
    t_range: [TBound; 2],
    low: Boundary,
    high: Boundary,
}

fn envelope(region: SweepRegion, x: Q, t: Q) -> Result<(Q, Q), Error> {
    match region {
        SweepRegion::BeforeCross => Ok((a_lower(x, t)?, b_upper(t)?)),
        SweepRegion::AfterCross => Ok((b_lower(t)?, a_upper(x, t)?)),
        _ => Err("oblique_pocket_region_has_no_crossed_envelope"),
    }
}

fn wall_t(region: SweepRegion, t: Q) -> Q {
    match region {
        SweepRegion::FullMaterial => int(7),
        _ => t,
    }
}

fn boundary_value(boundary: Boundary, region: SweepRegion, x: Q, t: Q) -> Result<Q, Error> {
    let tw = wall_t(region, t);
    match boundary {
        Boundary::Zero => Ok(Q::zero()),
        Boundary::Ten => Ok(int(10)),
        Boundary::EnvelopeLow => Ok(envelope(region, x, tw)?.0),
        Boundary::EnvelopeHigh => Ok(envelope(region, x, tw)?.1),
        Boundary::BLow => b_lower(tw),
        Boundary::ALow => a_lower(x, tw),
        Boundary::BHigh => b_upper(tw),
        Boundary::AHigh => a_upper(x, tw),
        Boundary::MiddleOne | Boundary::MiddleTwo => {
            // These polynomial separators reach the first pocket's lower wall
            // and the second pocket's upper wall exactly at the second floor.
            // They remain strictly inside the envelope at the wall crossing.
            let delta_at_floor = sub(profile(x)?, q(3, 2)?)?;
            let fraction_numerator = match boundary {
                Boundary::MiddleOne => {
                    add(int(1), div(mul(tw, sub(delta_at_floor, int(1))?)?, int(5))?)?
                }
                Boundary::MiddleTwo => add(q(3, 2)?, div(tw, int(10))?)?,
                _ => unreachable!(),
            };
            add(b_lower(tw)?, fraction_numerator)
        }
    }
}

fn sample_point(spec: ChartSpec, u: Q, w: Q, v: Q) -> Result<[Q; 3], Error> {
    let px = add(
        int(spec.physical_x[0]),
        mul(int(spec.physical_x[1] - spec.physical_x[0]), u)?,
    )?;
    let sx = if spec.profile_x[0] == spec.profile_x[1] {
        int(spec.profile_x[0])
    } else {
        add(
            int(spec.profile_x[0]),
            mul(int(spec.profile_x[1] - spec.profile_x[0]), u)?,
        )?
    };
    let t_bound = |bound: TBound| match bound {
        TBound::Zero => Ok(Q::zero()),
        TBound::Crossing => crossing_depth(sx),
        TBound::Five => Ok(int(5)),
        TBound::Seven => Ok(int(7)),
        TBound::Ten => Ok(int(10)),
    };
    let t0 = t_bound(spec.t_range[0])?;
    let t1 = t_bound(spec.t_range[1])?;
    let t = add(t0, mul(sub(t1, t0)?, w)?)?;
    let low = boundary_value(spec.low, spec.region, sx, t)?;
    let high = boundary_value(spec.high, spec.region, sx, t)?;
    let y = add(low, mul(v, sub(high, low)?)?)?;
    Ok([px, y, sub(int(10), t)?])
}

fn bernstein_controls(values: &[Q; 3]) -> Result<[Q; 3], Error> {
    // Convert the values at 0, 1/2 and 1 to quadratic Bernstein controls.
    let p0 = values[0];
    let p2 = values[2];
    let p1 = sub(mul(int(2), values[1])?, div(add(p0, p2)?, int(2))?)?;
    Ok([p0, p1, p2])
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.max(1)
}

fn common_scale(values: &[Q]) -> Result<i128, Error> {
    let mut scale = 1i128;
    for value in values {
        let den = value.denominator();
        let common = i128::try_from(gcd(scale as u128, den as u128))
            .map_err(|_| "oblique_pocket_scale_exceeds_i128")?;
        scale = (scale / common)
            .checked_mul(den)
            .ok_or("oblique_pocket_scale_exceeds_i128")?;
    }
    Ok(scale)
}

fn block(id: String, spec: ChartSpec) -> Result<BracketBlock, Error> {
    let u_samples = [Q::zero(), q(1, 2)?, Q::one()];
    let mut controls = vec![[Q::zero(); 4]; 12];
    for j in 0..2 {
        for k in 0..2 {
            let w = int(j as i64);
            let v = int(k as i64);
            let samples = u_samples
                .map(|u| sample_point(spec, u, w, v))
                .into_iter()
                .collect::<Result<Vec<_>, _>>()?;
            for axis in 0..3 {
                let converted =
                    bernstein_controls(&[samples[0][axis], samples[1][axis], samples[2][axis]])?;
                for (i, value) in converted.into_iter().enumerate() {
                    let index = (i * 2 + j) * 2 + k;
                    controls[index][axis] = value;
                }
            }
        }
    }
    for p in &mut controls {
        p[3] = Q::one();
    }
    let scale = common_scale(
        &controls
            .iter()
            .flat_map(|p| p.iter().copied())
            .collect::<Vec<_>>(),
    )?;
    let mut integer_controls = Vec::with_capacity(12);
    let mut divisor = 0u128;
    let mut scaled = Vec::with_capacity(12);
    for p in controls {
        let mut h = [0i128; 4];
        for axis in 0..4 {
            h[axis] = p[axis]
                .numerator()
                .checked_mul(scale / p[axis].denominator())
                .ok_or("oblique_pocket_control_overflow")?;
            divisor = gcd(divisor, h[axis].unsigned_abs());
        }
        scaled.push(h);
    }
    let divisor = i128::try_from(divisor).map_err(|_| "oblique_pocket_scale_exceeds_i128")?;
    for h in scaled {
        integer_controls.push([
            i64::try_from(h[0] / divisor).map_err(|_| "oblique_pocket_control_overflow")?,
            i64::try_from(h[1] / divisor).map_err(|_| "oblique_pocket_control_overflow")?,
            i64::try_from(h[2] / divisor).map_err(|_| "oblique_pocket_control_overflow")?,
            i64::try_from(h[3] / divisor).map_err(|_| "oblique_pocket_weight_overflow")?,
        ]);
    }
    Ok(BracketBlock {
        id,
        degrees: DEGREES,
        controls: integer_controls,
    })
}

fn add_band(out: &mut Vec<BracketBlock>, id: &str, spec: ChartSpec) -> Result<(), Error> {
    out.push(block(id.to_owned(), spec)?);
    Ok(())
}

/// Build the fixed `oblique_blind_pockets_v1` fixture as 64 whole rational
/// Bézier charts. Parameters and exact source profiles are recorded in the
/// companion reference recipe.
pub fn build_oblique_blind_pockets() -> Result<Vec<BracketBlock>, Error> {
    let mut out = Vec::new();
    let zones = [
        ([0, 1], [1, 1], true),
        ([1, 5], [1, 5], false),
        ([5, 9], [5, 9], false),
        ([9, 10], [9, 9], true),
    ];
    for (zone_id, (physical_x, profile_x, wing)) in zones.into_iter().enumerate() {
        let prefix = ["left_wing", "interior_left", "interior_right", "right_wing"][zone_id];
        for (region, trange) in [
            (SweepRegion::BeforeCross, [TBound::Zero, TBound::Crossing]),
            (SweepRegion::AfterCross, [TBound::Crossing, TBound::Five]),
        ] {
            let lower = Boundary::EnvelopeLow;
            let upper = Boundary::EnvelopeHigh;
            let bands = if wing {
                vec![
                    (Boundary::Zero, lower),
                    (lower, Boundary::MiddleOne),
                    (Boundary::MiddleOne, Boundary::MiddleTwo),
                    (Boundary::MiddleTwo, upper),
                    (upper, Boundary::Ten),
                ]
            } else {
                vec![(Boundary::Zero, lower), (upper, Boundary::Ten)]
            };
            let region_name = if matches!(region, SweepRegion::BeforeCross) {
                "before_cross"
            } else {
                "after_cross"
            };
            for (index, (lo, hi)) in bands.into_iter().enumerate() {
                add_band(
                    &mut out,
                    &format!("{prefix}_{region_name}_band{index}"),
                    ChartSpec {
                        physical_x,
                        profile_x,
                        region,
                        t_range: trange,
                        low: lo,
                        high: hi,
                    },
                )?;
            }
        }

        let floor_bands = [
            (Boundary::Zero, Boundary::BLow),
            (Boundary::BLow, Boundary::ALow),
            (Boundary::ALow, Boundary::BHigh),
            (Boundary::BHigh, Boundary::AHigh),
            (Boundary::AHigh, Boundary::Ten),
        ];
        for (index, (lo, hi)) in floor_bands.into_iter().enumerate() {
            if !wing && (index == 2 || index == 3) {
                continue; // Both intervals are still inside the first cutter.
            }
            add_band(
                &mut out,
                &format!("{prefix}_first_floor_band{index}"),
                ChartSpec {
                    physical_x,
                    profile_x,
                    region: SweepRegion::FirstPocketFloor,
                    t_range: [TBound::Five, TBound::Seven],
                    low: lo,
                    high: hi,
                },
            )?;
        }

        for (index, (lo, hi)) in floor_bands.into_iter().enumerate() {
            add_band(
                &mut out,
                &format!("{prefix}_full_band{index}"),
                ChartSpec {
                    physical_x,
                    profile_x,
                    region: SweepRegion::FullMaterial,
                    t_range: [TBound::Seven, TBound::Ten],
                    low: lo,
                    high: hi,
                },
            )?;
        }
    }
    Ok(out)
}

/// Rebuild and compare the complete fixed family. This checks coefficients,
/// whole-face conformity and closed-domain Jacobian certificates on replay.
pub fn audit_oblique_blind_pockets(blocks: &[BracketBlock]) -> Result<(), Error> {
    let expected = build_oblique_blind_pockets()?;
    if blocks != expected {
        return Err("oblique_blind_pocket_coefficients_mismatch");
    }
    crate::incidence::enumerate_faces(blocks)?;
    for block in blocks {
        let certificate = crate::spatial_certificate::certify_spatial(
            block.degrees,
            &block.controls,
            crate::spatial_certificate::CertificateBudget::default(),
        );
        if certificate.status != crate::spatial_certificate::CertificateStatus::Valid {
            return Err("oblique_blind_pocket_jacobian_not_certified");
        }
    }
    Ok(())
}
