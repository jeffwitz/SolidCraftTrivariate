//! Geometric continuity at an oriented planar rational Bézier junction.
//! The left span ends at the join and the right span starts there.

use crate::constraint_system::Ratio;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Continuity {
    Disconnected,
    Singular,
    Corner,
    Tangent,
    Curvature,
}

type Pair = [Ratio; 2];
type Error = &'static str;

struct Jet {
    position: Pair,
    velocity: Pair,
    acceleration: Pair,
}

fn quotient(a: Ratio, b: Ratio) -> Result<Ratio, Error> {
    a.checked_mul(Ratio::new(b.denominator(), b.numerator())?)
}

fn endpoint(controls: &[[i64; 3]], end: bool) -> Result<Jet, Error> {
    if controls.len() < 2 || controls.iter().any(|p| p[2] <= 0) {
        return Err("a junction requires positive-weight spans of degree at least one");
    }
    let n = controls.len() - 1;
    let degree = Ratio::new(n as i128, 1)?;
    let second_factor = degree.checked_mul(Ratio::new((n - 1) as i128, 1)?)?;
    let point = |i: usize| controls[if end { n - i } else { i }].map(i128::from);
    let h = point(0);
    let near = point(1);
    let far = point(n.min(2));
    let inv_w = Ratio::new(1, h[2])?;
    // Differences first point inward at both ends. Reversing the parameter
    // changes the first derivative's sign and preserves the second derivative.
    let first = |d: usize| Ratio::new(near[d] - h[d], 1)?.checked_mul(degree);
    let second = |d: usize| Ratio::new(far[d] - 2 * near[d] + h[d], 1)?.checked_mul(second_factor);
    let dw = first(2)?;
    let ddw = second(2)?;
    let mut jet = Jet {
        position: [Ratio::zero(); 2],
        velocity: [Ratio::zero(); 2],
        acceleration: [Ratio::zero(); 2],
    };
    for (d, x) in jet.position.iter_mut().enumerate() {
        *x = Ratio::new(h[d], h[2])?;
        let v = first(d)?
            .checked_sub(x.checked_mul(dw)?)?
            .checked_mul(inv_w)?;
        jet.acceleration[d] = second(d)?
            .checked_sub(x.checked_mul(ddw)?)?
            .checked_sub(v.checked_mul(dw)?.checked_mul(Ratio::new(2, 1)?)?)?
            .checked_mul(inv_w)?;
        jet.velocity[d] = if end {
            Ratio::zero().checked_sub(v)?
        } else {
            v
        };
    }
    Ok(jet)
}

fn cross(a: Pair, b: Pair) -> Result<Ratio, Error> {
    a[0].checked_mul(b[1])?.checked_sub(a[1].checked_mul(b[0])?)
}

/// Classify using exact endpoint jets of the stored homogeneous (wx, wy, w)
/// controls. Curvature means at least G2. Singular means a stationary endpoint
/// prevents this regular-jet test from deciding geometric continuity.
pub fn classify(left: &[[i64; 3]], right: &[[i64; 3]]) -> Result<Continuity, Error> {
    let a = endpoint(left, true)?;
    let b = endpoint(right, false)?;
    if a.position != b.position {
        return Ok(Continuity::Disconnected);
    }
    let Some(axis) = a.velocity.iter().position(|v| v.numerator() != 0) else {
        return Ok(Continuity::Singular);
    };
    if b.velocity == [Ratio::zero(); 2] {
        return Ok(Continuity::Singular);
    }
    if cross(a.velocity, b.velocity)? != Ratio::zero()
        || a.velocity[axis].numerator().signum() != b.velocity[axis].numerator().signum()
    {
        return Ok(Continuity::Corner);
    }
    // If v_b = speed * v_a, G2 permits arbitrary tangential acceleration.
    // The normal part of a_b must equal that of speed^2 * a_a.
    let speed = quotient(b.velocity[axis], a.velocity[axis])?;
    let right_bend = cross(a.velocity, b.acceleration)?;
    let left_bend = cross(a.velocity, a.acceleration)?
        .checked_mul(speed)?
        .checked_mul(speed)?;
    Ok(if right_bend == left_bend {
        Continuity::Curvature
    } else {
        Continuity::Tangent
    })
}

/// A candidate partition event at the start of the named right-hand span.
/// At closure the incoming span ends at the upper limit of the active domain.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SketchJunction {
    pub right_span: usize,
    pub parameter: i64,
    pub multiplicity: usize,
    pub closure: bool,
    pub continuity: Continuity,
}

/// Inspect every junction, including closure, before building any volume.
/// All events are preserved; callers decide which ones require partition lines.
/// Disconnected and singular results are diagnostics, not admissions.
pub fn plan_loop(
    degree: usize,
    knots: &[i64],
    controls: &[[i64; 3]],
) -> Result<Vec<SketchJunction>, Error> {
    let spans = crate::nurbs_curve::extract_active_bezier(degree, knots, controls)?;
    let mut events = Vec::with_capacity(spans.len());
    for (index, right) in spans.iter().enumerate() {
        let left = &spans[if index == 0 {
            spans.len() - 1
        } else {
            index - 1
        }];
        let parameter = right.interval[0];
        let first = knots.partition_point(|k| *k < parameter);
        let after = knots.partition_point(|k| *k <= parameter);
        events.push(SketchJunction {
            right_span: index,
            parameter,
            multiplicity: after - first,
            closure: index == 0,
            continuity: classify(&left.controls, &right.controls)?,
        });
    }
    Ok(events)
}
