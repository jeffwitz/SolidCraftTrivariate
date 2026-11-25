//! Exact recognition of quarter-cylinder boundary patches.
//!
//! This module never fits a surface. It accepts only a quadratic rational
//! quarter-circle extruded along one Cartesian axis and returns its exact
//! STEP placement data. All predicates use checked rational arithmetic.

use crate::{boundary_shell::ShellFace, constraint_system::Ratio};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct CylinderSurface {
    /// Exact axis-line point in millimetres.
    pub center: [Ratio; 3],
    /// Signed Cartesian unit vector, oriented along the extrusion parameter.
    pub axis: [i8; 3],
    /// Exact positive radius in millimetres.
    pub radius: Ratio,
    /// Exact unit radial vector at the first angular pole.
    pub reference: [Ratio; 3],
}

type V3 = [Ratio; 3];

fn zero3() -> V3 {
    [Ratio::zero(); 3]
}

fn add(a: V3, b: V3) -> Result<V3, &'static str> {
    Ok([
        a[0].checked_add(b[0])?,
        a[1].checked_add(b[1])?,
        a[2].checked_add(b[2])?,
    ])
}

fn sub(a: V3, b: V3) -> Result<V3, &'static str> {
    Ok([
        a[0].checked_sub(b[0])?,
        a[1].checked_sub(b[1])?,
        a[2].checked_sub(b[2])?,
    ])
}

fn mul(a: V3, scale: Ratio) -> Result<V3, &'static str> {
    Ok([
        a[0].checked_mul(scale)?,
        a[1].checked_mul(scale)?,
        a[2].checked_mul(scale)?,
    ])
}

fn div(a: V3, scale: Ratio) -> Result<V3, &'static str> {
    if scale.numerator() == 0 {
        return Err("cylinder_division_by_zero");
    }
    let reciprocal = Ratio::new(scale.denominator(), scale.numerator())?;
    mul(a, reciprocal)
}

fn dot(a: V3, b: V3) -> Result<Ratio, &'static str> {
    a[0].checked_mul(b[0])?
        .checked_add(a[1].checked_mul(b[1])?)?
        .checked_add(a[2].checked_mul(b[2])?)
}

fn cross(a: V3, b: V3) -> Result<V3, &'static str> {
    Ok([
        a[1].checked_mul(b[2])?
            .checked_sub(a[2].checked_mul(b[1])?)?,
        a[2].checked_mul(b[0])?
            .checked_sub(a[0].checked_mul(b[2])?)?,
        a[0].checked_mul(b[1])?
            .checked_sub(a[1].checked_mul(b[0])?)?,
    ])
}

fn point(h: [i64; 4]) -> Result<V3, &'static str> {
    if h[3] <= 0 {
        return Err("invalid_boundary_weight");
    }
    let denominator = i128::from(h[3]);
    Ok([
        Ratio::new(i128::from(h[0]), denominator)?,
        Ratio::new(i128::from(h[1]), denominator)?,
        Ratio::new(i128::from(h[2]), denominator)?,
    ])
}

fn is_zero(value: Ratio) -> bool {
    value.numerator() == 0
}

fn integer_sqrt(n: u128) -> Option<u128> {
    let root = n.isqrt();
    (root * root == n).then_some(root)
}

fn rational_sqrt(value: Ratio) -> Result<Option<Ratio>, &'static str> {
    if value.numerator() <= 0 {
        return Ok(None);
    }
    let Some(n) = integer_sqrt(value.numerator() as u128) else {
        return Ok(None);
    };
    let Some(d) = integer_sqrt(value.denominator() as u128) else {
        return Ok(None);
    };
    let n = i128::try_from(n).map_err(|_| "cylinder_radius_exceeds_i128")?;
    let d = i128::try_from(d).map_err(|_| "cylinder_radius_exceeds_i128")?;
    Ok(Some(Ratio::new(n, d)?))
}

/// Recognize one exact quarter-cylinder patch from its homogeneous face net.
///
/// The supported net has degrees `(2, 1)`: the first parameter is the
/// circular arc, the second is a straight extrusion. Unsupported geometry,
/// an inward angular normal, or a non-rational radius returns `Ok(None)`;
/// malformed weights and arithmetic overflow return an error.
pub(crate) fn cylinder_surface(face: &ShellFace) -> Result<Option<CylinderSurface>, &'static str> {
    if face.side != "v1" || face.degrees != [2, 1] || face.controls.len() != 6 {
        return Ok(None);
    }
    if face.controls.iter().any(|control| control[3] <= 0) {
        return Err("invalid_boundary_weight");
    }

    // Face nets are stored with the second index varying fastest:
    // (u0,v0),(u0,v1),(u1,v0),... .
    let mut layer = [zero3(); 3];
    let mut translated_layer = [zero3(); 3];
    let mut weights = [[0_i64; 2]; 3];
    for i in 0..3 {
        let h0 = face.controls[2 * i];
        let h1 = face.controls[2 * i + 1];
        weights[i] = [h0[3], h1[3]];
        if h0[3] != h1[3] {
            return Ok(None);
        }
        layer[i] = point(h0)?;
        translated_layer[i] = point(h1)?;
    }
    let extrusion = sub(translated_layer[0], layer[0])?;
    if extrusion == zero3() {
        return Ok(None);
    }
    for i in 1..3 {
        if sub(translated_layer[i], layer[i])? != extrusion {
            return Ok(None);
        }
    }
    let nonzero_axes = (0..3)
        .filter(|&i| !is_zero(extrusion[i]))
        .collect::<Vec<_>>();
    if nonzero_axes.len() != 1 {
        return Ok(None);
    }
    let axis_index = nonzero_axes[0];
    let mut axis = [0_i8; 3];
    axis[axis_index] = if extrusion[axis_index].numerator() > 0 {
        1
    } else {
        -1
    };

    let center = sub(add(layer[0], layer[2])?, layer[1])?;
    let radial_start = sub(layer[0], center)?;
    let radial_middle = sub(layer[1], center)?;
    let radial_end = sub(layer[2], center)?;
    if !is_zero(radial_start[axis_index])
        || !is_zero(radial_middle[axis_index])
        || !is_zero(radial_end[axis_index])
    {
        return Ok(None);
    }
    let start_sq = dot(radial_start, radial_start)?;
    let end_sq = dot(radial_end, radial_end)?;
    if start_sq != end_sq
        || dot(radial_start, radial_end)? != Ratio::zero()
        || start_sq.numerator() <= 0
    {
        return Ok(None);
    }

    let w0 = i128::from(weights[0][0]);
    let w1 = i128::from(weights[1][0]);
    let w2 = i128::from(weights[2][0]);
    if w0.checked_mul(w2).ok_or("cylinder_weight_overflow")?
        != 2_i128
            .checked_mul(w1.checked_mul(w1).ok_or("cylinder_weight_overflow")?)
            .ok_or("cylinder_weight_overflow")?
    {
        return Ok(None);
    }
    // Require the exact rational conic control to be the tangent intersection
    // already implied by center = P0 + P2 - P1; this also rejects malformed
    // degenerate nets before constructing an axis placement.
    if radial_middle == zero3() {
        return Ok(None);
    }
    // The natural parameter normal is du x dv. Keep only outward cylinder
    // normals so the caller can retain the source face's STEP sense unchanged.
    let angular_tangent = sub(layer[1], layer[0])?;
    let orientation = dot(cross(angular_tangent, extrusion)?, radial_start)?;
    if orientation.numerator() <= 0 {
        return Ok(None);
    }
    let Some(radius) = rational_sqrt(start_sq)? else {
        return Ok(None);
    };
    let reference = div(radial_start, radius)?;
    Ok(Some(CylinderSurface {
        center,
        axis,
        radius,
        reference,
    }))
}
