//! Exact classification and analytic sensitivities for two planar circles.
//!
//! Inputs are integer micrometres `[center_x, center_y, radius]`. Exact i128
//! predicates classify the topology before floating-point coordinates are
//! formed. Returned points are in millimetres and point jets are derivatives
//! in millimetres per millimetre, with parameter order
//! `[c1x,c1y,r1,c2x,c2y,r2]`.

const MAX_ABS_CENTER_UM: i64 = 1_000_000; // 1 m
const MAX_RADIUS_UM: i64 = 1_000_000; // 1 m

/// Topological result for the intersection of two circles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SketchRelation {
    /// Two isolated transverse solutions, with a differentiable local branch.
    Transverse,
    /// One coincident solution at external or internal tangency.
    TangentExternal,
    TangentInternal,
    /// Equal circles: the intersection is a one-dimensional set.
    Coincident,
    /// Concentric circles of different radii have no intersection.
    ConcentricDisjoint,
    /// The circles are separated by a positive gap.
    Disjoint,
    /// One circle is strictly contained inside the other.
    StrictContainment,
}

/// Exact cardinality class, distinguishing a continuum from zero points.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntersectionCardinality {
    Zero,
    One,
    Two,
    Infinite,
}

/// One transverse intersection, ordered relative to the directed center axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SketchPoint {
    /// Physical point in mm.
    pub point_mm: [f64; 2],
    /// `d(point_mm)/d([c1x,c1y,r1,c2x,c2y,r2])`, all in mm/mm.
    pub jet_mm_per_mm: [[f64; 2]; 6],
    /// Condition estimate for solving the two circle constraints.
    pub condition_estimate: f64,
    /// Maximum absolute circle-constraint residual in mm².
    pub constraint_residual_mm2: f64,
}

/// Exact relation and (when transverse) the two ordered points.
#[derive(Debug, Clone, PartialEq)]
pub struct SketchIntersection {
    pub relation: SketchRelation,
    pub cardinality: IntersectionCardinality,
    /// Finite count of isolated solutions; zero for the coincident continuum.
    pub intersection_count: usize,
    /// Jacobian rank of the two circle equations at an isolated solution.
    /// `None` means the solution is not an isolated point.
    pub rank: Option<usize>,
    pub points: Vec<SketchPoint>,
    /// Unique geometric point at a tangency; its derivative is unavailable.
    pub tangent_point_mm: Option<[f64; 2]>,
    /// Explanation when an isolated differentiable solution is unavailable.
    pub derivative_unavailable: Option<&'static str>,
}

/// Intersect circles `[cx_um, cy_um, radius_um]` using exact integer topology.
///
/// The squared-height discriminant is evaluated as
/// `((r1+r2)^2-D)*(D-(r1-r2)^2)` in checked i128 arithmetic. This avoids the
/// cancellation in `4*r1^2*D-(r1^2-r2^2+D)^2` near tangencies. Floating point
/// is used only after exact classification to return coordinates and jets.
pub fn intersect_circles(
    first: [i64; 3],
    second: [i64; 3],
) -> Result<SketchIntersection, &'static str> {
    validate_circle(first)?;
    validate_circle(second)?;
    let dx_i = second[0] as i128 - first[0] as i128;
    let dy_i = second[1] as i128 - first[1] as i128;
    let r1_i = first[2] as i128;
    let r2_i = second[2] as i128;
    let d2 = dx_i
        .checked_mul(dx_i)
        .and_then(|v| dy_i.checked_mul(dy_i).and_then(|w| v.checked_add(w)))
        .ok_or("integer distance overflow")?;
    let sum = r1_i.checked_add(r2_i).ok_or("integer radius overflow")?;
    let diff = r1_i.checked_sub(r2_i).ok_or("integer radius overflow")?;
    let outer = sum
        .checked_mul(sum)
        .and_then(|v| v.checked_sub(d2))
        .ok_or("integer discriminant overflow")?;
    let inner = d2
        .checked_sub(
            diff.checked_mul(diff)
                .ok_or("integer discriminant overflow")?,
        )
        .ok_or("integer discriminant overflow")?;

    let (relation, cardinality, count) = if d2 == 0 {
        if r1_i == r2_i {
            (
                SketchRelation::Coincident,
                IntersectionCardinality::Infinite,
                0,
            )
        } else {
            (
                SketchRelation::ConcentricDisjoint,
                IntersectionCardinality::Zero,
                0,
            )
        }
    } else if outer < 0 {
        (SketchRelation::Disjoint, IntersectionCardinality::Zero, 0)
    } else if inner < 0 {
        (
            SketchRelation::StrictContainment,
            IntersectionCardinality::Zero,
            0,
        )
    } else if outer == 0 {
        (
            SketchRelation::TangentExternal,
            IntersectionCardinality::One,
            1,
        )
    } else if inner == 0 {
        (
            SketchRelation::TangentInternal,
            IntersectionCardinality::One,
            1,
        )
    } else {
        (SketchRelation::Transverse, IntersectionCardinality::Two, 2)
    };

    if relation != SketchRelation::Transverse {
        let reason = match relation {
            SketchRelation::TangentExternal | SketchRelation::TangentInternal => {
                "tangent constraints have rank one; branch derivative is singular"
            }
            SketchRelation::Coincident => "coincident circles form a continuum of solutions",
            SketchRelation::ConcentricDisjoint => "concentric circles have no isolated solution",
            SketchRelation::Disjoint => "circles have no real intersection",
            SketchRelation::StrictContainment => "strictly contained circles have no intersection",
            SketchRelation::Transverse => unreachable!(),
        };
        return Ok(SketchIntersection {
            relation,
            cardinality,
            intersection_count: count,
            rank: if count == 1 { Some(1) } else { None },
            points: Vec::new(),
            tangent_point_mm: if count == 1 {
                let dx = dx_i as f64 / 1000.0;
                let dy = dy_i as f64 / 1000.0;
                let r1 = r1_i as f64 / 1000.0;
                let r2 = r2_i as f64 / 1000.0;
                let d2 = dx.mul_add(dx, dy * dy);
                let alpha = 0.5 + (r1 * r1 - r2 * r2) / (2.0 * d2);
                Some([
                    first[0] as f64 / 1000.0 + alpha * dx,
                    first[1] as f64 / 1000.0 + alpha * dy,
                ])
            } else {
                None
            },
            derivative_unavailable: Some(reason),
        });
    }

    // Convert to mm only after exact classification. The bounded inputs keep
    // these conversions and all subsequent products well inside f64 range.
    let dx = dx_i as f64 / 1000.0;
    let dy = dy_i as f64 / 1000.0;
    let r1 = r1_i as f64 / 1000.0;
    let r2 = r2_i as f64 / 1000.0;
    let d2f = dx.mul_add(dx, dy * dy);
    let outer_exact = outer as f64 / 1_000_000.0;
    let inner_exact = inner as f64 / 1_000_000.0;
    let alpha = 0.5 + (r1 * r1 - r2 * r2) / (2.0 * d2f);
    let height = (outer_exact * inner_exact).sqrt() / (2.0 * d2f.sqrt());
    if !height.is_finite() || height <= 0.0 {
        return Err("transverse height is outside the resolvable floating-point range");
    }

    let beta = height / d2f.sqrt();
    let mut points = Vec::with_capacity(2);
    // Positive sign is the left-hand side of the directed c1 -> c2 axis.
    for sign in [1.0, -1.0] {
        let x = first[0] as f64 / 1000.0 + alpha * dx - sign * beta * dy;
        let y = first[1] as f64 / 1000.0 + alpha * dy + sign * beta * dx;
        let mut jet = [[0.0; 2]; 6];
        for (p, tangent) in jet.iter_mut().enumerate() {
            let (ddx, ddy, dr1, dr2, dc1x, dc1y) = match p {
                0 => (-1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
                1 => (0.0, -1.0, 0.0, 0.0, 0.0, 1.0),
                2 => (0.0, 0.0, 1.0, 0.0, 0.0, 0.0),
                3 => (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                4 => (0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
                _ => (0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
            };
            let d_d2 = 2.0 * (dx * ddx + dy * ddy);
            let d_alpha =
                (r1 * dr1 - r2 * dr2) / d2f - (r1 * r1 - r2 * r2) * d_d2 / (2.0 * d2f * d2f);
            // Factor derivatives stay well-scaled close to a tangency: the
            // exact positive factors are used before conversion to f64.
            let sum_r = r1 + r2;
            let diff_r = r1 - r2;
            let l = outer_exact;
            let u = inner_exact;
            let d_l = 2.0 * sum_r * (dr1 + dr2) - d_d2;
            let d_u = d_d2 - 2.0 * diff_r * (dr1 - dr2);
            let d_beta = beta * 0.5 * (d_l / l + d_u / u) - beta * d_d2 / d2f;
            tangent[0] = dc1x + d_alpha * dx + alpha * ddx - sign * (d_beta * dy + beta * ddy);
            tangent[1] = dc1y + d_alpha * dy + alpha * ddy + sign * (d_beta * dx + beta * ddx);
        }
        let fdx = x - first[0] as f64 / 1000.0;
        let fdy = y - first[1] as f64 / 1000.0;
        let sdx = x - second[0] as f64 / 1000.0;
        let sdy = y - second[1] as f64 / 1000.0;
        let residual = (fdx.mul_add(fdx, fdy * fdy) - r1 * r1)
            .abs()
            .max((sdx.mul_add(sdx, sdy * sdy) - r2 * r2).abs());
        // Singular values of the two independent circle-constraint rows.
        // The product is d*h for the rows p-c_i, so the small singular
        // value is never obtained by subtracting nearly equal eigenvalues.
        let gram_trace = r1 * r1 + r2 * r2;
        let twice_dot = gram_trace - d2f;
        let gram_gap = (r1 * r1 - r2 * r2).hypot(twice_dot);
        let condition = (gram_trace + gram_gap) / (2. * d2f.sqrt() * height);
        if jet.iter().flatten().any(|v| !v.is_finite()) || !condition.is_finite() {
            return Err("intersection jet exceeds finite floating-point range");
        }
        points.push(SketchPoint {
            point_mm: [x, y],
            jet_mm_per_mm: jet,
            condition_estimate: condition,
            constraint_residual_mm2: residual,
        });
    }
    Ok(SketchIntersection {
        relation,
        cardinality,
        intersection_count: count,
        rank: Some(2),
        points,
        tangent_point_mm: None,
        derivative_unavailable: None,
    })
}

fn validate_circle(circle: [i64; 3]) -> Result<(), &'static str> {
    if circle[0].unsigned_abs() > MAX_ABS_CENTER_UM as u64
        || circle[1].unsigned_abs() > MAX_ABS_CENTER_UM as u64
        || circle[2] <= 0
        || circle[2] > MAX_RADIUS_UM
    {
        return Err("circle coordinates or radius are outside the supported domain");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{a} != {b} within {tol}");
    }

    #[test]
    fn almost_coincident_circles_have_ill_conditioned_constraint_rows() {
        let result = intersect_circles([0, 0, 1000000], [1, 0, 1000000]).unwrap();
        for p in result.points {
            assert!((p.condition_estimate / 2_000_000. - 1.).abs() < 1e-10);
        }
    }

    #[test]
    fn classifies_exact_transverse_tangent_and_empty_cases() {
        assert_eq!(
            intersect_circles([0, 0, 5], [6, 0, 5])
                .unwrap()
                .intersection_count,
            2
        );
        assert_eq!(
            intersect_circles([0, 0, 5], [10, 0, 5]).unwrap().relation,
            SketchRelation::TangentExternal
        );
        let tangent = intersect_circles([0, 0, 5], [10, 0, 5]).unwrap();
        assert_eq!(tangent.cardinality, IntersectionCardinality::One);
        assert_eq!(tangent.tangent_point_mm, Some([0.005, 0.0]));
        assert!(tangent.derivative_unavailable.is_some());
        let internal = intersect_circles([0, 0, 5], [4, 0, 1]).unwrap();
        assert_eq!(internal.tangent_point_mm, Some([0.005, 0.0]));
        assert_eq!(
            intersect_circles([0, 0, 5], [4, 0, 1]).unwrap().relation,
            SketchRelation::TangentInternal
        );
        assert_eq!(
            intersect_circles([0, 0, 5], [11, 0, 5]).unwrap().relation,
            SketchRelation::Disjoint
        );
        assert_eq!(
            intersect_circles([0, 0, 5], [2, 0, 1]).unwrap().relation,
            SketchRelation::StrictContainment
        );
        let coincident = intersect_circles([0, 0, 5], [0, 0, 5]).unwrap();
        assert_eq!(coincident.relation, SketchRelation::Coincident);
        assert_eq!(coincident.cardinality, IntersectionCardinality::Infinite);
        assert_eq!(
            intersect_circles([0, 0, 5], [0, 0, 4]).unwrap().relation,
            SketchRelation::ConcentricDisjoint
        );
        assert!(intersect_circles([0, 0, 0], [1, 0, 1]).is_err());
        assert!(intersect_circles([1_000_001, 0, 1], [0, 0, 1]).is_err());
    }

    #[test]
    fn orders_points_and_jets_covariantly_under_translation_and_rotation() {
        let base = intersect_circles([0, 0, 5_000], [6_000, 0, 5_000]).unwrap();
        near(base.points[0].point_mm[0], 3.0, 1e-14);
        near(base.points[0].point_mm[1], 4.0, 1e-14);
        near(base.points[1].point_mm[1], -4.0, 1e-14);
        let translated = intersect_circles([7000, -3000, 5000], [13000, -3000, 5000]).unwrap();
        near(
            translated.points[0].point_mm[0],
            base.points[0].point_mm[0] + 7.0,
            1e-14,
        );
        near(
            translated.points[0].point_mm[1],
            base.points[0].point_mm[1] - 3.0,
            1e-14,
        );
        // Rotate the center axis by +90 degrees; the left-normal ordering rotates too.
        let rotated = intersect_circles([0, 0, 5000], [0, 6000, 5000]).unwrap();
        near(rotated.points[0].point_mm[0], -4.0, 1e-14);
        near(rotated.points[0].point_mm[1], 3.0, 1e-14);
        let swapped = intersect_circles([6_000, 0, 5_000], [0, 0, 5_000]).unwrap();
        near(
            swapped.points[0].point_mm[1],
            base.points[1].point_mm[1],
            1e-14,
        );
        near(
            swapped.points[1].point_mm[1],
            base.points[0].point_mm[1],
            1e-14,
        );
    }

    #[test]
    fn analytic_jets_match_centered_differences_and_taylor_remainder() {
        let circles = [[0_i64, 0, 5000], [6000, 1000, 5000]];
        let result = intersect_circles(circles[0], circles[1]).unwrap();
        let point = result.points[0];
        let h_um = 1_i64;
        for p in 0..6 {
            let mut plus = circles;
            let mut minus = circles;
            let circle = p / 3;
            let field = p % 3;
            plus[circle][field] += h_um;
            minus[circle][field] -= h_um;
            let rp = intersect_circles(plus[0], plus[1]).unwrap().points[0];
            let rm = intersect_circles(minus[0], minus[1]).unwrap().points[0];
            for axis in 0..2 {
                let fd = (rp.point_mm[axis] - rm.point_mm[axis]) / 0.002;
                near(point.jet_mm_per_mm[p][axis], fd, 2e-6);
                let predicted = point.point_mm[axis] + point.jet_mm_per_mm[p][axis] * 0.001;
                let err = (rp.point_mm[axis] - predicted).abs();
                assert!(err < 2e-6, "first-order Taylor residual {err}");
            }
        }
        assert!(point.constraint_residual_mm2 <= 1e-12);
        assert!(point.condition_estimate.is_finite());
        assert_eq!(result.rank, Some(2));
        assert!(intersect_circles([0, 0, 5], [10, 0, 5])
            .unwrap()
            .points
            .is_empty());
    }

    #[test]
    fn factored_discriminant_resolves_one_micrometre_from_external_tangency() {
        let result = intersect_circles([0, 0, 5000], [9999, 0, 5000]).unwrap();
        assert_eq!(result.cardinality, IntersectionCardinality::Two);
        let expected =
            (5000.0_f64.powi(2) * 9999.0_f64.powi(2) - ((9999.0_f64.powi(2)) / 2.0).powi(2)).sqrt()
                / 9999.0
                / 1000.0;
        near(result.points[0].point_mm[1], expected, 1e-12);
        near(result.points[1].point_mm[1], -expected, 1e-12);
        assert!(result.points[0].condition_estimate > 1.0);
    }

    #[test]
    fn second_order_taylor_remainder_scales_quadratically() {
        let circles = [[0_i64, 0, 5000], [6000, 1000, 5000]];
        let base = intersect_circles(circles[0], circles[1]).unwrap().points[0];
        let residual = |step: i64| {
            let perturbed = intersect_circles(
                [circles[0][0], circles[0][1], circles[0][2] + step],
                circles[1],
            )
            .unwrap();
            let point = perturbed.points[0].point_mm;
            let h_mm = step as f64 / 1000.0;
            ((point[0] - base.point_mm[0] - base.jet_mm_per_mm[2][0] * h_mm).powi(2)
                + (point[1] - base.point_mm[1] - base.jet_mm_per_mm[2][1] * h_mm).powi(2))
            .sqrt()
        };
        let e1 = residual(100);
        let e2 = residual(200);
        let ratio = e2 / e1;
        assert!((ratio - 4.0).abs() < 0.15, "Taylor remainder ratio {ratio}");
    }
}
