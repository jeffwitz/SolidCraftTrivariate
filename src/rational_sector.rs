//! Exact integer-homogeneous charts for an annular sector prism.
use crate::bracket::BracketBlock;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RationalSector {
    /// Numerator and denominator of `tan(theta/2)`.
    pub half_angle_tangent: [i32; 2],
    /// Inner and outer radii in millimetres.
    pub radii_mm: [i32; 2],
    pub height_mm: i32,
}

type HPoint = [i128; 4];

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn normalize(id: &'static str, mut controls: Vec<HPoint>) -> Result<BracketBlock, &'static str> {
    if controls.len() != 12 || controls.iter().any(|p| p[3] <= 0) {
        return Err("sector_invalid_control_net");
    }
    let common = controls
        .iter()
        .flatten()
        .fold(0_u128, |g, &x| gcd(g, x.unsigned_abs()));
    if common == 0 {
        return Err("sector_zero_control_net");
    }
    let divisor = i128::try_from(common).map_err(|_| "sector_integer_overflow")?;
    for p in &mut controls {
        for x in p {
            *x /= divisor;
        }
    }
    let controls = controls
        .into_iter()
        .map(|p| {
            Ok([
                i64::try_from(p[0]).map_err(|_| "sector_control_exceeds_i64")?,
                i64::try_from(p[1]).map_err(|_| "sector_control_exceeds_i64")?,
                i64::try_from(p[2]).map_err(|_| "sector_control_exceeds_i64")?,
                i64::try_from(p[3]).map_err(|_| "sector_control_exceeds_i64")?,
            ])
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(BracketBlock {
        id: id.into(),
        degrees: [2, 1, 1],
        controls,
    })
}

fn restrict(source: &[HPoint], second: bool) -> Vec<HPoint> {
    let mut result = Vec::with_capacity(12);
    for u in 0..3 {
        for j in 0..4 {
            let a = source[j];
            let b = source[4 + j];
            let c = source[8 + j];
            let middle = std::array::from_fn(|k| a[k] + 2 * b[k] + c[k]);
            let p = match (second, u) {
                (false, 0) => a.map(|x| 4 * x),
                (false, 1) => std::array::from_fn(|k| 2 * (a[k] + b[k])),
                (false, 2) | (true, 0) => middle,
                (true, 1) => std::array::from_fn(|k| 2 * (b[k] + c[k])),
                (true, 2) => c.map(|x| 4 * x),
                _ => unreachable!(),
            };
            result.push(p);
        }
    }
    result
}

fn validate(recipe: &RationalSector) -> Result<(i128, i128, i128, i128, bool), &'static str> {
    let [n, d] = recipe.half_angle_tangent.map(i128::from);
    let [ri, ro] = recipe.radii_mm.map(i128::from);
    let h = i128::from(recipe.height_mm);
    if !(1..=256).contains(&n)
        || !(1..=256).contains(&d)
        || 4 * n < d
        || n > 4 * d
        || !(1..ro).contains(&ri)
        || ro > 64
        || !(1..=64).contains(&h)
    {
        return Err("sector_parameters_outside_domain");
    }
    // t > 1 + sqrt(2), compared exactly without approximating sqrt(2).
    let split = n > d && (n - d) * (n - d) > 2 * d * d;
    Ok((n, d, ri, ro, split))
}

/// Construct one or two exact polynomial charts of the rational polar sector.
pub fn build_rational_sector(recipe: &RationalSector) -> Result<Vec<BracketBlock>, &'static str> {
    let (n, d, ri, ro, split) = validate(recipe)?;
    let q = 4 * d;
    let tangent = 4 * n;
    let q2 = q * q;
    let t2 = tangent * tangent;
    let h = i128::from(recipe.height_mm);
    let mut base = Vec::with_capacity(12);
    // Bernstein coefficients of (Q²-S², 2QS, Q²+S²), S=T*u.
    let arc = [
        [q2, 0, q2],
        [q2, 4 * q * n, q2],
        [q2 - t2, 8 * q * n, q2 + t2],
    ];
    for angle in arc {
        for radius in [ro, ri] {
            for z in [0, h] {
                base.push([radius * angle[0], radius * angle[1], z * angle[2], angle[2]]);
            }
        }
    }
    if !split {
        return Ok(vec![normalize("whole", base)?]);
    }
    let first = restrict(&base, false);
    let second = restrict(&base, true);
    Ok(vec![
        normalize("first", first)?,
        normalize("second", second)?,
    ])
}

/// Exactly restrict the accepted whole chart in the global angular parameter.
pub(crate) fn split_for_analysis(block: &BracketBlock) -> Result<Vec<BracketBlock>, &'static str> {
    if block.id != "whole" || block.degrees != [2, 1, 1] || block.controls.len() != 12 {
        return Err("angular subdivision requires one quadratic sector chart");
    }
    let controls: Vec<_> = block.controls.iter().map(|p| p.map(i128::from)).collect();
    Ok(vec![
        normalize("first", restrict(&controls, false))?,
        normalize("second", restrict(&controls, true))?,
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rational_sector_audit::audit_rational_sector;

    fn recipe(t: [i32; 2]) -> RationalSector {
        RationalSector {
            half_angle_tangent: t,
            radii_mm: [3, 6],
            height_mm: 8,
        }
    }

    #[test]
    fn charts_switch_exactly_across_the_135_degree_threshold() {
        let below = build_rational_sector(&recipe([70, 29])).unwrap();
        let above = build_rational_sector(&recipe([169, 70])).unwrap();
        assert_eq!(
            below.iter().map(|b| b.id.as_str()).collect::<Vec<_>>(),
            ["whole"]
        );
        assert_eq!(
            above.iter().map(|b| b.id.as_str()).collect::<Vec<_>>(),
            ["first", "second"]
        );
        audit_rational_sector(&recipe([70, 29]), &below).unwrap();
        audit_rational_sector(&recipe([169, 70]), &above).unwrap();
    }

    #[test]
    fn exact_chart_controls_are_integer_positive_and_oracle_accepted() {
        for t in [[1, 4], [1, 1], [3, 2], [7, 5], [4, 1]] {
            let r = recipe(t);
            let blocks = build_rational_sector(&r).unwrap();
            assert!(blocks.iter().all(|b| b.controls.iter().all(|p| p[3] > 0)));
            audit_rational_sector(&r, &blocks).unwrap();
        }
    }

    #[test]
    fn invalid_rational_domains_are_rejected() {
        for r in [
            RationalSector {
                half_angle_tangent: [0, 1],
                ..recipe([1, 1])
            },
            RationalSector {
                half_angle_tangent: [1, 0],
                ..recipe([1, 1])
            },
            RationalSector {
                half_angle_tangent: [1, 5],
                ..recipe([1, 1])
            },
            RationalSector {
                radii_mm: [6, 6],
                ..recipe([1, 1])
            },
            RationalSector {
                height_mm: 0,
                ..recipe([1, 1])
            },
        ] {
            assert!(build_rational_sector(&r).is_err());
        }
    }

    #[test]
    fn regression_old_endpoint_y_coefficient_is_rejected() {
        let r = RationalSector {
            half_angle_tangent: [1, 4],
            radii_mm: [1, 2],
            height_mm: 1,
        };
        let mut old = build_rational_sector(&r).unwrap();
        assert_eq!(old[0].controls[8], [30, 16, 0, 17]);
        old[0].controls[8][1] = 8;
        assert_eq!(
            audit_rational_sector(&r, &old),
            Err("sector_homogeneous_polynomial_identity")
        );
    }
}
