//! Exact radial subdivision guide for the pocket attached to the rectangular channel.
use crate::constraint_system::Ratio;

type Point = [Ratio; 3];
type Failure = &'static str;

fn divide(a: Ratio, b: Ratio) -> Result<Ratio, Failure> {
    a.checked_mul(Ratio::new(b.denominator(), b.numerator())?)
}

fn cartesian(row: &Point) -> Result<[Ratio; 2], Failure> {
    if row[2].numerator() <= 0 {
        return Err("radial guide requires positive weights");
    }
    Ok([divide(row[0], row[2])?, divide(row[1], row[2])?])
}

fn mix(a: Ratio, b: Ratio, fraction: Ratio) -> Result<Ratio, Failure> {
    a.checked_add(b.checked_sub(a)?.checked_mul(fraction)?)
}

fn horizontal_fraction(
    outer: [Ratio; 2],
    inner: [Ratio; 2],
    stock_xmin: Ratio,
    pocket_xmin: Ratio,
) -> Result<Option<Ratio>, Failure> {
    if outer[0] != stock_xmin || outer[1] != inner[1] {
        return Ok(None);
    }
    if !stock_xmin.cmp_value(&pocket_xmin).is_lt() || !pocket_xmin.cmp_value(&inner[0]).is_lt() {
        return Err("pocket guide must lie strictly inside the horizontal ray");
    }
    Ok(Some(divide(
        pocket_xmin.checked_sub(stock_xmin)?,
        inner[0].checked_sub(stock_xmin)?,
    )?))
}

/// At each horizontal ray on the left wing the guide is exactly x=pocket_xmin.
/// Other sectors interpolate the two endpoint fractions. Endpoint construction
/// uses only that ray so reversing or rotating the ring preserves its trace.
pub(crate) fn radial_guide(
    outer: &[Point],
    inner: &[Point],
    stock_xmin: Ratio,
    pocket_xmin: Ratio,
) -> Result<Vec<Point>, Failure> {
    if outer.len() != inner.len() || !(2..=4).contains(&outer.len()) {
        return Err("radial guide requires matching linear to cubic rows");
    }
    let last = outer.len() - 1;
    let anchor = |i: usize| {
        horizontal_fraction(
            cartesian(&outer[i])?,
            cartesian(&inner[i])?,
            stock_xmin,
            pocket_xmin,
        )
    };
    let half = Ratio::new(1, 2)?;
    let first_fraction = anchor(0)?.unwrap_or(half);
    let last_fraction = anchor(last)?.unwrap_or(half);
    (0..=last)
        .map(|i| {
            let a = cartesian(&outer[i])?;
            let b = cartesian(&inner[i])?;
            let fraction = match horizontal_fraction(a, b, stock_xmin, pocket_xmin)? {
                Some(value) => value,
                None => mix(
                    first_fraction,
                    last_fraction,
                    Ratio::new(i as i128, last as i128)?,
                )?,
            };
            let w = inner[i][2];
            Ok([
                mix(a[0], b[0], fraction)?.checked_mul(w)?,
                mix(a[1], b[1], fraction)?.checked_mul(w)?,
                w,
            ])
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn n(i: i128) -> Ratio {
        Ratio::new(i, 1).unwrap()
    }

    #[test]
    fn fixes_the_left_cut_and_preserves_shared_endpoints_under_reversal() {
        let outer = [[n(0), n(36), n(6)], [n(0), n(63), n(9)]];
        let inner = [[n(132), n(12), n(2)], [n(201), n(21), n(3)]];
        let guide = radial_guide(&outer, &inner, n(0), n(55)).unwrap();
        assert_eq!(guide, [[n(110), n(12), n(2)], [n(165), n(21), n(3)]]);
        let neighbor = radial_guide(
            &[[n(20), n(0), n(2)], outer[1]],
            &[[n(100), n(9), n(2)], inner[1]],
            n(0),
            n(55),
        )
        .unwrap();
        assert_eq!(neighbor[1], guide[1]);
        assert!(radial_guide(&outer, &inner, n(0), n(67)).is_err());
    }

    #[test]
    fn splitting_the_source_ring_preserves_every_angular_trace() {
        use crate::bracket::BracketBlock;
        use crate::spatial_certificate::{certify_spatial, CertificateStatus};

        let model =
            crate::json::parse(include_str!("../docs/results/g2-cross-plate.model.json")).unwrap();
        let mut children = Vec::new();
        for block in crate::bracket_json::decode_blocks(&model).unwrap() {
            if block.controls[0][2] != 0 {
                continue;
            }
            let row = |v: usize| -> Vec<Point> {
                block
                    .controls
                    .as_chunks::<4>()
                    .0
                    .iter()
                    .map(|slab| {
                        let [x, y, _, w] = slab[v * 2];
                        [n(x.into()), n(y.into()), n(w.into())]
                    })
                    .collect()
            };
            let outer = row(0);
            let inner = row(1);
            let guide = radial_guide(&outer, &inner, n(0), n(55)).unwrap();
            for (side, (a, b)) in [(&outer, &guide), (&guide, &inner)].into_iter().enumerate() {
                let mut controls = Vec::new();
                for (first, second) in a.iter().zip(b) {
                    for &[x, y, w] in [first, second] {
                        controls.push([x, y, n(0), w]);
                        controls.push([x, y, w, w]);
                    }
                }
                let controls = crate::nurbs_curve::normalize(&controls).unwrap();
                let cert = certify_spatial([3, 1, 1], &controls, Default::default());
                assert_eq!(
                    cert.status,
                    CertificateStatus::Valid,
                    "{} side {side}",
                    block.id
                );
                children.push(BracketBlock {
                    id: format!("{}_{side}", block.id),
                    degrees: [3, 1, 1],
                    controls,
                });
            }
        }
        let faces = crate::incidence::enumerate_faces(&children).unwrap();
        assert_eq!(
            (children.len(), faces.interfaces.len(), faces.exterior.len()),
            (24, 36, 72)
        );
    }
}
