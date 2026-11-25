use solidcraft_trivariate::{
    bracket::BracketBlock, steinmetz_frame::intersecting_cylinders, RationalBezierMap,
};

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}

fn sample(block: &BracketBlock, uvw: [f64; 3]) -> ([f64; 3], f64) {
    let controls = block.controls.iter().map(|p| p.map(|x| x as f64)).collect();
    let evaluation = RationalBezierMap::new(block.degrees, controls)
        .unwrap()
        .evaluate(uvw)
        .unwrap();
    let determinant = dot(
        evaluation.jacobian_columns[0],
        cross(
            evaluation.jacobian_columns[1],
            evaluation.jacobian_columns[2],
        ),
    );
    (evaluation.position, determinant)
}

fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    let mut rule = Vec::with_capacity(n);
    for i in 0..n {
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..12 {
            let (mut p0, mut p1) = (1., x);
            for degree in 2..=n {
                let p =
                    ((2 * degree - 1) as f64 * x * p1 - (degree - 1) as f64 * p0) / degree as f64;
                p0 = p1;
                p1 = p;
            }
            let derivative = n as f64 * (x * p1 - p0) / (x * x - 1.);
            let update = p1 / derivative;
            x -= update;
            if update.abs() < 2e-16 {
                break;
            }
        }
        let (mut p0, mut p1) = (1., x);
        for degree in 2..=n {
            let p = ((2 * degree - 1) as f64 * x * p1 - (degree - 1) as f64 * p0) / degree as f64;
            p0 = p1;
            p1 = p;
        }
        let derivative = n as f64 * (x * p1 - p0) / (x * x - 1.);
        let weight = 1. / ((1. - x * x) * derivative * derivative);
        rule.push(((x + 1.) * 0.5, weight));
    }
    rule
}

fn volume(blocks: &[BracketBlock], order: usize) -> f64 {
    let rule = gauss_legendre(order);
    blocks
        .iter()
        .map(|block| {
            let controls = block.controls.iter().map(|p| p.map(|x| x as f64)).collect();
            let map = RationalBezierMap::new(block.degrees, controls).unwrap();
            let mut sum = 0.;
            for &(u, wu) in &rule {
                for &(v, wv) in &rule {
                    for &(w, ww) in &rule {
                        let jacobian = map.evaluate([u, v, w]).unwrap().jacobian_columns;
                        sum += wu * wv * ww * dot(jacobian[0], cross(jacobian[1], jacobian[2]));
                    }
                }
            }
            sum
        })
        .sum()
}

fn outside_cylinder(point: [f64; 3], center: [i32; 3], axis: usize, radius: i32) -> bool {
    (0..3)
        .filter(|&k| k != axis)
        .map(|k| (point[k] - f64::from(center[k])).powi(2))
        .sum::<f64>()
        >= f64::from(radius * radius) - 1e-9
}

fn check_pair(axes: [usize; 2]) -> Vec<BracketBlock> {
    let bounds = [[-21, 27], [-16, 19], [-12, 24]];
    let center = [3, -2, 5];
    let radius = 3;
    let blocks = intersecting_cylinders(bounds, center, radius, 7, axes, "pair").unwrap();
    assert_eq!(blocks.len(), 40);
    let mut ids = blocks.iter().map(|b| b.id.as_str()).collect::<Vec<_>>();
    ids.sort_unstable();
    ids.dedup();
    assert_eq!(ids.len(), blocks.len());
    for block in &blocks {
        assert!(block.controls.iter().all(|p| p[3] > 0));
        for u in [0., 0.5, 1.] {
            for v in [0., 0.5, 1.] {
                for w in [0., 0.5, 1.] {
                    let (p, determinant) = sample(block, [u, v, w]);
                    assert!(
                        determinant > 0.,
                        "{} at {u},{v},{w}: {determinant}",
                        block.id
                    );
                    for axis in 0..3 {
                        assert!(p[axis] >= f64::from(bounds[axis][0]) - 1e-10);
                        assert!(p[axis] <= f64::from(bounds[axis][1]) + 1e-10);
                    }
                    assert!(
                        outside_cylinder(p, center, axes[0], radius),
                        "{}: {p:?}",
                        block.id
                    );
                    assert!(
                        outside_cylinder(p, center, axes[1], radius),
                        "{}: {p:?}",
                        block.id
                    );
                }
            }
        }
    }
    blocks
}

#[test]
fn all_axis_pairs_build_in_non_cubic_off_origin_boxes() {
    for axes in [[0, 2], [0, 1], [1, 2], [2, 0], [1, 0], [2, 1]] {
        check_pair(axes);
    }
}

#[test]
fn reversing_the_axis_order_keeps_the_same_analytic_material_region() {
    let forward = check_pair([0, 2]);
    let reversed = check_pair([2, 0]);
    assert_eq!(forward.len(), reversed.len());
    let bounds = [[-21, 27], [-16, 19], [-12, 24]];
    let center = [3, -2, 5];
    let radius = 3;
    for (left, right) in forward.iter().zip(&reversed) {
        for uvw in [[0.23, 0.37, 0.61], [0.73, 0.5, 0.19]] {
            let (a, _) = sample(left, uvw);
            let (b, _) = sample(right, uvw);
            for p in [a, b] {
                for axis in 0..3 {
                    assert!(p[axis] >= f64::from(bounds[axis][0]) - 1e-10);
                    assert!(p[axis] <= f64::from(bounds[axis][1]) + 1e-10);
                }
                assert!(outside_cylinder(p, center, 0, radius));
                assert!(outside_cylinder(p, center, 2, radius));
            }
        }
    }
}

#[test]
fn invalid_axis_pairs_and_ranges_refuse() {
    let bounds = [[-10, 10]; 3];
    assert!(intersecting_cylinders(bounds, [0; 3], 2, 4, [0, 0], "x").is_err());
    assert!(intersecting_cylinders(bounds, [0; 3], 2, 4, [0, 3], "x").is_err());
    assert!(intersecting_cylinders(bounds, [0; 3], 2, 4, [0, 1], "").is_err());
}

#[test]
fn radius_one_asymmetric_frame_preserves_the_steinmetz_thirds_term() {
    let bounds = [[-8, 9], [-6, 5], [-9, 13]];
    let center = [1, 0, 2];
    let radius = 1;
    let blocks = intersecting_cylinders(bounds, center, radius, 4, [1, 2], "unit").unwrap();
    let lengths = bounds.map(|[lo, hi]| f64::from(hi - lo));
    let expected = lengths[0] * lengths[1] * lengths[2]
        - std::f64::consts::PI * f64::from(radius * radius) * (lengths[1] + lengths[2])
        + 16. * f64::from(radius.pow(3)) / 3.;
    assert!((volume(&blocks, 24) - expected).abs() < 1e-8);
}
