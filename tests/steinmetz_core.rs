use solidcraft_trivariate::elasticity::gauss;
use solidcraft_trivariate::incidence::enumerate_faces;
use solidcraft_trivariate::spatial_certificate::{
    certify_spatial, CertificateBudget, CertificateStatus,
};
use solidcraft_trivariate::steinmetz_core::steinmetz_box;
use solidcraft_trivariate::RationalBezierMap;

fn det([a, b, c]: [[f64; 3]; 3]) -> f64 {
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

#[test]
fn exact_circular_walls_positive_closed_cubes_and_union_volume() {
    let bounds = [[0, 60], [0, 40], [0, 30]];
    let center = [30, 20, 15];
    let blocks = steinmetz_box(bounds, center, 3, 4, "cross").unwrap();
    let incidence = enumerate_faces(&blocks).unwrap();
    assert_eq!(
        (
            blocks.len(),
            incidence.interfaces.len(),
            incidence.exterior.len()
        ),
        (40, 68, 104)
    );
    let quadrature = gauss(20).unwrap();
    let mut volume = 0.0;
    let mut lower_bound = f64::INFINITY;
    for block in blocks {
        let certificate =
            certify_spatial(block.degrees, &block.controls, CertificateBudget::default());
        assert_eq!(certificate.status, CertificateStatus::Valid, "{}", block.id);
        lower_bound = lower_bound.min(certificate.determinant_lower_bound_mm3.unwrap());
        let map = RationalBezierMap::new(
            block.degrees,
            block.controls.iter().map(|p| p.map(|x| x as f64)).collect(),
        )
        .unwrap();
        for &(u, wu) in &quadrature {
            for &(v, wv) in &quadrature {
                // All maps are affine along their physical straight rulings;
                // the determinant is constant along that parameter.
                volume += wu * wv * det(map.evaluate([u, v, 0.5]).unwrap().jacobian_columns);
            }
        }
        if block.degrees == [2, 1, 1] {
            for u in [0.0, 0.125, 0.5, 0.875, 1.0] {
                for s in [0.0, 0.5, 1.0] {
                    let [x, y, z] = map.evaluate([u, 0.0, s]).unwrap().position;
                    let wall = (y - 20.).powi(2) + (x - 30.).abs().min((z - 15.).abs()).powi(2);
                    assert!((wall - 9.).abs() < 1e-11);
                }
            }
        }
    }
    let expected = 72_000. - 810. * std::f64::consts::PI + 144.;
    eprintln!(
        "integrated volume={volume:.12}, residual={:.3e}, certified min det >= {lower_bound:.9}",
        volume - expected
    );
    assert!((volume - expected).abs() < 1e-7);
}

#[test]
fn off_center_axes_and_guard_limits() {
    let bounds = [[-40, 18], [-10, 12], [-15, 15]];
    let center = [-2, 2, 4];
    for g in [4, 8] {
        let blocks = steinmetz_box(bounds, center, 3, g, "offset").unwrap();
        assert_eq!(enumerate_faces(&blocks).unwrap().interfaces.len(), 68);
    }
    for g in [2, 3, 11] {
        assert!(steinmetz_box(bounds, center, 3, g, "offset").is_err());
    }
}
