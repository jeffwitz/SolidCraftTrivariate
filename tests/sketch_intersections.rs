use solidcraft_trivariate::sketch_intersections::intersect_circles;

#[test]
fn point_jets_satisfy_the_two_independent_linearized_circle_constraints() {
    let circles = [[-3000, 2000, 5000], [2000, 6000, 3000]];
    let result = intersect_circles(circles[0], circles[1]).unwrap();
    assert_eq!(result.points.len(), 2);
    for point in result.points {
        for (i, &[x, y, r]) in circles.iter().enumerate() {
            let dx = point.point_mm[0] - x as f64 / 1000.;
            let dy = point.point_mm[1] - y as f64 / 1000.;
            for (p, [px, py]) in point.jet_mm_per_mm.into_iter().enumerate() {
                let dcx = f64::from(u8::from(p == 3 * i));
                let dcy = f64::from(u8::from(p == 3 * i + 1));
                let dr = f64::from(u8::from(p == 3 * i + 2));
                assert!((dx * (px - dcx) + dy * (py - dcy) - r as f64 / 1000. * dr).abs() < 1e-12);
            }
        }
    }
}
