use solidcraft_trivariate::constraint_system::Ratio;
use solidcraft_trivariate::incidence::enumerate_faces;
use solidcraft_trivariate::nurbs_cross::{build_cross, build_cross_bounded, NurbsCurve2};

const BOX: [[i64; 2]; 3] = [[0, 20], [0, 15], [0, 22]];

fn line(x: i64, a: i64, b: i64, y0: i64, y1: i64) -> NurbsCurve2 {
    NurbsCurve2 {
        degree: 1,
        knots: vec![a, a, b, b],
        controls: vec![[x, y0, 1], [x, y1, 1]],
    }
}

#[test]
fn unequal_degrees_knot_vectors_and_cap_levels() {
    let xmin = NurbsCurve2 {
        degree: 2,
        knots: vec![0, 0, 0, 4, 4, 4],
        controls: vec![[2, 2, 1], [4, 6, 1], [3, 10, 1]],
    };
    let xmax = NurbsCurve2 {
        degree: 1,
        knots: vec![0, 0, 2, 4, 4],
        controls: vec![[12, 2, 1], [14, 6, 1], [12, 10, 1]],
    };
    let zmin = line(3, 1, 3, 4, 8);
    let zmax = line(10, 1, 3, 4, 8);
    let blocks = build_cross(BOX, [&xmin, &xmax, &zmin, &zmax]).unwrap();
    assert_eq!(blocks.len(), 38);
    // Eighteen stock end-cap blocks, two bands with one active feature,
    // two bands with both features active. Check normalized and elevated traces.
    let faces = enumerate_faces(&blocks).unwrap();
    assert_eq!(faces.exterior.len() + 2 * faces.interfaces.len(), 6 * 38);
    assert!(blocks.iter().any(|b| b.degrees == [1, 2, 1]));
}

#[test]
fn separated_parameter_intervals_fill_material_between_bores() {
    let curves = [
        line(3, 0, 1, 2, 4),
        line(9, 0, 1, 2, 4),
        line(5, 2, 3, 7, 11),
        line(13, 2, 3, 7, 11),
    ];
    let blocks = build_cross(BOX, curves.each_ref()).unwrap();
    assert_eq!(blocks.len(), 9 + 6 + 9 + 6 + 9);
    let mut reversed = curves;
    reversed.swap(0, 1);
    assert!(build_cross(BOX, reversed.each_ref()).is_err());
}

#[test]
fn opt_in_approximation_reconciles_independent_ordinates() {
    let simple = |x| NurbsCurve2 {
        degree: 2,
        knots: vec![0, 0, 0, 4, 4, 4],
        controls: vec![[x, 2, 1], [x + 2, 5, 1], [x, 9, 1]],
    };
    let mut curves = [simple(3), simple(12), simple(5), simple(16)];
    curves[2].controls[1][1] = 6;
    assert!(build_cross(BOX, curves.each_ref()).is_err());
    let tol = Ratio::new(1, 4).unwrap();
    assert!(build_cross_bounded(BOX, curves.each_ref(), tol, 0).is_err());
    let blocks = build_cross_bounded(BOX, curves.each_ref(), tol, 8).unwrap();
    assert!(blocks.iter().all(|b| b.degrees == [1; 3]));
    enumerate_faces(&blocks).unwrap();
}

#[test]
fn reflected_meridians_have_regular_poles_and_conforming_interfaces() {
    use solidcraft_trivariate::nurbs_cross::build_oval;
    let x = NurbsCurve2 {
        degree: 3,
        knots: vec![0, 0, 0, 0, 3, 3, 3, 3],
        controls: vec![[6, 0, 1], [12, 4, 2], [2, 5, 1], [0, 5, 1]],
    };
    let center = [10, 7, 11];
    let blocks = build_oval(BOX, center, &x, &x).unwrap();
    let faces = enumerate_faces(&blocks).unwrap();
    assert_eq!(
        (blocks.len(), faces.interfaces.len(), faces.exterior.len()),
        (40, 68, 104)
    );
    let mut unsynchronized = x.clone();
    unsynchronized.controls[1][1] += 1;
    assert!(build_oval(BOX, center, &x, &unsynchronized).is_err());
    let mut corner = x.clone();
    corner.controls[1][0] -= 1;
    assert!(build_oval(BOX, center, &corner, &corner).is_err());
}
