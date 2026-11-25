use solidcraft_trivariate::{
    constraint_system::Ratio,
    incidence::enumerate_faces,
    primitive_template::{Operation, Primitive},
    sketch_junction::plan_loop,
    spatial_certificate::{certify_spatial, CertificateStatus},
};

fn r(n: i128) -> Ratio {
    Ratio::new(n, 1).unwrap()
}

#[test]
fn templates_predict_cells_and_every_internal_interface_before_sweeping() {
    let primitives = [
        (
            Primitive::Rectangle {
                min: [r(4), r(-1)],
                max: [r(8), r(1)],
            },
            1,
            0,
        ),
        (
            Primitive::Circle {
                center: [r(6), r(0)],
                radius: r(2),
            },
            5,
            8,
        ),
        (
            Primitive::Ellipse {
                center: [r(6), r(0)],
                radii: [r(2), r(1)],
            },
            5,
            8,
        ),
    ];
    for (primitive, sections, planar_interfaces) in primitives {
        for (operation, sectors) in [
            (Operation::Extrude { length: 3 }, 1),
            (Operation::Revolve, 4),
        ] {
            let plan = primitive.plan(operation).unwrap();
            assert_eq!(plan.section_patches.len(), sections);
            assert_eq!(plan.block_count(), sections * sectors);
            let c = &plan.contour;
            let observed = plan_loop(c.degree, &c.knots, &c.controls).unwrap();
            assert_eq!(
                observed.iter().map(|j| j.continuity).collect::<Vec<_>>(),
                plan.junctions
            );
            let blocks = plan.build().unwrap();
            assert_eq!(blocks.len(), plan.block_count());
            let joins = enumerate_faces(&blocks).unwrap();
            assert_eq!(
                joins.interfaces.len(),
                sectors * planar_interfaces + if sectors == 4 { sectors * sections } else { 0 }
            );
            for b in blocks {
                assert_eq!(
                    certify_spatial(b.degrees, &b.controls, Default::default()).status,
                    CertificateStatus::Valid,
                    "{}",
                    b.id
                );
            }
        }
    }
}

#[test]
fn a_template_requires_a_nondegenerate_section_and_positive_axis_clearance() {
    let touching = Primitive::Circle {
        center: [r(2), r(0)],
        radius: r(2),
    };
    assert!(touching.plan(Operation::Revolve).is_err());
    assert!(touching.plan(Operation::Extrude { length: 3 }).is_ok());
    let rectangle = Primitive::Rectangle {
        min: [r(0), r(0)],
        max: [r(0), r(10)],
    };
    assert!(rectangle.plan(Operation::Extrude { length: 3 }).is_err());
}

#[test]
fn swept_volumes_match_the_area_and_centroid_of_each_source_primitive() {
    use solidcraft_trivariate::{elasticity::gauss, RationalBezierMap};
    use std::f64::consts::PI;
    let quadrature = gauss(12).unwrap();
    let cases = [
        (
            "rectangle",
            Primitive::Rectangle {
                min: [r(4), r(-1)],
                max: [r(8), r(1)],
            },
            8.0,
        ),
        (
            "circle",
            Primitive::Circle {
                center: [r(6), r(0)],
                radius: r(2),
            },
            4.0 * PI,
        ),
        (
            "ellipse",
            Primitive::Ellipse {
                center: [r(6), r(0)],
                radii: [r(2), r(1)],
            },
            2.0 * PI,
        ),
    ];
    for (name, primitive, area) in cases {
        for (operation, length) in [
            (Operation::Extrude { length: 3 }, 3.0),
            (Operation::Revolve, 12.0 * PI),
        ] {
            let mut volume = 0.0;
            for b in primitive.plan(operation).unwrap().build().unwrap() {
                let map = RationalBezierMap::new(
                    b.degrees,
                    b.controls.iter().map(|p| p.map(|a| a as f64)).collect(),
                )
                .unwrap();
                for &(u, wu) in &quadrature {
                    for &(v, wv) in &quadrature {
                        for &(w, ww) in &quadrature {
                            let [a, b, c] = map.evaluate([u, v, w]).unwrap().jacobian_columns;
                            let det = a[0] * (b[1] * c[2] - b[2] * c[1])
                                + a[1] * (b[2] * c[0] - b[0] * c[2])
                                + a[2] * (b[0] * c[1] - b[1] * c[0]);
                            volume += wu * wv * ww * det;
                        }
                    }
                }
            }
            let expected = area * length;
            let error = (volume - expected).abs() / expected;
            println!("{name} {operation:?}: volume={volume:.12}, reference={expected:.12}, relative={error:.3e}");
            assert!(error < 1e-9);
        }
    }
}
