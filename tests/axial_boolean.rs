use solidcraft_trivariate::{
    axial_boolean::{
        build_axial_boolean, AxialAxis, AxialBooleanOperation as Op, AxialOperand,
        AxialSectionFamily as Family, OperandLabel,
    },
    axial_boolean_audit::audit_axial_boolean,
    RationalBezierMap,
};

fn operand(section: Family, interval_mm: [i32; 2]) -> AxialOperand {
    AxialOperand {
        section,
        axis: AxialAxis::Z,
        interval_mm,
    }
}

#[test]
fn selected_geometry_and_components_match_unit_slice_oracle() {
    for family in [Family::ReferenceRectangleV1, Family::PiercedSquareV1] {
        let area = if family == Family::ReferenceRectangleV1 {
            400.
        } else {
            400. - 16. * std::f64::consts::PI
        };
        for (first, second) in [
            ([0, 8], [0, 8]),
            ([-8, 4], [-4, 0]),
            ([-4, 0], [-8, 4]),
            ([-7, 2], [-3, 4]),
            ([-3, 4], [-7, 2]),
            ([-8, -3], [2, 5]),
            ([-8, 0], [0, 5]),
            ([0, 5], [-8, 0]),
        ] {
            let a = operand(family, first);
            let b = operand(family, second);
            for op in [Op::Union, Op::Intersection, Op::Difference] {
                let model = build_axial_boolean(op, a, b).unwrap();
                audit_axial_boolean(op, a, b, &model).unwrap();
                let mut length = 0;
                let mut components = 0;
                let mut previous = false;
                for z in -9..10 {
                    let in_a = (first[0]..first[1]).contains(&z);
                    let in_b = (second[0]..second[1]).contains(&z);
                    let expected = match op {
                        Op::Union => in_a || in_b,
                        Op::Difference => in_a && !in_b,
                        Op::Intersection => in_a && in_b,
                    };
                    length += i32::from(expected);
                    components += usize::from(expected && !previous);
                    previous = expected;
                    let count = model
                        .blocks
                        .iter()
                        .filter(|blk| {
                            let z_of =
                                |i: usize| blk.controls[i][2] as f64 / blk.controls[i][3] as f64;
                            (z_of(0)..z_of(1)).contains(&(f64::from(z) + 0.5))
                        })
                        .count();
                    assert_eq!(
                        count,
                        usize::from(expected)
                            * if family == Family::ReferenceRectangleV1 {
                                1
                            } else {
                                4
                            }
                    );
                }
                assert_eq!(model.component_count, components);
                assert!((model.volume_mm3 - f64::from(length) * area).abs() < 1e-11 * area);
                // At every block the physical z-coordinate agrees with the affine
                // elevation interpolation and the axial Jacobian is independent of u,v.
                for blk in &model.blocks {
                    let p = |i: usize| blk.controls[i][2] as f64 / blk.controls[i][3] as f64;
                    let map = RationalBezierMap::new(
                        blk.degrees,
                        blk.controls.iter().map(|v| v.map(|x| x as f64)).collect(),
                    )
                    .unwrap();
                    for w in [0., 0.375, 1.] {
                        let eval = map.evaluate([0.3, 0.7, w]).unwrap();
                        assert!((eval.position[2] - (p(0) + (p(1) - p(0)) * w)).abs() < 1e-12);
                        assert!((eval.jacobian_columns[2][2] - (p(1) - p(0))).abs() < 1e-12);
                    }
                }
            }
        }
    }
}

#[test]
fn tool_caps_keep_negative_orientation_and_damaged_models_are_rejected() {
    let a = operand(Family::PiercedSquareV1, [0, 10]);
    let b = operand(Family::PiercedSquareV1, [4, 6]);
    let model = build_axial_boolean(Op::Difference, a, b).unwrap();
    let caps = model
        .faces
        .iter()
        .flat_map(|f| &f.provenance)
        .filter(|p| p.operand == OperandLabel::B)
        .collect::<Vec<_>>();
    assert_eq!(caps.len(), 8);
    assert!(caps.iter().all(|p| p.orientation_reversed));
    assert_eq!(caps.iter().filter(|p| p.source_face == "z0").count(), 4);
    for k in 0..6 {
        let mut damaged = model.clone();
        match k {
            0 => damaged.blocks[0].controls[2][0] += 1,
            1 => damaged.blocks[0].controls[5][2] += 1,
            2 => damaged.faces[0].provenance[0].orientation_reversed ^= true,
            3 => {
                damaged.incidence.exterior.pop();
            }
            4 => damaged.volume_mm3 = f64::NAN,
            _ => {
                damaged.selected_intervals[0].operands.clear();
            }
        }
        assert!(audit_axial_boolean(Op::Difference, a, b, &damaged).is_err());
    }
}

#[test]
fn canonical_json_round_trip_keeps_split_caps_and_audits_empty_solids() {
    use solidcraft_trivariate::{
        interchange::{audit, build},
        json::{parse, Json},
    };
    let request = |operation: &str| {
        parse(&format!(r#"{{"family":"axial_boolean_v1","section":"pierced_square_v1","operation":"{operation}","a_z_mm":[-8,0],"b_z_mm":[0,8]}}"#)).unwrap()
    };
    let empty = build(&request("intersection")).unwrap();
    assert!(empty.field("blocks").unwrap().items().unwrap().is_empty());
    assert_eq!(
        audit(&empty)
            .unwrap()
            .field("components")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        0
    );
    let mut joined = build(&request("union")).unwrap();
    assert_eq!(
        audit(&joined)
            .unwrap()
            .field("components")
            .unwrap()
            .number::<usize>()
            .unwrap(),
        1
    );
    if let Json::Object(map) = &mut joined {
        map.insert("intervals".into(), Json::Array(Vec::new()));
    }
    assert!(audit(&joined).is_err());
}

#[test]
fn invalid_geometry_and_unsupported_recipes_have_specific_refusals() {
    use solidcraft_trivariate::axial_boolean::AxialBooleanError as Error;
    let a = operand(Family::PiercedSquareV1, [-1000, 1000]);
    for (b, error) in [
        (operand(a.section, [1, 1]), Error::NullInterval),
        (operand(a.section, [2, 1]), Error::InvertedInterval),
        (operand(a.section, [-1001, 0]), Error::IntervalOutOfRange),
        (
            operand(Family::ReferenceRectangleV1, [0, 1]),
            Error::DifferentSections,
        ),
        (
            AxialOperand {
                axis: AxialAxis::X,
                ..a
            },
            Error::IncompatibleAxes,
        ),
    ] {
        assert_eq!(build_axial_boolean(Op::Union, a, b), Err(error));
    }
    let bad = operand(Family::CoreDiskV1, [0, 1]);
    assert_eq!(
        build_axial_boolean(Op::Difference, bad, bad),
        Err(Error::UnsupportedSection)
    );
    let bound = build_axial_boolean(Op::Union, a, a).unwrap();
    audit_axial_boolean(Op::Union, a, a, &bound).unwrap();
}
