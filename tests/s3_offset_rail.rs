use solidcraft_trivariate::s3_offset_rail;

#[test]
fn named_lip_exports_two_certified_continuous_whole_span_rails() {
    for side in [-1, 1] {
        let rail = s3_offset_rail::build(side).unwrap();
        assert!(rail.spans.len() <= 32);
        assert_eq!(rail.source_block, "s3_i_u_zbase_g2_s0_p0_z2_src0");
        assert_eq!(rail.continuity, "C0");
        assert_eq!(rail.spans.first().unwrap().interval[0], 0.0);
        assert_eq!(rail.spans.last().unwrap().interval[1], 1.0);
        for span in &rail.spans {
            assert!(span.error_bound_mm <= 1.0e-3);
            assert!(span.coefficient_error_bound_mm < 1.0e-9);
            assert!(span.controls.iter().flatten().all(|x| x.is_finite()));
            assert!(span.controls.iter().all(|p| p[2] == 1.0));
        }
        for pair in rail.spans.windows(2) {
            assert_eq!(pair[0].interval[1], pair[1].interval[0]);
            assert_eq!(pair[0].controls[3], pair[1].controls[0]);
        }
        let expected_y = if side < 0 { 34.25 } else { 33.75 };
        assert_eq!(
            rail.spans.last().unwrap().controls[3],
            [75.0, expected_y, 1.0]
        );
    }
}
