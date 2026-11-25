use solidcraft_trivariate::{
    elasticity::assemble,
    row_session::{RowAnalysisSession, RowRequest},
};
fn request(n: i32) -> RowRequest {
    RowRequest {
        half_widths: [10 * n, 20],
        centers: (0..n).map(|i| [20 * i + 10 - 10 * n, 0]).collect(),
        radii: vec![4; n as usize],
        young_mpa: 210000.,
        poisson: 0.3,
        gauss_order: 4,
    }
}
#[test]
fn an_accepted_local_edit_reuses_the_other_cells_and_matches_uncached_assembly() {
    for n in [1, 2, 4, 8, 16] {
        let mut session = RowAnalysisSession::default();
        let mut r = request(n);
        assert_eq!(
            session.update(r.clone()).unwrap().update.integrated,
            4 * n as usize
        );
        r.radii[n as usize / 2] = 5;
        let state = session.update(r.clone()).unwrap();
        assert_eq!(state.update.integrated, 4);
        assert_eq!(state.update.reused, 4 * (n as usize - 1));
        let direct = assemble(&state.space, r.young_mpa, r.poisson, r.gauss_order).unwrap();
        assert_eq!(state.stiffness.columns, direct.columns);
        assert_eq!(state.stiffness.offsets, direct.offsets);
        assert_eq!(state.stiffness.values, direct.values);
        assert!((state.stiffness.integrated_volume / direct.integrated_volume - 1.).abs() < 1e-11);
    }
}
#[test]
fn material_quadrature_and_geometry_faults_preserve_the_accepted_result() {
    let mut session = RowAnalysisSession::default();
    let mut r = request(2);
    session.update(r.clone()).unwrap();
    r.young_mpa *= 2.;
    assert_eq!(session.update(r.clone()).unwrap().update.integrated, 8);
    r.gauss_order = 6;
    assert_eq!(session.update(r.clone()).unwrap().update.integrated, 8);
    assert_eq!(session.update(r.clone()).unwrap().update.reused, 8);
    let saved = session.state().unwrap().stiffness.values.clone();
    for i in 0..4 {
        let mut bad = r.clone();
        match i {
            0 => bad.radii[1] = 0,
            1 => bad.gauss_order = 0,
            2 => bad.poisson = 0.5,
            _ => bad.young_mpa = f64::MAX,
        }
        assert!(session.update(bad).is_err());
        assert_eq!(session.state().unwrap().request, r);
        assert_eq!(session.state().unwrap().stiffness.values, saved);
    }
}
