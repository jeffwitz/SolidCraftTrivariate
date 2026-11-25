use solidcraft_trivariate::bracket_mechanics_session::{
    cache_preflight, BracketMechanicsConfig, BracketMechanicsSession,
};
use std::sync::Arc;

#[test]
fn edit_reject_return_and_single_axis_change_are_transactional() {
    let config = BracketMechanicsConfig {
        uniform_level: 0,
        ..BracketMechanicsConfig::campaign(0)
    };
    let mut session = BracketMechanicsSession::new(config).unwrap();
    session.update([50, 74]).unwrap();
    let base = Arc::clone(session.state().unwrap());
    let edited = session.update([51, 74]).unwrap();
    assert_eq!(edited.radii_sixteenths_mm, [51, 74]);
    assert!(edited.timing.reused_contributions > 0);
    assert!(edited.timing.stiffness_action_relative_error <= 1e-11);
    assert!(edited.timing.warm_cold_displacement_relative_error.unwrap() <= 1e-8);

    let before_rejection = Arc::clone(session.state().unwrap());
    let cache_before = session.cache_bytes_estimated();
    assert!(session.update([81, 74]).is_err());
    assert!(Arc::ptr_eq(session.state().unwrap(), &before_rejection));
    assert_eq!(session.cache_bytes_estimated(), cache_before);

    let returned = session.update([50, 74]).unwrap();
    assert_eq!(returned.radii_sixteenths_mm, [50, 74]);
    assert!(returned.timing.reused_contributions >= 36);
    assert_eq!(returned.blocks, base.blocks);
    let y_only = session.update([50, 73]).unwrap();
    assert_eq!(y_only.radii_sixteenths_mm, [50, 73]);
    assert!(y_only.timing.reused_contributions > 0);
    assert!(y_only.timing.integrated_contributions > 0);
}

#[test]
fn level_three_is_refused_by_cache_preflight() {
    let estimate = cache_preflight(3, 384 * 1024 * 1024).unwrap();
    assert_eq!(estimate.cells, 18_432);
    assert!(!estimate.fits);
}
