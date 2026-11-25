use solidcraft_trivariate::{
    analysis_space::{AnalysisBlock, AnalysisSpace},
    boundary_traction::{integrate_traction, LoadPatch},
    incidence::BlockSide,
};
use std::collections::{BTreeMap, BTreeSet};

fn box_space() -> (AnalysisSpace, BTreeSet<BlockSide>) {
    let mut block = AnalysisBlock {
        id: "box".into(),
        controls: Vec::new(),
        dofs: (0..27).collect(),
    };
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                block.controls.push([50 + 5 * i, 4 * j, 20 * k, 1]);
            }
        }
    }
    let sides = (0..3)
        .flat_map(|axis| {
            [false, true].map(move |upper| BlockSide {
                block: "box".into(),
                axis,
                upper,
            })
        })
        .collect();
    (
        AnalysisSpace {
            blocks: vec![block],
            ndof: 27,
            boundary_dofs: BTreeMap::new(),
        },
        sides,
    )
}
fn face() -> BlockSide {
    BlockSide {
        block: "box".into(),
        axis: 0,
        upper: true,
    }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-10 * b.abs().max(1.), "{a} != {b}");
}

#[test]
fn integrates_affine_traction_resultant_and_moment_on_a_half_and_a_quarter() {
    let (s, exterior) = box_space();
    for (yz, y, z) in [([[0, 1], [0, 2]], 4., 40.), ([[0, 1], [0, 1]], 4., 20.)] {
        let p = LoadPatch {
            side: face(),
            depth: 1,
            intervals: yz,
        };
        let r = integrate_traction(&s, &exterior, &[p], 4, |[_, y, z], normal| {
            close(normal[0], 1.);
            close(normal[1], 0.);
            close(normal[2], 0.);
            [0., -3. - y / 8. - z / 40., 0.]
        })
        .unwrap();
        let area = y * z;
        let fy = -area * (3. + y / 16. + z / 80.);
        let mx = area * ((3. + y / 16.) * z / 2. + z * z / 120.);
        close(r.area_mm2, area);
        close(r.resultant_n[1], fy);
        close(r.moment_n_mm[0], mx);
        close(r.moment_n_mm[2], 60. * fy);
        close(r.nodal_force_n.iter().skip(1).step_by(3).sum(), fy);
        let nodal_mx = s.blocks[0]
            .controls
            .iter()
            .zip(r.nodal_force_n.as_chunks::<3>().0)
            .map(|(h, f)| -(h[2] as f64) * f[1])
            .sum();
        close(nodal_mx, mx);
    }
}

#[test]
fn touching_quadrants_reproduce_full_load_and_oriented_closed_surface_flux() {
    let (s, exterior) = box_space();
    let q: Vec<_> = (0..2)
        .flat_map(|i| {
            (0..2).map(move |j| LoadPatch {
                side: face(),
                depth: 1,
                intervals: [[i, i + 1], [j, j + 1]],
            })
        })
        .collect();
    let traction = |[_, y, z]: [f64; 3], _| [y / 2., -3. - z / 2., 1.];
    let part = integrate_traction(&s, &exterior, &q, 4, traction).unwrap();
    let full = integrate_traction(&s, &exterior, &[LoadPatch::whole(face())], 4, traction).unwrap();
    for (a, b) in part.nodal_force_n.iter().zip(full.nodal_force_n) {
        close(*a, b);
    }
    let all = exterior
        .iter()
        .cloned()
        .map(LoadPatch::whole)
        .collect::<Vec<_>>();
    let flux = integrate_traction(&s, &exterior, &all, 4, |_, n| n).unwrap();
    for val in flux.resultant_n.into_iter().chain(flux.moment_n_mm) {
        close(val, 0.);
    }
    close(flux.area_mm2, 2. * (10. * 8. + 10. * 40. + 8. * 40.));
}

#[test]
fn rejects_nonexterior_overlapping_empty_invalid_and_nonfinite_loads() {
    let (s, exterior) = box_space();
    let p = LoadPatch::whole(face());
    for (d, iv) in [
        (0, [[0, 0], [0, 1]]),
        (1, [[0, 3], [0, 2]]),
        (31, [[0, 1]; 2]),
    ] {
        let bad = LoadPatch {
            side: face(),
            depth: d,
            intervals: iv,
        };
        assert!(integrate_traction(&s, &exterior, &[bad], 4, |_, _| [0.; 3]).is_err());
    }
    assert!(integrate_traction(
        &s,
        &BTreeSet::new(),
        std::slice::from_ref(&p),
        4,
        |_, _| [0.; 3]
    )
    .is_err());
    assert!(integrate_traction(&s, &exterior, &[p.clone(), p.clone()], 4, |_, _| [0.; 3]).is_err());
    assert!(integrate_traction(&s, &exterior, &[], 4, |_, _| [0.; 3]).is_err());
    assert!(
        integrate_traction(&s, &exterior, std::slice::from_ref(&p), 0, |_, _| [0.; 3]).is_err()
    );
    assert!(integrate_traction(&s, &exterior, &[p], 4, |_, _| [f64::NAN; 3]).is_err());
}
