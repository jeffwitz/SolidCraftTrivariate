use solidcraft_trivariate::{
    json::{object, Json},
    reentrant_ritz,
    thb_reentrant_edge::{build_hierarchy, build_hierarchy_with_widths, build_uniform},
};
use std::{fs, path::Path};
#[test]
#[ignore = "release multilevel scalar corner convergence and quadrature sequence"]
fn singular_error_per_dof() {
    let dir = Path::new("build/reentrant-ritz");
    fs::create_dir_all(dir).unwrap();
    let mut local = Vec::new();
    let mut uniform = Vec::new();
    for level in 0..=3 {
        for (name, curve) in [("thb", &mut local), ("uniform", &mut uniform)] {
            let h = if name == "thb" {
                build_hierarchy(level)
            } else {
                build_uniform(level)
            }
            .unwrap();
            let row = reentrant_ritz::solve(&h).unwrap();
            let probes = (0..3)
                .flat_map(|patch| {
                    [[0.013, 0.029], [0.17, 0.19], [0.43, 0.62], [0.87, 0.91]]
                        .map(move |uv| (patch, uv))
                })
                .map(|(patch, uv)| {
                    let n = |x: f64| Json::Number(x.to_string());
                    let values = h
                        .sample_patch(patch, level, uv)
                        .unwrap()
                        .into_iter()
                        .filter(|(_, v, _)| *v != 0.)
                        .map(|(_, v, g)| Json::Array(vec![n(v), n(g[0]), n(g[1])]))
                        .collect();
                    object([
                        ("patch", Json::Number(patch.to_string())),
                        ("uv", Json::Array(uv.map(n).to_vec())),
                        ("values_and_xy_gradients", Json::Array(values)),
                    ])
                })
                .collect();
            let value = object([
                ("curve", Json::String(name.into())),
                ("solve", row.clone()),
                ("basis_probes", Json::Array(probes)),
            ]);
            fs::write(dir.join(format!("{name}-{level}.json")), value.to_string()).unwrap();
            let get = |k| row.field(k).unwrap().number::<f64>().unwrap();
            curve.push((get("active_scalar_dofs"), get("relative_H1_seminorm")));
        }
    }
    // Freeze the efficiency verdict separately from basis validity.
    let efficient = local
        .iter()
        .any(|&(n, e)| uniform.iter().any(|&(m, f)| m >= 1.2 * n && e < f));
    fs::write(
        dir.join("efficiency.json"),
        Json::Bool(efficient).to_string(),
    )
    .unwrap();
}

#[test]
#[ignore = "supplementary wider mark on the same corner problem"]
fn wider_mark() {
    let dir = Path::new("build/reentrant-ritz");
    fs::create_dir_all(dir).unwrap();
    for level in 1..=3 {
        let space = build_hierarchy_with_widths(level, &vec![6; level]).unwrap();
        let row = reentrant_ritz::solve(&space).unwrap();
        fs::write(dir.join(format!("wide-{level}.json")), row.to_string()).unwrap();
    }
}
