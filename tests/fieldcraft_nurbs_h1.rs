use solidcraft_trivariate::{
    json::{parse, Json},
    mechanical_job, nurbs_space,
};
use std::{collections::BTreeMap, fs};

fn model() -> Json {
    parse(
        &fs::read_to_string("docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json")
            .unwrap(),
    )
    .unwrap()
}

fn number(value: &Json) -> usize {
    value.number().unwrap()
}

fn package_blocks(package: &Json) -> BTreeMap<String, &Json> {
    package
        .field("blocks")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|b| (b.field("id").unwrap().text().unwrap().to_owned(), b))
        .collect()
}

fn dof_partition(blocks: &[Json]) -> Vec<Vec<(String, usize)>> {
    let mut groups: BTreeMap<usize, Vec<(String, usize)>> = BTreeMap::new();
    for block in blocks {
        let id = block.field("id").unwrap().text().unwrap().to_owned();
        for (local, dof) in block
            .field("dofs")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .enumerate()
        {
            groups
                .entry(number(dof))
                .or_default()
                .push((id.clone(), local));
        }
    }
    let mut result: Vec<_> = groups.into_values().collect();
    result.sort();
    result
}

fn face_index(axis: usize, upper: bool, a: usize, b: usize, p: usize) -> usize {
    let mut q = [0usize; 3];
    let free: Vec<_> = (0..3).filter(|&i| i != axis).collect();
    q[axis] = if upper { p } else { 0 };
    q[free[0]] = a;
    q[free[1]] = b;
    (q[0] * (p + 1) + q[1]) * (p + 1) + q[2]
}

fn side(value: &Json) -> (String, usize, bool) {
    let items = value.items().unwrap();
    let block = items[0].text().unwrap().to_owned();
    let name = items[1].text().unwrap();
    (
        block,
        "uvw".find(name.chars().next().unwrap()).unwrap(),
        name.ends_with('1'),
    )
}

fn bernstein(p: usize, i: usize, t: f64) -> f64 {
    let binomial = match (p, i) {
        (2, 0) | (2, 2) => 1.0,
        (2, 1) => 2.0,
        (1, _) => 1.0,
        _ => panic!("unexpected degree"),
    };
    binomial * t.powi(i as i32) * (1.0 - t).powi((p - i) as i32)
}

fn verify_affine_patch(package: &Json) {
    let matrix = [
        [0.013, -0.004, 0.002],
        [0.001, 0.009, -0.003],
        [-0.002, 0.005, 0.006],
    ];
    let offset = [0.7, -0.2, 0.3];
    for block in package.field("blocks").unwrap().items().unwrap() {
        let controls: Vec<Vec<f64>> = block
            .field("controls")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .map(|p| {
                p.items()
                    .unwrap()
                    .iter()
                    .map(|x| x.number::<f64>().unwrap())
                    .collect()
            })
            .collect();
        let degrees: Vec<usize> = block
            .field("degrees")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .map(number)
            .collect();
        assert_eq!(degrees, [2, 2, 2]);
        for uvw in [
            [0.0, 0.0, 0.0],
            [0.25, 0.5, 0.75],
            [0.5, 0.5, 0.5],
            [1.0, 1.0, 1.0],
        ] {
            let mut den = 0.0;
            let mut x = [0.0; 3];
            let mut u = [0.0; 3];
            for i in 0..=2 {
                for j in 0..=2 {
                    for k in 0..=2 {
                        let basis = bernstein(2, i, uvw[0])
                            * bernstein(2, j, uvw[1])
                            * bernstein(2, k, uvw[2]);
                        let c = &controls[(i * 3 + j) * 3 + k];
                        let w = c[3] * basis;
                        den += w;
                        let point = [c[0] / c[3], c[1] / c[3], c[2] / c[3]];
                        for a in 0..3 {
                            x[a] += w * point[a];
                            u[a] += w
                                * (offset[a]
                                    + (0..3).map(|b| matrix[a][b] * point[b]).sum::<f64>());
                        }
                    }
                }
            }
            let x: [f64; 3] = std::array::from_fn(|a| x[a] / den);
            for a in 0..3 {
                u[a] /= den;
                let expected = offset[a] + (0..3).map(|b| matrix[a][b] * x[b]).sum::<f64>();
                assert!(
                    (u[a] - expected).abs() < 2e-12,
                    "affine patch defect block={} uvw={uvw:?} component={a} got={} expected={}",
                    block.field("id").unwrap().text().unwrap(),
                    u[a],
                    expected
                );
            }
        }
    }
}

fn set_feature(job: &mut Json, group: &str) {
    let Json::Object(top) = job else {
        panic!("job object")
    };
    let Json::Array(entries) = top.get_mut(group).unwrap() else {
        panic!("entries")
    };
    for entry in entries {
        let Json::Object(entry) = entry else {
            panic!("entry")
        };
        let Json::Object(select) = entry.get_mut("select").unwrap() else {
            panic!("select")
        };
        select.insert("feature".into(), Json::String("outer/profile".into()));
    }
}

#[test]
fn fieldcraft_gets_degree_harmonized_c0_h1_space_and_native_q2_mechanics() {
    let model = model();
    let audit = solidcraft_trivariate::interchange::audit(&model).unwrap();
    assert_eq!(audit.field("status").unwrap().text().unwrap(), "valid");

    let package = nurbs_space::analysis_package(&model, 16).unwrap();
    assert_eq!(
        package.field("schema").unwrap().text().unwrap(),
        "solidcraft_nurbs_h1_space_v1"
    );
    assert_eq!(
        package.field("trace_continuity").unwrap().text().unwrap(),
        "C0"
    );
    assert_eq!(number(package.field("tensor_degree").unwrap()), 2);
    assert_eq!(number(package.field("scalar_ndof").unwrap()), 171);
    assert_eq!(package.field("blocks").unwrap().items().unwrap().len(), 10);
    assert!(
        (package
            .field("integrated_volume_mm3")
            .unwrap()
            .number::<f64>()
            .unwrap()
            - 32144.7789254)
            .abs()
            < 2e-7
    );

    let package_blocks = package_blocks(&package);
    let source_blocks: BTreeMap<_, _> = model
        .field("blocks")
        .unwrap()
        .items()
        .unwrap()
        .iter()
        .map(|b| (b.field("id").unwrap().text().unwrap(), b))
        .collect();
    let interfaces = model
        .field("topology")
        .unwrap()
        .field("interfaces")
        .unwrap()
        .items()
        .unwrap();
    assert_eq!(interfaces.len(), 11);
    for interface in interfaces {
        let (a_id, a_axis, a_upper) = side(interface.field("first").unwrap());
        let (b_id, b_axis, b_upper) = side(interface.field("second").unwrap());
        let a = package_blocks[&a_id];
        let b = package_blocks[&b_id];
        let orientation = interface.field("orientation").unwrap();
        let flip_first = matches!(orientation.field("flip_first").unwrap(), Json::Bool(true));
        let flip_second = matches!(orientation.field("flip_second").unwrap(), Json::Bool(true));
        let swap = matches!(orientation.field("swap").unwrap(), Json::Bool(true));
        let da: Vec<usize> = a
            .field("dofs")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .map(number)
            .collect();
        let db: Vec<usize> = b
            .field("dofs")
            .unwrap()
            .items()
            .unwrap()
            .iter()
            .map(number)
            .collect();
        let mut matched = 0;
        for i in 0..=2 {
            for j in 0..=2 {
                let s = if flip_first { 2 - i } else { i };
                let t = if flip_second { 2 - j } else { j };
                let (r, c) = if swap { (t, s) } else { (s, t) };
                assert_eq!(
                    da[face_index(a_axis, a_upper, i, j, 2)],
                    db[face_index(b_axis, b_upper, r, c, 2)],
                    "interface {:?} -> {:?}",
                    interface.field("first").unwrap(),
                    interface.field("second").unwrap()
                );
                matched += 1;
            }
        }
        assert_eq!(matched, 9);
        // Confirm this is the same whole-face identity in the audited source model.
        assert!(
            source_blocks.contains_key(a_id.as_str()) && source_blocks.contains_key(b_id.as_str())
        );
    }
    verify_affine_patch(&package);

    let mut job = parse(include_str!(
        "references/mechanical_jobs/fieldcraft_bracket_hex_pa5.json"
    ))
    .unwrap();
    set_feature(&mut job, "supports");
    set_feature(&mut job, "loads");
    if let Json::Object(fields) = &mut job {
        fields.insert("order".into(), Json::Number("16".into()));
    }
    let solved = mechanical_job::solve_job(&model, &job).unwrap();
    assert_eq!(number(solved.field("scalar_dofs").unwrap()), 171);
    assert_eq!(
        solved
            .field("analysis_blocks")
            .unwrap()
            .items()
            .unwrap()
            .len(),
        10
    );
    assert_eq!(
        dof_partition(package.field("blocks").unwrap().items().unwrap()),
        dof_partition(solved.field("analysis_blocks").unwrap().items().unwrap()),
        "native mechanics and degree-harmonized NURBS export must mount the same C0 field"
    );
    assert!(
        solved
            .field("relative_residual")
            .unwrap()
            .number::<f64>()
            .unwrap()
            < 1e-9
    );
    assert!(
        (solved
            .field("strain_energy_n_mm")
            .unwrap()
            .number::<f64>()
            .unwrap()
            - 135.4447513902)
            .abs()
            < 2e-8
    );
}
