//! Export verification inputs and both operators for independent mixed-system assembly.
use solidcraft_trivariate::{
    analysis_space::{section_analysis_space, AnalysisBlock, AnalysisSpace},
    elasticity::{assemble, Stiffness},
    json::{object, Json},
    projected_bulk::assemble_projected_bulk,
    sections::section_blocks,
};
use std::collections::BTreeMap;
fn n(x: impl ToString) -> Json {
    Json::Number(x.to_string())
}
fn cube(warp: bool) -> AnalysisSpace {
    let mut controls = Vec::new();
    for i in 0..3_i128 {
        for j in 0..3_i128 {
            for k in 0..3_i128 {
                let mut p = [40 * i, 40 * j, 40 * k, 80];
                if warp {
                    p[0] += if i == 1 { 4 * j + 2 * k } else { 0 };
                    p[1] += if j == 1 { 3 * k - 2 * i } else { 0 };
                    p[2] += i * j * (2 - k);
                    let weight = 5 + (i + j + 2 * k) % 3;
                    p = p.map(|v| v * weight);
                }
                controls.push(p);
            }
        }
    }
    AnalysisSpace {
        ndof: 27,
        blocks: vec![AnalysisBlock {
            id: "cell".into(),
            controls,
            dofs: (0..27).collect(),
        }],
        boundary_dofs: BTreeMap::new(),
    }
}
fn dense(k: &Stiffness) -> Json {
    let mut dense = vec![vec![0.; k.size()]; k.size()];
    for i in 0..k.offsets.len() - 1 {
        for p in k.offsets[i]..k.offsets[i + 1] {
            for a in 0..3 {
                for b in 0..3 {
                    dense[3 * i + a][3 * k.columns[p] + b] = k.values[p][a][b];
                }
            }
        }
    }
    Json::Array(
        dense
            .into_iter()
            .map(|row| Json::Array(row.into_iter().map(n).collect()))
            .collect(),
    )
}
fn run(name: &str, space: &AnalysisSpace, nu: f64) -> Result<Json, &'static str> {
    let full = assemble(space, 210000., nu, 8)?;
    let mixed = assemble_projected_bulk(space, 210000., nu, 8)?;
    let u: Vec<_> = (0..full.size())
        .map(|i| 0.013 * ((2 * i + 1) as f64).sin() + 0.004 * ((13 * i + 3) as f64).cos())
        .collect();
    let mut pressures = BTreeMap::new();
    let mut blocks = Vec::new();
    for b in &space.blocks {
        pressures.insert(
            b.id.clone(),
            Json::Array(mixed.mean_stress_coefficients(&b.id, &u)?.map(n).into()),
        );
        blocks.push(object([
            ("id", Json::String(b.id.clone())),
            (
                "controls",
                Json::Array(
                    b.controls
                        .iter()
                        .map(|p| Json::Array(p.map(n).into()))
                        .collect(),
                ),
            ),
            ("dofs", Json::Array(b.dofs.iter().map(n).collect())),
        ]));
    }
    Ok(object([
        ("name", Json::String(name.into())),
        ("ndof", n(space.ndof)),
        ("young_mpa", n(210000)),
        ("poisson", n(nu)),
        ("order", n(8)),
        ("blocks", Json::Array(blocks)),
        ("projected_stiffness", dense(&mixed.stiffness)),
        ("full_stiffness", dense(&full)),
        (
            "displacement_test",
            Json::Array(u.into_iter().map(n).collect()),
        ),
        ("mean_stress_coefficients", Json::Object(pressures)),
    ]))
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("supply fixture output path")?;
    let ring = section_analysis_space("full_annulus_v1", &section_blocks("full_annulus_v1")?)?;
    let cases = vec![
        run("cube-nu03", &cube(false), 0.3)?,
        run("warped-nu03", &cube(true), 0.3)?,
        run("warped-nu04999", &cube(true), 0.4999)?,
        run("annulus-nu03", &ring, 0.3)?,
    ];
    std::fs::write(
        path,
        object([
            (
                "schema",
                Json::String("projected_bulk_verification_v1".into()),
            ),
            ("cases", Json::Array(cases)),
        ])
        .to_string(),
    )?;
    Ok(())
}
