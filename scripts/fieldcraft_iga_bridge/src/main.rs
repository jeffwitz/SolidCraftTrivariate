use fieldcraft_io::{
    results::{self, SolvedRef, StoredRef},
    Extras,
};
use fieldcraft_solver::{StaticResult, Stats, StepResult};
use serde_json::Value;
use std::{collections::HashMap, env, fs, path::Path};

fn run(sample_path: &Path, base_path: &Path, output_path: &Path, step: &str) -> Result<(), String> {
    let sample: Value = serde_json::from_slice(&fs::read(sample_path).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    let nodes = sample["nodes"].as_array().ok_or("sample has no nodes")?;
    let base = fs::read(base_path).map_err(|e| e.to_string())?;
    let file = fieldcraft_io::read(&base).map_err(|e| e.to_string())?;

    let mut displacement = Vec::with_capacity(file.model.nodes().len());
    let mut stress = Vec::with_capacity(file.model.nodes().len());
    let mut von_mises = Vec::with_capacity(file.model.nodes().len());
    let labels: Vec<u32> = file.model.nodes().keys().copied().collect();
    if labels.len() != nodes.len() {
        return Err(format!(
            "FieldCraft archive has {} nodes but IGA sample has {}",
            labels.len(),
            nodes.len()
        ));
    }
    let mut by_label = HashMap::with_capacity(nodes.len());
    for node in nodes {
        let label = u32::try_from(node["label"].as_u64().ok_or("invalid IGA node label")?)
            .map_err(|_| "IGA node label is out of range")?;
        if by_label.insert(label, node).is_some() {
            return Err(format!("IGA sample repeats node {label}"));
        }
    }
    for label in &labels {
        let value = by_label
            .get(label)
            .ok_or_else(|| format!("IGA sample is missing node {label}"))?;
        let xyz = value["xyz_m"].as_array().ok_or("invalid xyz_m")?;
        let u = value["displacement_m"]
            .as_array()
            .ok_or("invalid displacement_m")?;
        let s = value["stress_pa"].as_array().ok_or("invalid stress_pa")?;
        if xyz.len() != 3 || u.len() != 3 || s.len() != 6 {
            return Err(format!("IGA node {label} has an invalid vector dimension"));
        }
        let to_finite = |v: &Value| -> Result<f64, String> {
            v.as_f64()
                .filter(|x| x.is_finite())
                .ok_or_else(|| format!("IGA node {label} has a nonfinite value"))
        };
        let xyz: [f64; 3] = xyz
            .iter()
            .map(to_finite)
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .map_err(|_| "invalid xyz_m")?;
        let u: [f64; 3] = u
            .iter()
            .map(to_finite)
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .map_err(|_| "invalid displacement_m")?;
        let s: [f64; 6] = s
            .iter()
            .map(to_finite)
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .map_err(|_| "invalid stress_pa")?;
        let vm = to_finite(&value["von_mises_pa"])?;
        displacement.push(u);
        stress.push(s);
        von_mises.push(vm);
        let point = &file.model.nodes()[label].xyz;
        for i in 0..3 {
            let expected = xyz[i];
            if (point[i] - expected).abs() > 2e-11 {
                return Err(format!(
                    "node {label} does not match the imported IGA sample"
                ));
            }
        }
    }

    let static_result = StaticResult {
        step: step.to_owned(),
        nodes: labels,
        displacement,
        rotation: Vec::new(),
        reaction: Vec::new(),
        reaction_moment: Vec::new(),
        stress,
        von_mises,
        elements: Vec::new(),
        shells: Vec::new(),
        beams: Vec::new(),
        reference_points: Vec::new(),
        applied_total: vector3(&sample, "applied_total_n")?,
        reaction_total: vector3(&sample, "reaction_total_n")?,
        stats: Stats::default(),
        warnings: vec!["This is a SolidCraft Trivariate IGA result sampled onto a linear display mesh; FieldCraft did not solve it. Total reactions are retained, but native per-node reactions are unavailable and omitted.".to_owned()],
        temperature: None,
        nonlinear: None,
        contacts: Vec::new(),
    };
    let solved = StepResult::Static(static_result);
    let definition = results::definition_hash(&file.model, step, &file.geometry)
        .ok_or_else(|| format!("FieldCraft model has no step named `{step}`"))?;
    let extras = Extras {
        results: vec![(
            step,
            StoredRef {
                definition,
                result: SolvedRef::Implicit(&solved),
            },
        )],
        ..Extras::default()
    };
    let bytes = fieldcraft_io::write_with(
        &file.model,
        &file.geometry,
        &file.journal,
        &file.journal_header,
        &file.view,
        &extras,
    )
    .map_err(|e| e.to_string())?;
    fs::write(output_path, bytes).map_err(|e| e.to_string())?;
    // Re-open through the production reader to validate ZIP, hashes, mesh and result codec.
    let check = fieldcraft_io::read(&fs::read(output_path).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    if check.results.len() != 1 || results::check(&check.results[0].result).is_err() {
        return Err("written FieldCraft archive did not pass its result reader".into());
    }
    println!("{{\"status\":\"written\",\"nodes\":{},\"elements\":{},\"result\":\"SolidCraft Trivariate IGA sampled field\"}}", check.model.nodes().len(), check.model.elements().len());
    Ok(())
}

fn vector3(sample: &Value, name: &str) -> Result<[f64; 3], String> {
    let values = sample[name]
        .as_array()
        .ok_or_else(|| format!("sample has no {name}"))?;
    if values.len() != 3 {
        return Err(format!("sample {name} is not a vector3"));
    }
    values
        .iter()
        .map(|v| {
            v.as_f64()
                .filter(|x| x.is_finite())
                .ok_or_else(|| format!("sample {name} has a nonfinite value"))
        })
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| format!("sample {name} is not a vector3"))
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!(
            "usage: fieldcraft-iga-archive-bridge SAMPLE.json BASE.fieldcraft OUT.fieldcraft STEP"
        );
        std::process::exit(2);
    }
    if let Err(e) = run(
        Path::new(&args[1]),
        Path::new(&args[2]),
        Path::new(&args[3]),
        &args[4],
    ) {
        eprintln!("fieldcraft IGA archive bridge: {e}");
        std::process::exit(1);
    }
}
