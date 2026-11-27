use solidcraft_trivariate::{
    cylinder_tensor_covers::build_cover,
    json::{parse, Json},
};
use std::{env, fs};
fn sha256(path: &str) -> Result<String, String> {
    let output = std::process::Command::new("sha256sum")
        .arg(path)
        .output()
        .map_err(|e| format!("sha256sum: {e}"))?;
    if !output.status.success() {
        return Err("sha256sum failed".into());
    }
    String::from_utf8(output.stdout)
        .map_err(|e| e.to_string())
        .map(|x| x.split_whitespace().next().unwrap().to_owned())
}

fn run() -> Result<(), String> {
    let args: Vec<_> = env::args().collect();
    if args.len() != 3 {
        return Err("usage: cylinder_tensor_covers MODEL OUTDIR".into());
    }
    let model_text = fs::read_to_string(&args[1]).map_err(|e| e.to_string())?;
    let model = parse(&model_text).map_err(str::to_owned)?;
    let model_sha = sha256(&args[1])?;
    let recipe = std::path::Path::new(&args[1]).with_file_name("cylinder.recipe.json");
    let recipe_sha = sha256(recipe.to_str().ok_or("invalid recipe path")?)?;
    for count in [2, 3] {
        let mut result = build_cover(&model, count)?;
        if let Json::Object(fields) = &mut result {
            fields.insert(
                "source_model_sha256".into(),
                Json::String(model_sha.clone()),
            );
            fields.insert(
                "source_recipe_sha256".into(),
                Json::String(recipe_sha.clone()),
            );
        }
        let file = format!("{}/cylinder-{count}.patches.json", args[2]);
        fs::write(&file, format!("{result}\n")).map_err(|e| e.to_string())?;
        println!("wrote {file}");
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        std::process::exit(1)
    }
}
