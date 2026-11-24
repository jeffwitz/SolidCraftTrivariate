use solidcraft_trivariate::json::parse;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

static NEXT_DIR: AtomicUsize = AtomicUsize::new(0);

struct TempDir(PathBuf);

impl TempDir {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "solidcraft-trivariate-cli-{}-{}",
            std::process::id(),
            NEXT_DIR.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        Self(path)
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn cli(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
        .args(args)
        .output()
        .expect("launch CLI")
}

fn status(output: &Output) -> String {
    let text = std::str::from_utf8(&output.stdout).expect("CLI JSON is UTF-8");
    let value = parse(text).expect("CLI emits JSON");
    value.field("status").unwrap().text().unwrap().to_owned()
}

fn write(path: &Path, text: &str) {
    fs::write(path, text).unwrap();
}

#[test]
fn historical_affine_and_annular_recipes_build_and_audit() {
    for (name, recipe) in [
        (
            "affine",
            include_str!("references/affine_pair.json"),
        ),
        (
            "annular",
            include_str!("references/annular_pair.json"),
        ),
    ] {
        let dir = TempDir::new();
        let recipe_path = dir.path().join(format!("{name}-recipe.json"));
        let model_path = dir.path().join(format!("{name}-model.json"));
        write(&recipe_path, recipe);
        let recipe_path = recipe_path.to_str().unwrap();
        let model_path_str = model_path.to_str().unwrap();
        let built = cli(&["build", recipe_path, model_path_str]);
        assert!(built.status.success(), "{}", String::from_utf8_lossy(&built.stdout));
        assert_eq!(status(&built), "valid");
        let audited = cli(&["audit", model_path_str]);
        assert!(audited.status.success(), "{}", String::from_utf8_lossy(&audited.stdout));
        assert_eq!(status(&audited), "valid");
    }
}

#[test]
fn affine_decimal_integer_widths_are_not_limited_to_machine_numbers() {
    let dir = TempDir::new();
    let recipe = dir.path().join("huge.json");
    let model = dir.path().join("huge-model.json");
    write(
        &recipe,
        r#"{"family":"affine_pair_v1","widths_mm":[123456789012345678901234567890,987654321098765432109876543210]}"#,
    );
    let built = cli(&["build", recipe.to_str().unwrap(), model.to_str().unwrap()]);
    assert!(built.status.success(), "{}", String::from_utf8_lossy(&built.stdout));
    assert_eq!(status(&built), "valid");
    let saved = fs::read_to_string(&model).unwrap();
    assert!(saved.contains("1111111110111111111011111111100"), "{saved}");
    let audited = cli(&["audit", model.to_str().unwrap()]);
    assert!(audited.status.success());
    assert_eq!(status(&audited), "valid");
}

#[test]
fn refusal_preserves_existing_output_without_leaving_temporary_files() {
    let dir = TempDir::new();
    let recipe = dir.path().join("invalid.json");
    let output = dir.path().join("existing.json");
    write(&output, "keep this file\n");
    write(&recipe, r#"{"family":"annular_pair_v1","heights_mm":[1,0]}"#);
    let result = cli(&["build", recipe.to_str().unwrap(), output.to_str().unwrap()]);
    assert!(!result.status.success());
    assert_eq!(fs::read_to_string(&output).unwrap(), "keep this file\n");
    let files: Vec<_> = fs::read_dir(dir.path()).unwrap().collect();
    assert_eq!(files.len(), 2);
}

#[test]
fn concurrent_valid_builds_publish_one_complete_model_atomically() {
    let dir = TempDir::new();
    let output = dir.path().join("shared-model.json");
    let recipes = [
        r#"{"family":"affine_pair_v1","widths_mm":[2,3]}"#,
        r#"{"family":"affine_pair_v1","widths_mm":[5,7]}"#,
        r#"{"family":"annular_pair_v1","heights_mm":[1,2]}"#,
        r#"{"family":"annular_pair_v1","heights_mm":[3,4]}"#,
    ];
    let handles: Vec<_> = recipes
        .iter()
        .enumerate()
        .map(|(index, source)| {
            let recipe = dir.path().join(format!("recipe-{index}.json"));
            write(&recipe, source);
            let recipe = recipe.to_str().unwrap().to_owned();
            let output = output.to_str().unwrap().to_owned();
            thread::spawn(move || {
                let result = cli(&["build", &recipe, &output]);
                assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stdout));
            })
        })
        .collect();
    for handle in handles {
        handle.join().unwrap();
    }
    let audited = cli(&["audit", output.to_str().unwrap()]);
    assert!(audited.status.success(), "{}", String::from_utf8_lossy(&audited.stdout));
    assert_eq!(status(&audited), "valid");
    let remaining: Vec<_> = fs::read_dir(dir.path()).unwrap().collect();
    assert_eq!(remaining.len(), 5, "no temporary publication files remain");
}

#[test]
fn invalid_inputs_and_unknown_schema_have_distinct_cli_statuses() {
    let dir = TempDir::new();
    let recipe = dir.path().join("recipe.json");
    let model = dir.path().join("unknown.json");

    for invalid in [
        r#"{"family":"annular_pair_v1","heights_mm":[0,2]}"#,
        r#"{"family":"affine_pair_v1","widths_mm":[0,2]}"#,
        r#"{"family":"future_family"}"#,
        "{malformed",
    ] {
        write(&recipe, invalid);
        let result = cli(&["build", recipe.to_str().unwrap(), model.to_str().unwrap()]);
        assert!(!result.status.success());
        assert_eq!(status(&result), "error");
    }

    write(&model, r#"{"schema":"future_family"}"#);
    let result = cli(&["audit", model.to_str().unwrap()]);
    assert!(!result.status.success());
    assert_eq!(status(&result), "indeterminate");
}
