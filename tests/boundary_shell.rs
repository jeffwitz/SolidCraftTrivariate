use solidcraft_trivariate::{
    boundary_shell::boundary_shell, interchange::build, json::parse, step_boundary::export_step,
};
use std::{fs, process::Command};

#[test]
fn euler_and_step_shell_inventory_agree_with_the_analytic_topology() {
    for (family, euler) in [
        ("reference_rectangle_v1", 2),
        ("core_disk_v1", 2),
        ("full_annulus_v1", 0),
        ("pierced_square_v1", 0),
        ("rounded_corner_v1", 2),
        ("layered_pocket_v1", 2),
        ("filleted_bracket_v1", -2),
        ("overlap_bore_v1", 0),
    ] {
        let recipe = parse(&format!(r#"{{"family":"{family}"}}"#)).unwrap();
        let model = build(&recipe).unwrap();
        let shell = boundary_shell(&model).unwrap();
        let count =
            shell.vertices.len() as i64 - shell.edges.len() as i64 + shell.faces.len() as i64;
        assert_eq!(count, euler);
        let step = export_step(&model).unwrap();
        assert_eq!(step.matches("MANIFOLD_SOLID_BREP(").count(), 1);
        assert_eq!(step.matches(" = EDGE_CURVE(").count(), shell.edges.len());
        assert_eq!(step.matches(" = ADVANCED_FACE(").count(), shell.faces.len());
    }
}

#[test]
fn rejected_export_keeps_the_existing_step_bytes() {
    let dir = std::env::temp_dir().join(format!("solidcraft-step-refusal-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let input = dir.join("invalid.json");
    let output = dir.join("retained.stp");
    fs::write(&input, br#"{"schema":"pierced_square_v1","blocks":[]}"#).unwrap();
    fs::write(&output, b"prior accepted STEP bytes").unwrap();
    let res = Command::new(env!("CARGO_BIN_EXE_solidcraft-trivariate"))
        .arg("export-step")
        .arg(&input)
        .arg(&output)
        .output()
        .unwrap();
    assert!(!res.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"prior accepted STEP bytes");
    fs::remove_file(input).unwrap();
    fs::remove_file(output).unwrap();
    fs::remove_dir(dir).unwrap();
}
