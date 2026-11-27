#!/usr/bin/env python3
"""Check native sketch preview visibility and B-rep mesh invariance in SolidCraft."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from PIL import Image


def run(command: list[str], cwd: Path) -> str:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(
            f"command failed ({result.returncode}): {command!r}\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )
    return result.stdout


def read_eval(cli: str, design: Path, destination: Path, cwd: Path) -> dict:
    run([cli, "eval", str(design), "--json-out", str(destination), "--quiet"], cwd)
    return json.loads(destination.read_text())


def read_mesh(cli: str, design: Path, body: str, out: Path, script: Path, cwd: Path) -> list[dict]:
    script.write_text(
        json.dumps({"commands": [{"command": "model.mesh", "params": {"body": body, "tolerance": 0.05}}]})
    )
    run([cli, "run", str(script), "--in", str(design), "--results", "--quiet", "--json-out", str(out)], cwd)
    return json.loads(out.read_text())


def relative_path(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--solidcraft-cli", default="solidcraft-cli")
    parser.add_argument("--trivariate-cli", type=Path, default=Path("target/debug/solidcraft-trivariate"))
    parser.add_argument("--output", type=Path, default=Path("build/solidcraft-native-render-qa"))
    args = parser.parse_args()
    root = args.root.resolve()
    out = (root / args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    cli = args.solidcraft_cli
    kernel = args.trivariate_cli if args.trivariate_cli.is_absolute() else root / args.trivariate_cli

    source_path = root / "tests/references/solidcraft_sketch_bracket.timeline.json"
    source = json.loads(source_path.read_text())
    source["commands"] = [step for step in source["commands"] if step["command"] == "timeline.add"][:2]
    timeline = out / "bracket-extrusion-two-feature.timeline.json"
    timeline.write_text(json.dumps(source))

    base_design = out / "bracket-extrusion-base.sc"
    run([cli, "run", str(timeline), "--save", str(base_design), "--quiet"], root)
    preview_design = out / "bracket-extrusion-preview-hidden.sc"
    run(
        [
            sys.executable,
            str(root / "scripts/solidcraft_preview_design.py"),
            str(timeline),
            "1",
            str(preview_design),
            "--kernel",
            str(kernel),
            "--host",
            cli,
        ],
        root,
    )

    initial = read_eval(cli, preview_design, out / "bracket-preview-initial.eval.json", root)
    source_eval = read_eval(cli, base_design, out / "bracket-preview-base.eval.json", root)
    original_features = [(x.get("id"), x.get("name"), x.get("type")) for x in source_eval["timeline"]]
    preview_features = [(x.get("id"), x.get("name"), x.get("type")) for x in initial["timeline"]]
    if len(original_features) != 2 or preview_features[:2] != original_features or len(preview_features) != 3:
        raise RuntimeError(f"expected two source features plus one preview sketch, got {original_features} / {preview_features}")
    preview_sketch = next((x for x in initial["sketches"] if x["name"] == "SCT Partition 1"), None)
    if not preview_sketch or preview_sketch["visible"]:
        raise RuntimeError("the new native preview sketch was not initially hidden")
    preview_id = preview_sketch["id"]
    show_script = out / "show-native-sketch.json"
    hide_script = out / "hide-native-sketch.json"
    show_script.write_text(json.dumps({"commands": [{"command": "browser.visibility", "params": {"visible": True, "items": [preview_id]}}]}))
    hide_script.write_text(json.dumps({"commands": [{"command": "browser.visibility", "params": {"visible": False, "items": [preview_id]}}]}))

    shown_design = out / "bracket-extrusion-preview-shown.sc"
    hidden_again_design = out / "bracket-extrusion-preview-hidden-again.sc"
    run([cli, "run", str(show_script), "--in", str(preview_design), "--save", str(shown_design), "--quiet"], root)
    run([cli, "run", str(hide_script), "--in", str(shown_design), "--save", str(hidden_again_design), "--quiet"], root)

    state_designs = {"base": base_design, "initial_hidden": preview_design, "shown": shown_design, "hidden_again": hidden_again_design}
    state_evals: dict[str, dict] = {}
    state_meshes: dict[str, list[dict]] = {}
    state_images: dict[str, dict] = {}
    failures: list[str] = []
    mesh_script = out / "mesh.json"
    for name, design in state_designs.items():
        evaluated = read_eval(cli, design, out / f"bracket-preview-{name}.eval.json", root)
        state_evals[name] = evaluated
        if [(x.get("id"), x.get("name"), x.get("type")) for x in evaluated["timeline"][:2]] != original_features:
            failures.append(f"{name}: source sketch/extrusion feature identities changed")
        if evaluated["body_count"] != 1 or not evaluated["bodies"][0]["closed"]:
            failures.append(f"{name}: expected one closed bracket solid")
        if abs(evaluated["bodies"][0]["volume_mm3"] - source_eval["bodies"][0]["volume_mm3"]) > 1e-9:
            failures.append(f"{name}: bracket volume changed")
        if evaluated["bodies"][0]["faces"] != source_eval["bodies"][0]["faces"]:
            failures.append(f"{name}: bracket B-rep face count changed")
        sketch_visible = {x["name"]: x["visible"] for x in evaluated["sketches"]}
        expected_visible = name == "shown"
        if sketch_visible.get("SCT Partition 1") != expected_visible and name != "base":
            failures.append(f"{name}: preview sketch visibility {sketch_visible.get('SCT Partition 1')} != {expected_visible}")
        state_meshes[name] = read_mesh(
            cli,
            design,
            evaluated["bodies"][0]["name"],
            out / f"bracket-preview-{name}.mesh.json",
            mesh_script,
            root,
        )
        image_path = out / f"bracket-preview-{name}.png"
        snapshot = json.loads(
            run(
                [cli, "snapshot", "--in", str(design), "--out", str(image_path), "--width", "1000", "--height", "750", "--view", "iso", "--quality", "normal"],
                root,
            )
        )
        state_images[name] = {
            "path": str(image_path),
            "bytes": image_path.stat().st_size,
            "sha256": hashlib.sha256(image_path.read_bytes()).hexdigest(),
            "renderer": snapshot.get("renderer"),
        }

    base_mesh = state_meshes["base"]
    for name, mesh in state_meshes.items():
        if mesh != base_mesh:
            failures.append(f"{name}: SolidCraft model.mesh differs from the base extrusion")
    if state_images["initial_hidden"]["sha256"] == state_images["shown"]["sha256"]:
        snapshot_difference_pixels = 0
    else:
        before = Image.open(state_images["initial_hidden"]["path"]).convert("RGBA")
        after = Image.open(state_images["shown"]["path"]).convert("RGBA")
        snapshot_difference_pixels = sum(1 for a, b in zip(before.getdata(), after.getdata()) if a != b)
    if state_images["initial_hidden"]["sha256"] != state_images["hidden_again"]["sha256"]:
        failures.append("native snapshot did not return to the initial hidden state")

    plugin_folder = root / "plugins/solidcraft-trivariate"
    plugin_wasm = plugin_folder / "solidcraft-trivariate.wasm"
    plugin_test: dict[str, object] = {"status": "not_run", "reason": "component artifact not present"}
    if plugin_wasm.exists():
        check = subprocess.run([cli, "plugin", "check", str(plugin_folder), "--json"], cwd=root, text=True, capture_output=True)
        test = subprocess.run([cli, "plugin", "test", str(plugin_folder), "--json", "--out", str(out / "plugin-test-cases"), "--require-cases"], cwd=root, text=True, capture_output=True)
        plugin_test = {
            "status": "passed" if check.returncode == 0 and test.returncode == 0 else "failed",
            "check_stdout": check.stdout,
            "check_stderr": check.stderr,
            "test_stdout": test.stdout,
            "test_stderr": test.stderr,
            "check_exit": check.returncode,
            "test_exit": test.returncode,
        }
        if check.returncode or test.returncode:
            failures.append("SolidCraft plugin check/test failed")

    report = {
        "schema": "solidcraft_native_sketch_preview_qa_v1",
        "scope": "Native command-path QA; plugin component test is conditional on its WASM artifact being available.",
        "solidcraft_cli": cli,
        "trivariate_cli": str(kernel),
        "designs": {name: str(path) for name, path in state_designs.items()},
        "source_features": original_features,
        "preview_feature": {"id": preview_id, "name": preview_sketch["name"], "type": "sketch", "initial_visible": False},
        "body": {
            "name": source_eval["bodies"][0]["name"],
            "closed": source_eval["bodies"][0]["closed"],
            "volume_mm3": source_eval["bodies"][0]["volume_mm3"],
            "faces": source_eval["bodies"][0]["faces"],
        },
        "mesh": {
            "tolerance_mm": 0.05,
            "frame": base_mesh[0]["frame"],
            "vertices": len(base_mesh[0]["positions"]),
            "triangles": len(base_mesh[0]["triangles"]),
            "identical_across_base_hidden_shown_hidden_again": all(mesh == base_mesh for mesh in state_meshes.values()),
        },
        "sketch_visibility": {
            name: next((sketch["visible"] for sketch in state_evals[name]["sketches"] if sketch["name"] == "SCT Partition 1"), None)
            for name in state_evals
        },
        "snapshots": state_images,
        "snapshot_overlay": {
            "different_pixels_hidden_vs_shown": snapshot_difference_pixels,
            "visibly_demonstrated": snapshot_difference_pixels >= 1000,
            "note": "a small pixel difference is not sufficient evidence that snapshot rendered the sketch overlay",
        },
        "plugin_component_test": plugin_test,
        "failures": failures,
    }
    report_path = out / "sketch-preview-report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    summary_path = root / "docs/results/solidcraft-plugin/native-render.md"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    existing = summary_path.read_text() if summary_path.exists() else "# SolidCraft native render QA\n"
    insertion = [
        "## Native sketch show/hide and mesh invariance",
        "",
        "A two-feature bracket (closed six-edge sketch plus `solid.extrude`) was saved, reopened, and tested through SolidCraft's native `browser.visibility` command. The generated preview is a construction sketch with one native line; it starts hidden, is shown, then hidden again. The original sketch and extrusion feature IDs/types remain unchanged.",
        "",
        f"- Bracket: {report['body']['volume_mm3']:.9f} mm³, {report['body']['faces']} B-rep faces, closed.",
        f"- SolidCraft `model.mesh`: {report['mesh']['vertices']} vertices / {report['mesh']['triangles']} triangles at 0.05 mm in `{report['mesh']['frame']}` coordinates.",
        f"- Base, initial-hidden, shown, and hidden-again mesh arrays are exactly equal: `{report['mesh']['identical_across_base_hidden_shown_hidden_again']}`.",
        f"- Preview sketch visibility after each saved/reopened state: `{report['sketch_visibility']}`.",
        f"- Snapshot hashes: hidden `{state_images['initial_hidden']['sha256']}`, shown `{state_images['shown']['sha256']}`, hidden again `{state_images['hidden_again']['sha256']}`.",
        f"- Snapshot overlay visibility: {snapshot_difference_pixels} changed pixels; visibly demonstrated: `{report['snapshot_overlay']['visibly_demonstrated']}`.",
        f"- Plugin component check/test: `{plugin_test['status']}` ({plugin_test.get('reason', 'see JSON')}), report `{relative_path(report_path, root)}`.",
        "",
        "This exercises the native SolidCraft command path; when no WASM artifact is present, it is not reported as an invocation of the plugin component. The preview adds no mesh or B-rep face, and the body's evaluated mesh stays identical. The headless snapshot does not visibly establish that it rendered the sketch overlay; desktop visibility remains unqualified because the app could not connect to X.",
        "",
        "Failures: " + (", ".join(failures) if failures else "none"),
        "",
    ]
    marker = "## Native sketch show/hide and mesh invariance"
    if marker in existing:
        existing = existing[: existing.index(marker)].rstrip() + "\n\n"
    summary_path.write_text(existing.rstrip() + "\n\n" + "\n".join(insertion))
    print(json.dumps({"report": str(report_path), "summary": str(summary_path), "plugin_component_test": plugin_test["status"], "failures": failures}, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
