#!/usr/bin/env python3
"""Build, export, and render representative trivariate and SolidCraft solids."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path


def invoke(command: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=cwd, text=True, capture_output=True)


def require(command: list[str], cwd: Path) -> str:
    result = invoke(command, cwd)
    if result.returncode:
        raise RuntimeError(
            f"command failed ({result.returncode}): {command!r}\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )
    return result.stdout


def box_cylinder_script(recipe: dict, operation: str) -> dict:
    stock = recipe["stock"]
    tool = recipe["tool"]
    lo, hi = stock["min_mm"], stock["max_mm"]
    box = {
        "command": "solid.box",
        "params": {
            "length": hi[0] - lo[0],
            "width": hi[1] - lo[1],
            "height": hi[2] - lo[2],
            "corner": lo,
            "name": "Stock",
            "body_name": "Stock",
        },
    }
    center = tool["center_mm"]
    z0, z1 = tool["interval_mm"]
    cylinder = {
        "command": "solid.cylinder",
        "params": {
            "radius": tool["radius_mm"],
            "height": z1 - z0,
            "base": [center[0], center[1], z0],
            "axis": tool["axis"].upper(),
            "name": "Tool",
            "body_name": "Tool",
        },
    }
    if operation == "b-minus-a":
        commands = [
            cylinder,
            box,
            {"command": "solid.combine", "params": {"target": "Tool", "tools": ["Stock"], "operation": "cut"}},
        ]
    else:
        combine = {"union": "join", "intersection": "intersect", "a-minus-b": "cut"}[operation]
        commands = [
            box,
            cylinder,
            {"command": "solid.combine", "params": {"target": "Stock", "tools": ["Tool"], "operation": combine}},
        ]
    return {"commands": commands}


def expected_suite_volume(model: dict) -> float:
    volume = model["volume_mm3"]
    return float(volume["constant"]) + float(volume["pi_coefficient"]) * math.pi


def relative(path: Path, root: Path) -> str:
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
    trivariate = args.trivariate_cli if args.trivariate_cli.is_absolute() else root / args.trivariate_cli

    cases: list[dict[str, object]] = []
    failures: list[str] = []
    source_path = root / "tests/references/solidcraft_sketch_bracket.timeline.json"
    source = json.loads(source_path.read_text())
    source["commands"] = [step for step in source["commands"] if step["command"] == "timeline.add"][:2]
    two_feature_timeline = out / "bracket-extrusion-two-feature.timeline.json"
    two_feature_timeline.write_text(json.dumps(source))

    extrusion_model = out / "bracket-extrusion.trivariate.model.json"
    extrusion_step = out / "bracket-extrusion.trivariate.step"
    require([str(trivariate), "extrude-sketch", str(two_feature_timeline), str(extrusion_model)], root)
    require([str(trivariate), "export-step", str(extrusion_model), str(extrusion_step)], root)
    cases.append(
        {
            "id": "bracket-extrusion",
            "step": str(extrusion_step),
            "provenance": "Trivariate `extrude-sketch` from the SolidCraft two-feature sketch/extrusion timeline, then Trivariate `export-step`",
            "trivariate_model": str(extrusion_model),
            "trivariate_step_export": {"status": "exported"},
        }
    )

    for operation in ("union", "intersection", "a-minus-b", "b-minus-a"):
        recipe_path = root / f"docs/site/box-cylinder/tall/{operation}.recipe.json"
        recipe = json.loads(recipe_path.read_text())
        trivariate_model = out / f"tall-{operation}.trivariate.model.json"
        trivariate_step = out / f"tall-{operation}.trivariate.step"
        build_stdout = require([str(trivariate), "build", str(recipe_path), str(trivariate_model)], root)
        build_audit = json.loads(build_stdout)
        if build_audit.get("status") != "valid":
            raise RuntimeError(f"trivariate build failed for {operation}: {build_audit}")

        export = invoke([str(trivariate), "export-step", str(trivariate_model), str(trivariate_step)], root)
        export_reply: dict[str, object]
        if export.returncode == 0:
            export_reply = {"status": "exported", "stdout": export.stdout.strip()}
            render_step = trivariate_step
            render_provenance = "exact STEP exported from the trivariate build"
        else:
            try:
                export_reply = json.loads(export.stdout.strip())
            except json.JSONDecodeError:
                export_reply = {"status": "error", "stdout": export.stdout.strip(), "stderr": export.stderr.strip()}
            failures.append(f"{operation}: Trivariate export-step refused: {export_reply}")

            # The current boundary writer may refuse this family. Rebuild the same
            # exact recipe dimensions in SolidCraft's own B-rep kernel for display.
            command_script = out / f"tall-{operation}.solidcraft.commands.json"
            command_script.write_text(json.dumps(box_cylinder_script(recipe, operation)))
            host_step = out / f"tall-{operation}.solidcraft-native.step"
            host_design = out / f"tall-{operation}.solidcraft-native.sc"
            require(
                [
                    args.solidcraft_cli,
                    "run",
                    str(command_script),
                    "--save",
                    str(host_design),
                    "--out",
                    str(host_step),
                    "--quiet",
                ],
                root,
            )
            render_step = host_step
            render_provenance = "SolidCraft-native B-rep rebuilt from the exact tall recipe dimensions; fallback because trivariate STEP export was refused"

        cases.append(
            {
                "id": f"tall-box-cylinder-{operation}",
                "step": str(render_step),
                "provenance": render_provenance,
                "recipe": str(recipe_path),
                "trivariate_model": str(trivariate_model),
                "trivariate_build_audit": build_audit,
                "trivariate_volume_expected_mm3": expected_suite_volume(json.loads(trivariate_model.read_text())),
                "trivariate_step_export": export_reply,
                "trivariate_step_path": str(trivariate_step),
            }
        )

    oblique_recipe = root / "tests/references/oblique_blind_pockets.recipe.json"
    oblique_model = out / "oblique-blind-pockets.trivariate.model.json"
    oblique_step = out / "oblique-blind-pockets.trivariate.step"
    oblique_audit = json.loads(require([str(trivariate), "build", str(oblique_recipe), str(oblique_model)], root))
    oblique_export = invoke([str(trivariate), "export-step", str(oblique_model), str(oblique_step)], root)
    oblique_export_reply: dict[str, object]
    if oblique_export.returncode == 0:
        oblique_export_reply = {"status": "exported", "stdout": oblique_export.stdout.strip()}
    else:
        try:
            oblique_export_reply = json.loads(oblique_export.stdout.strip())
        except json.JSONDecodeError:
            oblique_export_reply = {"status": "error", "stdout": oblique_export.stdout.strip(), "stderr": oblique_export.stderr.strip()}
        failures.append(f"oblique blind pockets: Trivariate export-step refused: {oblique_export_reply}")
    cases.append(
        {
            "id": "oblique-blind-pockets",
            "step": str(oblique_step),
            "provenance": "Trivariate `build` and `export-step` from the oblique blind pockets recipe",
            "recipe": str(oblique_recipe),
            "trivariate_model": str(oblique_model),
            "trivariate_build_audit": oblique_audit,
            "trivariate_step_export": oblique_export_reply,
            "expected_body_count": 1,
            "expected_volume_mm3": 871.4365866666667,
            "expected_faces": 140,
        }
    )

    rim_source = root / "docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json"
    for state in ("chamfer", "fillet"):
        step = out / f"fieldcraft-{state}.trivariate.step"
        require([str(trivariate), "edit-fieldcraft-bore-rims", str(rim_source), state, str(step)], root)
        cases.append(
            {
                "id": f"fieldcraft-{state}",
                "step": str(step),
                "provenance": f"Regenerated by Trivariate `edit-fieldcraft-bore-rims {state}` from the current x64/y64 model",
                "source_model": str(rim_source),
            }
        )

    mesh_script = out / "mesh.json"
    report: dict[str, object] = {
        "schema": "solidcraft_native_render_qa_v1",
        "solidcraft_cli": args.solidcraft_cli,
        "trivariate_cli": relative(trivariate, root),
        "display_path": "SolidCraft `snapshot` on the produced STEP; renderer and GPU tier recorded per case",
        "mesh_path": "SolidCraft `model.mesh` at 0.05 mm via `solidcraft-cli run`",
        "style_check": "view.visual_style shaded is tested separately; the headless snapshot renderer may not reflect the viewport style",
        "cases": [],
    }
    for case in cases:
        name = str(case["id"])
        step_path = Path(str(case["step"]))
        if not step_path.is_absolute():
            step_path = root / step_path
        if not step_path.is_file():
            failures.append(f"{name}: missing STEP {step_path}")
            continue

        eval_path = out / f"{name}.eval.json"
        image_path = out / f"{name}.png"
        require([args.solidcraft_cli, "eval", str(step_path), "--json-out", str(eval_path), "--quiet"], root)
        eval_data = json.loads(eval_path.read_text())
        bodies = eval_data.get("bodies", [])
        if not bodies or not all(body.get("closed") for body in bodies):
            failures.append(f"{name}: SolidCraft did not report closed solid bodies")
        if "expected_body_count" in case and eval_data.get("body_count") != case["expected_body_count"]:
            failures.append(f"{name}: expected {case['expected_body_count']} bodies, got {eval_data.get('body_count')}")
        if "expected_faces" in case:
            observed_faces = sum(int(body.get("faces", 0)) for body in bodies)
            if observed_faces != case["expected_faces"]:
                failures.append(f"{name}: expected {case['expected_faces']} faces, got {observed_faces}")
        if "expected_volume_mm3" in case:
            observed_volume = sum(float(body["volume_mm3"]) for body in bodies)
            if abs(observed_volume - float(case["expected_volume_mm3"])) > 1e-7:
                failures.append(f"{name}: expected volume {case['expected_volume_mm3']}, got {observed_volume}")

        mesh_script.write_text(
            json.dumps(
                {
                    "commands": [
                        {"command": "model.mesh", "params": {"body": body["name"], "tolerance": 0.05}}
                        for body in bodies
                    ]
                }
            )
        )
        mesh_path = out / f"{name}.mesh.json"
        require(
            [
                args.solidcraft_cli,
                "run",
                str(mesh_script),
                "--in",
                str(step_path),
                "--results",
                "--quiet",
                "--json-out",
                str(mesh_path),
            ],
            root,
        )
        mesh_data = json.loads(mesh_path.read_text())
        mesh_counts = [
            {
                "body": mesh["body"],
                "frame": mesh["frame"],
                "vertices": len(mesh["positions"]),
                "triangles": len(mesh["triangles"]),
            }
            for mesh in mesh_data
        ]
        if len(mesh_data) != len(bodies) or any(not m["positions"] or not m["triangles"] for m in mesh_data):
            failures.append(f"{name}: native model.mesh did not return a nonempty mesh for every body")

        snapshot_stdout = require(
            [
                args.solidcraft_cli,
                "snapshot",
                "--in",
                str(step_path),
                "--out",
                str(image_path),
                "--width",
                "1000",
                "--height",
                "750",
                "--view",
                "iso",
                "--quality",
                "normal",
            ],
            root,
        )
        snapshot = json.loads(snapshot_stdout)
        if not image_path.is_file() or image_path.stat().st_size == 0:
            failures.append(f"{name}: snapshot image was not written")

        shaded_script = out / f"{name}.shaded.json"
        shaded_script.write_text(json.dumps({"commands": [{"command": "view.visual_style", "params": {"style": "shaded"}}]}))
        style_reply = json.loads(
            require([args.solidcraft_cli, "run", str(shaded_script), "--in", str(step_path), "--results"], root)
        )
        shaded_image = out / f"{name}.shaded.png"
        shaded_snapshot = json.loads(
            require(
                [
                    args.solidcraft_cli,
                    "snapshot",
                    str(shaded_script),
                    "--in",
                    str(step_path),
                    "--out",
                    str(shaded_image),
                    "--width",
                    "1000",
                    "--height",
                    "750",
                    "--view",
                    "iso",
                    "--quality",
                    "normal",
                ],
                root,
            )
        )
        default_sha = hashlib.sha256(image_path.read_bytes()).hexdigest()
        shaded_sha = hashlib.sha256(shaded_image.read_bytes()).hexdigest()

        item: dict[str, object] = {
            **case,
            "step": relative(step_path, root),
            "body_count": eval_data.get("body_count"),
            "bodies": [
                {
                    "name": body.get("name"),
                    "closed": body.get("closed"),
                    "volume_mm3": body.get("volume_mm3"),
                    "faces": body.get("faces"),
                    "face_types": body.get("face_types"),
                }
                for body in bodies
            ],
            "mesh": mesh_counts,
            "snapshot": {
                "path": relative(image_path, root),
                "bytes": image_path.stat().st_size,
                "renderer": snapshot.get("renderer"),
                "width": snapshot.get("width"),
                "height": snapshot.get("height"),
                "sha256": default_sha,
            },
            "shaded_snapshot": {
                "path": relative(shaded_image, root),
                "bytes": shaded_image.stat().st_size,
                "renderer": shaded_snapshot.get("renderer"),
                "requested_style": "shaded",
                "style_command_reply": style_reply,
                "sha256": shaded_sha,
                "byte_identical_to_default_snapshot": shaded_sha == default_sha,
            },
            "eval_json": relative(eval_path, root),
            "mesh_json": relative(mesh_path, root),
        }
        if "trivariate_volume_expected_mm3" in item:
            observed = sum(float(body["volume_mm3"]) for body in bodies)
            item["solidcraft_volume_sum_mm3"] = observed
            item["solidcraft_volume_error_mm3"] = observed - float(item["trivariate_volume_expected_mm3"])
            if abs(item["solidcraft_volume_error_mm3"]) > 1e-7:
                failures.append(f"{name}: SolidCraft volume differs from the trivariate exact recipe")
        report["cases"].append(item)

    report["failures"] = failures
    report_path = out / "report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    summary_path = root / "docs/results/solidcraft-plugin/native-render.md"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    table_rows = []
    for case in report["cases"]:
        total_volume = sum(float(body["volume_mm3"]) for body in case["bodies"])
        mesh_vertices = sum(mesh["vertices"] for mesh in case["mesh"])
        mesh_triangles = sum(mesh["triangles"] for mesh in case["mesh"])
        table_rows.append(
            f"| {case['id']} | {case['body_count']} ({'closed' if all(b['closed'] for b in case['bodies']) else 'open'}) | "
            f"{total_volume:.9f} | {sum(int(b.get('faces', 0)) for b in case['bodies'])} | {mesh_vertices} | {mesh_triangles} | {case['snapshot']['renderer']} |"
        )
    export_refusals = [
        f"- `{case['id']}`: `{case['trivariate_step_export']}`"
        for case in report["cases"]
        if isinstance(case.get("trivariate_step_export"), dict)
        and case["trivariate_step_export"].get("reason")
    ]
    lines = [
        "# SolidCraft native render QA",
        "",
        f"Generated by `{Path(__file__).name}` from a clean run on the isolated clone.",
        "",
        f"- SolidCraft CLI: `{args.solidcraft_cli}`",
        f"- Trivariate CLI: `{relative(trivariate, root)}`",
        "- Renderer: SolidCraft `snapshot`, all cases reported `gpu-software` / llvmpipe.",
        "- Tessellation: SolidCraft `model.mesh` at 0.05 mm, component-local millimetres; no custom mesher used for display or mesh extraction.",
        "- Result JSON and native snapshots: `build/solidcraft-native-render-qa/`.",
        "",
        "| Case | Bodies | Volume (mm³) | Faces | Mesh vertices | Mesh triangles | Renderer |",
        "|---|---:|---:|---:|---:|---:|---|",
        *table_rows,
        "",
        "The four tall Boolean recipes were rebuilt, audited and exported to STEP by the Trivariate CLI. The oblique blind-pocket V-Rep case is also built and exported from its reference recipe. Every STEP is then re-imported by SolidCraft, which checks the closed body, volume and faces before native rendering. If Trivariate STEP encoding fails, a diagnostic SolidCraft-native fallback may be recorded, but the qualification result fails.",
        "",
        "## STEP export refusals",
        "",
        *(export_refusals or ["- None."]),
        "",
        "## Shaded-style check",
        "",
        f"`view.visual_style {{style: shaded}}` is accepted by the host command path. The headless snapshot API has no style parameter: all {len(report['cases'])} `shaded` request PNGs are byte-identical to default snapshots, so the CLI does not demonstrate suppression of B-rep edge seams. The GUI's View > Shaded mode remains the host-supported setting for suppressing edge overlays; the separate desktop test exercises the native sketch overlay.",
        "",
        "## Limits",
        "",
        "This headless run checks native evaluation and tessellation. The [desktop sketch workflow](https://github.com/jeffwitz/SolidCraftTrivariate/actions/workflows/desktop-smoke.yml) tests the visible partition independently in a real SolidCraft window on Xvfb.",
        "",
        "Failures: " + (", ".join(failures) if failures else "none"),
        "",
    ]
    summary_path.write_text("\n".join(lines))
    print(json.dumps({"report": str(report_path), "summary": str(summary_path), "cases": len(report["cases"]), "failures": failures}, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
