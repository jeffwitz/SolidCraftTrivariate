#!/usr/bin/env python3
"""Replay native primitive planning, audited build and both publication formats."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def run(cli, work):
    work.mkdir(parents=True, exist_ok=True)
    def command(*args):
        p = subprocess.run([cli, *args], capture_output=True, text=True, check=True)
        return json.loads(p.stdout)
    def point(x, y):
        return [[x, 1], [y, 1]]
    shapes = [
        {"kind": "rectangle", "min": point(3, -1), "max": point(7, 1)},
        {"kind": "circle", "center": point(5, 0), "radius": [1, 1]},
        {"kind": "ellipse", "center": point(5, 0), "radii": point(2, 1)},
    ]
    results = []
    for shape in shapes:
        for operation in [{"kind": "extrude", "length_mm": 3}, {"kind": "revolve"}]:
            name = shape["kind"] + "-" + operation["kind"]
            recipe = {"family": "primitive_sweep_v1", "primitive": shape, "operation": operation}
            path = work / (name + ".recipe.json")
            path.write_text(json.dumps(recipe))
            model, step, space = [work / (name + ext) for ext in [".model.json", ".step", ".space.json"]]
            plan = command("plan-sketch", path)
            built = command("build", path, model)
            audit = command("audit", model)
            assert built == audit and audit["status"] == "valid"
            count = (1 if shape["kind"] == "rectangle" else 5) * (4 if operation["kind"] == "revolve" else 1)
            assert count == plan["block_count"] == audit["blocks"]
            assert audit["determinant_lower_bound_mm3"] > 0
            command("export-step", model, step)
            assert step.read_text().strip().endswith("END-ISO-10303-21;")
            command("export-analysis", model, "4", space)
            package = json.loads(space.read_text())
            assert len(package["blocks"]) == count
            results.append({"recipe": recipe, "blocks": count, "audit": audit,
                            "step_bytes": step.stat().st_size, "space_schema": package["schema"],
                            "scalar_dofs": package["scalar_ndof"]})
    return {"schema": "primitive_sweep_cli_checks_v1", "cli_sha256": hashlib.sha256(cli.read_bytes()).hexdigest(), "cases": results}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cli", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(args.cli.resolve(), args.work), indent=2) + "\n")
