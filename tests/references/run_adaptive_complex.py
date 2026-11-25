#!/usr/bin/env python3
"""Build six release adaptive packages and run the independent X3 oracle."""

import json
import os
from pathlib import Path
import subprocess
import sys
import time
import platform

import numpy
import scipy
from adaptive_homology import analyze

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "build/adaptive-complex"
OUTPUT = ROOT / "docs/results/adaptive-complex-20261125.json"


def run(args, *, env=None):
    return subprocess.run(list(map(str, args)), cwd=ROOT, env=env, check=True,
                          capture_output=True, text=True)


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    run(["rtk", "cargo", "build", "--release", "--locked", "--offline"])
    cli = ROOT / "target/release/solidcraft-trivariate"
    packages = []
    package_rows = []
    for radius in (64, 56):
        recipe = WORK / f"x{radius}.recipe.json"
        model = WORK / f"x{radius}.model.json"
        recipe.write_text(json.dumps({"family": "parameterized_bracket_v1",
                                      "radii_sixteenths_mm": [radius, 64]}))
        run([cli, "build", recipe, model])
        for depth in range(3):
            focus = "hole_x_sector_east_z10_30"
            steps = [
                {"refine": [{"root": focus, "level": level,
                             "index": [0, 0, 0] if level == 0 else [1, 0, 0]}]}
                for level in range(depth)
            ]
            plan = WORK / f"x{radius}-depth{depth}.plan.json"
            package = WORK / f"x{radius}-depth{depth}.package.json"
            plan.write_text(json.dumps({"schema": "adaptive_refinement_plan_v1",
                                        "max_depth": 4, "max_cells": 1500,
                                        "steps": steps}))
            start = time.monotonic()
            run([cli, "export-adaptive", model, plan, 8, package])
            export_seconds = time.monotonic() - start
            packages.append(package)
            package_rows.append({"radius_ticks": [radius, 64], "depth": depth,
                                 "export_wall_s": export_seconds,
                                 "package_bytes": package.stat().st_size})

    d4 = WORK / "two-cubes-d4.json"
    env = os.environ.copy()
    env["ADAPTIVE_COMPLEX_D4_OUTPUT"] = str(d4.relative_to(ROOT))
    d4_start = time.monotonic()
    run(["rtk", "cargo", "test", "--release", "--test", "adaptive_complex",
         "optionally_write_d4_complex_for_independent_oracle"], env=env)
    d4_generation_seconds = time.monotonic() - d4_start

    cases = []
    for package in packages:
        data = json.loads(package.read_text())
        roots = {block["id"]: block for block in data["source_model"]["blocks"]}
        complex_data = data["adaptive_mesh"]["oriented_complex"]
        result = analyze(complex_data, roots, 2)
        cases.append({"case": package.stem, "source": "adaptive_package", **result})
    d4_data = json.loads(d4.read_text())
    d4_result = analyze(d4_data, None, 0)
    cases.append({"case": "two-cubes-d4-refined-a", "source": "standalone_complex",
                  **d4_result})
    report = {
        "schema": "adaptive_oriented_complex_qualification_v1",
        "revision": run(["rtk", "git", "rev-parse", "HEAD"]).stdout.strip(),
        "commands": [
            "rtk cargo build --release --locked --offline",
            "rtk cargo test --release --test adaptive_complex optionally_write_d4_complex_for_independent_oracle",
            "rtk python3 tests/references/run_adaptive_complex.py",
        ],
        "oracle": "tests/references/adaptive_homology.py using exact GF(2) rank and j1_geometry chart evaluation",
        "environment": {"python": platform.python_version(), "numpy": numpy.__version__,
                        "scipy": scipy.__version__, "build": "Cargo release, locked offline"},
        "provenance_parameters": [0.25, 0.5, 0.75],
        "cases": cases,
        "package_costs": package_rows,
        "d4_complex_generation_wall_s": d4_generation_seconds,
        "limits": ["GF(2) homology is a diagnostic, not a general proof over all coefficient fields."],
    }
    OUTPUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"report": str(OUTPUT), "cases": len(cases),
                      "max_export_wall_s": max(row["export_wall_s"] for row in package_rows)}, indent=2))


if __name__ == "__main__":
    main()
