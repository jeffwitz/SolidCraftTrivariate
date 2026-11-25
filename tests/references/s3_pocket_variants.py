#!/usr/bin/env python3
"""Replay changes of pocket width, height and depth against independent section integration."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess

from g2_cross_model_volume import integrate
from g2_face_pocket_iga import pocket_volume_reference

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cli", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--work", type=Path, default=ROOT / "build/s3-pocket-variants")
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    cli = args.cli.resolve()
    base = json.loads((ROOT / "tests/references/s3_face_pocket.recipe.json").read_text())
    parent_volume = sum(integrate(base["parent_model"], 32).values())
    source = args.work / "parent.recipe.json"
    source.write_text(json.dumps(base["parent_model"]["source_recipe"]))
    results = []
    for name, xmin, xmax, ymin, depth in [
        ("wide", 54, 75, [5893, 266], 2),
        ("narrow", 56, 75, [5893, 266], 2),
        ("shallow", 55, 75, [5893, 266], 1),
        ("new_rational_cut", 54, 74, [48689, 2070], 2),
    ]:
        recipe = deepcopy(base)
        op = recipe["operation"]
        rect = op["loop"]["primitive"]
        rect["min"], rect["max"][0] = [[xmin, 1], ymin], [xmax, 1]
        op["depth_mm"] = depth
        path = args.work / (name + ".recipe.json")
        target = args.work / (name + ".model.json")
        path.write_text(json.dumps(recipe))
        proc = subprocess.run([cli, "build", path, target], capture_output=True, text=True)
        row = {"name": name, "operation": op, "returncode": proc.returncode, "passed": False}
        if proc.returncode:
            row["diagnostic"] = (proc.stdout + proc.stderr).strip()
        else:
            row["audit"] = json.loads(subprocess.check_output([cli, "audit", target]))
            volume = sum(integrate(json.loads(target.read_text()), 32).values())
            expected = pocket_volume_reference(source, xmin, xmax, float(Fraction(*ymin)), 25, depth)
            expected_removed = expected["quadrature"][-1]["removed_volume_mm3"]
            removed = parent_volume - volume
            error = abs(removed - expected_removed) / expected_removed
            row.update(removed_volume_mm3=removed, source_volume_mm3=expected_removed,
                       relative_error=error, passed=error <= 1.0e-8)
        results.append(row)
    args.output.write_text(json.dumps({
        "cli_sha256": hashlib.sha256(cli.read_bytes()).hexdigest(),
        "parent_q32_volume_mm3": parent_volume, "results": results,
    }, indent=2) + "\n")
    return int(any(not row["passed"] for row in results))


if __name__ == "__main__":
    raise SystemExit(main())
