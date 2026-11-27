#!/usr/bin/env python3
"""Run the audited pair-parent CLI and an independent Fraction oracle corpus."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
FIELDCRAFT = [
    "fieldcraft_x40_y80",
    "fieldcraft_x50_y74",
    "fieldcraft_x51_y74",
    "fieldcraft_x64_y64",
    "fieldcraft_x80_y40",
]


def run(command):
    completed = subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=True)
    return completed.stdout


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def check_case(cli, name, model, work, expected_interfaces, expected_lateral):
    model = Path(model)
    if not model.is_absolute():
        model = ROOT / model
    work = Path(work)
    if not work.is_absolute():
        work = ROOT / work
    work.mkdir(parents=True, exist_ok=True)
    report = work / f"{name}.pairs.json"
    oracle = work / f"{name}.oracle.json"
    run([str(cli), "recover-pair-parents", str(model), str(report)])
    run([
        sys.executable,
        str(ROOT / "scripts/check_pair_parent_report.py"),
        str(model),
        str(report),
        "--output",
        str(oracle),
    ])
    native = json.loads(report.read_text())
    independent = json.loads(oracle.read_text())
    if native["interface_count"] != expected_interfaces:
        raise RuntimeError(f"{name}: source interface coverage changed")
    if independent["status"] != "pass" or independent["interfaces"] != expected_interfaces:
        raise RuntimeError(f"{name}: independent replay failed")
    axis_rows = {
        row["axis"]: row for row in native["counts_by_first_interface_axis"]
    }
    lateral_recovered = sum(axis_rows[axis]["recovered_count"] for axis in "uv")
    if lateral_recovered != expected_lateral:
        raise RuntimeError(f"{name}: expected {expected_lateral} u/v parent recoveries, got {lateral_recovered}")
    if independent["full_jet_candidates_refused_by_exact_parent_reconstruction"]:
        raise RuntimeError(f"{name}: a full normal-jet candidate failed exact parent reconstruction")
    return {
        "name": name,
        "model": model.name,
        "model_sha256": digest(model),
        "report": report.name,
        "report_sha256": digest(report),
        "oracle": oracle.name,
        "oracle_sha256": digest(oracle),
        "interfaces": native["interface_count"],
        "recovered": native["recovered_count"],
        "refused": native["refused_count"],
        "indeterminate": native["indeterminate_count"],
        "by_first_axis": axis_rows,
        "full_normal_jet_candidates_by_axis": independent["independent_full_normal_jet_candidates_by_axis"],
        "exact_parent_replay": independent["status"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate")
    parser.add_argument("--work", type=Path, default=ROOT / "build/pair-parent-corpus")
    parser.add_argument("--output", type=Path, default=ROOT / "build/pair-parent-corpus/summary.json")
    args = parser.parse_args()
    cli = args.cli.resolve()
    work = args.work.resolve()
    # The first native build writes into this directory before any report exists.
    work.mkdir(parents=True, exist_ok=True)
    if not cli.is_file():
        raise SystemExit(f"missing CLI executable: {cli}")

    multiloop_model = work / "nurbs-multiloop.model.json"
    run([str(cli), "build", str(ROOT / "tests/references/nurbs_multiloop.recipe.json"), str(multiloop_model)])
    cases = [
        ("s3-face-pocket", "docs/results/s3-face-pocket.model.json", 191, 5),
        ("nurbs-multiloop", multiloop_model, 75, 8),
        ("g2-cross-plate", "docs/results/g2-cross-plate.model.json", 46, 0),
    ]
    results = [check_case(cli, name, model, work, nfaces, lateral) for name, model, nfaces, lateral in cases]
    for name in FIELDCRAFT:
        result = check_case(cli, name, f"docs/site/data/{name}/{name}.model.json", work, 11, 0)
        if result["recovered"] != 0 or result["full_normal_jet_candidates_by_axis"] != {"u": 0, "v": 0, "w": 0}:
            raise RuntimeError(f"{name}: expected the retained FieldCraft interfaces to refuse same-parent recovery")
        results.append(result)
    summary = {
        "schema": "pair_bezier_parent_corpus_v1",
        "method": "native exact Fraction-ratio recovery with independent full-net de Casteljau and normal-jet replay",
        "claim_scope": "pairwise mathematical parent support only; no feature-history provenance or global geometric certification",
        "cli_sha256": digest(cli),
        "oracle_script": "scripts/check_pair_parent_report.py",
        "source_models": results,
        "totals": {
            "models": len(results),
            "interfaces": sum(row["interfaces"] for row in results),
            "recovered": sum(row["recovered"] for row in results),
            "refused": sum(row["refused"] for row in results),
            "indeterminate": sum(row["indeterminate"] for row in results),
            "status": "pass",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary["totals"], sort_keys=True))


if __name__ == "__main__":
    main()
