#!/usr/bin/env python3
"""Independent oracle for the prospective row-bore CLI corpus.

The frozen manifest stores millimeters. This adapter converts its inputs to
the recipe's half-millimeter integer ticks and rejects values not exactly
representable in that unit. Geometry measurements use the generic rational
Bernstein evaluator in ``j1_geometry.py`` and analytic plate-minus-circles
volumes; no constructor geometry is generated here.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import j1_geometry


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "tests/references/row_bore_patterns.json"
DEFAULT_REPORT = ROOT / "build/row_bores_oracle.json"
SECTORS = ("east", "north", "west", "south")
SAMPLE = (0.0, 0.21132486540518713, 0.5, 0.7886751345948129, 1.0)


def half_ticks(value: float, label: str) -> int:
    ticks = 2.0 * value
    integer = round(ticks)
    if not math.isfinite(ticks) or abs(ticks - integer) > 1e-12:
        raise ValueError(f"{label}={value!r} is not an exact half-mm value")
    return int(integer)


def recipe(case: dict) -> dict:
    return {
        "family": "row_bores_v1",
        "half_widths_half_mm": [
            half_ticks(case["half_width_x"], f"{case['id']}.half_width_x"),
            half_ticks(case["half_width_y"], f"{case['id']}.half_width_y"),
        ],
        "centers_half_mm": [
            [half_ticks(x, f"{case['id']}.center[{i}].x"),
             half_ticks(y, f"{case['id']}.center[{i}].y")]
            for i, (x, y) in enumerate(case["centers"])
        ],
        "radii_half_mm": [
            half_ticks(radius, f"{case['id']}.radius[{i}]")
            for i, radius in enumerate(case["radii"])
        ],
    }


def analytic_volume(case: dict) -> float:
    a, b = case["half_width_x"], case["half_width_y"]
    removed = math.pi * sum(radius * radius for radius in case["radii"])
    return case["height"] * (4.0 * a * b - removed)


def radial_errors(blocks: list[dict], case: dict) -> dict:
    errors = []
    centers = case["centers"]
    radii = case["radii"]
    for cell, (center, radius) in enumerate(zip(centers, radii)):
        cx, cy = center
        for sector in SECTORS:
            block_id = f"row_bore_cell_{cell:02}_sector_{sector}_z0_8"
            block = next((b for b in blocks if b.get("id") == block_id), None)
            if block is None:
                errors.append({"block": block_id, "missing": True})
                continue
            for u in SAMPLE:
                point, _ = j1_geometry.evaluate(block, (u, 1.0, 0.5))
                error = abs(math.hypot(point[0] - cx, point[1] - cy) - radius)
                errors.append({"block": block_id, "u": u, "error_mm": error})
    return {
        "max_error_mm": max((item.get("error_mm", math.inf) for item in errors), default=math.inf),
        "samples": len(errors),
        "missing_blocks": sum("missing" in item for item in errors),
    }


def orientation(blocks: list[dict]) -> dict:
    determinants = []
    bad = []
    for block in blocks:
        for u in SAMPLE[1:-1]:
            for v in SAMPLE[1:-1]:
                for w in SAMPLE[1:-1]:
                    point, jac = j1_geometry.evaluate(block, (u, v, w))
                    det = j1_geometry.determinant(jac)
                    determinants.append(det)
                    if not math.isfinite(det) or det <= 0:
                        bad.append({"block": block.get("id"), "uvw": [u, v, w], "det": det})
    return {
        "sample_count": len(determinants),
        "minimum_sampled_det_mm3": min(determinants, default=None),
        "nonpositive_samples": bad,
    }


def analyze(case: dict, model: dict) -> dict:
    blocks = model.get("blocks", [])
    n = len(case["centers"])
    expected_ids = [f"row_bore_cell_{cell:02}_sector_{sector}_z0_8"
                    for cell in range(n) for sector in SECTORS]
    actual_ids = [block.get("id") for block in blocks]
    orders = {str(order): j1_geometry.volume(blocks, order) for order in (4, 8, 12, 16)}
    topology = j1_geometry.topology(blocks)
    volume_ref = analytic_volume(case)
    circles = radial_errors(blocks, case)
    orient = orientation(blocks)
    expected_interfaces = 5 * n - 1
    expected_exterior = 14 * n + 2
    expected_euler = 2 - 2 * n
    interfaces_match = topology["interfaces"] == expected_interfaces
    exterior_match = topology["exterior_faces"] == expected_exterior
    euler_match = topology["boundary_euler_sampled"] == expected_euler
    volume_error = abs(orders["16"] - volume_ref)
    qualified = (
        model.get("schema") == "row_bores_v1"
        and len(blocks) == 4 * n
        and actual_ids == expected_ids
        and interfaces_match and exterior_match and euler_match
        and topology["nonmanifold_faces"] == 0
        and circles["missing_blocks"] == 0
        and circles["max_error_mm"] <= 1e-12
        and not orient["nonpositive_samples"]
        and volume_error <= 2e-10 * max(volume_ref, 1.0)
    )
    return {
        "status": "measured",
        "qualified_by_oracle": qualified,
        "schema": model.get("schema"),
        "block_count": len(blocks),
        "expected_block_count": 4 * n,
        "ids_in_left_to_right_sector_order": actual_ids == expected_ids,
        "volume_by_gauss_order_mm3": orders,
        "volume_analytic_mm3": volume_ref,
        "volume_error_order16_mm3": volume_error,
        "quadrature_delta_12_16_mm3": abs(orders["16"] - orders["12"]),
        "interfaces_expected": expected_interfaces,
        "exterior_faces_expected": expected_exterior,
        "boundary_euler_expected": expected_euler,
        **topology,
        "topology_matches": interfaces_match and exterior_match and euler_match,
        "radial_error": circles,
        "orientation": orient,
    }


def run_cli(cli: str, recipe_data: dict, output: Path) -> dict:
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="row-bores-") as temp:
        recipe_path = Path(temp) / "recipe.json"
        recipe_path.write_text(json.dumps(recipe_data, indent=2) + "\n")
        start = time.perf_counter_ns()
        built = subprocess.run([cli, "build", str(recipe_path), str(output)],
                               capture_output=True, text=True, check=False)
        build_ms = (time.perf_counter_ns() - start) / 1e6
        try:
            build_status = json.loads(built.stdout)
        except json.JSONDecodeError:
            build_status = {"unparsed_stdout": built.stdout.strip()}
        result = {"build_returncode": built.returncode,
                  "build_wall_ms": build_ms,
                  "build_status": build_status,
                  "stderr": built.stderr.strip()}
        if built.returncode != 0 or not output.exists():
            return result
        start = time.perf_counter_ns()
        audited = subprocess.run([cli, "audit", str(output)], capture_output=True,
                                 text=True, check=False)
        result["audit_wall_ms"] = (time.perf_counter_ns() - start) / 1e6
        result["audit_returncode"] = audited.returncode
        try:
            result["audit_status"] = json.loads(audited.stdout)
        except json.JSONDecodeError:
            result["audit_status"] = {"unparsed_stdout": audited.stdout.strip()}
        result["audit_stderr"] = audited.stderr.strip()
        return result


def run_local_radius_change(cli: str, case: dict, baseline_model: dict, model_dir: Path) -> dict:
    if case["id"] != "row_4":
        return {"status": "not_applicable"}
    changed = dict(case)
    changed["radii"] = list(case["radii"])
    changed["radii"][0] += 0.5
    target = model_dir / "radius_perturbation_row_4.json"
    cli_result = run_cli(cli, recipe(changed), target)
    if cli_result["build_returncode"] != 0 or not target.exists():
        return {"status": "refused", **cli_result}
    perturbed = json.loads(target.read_text())
    before = {block["id"]: block["controls"] for block in baseline_model.get("blocks", [])}
    after = {block["id"]: block["controls"] for block in perturbed.get("blocks", [])}
    changed_ids = sorted(key for key in before.keys() & after.keys() if before[key] != after[key])
    expected_ids = sorted(f"row_bore_cell_00_sector_{sector}_z0_8" for sector in SECTORS)
    return {
        "status": "measured",
        "changed_block_ids": changed_ids,
        "expected_changed_block_ids": expected_ids,
        "only_four_networks_changed": changed_ids == expected_ids,
        "build_wall_ms": cli_result["build_wall_ms"],
        "audit_wall_ms": cli_result.get("audit_wall_ms"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", default="target/release/solidcraft-trivariate")
    parser.add_argument("--models", help="directory containing <case-id>.json model exports")
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--require-qualified", action="store_true")
    args = parser.parse_args()

    manifest = json.loads(MANIFEST.read_text())
    cases = manifest["regular_rows"] + manifest["interacting_pair"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    model_dir = args.output.parent / "row_bores_models"
    model_dir.mkdir(parents=True, exist_ok=True)
    report = {"manifest_version": manifest["version"], "frozen_base": manifest["frozen_base"],
              "unit_conversion": "manifest mm to exact half-mm integer recipe ticks",
              "families": {"regular_rows": [case["id"] for case in manifest["regular_rows"]],
                           "interacting_pair": [case["id"] for case in manifest["interacting_pair"]]},
              "cases": {}}
    models = {}
    for case in cases:
        cid = case["id"]
        target = Path(args.models) / f"{cid}.json" if args.models else model_dir / f"{cid}.json"
        timings = None
        if args.models:
            if not target.exists():
                report["cases"][cid] = {"status": "refused", "reason": f"missing model {target}"}
                continue
            model = json.loads(target.read_text())
        else:
            timings = run_cli(args.cli, recipe(case), target)
            if timings["build_returncode"] != 0 or not target.exists():
                report["cases"][cid] = {"status": "refused", "expected_state": case.get("geometry_state"),
                                         **timings}
                continue
            model = json.loads(target.read_text())
            models[cid] = model
        measured = analyze(case, model)
        if timings is not None:
            measured["cli"] = timings
        measured["geometry_state_manifest"] = case.get("geometry_state", "disjoint")
        report["cases"][cid] = measured

    regular_counts = [report["cases"].get(case["id"], {}).get("block_count")
                      for case in manifest["regular_rows"]]
    report["regular_series_block_counts"] = regular_counts
    report["regular_series_scales_with_n"] = regular_counts == [4 * len(case["centers"])
                                                                 for case in manifest["regular_rows"]]
    if not args.models and "row_4" in models:
        report["local_radius_perturbation"] = run_local_radius_change(
            args.cli, next(case for case in manifest["regular_rows"] if case["id"] == "row_4"),
            models["row_4"], model_dir)
    report["coverage"] = {
        "measured": sum(value.get("status") == "measured" for value in report["cases"].values()),
        "required": len(cases),
        "refused": sum(value.get("status") == "refused" for value in report["cases"].values()),
        "qualified_by_oracle": sum(value.get("qualified_by_oracle", False)
                                   for value in report["cases"].values()),
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.require_qualified:
        failures = report["coverage"]["qualified_by_oracle"] != len(cases)
        failures |= not report["regular_series_scales_with_n"]
        perturbation = report.get("local_radius_perturbation", {})
        failures |= perturbation.get("only_four_networks_changed") is False
        if failures:
            raise SystemExit(1)


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
