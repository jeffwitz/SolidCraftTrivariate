#!/usr/bin/env python3
"""Independent oracle and budget campaign for directional face quadrature."""

from __future__ import annotations

import argparse
import copy
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import adaptive_surface_load as surface_oracle

CORPUS = ROOT / "tests/references/circular_polygon_corpus.json"
WORK = ROOT / "build/directional-surface-load"
REPORT = ROOT / "docs/results/directional-surface-load-20261125.json"
SENTINEL = b"previous accepted result: directional surface quadrature\n"
TOLERANCE = 1e-9
DEPTH = 8
BUDGETS = (4096, 65536)


def build_cases():
    frozen = {item["id"]: item["recipe"] for item in json.loads(CORPUS.read_text())}
    x_shear = copy.deepcopy(frozen["extreme_shear"])
    x_shear["top_shift_mm"] = [32, 0]
    y_shear = copy.deepcopy(frozen["extreme_shear"])
    y_shear["top_shift_mm"] = [0, 32]
    return [
        {"id": "right_circle_pressure", "recipe_id": "one_round", "recipe": frozen["one_round"],
         "feature": "bore/pilot", "law": {"pressure_mpa": 0.05}, "map_direction": "right_circular"},
        {"id": "oblique_circle_stress", "recipe_id": "oblique_profile", "recipe": frozen["oblique_profile"],
         "feature": "bore/left", "law": {"stress_mpa": [[0, .3, .1], [.3, 0, -.2], [.1, -.2, 0]]},
         "map_direction": "oblique_circular"},
        {"id": "shear_x_vector", "recipe_id": "shear_x", "recipe": x_shear,
         "feature": "bore/micro", "law": {"vector_mpa": [.2, -.4, .7]}, "map_direction": "x_shear"},
        {"id": "shear_y_pressure", "recipe_id": "shear_y", "recipe": y_shear,
         "feature": "bore/micro", "law": {"pressure_mpa": .05}, "map_direction": "y_shear"},
        {"id": "shear_xy_mixed_stress", "recipe_id": "extreme_shear",
         "recipe": frozen["extreme_shear"], "feature": "bore/micro",
         "law": {"stress_mpa": [[0, .3, .1], [.3, 0, -.2], [.1, -.2, 0]]},
         "map_direction": "mixed_xy_shear"},
        {"id": "shear_xy_end_plane_pressure", "recipe_id": "extreme_shear",
         "recipe": frozen["extreme_shear"],
         "feature": {"feature": "end", "plane": {"axis": 2, "coordinate_mm": 1}},
         "law": {"pressure_mpa": .05}, "map_direction": "angular_radial_chart"},
        {"id": "shear_xy_rectangle_vector", "recipe_id": "extreme_shear",
         "recipe": frozen["extreme_shear"], "feature": "edge/outer/v00/v01",
         "law": {"vector_mpa": [.2, .4, -.1]}, "map_direction": "mixed_xy_shear_rectangle"},
    ]


def write_json(path, value):
    surface_oracle.write_json(path, value)


def independent_comparison(result, law):
    q96 = surface_oracle.reference_result_load(result, law, 96)
    q128 = surface_oracle.reference_result_load(result, law, 128)
    load = result["loads"][0]
    force_scale = max(q128["absolute_force_scale_n"], 1.0)
    moment_scale = max(q128["absolute_moment_scale_n_mm"], 1.0)
    nodal = np.asarray(result["force_n_xyz"], dtype=float)
    expected = q128["force_n_xyz"]
    displacement = np.asarray(result["displacement_mm_xyz"], dtype=float)
    reference_work = float(expected @ displacement)
    return {
        "reference_order": 128,
        "q96_q128_force_delta": float(np.linalg.norm(q96["force_n_xyz"] - expected) / force_scale),
        "area_reference_mm2": q128["area_mm2"],
        "area_relative_error": abs(float(load["area_mm2"]) / q128["area_mm2"] - 1.0),
        "full_nodal_force_relative_error": float(np.linalg.norm(nodal - expected) / force_scale),
        "resultant_relative_error": float(np.linalg.norm(np.asarray(load["resultant_n"]) - q128["resultant_n"]) / force_scale),
        "moment_relative_error": float(np.linalg.norm(np.asarray(load["moment_n_mm"]) - q128["moment_n_mm"]) / moment_scale),
        "reference_work_n_mm": reference_work,
        "work_relative_error": abs(reference_work-float(result["traction_work_n_mm"])) / max(abs(reference_work), 1.0),
        "faces": q128["face_count"],
        "panels": load.get("integration", {}).get("panels"),
        "deepest": load.get("integration", {}).get("deepest"),
        "maximum_estimated_relative_error": load.get("integration", {}).get("maximum_estimated_relative_error"),
        "split": load.get("integration", {}).get("split"),
        "patch_plans": load.get("integration", {}).get("patch_plans"),
    }


def run(binary):
    WORK.mkdir(parents=True, exist_ok=True)
    models = {}
    rows = []
    refusals = []
    for case in build_cases():
        case_root = WORK / case["id"]
        case_root.mkdir(parents=True, exist_ok=True)
        if case["recipe_id"] not in models:
            recipe_path = case_root / "recipe.json"; model_path = case_root / "model.json"
            write_json(recipe_path, case["recipe"])
            built = subprocess.run([str(binary), "build", str(recipe_path), str(model_path)],
                                   capture_output=True, text=True, timeout=180)
            if built.returncode:
                raise RuntimeError(f"{case['id']} model build failed: {built.stdout} {built.stderr}")
            models[case["recipe_id"]] = model_path
        for budget in BUDGETS:
            attempt = case_root / f"panels-{budget}"
            attempt.mkdir(parents=True, exist_ok=True)
            config = {"relative_tolerance": TOLERANCE, "max_depth": DEPTH,
                      "max_panels": budget, "split": "directional"}
            job = surface_oracle.make_job(case["feature"], case["law"], config)
            job_path = attempt / "job.json"; result_path = attempt / "result.json"
            write_json(job_path, job); result_path.write_bytes(SENTINEL)
            started = time.monotonic()
            proc = subprocess.run([str(binary), "solve-job", str(models[case["recipe_id"]]),
                                   str(job_path), str(result_path)], capture_output=True,
                                  text=True, timeout=1800)
            wall = time.monotonic() - started
            diagnostic = (proc.stdout + proc.stderr).strip()
            preserved = result_path.read_bytes() == SENTINEL
            row = {"id": case["id"], "map_direction": case["map_direction"],
                   "budget_panels": budget, "relative_tolerance": TOLERANCE,
                   "max_depth": DEPTH, "law": case["law"],
                   "command_status": "solved" if proc.returncode == 0 else "refused",
                   "wall_s": wall, "diagnostic": diagnostic,
                   "refused_target_preserved": preserved}
            if proc.returncode == 0:
                result = json.loads(result_path.read_text())
                row["scalar_dofs"] = result["scalar_dofs"]
                row["relative_residual"] = result["relative_residual"]
                row["relative_force_balance"] = result["relative_force_balance"]
                row["relative_moment_balance"] = result["relative_moment_balance"]
                row["comparison"] = independent_comparison(result, case["law"])
                plans = row["comparison"].get("patch_plans")
                row["directional_plan_valid"] = (
                    row["comparison"].get("split") == "directional"
                    and isinstance(plans, list)
                    and len(plans) == row["comparison"]["faces"]
                    and all(isinstance(p.get("depths"), list) and len(p["depths"]) == 2
                            and all(isinstance(v, int) and v >= 0 for v in p["depths"])
                            and sum(p["depths"]) > 0 for p in plans))
                row["qualified"] = (all(row["comparison"][k] <= 1e-8 for k in (
                    "q96_q128_force_delta", "area_relative_error", "full_nodal_force_relative_error",
                    "resultant_relative_error", "moment_relative_error", "work_relative_error"))
                    and row["relative_residual"] <= 1e-10
                    and row["directional_plan_valid"])
                row["refused_target_preserved"] = None
            else:
                row["qualified"] = False
                row["refused"] = True
                row["budget_refusal_verified"] = preserved and "budget exceeded" in diagnostic
                refusals.append(row.copy())
            rows.append(row)
            write_json(REPORT, {"protocol": "M12 directional surface load oracle v1",
                                "status": "running", "cases": rows, "refusals": refusals})
            print(case["id"], budget, row["command_status"], row.get("qualified"), flush=True)
    return rows, refusals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=ROOT / "target/release/solidcraft-trivariate")
    args = parser.parse_args()
    rows, refusals = run(args.binary.resolve())
    high_by_case = {row["id"]: row for row in rows if row["budget_panels"] == BUDGETS[-1]}
    low_refusals_valid = all(row.get("budget_refusal_verified") for row in refusals
                             if row["budget_panels"] == BUDGETS[0])
    directional_axis_counts = [0, 0]
    for row in high_by_case.values():
        for plan in (row.get("comparison") or {}).get("patch_plans") or []:
            for axis, depth in enumerate(plan["depths"]):
                if depth:
                    directional_axis_counts[axis] += 1
    axis_coverage = [i for i, count in enumerate(directional_axis_counts) if count]
    corpus_qualified = (len(high_by_case) == len(build_cases())
                        and all(row.get("qualified", False) for row in high_by_case.values())
                        and low_refusals_valid and axis_coverage == [0, 1])
    report = {"protocol": "M12 directional surface load oracle v1",
              "split_mode": "directional", "relative_tolerance": TOLERANCE,
              "max_depth": DEPTH, "budgets_panels": list(BUDGETS),
              "budget_refusals_verified": low_refusals_valid,
              "directional_axis_patch_counts": directional_axis_counts,
              "directional_axis_coverage": axis_coverage,
              "high_budget_case_coverage": sum(bool(r.get("qualified")) for r in high_by_case.values()),
              "high_budget_case_count": len(high_by_case),
              "cases": rows, "refusals": refusals,
              "status": "qualified" if corpus_qualified else "not_qualified",
              "limits": ["directional quadrature estimates are diagnostics, not mathematical error bounds",
                         "q128 reference is qualified by q96/q128 agreement, not a formal certificate"]}
    write_json(REPORT, report)


if __name__ == "__main__":
    main()
