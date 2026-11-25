#!/usr/bin/env python3
"""Independent finite-difference oracle for frozen-plan surface derivatives."""

from __future__ import annotations

import argparse
import copy
from fractions import Fraction
import json
import math
import subprocess
import sys
import time
import platform
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import adaptive_surface_load as surface_oracle
import affine_job_derivative as affine_oracle

CORPUS = ROOT / "tests/references/circular_polygon_corpus.json"
WORK = ROOT / "build/frozen-surface-derivative"
REPORT = ROOT / "docs/results/frozen-surface-derivative-20261125.json"
STEPS = (Fraction(1, 16), Fraction(1, 32), Fraction(1, 64))
ABS_TOL = 1e-9
REL_TOL = 1e-5
ADJOINT_TOL = 1e-8


def recipes_and_cases():
    corpus = {x["id"]: x["recipe"] for x in json.loads(CORPUS.read_text())}
    return [
        {"id": "rectangle_vector_isotropic", "recipe_id": "one_round",
         "recipe": corpus["one_round"], "feature": "edge/outer/v00/v01",
         "law": {"vector_mpa": [0.2, 0.4, -0.1]}, "split": "isotropic",
         "direction": [[0, 0, 1], [0, 0, 0], [0, 0, 0]]},
        {"id": "inclined_bore_pressure_directional", "recipe_id": "extreme_shear",
         "recipe": corpus["extreme_shear"], "feature": "bore/micro",
         "law": {"pressure_mpa": 0.05}, "split": "directional",
         "direction": [[0, 0, 0], [0, 0, 1], [0, 0, 0]]},
        {"id": "oblique_bore_stress_isotropic", "recipe_id": "oblique_profile",
         "recipe": corpus["oblique_profile"], "feature": "bore/left",
         "law": {"stress_mpa": [[0, .3, .1], [.3, 0, -.2], [.1, -.2, 0]]},
         "split": "isotropic", "direction": [[0, 0, 1], [0, 0, 0], [0, 0, 0]]},
    ]


def write_json(path, value):
    surface_oracle.write_json(path, value)


def job_for(case, *, freeze=True, budget=262144):
    config = {"relative_tolerance": 1e-9, "max_depth": 8,
              "max_panels": budget, "split": case["split"]}
    if freeze:
        config["derivative"] = "freeze_primal"
    return surface_oracle.make_job(case["feature"], case["law"], config)


def run_cli(binary, *args, timeout=1800):
    return subprocess.run([str(binary), *map(str, args)], capture_output=True,
                          text=True, timeout=timeout)


def perturbation_recipe(source, matrix_rate, step):
    return affine_oracle.image(source, 1, matrix_rate, [0, 0, 0], step)


def accepted_plan(result):
    return [load.get("integration", {}).get("patch_plans") for load in result.get("loads", [])]


def error_measures(finite_difference, derivative):
    fd = np.asarray(finite_difference, dtype=float)
    exact = np.asarray(derivative, dtype=float)
    absolute = float(np.linalg.norm(fd-exact))
    relative = absolute / max(float(np.linalg.norm(exact)), 1e-30)
    return {"absolute": absolute, "relative": relative,
            "accepted": absolute <= ABS_TOL or relative <= REL_TOL}


def finite_difference_pair(plus, minus, denominator):
    return (np.asarray(plus, dtype=float) - np.asarray(minus, dtype=float)) / (2.0 * float(denominator))


def run_case(binary, case):
    folder = WORK / case["id"]; folder.mkdir(parents=True, exist_ok=True)
    source_path = folder / "source-recipe.json"; model_path = folder / "base-model.json"
    write_json(source_path, case["recipe"])
    built = run_cli(binary, "build", source_path, model_path)
    if built.returncode:
        raise RuntimeError(f"{case['id']} base build refused: {built.stdout} {built.stderr}")
    model = json.loads(model_path.read_text())
    job = job_for(case)
    job_path = folder / "job.json"; write_json(job_path, job)
    primal_path = folder / "primal.result.json"
    primal_run = run_cli(binary, "solve-job", model_path, job_path, primal_path)
    if primal_run.returncode:
        raise RuntimeError(f"{case['id']} primal solve failed: {primal_run.stdout} {primal_run.stderr}")
    primal = json.loads(primal_path.read_text())
    direction_path = folder / "direction.json"
    direction = {"schema": "solidcraft_affine_velocity_v1",
                 "linear": case["direction"], "offset_mm": [0, 0, 0]}
    write_json(direction_path, direction)
    derivative_path = folder / "derivative.result.json"
    derivative_start = time.monotonic()
    derivative_run = run_cli(binary, "differentiate-job", model_path, job_path,
                             direction_path, derivative_path)
    derivative_wall = time.monotonic() - derivative_start
    if derivative_run.returncode:
        raise RuntimeError(f"{case['id']} derivative failed: {derivative_run.stdout} {derivative_run.stderr}")
    derivative = json.loads(derivative_path.read_text())
    if derivative.get("surface_derivative_plan") != "freeze_primal":
        raise AssertionError(f"{case['id']} did not report frozen primal surface plan")
    primal_plans = accepted_plan(primal)
    derivative_surface_loads = derivative.get("surface_loads", [])
    replay_plans = [load.get("integration", {}).get("patch_plans") for load in derivative_surface_loads]
    same_plan = primal_plans == replay_plans
    primal_panels = sum(int(load.get("integration", {}).get("panels", 0)) for load in primal.get("loads", []))
    replay_panels = sum(int(load.get("integration", {}).get("replay_panels", 0)) for load in derivative_surface_loads)
    job_budget = int(job["surface_quadrature"]["max_panels"])
    if not same_plan or primal_panels + replay_panels > job_budget:
        raise AssertionError(f"{case['id']} frozen plan/budget mismatch: {primal_panels}+{replay_panels}>{job_budget}")

    rows = []
    for step in STEPS:
        tagged = str(step).replace("/", "-")
        artifacts = {}
        for sign, signed_step in (("plus", step), ("minus", -step)):
            recipe = perturbation_recipe(case["recipe"], case["direction"], signed_step)
            recipe_path = folder / f"{tagged}-{sign}.recipe.json"
            perturbed_model_path = folder / f"{tagged}-{sign}.model.json"
            perturbed_result_path = folder / f"{tagged}-{sign}.result.json"
            write_json(recipe_path, recipe)
            build = run_cli(binary, "build", recipe_path, perturbed_model_path)
            if build.returncode:
                raise RuntimeError(f"{case['id']} {step} {sign} build failed: {build.stdout}")
            solve = run_cli(binary, "solve-job", perturbed_model_path, job_path, perturbed_result_path)
            if solve.returncode:
                raise RuntimeError(f"{case['id']} {step} {sign} solve failed: {solve.stdout}")
            artifacts[sign] = json.loads(perturbed_result_path.read_text())
        plus, minus = artifacts["plus"], artifacts["minus"]
        h = float(step)
        fd_u = finite_difference_pair(plus["displacement_mm_xyz"], minus["displacement_mm_xyz"], h)
        fd_f = finite_difference_pair(plus["force_n_xyz"], minus["force_n_xyz"], h)
        fd_work = (float(plus["external_work_n_mm"])-float(minus["external_work_n_mm"]))/(2*h)
        plan_plus, plan_minus = accepted_plan(plus), accepted_plan(minus)
        rows.append({"step": str(step), "displacement_rate_fd": fd_u.tolist(),
                     "force_rate_fd": fd_f.tolist(), "work_rate_fd_n_mm": fd_work,
                     "plus_plan_changed_faces": sum(a != b for group_a, group_b in zip(primal_plans, plan_plus)
                                                     for a, b in zip(group_a or [], group_b or [])),
                     "minus_plan_changed_faces": sum(a != b for group_a, group_b in zip(primal_plans, plan_minus)
                                                      for a, b in zip(group_a or [], group_b or [])),
                     "plus_panels": sum(int(x.get("integration", {}).get("panels", 0)) for x in plus.get("loads", [])),
                     "minus_panels": sum(int(x.get("integration", {}).get("panels", 0)) for x in minus.get("loads", []))})
    extrapolated = []
    for i in (0, 1):
        a = rows[i]["displacement_rate_fd"]; b = rows[i+1]["displacement_rate_fd"]
        fa = rows[i]["force_rate_fd"]; fb = rows[i+1]["force_rate_fd"]
        wa = rows[i]["work_rate_fd_n_mm"]; wb = rows[i+1]["work_rate_fd_n_mm"]
        extrapolated.append({"displacement": (4*np.asarray(b)-np.asarray(a))/3,
                             "force": (4*np.asarray(fb)-np.asarray(fa))/3,
                             "work": (4*wb-wa)/3})
    reported_u = np.asarray(derivative["physical_displacement_rate_mm_xyz"])
    reported_f = np.asarray(derivative["nodal_force_rate_n_xyz"])
    reported_work = float(derivative["compliance_rate_tangent_n_mm"])
    checks = []
    for estimate in extrapolated:
        checks.append({"displacement": error_measures(estimate["displacement"], reported_u),
                       "force": error_measures(estimate["force"], reported_f),
                       "work": error_measures([estimate["work"]], [reported_work])})
    adjoint_relative = abs(reported_work-float(derivative["compliance_rate_adjoint_n_mm"])) / max(abs(reported_work), 1.0)
    qualified = (same_plan and primal_panels+replay_panels <= job_budget
                 and adjoint_relative <= ADJOINT_TOL
                 and all(checks[1][key]["accepted"] for key in ("displacement", "force", "work")))
    step_summary = [{"step": item["step"],
                     "displacement_rate_fd_norm_mm_per_unit": float(np.linalg.norm(item["displacement_rate_fd"])),
                     "nodal_force_rate_fd_norm_n_per_unit": float(np.linalg.norm(item["force_rate_fd"])),
                     "work_rate_fd_n_mm": item["work_rate_fd_n_mm"],
                     "plus_plan_changed_faces": item["plus_plan_changed_faces"],
                     "minus_plan_changed_faces": item["minus_plan_changed_faces"],
                     "plus_panels": item["plus_panels"], "minus_panels": item["minus_panels"]}
                    for item in rows]
    extrapolated_summary = [{"from_steps": [str(STEPS[i]), str(STEPS[i+1])],
                             "displacement_rate_norm_mm_per_unit": float(np.linalg.norm(item["displacement"])),
                             "nodal_force_rate_norm_n_per_unit": float(np.linalg.norm(item["force"])),
                             "work_rate_n_mm": float(item["work"])}
                            for i, item in enumerate(extrapolated)]
    return {"id": case["id"], "recipe_id": case["recipe_id"], "feature": case["feature"],
            "law": case["law"], "split": case["split"], "direction": direction,
            "qualified": qualified, "primal_cells": len(primal["analysis_blocks"]),
            "scalar_dofs": primal["scalar_dofs"], "primal_residual": primal["relative_residual"],
            "derivative_residuals": derivative.get("relative_residuals"),
            "primal_panels": primal_panels, "replay_panels": replay_panels,
            "total_panels": primal_panels+replay_panels, "panel_budget": job_budget,
            "frozen_plan_matches_primal": same_plan,
            "derivative_wall_s": derivative_wall, "fd_steps": step_summary,
            "richardson_extrapolations": extrapolated_summary,
            "errors_vs_extrapolation": checks, "tangent_adjoint_relative_error": adjoint_relative}


def refusal_checks(binary, case, model_path, base_job, primal, derivative):
    folder = WORK / "refusals"; folder.mkdir(parents=True, exist_ok=True)
    direction_path = folder / "direction.json"
    write_json(direction_path, {"schema": "solidcraft_affine_velocity_v1",
                                "linear": case["direction"], "offset_mm": [0, 0, 0]})
    checks = []
    jobs = []
    no_opt = copy.deepcopy(base_job); del no_opt["surface_quadrature"]["derivative"]
    jobs.append(("missing_opt_in", no_opt, None))
    unknown = copy.deepcopy(base_job); unknown["surface_quadrature"]["derivative"] = "free_reselect"
    jobs.append(("unknown_mode", unknown, None))
    panels = sum(int(load.get("integration", {}).get("panels", 0)) for load in primal["loads"])
    replay_panels = sum(int(load.get("integration", {}).get("replay_panels", 0))
                        for load in derivative.get("surface_loads", []))
    replay_budget = max(panels, panels + replay_panels - 1)
    budget_job = copy.deepcopy(base_job); budget_job["surface_quadrature"]["max_panels"] = replay_budget
    jobs.append(("replay_budget", budget_job, replay_budget))
    for name, job, budget in jobs:
        job_path = folder / f"{name}.job.json"; output = folder / f"{name}.result.json"
        write_json(job_path, job); marker = b"previous accepted derivative\n"; output.write_bytes(marker)
        run = run_cli(binary, "differentiate-job", model_path, job_path, direction_path, output)
        preserved = output.read_bytes() == marker
        diagnostic = (run.stdout+run.stderr).strip()
        row = {"id": name, "max_panels": budget, "exit_code": run.returncode,
               "diagnostic": diagnostic, "previous_output_preserved": preserved}
        if run.returncode == 0 or not preserved:
            raise AssertionError(f"S4 refusal/atomicity check failed: {row}")
        checks.append(row)
    return checks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=ROOT / "target/release/solidcraft-trivariate")
    args = parser.parse_args(); binary = args.binary.resolve()
    WORK.mkdir(parents=True, exist_ok=True)
    rows = []
    for case in recipes_and_cases():
        rows.append(run_case(binary, case))
        write_json(REPORT, {"protocol": "S4 frozen surface derivative oracle v1",
                            "status": "running", "cases": rows})
    base_case = recipes_and_cases()[1]
    base_folder = WORK / base_case["id"]
    base_job = job_for(base_case)
    base_job_path = base_folder / "job.json"
    model_path = base_folder / "base-model.json"
    primal = json.loads((base_folder / "primal.result.json").read_text())
    derivative = json.loads((base_folder / "derivative.result.json").read_text())
    refusals = refusal_checks(binary, base_case, model_path, base_job, primal, derivative)
    report = {"protocol": "S4 frozen surface derivative oracle v1",
              "revision": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                         capture_output=True, text=True, check=True).stdout.strip(),
              "environment": {"python": platform.python_version(), "numpy": np.__version__},
              "commands": ["rtk python3 tests/references/frozen_surface_derivative.py",
                           "target/release/solidcraft-trivariate build",
                           "target/release/solidcraft-trivariate solve-job",
                           "target/release/solidcraft-trivariate differentiate-job"],
              "steps": [str(x) for x in STEPS], "relative_tolerance": REL_TOL,
              "absolute_tolerance": ABS_TOL, "adjoint_tolerance": ADJOINT_TOL,
              "cases": rows, "refusals": refusals,
              "status": "qualified" if all(row["qualified"] for row in rows)
                       and all(row["previous_output_preserved"] for row in refusals)
                       else "not_qualified",
              "limits": ["finite differences compare perturbed geometries whose adaptively selected plans may change",
                         "the sensitivity is that of the primal frozen quadrature plan, not an error bound for the continuous derivative"]}
    write_json(REPORT, report)
    print(json.dumps({"status": report["status"], "cases": rows, "refusals": refusals}, indent=2))


if __name__ == "__main__":
    main()
