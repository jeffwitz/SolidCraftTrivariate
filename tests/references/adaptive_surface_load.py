#!/usr/bin/env python3
"""Independent surface-load integration oracle for mechanical job v5.

The reference evaluates tensor-product rational Bernstein maps directly from
serialized homogeneous controls. It does not call the Rust load assembler.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "tests/references/circular_polygon_corpus.json"
WORK = ROOT / "build/adaptive-surface-load"
REPORT = ROOT / "docs/results/adaptive-surface-load-20261125.json"


def bernstein(degree, x):
    values = np.array([math.comb(degree, i) * x**i * (1-x)**(degree-i)
                       for i in range(degree+1)], dtype=float)
    if degree == 0:
        return values, np.zeros(1)
    lower = np.array([math.comb(degree-1, i) * x**i * (1-x)**(degree-1-i)
                      for i in range(degree)], dtype=float)
    derivative = np.zeros(degree+1)
    for i in range(degree+1):
        derivative[i] = degree * ((lower[i-1] if i else 0.0) - (lower[i] if i < degree else 0.0))
    return values, derivative


def map_sample(block, uvw):
    degrees = [int(x) for x in block["degrees"]]
    controls = np.asarray(block["controls"], dtype=float)
    expected = math.prod(d+1 for d in degrees)
    if len(controls) != expected or controls.shape[1] != 4:
        raise ValueError(f"malformed control net in {block.get('id')}")
    bases, derivatives = zip(*(bernstein(d, float(x)) for d, x in zip(degrees, uvw)))
    phi = []
    dphi = []
    for i in range(degrees[0]+1):
        for j in range(degrees[1]+1):
            for k in range(degrees[2]+1):
                p = bases[0][i] * bases[1][j] * bases[2][k]
                phi.append(p)
                dphi.append([derivatives[0][i]*bases[1][j]*bases[2][k],
                             bases[0][i]*derivatives[1][j]*bases[2][k],
                             bases[0][i]*bases[1][j]*derivatives[2][k]])
    phi = np.asarray(phi); dphi = np.asarray(dphi)
    weights = controls[:, 3]
    if np.any(weights <= 0):
        raise ValueError("nonpositive homogeneous weight")
    den = phi @ weights
    shape = phi * weights / den
    dden = dphi.T @ weights
    dshape = (dphi * weights[:, None] - shape[:, None] * dden[None, :]) / den
    points = controls[:, :3] / weights[:, None]
    position = shape @ points
    jacobian = points.T @ dshape
    return shape, position, jacobian


def face_sign(axis, upper):
    # Parametric outward orientation for increasing free axes.
    return {(0, 0): -1, (0, 1): 1, (1, 0): 1,
            (1, 1): -1, (2, 0): -1, (2, 1): 1}[(axis, upper)]


def integrate_face(block, side, order, law=None):
    axis = "uvw".index(side[0]); upper = int(side[1])
    free = [q for q in range(3) if q != axis]
    x, weights = np.polynomial.legendre.leggauss(order)
    nodes = (x + 1.0) * 0.5; weights = weights * 0.5
    nd = math.prod(int(d)+1 for d in block["degrees"])
    nodal_force = np.zeros((nd, 3)); area = 0.0
    resultant = np.zeros(3); moment = np.zeros(3)
    absolute_force = 0.0; absolute_moment = 0.0
    for ia, a in enumerate(nodes):
        for ib, b in enumerate(nodes):
            uvw = [0.0, 0.0, 0.0]; uvw[axis] = float(upper)
            uvw[free[0]] = float(a); uvw[free[1]] = float(b)
            shape, point, jac = map_sample(block, uvw)
            raw_normal = face_sign(axis, upper) * np.cross(jac[:, free[0]], jac[:, free[1]])
            jacobian = float(np.linalg.norm(raw_normal))
            measure = jacobian * weights[ia] * weights[ib]
            area += measure
            if law is None:
                continue
            normal = raw_normal / jacobian
            if "vector_mpa" in law:
                traction = np.asarray(law["vector_mpa"], dtype=float)
            elif "pressure_mpa" in law:
                traction = -float(law["pressure_mpa"]) * normal
            elif "stress_mpa" in law:
                traction = np.asarray(law["stress_mpa"], dtype=float) @ normal
            else:
                raise ValueError(f"unsupported surface load law {law}")
            nodal_force += shape[:, None] * (traction * measure)[None, :]
            resultant += traction * measure
            moment += np.cross(point, traction) * measure
            absolute_force += float(np.linalg.norm(traction)) * measure
            absolute_moment += float(np.linalg.norm(np.cross(point, traction))) * measure
    return {"area_mm2": area, "nodal_force_n": nodal_force,
            "resultant_n": resultant, "moment_n_mm": moment,
            "absolute_force_scale_n": absolute_force,
            "absolute_moment_scale_n_mm": absolute_moment}


def area_reference(model, feature, order):
    blocks = {b["id"]: b for b in model["blocks"]}
    faces = [f for f in model["design_faces"] if f["feature"] == feature]
    return sum(integrate_face(blocks[f["block"]], f["side"], order)["area_mm2"] for f in faces)


def integrate_body(blocks, specification, order):
    x, weights = np.polynomial.legendre.leggauss(order)
    nodes = (x + 1.0) * 0.5; weights = weights * 0.5
    ndof = max(int(d) for block in blocks for d in block["dofs"]) + 1
    nodal = np.zeros((ndof, 3)); volume = 0.0
    resultant = np.zeros(3); moment = np.zeros(3)
    constant = np.asarray(specification["constant_n_per_mm3"], dtype=float)
    gradient = np.asarray(specification["gradient_n_per_mm4"], dtype=float)
    for block in blocks:
        ids = np.asarray(block["dofs"], dtype=int)
        for iu, u in enumerate(nodes):
            for iv, v in enumerate(nodes):
                for iw, w in enumerate(nodes):
                    shape, point, jac = map_sample(block, (u, v, w))
                    det = float(np.linalg.det(jac))
                    if det <= 0 or not math.isfinite(det):
                        raise ValueError(f"nonpositive independent volume Jacobian in {block['id']}")
                    measure = det * weights[iu] * weights[iv] * weights[iw]
                    density = constant + gradient @ point
                    force = density * measure
                    volume += measure; resultant += force; moment += np.cross(point, force)
                    for n, dof in enumerate(ids):
                        nodal[dof] += shape[n] * force
    return {"force_n_xyz": nodal.reshape(-1), "volume_mm3": volume,
            "resultant_n": resultant, "moment_n_mm": moment}


def compute_q12_reference(binary):
    corpus = json.loads(CORPUS.read_text())
    case = next(x for x in corpus if x["id"] == "extreme_shear")
    WORK.mkdir(parents=True, exist_ok=True)
    recipe = WORK / "g20-extreme-shear.recipe.json"
    model_path = WORK / "g20-extreme-shear.model.json"
    recipe.write_text(json.dumps(case["recipe"], separators=(",", ":")))
    subprocess.run([str(binary), "build", str(recipe), str(model_path)], check=True,
                   capture_output=True, text=True, timeout=180)
    model = json.loads(model_path.read_text())
    start = time.monotonic()
    area_q12 = area_reference(model, "bore/micro", 12)
    area_q24 = area_reference(model, "bore/micro", 24)
    area_q48 = area_reference(model, "bore/micro", 48)
    area_q64 = area_reference(model, "bore/micro", 64)
    area_q96 = area_reference(model, "bore/micro", 96)
    area_q128 = area_reference(model, "bore/micro", 128)
    # q96/q128 agreement supplies the numerical angular-integral reference.
    high = area_q128
    record = {"case": "G20 extreme_shear", "feature": "bore/micro", "face_count": 8,
              "quadrature_order": 12, "area_q12_mm2": area_q12,
              "area_q24_mm2": area_q24, "area_q48_mm2": area_q48,
              "area_q64_mm2": area_q64, "area_q96_mm2": area_q96,
              "area_q128_mm2": area_q128,
              "q12_relative_to_q128": abs(area_q12/high-1.0),
              "q96_relative_to_q128": abs(area_q96/high-1.0),
              "q128_q96_relative": abs(area_q128/area_q96-1.0),
              "reference_method": "independent rational Bernstein surface Jacobian, tensor Gauss-Legendre",
              "reference_wall_s": time.monotonic()-start}
    return model, record


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False, sort_keys=True) + "\n")


def make_job(feature, law, surface=None):
    selector = feature if isinstance(feature, dict) else {"feature": feature}
    return {"schema": "solidcraft_mechanical_job_v5",
            "material": {"young_mpa": 1000, "poisson": 0.25},
            "formulation": "displacement", "order": 4,
            "plan": {"schema": "adaptive_refinement_plan_v1", "max_depth": 0,
                     "max_cells": 1024, "steps": []},
            "supports": [{"select": {"feature": "start"}, "components": [0, 1, 2]}],
            "loads": [{"select": selector, **law}],
            "surface_quadrature": surface or {"relative_tolerance": 1e-9,
                                                "max_depth": 8, "max_panels": 65536}}


def reference_result_load(result, law, order):
    blocks = {b["id"]: b for b in result["analysis_blocks"]}
    ndof = max(int(d) for b in blocks.values() for d in b["dofs"]) + 1
    nodal = np.zeros((ndof, 3)); area = 0.0
    resultant = np.zeros(3); moment = np.zeros(3)
    force_scale = 0.0; moment_scale = 0.0
    faces = result["loads"][0]["faces"]
    for block_id, side in faces:
        block = blocks[block_id]
        local = integrate_face(block, side, order, law)
        ids = np.asarray(block["dofs"], dtype=int)
        for i, dof in enumerate(ids):
            nodal[dof] += local["nodal_force_n"][i]
        area += local["area_mm2"]
        resultant += local["resultant_n"]
        moment += local["moment_n_mm"]
        force_scale += local["absolute_force_scale_n"]
        moment_scale += local["absolute_moment_scale_n_mm"]
    return {"force_n_xyz": nodal.reshape(-1), "area_mm2": area,
            "resultant_n": resultant, "moment_n_mm": moment,
            "absolute_force_scale_n": force_scale,
            "absolute_moment_scale_n_mm": moment_scale,
            "face_count": len(faces)}


def run_cases(binary, body_extension=False):
    corpus = json.loads(CORPUS.read_text())
    recipe_by_id = {x["id"]: x["recipe"] for x in corpus}
    specs = [
        ("right_circle_pressure", "one_round", "bore/pilot", {"pressure_mpa": 0.05}),
        ("right_circle_vector", "one_round", "bore/pilot", {"vector_mpa": [0.2, -0.4, 0.7]}),
        ("oblique_circle_stress", "oblique_profile", "bore/left",
         {"stress_mpa": [[0.0, 0.3, 0.1], [0.3, 0.0, -0.2], [0.1, -0.2, 0.0]]}),
        ("extreme_shear_circle_pressure", "extreme_shear", "bore/micro", {"pressure_mpa": 0.05}),
        ("extreme_shear_circle_vector", "extreme_shear", "bore/micro", {"vector_mpa": [0.2, -0.4, 0.7]}),
        ("extreme_shear_circle_stress", "extreme_shear", "bore/micro",
         {"stress_mpa": [[0.0, 0.3, 0.1], [0.3, 0.0, -0.2], [0.1, -0.2, 0.0]]}),
        ("extreme_shear_rectangle_vector", "extreme_shear", "edge/outer/v00/v01",
         {"vector_mpa": [0.2, 0.4, -0.1]}),
        ("extreme_shear_rectangle_body_and_surface", "extreme_shear", "edge/outer/v00/v01",
         {"vector_mpa": [0.2, 0.4, -0.1]}),
    ]
    models = {}
    rows = []
    for label, recipe_id, feature, law in specs:
        if label.endswith("body_and_surface") != body_extension:
            continue
        folder = WORK / label; folder.mkdir(parents=True, exist_ok=True)
        if recipe_id not in models:
            recipe_path = folder / "recipe.json"; model_path = folder / "model.json"
            write_json(recipe_path, recipe_by_id[recipe_id])
            built = subprocess.run([str(binary), "build", str(recipe_path), str(model_path)],
                                   capture_output=True, text=True, timeout=180)
            if built.returncode:
                raise RuntimeError(f"{label}: build failed: {built.stdout} {built.stderr}")
            models[recipe_id] = model_path
        surface = ({"relative_tolerance": 1e-9, "max_depth": 8, "max_panels": 262144}
                   if recipe_id in ("extreme_shear", "oblique_profile") else None)
        job = make_job(feature, law, surface)
        body_specification = None
        if label == "extreme_shear_rectangle_body_and_surface":
            body_specification = {"constant_n_per_mm3": [0.01, -0.02, 0.005],
                                  "gradient_n_per_mm4": [[0.001, 0.0, 0.0],
                                                         [0.0, -0.001, 0.0],
                                                         [0.0, 0.0, 0.0005]]}
            job["body_load"] = body_specification
        job_path = folder / "job.json"; result_path = folder / "result.json"
        write_json(job_path, job)
        started = time.monotonic()
        solved = subprocess.run([str(binary), "solve-job", str(models[recipe_id]),
                                 str(job_path), str(result_path)], capture_output=True,
                                text=True, timeout=1800)
        row = {"id": label, "recipe": recipe_id, "load": {"select": feature, **law},
               "command_status": "solved" if solved.returncode == 0 else "refused",
               "wall_s": time.monotonic()-started, "diagnostic": (solved.stdout+solved.stderr)[-2000:]}
        if solved.returncode:
            rows.append(row); continue
        result = json.loads(result_path.read_text())
        q96 = reference_result_load(result, law, 96)
        q128 = reference_result_load(result, law, 128)
        rust_load = result["loads"][0]
        rust_force = np.asarray(result["force_n_xyz"], dtype=float)
        body_result = None
        if body_specification is not None:
            body_result = integrate_body(result["analysis_blocks"], body_specification, 16)
        ref_force = q128["force_n_xyz"] + (body_result["force_n_xyz"] if body_result else 0.0)
        force_scale = max(q128["absolute_force_scale_n"], 1.0)
        moment_scale = max(q128["absolute_moment_scale_n_mm"], 1.0)
        vector_error = float(np.linalg.norm(rust_force-ref_force)/force_scale)
        vector_q_delta = float(np.linalg.norm(q96["force_n_xyz"]-q128["force_n_xyz"])/force_scale)
        area_error = abs(float(rust_load["area_mm2"])/q128["area_mm2"]-1.0)
        resultant_error = float(np.linalg.norm(np.asarray(rust_load["resultant_n"])-q128["resultant_n"])/force_scale)
        moment_error = float(np.linalg.norm(np.asarray(rust_load["moment_n_mm"])-q128["moment_n_mm"])/moment_scale)
        displacement = np.asarray(result["displacement_mm_xyz"], dtype=float)
        reference_work = float(q128["force_n_xyz"] @ displacement)
        work_error = abs(reference_work-float(result["traction_work_n_mm"])) / max(abs(reference_work), 1.0)
        row.update({"cells": len(result["analysis_blocks"]), "scalar_dofs": result["scalar_dofs"],
                    "faces": q128["face_count"], "area_reference_q128_mm2": q128["area_mm2"],
                    "area_rust_mm2": rust_load["area_mm2"], "area_relative_error": area_error,
                    "full_nodal_force_relative_error": vector_error,
                    "independent_q96_q128_force_delta": vector_q_delta,
                    "resultant_relative_error": resultant_error, "moment_relative_error": moment_error,
                    "independent_traction_work_n_mm": reference_work,
                    "traction_work_relative_error": work_error,
                    "rust_quadrature": rust_load.get("integration"),
                    "relative_residual": result["relative_residual"],
                    "relative_force_balance": result["relative_force_balance"],
                    "relative_moment_balance": result["relative_moment_balance"]})
        if body_result:
            body_report = result["body_load"]
            row["body_load"] = {"independent_volume_mm3": body_result["volume_mm3"],
                                 "reported_volume_mm3": body_report["volume_mm3"],
                                 "independent_resultant_n": body_result["resultant_n"].tolist(),
                                 "reported_resultant_n": body_report["resultant_n"],
                                 "independent_moment_n_mm": body_result["moment_n_mm"].tolist(),
                                 "reported_moment_n_mm": body_report["moment_n_mm"],
                                 "body_work_n_mm": body_report["work_n_mm"]}
            row["body_load"]["volume_relative_error"] = abs(body_report["volume_mm3"]/body_result["volume_mm3"]-1)
            row["body_load"]["resultant_relative_error"] = float(np.linalg.norm(np.asarray(body_report["resultant_n"])-body_result["resultant_n"])/max(np.linalg.norm(body_result["resultant_n"]),1.0))
            row["body_load"]["moment_relative_error"] = float(np.linalg.norm(np.asarray(body_report["moment_n_mm"])-body_result["moment_n_mm"])/max(np.linalg.norm(body_result["moment_n_mm"]),1.0))
            row["body_load"]["nodal_body_work_error"] = abs(body_report["work_n_mm"] - float(body_result["force_n_xyz"] @ displacement)) / max(1., abs(body_report["work_n_mm"]))
        row["qualified"] = all(row[key] <= 1e-8 for key in (
            "area_relative_error", "full_nodal_force_relative_error", "resultant_relative_error",
            "moment_relative_error", "traction_work_relative_error", "independent_q96_q128_force_delta"))
        if body_result:
            row["qualified"] &= max(row["body_load"][k] for k in ("volume_relative_error", "resultant_relative_error", "moment_relative_error", "nodal_body_work_error")) <= 1e-8
        rows.append(row)
        write_json(REPORT, {"protocol": "M11 adaptive surface load oracle v1", "cases": rows,
                            "status": "running", "limits": ["surface quadrature convergence is an estimator, not a mathematical error bound"]})
        print(label, row["qualified"], row.get("full_nodal_force_relative_error"), flush=True)
    return rows


def refusal_checks(binary):
    model_path = WORK / "extreme_shear_circle_pressure" / "model.json"
    model = json.loads(model_path.read_text())
    job = make_job("bore/micro", {"pressure_mpa": 0.05})
    direction_path = WORK / "identity-velocity.json"
    write_json(direction_path, {"schema": "solidcraft_affine_velocity_v1",
                                "linear": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                                "offset_mm": [0, 0, 0]})
    checks = []
    for name, mutate, derivative in (
        ("panel_budget", lambda j: j["surface_quadrature"].update(max_panels=1), False),
        ("invalid_budget", lambda j: j["surface_quadrature"].update(max_depth=0), False),
        ("shape_derivative_fixed_plan", lambda _j: None, True),
    ):
        candidate = json.loads(json.dumps(job)); mutate(candidate)
        candidate_path = WORK / f"refusal-{name}.job.json"
        output_path = WORK / f"refusal-{name}.result.json"
        marker = b"previous accepted result\n"; output_path.write_bytes(marker)
        write_json(candidate_path, candidate)
        if derivative:
            command = [str(binary), "differentiate-job", str(model_path), str(candidate_path),
                       str(direction_path), str(output_path)]
        else:
            command = [str(binary), "solve-job", str(model_path), str(candidate_path), str(output_path)]
        run = subprocess.run(command, capture_output=True, text=True, timeout=180)
        preserved = output_path.read_bytes() == marker
        diagnostic = (run.stdout + run.stderr).strip()
        record = {"id": name, "exit_code": run.returncode, "diagnostic": diagnostic,
                  "previous_output_preserved": preserved}
        if not preserved or run.returncode == 0:
            raise AssertionError(f"refusal/atomicity check failed: {record}")
        checks.append(record)
    return checks


def measured_panel_budgets(binary,rows):
    """Rerun each original job without changing its geometry, law or tolerance."""
    frozen = {r["id"]:r["recipe"] for r in json.loads(CORPUS.read_text())}
    trials=[]
    for row in rows:
        if "body_load" in row or row["id"].endswith("body_and_surface"):
            continue
        folder=WORK/"budget-replay"/row["id"];folder.mkdir(exist_ok=True,parents=True)
        rp,mp,jp,out=[folder/(x+".json") for x in ("recipe","model","job","result")]
        write_json(rp,frozen[row["recipe"]])
        subprocess.run([str(binary),"build",str(rp),str(mp)],capture_output=True,check=True)
        baseline=json.loads((WORK/row["id"]/"result.json").read_text())
        for strategy,budget in (("isotropic",65536),("directional",4096),("directional",65536)):
            job=json.loads((WORK/row["id"]/"job.json").read_text())
            job["surface_quadrature"].update(split=strategy,max_panels=budget)
            write_json(jp,job)
            marker=b"previous qualified result";out.write_bytes(marker)
            start=time.monotonic()
            proc=subprocess.run([str(binary),"solve-job",str(mp),str(jp),str(out)],capture_output=True,text=True)
            item=dict(id=row["id"],split=strategy,max_panels=budget,
                      elapsed_s=time.monotonic()-start,exit_code=proc.returncode,
                      command_status="refused" if proc.returncode else "solved",
                      diagnostic=proc.stdout.strip(),job=job)
            if proc.returncode:
                item["target_preserved"]=out.read_bytes()==marker
            else:
                result=json.loads(out.read_text())
                lhs=np.asarray(result["force_n_xyz"]);rhs=np.asarray(baseline["force_n_xyz"])
                delta=float(np.linalg.norm(lhs-rhs)/max(1.,np.linalg.norm(rhs)))
                item.update(panels=result["loads"][0]["integration"]["panels"],force_delta_to_qualified_isotropic=delta,
                            full_field_relative=float(np.linalg.norm(np.asarray(result["displacement_mm_xyz"])-np.asarray(baseline["displacement_mm_xyz"]))/max(np.linalg.norm(baseline["displacement_mm_xyz"]),1e-30)),
                            reproduced_force=delta<1e-8)
            trials.append(item)
    return trials


def main():
    global REPORT
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=ROOT / "target/release/solidcraft-trivariate")
    parser.add_argument("--reference-only", action="store_true")
    parser.add_argument("--output", type=Path, default=REPORT)
    parser.add_argument("--body-extension", action="store_true")
    parser.add_argument("--budget-replay-only", action="store_true")
    args = parser.parse_args()
    REPORT = args.output
    binary = args.binary.resolve()
    if args.budget_replay_only:
        previous=json.loads((ROOT/"docs/results/adaptive-surface-load-20261125.json").read_text())
        trials=measured_panel_budgets(binary,previous["cases"])
        valid=all((t["command_status"]=="solved" and t["reproduced_force"] and t["full_field_relative"]<1e-8)
                  or (t["split"]=="isotropic" and t.get("target_preserved") and "panel budget exceeded" in t["diagnostic"])
                  for t in trials)
        write_json(REPORT,dict(schema="matched_surface_strategy_replay_v1",cases=trials,qualified=valid))
        if not valid: raise SystemExit(1)
        return
    _model, q12 = compute_q12_reference(binary)
    rows = [] if args.reference_only else run_cases(binary,args.body_extension)
    refusals = [] if args.reference_only or args.body_extension else refusal_checks(binary)
    report = {"protocol": "M11 adaptive surface load oracle v1",
              "independent_reference": q12, "cases": rows, "refusals": refusals,
              "budget_trials": [] if args.reference_only else measured_panel_budgets(binary,rows),
              "status": "reference_only" if args.reference_only else
                       ("qualified" if all(x.get("qualified", False) for x in rows) and
                        all(x["previous_output_preserved"] for x in refusals) else "not_qualified"),
              "limits": ["surface quadrature convergence is an estimator, not a mathematical error bound",
                         "the q96/q128 difference is the numerical reference diagnostic, not a certified bound"]}
    write_json(REPORT, report)
    print(json.dumps(q12, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
