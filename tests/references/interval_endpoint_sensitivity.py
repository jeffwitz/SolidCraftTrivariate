#!/usr/bin/env python3
"""Independent oracle for G21 material interval-endpoint sensitivities."""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
import re
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HEIGHT = 64
LOWER, UPPER = 16, 48
RADIUS = 6
TOP_SHIFT = np.array([12.0, -8.0])
CENTER = np.array([2.0, -1.0])
AREA = math.pi * RADIUS**2
OUTER_AREA = 48 * 32
PLAN = {"schema": "adaptive_refinement_plan_v1", "max_depth": 0,
        "max_cells": 4096, "steps": []}


def write(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def run_cli(cli: Path, *args: object, expect_success: bool = True) -> subprocess.CompletedProcess:
    proc = subprocess.run([str(cli), *map(str, args)], capture_output=True, text=True)
    if expect_success and proc.returncode:
        raise RuntimeError(f"{[str(cli), *map(str, args)]} failed ({proc.returncode}): {proc.stdout}\n{proc.stderr}")
    return proc


def recipe(lo: int = LOWER, hi: int = UPPER, height: int = HEIGHT,
           top_shift=(12, -8), bores=None, wide=False):
    if bores is None:
        bores = [{"id": "pilot", "center_mm": [2, -1], "radius_mm": RADIUS,
                  "guard_mm": 12, "z_interval_mm": [lo, hi]}]
    half_x, half_y = (48, 32) if wide else (24, 16)
    return {
        "family": "polygon_interval_bores_v1",
        "outer": {"id": "outer", "vertices": [
            {"id": "v00", "xy_mm": [-half_x, -half_y]},
            {"id": "v01", "xy_mm": [half_x, -half_y]},
            {"id": "v02", "xy_mm": [half_x, half_y]},
            {"id": "v03", "xy_mm": [-half_x, half_y]},
        ]},
        "bores": bores, "height_mm": height,
        "top_shift_mm": list(top_shift),
    }


def job_for(law: str):
    common = {
        "schema": "solidcraft_mechanical_job_v4" if law == "body_affine" else "solidcraft_mechanical_job_v3",
        "material": {"young_mpa": 1000.0, "poisson": 0.25},
        "formulation": "displacement", "order": 8, "plan": copy.deepcopy(PLAN),
        "supports": [{"select": {"feature": "start"}, "components": [0, 1, 2]}],
    }
    if law == "pressure_lower":
        common["loads"] = [{"select": {"feature": "bore/pilot/lower"}, "pressure_mpa": 0.25}]
    elif law == "wall_vector":
        common["loads"] = [{"select": {"feature": "bore/pilot/wall"},
                            "vector_mpa": [0.04, -0.025, 0.03]}]
    elif law == "body_affine":
        common["loads"] = []
        common["body_load"] = {
            "constant_n_per_mm3": [0.2, -0.15, 0.3],
            "gradient_n_per_mm4": [[0.001, 0.0005, 0.0002],
                                    [0.0002, -0.0008, 0.0001],
                                    [-0.0004, 0.0007, 0.0003]],
        }
    else:
        raise ValueError(law)
    return common


def gauss_legendre(order: int):
    pairs = []
    for i in range(1, (order + 1) // 2 + 1):
        x = math.cos(math.pi * (i - 0.25) / (order + 0.5))
        for _ in range(64):
            p0, p1 = 1.0, x
            for k in range(2, order + 1):
                p0, p1 = p1, ((2*k-1)*x*p1-(k-1)*p0)/k
            derivative = order * (x*p1-p0)/(x*x-1.0)
            step = p1 / derivative
            x -= step
            if abs(step) < 2e-16:
                break
        weight = 2.0 / ((1.0-x*x)*derivative*derivative)
        pairs.extend(((-x, weight), (x, weight)) if 2*i <= order else ((x, weight),))
    return [(0.5*(x+1.0), 0.5*w) for x, w in sorted(pairs)]


def wall_resultant_rate(endpoint: str, traction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Independent area and moment rate of the cylindrical wall endpoint strip."""
    slope = TOP_SHIFT / HEIGHT
    z = float(UPPER if endpoint == "upper" else LOWER)
    sign = 1.0 if endpoint == "upper" else -1.0
    force = np.zeros(3)
    moment = np.zeros(3)
    for theta, weight in gauss_legendre(256):
        theta *= 2.0 * math.pi
        weight *= 2.0 * math.pi
        radial = np.array([math.cos(theta), math.sin(theta)])
        area_rate = RADIUS * math.sqrt(1.0 + float(slope @ radial)**2) * weight
        position = np.array([*(CENTER + slope*z + RADIUS*radial), z])
        df = sign * area_rate * traction
        force += df
        moment += np.cross(position, df)
    return force, moment


def body_resultant_rate(endpoint: str) -> tuple[np.ndarray, np.ndarray]:
    """Exact disk-area integration of the affine density at the moving cap."""
    sign = 1.0 if endpoint == "lower" else -1.0
    z = float(LOWER if endpoint == "lower" else UPPER)
    center = np.r_[CENTER + TOP_SHIFT*z/HEIGHT, z]
    density = np.array([0.2, -0.15, 0.3])
    gradient = np.array([[0.001, 0.0005, 0.0002],
                         [0.0002, -0.0008, 0.0001],
                         [-0.0004, 0.0007, 0.0003]])
    value = density + gradient @ center
    force = sign * AREA * value
    moment = sign * AREA * np.cross(center, value)
    # The affine in-plane density contributes the disk's second moment.
    m2 = AREA * RADIUS**2 / 4.0
    moment += sign * np.array([gradient[2, 1]*m2,
                               -gradient[2, 0]*m2,
                               (gradient[1, 0]-gradient[0, 1])*m2])
    return force, moment


def analytic_rates(spec: dict) -> dict:
    endpoint, law = spec["endpoint"], spec["law"]
    volume = AREA * (1.0 if endpoint == "lower" else -1.0)
    if law == "pressure_lower":
        pressure = 0.25
        force = np.zeros(3)
        f = np.array([0.0, 0.0, -pressure*AREA])
        axis = np.r_[TOP_SHIFT/HEIGHT, 1.0]
        moment = np.cross(axis, f)
    elif law == "wall_vector":
        force, moment = wall_resultant_rate(endpoint, np.array([0.04, -0.025, 0.03]))
    else:
        force, moment = body_resultant_rate(endpoint)
    return {"volume_rate_mm3": volume, "resultant_rate_n": force.tolist(),
            "moment_rate_n_mm": moment.tolist()}


def compliance(result: dict) -> float:
    return float(np.dot(np.asarray(result["force_n_xyz"]),
                        np.asarray(result["displacement_mm_xyz"])))


def root_key(block_id: str) -> str:
    return re.sub(r"_z-?\d+_-?\d+$", "", block_id)


def control_map(result: dict, field: str) -> dict[tuple[str, int, int], np.ndarray]:
    """Remap control values across layer IDs using root identity and z ordering."""
    source = np.asarray(result[field], dtype=float)
    grouped: dict[str, list[tuple[float, dict]]] = {}
    for block in result["analysis_blocks"]:
        controls = np.asarray(block["controls"], dtype=float).reshape(-1, 4)
        z = controls[:, 2] / controls[:, 3]
        grouped.setdefault(root_key(block["id"]), []).append((float(z.mean()), block))
    out = {}
    for ident, entries in grouped.items():
        entries.sort(key=lambda item: item[0])
        for occurrence, (_, block) in enumerate(entries):
            for local_index, dof in enumerate(block["dofs"]):
                out[(ident, occurrence, local_index)] = source[3*dof:3*dof+3]
    return out


def mapped_fd(reference: dict, plus: dict, minus: dict, field: str, step: float):
    ref = control_map(reference, field)
    pos = control_map(plus, field)
    neg = control_map(minus, field)
    common = sorted(set(ref) & set(pos) & set(neg))
    if not common:
        raise RuntimeError(f"no material identity matches for {field}")
    tangent = np.concatenate([ref[key] for key in common])
    centered = np.concatenate([(pos[key]-neg[key])/(2*step) for key in common])
    return centered, tangent, {"matched_control_vectors": len(common),
                               "reference_control_vectors": len(ref),
                               "plus_control_vectors": len(pos),
                               "minus_control_vectors": len(neg)}


def richardson(values: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    # Input order is h=4,2,1; return fourth- and sixth-order estimates.
    d4, d2, d1 = values
    r2 = (4*d2-d4)/3.0
    r1 = (4*d1-d2)/3.0
    return r1, (16*r1-r2)/15.0


def tolerance_fraction(error: float, reference: float, rel: float) -> float:
    return error / max(1e-8, rel*abs(reference))


def solve(cli: Path, folder: Path, name: str, rec: dict, job: dict):
    rp, mp, jp, op = (folder/f"{name}.{ext}" for ext in ("recipe.json", "model.json", "job.json", "solution.json"))
    write(rp, rec); write(jp, job)
    run_cli(cli, "build", rp, mp)
    run_cli(cli, "audit", mp)
    run_cli(cli, "solve-job", mp, jp, op)
    return json.loads(mp.read_text()), json.loads(op.read_text()), mp, jp, op


def case_spec():
    return [
        {"name": "lower_cap_pressure", "endpoint": "lower", "law": "pressure_lower"},
        {"name": "upper_wall_vector", "endpoint": "upper", "law": "wall_vector"},
        {"name": "upper_cavity_body", "endpoint": "upper", "law": "body_affine"},
    ]


def run(cli: Path, work: Path) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    rows=[]
    all_pass=True
    for spec in case_spec():
        folder=work/spec["name"]; folder.mkdir(exist_ok=True)
        job=job_for(spec["law"])
        base_model, base, mp, jp, op=solve(cli,folder,"base",recipe(),job)
        direction={"schema":"solidcraft_interval_endpoint_velocity_v1",
                   "bore":"pilot","endpoint":spec["endpoint"]}
        dp=folder/"direction.json"; derivative_path=folder/"derivative.json"
        write(dp,direction)
        proc=run_cli(cli,"differentiate-job",mp,jp,dp,derivative_path)
        derivative=json.loads(derivative_path.read_text())
        analytic=analytic_rates(spec)
        assert derivative["schema"]=="solidcraft_interval_endpoint_derivative_v1"
        assert derivative["active_planes_mm"]==([0,16,48] if spec["endpoint"]=="lower" else [16,48,64])

        sensitivity={}
        for key, expected in (("volume_rate_mm3",analytic["volume_rate_mm3"]),):
            got=float(derivative[key]); err=abs(got-expected)
            sensitivity[key]={"native":got,"independent":expected,"absolute_error":err,
                              "relative_error":err/max(abs(expected),1e-30)}
        for key, expected_key in (("applied_resultant_rate_n","resultant_rate_n"),):
            got=np.asarray(derivative[key],dtype=float); expected=np.asarray(analytic[expected_key])
            err=float(np.linalg.norm(got-expected)); denom=float(np.linalg.norm(expected))
            sensitivity[key]={"native":got.tolist(),"independent":expected.tolist(),
                              "absolute_error":err,"relative_error":err/max(denom,1e-30)}

        plus_results={}; minus_results={}; plus_models={}; minus_models={}
        compliances=[]; force_rates=[]; displacement_rates=[]; scalar_rates=[]
        volume_fd=[]; resultant_fd=[]; moment_fd=[]; map_records=[]
        for h in (4,2,1):
            lo,hi=LOWER,UPPER
            idx=0 if spec["endpoint"]=="lower" else 1
            plus_int=[lo,hi]; minus_int=[lo,hi]
            plus_int[idx]+=h; minus_int[idx]-=h
            plus_rec=recipe(plus_int[0],plus_int[1]); minus_rec=recipe(minus_int[0],minus_int[1])
            _, plus, plus_mp, plus_jp, plus_op=solve(cli,folder,f"plus_{h}",plus_rec,job)
            _, minus, minus_mp, minus_jp, minus_op=solve(cli,folder,f"minus_{h}",minus_rec,job)
            plus_results[h],minus_results[h]=plus,minus
            plus_models[h],minus_models[h]=plus_mp,minus_mp
            volume_fd.append(np.array([(plus["volume_mm3"]-minus["volume_mm3"])/(2*h)]))
            resultant_fd.append((np.asarray(plus["applied_resultant_n"])-np.asarray(minus["applied_resultant_n"]))/(2*h))
            moment_fd.append((np.asarray(plus["applied_moment_n_mm"])-np.asarray(minus["applied_moment_n_mm"]))/(2*h))
            scalar_rates.append(np.array([(compliance(plus)-compliance(minus))/(2*h)]))
            fd_force, tangent_force, force_ids=mapped_fd(base,plus,minus,"force_n_xyz",h)
            fd_disp, tangent_disp, disp_ids=mapped_fd(base,plus,minus,"displacement_mm_xyz",h)
            force_rates.append(fd_force); displacement_rates.append(fd_disp)
            map_records.append({"step_mm":h,"force_identity_map":force_ids,
                                "displacement_identity_map":disp_ids})

        vol6=richardson(volume_fd)[1][0]
        resultant6=richardson(resultant_fd)[1]
        moment6=richardson(moment_fd)[1]
        work6=richardson(scalar_rates)[1][0]
        force6=richardson(force_rates)[1]
        disp6=richardson(displacement_rates)[1]
        native_force=np.asarray(derivative["nodal_force_rate_n_xyz"],dtype=float)
        native_disp=np.asarray(derivative["physical_displacement_rate_mm_xyz"],dtype=float)
        # Keep reference alignment by the exact block/control map built from
        # the base result; the derivative vector uses that same DOF numbering.
        base_maps={}
        grouped={}
        for block in base["analysis_blocks"]:
            ctrl=np.asarray(block["controls"],dtype=float).reshape(-1,4)
            grouped.setdefault(root_key(block["id"]),[]).append((float((ctrl[:,2]/ctrl[:,3]).mean()),block))
        for ident, entries in grouped.items():
            entries.sort(key=lambda a:a[0])
            for occurrence,(_,block) in enumerate(entries):
                for local,dof in enumerate(block["dofs"]):
                    base_maps[(ident,occurrence,local)]=dof
        force_keys=sorted(control_map(base,"force_n_xyz"))
        disp_keys=sorted(control_map(base,"displacement_mm_xyz"))
        direct_force=np.concatenate([native_force[3*base_maps[k]:3*base_maps[k]+3]
                                     for k in force_keys])
        direct_disp=np.concatenate([native_disp[3*base_maps[k]:3*base_maps[k]+3]
                                    for k in disp_keys])
        force_err=float(np.linalg.norm(force6-direct_force))
        disp_err=float(np.linalg.norm(disp6-direct_disp))
        force_rel=force_err/max(float(np.linalg.norm(direct_force)),1e-30)
        disp_rel=disp_err/max(float(np.linalg.norm(direct_disp)),1e-30)
        tangent=float(derivative["compliance_rate_tangent_n_mm"])
        adjoint=float(derivative["compliance_rate_adjoint_n_mm"])
        adjoint_err=abs(tangent-adjoint)
        adjoint_rel=adjoint_err/max(abs(tangent),1e-30)
        expected_force=np.asarray(analytic["resultant_rate_n"])
        expected_moment=np.asarray(analytic["moment_rate_n_mm"])
        sensitivity.update({
            "resultant_fd_sixth": {"value":resultant6.tolist(),"independent":expected_force.tolist(),
                                    "relative_error":float(np.linalg.norm(resultant6-expected_force)/max(np.linalg.norm(expected_force),1e-30))},
            "moment_fd_sixth": {"value":moment6.tolist(),"independent":expected_moment.tolist(),
                                 "relative_error":float(np.linalg.norm(moment6-expected_moment)/max(np.linalg.norm(expected_moment),1e-30))},
            "volume_fd_sixth_mm3":float(vol6),
            "compliance_rate_tangent_n_mm":tangent,
            "compliance_rate_adjoint_n_mm":adjoint,
            "compliance_fd_sixth_n_mm":float(work6),
            "compliance_fd_tangent_relative_error":abs(work6-tangent)/max(abs(tangent),1e-30),
            "adjoint_absolute_error":adjoint_err,"adjoint_relative_error":adjoint_rel,
            "nodal_force_fd_sixth_relative_error":force_rel,
            "displacement_fd_sixth_relative_error":disp_rel,
            "fd_steps_mm":[4,2,1],"identity_remap":map_records,
        })
        volume_rel=abs(float(derivative["volume_rate_mm3"])-analytic["volume_rate_mm3"])/abs(analytic["volume_rate_mm3"])
        result_rel=float(sensitivity["applied_resultant_rate_n"]["relative_error"])
        # Published output currently exposes volume and applied force rate;
        # moment-rate accuracy is independently checked on cold differences.
        work_fraction=tolerance_fraction(abs(work6-tangent),tangent,1e-4)
        disp_fraction=tolerance_fraction(float(np.linalg.norm(disp6-direct_disp)),float(np.linalg.norm(direct_disp)),1e-4)
        force_fraction=tolerance_fraction(float(np.linalg.norm(force6-direct_force)),float(np.linalg.norm(direct_force)),1e-4)
        adjoint_fraction=adjoint_rel/1e-8
        checks={"volume_relative":volume_rel,"resultant_rate_relative":result_rel,
                "moment_rate_fd_relative":sensitivity["moment_fd_sixth"]["relative_error"],
                "compliance_fd_tolerance_fraction":work_fraction,
                "displacement_fd_tolerance_fraction":disp_fraction,
                "nodal_force_fd_tolerance_fraction":force_fraction,
                "adjoint_tolerance_fraction":adjoint_fraction,
                "compliance_value_matches_solve":abs(tangent-float(derivative["compliance_rate_tangent_n_mm"]))==0.0}
        passed=(volume_rel<=1e-7 and result_rel<=1e-5
                and sensitivity["moment_fd_sixth"]["relative_error"]<=1e-4
                and work_fraction<=1 and disp_fraction<=1 and force_fraction<=1 and adjoint_fraction<=1)
        row={"name":spec["name"],"law":spec["law"],"endpoint":spec["endpoint"],
             "recipe":recipe(),"job":job,"direction":direction,
             "model_blocks":len(base_model["blocks"]),"analysis_blocks":len(base["analysis_blocks"]),
             "scalar_dofs":base["scalar_dofs"],"base_volume_mm3":base["volume_mm3"],
             "base_compliance_n_mm":compliance(base),"derivative":derivative,
             "independent_rates":analytic,"measurements":sensitivity,"checks":checks,
             "passed":bool(passed)}
        write(folder/"measurement.json",row)
        rows.append(row); all_pass &= passed

    refusals=run_refusals(cli,work)
    all_pass &= all(r["refused"] and r["preserved"] for r in refusals)
    return {"schema":"interval_endpoint_sensitivity_oracle_v1",
            "provenance":{"base_revision":"4d5ae4c","checkout_revision":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
                          "python":sys.version.split()[0],"cli":str(cli),
                          "commands":["rtk cargo build --release --locked --offline --bin solidcraft-trivariate",
                                      "rtk python3 tests/references/interval_endpoint_sensitivity.py --cli target/release/solidcraft-trivariate --work build/interval-endpoint-sensitivity --output docs/results/interval-endpoint-sensitivity-20261125.json"]},
            "protocol":{"height_mm":HEIGHT,"interval_mm":[LOWER,UPPER],"radius_mm":RADIUS,
                        "top_shift_mm":TOP_SHIFT.astype(int).tolist(),"finite_difference_steps_mm":[4,2,1],
                        "volume_relative_tolerance":1e-7,"field_and_work_relative_tolerance":1e-4,
                        "absolute_floor":1e-8,"adjoint_relative_tolerance":1e-8},
            "cases":rows,"refusals":refusals,"passed":bool(all_pass)}


def run_refusals(cli: Path, work: Path) -> list[dict]:
    folder=work/"refusals"; folder.mkdir(exist_ok=True)
    base_rec=recipe(); base_job=job_for("pressure_lower")
    rec_path=folder/"base.recipe.json"; model_path=folder/"base.model.json"
    job_path=folder/"base.job.json"; direction_path=folder/"direction.json"
    output=folder/"preserved-output.json"
    write(rec_path,base_rec); write(job_path,base_job)
    write(direction_path,{"schema":"solidcraft_interval_endpoint_velocity_v1","bore":"pilot","endpoint":"lower"})
    run_cli(cli,"build",rec_path,model_path)
    run_cli(cli,"differentiate-job",model_path,job_path,direction_path,output)
    sentinel=output.read_bytes()
    cases=[]
    # Each run starts with a known accepted output and then checks byte preservation.
    invalid_specs=[]
    external=recipe(lo=0,hi=UPPER)
    invalid_specs.append(("external_endpoint",external,base_job,
                          {"schema":"solidcraft_interval_endpoint_velocity_v1","bore":"pilot","endpoint":"lower"},None))
    coincident_bores=[{"id":"pilot","center_mm":[2,-1],"radius_mm":6,"guard_mm":12,"z_interval_mm":[16,48]},
                      {"id":"second","center_mm":[-22,-1],"radius_mm":4,"guard_mm":8,"z_interval_mm":[16,40]}]
    invalid_specs.append(("coincident_endpoint",recipe(bores=coincident_bores,wide=True),base_job,
                          {"schema":"solidcraft_interval_endpoint_velocity_v1","bore":"pilot","endpoint":"lower"},None))
    invalid_specs.append(("unknown_bore",base_rec,base_job,
                          {"schema":"solidcraft_interval_endpoint_velocity_v1","bore":"missing","endpoint":"lower"},model_path))
    nonhom=copy.deepcopy(base_job)
    nonhom["supports"][0]["affine_mm"]={"denominator":1000,"numerators":[[1,0,0,0],[0,0,0,0],[0,0,0,0]]}
    invalid_specs.append(("nonhomogeneous_support",base_rec,nonhom,
                          {"schema":"solidcraft_interval_endpoint_velocity_v1","bore":"pilot","endpoint":"lower"},model_path))
    incompatible={"family":"reference_rectangle_v1"}
    invalid_specs.append(("incompatible_family",incompatible,base_job,
                          {"schema":"solidcraft_interval_endpoint_velocity_v1","bore":"pilot","endpoint":"lower"},None))
    for name,rec,job,spec,use_model in invalid_specs:
        rp=folder/f"{name}.recipe.json"; mp=folder/f"{name}.model.json"
        jp=folder/f"{name}.job.json"; dp=folder/f"{name}.direction.json"
        write(rp,rec);write(jp,job);write(dp,spec)
        if use_model is not None:
            mp=use_model
        else:
            run_cli(cli,"build",rp,mp)
        output.write_bytes(sentinel)
        proc=run_cli(cli,"differentiate-job",mp,jp,dp,output,expect_success=False)
        cases.append({"name":name,"refused":proc.returncode!=0,
                      "preserved":output.read_bytes()==sentinel,
                      "diagnostic":(proc.stdout+proc.stderr).strip()[:1000]})
    return cases


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli",type=Path,required=True)
    parser.add_argument("--work",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    result=run(args.cli.resolve(),args.work.resolve())
    write(args.output,result)
    print("S5 interval endpoint oracle",result["passed"])
    raise SystemExit(int(not result["passed"]))
