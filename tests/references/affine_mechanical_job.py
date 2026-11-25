#!/usr/bin/env python3
"""Independent SciPy oracle for affine mechanical-job supports (schema v2).

The oracle rebuilds K and boundary loads from serialized Q2 geometry, evaluates
the exact rational affine traces at every selected support control, forms
PᵀKP and Pᵀ(f−Kg) from the exported exact lifts, and solves the reduced sparse
system independently. It also checks work split, reactions, and balance.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.sparse import block_diag, coo_matrix
from scipy.sparse.linalg import spsolve

sys.path.insert(0, str(Path(__file__).parent))
import mechanical_job as m5


def exact_motion(motion, h, component):
    den = int(motion["denominator"])
    values = [int(v) for v in motion["numerators"][component]]
    homogeneous_numerator = values[0] * int(h[3]) + sum(values[i + 1] * int(h[i]) for i in range(3))
    return Fraction(homogeneous_numerator, den * int(h[3]))


def lift_matrix(result, ndof):
    maps=[]; offsets=[]
    for component, lift in enumerate(result["component_lifts"]):
        nfree=len(lift["free_roots"])
        rr=[];cc=[];vv=[];g=np.zeros(ndof)
        for line in lift["lines"]:
            row=int(line["dof"])
            for root,num,den in line["terms"]:
                rr.append(row);cc.append(int(root));vv.append(int(num)/int(den))
            g[row]=int(line["offset"][0])/int(line["offset"][1])
        maps.append(coo_matrix((vv,(rr,cc)),shape=(ndof,nfree)).tocsr())
        offsets.append(g)
    component_major=block_diag(maps,format="csr")
    physical_rows=np.asarray([c*ndof+d for d in range(ndof) for c in range(3)],dtype=np.int64)
    return component_major[physical_rows].tocsr(),np.concatenate(offsets)


def interleave_vector(component_vector, ndof):
    return np.asarray(component_vector).reshape(3,ndof).T.reshape(3*ndof)


def physical_points(result, ndof):
    points=np.full((ndof,3),np.nan)
    for block in result["analysis_blocks"]:
        for h,node in zip(block["controls"],block["dofs"]):
            node=int(node);x=np.asarray(h[:3],dtype=float)/int(h[3])
            if np.all(np.isfinite(points[node])) and np.linalg.norm(points[node]-x)>1e-10:
                raise ValueError(f"shared scalar DOF {node} has inconsistent position")
            points[node]=x
    if not np.all(np.isfinite(points)): raise ValueError("analysis space has missing coordinates")
    return points


def trace_error(job, result, ndof):
    blocks={b["id"]:b for b in result["analysis_blocks"]}
    field=np.asarray(result["displacement_mm_xyz"],dtype=float).reshape(ndof,3)
    largest=0.0
    for support,selected in zip(job["supports"],result["supports"]):
        motion=support.get("affine_mm")
        for block_id,side in selected["faces"]:
            block=blocks[block_id]
            axis,upper=m5.parse_side({"side":side})
            for local,h in enumerate(block["controls"]):
                ijk=(local//9,(local//3)%3,local%3)
                if ijk[axis]!=2*upper: continue
                node=int(block["dofs"][local])
                for component in support["components"]:
                    exact=Fraction(0) if motion is None else exact_motion(motion,h,component)
                    error=abs(field[node,component]-float(exact))
                    largest=max(largest,error)
    return largest


def run_valid(binary, fixture, work):
    data=json.loads(fixture.read_text());work.mkdir(parents=True,exist_ok=True)
    stem=data["id"]
    recipe=work/f"{stem}.recipe.json";model_path=work/f"{stem}.model.json"
    job_path=work/f"{stem}.job.json";result_path=work/f"{stem}.result.json"
    recipe.write_text(json.dumps(data["recipe"],separators=(",",":")))
    job_path.write_text(json.dumps(data["job"],separators=(",",":")))
    subprocess.run([str(binary),"build",str(recipe),str(model_path)],check=True,capture_output=True,text=True,timeout=1800)
    solved=subprocess.run([str(binary),"solve-job",str(model_path),str(job_path),str(result_path)],capture_output=True,text=True,timeout=1800)
    if solved.returncode: raise RuntimeError(f"job failed for {stem}: {solved.stdout} {solved.stderr}")
    model=json.loads(model_path.read_text());result=json.loads(result_path.read_text());job=data["job"]
    k,ndof=m5.assemble(result["analysis_blocks"],job["material"]["young_mpa"],job["material"]["poisson"],job["order"])
    f,applied,moment,area,load_reports,load_scale,max_radius=m5.independent_load(model,result,job,job["order"],ndof)
    u=np.asarray(result["displacement_mm_xyz"],dtype=float)
    P,g_block=lift_matrix(result,ndof)
    g=interleave_vector(g_block,ndof)
    K=P.T@k@P
    rhs=np.asarray(P.T@(f-k@g)).ravel()
    if K.shape[0]:
        reduced=spsolve(K.tocsc(),rhs)
        u_ind=np.asarray(P@reduced).ravel()+g
        reduced_res=float(np.linalg.norm(K@reduced-rhs)/max(np.linalg.norm(rhs),1.0))
    else:
        u_ind=g
        reduced_res=0.0
    coeff_error=float(np.linalg.norm(u-u_ind)/max(np.linalg.norm(u_ind),1e-30))
    coeff_max=float(np.max(np.abs(u-u_ind)))
    trace_max=trace_error(job,result,ndof)
    residual=k@u-f
    points=physical_points(result,ndof)
    reaction=residual.reshape(ndof,3).sum(axis=0)
    reaction_moment=np.cross(points,residual.reshape(ndof,3)).sum(axis=0)
    traction_work=float(f@u)
    support_work=float(u@residual)
    external_work=traction_work+support_work
    energy=float(0.5*u@(k@u))
    work_scale=max(abs(traction_work),abs(support_work),abs(external_work),abs(2*energy),1.0)
    work_error=abs(external_work-2*energy)/work_scale
    reaction_error=float(np.linalg.norm(reaction-np.asarray(result["support_reaction_n"])) / max(load_scale,1.0))
    reaction_moment_error=float(np.linalg.norm(reaction_moment-np.asarray(result["support_moment_n_mm"])) / max(load_scale*max(max_radius,1.0),1.0))
    force_balance=float(np.linalg.norm(reaction+applied)/max(load_scale,1.0))
    moment_balance=float(np.linalg.norm(reaction_moment+moment)/max(load_scale*max(max_radius,1.0),1.0))
    reported={"traction_work_n_mm":float(result["traction_work_n_mm"]),"support_work_n_mm":float(result["support_work_n_mm"]),"external_work_n_mm":float(result["external_work_n_mm"])}
    work_field_error=max(abs(reported["traction_work_n_mm"]-traction_work),abs(reported["support_work_n_mm"]-support_work),abs(reported["external_work_n_mm"]-external_work))/work_scale
    load_force_error=float(np.max(np.abs(f-np.asarray(result["force_n_xyz"]))))
    affine_error=None
    if data["oracle"]["kind"] in ("affine","rigid"):
        grad=np.asarray(data["oracle"]["gradient"],dtype=float);offset=np.asarray(data["oracle"]["offset"],dtype=float)
        expected=points@grad.T+offset
        affine_error=float(np.max(np.abs(u.reshape(ndof,3)-expected)))
    record={"id":stem,"status":"qualified","scalar_dofs":ndof,"vector_dofs":len(u),
        "sparse_nonzeros":int(k.nnz),"reduced_dofs":int(K.shape[0]),
        "independent_reduced_residual":reduced_res,"scipy_displacement_relative_error":coeff_error,
        "scipy_displacement_max_abs_error_mm":coeff_max,"affine_support_trace_max_error_mm":trace_max,
        "exact_affine_field_max_error_mm":affine_error,"reintegrated_resultant_n":applied.tolist(),
        "reintegrated_moment_n_mm":moment.tolist(),"load_rhs_max_abs_error_n":load_force_error,
        "reaction_n_independent":reaction.tolist(),"reaction_moment_n_mm_independent":reaction_moment.tolist(),
        "reaction_relative_error":reaction_error,"reaction_moment_relative_error":reaction_moment_error,
        "force_balance_relative":force_balance,"moment_balance_relative":moment_balance,
        "traction_work_n_mm":traction_work,"support_work_n_mm":support_work,"external_work_n_mm":external_work,
        "strain_energy_n_mm":energy,"work_energy_relative_error":work_error,
        "reported_work_relative_error":work_field_error,"rust_residual":result["relative_residual"],
        "rigid_constraint_rank":result["rigid_constraint_rank"],"quadrature_order":job["order"],
        "load_groups":load_reports}
    if reduced_res>1e-10 or coeff_error>1e-8 or trace_max>1e-10 or work_error>1e-8 or work_field_error>1e-8 or reaction_error>1e-8 or reaction_moment_error>1e-8 or force_balance>1e-8 or moment_balance>1e-8 or result["relative_residual"]>1e-10 or (affine_error is not None and affine_error>1e-10):
        raise AssertionError(f"M6 independent oracle failed {stem}: {record}")
    if data["oracle"]["kind"]=="rigid" and energy>1e-10:
        raise AssertionError(f"rigid motion stores nonzero strain energy: {energy}")
    return record


def run_refusal(binary,fixture,work):
    data=json.loads(fixture.read_text())
    with tempfile.TemporaryDirectory(dir=work) as temp:
        temp=Path(temp);recipe=temp/"recipe.json";model=temp/"model.json";job=temp/"job.json";out=temp/"result.json"
        recipe.write_text(json.dumps(data["recipe"]));job.write_text(json.dumps(data["job"]))
        subprocess.run([str(binary),"build",str(recipe),str(model)],check=True,capture_output=True,text=True,timeout=1800)
        out.write_text("previous-valid-result\n")
        response=subprocess.run([str(binary),"solve-job",str(model),str(job),str(out)],capture_output=True,text=True,timeout=1800)
        payload=json.loads(response.stdout)
        if response.returncode==0 or payload.get("reason")!=data["expected_error"] or out.read_text()!="previous-valid-result\n":
            raise AssertionError(f"M6 refusal mismatch {data['id']}: {payload}")
        return {"id":data["id"],"status":"refused","reason":payload["reason"],"previous_output_preserved":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary",type=Path)
    parser.add_argument("work",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args();args.work.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();records=[]
    for fixture in sorted(Path(__file__).parent.joinpath("affine_mechanical_jobs").glob("*.json")):
        data=json.loads(fixture.read_text())
        records.append(run_refusal(args.binary.resolve(),fixture,args.work) if "expected_error" in data else run_valid(args.binary.resolve(),fixture,args.work))
    if sum(x["status"]=="qualified" for x in records)!=8 or sum(x["status"]=="refused" for x in records)!=3:
        raise AssertionError("M6 corpus must contain eight positive cases and three refusals")
    output={"schema":"solidcraft_affine_mechanical_job_oracle_v1","oracle":"independent rational affine trace checks plus SciPy CSR congruence solve",
        "cases":records,"elapsed_seconds":time.monotonic()-started,
        "thresholds":{"reduced_residual":1e-10,"coefficient_relative_error":1e-8,"trace_and_affine_absolute_mm":1e-10,"work_and_balance_relative":1e-8}}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    print(json.dumps(output,indent=2,allow_nan=False))


if __name__=="__main__": main()
