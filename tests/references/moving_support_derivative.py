"""Independently perturb complete jobs with spatial affine Dirichlet data."""
import argparse
import copy
from fractions import Fraction
import json
from pathlib import Path
import subprocess

import numpy as np

from affine_job_derivative import image, read, write
from body_load import solve_command

HERE=Path(__file__).resolve().parent


def cases():
    folder=HERE/"affine_mechanical_jobs"
    bar=read(folder/"cube-normal-traction.json")
    zero=[[0,0,0]]*3
    yield "translated_bar",bar["recipe"],bar["job"],zero,[1,0,0],"translation"
    job=copy.deepcopy(bar["job"])
    job["loads"]=[]
    job["supports"]=[dict(select=dict(feature="*"),components=[0,1,2],
                          affine_mm=dict(denominator=1000,numerators=[[3,0,2,0],[0,0,0,0],[0,0,0,0]]))]
    yield "prescribed_shear",bar["recipe"],job,[[0,0,0],[0,1,0],[0,0,0]],[0,1,0],"volume"
    curved=read(folder/"curved-hanging.json")
    job=curved["job"];job["schema"]="solidcraft_mechanical_job_v4"
    job["supports"][0]["affine_mm"]=dict(denominator=1000,numerators=[[1,0,1,0],[0,0,0,1],[0,0,0,0]])
    job["body_load"]=dict(constant_n_per_mm3=[.005,-.002,.001],gradient_n_per_mm4=[[0,.0001,0],[.0001,0,.0001],[0,0,0]])
    yield "hanging_bend_body",curved["recipe"],job,[[0,1,0],[0,0,0],[0,0,0]],[0,0,1],"numeric"


def direction(B,b):return dict(schema="solidcraft_affine_velocity_v2",
                                  affine_mm=dict(denominator=1,numerators=[[o,*row] for o,row in zip(b,B)]))


def run(cli,work):
    rows=[]
    for name,source,job,B,b,law in cases():
        folder=work/name;folder.mkdir(exist_ok=True,parents=True)
        jp,mp,rp,out=[folder/(n+".json") for n in ("job","model","recipe","result")]
        write(jp,job)
        def solve(s):
            write(rp,image(source,1,B,b,s))
            solve_command(cli,"build",rp,mp);solve_command(cli,"solve-job",mp,jp,out)
            return read(out)
        base=solve(Fraction(0))
        dp=folder/"direction.json";write(dp,direction(B,b))
        solve_command(cli,"differentiate-job",mp,jp,dp,out)
        d=read(out)
        def values(result):
            traction=result["traction_work_n_mm"]
            body=result.get("body_load",{}).get("work_n_mm",0.)
            return (np.array(result["displacement_mm_xyz"]),np.array([traction+body,result["support_work_n_mm"],result["strain_energy_n_mm"]]))
        u,work0=values(base)
        du=np.array(d["physical_displacement_rate_mm_xyz"])
        rates=np.array([d["applied_work_rate_tangent_n_mm"],d["support_work_rate_n_mm"],d["strain_energy_rate_tangent_n_mm"]])
        base_agreement=abs(work0-np.array([d["applied_work_n_mm"],d["support_work_n_mm"],d["strain_energy_n_mm"]]))
        hist=[];previous=None;previous4=None
        def errors(a,ref):
            delta=abs(a-ref)
            return np.minimum(delta/1e-9,delta/(1e-5*np.maximum(abs(ref),1e-30)))
        for denom in (8,16,32):
            step=Fraction(1,denom)
            p,w=values(solve(step));m,z=values(solve(-step))
            grad=np.r_[(p-m)/(2*float(step)),(w-z)/(2*float(step))]
            expected=np.r_[du,rates]
            def discrepancy(g):
                delta=g-expected
                return dict(displacement_norm_mm=float(np.linalg.norm(delta[:-3])),
                            displacement_relative=float(np.linalg.norm(delta[:-3])/max(np.linalg.norm(du),1e-30)),
                            work_errors_n_mm=delta[-3:].tolist(),
                            tolerance_fraction=max(min(np.linalg.norm(delta[:-3])/1e-9,
                                np.linalg.norm(delta[:-3])/(1e-5*max(np.linalg.norm(du),1e-30))),max(errors(g[-3:],rates))))
            rec=dict(step=float(step),centered=discrepancy(grad))
            if previous is not None:
                fourth=(4*grad-previous)/3
                rec["fourth"]=discrepancy(fourth)
                if previous4 is not None:rec["sixth"]=discrepancy((16*fourth-previous4)/15)
                previous4=fourth
            previous=grad
            hist.append(rec)
        adjoints=[d["applied_work_rate_adjoint_n_mm"],d["strain_energy_rate_adjoint_n_mm"]]
        errs=abs(np.array(adjoints)-rates[[0,2]])
        adj_pass=np.all(errs<=np.maximum(1e-9,1e-8*abs(rates[[0,2]])))
        exact={}
        if law=="translation":
            fixed=np.tile([.001,0.,0.],len(du)//3)
            exact["max_displacement_mm"]=float(abs(fixed-du).max())
            exact["strain_energy_rate"]=abs(rates[2])
        elif law=="volume":
            velocity=np.array(d["geometry_velocity_mm_xyz"],dtype=float)
            fixed=np.zeros_like(velocity);fixed[:,0]=.002*velocity[:,1]
            exact["max_displacement_mm"]=float(abs(fixed.ravel()-du).max())
            exact["strain_energy_rate"]=abs(rates[2]-work0[2])
        rows.append(dict(name=name,scalar_dofs=base["scalar_dofs"],baseline_work=work0.tolist(),rates=rates.tolist(),
                         applied_and_strain_adjoint_errors=errs.tolist(),base_work_errors=base_agreement.tolist(),
                         work_rate_closure=d["work_rate_closure_n_mm"],independent_exact=exact,
                         centered_history=hist,relative_residuals=d["relative_residuals"],
                         passed=bool(hist[-1]["sixth"]["tolerance_fraction"]<=1 and adj_pass and
                                     max(abs(d["work_rate_closure_n_mm"]),*base_agreement,*exact.values())<1e-9)))
    folder=work/"refusal";folder.mkdir(exist_ok=True)
    bar=read(HERE/"affine_mechanical_jobs/cube-normal-traction.json");job=bar["job"]
    rp,mp,jp,dp,out=[folder/(n+".json") for n in ("recipe","model","job","direction","out")]
    write(rp,bar["recipe"]);solve_command(cli,"build",rp,mp)
    failed=[]
    for case in ("support_velocity_conflict","v1_nonhomogeneous","zero_velocity_denominator"):
        spec=copy.deepcopy(job);flow=direction([[0,0,0]]*3,[1,0,0])
        if case=="support_velocity_conflict":
            spec["supports"].append(dict(select=dict(feature="x-10"),components=[0]))
        if case=="v1_nonhomogeneous":flow=dict(schema="solidcraft_affine_velocity_v1",linear=[[0,0,0]]*3,offset_mm=[1,0,0])
        if case=="zero_velocity_denominator":flow["affine_mm"]["denominator"]=0
        write(jp,spec);write(dp,flow)
        solve_command(cli,"solve-job",mp,jp,out)
        original=out.read_bytes()
        proc=subprocess.run([cli,"differentiate-job",mp,jp,dp,out],capture_output=True,text=True)
        failed.append(dict(name=case,refused=proc.returncode!=0,preserved=original==out.read_bytes(),diagnostic=proc.stdout.strip()))
    return dict(schema="prescribed_motion_material_derivative_oracle_v1",cases=rows,refusals=failed,
                passed=all(r["passed"] for r in rows) and all(r["refused"] and r["preserved"] for r in failed))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","work","output"):p.add_argument(k,type=Path)
    a=p.parse_args();report=run(a.cli.resolve(),a.work)
    write(a.output,report);print("Moving supports",report["passed"]);raise SystemExit(int(not report["passed"]))
