"""Check sector angle jets using rebuilt mechanical jobs and exact surface loads."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import numpy as np


def write(p,value):
    p.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")


def call(cli,*args):
    p=subprocess.run([cli,*args],capture_output=True,text=True)
    if p.returncode:raise ValueError((args,p.stdout,p.stderr))


def recipe(t):
    return dict(family="rational_sector_prism_v1",half_angle_tangent=[t.numerator,t.denominator],
                radii_mm=[3,8],height_mm=7)


def job(law,order=16):
    return dict(schema="solidcraft_mechanical_job_v1",order=order,formulation="displacement",
                material=dict(young_mpa=1700,poisson=.27),
                supports=[dict(select=dict(feature="angle_start"),components=[0,1,2])],
                loads=[law],plan=dict(schema="adaptive_refinement_plan_v1",max_depth=0,max_cells=16,steps=[]))


def solve(cli,folder,t,j):
    rp,mp,jp,out=(folder/k for k in ["recipe.json","model.json","job.json","solution.json"])
    write(rp,recipe(t));write(jp,j);call(cli,"build",rp,mp);call(cli,"solve-job",mp,jp,out)
    return json.loads(out.read_text())


def exact_rate(t,law):
    theta=2*math.atan(t);speed=2/(1+t*t);er=np.array([math.cos(theta),math.sin(theta),0.])
    if "pressure_mpa" in law:return law["pressure_mpa"]*3*7*speed*er
    if "vector_mpa" in law:return np.array(law["vector_mpa"])*(8**2-3**2)/(1+t*t)
    return -np.array(law["stress_mpa"])@er*speed*7*(8-3)


def run(cli,folder):
    folder.mkdir(exist_ok=True,parents=True)
    out=dict(schema="rational_sector_job_derivative_oracle_v1",executable_sha256=hashlib.sha256(cli.read_bytes()).hexdigest(),cases=[])
    laws=[dict(select=dict(feature="inner"),pressure_mpa=.3),
          dict(select=dict(feature="end"),vector_mpa=[.13,-.2,.09]),
          dict(select=dict(feature="angle_end"),stress_mpa=[[.3,.15,0],[.15,-.2,.3],[0,.3,.4]])]
    direction={"schema":"solidcraft_sector_angle_velocity_v1"}
    for t,law in zip([Fraction(1,2),Fraction(3,2),Fraction(3)],laws):
        work=folder/f"t{t.numerator}_{t.denominator}";work.mkdir(exist_ok=True)
        j=job(law);base=solve(cli,work,t,j);write(work/"dir.json",direction)
        call(cli,"differentiate-job",work/"model.json",work/"job.json",work/"dir.json",work/"jet.json")
        d=json.loads((work/"jet.json").read_text())
        exact=exact_rate(float(t),law);actual=np.array(d["applied_resultant_rate_n"])
        volume_rate=7*(8**2-3**2)/(1+float(t)**2)
        checks=dict(volume_relative=abs(d["volume_rate_mm3"]/volume_rate-1),
                    resultant_relative=float(np.linalg.norm(actual-exact)/np.linalg.norm(exact)),
                    adjoint_threshold_fraction=abs(d["compliance_rate_tangent_n_mm"]-d["compliance_rate_adjoint_n_mm"])/
                        max(1e-9,1e-8*abs(d["compliance_rate_tangent_n_mm"])))
        targets=[np.array(d["nodal_force_rate_n_xyz"]),np.array(d["physical_displacement_rate_mm_xyz"]),
                 np.array([d["compliance_rate_tangent_n_mm"]])]
        history=[];previous=None;r4_prev=None
        def differences(approx):
            return [float(np.linalg.norm(a-b)/max(1e-8,1e-5*np.linalg.norm(b))) for a,b in zip(approx,targets)]
        for den in (16,32,64):
            step=Fraction(1,den);values=[]
            for sign in (1,-1):
                result=solve(cli,work,t+sign*step,j)
                assert [(b["id"],b["dofs"]) for b in result["analysis_blocks"]]==[(b["id"],b["dofs"]) for b in base["analysis_blocks"]]
                values.append([np.array(result["force_n_xyz"]),np.array(result["displacement_mm_xyz"]),
                               np.array([result["external_work_n_mm"]])])
            rate=[(a-b)*den/2 for a,b in zip(*values)]
            row=dict(step=1/den,centered_threshold_fractions=differences(rate))
            row["forward_work_remainder"]=float(abs(values[0][2][0]-base["external_work_n_mm"]-targets[2][0]/den))
            if previous is not None:
                r4=[(4*a-b)/3 for a,b in zip(rate,previous)]
                row["fourth_threshold_fractions"]=differences(r4)
                if r4_prev is not None:row["sixth_threshold_fractions"]=differences([(16*a-b)/15 for a,b in zip(r4,r4_prev)])
                r4_prev=r4
            previous=rate;history.append(row)
        passed=max(checks["volume_relative"],checks["resultant_relative"])<1e-8 and checks["adjoint_threshold_fraction"]<=1 and max(history[-1]["sixth_threshold_fractions"])<=1
        out["cases"].append(dict(recipe=recipe(t),job=j,volume_rate_mm3=d["volume_rate_mm3"],resultant_rate_n=actual.tolist(),
                          work_rate_n_mm=targets[2][0],residuals=d["relative_residuals"],checks=checks,history=history,passed=bool(passed)))
        print(t,passed,checks,flush=True)
    event=[]
    for t in (Fraction(70,29),Fraction(169,70)):
        v=solve(cli,folder,t,job(laws[1]));event.append(dict(t=str(t),charts=len(v["analysis_blocks"]),
             volume_mm3=v["volume_mm3"],analytic_volume_mm3=7*(64-9)*math.atan(t),work_n_mm=v["external_work_n_mm"]))
    out["across_chart_change"]=event
    refusals=[]
    for name,src,jj,spec in [("family",dict(family="reference_rectangle_v1"),job(laws[1]),direction),
            ("refinement",recipe(Fraction(1)),job(laws[1]),direction),
            ("spatial_selector",recipe(Fraction(1)),job(laws[1]),direction)]:
        if name=="refinement":jj["plan"].update(max_depth=1,steps=[dict(refine=[dict(root="whole",level=0,index=[0,0,0])])])
        if name=="spatial_selector":jj["supports"][0]["select"]=dict(plane=dict(axis=2,coordinate_mm=0))
        write(folder/"recipe.json",src);write(folder/"job.json",jj);write(folder/"direction.json",spec)
        call(cli,"build",folder/"recipe.json",folder/"model.json")
        output=folder/"refused.json";output.write_bytes(b"kept")
        proc=subprocess.run([cli,"differentiate-job",folder/"model.json",folder/"job.json",folder/"direction.json",output],capture_output=True,text=True)
        refusals.append(dict(name=name,diagnostic=proc.stdout.strip(),passed=proc.returncode!=0 and output.read_bytes()==b"kept"))
    out["refusals"]=refusals;out["qualified"]=all(i["passed"] for i in [*out["cases"],*refusals]) and [v["charts"] for v in event]==[1,2]
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("cli","work","output"):p.add_argument(name,type=Path)
    a=p.parse_args();r=run(a.cli.resolve(),a.work);write(a.output,r)
    raise SystemExit(int(not r["qualified"]))
