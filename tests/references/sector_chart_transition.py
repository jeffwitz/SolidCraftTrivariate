"""Measure a 135-degree representation change using geometric and fixed analysis charts."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import subprocess
import numpy as np

from sector_job_sensitivity import call,job,recipe,solve,write


def run(cli,work):
    work.mkdir(exist_ok=True,parents=True)
    left=[F(12,5),F(41,17),F(70,29),F(239,99)]
    right=[F(17,7),F(29,12),F(99,41),F(169,70)]
    t0=1+math.sqrt(2)
    cases=[]
    for name,law,refined in [
        ("cap_vector",dict(select=dict(feature="end"),vector_mpa=[.13,-.2,.09]),False),
        ("inner_pressure",dict(select=dict(feature="inner"),pressure_mpa=.3),False),
        ("local_cap_vector",dict(select=dict(feature="end"),vector_mpa=[.13,-.2,.09]),True)]:
        folder=work/name;folder.mkdir(exist_ok=True)
        base=job(law,12)
        if refined:
            base["plan"].update(max_depth=2,max_cells=64,steps=[
                dict(refine=[dict(root="first",level=0,index=[0,0,0])]),
                dict(refine=[dict(root="first",level=1,index=[0,0,0])])])
        wrapped=dict(schema="solidcraft_sector_analysis_job_v1",analysis_partition="angular_halves",job=base)
        spec=dict(schema="solidcraft_sector_angle_velocity_v2")
        write(folder/"direction.json",spec)
        intervals=[]
        for a,b in zip(left,right):
            assert float(a)<t0<float(b)
            value=[]
            for t in (a,b):
                sol=solve(cli,folder,t,wrapped)
                call(cli,"differentiate-job",folder/"model.json",folder/"job.json",folder/"direction.json",folder/"rate.json")
                rate=json.loads((folder/"rate.json").read_text())
                assert rate["analysis_partition"]==sol["analysis_partition"]
                obs=dict(t=str(t),assembly=sol["analysis_partition"],cells=len(sol["analysis_blocks"]),
                         dofs=sol["scalar_dofs"],work=sol["external_work_n_mm"],work_rate=rate["compliance_rate_tangent_n_mm"],
                         adjoint_work_rate=rate["compliance_rate_adjoint_n_mm"],volume=sol["volume_mm3"],
                         reference_volume=7*55*math.atan(t))
                u=np.array(sol["displacement_mm_xyz"]);du=np.array(rate["physical_displacement_rate_mm_xyz"])
                if t==b or not refined:
                    default=solve(cli,folder,t,base)
                    obs["source_partition_work"]=default["external_work_n_mm"]
                    if t==b:
                        assert default["analysis_blocks"]==sol["analysis_blocks"]
                        assert np.array_equal(default["displacement_mm_xyz"],u)
                value.append((obs,u,du))
            (lo,u0,d0),(hi,u1,d1)=value
            h=float(b-a);secant=(hi["work"]-lo["work"])/h;field=(u1-u0)/h
            avg=.5*(lo["work_rate"]+hi["work_rate"]);avgfield=(d0+d1)/2
            lwork=lo["work"]+(t0-float(a))*lo["work_rate"]
            rwork=hi["work"]+(t0-float(b))*hi["work_rate"]
            intervals.append(dict(left=lo,right=hi,gap_t=h,
                        extrapolated_work_gap_relative=abs(lwork-rwork)/max(abs(lwork),1e-12),
                        secant_work_rate_relative=abs(secant-avg)/max(abs(avg),1e-12),
                        secant_field_rate_relative=float(np.linalg.norm(field-avgfield)/max(np.linalg.norm(avgfield),1e-12))))
        last=intervals[-1]
        checks=dict(extrapolated_work=last["extrapolated_work_gap_relative"]<=1e-8,
                    work_rate=last["secant_work_rate_relative"]<=1e-4,
                    field_rate=last["secant_field_rate_relative"]<=1e-4,
                    charts=all(r["left"]["assembly"]["source_charts"]==1 and r["right"]["assembly"]["source_charts"]==2
                               and r["left"]["dofs"]==r["right"]["dofs"] for r in intervals),
                    adjoint=all(abs(v["work_rate"]-v["adjoint_work_rate"])<=max(1e-9,1e-8*abs(v["work_rate"]))
                                for row in intervals for v in [row["left"],row["right"]]),
                    analytic_volume=all(abs(v["volume"]/v["reference_volume"]-1)<=1e-8
                                for row in intervals for v in [row["left"],row["right"]]))
        if not refined:checks["scalar_dimension"]=last["left"]["dofs"]==45
        cases.append(dict(name=name,job=wrapped,intervals=intervals,checks=checks,passed=all(checks.values())))
        print(name,checks,flush=True)
    refusals=[]
    for label,rec,policy in [("wrong_family",dict(family="reference_rectangle_v1"),"angular_halves"),
                             ("unknown_policy",recipe(F(1)),"approximate")]:
        write(work/"recipe.json",rec);call(cli,"build",work/"recipe.json",work/"model.json")
        write(work/"job.json",dict(schema="solidcraft_sector_analysis_job_v1",analysis_partition=policy,job=job(law)))
        target=work/"protected.json";target.write_bytes(b"previous")
        proc=subprocess.run([cli,"solve-job",work/"model.json",work/"job.json",target],capture_output=True,text=True)
        refusals.append(dict(id=label,diagnostic=proc.stdout.strip(),passed=proc.returncode!=0 and target.read_bytes()==b"previous"))
    return dict(schema="sector_chart_transition_oracle_v1",executable_sha256=hashlib.sha256(cli.read_bytes()).hexdigest(),
                cases=cases,refusals=refusals,qualified=all(v["passed"] for v in [*cases,*refusals]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","work","output"):p.add_argument(k,type=Path)
    a=p.parse_args();out=run(a.cli.resolve(),a.work);write(a.output,out)
    raise SystemExit(int(not out["qualified"]))
