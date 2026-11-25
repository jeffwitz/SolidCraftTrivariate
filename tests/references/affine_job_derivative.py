"""Verify material shape derivatives against fresh rationally perturbed native jobs."""
import argparse
import copy
from fractions import Fraction
from functools import reduce
import json
import math
from pathlib import Path
import subprocess

import numpy as np

HERE=Path(__file__).resolve().parent


def read(p): return json.loads(p.read_text())
def write(p,v): p.write_text(json.dumps(v,allow_nan=False,indent=2)+"\n")


def cases():
    polygon=next(r["recipe"] for r in read(HERE/"polygon_prism_corpus.json") if r["id"]=="L")
    recipe=dict(polygon,top_shift_mm=[0,0])
    job=dict(schema="solidcraft_mechanical_job_v1",material=dict(young_mpa=10000,poisson=.25),
        formulation="displacement",order=4,plan=dict(schema="adaptive_refinement_plan_v1",max_depth=0,max_cells=1000,steps=[]),
        supports=[dict(select=dict(feature="start"),components=[0,1,2])],
        loads=[dict(select=dict(feature="end"),vector_mpa=[1.,-2.,.5])])
    ident=np.eye(3,dtype=int);zero=np.zeros((3,3),dtype=int)
    yield "material_translation",recipe,job,1,zero,[1,2,-1],"translation"
    yield "dilatation",recipe,job,Fraction(1,2),2*ident,[0,0,0],"dilatation"
    # Two components of traction and an interior reentrant corner.
    local=copy.deepcopy(job)
    root="tri/p00/p01/p02/p00"
    local["plan"].update(max_depth=1,steps=[dict(refine=[dict(root=root,level=0,index=[0,0,0])])])
    yield "hanging_shear",recipe,local,1,[[0,1,0],[0,0,0],[1,0,0]],[0,0,0],"nonlinear"
    bend=read(HERE/"mechanical_jobs/open-bend-affine.json")
    yield "curved_stress",bend["recipe"],bend["job"],1,[[0,0,0],[1,0,1],[0,0,0]],[0,0,0],"nonlinear"
    cavity=read(HERE/"mechanical_jobs/two-cavities.json")
    yield "cavity_pressure",cavity["recipe"],cavity["job"],1,[[0,1,0],[0,0,1],[0,0,0]],[0,0,0],"nonlinear"


def image(source,base,B,shift,s):
    M=(np.eye(3,dtype=int)+s*np.asarray(B,dtype=object))*base
    t=s*np.asarray(shift,dtype=object)
    values=[Fraction(x) for x in [*M.flat,*t]]
    den=reduce(math.lcm,[v.denominator for v in values],1)
    ints=[int(v*den) for v in values]
    return dict(family="affine_image_v1",source=source,transform=dict(
                matrix=np.array(ints[:9]).reshape(3,3).tolist(),translation=ints[9:],
                denominator=den,homogeneous_range="i64"))


def run(cli,work):
    work.mkdir(exist_ok=True,parents=True);cases_out=[]
    for name,source,job,base,B,shift,law in cases():
        folder=work/name;folder.mkdir(exist_ok=True);jp=folder/"job.json";write(jp,job)
        B=np.asarray(B,dtype=int);dp=folder/"direction.json"
        direction=dict(schema="solidcraft_affine_velocity_v1",linear=B.tolist(),offset_mm=shift);write(dp,direction)
        def call(*args):
            cmd=subprocess.run([str(cli),*map(str,args)],capture_output=True,text=True)
            if cmd.returncode:raise ValueError(f"{name} {args[0]}: {cmd.stdout}")
            return json.loads(cmd.stdout)
        def solve(s):
            tag=str(s).replace("/","_");rp=folder/(tag+".recipe");mp=folder/(tag+".model");result=folder/(tag+".result")
            write(rp,image(source,base,B,shift,s));call("build",rp,mp);call("solve-job",mp,jp,result)
            return read(result),mp
        original,mp=solve(Fraction(0));out=folder/"derivative.json"
        call("differentiate-job",mp,jp,dp,out);der=read(out)
        u=np.array(original["displacement_mm_xyz"]);du=np.array(der["physical_displacement_rate_mm_xyz"])
        f=np.array(original["force_n_xyz"]);df=np.array(der["nodal_force_rate_n_xyz"])
        C=original["external_work_n_mm"];Cd=der["compliance_rate_tangent_n_mm"];vol=original["volume_mm3"]
        adj_error=abs(Cd-der["compliance_rate_adjoint_n_mm"])/max(abs(C),abs(Cd),1e-20)
        verror=float(abs(der["volume_rate_mm3"]/vol-np.trace(B)))
        history=[]
        def relative(delta,scale):return float(np.linalg.norm(delta)/max(np.linalg.norm(scale),1e-30))
        for den in (16,32,64):
            s=Fraction(1,den);p,_=solve(s);m,_=solve(-s);h=float(s)
            up=np.array(p["displacement_mm_xyz"]);um=np.array(m["displacement_mm_xyz"])
            fp=np.array(p["force_n_xyz"]);fm=np.array(m["force_n_xyz"])
            cp=p["external_work_n_mm"];cm=m["external_work_n_mm"]
            history.append(dict(step=h,
                centered_u_error=relative((up-um)/(2*h)-du,u),
                centered_load_error=relative((fp-fm)/(2*h)-df,f),
                centered_compliance_error=abs((cp-cm)/(2*h)-Cd)/abs(C),
                remainder_u=relative(up-u-h*du,u),
                remainder_compliance=abs(cp-C-h*Cd)/abs(C)))
        ratios={key:[history[i][key]/history[i+1][key] if history[i+1][key]>0 else None for i in (0,1)]
                for key in ("remainder_u","remainder_compliance")}
        # Quadratic remainder is tested when it is resolved above the primal tolerance.
        taylor=all(3.1<=ratio<=5.1 for key in ratios for i,ratio in enumerate(ratios[key])
                   if history[i+1][key]>1e-8)
        final=history[-1]
        checks=dict(adjoint=adj_error<1e-8,volume=verror<1e-9,taylor=taylor,
                    centered_load=final["centered_load_error"]<.003,
                    centered_u=final["centered_u_error"]<.01,
                    centered_compliance=final["centered_compliance_error"]<.01,
                    residuals=max(der["relative_residuals"])<1e-9)
        if law=="translation":checks["exact_invariance"]=max(abs(Cd/C),relative(du,u),relative(df,f))<1e-10
        if law=="dilatation":checks["exact_scale_law"]=max(abs(Cd/C-6),relative(du-2*u,u),relative(df-4*f,f))<1e-8
        row=dict(id=name,source=source,base_scale=float(base),direction=direction,job=job,checks=checks,
                 qualified=all(checks.values()),taylor_ratios=ratios,history=history,
                 compliance_n_mm=C,compliance_rate_n_mm=Cd,adjoint_error=adj_error,volume_rate_error=float(verror),
                 residuals=der["relative_residuals"],cells=len(original["analysis_blocks"]))
        cases_out.append(row)
        print(name,checks,ratios,flush=True)
    # Verify atomic refusals for unsupported hypotheses using an already accepted target.
    mp=work/"material_translation/0.model";jp=work/"refusal.job";dp=work/"material_translation/direction.json"
    target=work/"refusal.result";target.write_bytes(b"previous accepted result\n");cases_refused=[]
    base=read(work/"material_translation/job.json")
    for name in ("plane_selection","projected_bulk","prescribed_motion"):
        job=copy.deepcopy(base)
        if name=="plane_selection":job["supports"][0]["select"]["plane"]=dict(axis=0,coordinate_mm=0)
        elif name=="projected_bulk":job["formulation"]="projected_bulk_q1"
        else:
            job["schema"]="solidcraft_mechanical_job_v2"
            job["supports"][0]["affine_mm"]=dict(denominator=1,numerators=[[0]*4 for _ in range(3)])
        write(jp,job);p=subprocess.run([str(cli),"differentiate-job",mp,jp,dp,target],capture_output=True,text=True)
        cases_refused.append(dict(id=name,diagnostic=json.loads(p.stdout),qualified=bool(p.returncode and target.read_bytes()==b"previous accepted result\n")))
    return dict(schema="affine_job_derivative_cold_rebuild_oracle_v1",cases=cases_out,refusals=cases_refused)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","work","output"):p.add_argument(k,type=Path)
    args=p.parse_args();report=run(args.cli.resolve(),args.work);write(args.output,report)
    if not all(r["qualified"] for r in report["cases"]+report["refusals"]):raise SystemExit(1)
