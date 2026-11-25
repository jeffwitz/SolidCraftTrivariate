"""Independent volume-force and material-rate checks for order-adaptive integration."""
import argparse
import copy
from fractions import Fraction
import json
from pathlib import Path
import subprocess

import numpy as np
from adaptive_surface_load import reference_result_load
from affine_job_derivative import image
from body_load import solve_command
from solver_package import shape

HERE=Path(__file__).resolve().parent


def read(p):return json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")


def oracle(blocks,density,order,B,b):
    x,w=np.polynomial.legendre.leggauss(order)
    x=(x+1)/2;w=w/2
    uvw=np.array(np.meshgrid(x,x,x,indexing="ij")).reshape(3,-1).T
    rule=np.einsum("i,j,k->ijk",w,w,w).ravel()
    ndof=max(max(block["dofs"]) for block in blocks)+1
    f=np.zeros((ndof,3));rate=np.zeros_like(f);V=0.;moment=np.zeros(3);scale=0.;mscale=0.
    grad=np.asarray(density["gradient_n_per_mm4"]);const=np.asarray(density["constant_n_per_mm3"])
    for block in blocks:
        ids=block["dofs"]
        for start in range(0,len(uvw),1024):
            selection=slice(start,start+1024)
            N,point,J,_=shape(uvw[selection],block["controls"])
            dw=rule[selection]*np.linalg.det(J)
            assert np.all(dw>0)
            V+=dw.sum()
            load=(const+point@grad.T)*dw[:,None]
            dload=((point@B.T+b)@grad.T)*dw[:,None]+np.trace(B)*load
            np.add.at(f,ids,N.T@load);np.add.at(rate,ids,N.T@dload)
            m=np.cross(point,load);moment+=m.sum(axis=0)
            scale+=np.linalg.norm(load,axis=1).sum();mscale+=np.linalg.norm(m,axis=1).sum()
    return dict(f=f.ravel(),rate=rate.ravel(),volume=V,moment=moment,force=f.sum(axis=0),scale=scale,mscale=mscale)


def candidates():
    circ={r["id"]:r["recipe"] for r in read(HERE/"circular_polygon_corpus.json")}
    finite={r["id"]:r["recipe"] for r in read(HERE/"interval_bores_corpus.json")}
    mixed=dict(constant_n_per_mm3=[.01,-.02,.005],gradient_n_per_mm4=[[.001,0,0],[0,-.001,0],[0,0,.0005]])
    yield "extreme_surface_body",circ["extreme_shear"],mixed,dict(select=dict(feature="edge/outer/v00/v01"),vector_mpa=[.2,.4,-.1])
    odd=dict(constant_n_per_mm3=[0,0,0],gradient_n_per_mm4=[[.015,0,0],[0,0,0],[0,0,0]])
    yield "zero_resultant_closed",finite["closed"],odd,None
    yield "two_blind_and_sealed",finite["mixed"],mixed,None


def run(cli,work):
    rows=[];refusals=[]
    for name,source,density,load in candidates():
        folder=work/name;folder.mkdir(exist_ok=True,parents=True)
        rp,mp,jp,op,dp=[folder/f"{v}.json" for v in ("recipe","model","job","result","direction")]
        write(rp,source);solve_command(cli,"build",rp,mp)
        job=dict(schema="solidcraft_mechanical_job_v6",material=dict(young_mpa=1000,poisson=.25),formulation="displacement",
                 order=4,plan=dict(schema="adaptive_refinement_plan_v1",max_depth=0,max_cells=2048,steps=[]),
                 supports=[dict(select=dict(feature="start"),components=[0,1,2])],
                 loads=[load] if load else [],body_load=density,
                 body_quadrature=dict(min_order=2,max_order=16,relative_tolerance=1e-9,max_points=8000000))
        if load:job["surface_quadrature"]=dict(split="directional",derivative="freeze_primal",
                                               max_depth=8,max_panels=4096,relative_tolerance=1e-9)
        write(jp,job);solve_command(cli,"solve-job",mp,jp,op);res=read(op)
        B=np.array([[0,1,0],[0,0,0],[0,0,0]],dtype=int);offset=np.array([1,0,0])
        a=oracle(res["analysis_blocks"],density,16,B,offset);b=oracle(res["analysis_blocks"],density,24,B,offset)
        sol=np.array(res["displacement_mm_xyz"]);sum_force=b["f"].copy()
        if load:
            law={k:v for k,v in load.items() if k!="select"}
            sum_force+=reference_result_load(res,law,96)["force_n_xyz"]
        def normalized(v,scale):return float(np.linalg.norm(v)/max(1.,scale))
        actual=res["body_load"];quad=res["body_quadrature"]
        errors=dict(full_nodal=normalized(np.array(res["force_n_xyz"])-sum_force,np.linalg.norm(sum_force)),
                    reference_convergence=normalized(a["f"]-b["f"],np.linalg.norm(b["f"])),
                    volume=abs(actual["volume_mm3"]/b["volume"]-1),
                    resultant=normalized(actual["resultant_n"]-b["force"],b["scale"]),
                    moment=normalized(actual["moment_n_mm"]-b["moment"],b["mscale"]),
                    work=normalized(actual["work_n_mm"]-b["f"]@sol,abs(b["f"]@sol)))
        direction=dict(schema="solidcraft_affine_velocity_v2",affine_mm=dict(denominator=1,
                    numerators=[[int(v),*row.tolist()] for v,row in zip(offset,B)]))
        write(dp,direction);solve_command(cli,"differentiate-job",mp,jp,dp,op);der=read(op)
        work_scale=max(1e-4,abs(der["applied_work_rate_tangent_n_mm"]))
        errors["applied_work_adjoint"]=abs(der["applied_work_rate_tangent_n_mm"]-der["applied_work_rate_adjoint_n_mm"])/work_scale
        errors["work_rate_balance"]=abs(der["work_rate_closure_n_mm"])/work_scale
        errors["explicit_work_rate"]=abs(der["applied_load_explicit_work_rate_n_mm"]-np.dot(der["nodal_force_rate_n_xyz"],sol))/work_scale
        errors["measure_rate"]=abs(der["volume_rate_mm3"]-np.trace(B)*der["volume_mm3"])/max(1.,b["volume"])
        same_rule=der["body_quadrature"]==quad
        if not load:errors["analytic_force_rate"]=normalized(np.asarray(der["nodal_force_rate_n_xyz"])-b["rate"],np.linalg.norm(b["rate"]))
        fields=np.r_[der["nodal_force_rate_n_xyz"],der["physical_displacement_rate_mm_xyz"],der["applied_work_rate_tangent_n_mm"]]
        sz=len(sol);hist=[];previous=None;fourth=None
        for den in (8,16,32):
            step=Fraction(1,den);vals=[];orders=[]
            for sign in (1,-1):
                write(rp,image(source,1,B,offset,sign*step));solve_command(cli,"build",rp,mp)
                solve_command(cli,"solve-job",mp,jp,op);v=read(op)
                vals.append(np.r_[v["force_n_xyz"],v["displacement_mm_xyz"],v["body_load"]["work_n_mm"]+v["traction_work_n_mm"]])
                orders.append(v["body_quadrature"]["selected_order"])
            rate=(vals[0]-vals[1])/(2*float(step))
            def err(g):
                return [float(np.linalg.norm(g[sl]-fields[sl])/max(1e-4,np.linalg.norm(fields[sl])))
                        for sl in (slice(0,sz),slice(sz,2*sz),slice(2*sz,None))]
            item=dict(step=1/den,orders=orders,centered=err(rate))
            if previous is not None:
                r4=(4*rate-previous)/3;item["fourth"]=err(r4)
                if fourth is not None:item["sixth"]=err((16*r4-fourth)/15)
                fourth=r4
            previous=rate;hist.append(item)
        accounting=sum(level["order"]**3*len(res["analysis_blocks"]) for level in quad["history"])==quad["points"]
        row=dict(id=name,blocks=len(res["analysis_blocks"]),body_quadrature=quad,independent=qinfo(b),
                 errors=errors,fixed_primal_rule=same_rule,point_accounting=accounting,derivative_history=hist,
                 passed=bool(same_rule and accounting and max(errors.values())<1e-8 and max(hist[-1]["sixth"])<1e-5))
        rows.append(row);print(name,row["passed"],quad["selected_order"],flush=True)
        if name=="extreme_surface_body":
            write(rp,source);solve_command(cli,"build",rp,mp)
            for label,overrides in (("point_budget",dict(max_points=quad["points"]-1)),
                                    ("order_budget",dict(max_order=3)),("invalid_range",dict(min_order=16))):
                bad=copy.deepcopy(job);bad["body_quadrature"].update(overrides)
                write(jp,bad);op.write_bytes(b"previous")
                proc=subprocess.run([cli,"solve-job",mp,jp,op],capture_output=True,text=True)
                refusals.append(dict(id=label,diagnostic=proc.stdout.strip(),passed=proc.returncode!=0 and op.read_bytes()==b"previous"))
    return dict(schema="adaptive_body_quadrature_oracle_v1",cases=rows,refusals=refusals,
                qualified=all(r["passed"] for r in [*rows,*refusals]))


def qinfo(b):
    return dict(order=24,volume_mm3=b["volume"],resultant_n=b["force"].tolist(),moment_n_mm=b["moment"].tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for v in ("cli","work","output"):p.add_argument(v,type=Path)
    a=p.parse_args();result=run(a.cli.resolve(),a.work);write(a.output,result)
    raise SystemExit(int(not result["qualified"]))
