"""Run the two-program, same-domain M7 elasticity study under fixed resource limits."""
import argparse
from itertools import product
import json
import os
from pathlib import Path
import sys

import numpy as np

from refine_tetra_fields import measure

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent


def save(path,data):
    tmp=path.with_suffix(".tmp");tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n");tmp.replace(path)


def read(path):
    return json.loads(path.read_text())


def recipes(case):
    if case=="L":
        shape=next(c["recipe"] for c in read(HERE/"polygon_prism_corpus.json") if c["id"]=="L")
        shape=dict(shape,height_mm=4,top_shift_mm=[0,0])
        box={"family":"orthogonal_feature_program_v1","base_mm":[[0,12],[0,11],[0,4]],
             "features":[{"id":"notch","operation":"cut","bounds_mm":[[4,12],[3,11],[0,4]]}]}
        return box,shape,(0,12),3.*4,np.array([12.,1.5,2.])
    shape=next(c["recipe"] for c in read(HERE/"perforated_prism_corpus.json") if c["id"]=="two_staggered")
    shape=dict(shape,height_mm=4,top_shift_mm=[0,0])
    box={"family":"orthogonal_feature_program_v1","base_mm":[[-16,16],[-12,12],[0,4]],"features":[]}
    for loop in shape["holes"]:
        pos=np.array([v["xy_mm"] for v in loop["vertices"]])
        box["features"].append(dict(id=loop["id"],operation="cut",bounds_mm=[*[list(map(int,q)) for q in zip(pos.min(axis=0),pos.max(axis=0))],[0,4]]))
    return box,shape,(-16,16),24.*4,np.array([16.,0.,2.])


def plan(roots,h):
    steps=[]
    for level in range(h):
        entries=[dict(root=r,level=level,index=ijk) for r in roots for ijk in product(range(1<<level),repeat=3)]
        steps.append(dict(refine=entries))
    return dict(schema="adaptive_refinement_plan_v1",max_depth=h,max_cells=len(roots)*8**h,steps=steps)


def compact(result,folder):
    data=read(result)
    blocks=data["analysis_blocks"]
    coeff=np.array(data["displacement_mm_xyz"]).reshape(-1,3)
    geom=np.array([b["controls"] for b in blocks],dtype=float)
    table=np.array([b["dofs"] for b in blocks],dtype=int)
    np.savez_compressed(folder/"field.npz",names=np.array([b["id"] for b in blocks]),geom=geom,field=coeff[table])
    info={k:data[k] for k in ("strain_energy_n_mm","external_work_n_mm","relative_residual","relative_force_balance",
                              "relative_moment_balance","applied_resultant_n","applied_moment_n_mm","volume_mm3","scalar_dofs","iterations")}
    info.update(cells=len(blocks),load_area_mm2=sum(r["area_mm2"] for r in data["loads"]))
    if not data["loads"]:
        points=geom[:,:,:3]/geom[:,:,3,None]
        gradient=np.array([[1,1,0],[0,-2,1],[2,0,1]],dtype=float)/1000
        shift=np.array([3,-2,1],dtype=float)/1000
        info["affine_control_max_error_mm"]=float(np.max(np.linalg.norm(coeff[table]-(points@gradient.T+shift),axis=2)))
    save(folder/"physical.json",info)


def campaign(cli,work,case):
    roots=work/case;roots.mkdir(exist_ok=True,parents=True)
    box,shape,(xmin,xmax),load_area,center=recipes(case)
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
    source={}
    for tag,recipe in (("boxes",box),("contour",shape)):
        fold=roots/tag;fold.mkdir(exist_ok=True);rp=fold/"recipe.json";mp=fold/"model.json";save(rp,recipe)
        m=measure(list(map(str,[cli,"build",rp,mp])),env,fold/"build.log")
        if m["status"]!="completed":raise ValueError(f"{tag} model refused: {(fold/'build.log').read_text()}")
        source[tag]=(mp,sorted(b["id"] for b in read(mp)["blocks"]))
    rows=[]
    for h,q,affine in [(0,4,True),(0,4,False),(1,4,False),(2,4,False),(2,6,False),(3,4,False)]:
        for tag,(mp,ids) in source.items():
            name="affine" if affine else f"h{h}-q{q}"
            folder=roots/tag/name;folder.mkdir(exist_ok=True)
            jp=folder/"job.json";result=folder/"result.json"
            support={"select":{"feature":"*","plane":{"axis":0,"coordinate_mm":xmin}},"components":[0,1,2]}
            traction=[2.,-1.,.5]
            loads=[{"select":{"feature":"*","plane":{"axis":0,"coordinate_mm":xmax}},"vector_mpa":traction}]
            if affine:
                support={"select":{"feature":"*"},"components":[0,1,2],"affine_mm":{
                    "denominator":1000,"numerators":[[3,1,1,0],[-2,0,-2,1],[1,2,0,1]]}}
                loads=[]
            job=dict(schema="solidcraft_mechanical_job_v2",material=dict(young_mpa=100000,poisson=.25),
                     formulation="displacement",order=q,plan=plan(ids,h),supports=[support],loads=loads)
            save(jp,job)
            run=measure(list(map(str,[cli,"solve-job",mp,jp,result])),env,folder/"solve.log")
            row=dict(case=case,method=tag,level=h,order=q,affine=affine,solve=run)
            if run["status"]=="completed":
                pack=measure([sys.executable,str(Path(__file__).resolve()),"compact",str(result),str(folder)],env,folder/"compact.log")
                row["compact"]=pack
                if pack["status"]=="completed":
                    physical=read(folder/"physical.json");row["physical"]=physical
                    load=np.array(traction)*load_area
                    if affine:ok=physical["affine_control_max_error_mm"]<=1e-9
                    else:ok=max(abs(physical["load_area_mm2"]/load_area-1),float(np.linalg.norm(np.array(physical["applied_resultant_n"])-load)/np.linalg.norm(load)),
                                float(np.linalg.norm(np.array(physical["applied_moment_n_mm"])-np.cross(center,load))/max(1.,np.linalg.norm(np.cross(center,load)))))<1e-9
                    row["qualified_equilibrium_and_loading"]=bool(ok and max(physical["relative_force_balance"],physical["relative_moment_balance"])<1e-8)
            rows.append(row);save(roots/"solves.json",rows)
            print(case,tag,name,run["status"],row.get("qualified_equilibrium_and_loading"),flush=True)
    return rows


if __name__=="__main__":
    if len(sys.argv)==4 and sys.argv[1]=="compact":compact(Path(sys.argv[2]),Path(sys.argv[3]))
    else:
        p=argparse.ArgumentParser(description=__doc__);p.add_argument("cli",type=Path);p.add_argument("work",type=Path);p.add_argument("case",choices=("L","two_holes"))
        a=p.parse_args();campaign(a.cli.resolve(),a.work.resolve(),a.case)
