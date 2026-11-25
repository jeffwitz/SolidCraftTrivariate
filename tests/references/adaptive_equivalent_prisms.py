"""Evaluate paired adaptive Q2 meshes driven by a cross-partition energy indicator."""
import argparse
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

from equivalent_prism_fields import Field, Geometry, clouds, compare, energy
from equivalent_prism_mechanics import read, save
from refine_tetra_fields import measure

HERE=Path(__file__).resolve().parent
METHODS=("boxes","contour")


def mark(fields,batches):
    indicators=[np.zeros(len(f.names)) for f in fields]
    for _,vol,x,locs in batches:
        (_,a),(_,b)=[f.evaluate(loc,x) for f,loc in zip(fields,locs)]
        delta=energy(a-b)*vol/(8*len(x))
        for scores,f,loc in zip(indicators,fields,locs):
            scores+=np.bincount(f.cell_ids(loc),weights=delta,minlength=len(scores))
    requests=[]
    stats=[]
    for scores,field in zip(indicators,fields):
        total=scores.sum()
        ranked=sorted((i for i,key in enumerate(field.keys) if key[1]<6 and scores[i]>0),
                      key=lambda i:(-scores[i],field.names[i]))
        selected=[]
        mass=0.
        for i in ranked:
            if len(selected)==256 or mass>=.5*total:break
            root,level,index=field.keys[i]
            selected.append(dict(root=root,level=level,index=index.tolist()))
            mass+=scores[i]
        requests.append(selected)
        stats.append(dict(selected=len(selected),captured_indicator_fraction=mass/total if total else 0.,
                          integrated_difference_n_mm=float(total),largest_leaf_fraction=float(scores.max()/total) if total else 0.))
    return requests,stats


def metrics(fields,batches,counts=(128,512,2048)):
    rows=compare(fields,batches,False,counts)["metrics"]
    final,previous=rows[-1],rows[-2]
    work=[f.physical["external_work_n_mm"] for f in fields]
    checks=dict(energy=final["energy_distance"]<=.05,work=abs(work[0]/work[1]-1)<=.01,
                recovery=max(map(abs,final["relative_energy_recovery"]))<=.015)
    for k,width in zip(("energy_distance","l2_distance"),final["empirical_half_width_95"]):
        tol=max(.001,.05*final[k])
        checks[k+"_stable"]=width<=tol and abs(previous[k]-final[k])<=tol
    return dict(levels=rows,checks=checks,qualified=all(checks.values()),work_n_mm=work)


def solve(cli,model,job,folder,base):
    folder.mkdir(parents=True,exist_ok=True);jp=folder/"job.json";out=folder/"result.json";save(jp,job)
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
    row=measure(list(map(str,[cli,"solve-job",model,jp,out])),env,folder/"solve.log")
    if row["status"]!="completed":return row
    pack=measure([sys.executable,str(HERE/"equivalent_prism_mechanics.py"),"compact",str(out),str(folder)],env,folder/"compact.log")
    row["compact"]=pack
    if pack["status"]!="completed":row["status"]="field_unavailable";return row
    physical=read(folder/"physical.json")
    errors={k:float(np.linalg.norm(np.array(physical[k])-np.array(base[k]))/max(np.linalg.norm(base[k]),1.))
            for k in ("volume_mm3","load_area_mm2","applied_resultant_n","applied_moment_n_mm")}
    row.update(physical=physical,relative_to_h2_geometry_and_loading=errors)
    if max(*errors.values(),physical["relative_force_balance"],physical["relative_moment_balance"])>1e-8:
        row["status"]="physical_checks_failed"
    return row


def campaign(cli,base,work,case,eliminate=False):
    ref=base/case;root=work/case;root.mkdir(parents=True,exist_ok=True)
    geoms=[Geometry(ref/tag) for tag in METHODS]
    held=clouds(geoms);train=clouds(geoms,seed=503430912,log2=9)
    locations=[ref/tag/"h2-q4" for tag in METHODS]
    jobs=[read(d/"job.json") for d in locations];bases=[read(d/"physical.json") for d in locations]
    for job in jobs:
        job["plan"].update(max_cells=24576,max_depth=6)
        if eliminate:job["schema"]="solidcraft_mechanical_job_v3"
    report=dict(case=case,rounds=[],qualified=False,training_seed=503430912,validation_seed=20261125,confirmation_seed=9031175)
    for count in range(7):
        fields=[Field(d,g,None) for d,g in zip(locations,geoms)]
        row=dict(round=count,leaves=[len(f.names) for f in fields],field=metrics(fields,held))
        report["rounds"].append(row);save(root/"report.json",report)
        print(case,count,row["leaves"],row["field"]["levels"][-1],flush=True)
        if row["field"]["qualified"]:
            confirmation=clouds(geoms,seed=9031175,log2=13)
            row["confirmation"]=metrics(fields,confirmation,(512,2048,8192))
            del confirmation
            if row["confirmation"]["qualified"]:
                report["qualified"]=True;report["stop"]="confirmed";break
        if count==6:
            report["stop"]="six_round_budget";break
        requests,row["marking"]=mark(fields,train)
        del fields
        if not all(requests):
            report["stop"]="no_admissible_marked_leaves";break
        row["next_solves"]=[]
        for i,tag in enumerate(METHODS):
            job=jobs[i];job["plan"]["steps"].append(dict(refine=requests[i]))
            locations[i]=root/tag/f"r{count+1}"
            status=solve(cli,ref/tag/"model.json",job,locations[i],bases[i])
            row["next_solves"].append(status);save(root/"report.json",report)
            print(case,count+1,tag,status["status"],round(status["elapsed_s"],1),flush=True)
        if any(s["status"]!="completed" for s in row["next_solves"]):
            report["stop"]="incomplete_paired_solve";break
    save(root/"report.json",report)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","base","work","output"):p.add_argument(k,type=Path)
    p.add_argument("--essential-elimination",action="store_true")
    a=p.parse_args();start=time.monotonic();reports=[]
    for case in ("L","two_holes"):
        reports.append(campaign(a.cli.resolve(),a.base,a.work,case,a.essential_elimination))
        save(a.output,dict(schema="paired_adaptive_physical_energy_discrepancy_v1",cases=reports,elapsed_s=time.monotonic()-start))
