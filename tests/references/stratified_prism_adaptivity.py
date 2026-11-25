"""Held-out integration on fixed physical strata for two independent L-prism meshes."""
import argparse
from pathlib import Path
import time

import numpy as np
from scipy.stats import qmc, t as student

from adaptive_equivalent_prisms import METHODS, solve
from equivalent_prism_fields import Field, Geometry, energy
from equivalent_prism_mechanics import read, save


def strata(field):
    xyz = field.geom[..., :3] / field.geom[..., 3, None]
    lo, hi = xyz.min(axis=1), xyz.max(axis=1)
    corners = xyz.reshape(-1, 3, 3, 3, 3)[:, ::2, ::2, ::2].reshape(-1, 8, 3)
    assert np.all((corners == lo[:,None,:]) | (corners == hi[:,None,:]))
    volume = np.prod(hi-lo, axis=1)
    assert np.all(volume>0) and abs(volume.sum()/field.physical["volume_mm3"]-1) < 1e-10
    # Field's per-root dyadic occupancy check guarantees a gap-free partition.
    return lo, hi, volume


def integrate(fields, geometries, boxes, seed, counts, marking=False):
    lo, hi, volumes = boxes
    out = np.zeros((len(counts),8,5))
    scores = [np.zeros(len(f.names)) for f in fields]
    n = counts[-1]
    start = time.monotonic()
    for batch in range(0,len(lo),16):
        last = min(batch+16, len(lo))
        weights = volumes[batch:last]
        for rep in range(8):
            sample = np.stack([qmc.Sobol(3, scramble=True, seed=seed+8*k+rep).random_base2(n.bit_length()-1)
                               for k in range(batch,last)])
            x = (lo[batch:last,None,:]+(hi-lo)[batch:last,None,:]*sample).reshape(-1,3)
            loc = [g.locate(x) for g in geometries]
            (u, g), (v, h) = [f.evaluate(t,x) for f,t in zip(fields,loc)]
            diff = energy(g-h)
            values = np.column_stack((energy(g),energy(h),diff,np.sum((u-v)**2,axis=1),np.sum(v*v,axis=1))).reshape(last-batch,n,5)
            for j,count in enumerate(counts):
                out[j,rep] += weights @ values[:,:count].mean(axis=1)
            if marking:
                mass = (diff.reshape(last-batch,n)*weights[:,None]/(8*n)).ravel()
                for score, f, t in zip(scores,fields,loc):
                    score += np.bincount(f.cell_ids(t), weights=mass, minlength=len(score))
    metrics = []
    for count, values in zip(counts,out):
        ratios = np.sqrt(values[:,[2,3]]/values[:,[1,4]])
        mean = ratios.mean(axis=0)
        ci = student.ppf(.975,7)*ratios.std(axis=0,ddof=1)/np.sqrt(8)
        recovery = values[:,:2].mean(axis=0)/[f.physical["strain_energy_n_mm"] for f in fields]-1
        metrics.append(dict(points_per_stratum_scramble=count, physical_points=count*8*len(lo),
                            energy_distance=float(mean[0]), l2_distance=float(mean[1]),
                            half_width_95=ci.tolist(), energy_recovery=recovery.tolist()))
    final, prev = metrics[-1], metrics[-2] if len(metrics)>1 else metrics[-1]
    stable = all(max(ci,abs(final[key]-prev[key]))<=max(.001,.05*final[key])
                 for ci,key in zip(final["half_width_95"],("energy_distance","l2_distance")))
    work = [f.physical["external_work_n_mm"] for f in fields]
    checks = dict(stable=stable,energy=final["energy_distance"]<=.05,work=abs(work[0]/work[1]-1)<=.01,
                  integrated_energy=max(map(abs,final["energy_recovery"]))<=.015)
    return dict(levels=metrics,work_n_mm=work,checks=checks,qualified=all(checks.values()),
                elapsed_s=time.monotonic()-start), scores


def requests(field, scores):
    ranked = sorted((i for i,key in enumerate(field.keys) if key[1]<8 and scores[i]>0),
                    key=lambda i:(-scores[i],field.names[i]))
    total = scores.sum()
    result = []
    mass = 0.
    for i in ranked:
        if len(result)==384 or mass>=.5*total: break
        root, level, ids = field.keys[i]
        result.append(dict(root=root,level=level,index=ids.tolist()))
        mass += scores[i]
    return result, dict(leaves=len(result),mass_fraction=float(mass/total))


def run(cli, base, previous, work, report):
    geometries = [Geometry(base/t) for t in METHODS]
    locations = [previous/t/"r6" for t in METHODS]
    physics = [read(p/"physical.json") for p in locations]
    jobs = [read(p/"job.json") for p in locations]
    for job in jobs:job["plan"].update(max_depth=8, max_cells=24576)
    anchor = Field(locations[0],geometries[0],None)
    boxes = strata(anchor)
    del anchor
    data = dict(schema="stratified_equivalent_prism_v1", seeds=dict(training=62503231,validation=183214505,confirmation=371663271),
                starting_meshes="M8c L r6", strata=len(boxes[0]), volume_mm3=float(boxes[2].sum()),rounds=[],qualified=False)
    for k in range(5):
        fields = [Field(p,g,None) for p,g in zip(locations,geometries)]
        row = dict(round=k, cells=[len(f.names) for f in fields])
        data["rounds"].append(row)
        row["validation"],_ = integrate(fields,geometries,boxes,data["seeds"]["validation"],(8,32,128))
        save(report,data)
        print("validation",k,row["cells"],row["validation"]["levels"][-1],flush=True)
        if row["validation"]["qualified"]:
            row["confirmation"],_ = integrate(fields,geometries,boxes,data["seeds"]["confirmation"],(16,64,256))
            if row["confirmation"]["qualified"]:
                data.update(qualified=True,stop="independently_confirmed");break
        if k==4:
            data["stop"]="four_new_round_budget";break
        row["training"],scores=integrate(fields,geometries,boxes,data["seeds"]["training"],(32,),True)
        chosen = [requests(f,s) for f,s in zip(fields,scores)]
        row["marking"]=[s[1] for s in chosen]
        del fields
        row["solves"]=[]
        for i,tag in enumerate(METHODS):
            jobs[i]["plan"]["steps"].append(dict(refine=chosen[i][0]))
            locations[i]=work/tag/f"r{k+1}"
            status=solve(cli,base/tag/"model.json",jobs[i],locations[i],physics[i])
            row["solves"].append(status)
            save(report,data)
            print("solve",k+1,tag,status["status"],round(status["elapsed_s"],2),flush=True)
        if any(s["status"]!="completed" for s in row["solves"]):
            data["stop"]="incomplete_paired_solve";break
    save(report,data)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("cli","base","previous","work","report"):p.add_argument(name,type=Path)
    a=p.parse_args()
    run(a.cli.resolve(),a.base,a.previous,a.work,a.report)
