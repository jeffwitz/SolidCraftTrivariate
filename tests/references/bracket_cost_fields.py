"""Compare J4r tiers to both archived and refined C3D10 fields on common physical points."""
import argparse
import json
import os
from pathlib import Path
import time

import numpy as np
from scipy.stats import qmc, t as student
from adaptive_field_distance import AdaptiveField
from bracket_field_distance import TetraField, ImmersedField, density, verify_affine_readers
from bracket_serial_cost import first_local_tiers
from refine_tetra_fields import measure
from solver_package import shape

ROOT=Path(__file__).resolve().parents[2]


def read(path):
    return json.loads(path.read_text())


def save(path,data):
    tmp=path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    tmp.replace(path)


def summarize(rows,counts,scales,pairs,nblock):
    out=[]
    for n in counts:
        per=rows[n]
        def array(key):
            return np.array([r[key] for r in per])
        energy={k:float(array(k+"_energy").mean()) for k in scales}
        distances={}
        for a,b in pairs:
            pair=f"{a}__{b}"
            entry={}
            for measure_name in ("energy","l2"):
                ratios=np.sqrt(array(pair+"_"+measure_name)/array(b+"_"+measure_name))
                avg=float(ratios.mean())
                ci=float(student.ppf(.975,7)*ratios.std(ddof=1)/np.sqrt(8))
                entry[measure_name]={"mean":avg,"half_width_95":ci}
            distances[pair]=entry
        out.append(dict(samples_per_scramble_per_block=n,points_total=n*8*nblock,
                        full_volume_mm3=float(array("volume").mean()),
                        excluded_volume_fraction=float(array("outside").sum()/array("volume").sum()),
                        measured_energy_n_mm=energy,recovered_energy_relative={k:abs(energy[k]/v-1) for k,v in scales.items()},
                        distances=distances))
    return out


def compare(blocks,fields,scales,log2):
    counts=[n for n in (128,512,2048,8192) if n<=2**log2]
    totals={n:[{} for _ in range(8)] for n in counts}
    tetra={k:v for k,v in fields.items() if isinstance(v,TetraField)}
    pairs=[(a,b) for b in ("tetra_h125","tetra_h1") if b in fields for a in fields if a!=b and (a!="tetra_h1")]
    pairs=list(dict.fromkeys(pairs))
    checks={k:{"not_located":0,"multiple":0,"max_affine_inverse_error":0.} for k in tetra}
    outside_witness=[]
    start=time.monotonic()
    for bi,block in enumerate(sorted(blocks,key=lambda b:b["id"])):
        for rep in range(8):
            uvw=qmc.Sobol(3,scramble=True,seed=20261125+1000*bi+rep).random_base2(log2)
            _,x,j,_=shape(uvw,block["controls"])
            det=np.linalg.det(j);inv=np.linalg.inv(j)
            available=np.ones(len(uvw),dtype=bool);values={}
            for k,f in fields.items():
                if k in tetra:
                    u,g,located,multi,err=f.evaluate(x)
                    info=checks[k];info["not_located"]+=int(np.count_nonzero(~located))
                    info["multiple"]+=int(multi);info["max_affine_inverse_error"]=max(info["max_affine_inverse_error"],float(err))
                    if len(outside_witness)<16:
                        outside_witness.extend(x[~located][:16-len(outside_witness)].tolist())
                    available &= located;values[k]=(u,g)
                elif isinstance(f,AdaptiveField):values[k]=f.evaluate(block["id"],uvw,inv)
                else:values[k]=f.evaluate(x)
            raw={"volume":det,"outside":det*(~available)}
            for k,(u,g) in values.items():
                raw[k+"_energy"]=np.where(available,det*density(g)/2,0.)
                raw[k+"_l2"]=np.where(available,det*np.einsum("ij,ij->i",u,u),0.)
            for a,b in pairs:
                u,g=values[a];v,h=values[b];key=f"{a}__{b}"
                raw[key+"_energy"]=np.where(available,det*density(g-h)/2,0.)
                raw[key+"_l2"]=np.where(available,det*np.sum((u-v)**2,axis=1),0.)
            for n in counts:
                row=totals[n][rep]
                for key,v in raw.items():row[key]=row.get(key,0.)+float(np.mean(v[:n]))
        if bi%8==7:print("block",bi+1,flush=True)
    levels=summarize(totals,counts,scales,pairs,len(blocks))
    previous,final=levels[-2:]
    stable=True
    for p,metrics in final["distances"].items():
        for name,detail in metrics.items():
            allowance=max(.001,.05*detail["mean"])
            change=abs(detail["mean"]-previous["distances"][p][name]["mean"])
            detail.update(previous_level_change=change,stable=change<=allowance and detail["half_width_95"]<=allowance)
            stable &= detail["stable"]
    return dict(elapsed_s=time.monotonic()-start,inverse_checks=checks,excluded_locations_mm=outside_witness,
                reader_checks=verify_affine_readers(blocks),levels=levels,
                recovered_energies=max(final["recovered_energy_relative"].values())<=.015,
                stable=bool(stable),all_points_located_once=all(c["not_located"]==0 and c["multiple"]==0 for c in checks.values()))


def one(work,rx,ry,tier,log2):
    folder=work/f"x{rx}-y{ry}"
    source=ROOT/f"build/bracket-cost/x{rx}-y{ry}/tetra"
    inputs=next(r for r in read(work/"preparation.json") if r["radii_ticks"]==[rx,ry])
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
    cli=ROOT/"target/release/solidcraft-trivariate"
    recipe=folder/"native.recipe.json";model=folder/"native.model.json";package=folder/"native.package.json"
    save(recipe,dict(family="parameterized_bracket_v1",radii_sixteenths_mm=[rx,ry]))
    for name,cmd in (("geometry",[cli,"build",recipe,model]),("package",[cli,"export-analysis",model,8,package])):
        m=measure(list(map(str,cmd)),env,folder/(name+".log"))
        if m["status"]!="completed":raise ValueError(f"{name}: {m}")
    blocks=read(package)["analysis_space"]["blocks"]
    fields={};scales={}
    if tier:
        result=folder/"native-tier.json";checkpoint=folder/"native-tier.bin"
        cmd=[ROOT/"target/release/examples/local_fillet_cost","--tier",rx,ry,tier["path"],tier["step"],result,checkpoint]
        m=measure(list(map(str,cmd)),env,folder/"local-checkpoint.log")
        if m["status"]!="completed":raise ValueError(m)
        fields["native_local"]=AdaptiveField(checkpoint)
        scales["native_local"]=read(result)["q8"]["strain_energy_n_mm"]
    fields["immersed_h125"]=ImmersedField(folder/"immersed.npz")
    scales["immersed_h125"]=read(folder/"immersed.json")["strain_energy_n_mm"]
    for name,loc in (("tetra_h15",source/"h1.5-qNone"),("tetra_h125",source/"h1.25-qNone"),("tetra_h1",folder/"h1")):
        if name=="tetra_h1" and not inputs.get("fine_qualified"):continue
        fields[name]=TetraField(loc)
        d=read(loc/"result.json");scales[name]=d.get("strain_energy_n_mm",-500.*d["mean_tip_y_mm"])
    return dict(radii_ticks=[rx,ry],local_tier=tier,reference_energy_n_mm=scales,preparation=inputs,**compare(blocks,fields,scales,log2))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("work",type=Path);p.add_argument("output",type=Path)
    p.add_argument("rx",type=int);p.add_argument("ry",type=int);p.add_argument("--log2",type=int,choices=(9,11,13),default=11)
    p.add_argument("--serial-report",type=Path,help="Read the local tier actually replayed by J4r")
    a=p.parse_args()
    if a.serial_report:
        serial=read(a.serial_report)
        rows=serial.get("schedule_rows",serial.get("completed_schedule_rows",[]))
        candidates=[t["j4d_tier"] for r in rows if r["radii_ticks"]==[a.rx,a.ry]
                    for t in r["trials"] if t["method"]=="local" and t.get("status")=="completed"]
        if not candidates or any(c!=candidates[0] for c in candidates):
            raise ValueError("Missing or inconsistent frozen local tier in serial trials")
        tier=candidates[0]
    else:
        tier=first_local_tiers(read(ROOT/"docs/results/local-fillet-cost-20261125.json")).get((a.rx,a.ry))
    save(a.output,one(a.work.resolve(),a.rx,a.ry,tier,a.log2))
