"""Analytic, rational and Open CASCADE oracles for counterbored through holes."""
import FreeCAD as App
import Part
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess
import time
import platform

import numpy as np

from circular_polygon import incidence, reference
from coaxial_program import compile_mass_oracle
from polygon_prism import interior

HERE=Path(__file__).resolve().parent


def extra(base):
    rng=random.Random(2026112523)
    for i in range(96):
        r=copy.deepcopy(base["four_alternating"]["recipe"])
        r["bores"]=r["bores"][:1+i%4]
        r["height_mm"]=h=rng.randint(2,24)
        r["top_shift_mm"]=[rng.randint(-14,14),rng.randint(-14,14)]
        for bore in r["bores"]:
            bore["depth_mm"]=rng.randrange(1,h)
            bore["opening"]=rng.choice(("start","end"))
        ok=i%6!=0
        if not ok:
            key,value=("depth_mm",h) if i%12 else ("counterbore_radius_mm",0)
            r["bores"][0][key]=value
        yield dict(id=f"seed-{i:03d}",expected=ok,recipe=r)


def spans(b,h):
    d=b["depth_mm"];r=b["radius_mm"];R=b["counterbore_radius_mm"]
    return [(0,d,R),(d,h,r)] if b["opening"]=="start" else [(0,h-d,r),(h-d,h,R)]


def interval_recipe(r):
    cuts=[dict(b,radius_mm=radius,z_interval_mm=[a,z]) for b in r["bores"]
          for a,z,radius in spans(b,r["height_mm"])]
    return dict(r,bores=cuts)


def closed_form(r):
    p,twice,area,center=reference(dict(r,bores=[]))
    h=r["height_mm"];volume=area*h;moment=volume*center
    coeff=0
    for b in r["bores"]:
        pilot=b["radius_mm"]**2;annulus=b["counterbore_radius_mm"]**2-pilot
        d=b["depth_mm"];z=d/2 if b["opening"]=="start" else h-d/2
        for c,zc in [(pilot*h,h/2),(annulus*d,z)]:
            xy=np.array(b["center_mm"])+np.array(r["top_shift_mm"])*zc/h
            coeff+=c;moment-=math.pi*c*np.r_[xy,zc];volume-=math.pi*c
    return p,twice*h,coeff,volume,moment/volume,0,len(r["bores"])


def independent_cad(p,r):
    from interval_bores import independent_cad as cut_cylinders
    return cut_cylinders(p,interval_recipe(r))


def classify(r,p,shapes):
    physical=r; r=interval_recipe(r)
    h=r["height_mm"];shift=np.asarray(r["top_shift_mm"])
    rng=np.random.default_rng(2026112523)
    pts=[(*rng.uniform(p.min(axis=0)-1,p.max(axis=0)+1),rng.uniform(-.1,h+.1)) for _ in range(384)]
    cuts=sorted({0,h,*[z for b in r["bores"] for z in b["z_interval_mm"]]})
    for b in r["bores"]:
        center=np.asarray(b["center_mm"]);radius=b["radius_mm"]
        for lo,hi in zip(cuts,cuts[1:]):
            for extra_r in (0,radius-.0001,radius+.0001):
                for t in (.31,1.47,4.39):
                    xy=center+extra_r*np.array([math.cos(t),math.sin(t)])
                    pts.append((*xy,(lo+hi)/2))
        for z in b["z_interval_mm"]:
            for dz in (-.0001,.0001):pts.append((*center,z+dz))
    for b in physical["bores"]:
        shoulder=spans(b,h)[0][1]
        for rad in [b["radius_mm"]-.0001,b["radius_mm"]+.0001,
                    (b["radius_mm"]+b["counterbore_radius_mm"])/2,
                    b["counterbore_radius_mm"]-.0001,b["counterbore_radius_mm"]+.0001]:
            for dz in [-.0001,.0001]:
                for angle in [.73,3.63]:
                    xy=np.array(b["center_mm"])+rad*np.array([math.cos(angle),math.sin(angle)])
                    pts.append((*xy,shoulder+dz))
    bad=[]
    for x,y,z in pts:
        good=bool(0<z<h and interior(p,x,y) and not any(
            b["z_interval_mm"][0]<z<b["z_interval_mm"][1] and math.dist([x,y],b["center_mm"])<b["radius_mm"] for b in r["bores"]))
        world=App.Vector(*(np.array([x,y])+shift*z/h),z)
        vals=[bool(s.isInside(world,1e-8,False)) for s in shapes]
        if any(v!=good for v in vals):bad.append(dict(xyz=list(world),expected=good,observed=vals))
    return len(pts),bad


def walls(model,r):
    by_id={b["id"]:b for b in model["blocks"]}
    bores={b["id"]:b for b in r["bores"]}
    count=0;tags={}
    for f in model["design_faces"]:
        feature=f["feature"]
        if not feature.startswith("bore/"):continue
        _,key,kind=feature.split("/")
        tags[feature]=tags.get(feature,0)+1
        b=bores[key]
        if kind=="shoulder":
            side=f["side"];assert side==("w0" if b["opening"]=="start" else "w1")
            controls=np.array(by_id[f["block"]]["controls"]).reshape(3,2,2,4)[:,:,int(side[1])]
            z=spans(b,r["height_mm"])[0][1]
            assert np.all(controls[:,:,2]==z*controls[:,:,3]);continue
        assert kind in ("counterbore_wall","pilot_wall")
        rad=b["radius_mm"] if kind=="pilot_wall" else b["counterbore_radius_mm"]
        lo,hi=next((lo,hi) for lo,hi,rr in spans(b,r["height_mm"]) if rr==rad)
        src=by_id[f["block"]]
        assert f["side"]=="v1" and src["degrees"]==[2,1,1]
        net=np.array(src["controls"],dtype=object).reshape(3,2,2,4)[:,1]
        for u in (Q(0),Q(1,5),Q(1,2),Q(2,3),Q(1)):
            for w in (Q(0),Q(1,2),Q(1)):
                bu=[(1-u)**2,2*u*(1-u),u*u]
                xyzw=sum(t*(layer[0]*(1-w)+layer[1]*w) for t,layer in zip(bu,net))
                x,y,z=[v/xyzw[3] for v in xyzw[:3]]
                dx,dy=r["top_shift_mm"];cx,cy=b["center_mm"]
                assert (x-cx-dx*z/r["height_mm"])**2+(y-cy-dy*z/r["height_mm"])**2==rad**2
                assert lo<=z<=hi
                count+=1
    for b in bores.values():
        prefix="bore/"+b["id"]+"/"
        assert tags.get(prefix+"pilot_wall",0)>=8 and tags.get(prefix+"counterbore_wall",0)>=8
        assert tags.get(prefix+"shoulder",0)==8
    return count,tags


def pressure(cli,mp,folder,r):
    bore=r["bores"][0]
    prefix="bore/"+bore["id"]+"/";pres=2.5
    rec=[]
    for names in (["shoulder"],):
        job=dict(schema="solidcraft_mechanical_job_v3",material=dict(young_mpa=1e5,poisson=.25),formulation="displacement",order=10,
                 plan=dict(schema="adaptive_refinement_plan_v1",max_depth=0,max_cells=2048,steps=[]),
                 supports=[dict(select=dict(feature="start"),components=[0,1,2])],
                 loads=[dict(select=dict(feature=prefix+k),pressure_mpa=pres) for k in names])
        jp=folder/"pressure-job.json";out=folder/"pressure-result.json";jp.write_text(json.dumps(job))
        p=subprocess.run([cli,"solve-job",mp,jp,out],capture_output=True,text=True)
        if p.returncode:
            rec.append(dict(selection=names,passed=False,diagnostic=p.stdout));continue
        d=json.loads(out.read_text());scale=pres*math.pi*(bore["counterbore_radius_mm"]**2-bore["radius_mm"]**2)
        z=spans(bore,r["height_mm"])[0][1];xy=np.array(bore["center_mm"])+np.asarray(r["top_shift_mm"])*z/r["height_mm"]
        expected=np.array([0.,0.,scale if bore["opening"]=="start" else -scale])
        m=np.cross(np.r_[xy,z],expected)
        ef=np.linalg.norm(d["applied_resultant_n"]-expected)/scale
        em=np.linalg.norm(d["applied_moment_n_mm"]-m)/(scale*max(1.,np.linalg.norm(np.r_[xy,z])))
        rec.append(dict(selection=names,resultant_n=d["applied_resultant_n"],moment_n_mm=d["applied_moment_n_mm"],
                        normalized_force_error=ef,normalized_moment_error=em,
                        strain_energy_n_mm=d["strain_energy_n_mm"],passed=bool(max(ef,em,d["relative_force_balance"],d["relative_moment_balance"])<1e-8 and d["strain_energy_n_mm"]>0)))
    return rec


def run(cli,work):
    start=time.monotonic();binary_sha=hashlib.sha256(cli.read_bytes()).hexdigest()
    work.mkdir(exist_ok=True,parents=True);helper=compile_mass_oracle(work)
    corpus=json.loads((HERE/"counterbores_corpus.json").read_text());base={i["id"]:i for i in corpus}
    rows=[];hashes={};rp=work/"recipe.json";mp=work/"model.json"
    for item in [*corpus,*extra(base)]:
        r=item["recipe"];name=item["id"];rp.write_text(json.dumps(r))
        old=mp.read_bytes() if mp.exists() else None
        before=time.monotonic()
        proc=subprocess.run([cli,"build",rp,mp],capture_output=True,text=True)
        info=json.loads(proc.stdout)
        row=dict(id=name,recipe=r,expected=item["expected"],audit=info,build_and_audit_s=time.monotonic()-before);rows.append(row)
        ok=(proc.returncode==0)==item["expected"]
        if proc.returncode:
            row["preserved"]=(mp.read_bytes() if mp.exists() else None)==old;ok &= row["preserved"]
        elif item["expected"]:
            model=json.loads(mp.read_text());poly,rational,coef,vol,cent,voids,tunnels=closed_form(r)
            chi,mult=incidence(model)
            digest=hashlib.sha256(json.dumps([model[k] for k in ("blocks","topology","design_faces")],sort_keys=True).encode()).hexdigest()
            hashes[name]=digest
            expected_volume=dict(twice_rational_part_mm3=rational,pi_coefficient_mm3=-coef)
            row.update(volume_mm3=vol,centroid_mm=cent.tolist(),euler=chi,voids=voids,tunnels=tunnels,partition_digest=digest)
            ok &= (info["volume_exact"]==expected_volume and info["cavities"]==voids and info["through_bores"]==tunnels and chi==1-tunnels+voids and mult==2)
            if "equivalent_to" in item:ok &= digest==hashes[item["equivalent_to"]]
            elif not name.startswith("seed-"):
                row["exact_wall_points"],row["bore_tags"]=walls(model,r)
                before=time.monotonic()
                sp=work/(name+".step");proc=subprocess.run([cli,"export-step",mp,sp],capture_output=True,text=True)
                row["step_export_s"]=time.monotonic()-before
                row["step"]=json.loads(proc.stdout);ok &= proc.returncode==0
                if proc.returncode==0:
                    before=time.monotonic()
                    shapes=[independent_cad(poly,r),Part.read(str(sp))]
                    masses=[]
                    for index,shape in enumerate(shapes):
                        bp=work/f"{name}-{index}.brep";shape.exportBrep(str(bp))
                        mass=json.loads(subprocess.check_output([helper,bp,"--centroid"],text=True))[-1];masses.append(mass)
                    n,bad=classify(r,poly,shapes)
                    ev=max(abs(m["volume_mm3"]/vol-1) for m in masses)
                    ec=max(float(np.linalg.norm(m["centroid_mm"]-cent)) for m in masses)
                    row.update(masses=masses,classifications=n,mismatches=len(bad),witnesses=bad[:6],volume_error=ev,centroid_error_mm=ec,
                               shell_counts=[len(s.Shells) for s in shapes],valid=[bool(s.isValid() and len(s.Solids)==1) for s in shapes])
                    bounds=np.r_[np.minimum(poly.min(axis=0),poly.min(axis=0)+r["top_shift_mm"]),0,
                                 np.maximum(poly.max(axis=0),poly.max(axis=0)+r["top_shift_mm"]),r["height_mm"]]
                    def box(shape):
                        bb=shape.optimalBoundingBox(False)
                        return np.array([bb.XMin,bb.YMin,bb.ZMin,bb.XMax,bb.YMax,bb.ZMax])
                    row["bounds_error_mm"]=max(float(np.max(np.abs(box(s)-bounds))) for s in shapes)
                    row["independent_cad_and_reimport_validation_s"]=time.monotonic()-before
                    ok &= row["bounds_error_mm"]<1e-6
                    ok &= (ev<1e-8 and ec<1e-7 and not bad and all(row["valid"]) and all(c==1+voids for c in row["shell_counts"]))
                if name in ("from_start","from_end"):
                    row["pressure"]=pressure(cli,mp,work,r);ok &= all(p["passed"] for p in row["pressure"])
        row["qualified"]=bool(ok)
        if not name.startswith("seed-") or len(rows)%16==0:print(name,bool(ok),flush=True)
    if binary_sha!=hashlib.sha256(cli.read_bytes()).hexdigest():raise ValueError("executable changed during qualification")
    return dict(schema="counterbore_prism_cad_oracle_v1",cli_sha256=binary_sha,python=platform.python_version(),elapsed_s=time.monotonic()-start,occ=Part.OCC_VERSION,seed=2026112523,cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","work","output"):p.add_argument(k,type=Path)
    args=p.parse_args();report=run(args.cli.resolve(),args.work)
    args.output.write_text(json.dumps(report,allow_nan=False,indent=2)+"\n")
    raise SystemExit(int(not all(r["qualified"] for r in report["cases"])))
