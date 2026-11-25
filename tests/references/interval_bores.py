"""Analytic, rational and Open CASCADE oracles for finite circular excavations."""
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

import numpy as np

from circular_polygon import incidence, reference
from coaxial_program import compile_mass_oracle
from polygon_prism import interior

HERE=Path(__file__).resolve().parent


def extra(base):
    rng=random.Random(2026112521)
    for i in range(96):
        r=copy.deepcopy(base["mixed"]["recipe"])
        r["bores"]=r["bores"][:1+i%4]
        r["height_mm"]=h=rng.randint(4,24)
        r["top_shift_mm"]=[rng.randint(-14,14),rng.randint(-14,14)]
        for b in r["bores"]:b["z_interval_mm"]=sorted(rng.sample(range(h+1),2))
        ok=i%6!=0
        if not ok:
            if i%12==0:r["bores"][0]["z_interval_mm"][1]=h+1
            else:r["bores"][0]["guard_mm"]=0
        yield dict(id=f"seed-{i:03d}",expected=ok,recipe=r)


def closed_form(r):
    p,twice,area,center=reference(dict(r,bores=[]))
    h=r["height_mm"];volume=area*h
    moment=volume*center;shift=np.array(r["top_shift_mm"])
    removed=0;sealed=through=0
    for b in r["bores"]:
        lo,hi=b["z_interval_mm"];z=(lo+hi)/2
        coef=(hi-lo)*b["radius_mm"]**2;removed+=coef
        v=math.pi*coef;volume-=v
        moment-=v*np.r_[np.array(b["center_mm"])+shift*z/h,z]
        sealed+=int(0<lo and hi<h);through+=int(lo==0 and hi==h)
    return p,twice*h,removed,volume,moment/volume,sealed,through


def independent_cad(p,r):
    h=r["height_mm"];shift=np.array(r["top_shift_mm"],dtype=float)
    pts=[App.Vector(*map(float,xy),0) for xy in p]
    solid=Part.Face(Part.makePolygon([*pts,pts[0]])).extrude(App.Vector(*shift,h))
    for b in r["bores"]:
        lo,hi=b["z_interval_mm"];v=App.Vector(*(shift*(hi-lo)/h),hi-lo)
        center=np.array(b["center_mm"])+shift*lo/h
        ring=Part.makeCircle(b["radius_mm"],App.Vector(*center,lo))
        solid=solid.cut(Part.Face(Part.Wire([ring])).extrude(v))
    return solid.removeSplitter()


def classify(r,p,shapes):
    h=r["height_mm"];shift=np.asarray(r["top_shift_mm"])
    rng=np.random.default_rng(2026112521)
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
        if kind!="wall":continue
        src=by_id[f["block"]]
        assert f["side"]=="v1" and src["degrees"]==[2,1,1]
        net=np.array(src["controls"],dtype=object).reshape(3,2,2,4)[:,1]
        for u in (Q(0),Q(1,5),Q(1,2),Q(2,3),Q(1)):
            for w in (Q(0),Q(1,2),Q(1)):
                bu=[(1-u)**2,2*u*(1-u),u*u]
                xyzw=sum(t*(layer[0]*(1-w)+layer[1]*w) for t,layer in zip(bu,net))
                x,y,z=[v/xyzw[3] for v in xyzw[:3]]
                dx,dy=r["top_shift_mm"];cx,cy=b["center_mm"]
                assert (x-cx-dx*z/r["height_mm"])**2+(y-cy-dy*z/r["height_mm"])**2==b["radius_mm"]**2
                assert b["z_interval_mm"][0]<=z<=b["z_interval_mm"][1]
                count+=1
    for b in bores.values():
        prefix="bore/"+b["id"]+"/"
        assert tags.get(prefix+"wall",0)>=8
        for name,present in [("lower",b["z_interval_mm"][0]>0),("upper",b["z_interval_mm"][1]<r["height_mm"])]:
            assert tags.get(prefix+name,0)==(12 if present else 0)
    return count,tags


def pressure(cli,mp,folder,r):
    bore=r["bores"][0]
    prefix="bore/"+bore["id"]+"/";pres=2.5
    rec=[]
    for names in (["lower"],["lower","upper","wall"]):
        job=dict(schema="solidcraft_mechanical_job_v3",material=dict(young_mpa=1e5,poisson=.25),formulation="displacement",order=10,
                 plan=dict(schema="adaptive_refinement_plan_v1",max_depth=0,max_cells=2048,steps=[]),
                 supports=[dict(select=dict(feature="start"),components=[0,1,2])],
                 loads=[dict(select=dict(feature=prefix+k),pressure_mpa=pres) for k in names])
        jp=folder/"pressure-job.json";out=folder/"pressure-result.json";jp.write_text(json.dumps(job))
        p=subprocess.run([cli,"solve-job",mp,jp,out],capture_output=True,text=True)
        if p.returncode:
            rec.append(dict(selection=names,passed=False,diagnostic=p.stdout));continue
        d=json.loads(out.read_text());scale=pres*math.pi*bore["radius_mm"]**2
        z=bore["z_interval_mm"][0];xy=np.array(bore["center_mm"])+np.asarray(r["top_shift_mm"])*z/r["height_mm"]
        expected=np.array([0.,0.,-scale]) if len(names)==1 else np.zeros(3)
        m=np.cross(np.r_[xy,z],expected)
        ef=np.linalg.norm(d["applied_resultant_n"]-expected)/scale
        em=np.linalg.norm(d["applied_moment_n_mm"]-m)/(scale*max(1.,np.linalg.norm(np.r_[xy,z])))
        rec.append(dict(selection=names,resultant_n=d["applied_resultant_n"],moment_n_mm=d["applied_moment_n_mm"],
                        normalized_force_error=ef,normalized_moment_error=em,
                        strain_energy_n_mm=d["strain_energy_n_mm"],passed=bool(max(ef,em,d["relative_force_balance"],d["relative_moment_balance"])<1e-8 and d["strain_energy_n_mm"]>0)))
    return rec


def run(cli,work):
    work.mkdir(exist_ok=True,parents=True);helper=compile_mass_oracle(work)
    corpus=json.loads((HERE/"interval_bores_corpus.json").read_text());base={i["id"]:i for i in corpus}
    rows=[];hashes={};rp=work/"recipe.json";mp=work/"model.json"
    off_axis=copy.deepcopy(base["closed"]["recipe"])
    off_axis["top_shift_mm"]=[3,-2]
    off_axis["bores"][0]["center_mm"]=[2,-1]
    supplement=dict(id="closed_off_axis",recipe=off_axis,expected=True)
    for item in [*corpus,*extra(base),supplement]:
        r=item["recipe"];name=item["id"];rp.write_text(json.dumps(r))
        old=mp.read_bytes() if mp.exists() else None
        proc=subprocess.run([cli,"build",rp,mp],capture_output=True,text=True)
        info=json.loads(proc.stdout)
        row=dict(id=name,recipe=r,expected=item["expected"],audit=info);rows.append(row)
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
                sp=work/(name+".step");proc=subprocess.run([cli,"export-step",mp,sp],capture_output=True,text=True)
                row["step"]=json.loads(proc.stdout);ok &= proc.returncode==0
                if proc.returncode==0:
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
                    ok &= (ev<1e-8 and ec<1e-7 and not bad and all(row["valid"]) and all(c==1+voids for c in row["shell_counts"]))
                if name in ("closed","closed_off_axis"):
                    row["pressure"]=pressure(cli,mp,work,r);ok &= all(p["passed"] for p in row["pressure"])
        row["qualified"]=bool(ok)
        if not name.startswith("seed-") or len(rows)%16==0:print(name,bool(ok),flush=True)
    return dict(schema="interval_bore_prism_cad_oracle_v1",occ=Part.OCC_VERSION,seed=2026112521,cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","work","output"):p.add_argument(k,type=Path)
    args=p.parse_args();report=run(args.cli.resolve(),args.work)
    args.output.write_text(json.dumps(report,allow_nan=False,indent=2)+"\n")
    raise SystemExit(int(not all(r["qualified"] for r in report["cases"])))
