"""Check G20 circles, named surfaces, shell topology and pressure against CAD and analysis."""
import FreeCAD as App
import Part

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import random
import subprocess

import numpy as np
from circular_bend import integrated
from coaxial_program import compile_mass_oracle
from polygon_prism import interior

HERE=Path(__file__).resolve().parent


def generated():
    rng=random.Random(2026112520)
    outer=dict(id="outline",vertices=[dict(id=f"p{i}",xy_mm=list(x)) for i,x in enumerate(((-40,-16),(40,-16),(40,16),(-40,16)))])
    for i in range(128):
        bs=[]
        for j in range(i%4+1):
            r=rng.randint(1,3)
            bs.append(dict(id=f"drill{j}",center_mm=[-29+18*j+rng.randint(-1,1),rng.randint(-4,4)],
                           radius_mm=r,guard_mm=2*r+rng.randint(0,1)))
        ok=i%5!=0
        if not ok:
            if len(bs)>1:bs[1]["center_mm"]=bs[0]["center_mm"].copy()
            else:bs[0]["center_mm"]=[40,0]
        yield dict(id=f"seed-{i:03d}",expected=ok,recipe=dict(family="polygon_circular_bores_v1",
            outer=outer,bores=bs,height_mm=rng.randint(1,32),top_shift_mm=[rng.randint(-12,12),rng.randint(-12,12)]))


def reference(r):
    p=np.array([v["xy_mm"] for v in r["outer"]["vertices"]]);q=np.roll(p,-1,axis=0)
    cross=p[:,0]*q[:,1]-q[:,0]*p[:,1]
    signed=int(cross.sum());area=abs(signed)/2.
    centre=np.sum((p+q)*cross[:,None],axis=0)/(3.*signed)
    mom=area*centre
    for bore in r["bores"]:
        hole=math.pi*bore["radius_mm"]**2;area-=hole;mom-=hole*np.array(bore["center_mm"])
    centre=np.r_[mom/area+np.array(r["top_shift_mm"])/2.,r["height_mm"]/2.]
    return p,abs(signed),area,centre


def cad(p,r):
    vector=App.Vector(*r["top_shift_mm"],r["height_mm"]);verts=[App.Vector(*map(float,xy),0) for xy in p]
    solid=Part.Face(Part.makePolygon([*verts,verts[0]])).extrude(vector)
    for bore in r["bores"]:
        boundary=Part.makeCircle(bore["radius_mm"],App.Vector(*bore["center_mm"],0))
        solid=solid.cut(Part.Face(Part.Wire([boundary])).extrude(vector))
    return solid.removeSplitter()


def incidence(d):
    vertices=set();edges=set();faces=Counter()
    for b in d["blocks"]:
        shape=[v+1 for v in b["degrees"]]
        net=np.array(b["controls"],dtype=object).reshape(*shape,4)
        corners=[tuple(F(a,w) for a in pt[:3]) for ix in product(*[(0,n-1) for n in shape])
                 for pt in [net[ix]] for w in [pt[3]]]
        vertices.update(corners)
        for bit in (1,2,4):
            edges.update(tuple(sorted((corners[i],corners[i|bit]))) for i in range(8) if not i&bit)
            for flag in (False,True):faces[tuple(sorted(corners[i] for i in range(8) if bool(i&bit)==flag))]+=1
    return len(vertices)-len(edges)+len(faces)-len(d["blocks"]),max(faces.values())


def check_circles(d,r):
    blocks={b["id"]:b for b in d["blocks"]};bores={b["id"]:b for b in r["bores"]}
    tested=0;faces=Counter()
    for f in d["design_faces"]:
        if not f["feature"].startswith("bore/"):continue
        bore=bores[f["feature"][5:]]
        b=blocks[f["block"]];assert f["side"]=="v1" and b["degrees"]==[2,1,1]
        net=np.array(b["controls"],dtype=object).reshape(3,2,2,4)
        for z in (0,1):
            for t in (F(0),F(1,7),F(1,3),F(3,4),F(1)):
                coeff=[(1-t)**2,2*t*(1-t),t*t]
                H=[sum(coeff[i]*net[i,1,z,k] for i in range(3)) for k in range(4)]
                xy=[H[k]/H[3]-bore["center_mm"][k]-z*r["top_shift_mm"][k] for k in range(2)]
                assert sum(v*v for v in xy)==bore["radius_mm"]**2
                assert H[2]/H[3]==z*r["height_mm"]
                tested+=1
        faces[bore["id"]]+=1
    assert faces=={k:8 for k in bores}
    return tested


def classifications(pts,r,shapes):
    rng=np.random.default_rng(2026112520)
    lo,hi=pts.min(axis=0)-1,pts.max(axis=0)+1
    probes=[(*rng.uniform(lo,hi),rng.uniform(-.03,1.03)) for _ in range(256)]
    for b in r["bores"]:
        for theta,delta,t in product((0.,.25,1.8,3.3,5.1),(-.0001,.0001),(.1,.63)):
            xy=np.array(b["center_mm"])+(b["radius_mm"]+delta)*np.array([np.cos(theta),np.sin(theta)])
            probes.append((*xy,t))
    failures=[]
    for x,y,t in probes:
        expected=0<t<1 and interior(pts,x,y) and all(math.dist((x,y),b["center_mm"])>b["radius_mm"] for b in r["bores"])
        xy=np.array([x,y])+t*np.array(r["top_shift_mm"]);p=App.Vector(*xy,t*r["height_mm"])
        obs=[s.isInside(p,1e-8,False) for s in shapes]
        if any(v!=expected for v in obs):failures.append(dict(xyz=list(p),expected=bool(expected),actual=obs))
    return len(probes),failures


def wall_area(bore,r):
    # Split into smooth half-quadrants before Gaussian quadrature; also handles high shear.
    nodes,weights=np.polynomial.legendre.leggauss(64)
    result=0.
    for k in range(16):
        t=(nodes+1)*math.pi/16+k*math.pi/8
        projected=r["top_shift_mm"][0]*np.cos(t)+r["top_shift_mm"][1]*np.sin(t)
        result+=weights@np.hypot(r["height_mm"],projected)*math.pi/16
    return float(bore["radius_mm"]*result)


def run(cli,work):
    work.mkdir(parents=True,exist_ok=True);helper=compile_mass_oracle(work)
    def call(*cmd):
        p=subprocess.run([str(cli),*map(str,cmd)],capture_output=True,text=True)
        return p.returncode,json.loads(p.stdout)
    hashes={};rows=[];mp=work/"adopted.json";rp=work/"recipe.json"
    for item in [*json.loads((HERE/"circular_polygon_corpus.json").read_text()),*generated()]:
        name=item["id"];r=item["recipe"];rp.write_text(json.dumps(r))
        old=mp.read_bytes() if mp.exists() else None
        rc,report=call("build",rp,mp)
        row=dict(id=name,expected=item["expected"],audit=report,recipe=r);rows.append(row)
        ok=(rc==0)==item["expected"]
        if rc:
            row["old_model_preserved"]=(mp.read_bytes() if mp.exists() else None)==old;ok &= row["old_model_preserved"]
        elif item["expected"]:
            d=json.loads(mp.read_text());pts,twice,area,ctr=reference(r);height=r["height_mm"]
            digest=hashlib.sha256(json.dumps([d[k] for k in ("blocks","design_faces","topology")],sort_keys=True).encode()).hexdigest()
            hashes[name]=digest;chi,n=incidence(d);circ=check_circles(d,r)
            row.update(cell_characteristic=chi,max_face_multiplicity=n,circle_rational_checks=circ,
                       ideal_volume_mm3=area*height,ideal_centroid_mm=ctr.tolist(),partition_sha256=digest)
            ok &= chi==1-len(r["bores"]) and n==2 and report["projected_area_exact"]==dict(outer_area_twice_mm2=twice,pi_coefficient_mm2=-sum(b["radius_mm"]**2 for b in r["bores"]))
            if "equivalent_to" in item:ok &= digest==hashes[item["equivalent_to"]]
            elif not name.startswith("seed-"):
                sp=work/f"{name}.step";rc,step=call("export-step",mp,sp);row["step"]=step
                ok &= rc==0
                if not rc:
                    shapes=(cad(pts,r),Part.read(str(sp)))
                    valid=all(s.isValid() and len(s.Solids)==1 and len(s.Shells)==1 for s in shapes)
                    masses=[]
                    for i,s in enumerate(shapes):
                        bp=work/f"{name}-{i}.brep";s.exportBrep(str(bp))
                        masses.append(json.loads(subprocess.check_output([str(helper),str(bp),"--centroid"],text=True))[-1])
                    ntest,fail=classifications(pts,r,shapes);ve=max(abs(m["volume_mm3"]/(area*height)-1) for m in masses)
                    ce=max(float(np.linalg.norm(m["centroid_mm"]-ctr)) for m in masses)
                    row.update(valid=valid,occ_masses=masses,classification_count=ntest,classification_failures=fail[:8],
                               classification_mismatches=len(fail),cad_volume_error=ve,cad_centroid_error_mm=ce)
                    ok &= valid and not fail and ve<1e-8 and ce<1e-7
                native=integrated(d,12)
                expected_walls={f'bore/{b["id"]}':wall_area(b,r) for b in r["bores"]}
                errs={key:abs(native["areas_by_feature_mm2"][key]/expected-1) for key,expected in expected_walls.items()}
                row.update(native_q12_volume_mm3=native["volume_mm3"],q12_bore_area_relative=errs,analytic_wall_areas_mm2=expected_walls)
                ok &= abs(native["volume_mm3"]/(area*height)-1)<1e-9
                if name=="staggered":
                    bore=r["bores"][0];pressure=2.5
                    job=dict(schema="solidcraft_mechanical_job_v1",material=dict(young_mpa=100000,poisson=.25),
                             formulation="displacement",order=12,plan=dict(schema="adaptive_refinement_plan_v1",max_depth=0,max_cells=1000,steps=[]),
                             supports=[dict(select=dict(feature="start"),components=[0,1,2])],
                             loads=[dict(select=dict(feature="bore/"+bore["id"]),pressure_mpa=pressure)])
                    jp=work/"pressure.job.json";out=work/"pressure.result.json";jp.write_text(json.dumps(job));rc,diagnostic=call("solve-job",mp,jp,out)
                    row["pressure_job"]=dict(status=diagnostic,qualified=False)
                    if not rc:
                        result=json.loads(out.read_text());dx,dy=r["top_shift_mm"]
                        torque=pressure*math.pi*bore["radius_mm"]**2*np.array([-dy,dx,0.])
                        terr=np.linalg.norm(torque-result["applied_moment_n_mm"])/np.linalg.norm(torque)
                        ferr=np.linalg.norm(result["applied_resultant_n"])/(pressure*expected_walls["bore/"+bore["id"]])
                        row["pressure_job"].update(analytic_moment_n_mm=torque.tolist(),measured_moment_n_mm=result["applied_moment_n_mm"],
                            relative_moment_error=terr,normalized_force_error=ferr,strain_energy_n_mm=result["strain_energy_n_mm"],
                            qualified=bool(max(terr,ferr,result["relative_force_balance"],result["relative_moment_balance"])<1e-8 and result["strain_energy_n_mm"]>0))
                    ok &= row["pressure_job"]["qualified"]
        row["qualified"]=bool(ok)
        if not name.startswith("seed-") or len(rows)%32==0:print(name,row["qualified"],flush=True)
    return dict(schema="exact_circular_bores_in_polygon_prism_oracle_v1",occ=Part.OCC_VERSION,seed=2026112520,cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","work","output"):p.add_argument(k,type=Path)
    a=p.parse_args();result=run(a.cli.resolve(),a.work)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    if not all(r["qualified"] for r in result["cases"]):raise SystemExit(1)
