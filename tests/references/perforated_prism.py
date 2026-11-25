"""Independent area, topology and CAD checks for multiply connected polygon extrusions."""
import FreeCAD as App
import Part

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import subprocess

import numpy as np

from polygon_prism import decoded_measure, interior, measures

HERE=Path(__file__).resolve().parent


def generated():
    rng=random.Random(2026112518)
    def rectangle(name,x,y,w,h):
        pts=((x,y),(x+w,y),(x+w,y+h),(x,y+h))
        return dict(id=name,vertices=[dict(id=name+str(i),xy_mm=xy) for i,xy in enumerate(pts)])
    for i in range(128):
        count=1+i%4
        holes=[]
        for j in range(count):
            x=-25+13*j+rng.randint(0,4)
            y=rng.randint(-9,3)
            holes.append(rectangle(f"h{j}",x,y,rng.randint(3,7),rng.randint(3,8)))
        ok=i%5!=4
        if not ok:
            if len(holes)>1:
                # Keep the second set of vertex names; make its region overlap the first.
                for v,w in zip(holes[1]["vertices"],holes[0]["vertices"]):
                    v["xy_mm"]=[w["xy_mm"][0]+1,w["xy_mm"][1]+1]
            else:
                for v in holes[0]["vertices"]:v["xy_mm"]=[v["xy_mm"][0]+60,v["xy_mm"][1]]
        yield dict(id=f"seed-{i:03d}",expected=ok,recipe=dict(family="polygonal_perforated_prism_v1",
                   outer=rectangle("outer",-30,-15,60,30),holes=holes,height_mm=rng.randint(1,32),
                   top_shift_mm=[rng.randint(-15,15),rng.randint(-15,15)]))


def exact_geometry(r):
    loops=[r["outer"],*r["holes"]]
    area2=0;vol6=0;moment=np.zeros(3);areas={};points=[]
    for index,loop in enumerate(loops):
        pts,a,v,c,faces=measures(dict(vertices=loop["vertices"],height_mm=r["height_mm"],top_shift_mm=r["top_shift_mm"]))
        sign=1 if index==0 else -1
        area2+=sign*a;vol6+=sign*v;moment+=sign*v*c;points.append(pts)
        for name,value in faces.items():
            if name in ("start","end"):
                areas[name]=areas.get(name,0.)+sign*value
            else:
                areas[name.replace("edge/","edge/"+loop["id"]+"/",1)]=value
    return points,area2,vol6,moment/vol6,areas


def topology(blocks):
    vertices=set();edges=set();faces=Counter()
    for b in blocks:
        h=b["controls"]
        # G17/G18 canonical storage has a common denominator 6.
        if any(p[3]!=6 for p in h):raise ValueError("Unexpected G18 Cartesian gauge")
        p=[tuple(v[:3]) for v in h]
        vertices.update(p)
        for i in range(8):
            for bit in (1,2,4):
                if not i&bit:edges.add(tuple(sorted((p[i],p[i|bit]))))
        for bit in (1,2,4):
            for value in (False,True):
                ids=tuple(sorted(p[i] for i in range(8) if bool(i&bit)==value))
                faces[ids]+=1
    outside=[f for f,n in faces.items() if n==1]
    # Each edge of a closed quadrilateral boundary is incident to exactly two boundary faces.
    bvertices={v for f in outside for v in f}
    bedges=set()
    for f in outside:
        for i in range(len(f)):
            for j in range(i):
                edge=tuple(sorted((f[i],f[j])))
                if edge in edges:bedges.add(edge)
    if max(faces.values())>2:raise ValueError("Nonmanifold cell face")
    return [len(vertices)-len(edges)+len(faces)-len(blocks),len(bvertices)-len(bedges)+len(outside)]


def occ_solid(pts,r):
    def body(loop):
        p=[App.Vector(float(x),float(y),0) for x,y in loop]
        return Part.Face(Part.makePolygon([*p,p[0]])).extrude(App.Vector(*r["top_shift_mm"],r["height_mm"]))
    shape=body(pts[0])
    for hole in pts[1:]:shape=shape.cut(body(hole))
    return shape.removeSplitter()


def probes(pts,r,shapes):
    lo,hi=pts[0].min(axis=0)-1,pts[0].max(axis=0)+1
    rng=random.Random(6065)
    samples=[(rng.uniform(lo[0],hi[0]),rng.uniform(lo[1],hi[1]),rng.uniform(-.05,1.05)) for _ in range(512)]
    for loop in pts:
        for i,a in enumerate(loop):
            b=loop[(i+1)%len(loop)];d=b-a;normal=np.array([-d[1],d[0]])/np.linalg.norm(d)
            for sign in (-1,1):
                x,y=(a+b)/2+sign*.0001*normal
                samples.extend((x,y,t) for t in (.1,.5,.9))
    fail=0;shift=np.array(r["top_shift_mm"])
    for x,y,t in samples:
        expected=(0<t<1 and interior(pts[0],x,y) and not any(interior(h,x,y) for h in pts[1:]))
        p=App.Vector(x+t*shift[0],y+t*shift[1],t*r["height_mm"])
        fail+=int(any(s.isInside(p,1e-8,False)!=expected for s in shapes))
    return len(samples),fail


def run(cli,work):
    work.mkdir(parents=True,exist_ok=True);hashes={};rows=[]
    rp=work/"recipe.json";mp=work/"adopted.json"
    for item in [*json.loads((HERE/"perforated_prism_corpus.json").read_text()),*generated()]:
        name=item["id"];r=item["recipe"];rp.write_text(json.dumps(r))
        prior=mp.read_bytes() if mp.exists() else None
        p=subprocess.run([cli,"build",rp,mp],capture_output=True,text=True)
        audit=json.loads(p.stdout);row={"id":name,"recipe":r,"expected":item["expected"],"audit":audit}
        ok=(p.returncode==0)==item["expected"]
        if p.returncode:
            row["previous_unchanged"]=prior==(mp.read_bytes() if mp.exists() else None)
            ok=ok and row["previous_unchanged"]
        elif item["expected"]:
            d=json.loads(mp.read_text());pts,area2,vol6,centroid,areas=exact_geometry(r)
            native=decoded_measure(d)
            err=max(abs(native[k]-v)/max(1.,v) for k,v in areas.items())
            n=sum(map(len,pts));h=len(r["holes"]);chi=topology(d["blocks"])
            row.update(vertices=n,holes=h,area_twice_mm2=area2,volume_sixths_mm3=vol6,
                       characteristic_solid_boundary=chi,max_named_area_relative=err)
            ok=ok and chi==[1-h,2-2*h] and len(d["blocks"])==3*(n+2*h-2) and area2==audit["projected_area_twice_mm2"] and vol6==audit["volume_sixths_mm3"] and native.keys()==areas.keys() and err<1e-12
            digest=hashlib.sha256(json.dumps([d[k] for k in ("blocks","topology","design_faces")],sort_keys=True).encode()).hexdigest()
            hashes[name]=digest
            if "equivalent_to" in item:
                row["same_partition"]=digest==hashes[item["equivalent_to"]];ok=ok and row["same_partition"]
            if not name.startswith("seed-") and "equivalent_to" not in item:
                step=work/(name+".step")
                p=subprocess.run([cli,"export-step",mp,step],capture_output=True,text=True)
                row["step"]=json.loads(p.stdout);ok=ok and p.returncode==0
                if not p.returncode:
                    shapes=[occ_solid(pts,r),Part.read(str(step))]
                    valid=all(s.isValid() and len(s.Solids)==1 and len(s.Shells)==1 for s in shapes)
                    ve=max(abs(s.Volume/(vol6/6)-1) for s in shapes)
                    ae=max(abs(s.Area/sum(areas.values())-1) for s in shapes)
                    ce=max(np.linalg.norm(np.array(list(s.Solids[0].CenterOfMass))-centroid) for s in shapes)
                    count,miss=probes(pts,r,shapes)
                    row.update(valid_solids=valid,cad_volume_relative=ve,cad_area_relative=ae,centroid_error_mm=ce,
                               classified_points=count,classification_mismatches=miss)
                    ok=ok and valid and ve<1e-9 and ae<1e-9 and ce<1e-8 and miss==0
        row["qualified"]=bool(ok);rows.append(row)
        if not name.startswith("seed-") or len(rows)%32==0:print(name,ok,flush=True)
    return dict(schema="perforated_polygon_prism_mass_topology_and_occ_v1",seed=2026112518,occ=Part.OCC_VERSION,cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","work","output"):p.add_argument(k,type=Path)
    a=p.parse_args();d=run(a.cli.resolve(),a.work)
    a.output.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    if not all(r["qualified"] for r in d["cases"]):raise SystemExit(1)
