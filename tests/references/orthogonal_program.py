"""Verify integer box histories by an independent unit-voxel boundary and exact cubical incidences."""
import FreeCAD as App
import Part

import argparse
from collections import Counter,defaultdict
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import time

import numpy as np
from scipy.ndimage import label

HERE=Path(__file__).resolve().parent


def voxel(recipe):
    features=recipe["features"]
    ids=[f["id"] for f in features]
    if len(ids)!=len(set(ids)) or "base" in ids:return {"admit":False,"reason":"parameter"}
    boxes=[recipe["base_mm"],*(f["bounds_mm"] for f in features)]
    if any(a>=b for box in boxes for a,b in box):return {"admit":False,"reason":"parameter"}
    lo=np.array([min(b[j][0] for b in boxes)-1 for j in range(3)])
    hi=np.array([max(b[j][1] for b in boxes)+1 for j in range(3)])
    solid=np.zeros(tuple(hi-lo),bool)
    changed=np.zeros(solid.shape,dtype=int)
    for i,(bounds,operation) in enumerate(zip(boxes,["add",*(f["operation"] for f in features)])):
        region=tuple(slice(a-off,b-off) for (a,b),off in zip(bounds,lo))
        target=operation=="add"
        changed[region][solid[region]!=target]=i
        solid[region]=target
    volume=int(solid.sum())
    if volume==0:return {"admit":False,"reason":"empty"}
    # Enumerate the square boundary complex independently of the compressed block mesh.
    edge_uses=Counter();links=defaultdict(lambda:defaultdict(set));areas=Counter()
    for axis in range(3):
        tangents=[j for j in range(3) if j!=axis]
        for sign in (-1,1):
            outside=np.roll(~solid,-sign,axis=axis)
            for cell in np.argwhere(solid & outside):
                other=cell.copy();other[axis]+=sign
                owner=max(changed[tuple(cell)],changed[tuple(other)])
                feature=ids[owner-1] if owner else "base"
                base=cell+lo;base[axis]+=int(sign==1)
                areas[feature,axis,sign,int(base[axis])]+=1
                vertices=[]
                for x,y in ((0,0),(1,0),(1,1),(0,1)):
                    v=base.copy();v[tangents[0]]+=x;v[tangents[1]]+=y
                    vertices.append(tuple(map(int,v)))
                for j,v in enumerate(vertices):
                    prev=vertices[(j-1)%4];nxt=vertices[(j+1)%4]
                    edge_uses[tuple(sorted([v,nxt]))]+=1
                    links[v][prev].add(nxt);links[v][nxt].add(prev)
    if any(count!=2 for count in edge_uses.values()):return {"admit":False,"reason":"edge"}
    for link in links.values():
        if any(len(n)!=2 for n in link.values()):return {"admit":False,"reason":"vertex"}
        found=set();todo=[next(iter(link))]
        while todo:
            n=todo.pop()
            if n not in found:found.add(n);todo.extend(link[n]-found)
        if len(found)!=len(link):return {"admit":False,"reason":"vertex"}
    if label(solid)[1]!=1:return {"admit":False,"reason":"disconnected"}
    return {"admit":True,"volume":volume,"cavities":label(~solid)[1]-1,"areas":areas,"solid":solid,"origin":lo}


def native_measures(model):
    byid={b["id"]:b for b in model["blocks"]};areas=Counter();volume=0
    for b in byid.values():
        pts=np.array(b["controls"],dtype=np.int64);xyz=pts[:,:3]//pts[:,3:]
        assert b["degrees"]==[1,1,1] and np.all(xyz*pts[:,3:]==pts[:,:3])
        volume+=int(np.prod(np.ptp(xyz,axis=0)))
    for f in model["design_faces"]:
        b=byid[f["block"]];pts=np.array(b["controls"],dtype=np.int64);xyz=(pts[:,:3]//pts[:,3:]).reshape(2,2,2,3)
        axis="uvw".index(f["side"][0]);upper=int(f["side"][1])
        face=np.take(xyz,upper,axis=axis).reshape(4,3);spans=np.ptp(face,axis=0)
        normal=np.flatnonzero(spans==0);assert len(normal)==1
        j=int(normal[0]);assert j==axis
        areas[f["feature"],j,2*upper-1,int(face[0,j])]+=int(np.prod(np.delete(spans,j)))
    return volume,areas


def occ(recipe):
    def box(coords):
        return Part.makeBox(*(b-a for a,b in coords),App.Vector(*(a for a,b in coords)))
    shape=box(recipe["base_mm"])
    for f in recipe["features"]:
        tool=box(f["bounds_mm"])
        shape=shape.fuse(tool) if f["operation"]=="add" else shape.cut(tool)
    return shape.removeSplitter()


def random_programs():
    rng=random.Random(20261125)
    for i in range(256):
        features=[]
        for j in range(rng.randrange(1,9)):
            bounds=[]
            for _ in range(3):
                a,b=sorted(rng.sample(range(-2,11),2));bounds.append([a,b])
            features.append(dict(id=f"op{j}",operation=rng.choice(("add","cut")),bounds_mm=bounds))
        yield {"id":f"seed-{i:03}","recipe":{"family":"orthogonal_feature_program_v1","base_mm":[[0,8]]*3,"features":features}}


def run(cli,work):
    work.mkdir(parents=True,exist_ok=True)
    frozen=json.loads((HERE/"orthogonal_corpus.json").read_text());rows=[];random_rows=[]
    current=work/"current.json";rp=work/"current.recipe"
    for item in itertools.chain(frozen,random_programs()):
        start=time.monotonic();recipe=item["recipe"];random_case=item["id"].startswith("seed-")
        expected=voxel(recipe);rp.write_text(json.dumps(recipe))
        prior=current.read_bytes() if current.exists() else None
        proc=subprocess.run([str(cli),"build",str(rp),str(current)],text=True,capture_output=True)
        diag=json.loads(proc.stdout)
        row={"id":item["id"],"recipe":recipe,"expected":expected["admit"],"build":diag}
        if "admitted" in item and expected["admit"]!=item["admitted"]:raise ValueError(f"Contradictory frozen case {item['id']}")
        ok=(proc.returncode==0)==expected["admit"]
        if proc.returncode:
            row.update(previous_unchanged=prior==(current.read_bytes() if current.exists() else None),voxel_reason=expected.get("reason"))
            ok=ok and row["previous_unchanged"]
        elif expected["admit"]:
            volume,areas=native_measures(json.loads(current.read_text()))
            row.update(volume_mm3=volume,expected_volume_mm3=expected["volume"],cavities=int(expected["cavities"]),
                       area_mm2=sum(areas.values()),feature_supports=len(areas),exact_feature_areas=(areas==expected["areas"]))
            ok=ok and volume==expected["volume"]==diag["volume_mm3"] and areas==expected["areas"] and diag["cavities"]==expected["cavities"]
            if not random_case:
                step=work/(item["id"]+".step");export=subprocess.run([str(cli),"export-step",str(current),str(step)],capture_output=True,text=True)
                row["export"]=json.loads(export.stdout);ok=ok and export.returncode==0
                if export.returncode==0:
                    cad=occ(recipe);reloaded=Part.read(str(step));shapes=[cad,reloaded]
                    valid=all(s.isValid() and len(s.Solids)==1 and len(s.Shells)==1+expected["cavities"] for s in shapes)
                    err=max(abs(s.Volume/volume-1) for s in shapes)
                    areaerr=max(abs(s.Area/sum(areas.values())-1) for s in shapes)
                    # Sample midpoints on every compressed interval, including the adjacent exterior.
                    boxes=[recipe["base_mm"],*(f["bounds_mm"] for f in recipe["features"])]
                    axes=[]
                    for axis in range(3):
                        cuts=sorted({n for bounds in boxes for n in bounds[axis]})
                        axes.append([cuts[0]-.5,*[(a+b)/2 for a,b in zip(cuts,cuts[1:])],cuts[-1]+.5])
                    mismatches=0;count=0
                    for xyz in itertools.product(*axes):
                        index=tuple(np.floor(np.array(xyz)-expected["origin"]).astype(int))
                        ref=bool(expected["solid"][index]);point=App.Vector(*xyz)
                        mismatches+=int(any(s.isInside(point,1e-8,False)!=ref for s in shapes));count+=1
                    row.update(valid_solids=valid,cad_volume_relative=err,cad_area_relative=areaerr,classified_points=count,mismatches=mismatches)
                    ok=ok and valid and err<1e-9 and areaerr<1e-9 and mismatches==0
        row.update(qualified=bool(ok),elapsed_s=time.monotonic()-start)
        (random_rows if random_case else rows).append(row)
        if not random_case or len(random_rows)%32==0:print(item["id"],ok,flush=True)
    return {"schema":"orthogonal_unit_voxel_and_occ_v1","seed":20261125,"occ":Part.OCC_VERSION,"required":rows,"randomized":random_rows}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("cli","work","output"):parser.add_argument(key,type=Path)
    a=parser.parse_args();d=run(a.cli.resolve(),a.work)
    a.output.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    if not all(r["qualified"] for r in d["required"]+d["randomized"]):raise SystemExit(1)
