"""Check polygon prisms against closed-form moments, physical probes, and OCC."""
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

HERE = Path(__file__).resolve().parent


def random_cases():
    rng = random.Random(2026112517)
    rays = [(1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1)]
    for i in range(128):
        while True:
            x0,y0 = (rng.randint(-12,12) for _ in range(2))
            points = []
            for x,y in rays:
                r = rng.randint(2,14)
                points.append((x0+x*r,y0+y*r))
            determinants = []
            for j,b in enumerate(points):
                a,c=points[j-1],points[(j+1)%8]
                determinants.append((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))
            if all(determinants):
                break
        valid = i%4 != 3
        if not valid:
            points[4]=points[0]
        ids = [{"id":f"r{j}","xy_mm":p} for j,p in enumerate(points)]
        if i%2:
            ids.reverse()
        rotation = rng.randrange(8)
        yield {"id":f"seed-{i:03d}","admitted":valid,"recipe":{
            "family":"polygon_prism_v1", "vertices":ids[rotation:]+ids[:rotation],
            "height_mm":rng.randint(1,64),"top_shift_mm":[rng.randint(-32,32),rng.randint(-32,32)]}}


def measures(recipe):
    v = recipe["vertices"]
    points = np.array([p["xy_mm"] for p in v])
    x,y=points.T
    cross=x*np.roll(y,-1)-y*np.roll(x,-1)
    area2=int(cross.sum())
    center=(points+np.roll(points,-1,axis=0)).T @ cross/(3*area2)
    h=recipe["height_mm"]
    shift=np.array(recipe["top_shift_mm"])
    translation=np.array([*shift,h])
    areas={"start":abs(area2)/2,"end":abs(area2)/2}
    for j,a in enumerate(v):
        b=v[(j+1)%len(v)]
        edge=np.array([*(points[(j+1)%len(v)]-points[j]),0])
        name="edge/"+"/".join(sorted((a["id"],b["id"])))
        areas[name]=float(np.linalg.norm(np.cross(edge,translation)))
    return points,abs(area2),6*h*abs(area2)//2,np.array([*(center+shift/2),h/2]),areas


def decoded_measure(model):
    lookup={b["id"]:b for b in model["blocks"]}
    areas=Counter()
    for face in model["design_faces"]:
        block=lookup[face["block"]]
        xyz=np.array(block["controls"],dtype=float)
        xyz=(xyz[:,:3]/xyz[:,3:]).reshape(2,2,2,3)
        axis="uvw".index(face["side"][0])
        plane=np.take(xyz,int(face["side"][1]),axis=axis)
        a,b,c,d=plane[0,0],plane[1,0],plane[1,1],plane[0,1]
        areas[face["feature"]]+=0.5*(np.linalg.norm(np.cross(b-a,c-a))+np.linalg.norm(np.cross(c-a,d-a)))
    return areas


def interior(pts,x,y):
    inside=False
    for i,(a,b) in enumerate(pts):
        c,d=pts[(i+1)%len(pts)]
        if (b>y)!=(d>y) and x<(c-a)*(y-b)/(d-b)+a:
            inside=not inside
    return inside


def classify(pts,recipe,shapes):
    rng=random.Random(11125)
    h=recipe["height_mm"];shift=np.array(recipe["top_shift_mm"])
    lo,hi=pts.min(axis=0)-1,pts.max(axis=0)+1
    queries=[(rng.uniform(lo[0],hi[0]),rng.uniform(lo[1],hi[1]),rng.uniform(-.2,1.2)) for _ in range(384)]
    for i,a in enumerate(pts):
        b=pts[(i+1)%len(pts)]
        d=b-a;normal=np.array([-d[1],d[0]])/np.linalg.norm(d)
        for sign in (-1,1):
            x,y=(a+b)*.5+sign*0.0001*normal
            queries.extend([(x,y,f) for f in (.12,.5,.88)])
    mismatch=0
    for x,y,f in queries:
        ref=0<f<1 and interior(pts,x,y)
        point=App.Vector(x+f*shift[0],y+f*shift[1],h*f)
        mismatch+=int(any(s.isInside(point,1e-8,False)!=ref for s in shapes))
    return len(queries),mismatch


def run(cli,work):
    work.mkdir(parents=True,exist_ok=True)
    frozen=json.loads((HERE/"polygon_prism_corpus.json").read_text())
    model=work/"adopted.json";recipe_path=work/"recipe.json"
    hashes={};rows=[]
    for item in [*frozen,*random_cases()]:
        rid=item["id"];recipe=item["recipe"]
        recipe_path.write_text(json.dumps(recipe))
        before=model.read_bytes() if model.exists() else None
        proc=subprocess.run([cli,"build",recipe_path,model],capture_output=True,text=True)
        report=json.loads(proc.stdout)
        row={"id":rid,"expected":item["admitted"],"native":report,"recipe":recipe}
        ok=(proc.returncode==0)==item["admitted"]
        if proc.returncode:
            row["previous_unchanged"]=before==(model.read_bytes() if model.exists() else None)
            ok=ok and row["previous_unchanged"]
        elif item["admitted"]:
            data=json.loads(model.read_text())
            pts,area2,vol6,center,areas=measures(recipe)
            actual=decoded_measure(data)
            area_err=max((abs(actual[f]-a)/max(1.,a) for f,a in areas.items()),default=0)
            row.update(expected_area_twice_mm2=area2,expected_volume_sixths_mm3=vol6,feature_area_relative=area_err)
            ok=ok and len(data["blocks"])==3*(len(pts)-2) and area2==report["projected_area_twice_mm2"] and vol6==report["volume_sixths_mm3"] and areas.keys()==actual.keys() and area_err<1e-12
            digest=hashlib.sha256(json.dumps([data[k] for k in ("blocks","topology","design_faces")],sort_keys=True).encode()).hexdigest()
            hashes[rid]=digest
            if "equivalent_to" in item:
                row["same_partition"]=digest==hashes[item["equivalent_to"]]
                ok=ok and row["same_partition"]
            if not rid.startswith("seed-") and "equivalent_to" not in item:
                step=work/(rid+".step")
                exp=subprocess.run([cli,"export-step",model,step],capture_output=True,text=True)
                row["step_export"]=json.loads(exp.stdout)
                ok=ok and exp.returncode==0
                if exp.returncode==0:
                    vs=[App.Vector(float(p[0]),float(p[1]),0) for p in pts]
                    reference=Part.Face(Part.makePolygon([*vs,vs[0]])).extrude(App.Vector(*recipe["top_shift_mm"],recipe["height_mm"]))
                    shapes=[reference,Part.read(str(step))]
                    valid=all(s.isValid() and len(s.Solids)==1 and len(s.Shells)==1 for s in shapes)
                    volume_error=max(abs(s.Volume/(vol6/6)-1) for s in shapes)
                    ca_error=max(abs(s.Area/sum(areas.values())-1) for s in shapes)
                    centroid_error=max(np.linalg.norm(np.array([s.Solids[0].CenterOfMass.x,s.Solids[0].CenterOfMass.y,s.Solids[0].CenterOfMass.z])-center) for s in shapes)
                    count,miss=classify(pts,recipe,shapes)
                    row.update(valid_solids=valid,cad_volume_relative=volume_error,cad_area_relative=ca_error,
                               centroid_error_mm=centroid_error,classified_points=count,mismatches=miss)
                    ok=ok and valid and max(volume_error,ca_error)<1e-9 and centroid_error<1e-8 and miss==0
        row["qualified"]=bool(ok)
        rows.append(row)
        if not rid.startswith("seed-") or len(rows)%32==0:print(rid,ok,flush=True)
    return {"schema":"polygon_prism_mass_and_occ_v1","seed":2026112517,"occ":Part.OCC_VERSION,"cases":rows}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","work","output"):
        parser.add_argument(k,type=Path)
    a=parser.parse_args();result=run(a.cli.resolve(),a.work)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    if not all(r["qualified"] for r in result["cases"]):
        raise SystemExit(1)
