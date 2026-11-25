"""Independent layer integrals, OCC Boolean history, and STEP replay for G14."""
import FreeCAD as App
import Part
import argparse
from collections import defaultdict
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import time

from j1_geometry import determinant,evaluate,gauss

HERE=Path(__file__).resolve().parent


def occupied(recipe,r,z):
    inside=z>0 and z<recipe["base_height_mm"]
    for feature in recipe["features"]:
        lo,hi=feature["z_mm"]
        if lo<z<hi and r<feature["radius_mm"]:
            inside=feature["operation"]=="add"
    return inside


def analytical(recipe):
    features=recipe["features"]
    zs=sorted({0,recipe["base_height_mm"],*(x for f in features for x in f["z_mm"])})
    radii=sorted({0,*(f["radius_mm"] for f in features)})
    rs=[*radii,15]
    area=[math.pi*(b*b-a*a) if b!=15 else 400-math.pi*a*a for a,b in zip(rs,rs[1:])]
    group=defaultdict(float)
    volume=0.
    for lo,hi in zip(zs,zs[1:]):
        middle=(hi+lo)/2
        for i,(a,b) in enumerate(zip(rs,rs[1:])):
            if not occupied(recipe,(a+b)/2,middle):continue
            volume+=area[i]*(hi-lo)
            for z,d in ((lo,-1),(hi,1)):
                if not occupied(recipe,(a+b)/2,z+d*.125):
                    group[f"z{z}:{d}"]+=area[i]
            for radius,sign in ((a,-1),(b,1)):
                if radius==0:continue
                if radius==15:
                    for axis in "xy":
                        for side in (-10,10):group[f"{axis}{side}"]+=20*(hi-lo)
                elif not occupied(recipe,radius+sign*.125,middle):
                    group[f"r{radius}:{sign}"]+=2*math.pi*radius*(hi-lo)
    return volume,dict(group)


def occ_reference(recipe):
    shape=Part.makeBox(20,20,recipe["base_height_mm"],App.Vector(-10,-10,0))
    for f in recipe["features"]:
        lo,hi=f["z_mm"]
        tool=Part.makeCylinder(f["radius_mm"],hi-lo,App.Vector(0,0,lo))
        shape=shape.cut(tool) if f["operation"]=="cut" else shape.fuse(tool)
    return shape.removeSplitter()


def material_history(recipe,point):
    x,y,z=point
    inside=abs(x)<10 and abs(y)<10 and 0<z<recipe["base_height_mm"]
    changed=-1
    for i,f in enumerate(recipe["features"]):
        if f["z_mm"][0]<z<f["z_mm"][1] and math.hypot(x,y)<f["radius_mm"]:
            next_inside=f["operation"]=="add"
            if next_inside!=inside:
                changed=i
                inside=next_inside
    return inside,changed


def native_integrals(model,recipe,order):
    qs=list(zip(*gauss(order)))
    zrule=list(zip(*gauss(2)))
    blocks={b["id"]:b for b in model["blocks"]}
    volume=0.
    for block in blocks.values():
        # Axial coordinates are affine and the complete section is constant in z.
        for (u,wu),(v,wv),(w,ww) in itertools.product(qs,qs,zrule):
            _,jac=evaluate(block,(u,v,w))
            volume+=wu*wv*ww*determinant(jac)
    group=defaultdict(float)
    provenance=defaultdict(lambda:defaultdict(float))
    bad_labels=[]
    for face in model["design_faces"]:
        b=blocks[face["block"]]
        ax="uvw".index(face["side"][0])
        sign=1 if face["side"][1]=="1" else -1
        d0,d1=(ax+1)%3,(ax+2)%3
        key=None
        checked_label=False
        for (u,wu),(v,wv) in itertools.product(qs,repeat=2):
            t=[0.]*3;t[ax]=(1+sign)/2;t[d0]=u;t[d1]=v
            x,j=evaluate(b,t);a,bb=j[d0],j[d1]
            normal=[sign*(a[(k+1)%3]*bb[(k+2)%3]-a[(k+2)%3]*bb[(k+1)%3]) for k in range(3)]
            if abs(normal[2])>math.hypot(*normal[:2]):
                current=f"z{round(x[2])}:{1 if normal[2]>0 else -1}"
            elif abs(abs(x[0])-10)<1e-9:current=f"x{round(x[0])}"
            elif abs(abs(x[1])-10)<1e-9:current=f"y{round(x[1])}"
            else:
                radius=math.hypot(*x[:2]);direction=normal[0]*x[0]+normal[1]*x[1]
                assert abs(radius-round(radius))<1e-10
                current=f"r{round(radius)}:{1 if direction>0 else -1}"
            assert key is None or key==current
            key=current
            surface=math.sqrt(sum(v*v for v in normal))*wu*wv
            if not checked_label:
                magnitude=math.sqrt(sum(v*v for v in normal))
                before,inside_index=material_history(recipe,[v-1e-6*n/magnitude for v,n in zip(x,normal)])
                after,outside_index=material_history(recipe,[v+1e-6*n/magnitude for v,n in zip(x,normal)])
                owner=max(inside_index,outside_index)
                expected=recipe["features"][owner]["id"] if owner>=0 else "base"
                if not before or after or face["feature"]!=expected:
                    bad_labels.append([face["block"],face["side"],face["feature"],expected,before,after])
                checked_label=True
            group[key]+=surface
            provenance[face["feature"]][key]+=surface
    return volume,dict(group),dict(provenance),bad_labels


def compile_mass_oracle(work):
    prefix=Path(sys.prefix)
    helper=work/"occ-volume"
    subprocess.run(["g++","-std=c++17","-O2","-I"+str(prefix/"include/opencascade"),
                    str(HERE/"occ_volume.cpp"),"-L"+str(prefix/"lib"),"-Wl,-rpath,"+str(prefix/"lib"),
                    "-lTKTopAlgo","-lTKBRep","-lTKG3d","-lTKMath","-lTKernel","-o",str(helper)],check=True)
    return helper


def replay(cli,work,corpus):
    work.mkdir(parents=True,exist_ok=True)
    mass_helper=compile_mass_oracle(work)
    def call(*args):
        p=subprocess.run([str(cli),*map(str,args)],capture_output=True,text=True)
        return p.returncode,json.loads(p.stdout)
    current=work/"current.json"
    rows=[]
    for item in corpus:
        start=time.monotonic()
        recipe=item["recipe"];name=item["id"]
        rp=work/f"{name}.recipe";rp.write_text(json.dumps(recipe))
        prior=current.read_bytes() if current.exists() else None
        rc,diag=call("build",rp,current)
        row=dict(id=name,expected_admission=item["admitted"],build_diagnostic=diag)
        if rc:
            row.update(previous_unchanged=(current.read_bytes() if current.exists() else None)==prior,
                       qualified=not item["admitted"])
            row["qualified"] &= row["previous_unchanged"]
        elif item["admitted"]:
            target,areas=analytical(recipe)
            model=json.loads(current.read_text())
            history={}
            for order in (8,16,24):
                volume,actual,provenance,bad_labels=native_integrals(model,recipe,order)
                history[order]={"volume_mm3":volume,"support_areas_mm2":actual}
            stable=(abs(volume-history[16]["volume_mm3"])<=1e-9*target
                    and all(abs(actual[k]-v)<=1e-8*v for k,v in history[16]["support_areas_mm2"].items()))
            area_error=max([abs(actual[k]/v-1) if k in actual else 1. for k,v in areas.items()]+[int(k not in areas) for k in actual])
            cad=occ_reference(recipe)
            ok_geom=(stable and not bad_labels and abs(volume/target-1)<=1e-9 and area_error<=1e-8 and abs(cad.Volume/target-1)<1e-9 and cad.isValid() and len(cad.Solids)==1)
            assert call("audit",current)[1]["status"]=="valid"
            step=work/f"{name}.step";code,export=call("export-step",current,step)
            row.update(ideal_volume_mm3=target,integrated_volume_mm3=volume,relative_volume_error=abs(volume/target-1),
                       reference_area_by_support=areas,native_area_by_support=actual,native_area_by_feature=provenance,provenance_mismatches=bad_labels,quadrature_by_order=history,quadrature_stable=stable,
                       maximum_relative_area_error=area_error,occ_reference_volume_mm3=cad.Volume,
                       occ_reference_shells=len(cad.Shells),native_qualified=ok_geom,step_diagnostic=export,step_qualified=False)
            if code==0:
                imported=Part.read(str(step))
                mismatch=0;checked=0
                rng=random.Random(2026112503)
                points=[(rng.uniform(-11,11),rng.uniform(-11,11),rng.uniform(-1,25)) for _ in range(400)]
                radii=sorted({0,10,*(f["radius_mm"] for f in recipe["features"])})
                heights=sorted({-1,0,25,recipe["base_height_mm"],*(z for f in recipe["features"] for z in f["z_mm"])})
                for lo,hi in zip(heights,heights[1:]):
                    for a,b in zip(radii,radii[1:]):
                        for angle in (0.31,1.74,3.52,5.68):
                            r=(a+b)/2
                            points.append((r*math.cos(angle),r*math.sin(angle),(hi+lo)/2))
                for x,y,z in points:
                    ideal=abs(x)<10 and abs(y)<10 and occupied(recipe,math.hypot(x,y),z)
                    point=App.Vector(x,y,z)
                    mismatch+=int(cad.isInside(point,1e-8,False)!=ideal or imported.isInside(point,1e-8,False)!=ideal)
                    checked+=1
                brep=work/f"{name}.brep"
                imported.exportBrep(str(brep))
                adaptive_mass=json.loads(subprocess.check_output([str(mass_helper),str(brep)],text=True))
                err=max(abs(row["volume_mm3"]/target-1) for row in adaptive_mass)
                row.update(step_default_volume_mm3=imported.Volume,step_default_relative_volume_error=abs(imported.Volume/target-1),
                           step_adaptive_mass=adaptive_mass)
                row.update(step_solids=len(imported.Solids),step_shells=len(imported.Shells),step_valid=imported.isValid(),
                           step_relative_volume_error=err,points_classified=checked,classification_mismatches=mismatch)
                row["step_qualified"]=(imported.isValid() and len(imported.Solids)==1 and len(imported.Shells)==len(cad.Shells) and err<1e-8 and mismatch==0)
            row["qualified"]=bool(ok_geom and row["step_qualified"])
        else:
            row["qualified"]=False
        row["wall_s"]=time.monotonic()-start;rows.append(row)
        print(name,row["qualified"],flush=True)
    return {"occ":Part.OCC_VERSION,"freecad":App.Version()[:3],"cases":rows}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("cli",type=Path);p.add_argument("work",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args()
    result=replay(a.cli.resolve(),a.work,json.loads((HERE/"coaxial_corpus.json").read_text()))
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    if not all(row["qualified"] for row in result["cases"]):
        raise SystemExit(1)
