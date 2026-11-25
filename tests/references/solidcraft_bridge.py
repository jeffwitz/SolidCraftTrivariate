"""Replay native and SolidCraft bracket histories, including saved host documents."""
import FreeCAD as App
import Part

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess
import time

from coaxial_program import compile_mass_oracle
ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"upstream/fieldcraft/tests/oracle/filleted-bracket-fillet-stress/design.json"


def shape_member(x,y,z,rx,ry):
    material=(0<x<60 and 0<y<50 and 0<z<40 and
              (x<8 or y<8 or (x<14 and y<14 and (14-x)**2+(14-y)**2>36)))
    return material and (x-40)**2+(z-20)**2>rx**2 and (y-32)**2+(z-20)**2>ry**2


def samples(rx,ry):
    rng=random.Random(2026112501)
    for _ in range(2000):
        yield rng.uniform(-1,61),rng.uniform(-1,51),rng.uniform(-1,41)
    for theta in (0.17,1.42,2.9,3.4,4.78,5.9):
        for sign in (-1,1):
            delta=sign*.001
            for x in (0.25,4,7.75):
                yield 40+(rx+delta)*math.cos(theta),x,20+(rx+delta)*math.sin(theta)
                yield x,32+(ry+delta)*math.cos(theta),20+(ry+delta)*math.sin(theta)
    for theta in (3.3,3.7,4.,4.5):
        for delta in (-.001,.001):
            for z in (1.,20.,39.):
                yield 14+(6+delta)*math.cos(theta),14+(6+delta)*math.sin(theta),z


def command(*args):
    p=subprocess.run(list(map(str,args)),text=True,capture_output=True)
    if p.returncode:
        raise RuntimeError(f"{args[0]} {args[1]}: {p.stderr[-500:]} {p.stdout[-500:]}")
    return p.stdout


def bounds(shape):
    box=shape.optimalBoundingBox(False,False)
    return [box.XMin,box.YMin,box.ZMin,box.XMax,box.YMax,box.ZMax]


def run(cli,host,work):
    work.mkdir(parents=True,exist_ok=True)
    helper=compile_mass_oracle(work)
    source=json.loads(SOURCE.read_text())
    rows=[]
    for index,(rx,ry) in enumerate(((64,64),(40,80),(80,40),(50,74),(51,74))):
        case=work/f"x{rx}-y{ry}";case.mkdir(exist_ok=True)
        script=copy.deepcopy(source);commands=script["commands"]
        commands[5]["params"]["diameter"]=rx/8
        commands[6]["params"]["diameter"]=ry/8
        points=commands[1]["params"]["points"]
        rotation=index%6
        commands[1]["params"]["points"]=points[rotation:]+points[:rotation]
        if index%2:commands[5],commands[6]=commands[6],commands[5]
        inp=case/"script.json";inp.write_text(json.dumps(script))
        model=case/"native.json";step=case/"native.step";sc=case/"editable.sc";host_step=case/"host.step";reloaded=case/"saved.step"
        start=time.monotonic()
        command(host,"run",inp,"--save",sc,"--out",host_step,"--param-check","strict","--timings","--json-out",case/"host-timing.json")
        host_cost=time.monotonic()-start
        command(host,"run","--in",sc,"--out",reloaded,"--quiet")
        start=time.monotonic();report=json.loads(command(cli,"import-solidcraft",inp,model))
        native_cost=time.monotonic()-start
        command(cli,"export-step",model,step)
        shapes=[Part.read(str(path)) for path in (step,host_step,reloaded)]
        mass=[]
        for name,s in zip(("native","host","reopened"),shapes):
            brep=case/(name+".brep");s.exportBrep(str(brep))
            mass.append(json.loads(command(helper,brep)))
        expected=34080-math.pi*(360+8*((rx/16)**2+(ry/16)**2))
        volume=max(abs(m["volume_mm3"]/expected-1) for measures in mass for m in measures)
        area=max(abs(measures[-1]["surface_mm2"]/mass[1][-1]["surface_mm2"]-1) for measures in mass)
        boxes=list(map(bounds,shapes))
        box_error=max(abs(a-b) for box in boxes for a,b in zip(box,[0,0,0,60,50,40]))
        bad=[];checked=0
        for x,y,z in samples(rx/16,ry/16):
            v=App.Vector(x,y,z);inside=shape_member(x,y,z,rx/16,ry/16)
            actual=[shape.isInside(v,1e-8,False) for shape in shapes]
            if any(result!=inside for result in actual):
                bad.append(dict(xyz=[x,y,z],expected=inside,statuses=actual))
            checked+=1
        valid=[s.isValid() and len(s.Solids)==1 for s in shapes]
        row=dict(radii_ticks=[rx,ry],rotated_contour_vertices=rotation,reversed_independent_holes=bool(index%2),
                 source_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),native=report,host_timing=json.loads((case/"host-timing.json").read_text()),
                 host_build_save_step_wall_s=host_cost,native_translate_build_wall_s=native_cost,
                 saved_document_bytes=sc.stat().st_size,shape_order=["native","host","reopened"],
                 valid_single_solid=valid,faces=[len(s.Faces) for s in shapes],edges=[len(s.Edges) for s in shapes],
                 mass_properties=mass,analytic_volume_mm3=expected,relative_volume_error=volume,
                 maximum_relative_surface_difference=area,optimal_bounds_mm=boxes,bounds_error_mm=box_error,
                 classified_points=checked,mismatches=bad,
                 qualified=all(valid) and volume<=1e-8 and area<=1e-8 and box_error<=1e-6 and not bad and json.loads(model.read_text())["radii_sixteenths_mm"]==[rx,ry])
        rows.append(row)
        print(rx,ry,row["qualified"],volume,area,box_error,len(bad),flush=True)
    # The same previously accepted file is used for all rejected interpretations.
    prior=model.read_bytes()
    invalid={}
    def rejected(label,edit):
        s=copy.deepcopy(source);edit(s["commands"])
        path=work/f"invalid-{label}.json";path.write_text(json.dumps(s))
        result=subprocess.run([str(cli),"import-solidcraft",str(path),str(model)],capture_output=True,text=True)
        invalid[label]={"code":result.returncode,"diagnostic":json.loads(result.stdout),"previous_unchanged":prior==model.read_bytes()}
    rejected("missing_fillet",lambda c:c.pop(4))
    rejected("unknown_command",lambda c:c[4].update(command="solid.chamfer"))
    rejected("different_fillet",lambda c:c[4]["params"].update(edges=[[8,7,20]]))
    rejected("blind_hole",lambda c:c[5]["params"].update(depth=4))
    rejected("duplicate_hole",lambda c:c[6]["params"].update(position=[40,8,20]))
    rejected("off_lattice",lambda c:c[5]["params"].update(diameter=8.000000000000002))
    rejected("extrusion",lambda c:c[3]["params"].update(distance=39))
    rejected("wrong_outline",lambda c:c[1]["params"]["points"][2].__setitem__(1,9))
    return {"schema":"solidcraft_host_bracket_parity_v1","solidcraft_cli":command(host,"--version").strip(),
            "source_fixture":str(SOURCE.relative_to(ROOT)),"fixture_sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            "fieldcraft_revision":command("git","-C",ROOT/"upstream/fieldcraft","rev-parse","HEAD").strip(),
            "occ":Part.OCC_VERSION,"freecad":App.Version()[:3],"cases":rows,"refusals":invalid}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("cli",type=Path);p.add_argument("work",type=Path);p.add_argument("output",type=Path)
    p.add_argument("--host",default="solidcraft-cli")
    a=p.parse_args();result=run(a.cli.resolve(),a.host,a.work.resolve())
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    if not all(x["qualified"] for x in result["cases"]) or not all(x["code"]!=0 and x["previous_unchanged"] for x in result["refusals"].values()):
        raise SystemExit(1)
