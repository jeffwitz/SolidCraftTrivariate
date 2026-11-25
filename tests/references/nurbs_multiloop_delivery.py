#!/usr/bin/env python3
"""End-to-end M1 geometry, STEP, adaptive OCC mass, walls, and refusal oracle."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import random
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import nurbs_multiloop_oracle as source


def cases(base):
    through = copy.deepcopy(base); through["loops"] = [through["loops"][0]]
    pocket = copy.deepcopy(through); pocket["loops"][0]["z_interval_mm"] = [0, 8]
    mixed = copy.deepcopy(base)
    permuted = copy.deepcopy(base); permuted["loops"].reverse()
    permuted["loops"][0]["id"] = "alpha"; permuted["loops"][1]["id"] = "beta"
    touch = copy.deepcopy(through)
    for curve_name in ("lower", "upper"):
        touch["loops"][0][curve_name]["controls"][0][0] = 0
    overlap = copy.deepcopy(base)
    for curve_name in ("lower", "upper"):
        for pole in overlap["loops"][1][curve_name]["controls"]:
            pole[0] -= 35 * pole[2]
    degenerate = copy.deepcopy(pocket); degenerate["loops"][0]["z_interval_mm"] = [8, 8]
    stationary = json.loads((ROOT / "tests/counterexamples/nurbs_multiloop_stationary.json").read_text())
    return [
        {"id":"single-through", "accepted":True, "recipe":through},
        {"id":"single-blind-pocket", "accepted":True, "recipe":pocket},
        {"id":"two-loops-mixed-depth", "accepted":True, "recipe":mixed},
        {"id":"permuted-renamed", "accepted":True, "recipe":permuted, "same_geometry_as":"two-loops-mixed-depth"},
        {"id":"touching-stock-boundary", "accepted":False, "recipe":touch},
        {"id":"overlapping-x-projections", "accepted":False, "recipe":overlap},
        {"id":"degenerate-z-interval", "accepted":False, "recipe":degenerate},
        {"id":"stationary-profile", "accepted":False, "recipe":stationary},
    ]


def bernstein(n, t):
    return [math.comb(n, i) * t**i * (1-t)**(n-i) for i in range(n+1)]


def bernstein_derivative(n, t):
    if n == 0: return [0.0]
    low = bernstein(n-1, t)
    return [n * ((low[i-1] if i else 0.0) - (low[i] if i<n else 0.0)) for i in range(n+1)]


def evaluate_block(block, uvw):
    degree = block["degrees"]
    net = block["controls"]
    basis = [bernstein(n,t) for n,t in zip(degree,uvw)]
    dbasis = [bernstein_derivative(n,t) for n,t in zip(degree,uvw)]
    h = [0.0]*4; dh = [[0.0]*4 for _ in range(3)]
    q = 0
    for i in range(degree[0]+1):
        for j in range(degree[1]+1):
            for k in range(degree[2]+1):
                c = net[q]; q += 1
                factors = [basis[0][i],basis[1][j],basis[2][k]]
                w = factors[0]*factors[1]*factors[2]
                for m in range(4):
                    h[m] += w*c[m]
                    for ax,idx in enumerate((i,j,k)):
                        f = dbasis[ax][idx]
                        for other in range(3):
                            if other != ax: f *= factors[other]
                        dh[ax][m] += f*c[m]
    p = [h[m]/h[3] for m in range(3)]
    cols = [[(dh[a][m]*h[3]-h[m]*dh[a][3])/(h[3]**2) for m in range(3)] for a in range(3)]
    return p, cols


def determinant(cols):
    a,b,c=cols
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - b[0]*(a[1]*c[2]-a[2]*c[1]) + c[0]*(a[1]*b[2]-a[2]*b[1]))


def cross(a,b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def affine_uv_w_structure(block):
    """Check exact homogeneous structure for v- and w-affine extrusion axes."""
    du,dv,dw=block["degrees"]
    if dv != 1 or dw != 1: return False
    net=block["controls"]
    at=lambda i,j,k:net[(i*(dv+1)+j)*(dw+1)+k]
    for i in range(du+1):
        ref=at(i,0,0); w=ref[3]
        for j in range(dv+1):
            for k in range(dw+1):
                p=at(i,j,k)
                if p[3]!=w: return False
                if p[0]*w!=ref[0]*p[3]: return False
                if p[1]*at(i,j,0)[3]!=at(i,j,0)[1]*p[3]: return False
                if p[2]*at(0,0,k)[3]!=at(0,0,k)[2]*p[3]: return False
    return True


def source_state(recipe, xyz, tol=1e-8):
    x,y,z=xyz
    (xmin,xmax),(ymin,ymax)=recipe["bounds_xy_mm"]
    h=recipe["height_mm"]
    if x < xmin-tol or x > xmax+tol or y < ymin-tol or y > ymax+tol or z < -tol or z > h+tol:
        return "OUT"
    on_stock = min(abs(x-xmin),abs(x-xmax),abs(y-ymin),abs(y-ymax),abs(z),abs(z-h)) <= tol
    for loop in recipe["loops"]:
        lo,hi=loop["z_interval_mm"]
        spec=loop["lower"]; knots=spec["knots"]; p=spec["degree"]
        t0,t1=knots[p],knots[len(spec["controls"])]
        x0=source.curve(spec,t0)[0][0]; x1=source.curve(spec,t1)[0][0]
        if x < x0-tol or x > x1+tol: continue
        a,b=t0,t1
        for _ in range(70):
            m=(a+b)/2
            if source.curve(spec,m)[0][0] < x: a=m
            else: b=m
        t=(a+b)/2
        yl=source.curve(loop["lower"],t)[0][1]
        yu=source.curve(loop["upper"],t)[0][1]
        on_curve=min(abs(y-yl),abs(y-yu)) <= tol
        on_end=min(abs(x-x0),abs(x-x1)) <= tol and yl-tol <= y <= yu+tol
        in_channel=yl-tol <= y <= yu+tol
        if on_curve or on_end:
            if lo+tol < z < hi-tol: return "ON"
        if in_channel:
            # Open through cuts meet the stock end planes; these are openings,
            # while a blind-cut termination is an exposed internal face.
            if lo < z < hi:
                return "OUT"
            if lo == 0 and abs(z) <= tol or hi == h and abs(z-h) <= tol:
                return "OUT"
            if abs(z-lo) <= tol or abs(z-hi) <= tol:
                return "ON"
    return "ON" if on_stock else "IN"


def gauss_integrals(model, order):
    rule=[((x+1.0)*0.5,w*0.5) for x,w in source.gauss_legendre(order)]; volume=0.0; face_rows=[]
    blocks={b["id"]:b for b in model["blocks"]}
    for b in blocks.values():
        if affine_uv_w_structure(b):
            for u,wu in rule:
                _,cols=evaluate_block(b,(u,0.5,0.5));volume += wu*determinant(cols)
        else:
            for (u,wu),(v,wv),(w,ww) in __import__("itertools").product(rule,repeat=3):
                _,cols=evaluate_block(b,(u,v,w)); volume += wu*wv*ww*determinant(cols)
    for face in model["design_faces"]:
        b=blocks[face["block"]]; axis="uvw".index(face["side"][0]); side=int(face["side"][1])
        other=[a for a in range(3) if a!=axis]; area=0.0; points=[]
        for (u,wu),(v,wv) in __import__("itertools").product(rule,repeat=2):
            t=[0.0]*3;t[axis]=float(side);t[other[0]]=u;t[other[1]]=v
            p,cols=evaluate_block(b,t); n=cross(cols[other[0]],cols[other[1]])
            area += wu*wv*math.sqrt(sum(q*q for q in n))
        for u in (0.11,0.5,0.89):
            for v in (0.11,0.5,0.89):
                t=[0.0]*3;t[axis]=float(side);t[other[0]]=u;t[other[1]]=v
                points.append(evaluate_block(b,t)[0])
        face_rows.append({"block":face["block"],"side":face["side"],"feature":face["feature"],"area_mm2":area,"sample_points":points})
    return volume,face_rows


def verify_face_label(recipe, face):
    label=face["feature"]
    tol=2e-8
    max_error=0.0
    for p in face["sample_points"]:
        if source_state(recipe,p,tol)!="ON": return {"valid":False,"max_error_mm":math.inf}
        if label.startswith("stock:"):
            coordinate={"xmin":p[0]-recipe["bounds_xy_mm"][0][0],"xmax":p[0]-recipe["bounds_xy_mm"][0][1],
                        "ymin":p[1]-recipe["bounds_xy_mm"][1][0],"ymax":p[1]-recipe["bounds_xy_mm"][1][1],
                        "zmin":p[2],"zmax":p[2]-recipe["height_mm"]}[label.split(":")[1]]
            max_error=max(max_error,abs(coordinate))
        else:
            loop_id,surface=label.split(":",1)
            loop=next((x for x in recipe["loops"] if x["id"]==loop_id),None)
            if loop is None:return False
            if surface.startswith("wall_"):
                spec=loop["lower" if surface=="wall_lower" else "upper"]
                tlo=spec["knots"][spec["degree"]];thi=spec["knots"][len(spec["controls"])]
                a,b=tlo,thi
                for _ in range(60):
                    m=(a+b)/2
                    if source.curve(spec,m)[0][0]<p[0]:a=m
                    else:b=m
                expected=source.curve(spec,(a+b)/2)[0]
                max_error=max(max_error,math.dist((p[0],p[1]),expected))
                if p[2]<loop["z_interval_mm"][0]-tol or p[2]>loop["z_interval_mm"][1]+tol:return {"valid":False,"max_error_mm":max_error}
            elif surface.startswith("end_"):
                xs=[source.curve(loop["lower"],loop["lower"]["knots"][loop["lower"]["degree"]])[0][0],source.curve(loop["lower"],loop["lower"]["knots"][len(loop["lower"]["controls"])])[0][0]]
                target_x=xs[0] if surface=="end_lower" else xs[1]
                max_error=max(max_error,abs(p[0]-target_x))
                t=loop["lower"]["knots"][loop["lower"]["degree"]] if surface=="end_lower" else loop["lower"]["knots"][len(loop["lower"]["controls"])]
                lo_y=source.curve(loop["lower"],t)[0][1];hi_y=source.curve(loop["upper"],t)[0][1]
                if not min(lo_y,hi_y)-tol<=p[1]<=max(lo_y,hi_y)+tol:return {"valid":False,"max_error_mm":max_error}
            elif surface in ("floor","ceiling"):
                z=loop["z_interval_mm"][1] if surface=="ceiling" else loop["z_interval_mm"][0]
                max_error=max(max_error,abs(p[2]-z))
    return {"valid":True,"max_error_mm":max_error}


def run(cli, work, output):
    import os
    saved={k:globals().get(k) for k in ("json","Path","sys","subprocess","time","math","random","copy","argparse","source","HERE","hashlib")}
    import FreeCAD as App
    import Part
    globals().update(saved)
    base=json.loads((HERE/"nurbs_multiloop.recipe.json").read_text())
    work.mkdir(parents=True,exist_ok=True)
    sys.path.insert(0,str(HERE))
    import coaxial_program
    import pathlib
    coaxial_program.Path=pathlib.Path
    mass_helper=coaxial_program.compile_mass_oracle(work)
    target_path=work/"candidate.json"
    target_path.write_bytes(b"seed-for-atomic-refusal-check\n")
    rows=[]; reference_volumes={}
    def cli_call(*args):
        p=subprocess.run([str(cli),*map(str,args)],capture_output=True,text=True)
        try: detail=json.loads(p.stdout)
        except json.JSONDecodeError: detail={"stdout":p.stdout,"stderr":p.stderr}
        return p.returncode,detail
    for item in cases(base):
        start=time.monotonic(); recipe=item["recipe"]; prior=target_path.read_bytes()
        rp=work/f"{item['id']}.recipe.json";rp.write_text(json.dumps(recipe,separators=(",",":")))
        rc,diag=cli_call("build",rp,target_path)
        row={"id":item["id"],"expected_accepted":item["accepted"],"build_exit":rc,"build_diagnostic":diag}
        if not item["accepted"]:
            row["previous_model_preserved"]=(target_path.read_bytes()==prior)
            row["qualified"]=(rc!=0 and row["previous_model_preserved"])
            row["wall_s"]=time.monotonic()-start;rows.append(row);continue
        if rc:
            row.update(qualified=False,required_case_refused=True)
            row["wall_s"]=time.monotonic()-start;rows.append(row);continue
        model=json.loads(target_path.read_text())
        arc,audit=cli_call("audit",target_path)
        vol32,faces32=gauss_integrals(model,32)
        vol16,_=gauss_integrals(model,16)
        target=source.audit(recipe)["remaining_volume_mm3"]
        ref=source.audit(recipe)
        row.update(audit_exit=arc,audit=audit,blocks=len(model["blocks"]),design_faces=len(model["design_faces"]),
                   source_quadrature=ref,volume_q16_mm3=vol16,volume_q32_mm3=vol32,
                   volume_relative_error=abs(vol32/target-1),volume_quadrature_relative_delta=abs(vol32-vol16)/target,
                   wall_face_count=len(faces32),
                   area_by_feature={})
        face_checks=[verify_face_label(recipe,f) for f in faces32]
        row["wall_faces_valid"]=all(x["valid"] for x in face_checks)
        row["wall_max_error_mm"]=max(x["max_error_mm"] for x in face_checks)
        for f in faces32: row["area_by_feature"][f["feature"]]=row["area_by_feature"].get(f["feature"],0.0)+f["area_mm2"]
        step=work/f"{item['id']}.step";step_rc,step_diag=cli_call("export-step",target_path,step)
        row["step_export"]={"exit":step_rc,"diagnostic":step_diag}
        if step_rc==0:
            imported=Part.read(str(step));brep=work/f"{item['id']}.brep";imported.exportBrep(str(brep))
            adaptive=json.loads(subprocess.check_output([str(mass_helper),str(brep)],text=True))
            shell_compound=Part.makeCompound(imported.Shells)
            row.update(step_valid=imported.isValid(),step_solids=len(imported.Solids),step_shells=len(imported.Shells),
                       step_closed=imported.isClosed(),
                       step_faces=len(imported.Faces),step_default_volume_mm3=imported.Volume,step_adaptive_mass=adaptive,
                       step_relative_volume_error=abs(adaptive[-1]["volume_mm3"]/target-1))
            rng=random.Random(2026112508)
            points=[]
            for loop in recipe["loops"]:
                for q in (0.11,0.29,0.53,0.77,0.91):
                    lo=loop["lower"];t0=lo["knots"][lo["degree"]];t1=lo["knots"][len(lo["controls"])];t=t0+(t1-t0)*q
                    low=source.curve(loop["lower"],t)[0];high=source.curve(loop["upper"],t)[0];z0,z1=loop["z_interval_mm"]
                    for z in (z0+0.37*(z1-z0), (z0+z1)/2):
                        mid=(low[1]+high[1])/2
                        points.extend([((low[0]+high[0])/2,mid,z,"OUT"),
                                       (low[0],low[1],z,"ON"),(high[0],high[1],z,"ON"),
                                       ((low[0]+high[0])/2,low[1]-1,z,"IN"),((low[0]+high[0])/2,high[1]+1,z,"IN")])
                    if z1 < recipe["height_mm"]:
                        points.append(((low[0]+high[0])/2,(low[1]+high[1])/2,z1,"ON"))
            bounds=recipe["bounds_xy_mm"]
            for _ in range(300):
                p=(rng.uniform(*bounds[0]),rng.uniform(*bounds[1]),rng.uniform(0,recipe["height_mm"]))
                points.append((*p,source_state(recipe,p)))
            failures=[]
            for x,y,z,want in points:
                point=App.Vector(x,y,z);inside=imported.isInside(point,1e-8,False)
                if want=="IN" and not inside: failures.append([x,y,z,want,"OUT"])
                elif want=="OUT" and inside: failures.append([x,y,z,want,"IN"])
                elif want=="ON":
                    dist=shell_compound.distToShape(Part.Vertex(point))[0]
                    if dist>2e-7: failures.append([x,y,z,want,f"DIST={dist}"])
            row.update(classifications_checked=len(points),classification_mismatches=failures)
        row["qualified"]=(arc==0 and audit.get("status")=="valid" and row["volume_relative_error"]<1e-8 and
                          row["volume_quadrature_relative_delta"]<1e-8 and row["wall_faces_valid"] and step_rc==0 and
                          row.get("wall_max_error_mm",math.inf)<1e-10 and row.get("step_valid",False) and row.get("step_closed",False) and row.get("step_solids")==1 and row.get("step_relative_volume_error",math.inf)<1e-8 and
                          not row.get("classification_mismatches",[1]))
        if item.get("same_geometry_as"):
            row["invariance_reference_volume_mm3"]=reference_volumes.get(item["same_geometry_as"])
            row["invariance_relative_delta"]=abs(target-row["invariance_reference_volume_mm3"])/target if item["same_geometry_as"] in reference_volumes else None
            row["qualified"] &= row["invariance_relative_delta"] is not None and row["invariance_relative_delta"]<1e-12
        reference_volumes[item["id"]]=target
        row["wall_s"]=time.monotonic()-start;rows.append(row)
        print(item["id"],row["qualified"],flush=True)
    result={"oracle":"nurbs_multiloop_delivery_step_adaptive_mass_walls_classification_v1","cli":str(cli),"cli_sha256":hashlib.sha256(cli.read_bytes()).hexdigest(),"occ":Part.OCC_VERSION,"freecad":App.Version()[:3],"cases":rows}
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cli",type=Path);ap.add_argument("work",type=Path);ap.add_argument("output",type=Path)
    a=ap.parse_args();result=run(a.cli.resolve(),a.work.resolve(),a.output.resolve())
    if not all(c["qualified"] for c in result["cases"]): raise SystemExit(1)


if __name__=="__main__": main()
