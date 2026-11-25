#!/usr/bin/env python3
"""Run the cross-cut plate through SolidCore's real B-rep sweeps and booleans."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import random
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOLIDCORE = ROOT / "upstream/solidcraft/solidcore"
DEFAULT_RECIPE = HERE / "g2_cross_sketches.recipe.json"


def deboor(degree, knots, controls, t, side):
    n=len(controls)-1
    spans=[i for i in range(degree,n+1) if knots[i]<knots[i+1] and
           (knots[i]<=t<knots[i+1] if side=="right" else knots[i]<t<=knots[i+1])]
    if not spans and side=="right" and t==knots[n+1]: spans=[n]
    if not spans: raise ValueError(f"no span for t={t} side={side}")
    k=spans[0] if side=="right" else spans[-1]
    work=[list(p) for p in controls[k-degree:k+1]]
    for r in range(1,degree+1):
        for j in range(degree,r-1,-1):
            i=k-degree+j; den=knots[i+degree-r+1]-knots[i]
            a=Fraction(0) if den==0 else (t-knots[i])/den
            work[j]=[(1-a)*work[j-1][c]+a*work[j][c] for c in range(len(work[j]))]
    return work[degree]


def derivative_spline(degree,knots,controls):
    out=[]
    for i in range(len(controls)-1):
        den=knots[i+degree+1]-knots[i+1]
        out.append([Fraction(0) if den==0 else degree*(controls[i+1][c]-controls[i][c])/den
                    for c in range(len(controls[i]))])
    return degree-1,knots[1:-1],out


def homogeneous_jets(curve,t,side):
    p=int(curve["degree"]);u=[Fraction(x) for x in curve["knots"]]
    c=[[Fraction(x) for x in row] for row in curve["controls"]]
    h=deboor(p,u,c,Fraction(t),side)
    p1,u1,c1=derivative_spline(p,u,c);dh=deboor(p1,u1,c1,Fraction(t),side)
    p2,u2,c2=derivative_spline(p1,u1,c1);ddh=deboor(p2,u2,c2,Fraction(t),side)
    return h,dh,ddh


def rational_jet(curve,t,side="right"):
    h,dh,ddh=homogeneous_jets(curve,t,side);w=h[-1];dw=dh[-1];ddw=ddh[-1]
    p=[h[0]/w,h[1]/w]
    d=[(dh[i]-p[i]*dw)/w for i in range(2)]
    dd=[(ddh[i]-2*d[i]*dw-p[i]*ddw)/w for i in range(2)]
    return [float(x) for x in p],[float(x) for x in d],[float(x) for x in dd]


def rational_value_fast(curve,t):
    """Floating source-curve evaluator used only to bracket oracle roots."""
    p=int(curve["degree"]);u=curve["knots"];c=curve["controls"];n=len(c)-1
    k=next((i for i in range(p,n+1) if u[i] <= t < u[i+1]), n)
    work=[list(map(float,row)) for row in c[k-p:k+1]]
    for r in range(1,p+1):
        for j in range(p,r-1,-1):
            i=k-p+j;den=u[i+p-r+1]-u[i]
            a=0.0 if den==0 else (t-u[i])/den
            work[j]=[(1-a)*work[j-1][d]+a*work[j][d] for d in range(len(work[j]))]
    h=work[p]
    return h[0]/h[2],h[1]/h[2]


def bezier_spans(curve, frame):
    p=int(curve["degree"]);u=[Fraction(x) for x in curve["knots"]]
    c=[[Fraction(x) for x in row] for row in curve["controls"]]
    n=len(c)-1;domains=[(u[i],u[i+1]) for i in range(p,n+1) if u[i]<u[i+1]]
    p1,u1,c1=derivative_spline(p,u,c)
    origin=list(map(float,frame["origin_mm"]));axis_u=list(map(float,frame["u"]));axis_v=list(map(float,frame["v"]))
    result=[]
    for a,b in domains:
        h0=deboor(p,u,c,a,"right");h1=deboor(p,u,c,b,"left")
        d0=deboor(p1,u1,c1,a,"right");d1=deboor(p1,u1,c1,b,"left")
        dt=b-a
        homogeneous=[h0,[h0[i]+dt*d0[i]/3 for i in range(len(h0))],
                     [h1[i]-dt*d1[i]/3 for i in range(len(h1))],h1]
        points=[]
        for hx,hy,w in homogeneous:
            xyz=[w*origin[k]+hx*axis_u[k]+hy*axis_v[k] for k in range(3)]
            points.append([float(xyz[k]/w) for k in range(3)]+[float(w)])
        result.append({"degree":3,"knots":[0,0,0,0,1,1,1,1],"points":points,
                       "source_interval":[int(a),int(b)]})
    return result


def gauss_legendre(n):
    pairs=[]
    for i in range(1,n+1):
        z=math.cos(math.pi*(i-.25)/(n+.5))
        for _ in range(100):
            p0,p1=1.0,z
            for j in range(2,n+1): p0,p1=p1,((2*j-1)*z*p1-(j-1)*p0)/j
            pn=p1 if n>1 else z;pnm1=p0 if n>1 else 1.0
            dp=n*(z*pn-pnm1)/(z*z-1);dz=pn/dp;z-=dz
            if abs(dz)<2e-16: break
        pairs.append((z,2/((1-z*z)*dp*dp)))
    return pairs


def curve_xy_jets(curve,t,side="right"):
    p,d,dd=rational_jet(curve,t,side)
    return p,d


def profile_area(curve,order):
    degree=int(curve["degree"]);u=curve["knots"];c=curve["controls"]
    total=0.0;rule=gauss_legendre(order)
    for i in range(degree,len(c)):
        a,b=u[i],u[i+1]
        if a>=b: continue
        for x,w in rule:
            t=(a+b)/2+(b-a)*x/2
            (px,py),(dx,dy)=curve_xy_jets(curve,t)
            total += w*(b-a)/2*.5*(px*dy-py*dx)
    return total


def horizontal_roots(curve,y):
    p=int(curve["degree"]);u=curve["knots"];c=curve["controls"]
    roots=[]
    def f(t): return rational_value_fast(curve,t)[1]-y
    for i in range(p,len(c)):
        a,b=float(u[i]),float(u[i+1])
        if a>=b: continue
        cuts=[a+(b-a)*j/16 for j in range(17)]
        vals=[f(t) for t in cuts]
        for left,right,fl,fr in zip(cuts,cuts[1:],vals,vals[1:]):
            if fl==0: roots.append(left)
            if fl*fr<0:
                lo,hi=left,right
                for _ in range(70):
                    mid=(lo+hi)/2
                    if f(lo)*f(mid)<=0: hi=mid
                    else: lo=mid
                roots.append((lo+hi)/2)
    unique=[]
    for t in sorted(roots):
        if not unique or abs(t-unique[-1])>1e-8: unique.append(t)
    return unique


def strip_width(curve,y):
    roots=horizontal_roots(curve,y)
    if len(roots)%2: raise ValueError(f"odd source-profile crossings at y={y}: {roots}")
    xs=sorted(rational_value_fast(curve,t)[0] for t in roots)
    return sum(xs[i+1]-xs[i] for i in range(0,len(xs),2))


def source_measure(recipe,order):
    bounds=recipe["bounds_xyz_mm"];sx=bounds[0][1]-bounds[0][0];sy=bounds[1][1]-bounds[1][0];sz=bounds[2][1]-bounds[2][0]
    free,cross=recipe["sketches"]
    profile=free["loops"][0]["curve"]
    rect=cross["loops"][0]["curve"]
    y0=rect["controls"][0][0]/rect["controls"][0][2];y1=rect["controls"][1][0]/rect["controls"][1][2]
    z0=rect["controls"][0][1]/rect["controls"][0][2];z1=rect["controls"][2][1]/rect["controls"][2][2]
    area=profile_area(profile,order); rule=gauss_legendre(order)
    common=(y1-y0)*0.0
    common=sum(w*(y1-y0)/2*strip_width(profile,(y0+y1)/2+(y1-y0)*x/2) for x,w in rule)*(z1-z0)
    vp=area*sz;vr=sx*(y1-y0)*(z1-z0);vu=vp+vr-common;stock=sx*sy*sz
    return {"order":order,"profile_area_mm2":area,"profile_tool_volume_mm3":vp,
            "rectangle_tool_volume_mm3":vr,"tool_intersection_volume_mm3":common,
            "tool_union_volume_mm3":vu,"remaining_solid_volume_mm3":stock-vu,
            "stock_volume_mm3":stock,"rectangle_y_mm":[y0,y1],"rectangle_z_mm":[z0,z1]}


def source_class(recipe,p,tol=1e-8):
    x,y,z=p;(xb,yb,zb)=recipe["bounds_xyz_mm"]
    if not(xb[0]-tol<=x<=xb[1]+tol and yb[0]-tol<=y<=yb[1]+tol and zb[0]-tol<=z<=zb[1]+tol): return "Outside"
    curve=recipe["sketches"][0]["loops"][0]["curve"]
    rect=recipe["sketches"][1]["loops"][0]["curve"]
    y0,y1=[rect["controls"][i][0]/rect["controls"][i][2] for i in (0,1)]
    z0,z1=[rect["controls"][i][1]/rect["controls"][i][2] for i in (0,2)]
    # The source boundary itself is a B-spline, not only the roots of horizontal
    # slices. Check its distance explicitly so exact knot/crossing points remain ON.
    degree=int(curve["degree"]);knots=curve["knots"];controls=curve["controls"]
    closest=float("inf")
    for i in range(degree,len(controls)):
        a,b=float(knots[i]),float(knots[i+1])
        if not a<b:continue
        vals=[a+(b-a)*j/24 for j in range(25)]
        ds=[sum((q-r)**2 for q,r in zip(rational_value_fast(curve,t),(x,y))) for t in vals]
        j=min(range(len(ds)),key=ds.__getitem__)
        if ds[j]<closest:closest=ds[j]
        if 0<j<len(vals)-1:
            lo,hi=vals[j-1],vals[j+1]
            for _ in range(32):
                m1=lo+(hi-lo)*0.3819660112501051;m2=lo+(hi-lo)*0.6180339887498949
                d1=sum((q-r)**2 for q,r in zip(rational_value_fast(curve,m1),(x,y)))
                d2=sum((q-r)**2 for q,r in zip(rational_value_fast(curve,m2),(x,y)))
                if d1<d2:hi=m2
                else:lo=m1
            closest=min(closest,d1,d2)
    on_profile=closest<=tol*tol*16
    # Find a ray crossing count along the source profile at this ordinate.
    roots=horizontal_roots(curve,y)
    xs=sorted(curve_xy_jets(curve,t)[0][0] for t in roots)
    in_profile=any(xs[i]-tol<=x<=xs[i+1]+tol for i in range(0,len(xs)-1,2))
    in_rect=y0-tol<=y<=y1+tol and z0-tol<=z<=z1+tol
    in_profile_strict=any(xs[i]+tol<x<xs[i+1]-tol for i in range(0,len(xs)-1,2))
    in_rect_strict=y0+tol<y<y1-tol and z0+tol<z<z1-tol
    if in_profile_strict and not on_profile and zb[0]+tol<z<zb[1]-tol:return "Outside"
    if in_rect_strict and xb[0]+tol<x<xb[1]-tol:return "Outside"
    if on_profile:return "OnBoundary"
    if in_profile and (min(abs(x-a) for a in xs)<=tol or min(abs(y-b) for b in [y0,y1])<=tol or min(abs(z-b) for b in zb)<=tol): return "OnBoundary"
    if in_rect and (min(abs(y-b) for b in [y0,y1])<=tol or min(abs(z-b) for b in [z0,z1])<=tol): return "OnBoundary"
    if min(abs(x-b) for b in xb+[0])<=tol or min(abs(y-b) for b in yb)<=tol or min(abs(z-b) for b in zb)<=tol: return "OnBoundary"
    return "Outside" if in_profile or in_rect else "Inside"


def classification_points(recipe):
    profile=recipe["sketches"][0]["loops"][0]["curve"]
    rect=recipe["sketches"][1]["loops"][0]["curve"]
    y0,y1=[rect["controls"][i][0]/rect["controls"][i][2] for i in (0,1)]
    z0,z1=[rect["controls"][i][1]/rect["controls"][i][2] for i in (0,2)]
    pts=[[10,10,2],[50,25,2],[75,25,2],[10,25,6],[75,25,6],[-1,10,2],
         [50,y0,6],[50,y1,6],[50,25,z0],[50,25,z1]]
    knots=profile["knots"];p=profile["degree"]
    for t in range(8):
        xy=curve_xy_jets(profile,t,"right")[0]
        pts.append([xy[0],xy[1],2.0]);pts.append([xy[0],xy[1],10.0])
        for side,tt in (("left",t-1e-5),("right",t+1e-5)):
            if tt < profile["knots"][profile["degree"]] or tt > profile["knots"][len(profile["controls"])]:
                continue
            q=curve_xy_jets(profile,tt,"right" if side=="right" else "left")[0]
            pts.extend([[q[0],q[1],2.0],[q[0],q[1],10.0]])
    # The four transverse crossings and all eight points where the rectangle's
    # horizontal edges meet the curved source wall. Probe both material sides and
    # the four nearby z layers around each horizontal seam.
    eps=1e-5
    for t in (2.5,3.5,6.5,7.5):
        (x,y),(dx,dy)=rational_jet(profile,t,"right")[:2]
        n=math.hypot(dx,dy);nx,ny=-dy/n,dx/n
        for z in (4-eps,4+eps,6.0,8-eps,8+eps):
            pts.extend([[x+nx*s*eps,y+ny*s*eps,z] for s in (-2,0,2)])
    rng=random.Random(2026112517)
    for _ in range(180):
        pts.append([rng.uniform(-1,101),rng.uniform(-1,51),rng.uniform(-1,13)])
    return [{"point_mm":p,"expected":source_class(recipe,p)} for p in pts]


def payload_for(recipe):
    sketches=recipe["sketches"];profile_sketch=sketches[0];rect_sketch=sketches[1]
    arcs=bezier_spans(profile_sketch["loops"][0]["curve"],profile_sketch["frame"])
    rect_curve=rect_sketch["loops"][0]["curve"]
    frame=rect_sketch["frame"];origin=list(map(float,frame["origin_mm"]));u=list(map(float,frame["u"]));v=list(map(float,frame["v"]))
    rect_points=[]
    for wx,wy,w in rect_curve["controls"]:
        rect_points.append([origin[k]+wx/w*u[k]+wy/w*v[k] for k in range(3)])
    # Degree one and equal weights make each source span exactly its endpoint segment.
    if rect_curve["degree"]!=1 or len(rect_curve["controls"])!=5 or rect_points[0]!=rect_points[-1]:
        raise ValueError("SolidCore line conversion requires the fixture's exact degree-one closed rectangle")
    payload={"bounds_xyz_mm":recipe["bounds_xyz_mm"],"profile_arcs":arcs,"rect_points_world_mm":rect_points,
             "free_direction":profile_sketch["direction"],"free_length_mm":profile_sketch["length_mm"],
             "cross_direction":rect_sketch["direction"],"cross_length_mm":rect_sketch["length_mm"],
             "classification_points":[x["point_mm"] for x in classification_points(recipe)]}
    return payload


def run(recipe_path,work,output):
    recipe=json.loads(recipe_path.read_text());work.mkdir(parents=True,exist_ok=True)
    payload=payload_for(recipe);payload_path=work/"solidcore-input.json";payload_path.write_text(json.dumps(payload,indent=2)+"\n")
    package=work/"solidcore-harness";(package/"src").mkdir(parents=True,exist_ok=True)
    (package/"src/main.rs").write_text((HERE/"g2_cross_brep_solidcore.rs").read_text())
    def fmt(x):
        value=format(float(x),".17g")
        return value if any(ch in value for ch in ".eE") else value+".0"
    arcs=payload["profile_arcs"]
    if any(p[3] <= 0.0 for arc in arcs for p in arc["points"]):
        raise ValueError("rational Bezier span has a nonpositive weight")
    arc_source="["+",".join("["+",".join("["+",".join(fmt(v) for v in p)+"]" for p in arc["points"])+"]" for arc in arcs)+"]"
    rect_source="["+",".join("["+",".join(fmt(v) for v in p)+"]" for p in payload["rect_points_world_mm"])+"]"
    point_source="["+",".join("["+",".join(fmt(v) for v in p)+"]" for p in payload["classification_points"])+"]"
    bounds_source="["+",".join("["+",".join(fmt(v) for v in axis)+"]" for axis in payload["bounds_xyz_mm"])+"]"
    input_source=(f"const PROFILE_ARCS: [[[f64;4];4];{len(arcs)}]={arc_source};\n"+
        f"const RECT_POINTS: [[f64;3];{len(payload['rect_points_world_mm'])}]={rect_source};\n"+
        f"const CLASS_POINTS: [[f64;3];{len(payload['classification_points'])}]={point_source};\n"+
        f"const BOUNDS: [[f64;2];3]={bounds_source};\n"+
        f"const FREE_DIRECTION: [f64;3]=["+",".join(fmt(v) for v in payload['free_direction'])+"];\n"+
        f"const CROSS_DIRECTION: [f64;3]=["+",".join(fmt(v) for v in payload['cross_direction'])+"];\n"+
        f"const FREE_LENGTH: f64={fmt(payload['free_length_mm'])};\n"+
        f"const CROSS_LENGTH: f64={fmt(payload['cross_length_mm'])};\n")
    profile=recipe["sketches"][0]["loops"][0]["curve"]
    jet_targets=[]
    for t in range(8):
        for side,tt,code in (("left",(8-1e-5 if t==0 else t-1e-5),-1.0),("right",t+1e-5,1.0)):
            if tt<0 or tt>8: continue
            x,y=rational_jet(profile,tt,side)[0]
            for z in (2.0,10.0):jet_targets.append([x,y,z,float(t),code])
    jt="["+",".join("["+",".join(fmt(v) for v in row)+"]" for row in jet_targets)+"]"
    input_source+=f"const JET_TARGETS: [[f64;5];{len(jet_targets)}]={jt};\n"
    seam_probes=[];eps=1e-5
    for t in (2.5,3.5,6.5,7.5):
        (x,y),(dx,dy)=rational_jet(profile,t,"right")[:2];norm=math.hypot(dx,dy);nx,ny=-dy/norm,dx/norm
        for z in (4+eps,6.0,8-eps):
            for sign in (-1.0,1.0):
                seam_probes.append([t,z,sign,x+sign*eps*nx,y+sign*eps*ny,z])
    st="["+",".join("["+",".join(fmt(v) for v in row)+"]" for row in seam_probes)+"]"
    input_source+=f"const SEAM_PROBES: [[f64;6];{len(seam_probes)}]={st};\n"
    tol_probes=[]
    for t in (2.5,3.5,6.5,7.5):
        (x,y),(dx,dy)=rational_jet(profile,t,"right")[:2];norm=math.hypot(dx,dy);nx,ny=-dy/norm,dx/norm
        for z in (4+eps,6.0,8-eps):
            for sign in (-1.0,1.0):
                for delta in (1e-5,1e-4,1e-3,1e-2):
                    tol_probes.append([t,z,sign,delta,x+sign*delta*nx,y+sign*delta*ny,z])
    tp="["+",".join("["+",".join(fmt(v) for v in row)+"]" for row in tol_probes)+"]"
    input_source+=f"const SEAM_TOLERANCE_PROBES: [[f64;7];{len(tol_probes)}]={tp};\n"
    (package/"src/input.rs").write_text(input_source)
    patch_path=ROOT/"build/solidcore-deps/patch.toml"
    patch_text=patch_path.read_text() if patch_path.exists() else ""
    manifest=f'''[package]\nname = "g2-cross-solidcore-harness"\nversion = "0.1.0"\nedition = "2024"\n\n[dependencies]\nsolidcore = {{ path = "{SOLIDCORE}/crates/solidcore" }}\nsolidcore-ops = {{ path = "{SOLIDCORE}/crates/solidcore-ops" }}\nsolidcore-step = {{ path = "{SOLIDCORE}/crates/solidcore-step" }}\n'''+patch_text
    (package/"Cargo.toml").write_text(manifest)
    report=work/"solidcore.json";step=work/"solidcore.step";target=work/"target"
    env=dict(__import__("os").environ);env["CARGO_TARGET_DIR"]=str(package/"target");env["RUST_MIN_STACK"]="8388608"
    cmd=["rtk","cargo","run"]
    if patch_text:cmd.append("--offline")
    cmd.extend(["--manifest-path",str(package/"Cargo.toml"),"--",str(step)])
    proc=subprocess.run(cmd,cwd=ROOT,env=env,capture_output=True,text=True)
    backend={}
    for line in proc.stdout.splitlines():
        if "=" in line:
            key,value=line.split("=",1)
            if key=="classifications": backend[key]=value.split(",")
            elif key=="boundary_jets":
                records=[]
                for raw in filter(None,value.split(";")):
                    f=raw.split(",")
                    nums=[float(x) for x in f[3:]]
                    records.append({"knot_parameter":float(f[0]),"side":int(float(f[1])),"face":f[2],"distance_mm":nums[0],
                        "normal":nums[1:4],"gaussian_curvature_per_mm2":nums[4],"principal_curvatures_per_mm":[nums[5],nums[6]],
                        "principal_curvature_vectors":[nums[7:10],nums[10:13]],"surface_du":nums[13:16],"surface_dv":nums[16:19],"z_mm":nums[19] if len(nums)>19 else None})
                backend[key]=records
            elif key=="seam_probes":
                backend[key]=[{"t":float(f[0]),"z_mm":float(f[1]),"side":int(float(f[2])),"native":f[3],"nearest_brep_face_distance_mm":float(f[4])}
                              for f in (raw.split(",") for raw in filter(None,value.split(";")))]
            elif key=="seam_tolerance_probes":
                backend[key]=[{"t":float(f[0]),"z_mm":float(f[1]),"side":int(float(f[2])),"offset_mm":float(f[3]),"native":f[4]}
                              for f in (raw.split(",") for raw in filter(None,value.split(";")))]
            else:
                try: backend[key]=int(value)
                except ValueError:
                    try: backend[key]=float(value)
                    except ValueError: backend[key]=value
    run_data={"command":cmd,"returncode":proc.returncode,"stdout":proc.stdout,"stderr":proc.stderr,
              "solidcore_git_revision":subprocess.check_output(["rtk","git","-C",str(ROOT/"upstream/solidcraft"),"rev-parse","HEAD"],text=True).strip(),
              "dependency_patch":str(patch_path) if patch_text else None,
              "dependencies_without_patch_note":None if patch_text else "Use standard Cargo dependencies pinned by upstream/solidcraft/solidcore/Cargo.lock; run `cargo fetch --locked --manifest-path <harness>/Cargo.toml` with registry access, then rerun this script.",
              "source_recipe":str(recipe_path),"backend_result":backend if proc.returncode==0 else None,
              "step_path":str(step) if step.exists() else None,
              "step_sha256":hashlib.sha256(step.read_bytes()).hexdigest() if step.exists() else None}
    try:
        oracle16=source_measure(recipe,16);oracle32=source_measure(recipe,32);oracle64=source_measure(recipe,64)
        run_data["independent_source_oracle"]={"q16":oracle16,"q32":oracle32,"q64":oracle64,
            "q32_q64_relative_delta":abs(oracle64["remaining_solid_volume_mm3"]-oracle32["remaining_solid_volume_mm3"])/oracle64["remaining_solid_volume_mm3"],
            "classifications":classification_points(recipe),
            "profile_knot_jets":[{"t":t,"left":rational_jet(recipe["sketches"][0]["loops"][0]["curve"],8 if t==0 else t,"left"),
                                   "right":rational_jet(recipe["sketches"][0]["loops"][0]["curve"],t,"right")} for t in range(8)]}
        run_data["source_integration_relative_change"]={k:abs(oracle64[k]-oracle32[k])/max(1,abs(oracle64[k])) for k in oracle64 if k.endswith("volume_mm3")}
        native=backend.get("classifications",[])
        expected=[x["expected"] for x in run_data["independent_source_oracle"]["classifications"]]
        run_data["native_classification_comparison"]={"point_count":len(expected),"native_count":len(native),
            "mismatch_count":sum(a!=b for a,b in zip(expected,native)),
            "mismatches":[{"index":i,"point_mm":run_data["independent_source_oracle"]["classifications"][i]["point_mm"],"expected":a,"native":b}
                          for i,(a,b) in enumerate(zip(expected,native)) if a!=b],
            "boundary_tolerance_mm":2e-7}
        jets=backend.get("boundary_jets",[]);profile=recipe["sketches"][0]["loops"][0]["curve"];jet_comparison=[]
        for rec in jets:
            t=rec["knot_parameter"]
            tt=(8-1e-5 if t==0 and rec["side"]<0 else t+1e-5 if rec["side"]>0 else t-1e-5)
            _,d,dd=rational_jet(profile,tt,"right" if rec["side"]>0 else "left")
            speed=math.hypot(*d);tan=[d[0]/speed,d[1]/speed,0.0]
            tangents=(rec["surface_du"],rec["surface_dv"])
            tangent=sorted(tangents,key=lambda q:math.hypot(q[0],q[1]),reverse=True)[0]
            tnorm=math.sqrt(sum(q*q for q in tangent));tangent=[q/tnorm for q in tangent]
            align=abs(sum(tan[i]*tangent[i] for i in range(3)))
            dot=(d[0]*dd[0]+d[1]*dd[1])/(speed*speed)
            kvec=[(dd[i]-tan[i]*(d[0]*dd[0]+d[1]*dd[1])/speed)/ (speed*speed) for i in range(2)]+[0.0]
            idx=max(range(2),key=lambda i:abs(rec["principal_curvatures_per_mm"][i]))
            # On an extrusion the nonzero principal direction is along the section
            # tangent; signed section curvature acts along the surface normal.
            kv=[rec["principal_curvatures_per_mm"][idx]*q for q in rec["normal"]]
            kv_error=min(math.dist(kvec,kv),math.dist(kvec,[-q for q in kv]))
            jet_comparison.append({**rec,"source_tangent":tan,"tangent_alignment_abs_dot":align,
                "tangent_alignment_signed_dot":sum(tan[i]*tangent[i] for i in range(3)),
                "source_curvature_vector_per_mm":kvec,"brep_curvature_vector_error_per_mm":kv_error,
                "source_curvature_magnitude_per_mm":math.sqrt(sum(q*q for q in kvec)),"brep_section_curvature_per_mm":rec["principal_curvatures_per_mm"][idx]})
        run_data["native_wall_jet_comparison"]={"method":"nearest B-rep face surface, bounded Newton projection in face UV domain; surface evaluate through second derivatives and EvaluateGeometric; sample curve parameters knot +/- 1e-5 mm-param on Z=2,10 walls","sample_count":len(jet_comparison),"maximum_distance_mm":max((x["distance_mm"] for x in jet_comparison),default=None),
            "minimum_tangent_alignment_abs_dot":min((x["tangent_alignment_abs_dot"] for x in jet_comparison),default=None),
            "minimum_tangent_alignment_signed_dot":min((x["tangent_alignment_signed_dot"] for x in jet_comparison),default=None),
            "maximum_curvature_vector_error_per_mm":max((x["brep_curvature_vector_error_per_mm"] for x in jet_comparison),default=None),"samples":jet_comparison}
        seam_records=[]
        for rec in backend.get("seam_probes",[]):
            (x,y),(dx,dy)=rational_jet(profile,rec["t"],"right")[:2];norm=math.hypot(dx,dy)
            px=x+rec["side"]*1e-5*(-dy/norm);py=y+rec["side"]*1e-5*(dx/norm)
            point=[px,py,rec["z_mm"]];seam_records.append({**rec,"point_mm":point,"source_expected":source_class(recipe,point)})
        run_data["native_seam_probe_audit"]={"source_normal_offset_mm":1e-5,"classifier_tolerance_mm":2e-7,"probes":seam_records,
            "mismatch_count":sum(p["source_expected"]!=p["native"] for p in seam_records),
            "on_boundary_max_distance_mm":max((p["nearest_brep_face_distance_mm"] for p in seam_records if p["native"]=="OnBoundary"),default=None)}
        tolerance_records=[]
        for rec in backend.get("seam_tolerance_probes",[]):
            (x,y),(dx,dy)=rational_jet(profile,rec["t"],"right")[:2];norm=math.hypot(dx,dy)
            point=[x+rec["side"]*rec["offset_mm"]*(-dy/norm),y+rec["side"]*rec["offset_mm"]*(dx/norm),rec["z_mm"]]
            tolerance_records.append({**rec,"point_mm":point,"source_expected":source_class(recipe,point)})
        run_data["native_seam_tolerance_sweep"]={"passed_classifier_tolerance_mm":2e-7,"offsets_mm":[1e-5,1e-4,1e-3,1e-2],
            "first_offset_where_all_crossing_probes_match_source":next((delta for delta in (1e-5,1e-4,1e-3,1e-2)
                if all(p["source_expected"]==p["native"] for p in tolerance_records if p["offset_mm"]==delta)),None),"probes":tolerance_records}
    except Exception as exc:
        run_data["source_oracle_error"]=f"{type(exc).__name__}: {exc}"
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(run_data,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print(json.dumps({"status":"solidcore_ok" if proc.returncode==0 else "solidcore_failed","output":str(output),"step":run_data["step_path"]}))
    return proc.returncode


def main():
    ap=argparse.ArgumentParser();ap.add_argument("recipe",type=Path,default=DEFAULT_RECIPE,nargs="?");ap.add_argument("--work",type=Path,default=ROOT/"build/g2-cross-brep");ap.add_argument("--output",type=Path,default=ROOT/"docs/results/g2-cross-brep-solidcore-20261125.json")
    a=ap.parse_args();raise SystemExit(run(a.recipe,a.work,a.output))


if __name__=="__main__": main()
