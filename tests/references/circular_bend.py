"""Independent ideal revolutions, physical classifications, and rational quadrature for G15."""
import FreeCAD as App
import Part
import numpy as np

import argparse
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import time

from coaxial_program import compile_mass_oracle

SECTIONS = ("reference_rectangle_v1", "core_disk_v1", "full_annulus_v1", "pierced_square_v1")


def section_properties(section):
    return {
        SECTIONS[0]: (400., 80.),
        SECTIONS[1]: (16*math.pi, 8*math.pi),
        SECTIONS[2]: (12*math.pi, 12*math.pi),
        SECTIONS[3]: (400.-16*math.pi, 80.+8*math.pi),
    }[section]


def inside_section(section, x, y):
    rad = math.hypot(x, y)
    if section == SECTIONS[0]:
        return max(abs(x), abs(y)) < 10
    if section == SECTIONS[1]:
        return rad < 4
    if section == SECTIONS[2]:
        return 2 < rad < 4
    return max(abs(x), abs(y)) < 10 and rad > 4


def reference_solid(section, radius, quadrants):
    axis, origin = App.Vector(0,0,1), App.Vector(0,0,0)
    def torus(r):
        return Part.makeTorus(radius,r,origin,axis,-180,180,90*quadrants)
    if section in (SECTIONS[1], SECTIONS[2]):
        body = torus(4)
        return body if section == SECTIONS[1] else body.cut(torus(2)).removeSplitter()
    corners = [App.Vector(radius+x,0,z) for x,z in ((-10,-10),(10,-10),(10,10),(-10,10),(-10,-10))]
    block = Part.Face(Part.makePolygon(corners)).revolve(origin,axis,90*quadrants)
    return block if section == SECTIONS[0] else block.cut(torus(4)).removeSplitter()


def basis(degree, t):
    if degree == 1:
        return np.column_stack((1-t,t)), np.column_stack((-np.ones_like(t),np.ones_like(t)))
    if degree != 2:
        raise ValueError("Unexpected degree")
    return np.column_stack(((1-t)**2,2*t*(1-t),t*t)), np.column_stack((-2*(1-t),2-4*t,2*t))


def geometry(block, nodes):
    degrees=block["degrees"]
    controls=np.asarray(block["controls"],dtype=float).reshape(*(d+1 for d in degrees),4)
    bz=[basis(d,x) for d,x in zip(degrees,nodes)]
    def evaluate(derivative):
        return np.einsum("ai,bj,ck,ijkd->abcd",*(b[int(derivative==i)] for i,b in enumerate(bz)),controls,optimize=True)
    h=evaluate(None)
    x=h[...,:3]/h[...,3:]
    derivatives=[evaluate(i) for i in range(3)]
    jac=np.stack([(dh[...,:3]-x*dh[...,3:])/h[...,3:] for dh in derivatives],axis=-1)
    return x,jac


def integrated(model,order):
    t,w=np.polynomial.legendre.leggauss(order);t=(1+t)/2;w=w/2
    blocks={b["id"]:b for b in model["blocks"]}
    volume=0.
    min_det=math.inf
    for b in blocks.values():
        _,j=geometry(b,[t,t,t])
        det=np.linalg.det(j)
        min_det=min(min_det,float(det.min()))
        volume+=float(np.einsum("i,j,k,ijk",w,w,w,det))
    areas={}
    for f in model["design_faces"]:
        axis="uvw".index(f["side"][0]); cap=int(f["side"][1])
        points=[t,t,t];points[axis]=np.array([float(cap)])
        _,j=geometry(blocks[f["block"]],points)
        a,b=(axis+1)%3,(axis+2)%3
        measure=np.linalg.norm(np.cross(j[...,a],j[...,b]),axis=-1).squeeze(axis=axis)
        areas[f["feature"]]=areas.get(f["feature"],0.)+float(np.einsum("i,j,ij",w,w,measure))
    return dict(volume_mm3=volume,area_mm2=sum(areas.values()),areas_by_feature_mm2=areas,sampled_det_min=min_det)


def classify(analytic, imported, section, radius, quadrants):
    rng=random.Random(202611250415)
    checked=0
    failures=[]
    for theta,xi,eta in itertools.chain(
        itertools.product((0.15,1.7,3.4,5.,6.1),(-9.,-3.,0.,3.,9.),(-9.,-2.5,0.,2.5,9.)),
        ((rng.uniform(0,2*math.pi),rng.uniform(-11,11),rng.uniform(-11,11)) for _ in range(384)),
    ):
        # Use the physical radius recovered from the Cartesian point. R+ξ
        # may be negative in test samples outside a small-radius disk bend.
        x=(radius+xi)*math.cos(theta);y=(radius+xi)*math.sin(theta)
        angle=math.atan2(y,x)%(2*math.pi);radial=math.hypot(x,y)-radius
        expected=angle<quadrants*math.pi/2 and inside_section(section,radial,eta)
        point=App.Vector(x,y,-eta)
        statuses=[s.isInside(point,1e-8,False) for s in (analytic,imported)]
        if any(v!=expected for v in statuses):
            failures.append(dict(xyz=[x,y,-eta],expected=expected,actual=statuses))
        checked+=1
    return checked,failures


def corpus():
    cases=[(s,r,q,True) for s in SECTIONS for r,q in ((16,1),(16,4),(24,2),(24,3))]
    for section,extent in ((SECTIONS[2],4),(SECTIONS[0],10)):
        cases.extend((section,r,q,True) for r,q in ((extent+1,4),(64,3)))
    for s in SECTIONS:
        extent=4 if s in (SECTIONS[1],SECTIONS[2]) else 10
        cases.extend((s,r,2,False) for r in (extent-1,extent))
    cases.extend((SECTIONS[0],r,q,False) for r,q in ((16,0),(16,5),(65,1)))
    return cases


def run(cli,work):
    work.mkdir(exist_ok=True,parents=True)
    helper=compile_mass_oracle(work)
    model=work/"current.json"
    def call(*args):
        p=subprocess.run([str(cli),*map(str,args)],capture_output=True,text=True)
        return p.returncode,json.loads(p.stdout)
    rows=[]
    for section,r,q,admit in corpus():
        start=time.monotonic()
        name=f"{section}-R{r}-Q{q}"
        recipe={"family":"circular_bend_v1","section":section,"major_radius_mm":r,"quadrants":q}
        path=work/(name+".recipe");path.write_text(json.dumps(recipe))
        old=model.read_bytes() if model.exists() else None
        code,status=call("build",path,model)
        row=dict(recipe=recipe,expected_admission=admit,build=status,qualified=False)
        if code:
            row["preserved"]=(model.read_bytes() if model.exists() else None)==old
            row["qualified"]=not admit and row["preserved"]
        elif admit:
            m=json.loads(model.read_text())
            A,L=section_properties(section)
            angle=q*math.pi/2
            volume=A*r*angle
            area=L*r*angle+(2*A if q<4 else 0)
            native={k:integrated(m,k) for k in (12,20)}
            ideal=reference_solid(section,r,q)
            row.update(target_volume_mm3=volume,target_area_mm2=area,native=native,
                       analytic_cad_volume_mm3=ideal.Volume,analytic_cad_area_mm2=ideal.Area)
            step=work/(name+".step")
            rc,export=call("export-step",model,step);row["step"]=export
            if rc==0:
                imported=Part.read(str(step))
                brep=work/(name+".brep");imported.exportBrep(str(brep))
                masses=json.loads(subprocess.check_output([str(helper),str(brep)],text=True))
                checked,bad=classify(ideal,imported,section,r,q)
                shells=2 if q==4 and section in (SECTIONS[2],SECTIONS[3]) else 1
                bound=status["determinant_lower_bound_mm3"]
                row.update(step_solids=len(imported.Solids),step_shells=len(imported.Shells),step_valid=imported.isValid(),
                           checked_points=checked,classifier_mismatches=bad,step_mass=masses)
                row["qualified"]=(ideal.isValid() and len(ideal.Solids)==1 and len(ideal.Shells)==shells
                    and abs(ideal.Volume/volume-1)<1e-9 and abs(ideal.Area/area-1)<1e-8
                    and imported.isValid() and len(imported.Solids)==1 and len(imported.Shells)==shells
                    and all(abs(t["volume_mm3"]/volume-1)<1e-8 and abs(t["surface_mm2"]/area-1)<1e-8 for t in masses)
                    and not bad and abs(native[20]["volume_mm3"]/volume-1)<1e-9
                    and abs(native[20]["area_mm2"]/area-1)<1e-8
                    and 0<bound<=native[20]["sampled_det_min"]*(1+1e-12))
        row["wall_s"]=time.monotonic()-start
        rows.append(row)
        print(name,row["qualified"],flush=True)
    return {"schema":"circular_bend_independent_revolution_v1","occ":Part.OCC_VERSION,
            "freecad":App.Version()[:3],"numpy":np.__version__,"cases":rows}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("cli","work","output"):parser.add_argument(key,type=Path)
    a=parser.parse_args()
    result=run(a.cli.resolve(),a.work.resolve())
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    if not all(row["qualified"] for row in result["cases"]):raise SystemExit(1)
