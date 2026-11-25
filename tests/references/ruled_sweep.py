"""Check actual interchange maps and intersect their reloaded STEP shells in OCC."""
import FreeCAD as App
import Part
import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import time
from j1_geometry import evaluate, determinant, gauss

CASES=[(8,0,8,0,0),(0,8,12,4,-8),(-8,1,7,-16,24),(-8,-8,16,8,4),
       (16,0,1,-32,32),(4,-12,32,32,-32),(-16,0,8,0,0),(0,0,8,0,0)]


def invoke(cli,*args):
    result=subprocess.run([str(cli),*map(str,args)],text=True,capture_output=True)
    return result.returncode,result.stdout,result.stderr


def contour_area(wire, order):
    nodes, weights = gauss(order)
    integral = 0.
    for edge in wire.OrderedEdges:
        lo, hi = edge.ParameterRange
        knots = getattr(edge.Curve, "getKnots", lambda: [])()
        cuts = [lo, *(k for k in knots if lo < k < hi), hi]
        orientation = 1 if edge.Orientation == "Forward" else -1
        for a, b in zip(cuts, cuts[1:]):
            for t, weight in zip(nodes, weights):
                at = a + (b-a)*t
                p, d = edge.valueAt(at), edge.derivative1At(at)
                integral += orientation * weight * (b-a) * (p.x*d.y-p.y*d.x)/2
    return abs(integral)


def sections(solid,h,pierced,q,offset):
    output=[]
    for w in [0.,.25,.5,.75,1.]:
        z=w*h
        face_area=None
        area_q12=None
        if w in (0,1):
            faces=[face for face in solid.Faces if all(abs(v.Point.z-z)<1e-7 for v in face.Vertexes)]
            area=sum(face.Area for face in faces)
            wires=None
        else:
            wires=solid.slice(App.Vector(0,0,1),z)
            areas=sorted((Part.Face(wire).Area for wire in wires),reverse=True)
            if len(areas)!=1+int(pierced) or not all(wire.isClosed() for wire in wires):
                raise ValueError(f"Unexpected section wires at {w}: {len(areas)}")
            face_area=areas[0]-sum(areas[1:])
            def section_area(order):
                contour=sorted((contour_area(wire,order) for wire in wires),reverse=True)
                return contour[0]-sum(contour[1:])
            area_q12=section_area(12)
            area=section_area(20)
        scale=1.+(q-1.)*w
        expected=(400.-16.*math.pi*int(pierced))*abs(scale)**2
        # Three radial and square boundaries in both directions, away from
        # the OCC classification tolerance and away from either cap.
        mismatch=0
        if 0<w<1:
            for sx,sy in itertools.product([-10.1,-9.9,-4.01,-3.99,0.,3.99,4.01,9.9,10.1],repeat=2):
                source=complex(sx,sy)
                point=scale*source+w*offset
                ideal=abs(sx)<10 and abs(sy)<10 and (not pierced or abs(source)>4)
                actual=solid.isInside(App.Vector(point.real,point.imag,z),1e-8,False)
                mismatch+=int(ideal!=actual)
        output.append({"fraction":w,"area_mm2":area,"ideal_mm2":expected,"relative_error":abs(area/expected-1.),
                       "occ_face_area_mm2":face_area,"contour_q12_mm2":area_q12,
                       "wires":len(wires) if wires is not None else None,"classification_mismatches":mismatch})
    return output


def run(cli,work):
    work.mkdir(parents=True,exist_ok=True)
    rows=[]
    transactions=[]
    for section in ["reference_rectangle_v1","pierced_square_v1"]:
        pierced=section=="pierced_square_v1"
        area=400.-16.*math.pi*int(pierced)
        src=work/(section+".json");base_recipe=work/(section+".recipe")
        base_recipe.write_text(json.dumps({"family":section}))
        assert invoke(cli,"build",base_recipe,src)[0]==0
        base={b["id"]:b for b in json.loads(src.read_text())["blocks"]}
        for idx,(a,b,h,tx,ty) in enumerate(CASES):
            started=time.monotonic()
            prefix=work/f"{section}-{idx}"
            recipe={"family":"ruled_section_sweep_v1","section":section,"height_mm":h,
                    "top_multiplier_eighths":[a,b],"top_translation_eighths_mm":[tx,ty]}
            rp=prefix.with_suffix(".recipe");rp.write_text(json.dumps(recipe))
            model=prefix.with_suffix(".json")
            ret,out,error=invoke(cli,"build",rp,model)
            error=error or out
            row={"recipe":recipe,"build_status":ret}
            if idx>=6:
                row.update(qualified=ret!=0 and "ruled_sweep_collapsed_section" in error,diagnostic=error.strip())
                rows.append(row);continue
            if ret!=0:
                row.update(qualified=False,diagnostic=error.strip());rows.append(row);continue
            ret,s,e=invoke(cli,"audit",model)
            if ret:
                row.update(qualified=False,diagnostic=e.strip());rows.append(row);continue
            audit=json.loads(s)
            record=json.loads(model.read_text())
            q=complex(a,b)/8;shift=complex(tx,ty)/8
            ideal=area*h*(1+q.real+abs(q)**2)/3
            max_x=max_dx=max_inv=0.;min_jac=math.inf; volume=0.
            bound_verified=True; exact_rechecks=0
            nodes=list(zip(*gauss(8)))
            for block in record["blocks"]:
                reference=base[block["id"]]
                rng=random.Random(11871)
                for t in [(rng.random(),rng.random(),rng.random()) for _ in range(100)]:
                    x,j=evaluate(reference,t);y,k=evaluate(block,t)
                    source=complex(*x[:2]);factor=1+(q-1)*t[2]
                    xy=factor*source+t[2]*shift
                    max_x=max(max_x,abs(complex(*y[:2])-xy),abs(y[2]-t[2]*h))
                    inv=(complex(*y[:2])-y[2]*shift/h)/(1+(q-1)*y[2]/h)
                    max_inv=max(max_inv,abs(inv-source))
                    for axis in range(3):
                        d=factor*complex(*j[axis][:2])+((q-1)*source+shift)*(axis==2)
                        max_dx=max(max_dx,abs(complex(*k[axis][:2])-d),abs(k[axis][2]-h*(axis==2)))
                    det=determinant(k)
                    min_jac=min(min_jac,det)
                    if det<audit["determinant_lower_bound_mm3"]:
                        exact_rechecks+=1
                        _,exact=evaluate(block,tuple(map(Fraction,t)))
                        bound_verified &= determinant(exact)>=Fraction(audit["determinant_lower_bound_mm3"])
                for (u,wu),(v,wv),(w,ww) in itertools.product(nodes,repeat=3):
                    volume+=wu*wv*ww*determinant(evaluate(block,(u,v,w))[1])
            step=prefix.with_suffix(".step");ret,_,e=invoke(cli,"export-step",model,step)
            cad={"export_status":ret}
            if ret==0:
                try:
                    shape=Part.read(str(step));solids=shape.Solids
                    cad.update(solids=len(solids),valid=shape.isValid(),closed=all(s.isClosed() for s in shape.Shells),
                               volume_mm3=shape.Volume,volume_relative_error=abs(shape.Volume/ideal-1.))
                    cad["sections"]=sections(solids[0],h,pierced,q,shift)
                except Exception as exc:
                    cad["exception"]=repr(exc)
            else:
                cad["diagnostic"]=e.strip()
            row.update(audit=audit,integrated_q8_mm3=volume,ideal_volume_mm3=ideal,volume_relative_error=abs(volume/ideal-1.),
                       position_max_mm=max_x,inverse_max_mm=max_inv,jacobian_columns_max_mm=max_dx,
                       sampled_minimum_det_mm3=min_jac,exact_determinant_rechecks=exact_rechecks,cad=cad,elapsed_s=time.monotonic()-started)
            row["map_qualified"]=max(max_x,max_dx,max_inv)<=1e-11 and bound_verified and abs(volume/ideal-1.)<=1e-8
            row["cad_qualified"]=(cad.get("solids")==1 and cad.get("valid") and cad.get("closed")
                and cad.get("volume_relative_error",math.inf)<=1e-7 and len(cad.get("sections",[]))==5
                and all(v["relative_error"]<=1e-7 and v["classification_mismatches"]==0 for v in cad["sections"]))
            row["qualified"]=bool(row["map_qualified"] and row["cad_qualified"])
            rows.append(row)
        live=work/(section+"-adopted.json")
        assert invoke(cli,"build",work/f"{section}-2.recipe",live)[0]==0
        old=live.read_bytes()
        for idx in (6,7):
            ret,out,_=invoke(cli,"build",work/f"{section}-{idx}.recipe",live)
            assert ret and json.loads(out)["reason"]=="ruled_sweep_collapsed_section" and live.read_bytes()==old
        rejected=work/(section+"-tampered.json")
        step=work/(section+"-retained.step")
        assert invoke(cli,"export-step",live,step)[0]==0
        step_before=step.read_bytes()
        for field in ("top_multiplier_eighths","design_faces","blocks"):
            changed=json.loads(old)
            if field=="top_multiplier_eighths":
                changed[field][0]+=1
            elif field=="design_faces":
                changed[field].pop()
            else:
                changed[field][0]["controls"][0][0]+=1
            rejected.write_text(json.dumps(changed))
            ret,out,_=invoke(cli,"audit",rejected)
            assert ret and json.loads(out)["status"]=="invalid"
            assert invoke(cli,"export-step",rejected,step)[0] and step.read_bytes()==step_before
        transactions.append({"section":section,"pinches_preserving_model":2,
                             "tampered_models_refused_preserving_step":3})
    return {"transactions":transactions,"occ_version":Part.OCC_VERSION,"freecad_version":App.Version()[:3],"cases":rows}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("cli",type=Path);p.add_argument("work",type=Path);p.add_argument("output",type=Path)
    args=p.parse_args()
    report=run(args.cli.resolve(),args.work)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(sum(row["qualified"] for row in report["cases"]),"of",len(report["cases"]))
