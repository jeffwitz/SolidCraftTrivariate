"""Independent geometric oracle for the fixed filleted bracket template."""

from fractions import Fraction as F
from itertools import product
import json

from .arc_strip import certify_arc_strip
from .incidence import enumerate_faces
from .trace import extract_trace


def _regions():
    """Partition the manufacturing dimensions without invoking the builder."""
    for z0,z1 in ((0,10),(10,30),(30,40)):
        suffix = f"_z{z0}_{z1}"
        for side in ("east","north"):
            yield "elbow_"+side+suffix, ("elbow",side,z0,z1)
        for axis,ticks in (("x",(14,16,30,50,60)),("y",(14,16,22,42,50))):
            for i,(a,b) in enumerate(zip(ticks,ticks[1:])):
                if i==2 and z0==10:
                    for turn,side in enumerate(("east","north","west","south")):
                        yield f"hole_{axis}_sector_{side}"+suffix, ("hole",axis,turn)
                else:
                    prefix = "connector" if i==0 else "arm"
                    bounds = ((a,b),(0,8)) if axis=="x" else ((0,8),(a,b))
                    yield f"{prefix}_{axis}_{a}_{b}"+suffix, ("box",*bounds,(z0,z1))


def _require(ok, reason):
    if not ok:
        raise ValueError(reason)


def _straight_box(block, bounds):
    _require(block["degrees"]==[1,1,1], "box_degrees")
    controls = block["controls"]
    weight = controls[0][3]
    for pt,corner in zip(controls,product(*bounds),strict=True):
        _require(pt==[weight*c for c in corner]+[weight], "box_partition_coordinates")
    return (bounds[0][1]-bounds[0][0])*(bounds[1][1]-bounds[1][0])*(bounds[2][1]-bounds[2][0])


def _section(block, region):
    """Undo the specified rigid placement and verify the axial product."""
    _require(block["degrees"]==[2,1,1], "curved_degrees")
    hole = region[0]=="hole"
    if hole:
        _,axis,turn = region
        height = 8
    else:
        _,side,z0,z1 = region
        height = z1-z0
    rows = []
    for x,y,z,q in block["controls"]:
        if hole:
            if axis=="x":
                a,b,c = x-40*q, z-20*q, 8*q-y
            else:
                a,b,c = y-32*q, z-20*q, x
            for _ in range(turn):
                a,b = b,-a
        else:
            a,b = (x,y) if side=="east" else (y,x)
            c = z-z0*q
        rows.append((a,b,c,q))
    if not hole and side=="north":
        rows = rows[8:12]+rows[4:8]+rows[0:4]
    profile = []
    for (a,b,c,q),(aa,bb,cc,qq) in zip(rows[::2],rows[1::2], strict=True):
        _require(c==0 and cc==q*height and (a,b,q)==(aa,bb,qq), "extrusion_coefficients")
        profile.append((a,b,q))
    points = [(F(x,q),F(y,q)) for x,y,q in profile]
    b0,a0,_,_,b1,a1 = points
    if hole:
        _require(b0==(10,-10) and b1==(10,10) and a0[0]==-a0[1]
                 and a1[0]==a1[1] and 0<a0[0]<10 and 0<a1[0]<10, "hole_sector_rays")
        result = certify_arc_strip(profile,(0,0),4)
    else:
        _require(b0==(0,0) and b1==(14,0) and a1==(14,8)
                 and a0[0]==a0[1] and 8<a0[0]<14, "fillet_endpoints")
        result = certify_arc_strip(profile,(14,14),6)
    if result["status"]!="valid":
        return result
    return {"status":"valid","det":height*result["determinant_lower_bound"],
            "error":result["radius_error_mm"]}


def design_faces(blocks, exterior):
    """Resolve every exterior side to exactly one fixed design feature."""
    result = []
    planes = (("x",0,(0,8,60)), ("y",1,(0,8,50)), ("z",2,(0,40)))
    for name,side in exterior:
        if side=="v1" and name.startswith(("hole_","elbow_")):
            feature = name[:6] if name.startswith("hole_") else "inside_fillet"
        else:
            trace = extract_trace(blocks[name],side)["controls"]
            hits = [f"{axis}{value}" for axis,index,levels in planes for value in levels
                    if all(pt[index]==value*pt[3] for pt in trace)]
            _require(len(hits)==1, "unclassified_boundary")
            feature = hits[0]
        result.append({"block":name,"side":side,"feature":feature})
    return result


def audit(model):
    """Certify each member and the disjoint tiling of the represented bracket."""
    try:
        _require(isinstance(model,dict) and model.get("schema")=="filleted_bracket_v1"
                 and set(model)=={"schema","blocks","topology","design_faces"}, "bracket_schema")
        specs = dict(_regions())
        sequence = model["blocks"]
        _require(isinstance(sequence,list) and len(sequence)==len(specs), "block_count")
        blocks = {}
        dets,errors = [],[]
        for block in sequence:
            _require(isinstance(block,dict) and set(block)=={"id","controls","degrees"}, "block_keys")
            name = block["id"]
            _require(type(name) is str and name in specs and name not in blocks,"block_identity")
            blocks[name]=block
            ds = block["degrees"]
            _require(isinstance(ds,list) and len(ds)==3 and all(type(d) is int for d in ds),
                     "block_degrees")
            pts = block["controls"]
            _require(isinstance(pts,list) and len(pts) in (8,12) and all(
                isinstance(pt,list) and len(pt)==4 and all(type(x) is int for x in pt)
                and pt[3]>0 for pt in pts), "integer_homogeneous_controls")
            region=specs[name]
            if region[0]=="box":
                dets.append(_straight_box(block,region[1:]))
            else:
                _require(len(pts)==12,"curved_control_count")
                cert = _section(block,region)
                if cert["status"]!="valid":
                    return cert
                dets.append(cert["det"])
                errors.append(cert["error"])
        topology = enumerate_faces(sequence)
        _require(json.dumps(topology, sort_keys=True)==json.dumps(model["topology"], sort_keys=True),
                 "incidence_mismatch")
        _require(len(topology["interfaces"])==55 and len(topology["exterior"])==106,
                 "partition_topology")
        _require(design_faces(blocks,topology["exterior"])==model["design_faces"], "design_faces")
    except (ValueError, TypeError, KeyError) as error:
        return {"status":"invalid","reason":str(error)}
    return {"status":"valid","method":"fixed_extrusion_tiling_with_certified_arc_boundaries",
            "blocks":36,"shared_faces":55,"exterior_sides":106,
            "determinant_lower_bound_mm3":str(min(dets)),
            "radial_error_upper_bound_mm":str(max(errors)),
            "domain":"closed cubes of the fixed 60 x 50 x 40 mm bracket with two through holes"}
