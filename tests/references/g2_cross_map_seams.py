#!/usr/bin/env python3
"""Classify offset samples around S2 seam loci against the exported volume maps."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import g2_cross_junctions as junctions

ROOT = Path(__file__).resolve().parents[2]


def bernstein(n, t):
    return [math.comb(n, i) * t**i * (1.0 - t)**(n-i) for i in range(n+1)]


def bernstein_d(n, t):
    if n == 0:
        return [0.0]
    lower = bernstein(n-1, t)
    return [n * ((lower[i-1] if i else 0.0) - (lower[i] if i < n else 0.0))
            for i in range(n+1)]


def map_xy(block, uvw):
    p, q, r = block["degrees"]
    u, v, w = uvw
    bu, bv, bw = bernstein(p, u), bernstein(q, v), bernstein(r, w)
    du, dv = bernstein_d(p, u), bernstein_d(q, v)
    h = [0.0] * 4
    hu, hv = [0.0] * 4, [0.0] * 4
    for i in range(p+1):
        for j in range(q+1):
            for k in range(r+1):
                c = block["controls"][(i*(q+1)+j)*(r+1)+k]
                b0, bu0, bv0 = bu[i]*bv[j]*bw[k], du[i]*bv[j]*bw[k], bu[i]*dv[j]*bw[k]
                for a in range(4):
                    h[a] += b0*c[a]
                    hu[a] += bu0*c[a]
                    hv[a] += bv0*c[a]
    den = h[3]
    xy = [h[a]/den for a in range(2)]
    ju = [(hu[a]*den-h[a]*hu[3])/(den*den) for a in range(2)]
    jv = [(hv[a]*den-h[a]*hv[3])/(den*den) for a in range(2)]
    return xy, ju, jv


def invert_xy(block, target, z, tol=1e-8):
    p, q, r = block["degrees"]
    if r != 1:
        return None
    controls = block["controls"]
    z0, z1 = controls[0][2]/controls[0][3], controls[1][2]/controls[1][3]
    if not min(z0, z1)-tol <= z <= max(z0, z1)+tol or z1 == z0:
        return None
    w = (z-z0)/(z1-z0)
    for su in (0.2, 0.5, 0.8):
        for sv in (0.2, 0.5, 0.8):
            u, v = su, sv
            for _ in range(32):
                xy, ju, jv = map_xy(block, (u, v, w))
                rx, ry = xy[0]-target[0], xy[1]-target[1]
                if math.hypot(rx, ry) <= tol:
                    if -tol <= u <= 1+tol and -tol <= v <= 1+tol:
                        return (min(1.0, max(0.0, u)), min(1.0, max(0.0, v)),
                                min(1.0, max(0.0, w)))
                    break
                det = ju[0]*jv[1]-jv[0]*ju[1]
                if abs(det) < 1e-16:
                    break
                du = (rx*jv[1]-ry*jv[0])/det
                dv = (ju[0]*ry-ju[1]*rx)/det
                # Keep every evaluation in the closed parameter square; a
                # positive-weight rational map need not have a nonzero
                # denominator for extrapolated parameters outside that square.
                scale = max(1.0, abs(du)/0.75, abs(dv)/0.75)
                accepted = False
                for _ in range(24):
                    next_u, next_v = u-du/scale, v-dv/scale
                    if 0.0 <= next_u <= 1.0 and 0.0 <= next_v <= 1.0:
                        u, v = next_u, next_v
                        accepted = True
                        break
                    scale *= 2.0
                if not accepted:
                    break
    return None


def model_classify(model, point):
    exterior = {(face["block"], face["side"]) for face in model["design_faces"]}
    candidates = []
    for block in model["blocks"]:
        uvw = invert_xy(block, point[:2], point[2])
        if uvw is not None:
            candidates.append((block["id"], uvw))
    if not candidates:
        return {"classification": "OUT", "candidate_maps": 0}
    epsilon = 2e-7
    has_interior_shared_face = False
    has_exterior_face = False
    for block_id, uvw in candidates:
        boundary_sides = []
        for axis, value in enumerate(uvw):
            if abs(value) <= epsilon:
                boundary_sides.append("uvw"[axis]+"0")
            elif abs(value-1.0) <= epsilon:
                boundary_sides.append("uvw"[axis]+"1")
        if any((block_id, side) not in exterior for side in boundary_sides):
            has_interior_shared_face = True
        if any((block_id, side) in exterior for side in boundary_sides):
            has_exterior_face = True
        if not boundary_sides:
            return {"classification": "IN", "candidate_maps": len(candidates),
                    "reason": "strict_block_interior"}
    if has_exterior_face:
        return {"classification": "ON", "candidate_maps": len(candidates),
                "reason": "exterior_design_face"}
    if has_interior_shared_face:
        return {"classification": "IN", "candidate_maps": len(candidates),
                "reason": "interior_shared_block_interface"}
    return {"classification": "INDETERMINATE", "candidate_maps": len(candidates)}


def inside_profile(curve, x, y):
    roots = junctions.find_level_roots(curve, y, axis=1)
    xs = sorted(junctions.jets.rational_jets(curve, t, "right")["position_mm"][0]
                for t in roots)
    return sum(root_x > x for root_x in xs) % 2 == 1


def independent_classify(point, bounds, curve, ylevels, zrange):
    x, y, z = point
    if not all(bounds[a][0] < point[a] < bounds[a][1] for a in range(3)):
        return "OUT"
    in_hole = inside_profile(curve, x, y)
    in_rectangle = ylevels[0] < y < ylevels[1] and zrange[0] < z < zrange[1]
    return "IN" if not in_hole and not in_rectangle else "OUT"


def norm(v):
    s = math.hypot(*v)
    return [a/s for a in v]


def run(model_path, recipe_path):
    model = json.loads(model_path.read_text())
    recipe = json.loads(recipe_path.read_text())
    free = next(s for s in recipe["sketches"] if s["direction"] == [0, 0, 1])
    curve = junctions.sketch_curve(free)
    bounds = model["source_recipe"]["bounds_xyz_mm"]
    _, ylevels, zrange, _ = junctions.sketch_geometry(recipe)
    cross = next(s for s in recipe["sketches"] if s["direction"] == [1, 0, 0])
    offsets = [0.02, -0.02]
    checks, indeterminate = [], 0
    for event in model["sketch_plan"]["crossings"]:
        t = event["parameter"][0]/event["parameter"][1]
        position = junctions.jets.rational_jets(curve, t, "right")["position_mm"]
        tangent = norm(junctions.jets.rational_jets(curve, t, "right")["first_derivative"])
        normal = [-tangent[1], tangent[0]]
        center = [position[0], position[1], 6.0]
        for dy_sign in (1, -1):
            for dn_sign in (1, -1):
                point = [position[0]+dn_sign*offsets[0]*normal[0],
                         position[1]+dy_sign*offsets[0]+dn_sign*offsets[0]*normal[1], 6.0]
                checks.append({"seam": "generatrix", "parameter": event["parameter"],
                               "point_mm": point,
                               "expected": independent_classify(point,bounds,curve,ylevels,zrange),
                               "actual": model_classify(model,point)})
        checks.append({"seam": "generatrix", "parameter": event["parameter"],
                       "point_mm": center, "expected": "ON", "actual": model_classify(model,center)})
    for arc in json.loads((ROOT/"docs/results/g2-cross-junctions.json").read_text())["horizontal_profile_arc_seam_loci"]:
        ta, tb = arc["parameter_interval"]
        t = (ta+tb)/2
        z = arc["sampled_curve_3d_mm"][4][2]
        jet = junctions.jets.rational_jets(curve,t,"right")
        tangent = norm(jet["first_derivative"])
        normal = [-tangent[1],tangent[0]]
        xy = jet["position_mm"]
        for dn_sign in (1,-1):
            for dz_sign in (1,-1):
                point=[xy[0]+dn_sign*offsets[0]*normal[0],
                       xy[1]+dn_sign*offsets[0]*normal[1],z+dz_sign*offsets[0]]
                checks.append({"seam":arc["branch"],"parameter_interval":[ta,tb],
                               "point_mm":point,
                               "expected":independent_classify(point,bounds,curve,ylevels,zrange),
                               "actual":model_classify(model,point)})
        checks.append({"seam":arc["branch"],"parameter_interval":[ta,tb],
                       "point_mm":[xy[0],xy[1],z],"expected":"ON",
                       "actual":model_classify(model,[xy[0],xy[1],z])})
    for check in checks:
        actual = check["actual"]["classification"]
        check["agreement"] = actual == check["expected"]
        if actual not in ("IN","OUT","ON"):
            indeterminate += 1
    agree = sum(row["agreement"] for row in checks)
    return {"schema":"g2_cross_map_seam_classification_v1",
            "model_sha256":hashlib.sha256(model_path.read_bytes()).hexdigest(),
            "source_curve":"raw periodic NURBS recipe, independent sketch parity classifier",
            "point_offset_mm":0.02,
            "query_count":len(checks),"agreement_count":agree,
            "disagreement_count":len(checks)-agree,"indeterminate_count":indeterminate,
            "status":"verified" if agree==len(checks) and indeterminate==0 else "failed",
            "checks":checks,
            "limits":["IN/OUT samples use numerical inverse maps over every volume block and an independent source-curve parity test.",
                      "ON samples are exact seam midpoints and rely on exterior design-face labels; this is not a B-rep or general point-classification certificate."]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model",type=Path)
    parser.add_argument("recipe",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    result=run(args.model,args.recipe)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:result[k] for k in ("status","query_count","agreement_count","disagreement_count","indeterminate_count")}))
    return 0 if result["status"]=="verified" else 1


if __name__=="__main__":
    raise SystemExit(main())
