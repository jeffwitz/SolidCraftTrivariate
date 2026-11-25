"""Check the stored ruled walls against both source one-sided jets in exact rationals."""
import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

from g2_cross_brep_solidcore import homogeneous_jets


def geometric_jet(h, dh, ddh):
    """Project velocity to a direction and acceleration to the curvature vector."""
    w = h[-1]
    p = tuple(x / w for x in h[:2])
    v = tuple((dh[i] - p[i] * dh[-1]) / w for i in range(2))
    a = tuple((ddh[i] - p[i] * ddh[-1] - 2*v[i]*dh[-1]) / w for i in range(2))
    speed2 = sum(x*x for x in v)
    if speed2 == 0:
        raise ValueError("Singular wall endpoint")
    tangent_accel = sum(x*y for x, y in zip(v, a)) / speed2
    curvature = tuple((a[i] - tangent_accel * v[i]) / speed2 for i in range(2))
    return p, v, curvature


def stored_jet(net, end):
    n = len(net) - 1
    h, near, far = net[-1:-4:-1] if end else net[:3]
    direction = -1 if end else 1
    dh = [direction * n * (x-y) for x,y in zip(near,h)]
    ddh = [n*(n-1)*(z-2*y+x) for x,y,z in zip(h,near,far)]
    return geometric_jet(h,dh,ddh)


def run(model):
    profile = next(s for s in model["source_recipe"]["sketches"]
                   if s["direction"] == [0,0,1])["loops"][0]["curve"]
    knots = profile["knots"]
    p = profile["degree"]
    lo, hi = knots[p], knots[-p-1]
    levels = set(map(F, knots[p:-p]))
    levels.update(F(*c["parameter"]) for c in model["sketch_plan"]["crossings"])
    reference = {}
    for t in levels - {hi}:
        a = geometric_jet(*homogeneous_jets(profile, hi if t == lo else t, "left"))
        b = geometric_jet(*homogeneous_jets(profile, t, "right"))
        assert a[0] == b[0] and a[2] == b[2]
        assert a[1][0]*b[1][1] == a[1][1]*b[1][0]
        assert sum(x*y for x,y in zip(a[1],b[1])) > 0
        reference[b[0]] = b

    blocks = {b["id"]:b for b in model["blocks"]}
    layers = defaultdict(Counter)
    checks = faces = 0
    covered = set()
    for face in model["design_faces"]:
        if face["feature"] != "nurbs_z_hole/wall":
            continue
        block = blocks[face["block"]]
        assert face["side"] == "v1" and block["degrees"][1:] == [1,1]
        pts = [tuple(map(F,row)) for row in block["controls"]]
        bottom, top = pts[2::4], pts[3::4]
        z = (bottom[0][2]/bottom[0][3], top[0][2]/top[0][3])
        assert z[1] > z[0]
        net = []
        for a,b in zip(bottom,top):
            assert (a[0],a[1],a[3]) == (b[0],b[1],b[3])
            assert a[2] == a[3]*z[0] and b[2] == b[3]*z[1]
            net.append((a[0],a[1],a[3]))
        orientations = []
        for end in (False,True):
            point,v,curvature = stored_jet(net,end)
            _, source_v, source_curvature = reference[point]
            assert v[0]*source_v[1] == v[1]*source_v[0]
            assert curvature == source_curvature
            orientations.append(sum(x*y for x,y in zip(v,source_v)) > 0)
            layers[z][point] += 1
            covered.add(point)
            checks += 1
        assert orientations[0] == orientations[1]
        faces += 1
    assert covered == reference.keys()
    counts = {}
    for z, endpoints in sorted(layers.items()):
        assert all(n in (1,2) for n in endpoints.values())
        counts[str(tuple(map(int,z)))] = {"smooth_paired_joins":sum(n==2 for n in endpoints.values()),
                                         "ends_at_transverse_cut":sum(n==1 for n in endpoints.values())}
    return {"method":"Exact Fraction endpoint jets of saved wall controls vs source homogeneous de Boor derivatives.",
            "curved_wall_faces":faces,"exact_wall_endpoint_checks":checks,"source_locations":len(reference),
            "layer_join_counts":counts,"curvature_vector_difference_inv_mm":0,
            "squared_normal_direction_cross_product":0,"generatrix_curvature_inv_mm":0,
            "scope":"The ruled wall is exactly a translation along Z. These checks apply to retained NURBS walls; intersections with cutter planes are sharp edges."}


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    report=run(json.loads(args.model.read_text()))
    args.output.write_text(json.dumps(report,indent=2)+"\n")
