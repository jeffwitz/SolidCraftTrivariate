#!/usr/bin/env python3
"""Independent exact-rational audit for the curved-sketch gallery models.

All geometric predicates in this file use :class:`fractions.Fraction`; no
constructor certificates or floating point sampling are trusted.  A positive
Bernstein subdivision certificate proves det(DX)>0 on each closed block.
"""
from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path


def q(value):
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return F(int(value[0]), int(value[1]))
    return F(str(value))


def poly_add(a, b, scale=F(1)):
    out = dict(a)
    for key, value in b.items():
        out[key] = out.get(key, F(0)) + scale * value
        if not out[key]:
            del out[key]
    return out


def poly_mul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            key = tuple(x + y for x, y in zip(ka, kb))
            out[key] = out.get(key, F(0)) + va * vb
    return {k: v for k, v in out.items() if v}


def derivative(a, axis):
    out = {}
    for key, value in a.items():
        if key[axis]:
            new = list(key)
            out[tuple(new[:axis] + [new[axis] - 1] + new[axis + 1:])] = value * key[axis]
    return out


def bernstein_to_power(degree, values):
    """Convert one Bernstein coefficient vector to the exact power basis."""
    out = [F(0)] * (degree + 1)
    for i, value in enumerate(values):
        for k in range(i, degree + 1):
            out[k] += value * math.comb(degree, i) * math.comb(degree - i, k - i) * (-1) ** (k - i)
    return out


def tensor_power(degrees, controls, component):
    """Tensor Bernstein net (u-major, then v, then w) to sparse power form."""
    dims = [d + 1 for d in degrees]
    grid = [[[F(0) for _ in range(dims[2])] for _ in range(dims[1])] for _ in range(dims[0])]
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                p = controls[(i * dims[1] + j) * dims[2] + k]
                weight = q(p[3])
                if weight <= 0:
                    raise ValueError("all homogeneous weights must be positive")
                grid[i][j][k] = q(p[component]) / weight
    # Successive exact basis changes along each tensor axis.
    for axis in range(3):
        n = degrees[axis]
        updated = [[[F(0) for _ in range(dims[2])] for _ in range(dims[1])] for _ in range(dims[0])]
        for i in range(dims[0]):
            for j in range(dims[1]):
                for k in range(dims[2]):
                    idx = [i, j, k]
                    if idx[axis] != 0:
                        continue
                    line = []
                    for s in range(n + 1):
                        here = idx.copy(); here[axis] = s
                        line.append(grid[here[0]][here[1]][here[2]])
                    coeffs = bernstein_to_power(n, line)
                    for s, value in enumerate(coeffs):
                        here = idx.copy(); here[axis] = s
                        updated[here[0]][here[1]][here[2]] = value
        grid = updated
    return {(i, j, k): grid[i][j][k]
            for i in range(dims[0]) for j in range(dims[1]) for k in range(dims[2])
            if grid[i][j][k]}


def determinant_poly(block):
    deg = list(map(int, block["degrees"]))
    coords = [tensor_power(deg, block["controls"], c) for c in range(3)]
    jac = [[derivative(coords[row], col) for col in range(3)] for row in range(3)]
    a = poly_mul(jac[0][0], poly_add(poly_mul(jac[1][1], jac[2][2]), poly_mul(jac[1][2], jac[2][1]), -1))
    b = poly_mul(jac[0][1], poly_add(poly_mul(jac[1][0], jac[2][2]), poly_mul(jac[1][2], jac[2][0]), -1))
    c = poly_mul(jac[0][2], poly_add(poly_mul(jac[1][0], jac[2][1]), poly_mul(jac[1][1], jac[2][0]), -1))
    return poly_add(poly_add(a, b, -1), c)


def power_to_bernstein(poly):
    deg = tuple(max((k[a] for k in poly), default=0) for a in range(3))
    net = {}
    for i in range(deg[0] + 1):
        for j in range(deg[1] + 1):
            for k in range(deg[2] + 1):
                index = (i, j, k)
                value = F(0)
                for power, coeff in poly.items():
                    if all(power[a] <= index[a] for a in range(3)):
                        factor = F(1)
                        for a in range(3):
                            factor *= F(math.comb(index[a], power[a]), math.comb(deg[a], power[a]))
                        value += coeff * factor
                net[index] = value
    return deg, net


def split_net(deg, net, axis):
    """Exact de Casteljau split at one half along a tensor direction."""
    dims = [d + 1 for d in deg]
    left, right = {}, {}
    other = [a for a in range(3) if a != axis]
    for x in range(dims[other[0]]):
        for y in range(dims[other[1]]):
            base = [0, 0, 0]; base[other[0]] = x; base[other[1]] = y
            row = []
            for t in range(dims[axis]):
                base[axis] = t; row.append(net[tuple(base)])
            levels = [row]
            while len(levels[-1]) > 1:
                prev = levels[-1]
                levels.append([(prev[i] + prev[i + 1]) / 2 for i in range(len(prev) - 1)])
            lrow = [levels[i][0] for i in range(len(row))]
            rrow = [levels[len(row) - 1 - i][i] for i in range(len(row))]
            for t, (lv, rv) in enumerate(zip(lrow, rrow)):
                base[axis] = t; left[tuple(base)] = lv; right[tuple(base)] = rv
    return left, right


def certify_positive(deg, net, depth=0, max_depth=12, boxes=None):
    vals = list(net.values())
    if min(vals) > 0:
        return depth, min(vals), boxes or 1
    if min(vals) <= 0 and max(vals) <= 0:
        raise ValueError("Jacobian determinant is nonpositive on a Bernstein control hull")
    if depth >= max_depth:
        raise ValueError("strict-positive Jacobian not certified within subdivision budget")
    dimensions = [d + 1 for d in deg]
    variation = []
    for axis in range(3):
        if deg[axis] == 0:
            variation.append(F(0))
            continue
        worst = F(0)
        other = [a for a in range(3) if a != axis]
        for i in range(dimensions[other[0]]):
            for j in range(dimensions[other[1]]):
                index = [0, 0, 0]
                index[other[0]], index[other[1]] = i, j
                for t in range(deg[axis]):
                    index[axis] = t; left = net[tuple(index)]
                    index[axis] = t + 1; right = net[tuple(index)]
                    worst = max(worst, abs(right - left))
        variation.append(worst)
    axis = max(range(3), key=lambda a: variation[a])
    if deg[axis] == 0 or variation[axis] == 0:
        raise ValueError("nonpositive constant determinant")
    lnet, rnet = split_net(deg, net, axis)
    dl, ml, bl = certify_positive(deg, lnet, depth + 1, max_depth)
    dr, mr, br = certify_positive(deg, rnet, depth + 1, max_depth)
    return max(dl, dr), min(ml, mr), bl + br


def exact_integral(poly):
    result = F(0)
    for powers, coeff in poly.items():
        result += coeff / math.prod(p + 1 for p in powers)
    return result


def face_net(block, side):
    axis, end = "uvw".index(side[0]), int(side[1])
    deg = list(map(int, block["degrees"]))
    dims = [d + 1 for d in deg]
    controls = []
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                idx = [i, j, k]
                if idx[axis] == (deg[axis] if end else 0):
                    controls.append(tuple(q(x) for x in block["controls"][(i*dims[1]+j)*dims[2]+k]))
    axes = [a for a in range(3) if a != axis]
    shape = [dims[a] for a in axes]
    return tuple(shape), controls


def design_wall_curve_rows(model, feature_id):
    """Return actual cut-wall quadratic traces at both axial faces."""
    blocks = {b["id"]: b for b in model["blocks"]}
    rows = set()
    feature_name = f"cut/{feature_id}/wall"
    for face in model.get("design_faces", []):
        if face.get("feature") != feature_name:
            continue
        block = blocks[face["block"]]
        shape, net = face_net(block, face["side"])
        if 3 not in shape or 2 not in shape:
            continue
        curve_axis = shape.index(3)
        sweep_axis = shape.index(2)
        for sweep_index in (0, 1):
            pts = []
            for curve_index in range(3):
                idx = [0, 0]
                idx[curve_axis] = curve_index
                idx[sweep_axis] = sweep_index
                p = net[idx[0] * shape[1] + idx[1]]
                pts.append((q(p[0]) / q(p[3]), q(p[1]) / q(p[3])))
            rows.add(tuple(pts))
    return rows


def orient_face(shape, points, orientation):
    rows, cols = shape
    net = [list(points[r*cols:(r+1)*cols]) for r in range(rows)]
    if orientation.get("flip_first", False): net.reverse()
    if orientation.get("flip_second", False): net = [list(reversed(row)) for row in net]
    if orientation.get("swap", False): net = [list(row) for row in zip(*net)]
    return (len(net), len(net[0])), [p for row in net for p in row]


def check_interfaces(model):
    blocks = {b["id"]: b for b in model["blocks"]}
    interfaces = model.get("topology", {}).get("interfaces", [])
    for itf in interfaces:
        a, sa = itf["first"]; b, sb = itf["second"]
        o = itf.get("orientation", {})
        sha, pa = face_net(blocks[a], sa); shb, pb = face_net(blocks[b], sb)
        sha, pa = orient_face(sha, pa, {"flip_first": o.get("flip_first", False), "flip_second": o.get("flip_second", False), "swap": o.get("swap", False)})
        if sha != shb or pa != pb:
            raise ValueError(f"oriented interface controls disagree: {a}:{sa} / {b}:{sb}")
    return len(interfaces)


def volume_audit(model):
    total = F(0)
    certs = []
    for block in model.get("blocks", []):
        ws = {q(p[3]) for p in block["controls"]}
        if len(ws) != 1 or next(iter(ws)) <= 0:
            raise ValueError(f"block {block['id']} does not have equal positive homogeneous weights")
        det = determinant_poly(block)
        deg, net = power_to_bernstein(det)
        levels, lower, boxes = certify_positive(deg, net)
        total += exact_integral(det)
        certs.append({"block": block["id"], "determinant_degree": deg,
                      "bernstein_subdivision_depth": levels, "certified_lower_bound": str(lower),
                      "subboxes": boxes})
    return total, certs


def check_g1_ovals(recipe, model):
    """Reconstruct four quadratic arcs from each source rectangle.

    The two exact construction rules are the current family contract: a
    through-hole scales corners by 2/3 and uses source-span/3 control sag;
    a pocket scales by 4/5 and uses source-span/5. Controls are checked
    against all model net rows/columns, while endpoint tangent cross products
    classify the actual junctions without assuming G1 for pocket corners.
    """
    scaffold = recipe["scaffold"]
    outer = scaffold["outer"]["vertices"]
    bounds = [[min(p[a] for p in outer), max(p[a] for p in outer)] for a in range(2)]
    height = q(scaffold.get("height_mm", model.get("height_mm", 0)))
    results = []
    for cut in scaffold.get("cuts", []):
        verts = cut["vertices"]
        minx, maxx = min(p[0] for p in verts), max(p[0] for p in verts)
        miny, maxy = min(p[1] for p in verts), max(p[1] for p in verts)
        if minx >= maxx or miny >= maxy:
            raise ValueError(f"feature {cut['id']} has nonpositive source widths")
        zi = cut.get("z_interval_mm", [0, height])
        through = q(zi[0]) == 0 and q(zi[1]) == height
        corner_scale = F(2, 3) if through else F(4, 5)
        sag_scale = F(1, 3) if through else F(1, 5)
        cx, cy = F(minx + maxx, 2), F(miny + maxy, 2)
        corners = [(cx + corner_scale * (q(x) - cx), cy + corner_scale * (q(y) - cy)) for x, y in verts]
        expected = []
        wall_rows = design_wall_curve_rows(model, cut["id"])
        if not wall_rows:
            raise ValueError(f"model has no design wall faces for cut {cut['id']}")
        shift = list(map(q, scaffold.get("top_shift_mm", [0, 0])))
        for i, a in enumerate(corners):
            b = corners[(i + 1) % len(corners)]
            if a[1] == b[1]:
                sag = (q(maxy) - q(miny)) * sag_scale
                side = 1 if a[1] > cy else -1
                control = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + side * sag)
            elif a[0] == b[0]:
                sag = (q(maxx) - q(minx)) * sag_scale
                side = 1 if a[0] > cx else -1
                control = ((a[0] + b[0]) / 2 + side * sag, (a[1] + b[1]) / 2)
            else:
                raise ValueError(f"cut {cut['id']} is not an orthogonal rectangle")
            curve = (a, control, b)
            for axial_fraction in (F(0), F(1)):
                shifted = tuple((p[0] + shift[0]*axial_fraction, p[1] + shift[1]*axial_fraction) for p in curve)
                if not any(shifted == row or shifted == tuple(reversed(row)) for row in wall_rows):
                    raise ValueError(f"model wall controls do not reproduce source arc {cut['id']} edge {i} at w={axial_fraction}")
            expected.append(curve)
        joints = []
        for i in range(4):
            prev, cur = expected[i - 1], expected[i]
            # Orient both arcs along the source loop and compare incoming and
            # outgoing tangent vectors at their shared source corner.
            vin = (prev[2][0] - prev[1][0], prev[2][1] - prev[1][1])
            vout = (cur[1][0] - cur[0][0], cur[1][1] - cur[0][1])
            cross = vin[0] * vout[1] - vin[1] * vout[0]
            dot = vin[0] * vout[0] + vin[1] * vout[1]
            c1 = vin == vout
            joints.append({"cross_product": str(cross), "dot_product": str(dot),
                           "g1": cross == 0 and dot > 0, "c1": c1})
        expected_class = "C1" if through else "G0"
        if through and not all(j["g1"] for j in joints):
            raise ValueError(f"through-hole {cut['id']} arc junctions are not oriented G1")
        if through and not all(j["c1"] for j in joints):
            raise ValueError(f"through-hole {cut['id']} arc endpoint derivatives do not match exactly")
        results.append({"feature": cut["id"], "expected_continuity": expected_class,
                        "four_source_arcs_matched": len(expected), "junctions": joints,
                        "source_widths_mm": [maxx-minx, maxy-miny],
                        "quadratic_arcs_mm": [[[float(c) for c in point] for point in arc]
                                               for arc in expected]})
    # The bowed outer stock side is checked in the actual outer design face.
    ymin = min(q(p[1]) for p in outer)
    for i, a0 in enumerate(outer):
        b0 = outer[(i + 1) % len(outer)]
        if q(a0[1]) == q(b0[1]) == ymin:
            a, b = tuple(map(q, a0)), tuple(map(q, b0))
            # Constructor contract: quadratic sag is one eighth of stock height.
            control = ((a[0] + b[0]) / 2, ymin - (bounds[1][1] - bounds[1][0]) / 8)
            stock_rows = set()
            for face in model.get("design_faces", []):
                if face.get("feature") != "outer/body":
                    continue
                block = {b["id"]: b for b in model["blocks"]}[face["block"]]
                shape, net = face_net(block, face["side"])
                if 3 in shape and 2 in shape:
                    ca, sa = shape.index(3), shape.index(2)
                    for si in (0, 1):
                        row=[]
                        for ci in range(3):
                            idx=[0,0]; idx[ca]=ci; idx[sa]=si
                            p=net[idx[0]*shape[1]+idx[1]]
                            row.append((q(p[0])/q(p[3]),q(p[1])/q(p[3])))
                        stock_rows.add(tuple(row))
            if not any((a, control, b) == row or (a, control, b) == tuple(reversed(row)) for row in stock_rows):
                raise ValueError("model is missing the independently predicted bowed outer stock side")
    return results


def green_volume(recipe):
    """Integrate source contours by Green's theorem, layer by layer."""
    scaffold = recipe["scaffold"]
    outer = scaffold["outer"]["vertices"]
    cuts = scaffold.get("cuts", [])
    height = q(scaffold.get("height_mm"))
    zcuts = sorted({F(0), height, *(q(z) for cut in cuts for z in cut.get("z_interval_mm", []))})
    total = F(0)
    for z0, z1 in zip(zcuts, zcuts[1:]):
        if z1 <= z0:
            continue
        mid = (z0 + z1) / 2
        area = green_loop(outer, outer_side=True)
        for cut in cuts:
            zlo, zhi = map(q, cut.get("z_interval_mm", [0, height]))
            if zlo < mid < zhi:
                area -= green_loop(cut["vertices"], cut=True,
                                   through=(zlo == 0 and zhi == height))
        if area <= 0:
            raise ValueError("source contour has nonpositive material area")
        total += area * (z1 - z0)
    return total


def green_loop(vertices, outer_side=False, cut=False, through=False):
    """Exact 1/2 integral (x dy-y dx) for lines and quadratic Bezier arcs."""
    pts = [tuple(map(q, p)) for p in vertices]
    minx, maxx = min(p[0] for p in pts), max(p[0] for p in pts)
    miny, maxy = min(p[1] for p in pts), max(p[1] for p in pts)
    cx, cy = (minx + maxx) / 2, (miny + maxy) / 2
    height = maxy - miny
    scale = F(2, 3) if through else F(4, 5)
    sag_scale = F(1, 3) if through else F(1, 5)
    contour = [(cx + scale*(x-cx), cy + scale*(y-cy)) for x,y in pts] if cut else pts
    total = F(0)
    for i, a in enumerate(contour):
        b = contour[(i + 1) % len(contour)]
        control = None
        if cut and a[1] == b[1]:
            side = 1 if a[1] > cy else -1
            sag = (maxy - miny) * sag_scale
            control = ((a[0]+b[0])/2, (a[1]+b[1])/2 + side*sag)
        elif cut and a[0] == b[0]:
            side = 1 if a[0] > cx else -1
            sag = (maxx - minx) * sag_scale
            control = ((a[0]+b[0])/2 + side*sag, (a[1]+b[1])/2)
        elif outer_side and a[1] == b[1] == miny:
            control = ((a[0]+b[0])/2, miny-height/8)
        if control is None:
            # For a line B(t)=a+(b-a)t, Green contribution is cross(a,b)/2.
            total += (a[0]*b[1] - a[1]*b[0]) / 2
        else:
            # Convert quadratic Bezier coordinates to power coefficients.
            x = [a[0], 2*(control[0]-a[0]), a[0]-2*control[0]+b[0]]
            y = [a[1], 2*(control[1]-a[1]), a[1]-2*control[1]+b[1]]
            dx, dy = [x[1], 2*x[2]], [y[1], 2*y[2]]
            integ = F(0)
            for j, xv in enumerate(x):
                for k, ydv in enumerate(dy): integ += xv*ydv/F(j+k+1)
            for j, yv in enumerate(y):
                for k, xdv in enumerate(dx): integ -= yv*xdv/F(j+k+1)
            total += integ/2
    return total


def layer_minima(model, scaffold):
    """Count required arc junctions on a common footprint with whole-face joins."""
    cuts = scaffold["cuts"]
    levels = sorted({0, scaffold["height_mm"], *(v for c in cuts for v in c["z_interval_mm"])})
    intervals = []
    for block in model["blocks"]:
        zs = [q(p[2]) / q(p[3]) for p in block["controls"]]
        intervals.append((min(zs), max(zs)))
    result = []
    for lo, hi in zip(levels, levels[1:]):
        holes = [c for c in cuts if c["z_interval_mm"][0] <= lo and hi <= c["z_interval_mm"][1]]
        boundary = len(scaffold["outer"]["vertices"]) + 4 * len(holes)
        interior = 4 * (len(cuts) - len(holes))
        minimum = F(boundary, 2) + interior + len(holes) - 1
        observed = intervals.count((lo, hi))
        if observed != minimum:
            raise ValueError("layer does not attain its marked-junction lower bound")
        result.append({"z_mm": [lo, hi], "required_boundary_vertices": boundary,
                       "required_interior_vertices": interior, "holes": len(holes),
                       "minimum_blocks": observed})
    if sum(row["minimum_blocks"] for row in result) != len(intervals):
        raise ValueError("found a block spanning a layer event")
    return result


def read_step(path, target_volume):
    """Optional independent STEP qualification using FreeCAD's OCCT reader."""
    import FreeCAD  # initialize the geometry kernel before importing Part
    import Part
    shape = Part.read(str(path))
    if len(shape.Solids) != 1 or not shape.isValid() or not all(s.isClosed() for s in shape.Shells):
        raise ValueError("STEP does not describe a valid solid with closed shells")
    volume = shape.Volume
    error = abs(volume - target_volume) / target_volume
    if error > 1e-7:
        raise ValueError(f"STEP relative volume error {error} exceeds 1e-7")
    return {"reader": "OCCT " + Part.OCC_VERSION, "closed_shells": len(shape.Shells),
            "valid_solids": len(shape.Solids), "volume_mm3": volume, "relative_volume_error": error}


def run(recipe_path, model_path, expected_block_count=None):
    recipe, model = json.loads(recipe_path.read_text()), json.loads(model_path.read_text())
    if recipe.get("family") != "curved_sketch_layers_v1":
        raise ValueError("unexpected recipe family")
    scaffold = recipe.get("scaffold")
    if not isinstance(scaffold, dict):
        raise ValueError("curved recipe must wrap a scaffold recipe")
    count = len(model.get("blocks", []))
    expected = (expected_block_count or model.get("expected_block_count")
                or recipe.get("expected_block_count"))
    if expected is not None and count != int(expected):
        raise ValueError(f"block count mismatch: {count} != {expected}")
    common_faces = check_interfaces(model)
    exact_volume, certificates = volume_audit(model)
    g1 = check_g1_ovals(recipe, model)
    reference_volume = green_volume(recipe)
    if exact_volume != reference_volume:
        raise ValueError(f"analytic Green volume mismatch: map={exact_volume}, source={reference_volume}")
    return {"schema": "curved_sketch_gallery_audit_v1", "status": "pass",
            "block_count": count, "common_faces_checked": common_faces,
            "layer_lower_bounds": layer_minima(model, scaffold),
            "exact_map_volume_mm3": str(exact_volume), "green_volume_mm3": str(reference_volume),
            "g1_oval_arcs": g1, "jacobian_certificates": certificates}


def main():
    # Some FreeCAD builds initialize the main Python namespace during import.
    serialize = json.dumps
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-block-count", type=int)
    parser.add_argument("--step", type=Path, help="also read the STEP solid with FreeCAD/OCCT")
    args = parser.parse_args()
    audit = run(args.recipe, args.model, args.expected_block_count)
    if args.step:
        audit["step_readback"] = read_step(args.step, float(F(audit["exact_map_volume_mm3"])))
    args.output.write_text(serialize(audit, indent=2) + "\n")


if __name__ == "__main__":
    main()
