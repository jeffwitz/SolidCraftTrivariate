#!/usr/bin/env python3
"""Independent source-curve oracle for the S3 face-pocket prospective case.

No V-Rep blocks, STEP, B-Rep, or constructor output is read. The original
periodic rational cubic is evaluated directly from its Cox/de Boor data.
"""
from __future__ import annotations

import json
import hashlib
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
RECIPE = ROOT / "tests/references/g2_cross_sketches.recipe.json"
S3_RECIPES = {
    "interior_principal": ROOT / "tests/references/s3_face_pocket.recipe.json",
    "open_stock_edge_witness": ROOT / "tests/references/s3_face_pocket_open_side.recipe.json",
}
REPORT = ROOT / "docs/results/g2-face-pocket-geometry-20261125.json"
Z_SUPPORT, Z_FLOOR = 4.0, 2.0


def curve_data():
    data = json.loads(RECIPE.read_text())
    sketch = data["sketches"][0]
    loop = sketch["loops"][0]
    curve = loop["curve"]
    return curve["degree"], curve["knots"], curve["controls"], sketch["id"], loop["id"]


def fixture_cases():
    cases, inputs = {}, {}
    for name, path in S3_RECIPES.items():
        recipe = json.loads(path.read_text())
        op = recipe["operation"]
        loop = op["loop"]
        primitive = loop["primitive"]
        frac = lambda pair: float(pair[0]) / float(pair[1])
        x0, y0 = (frac(v) for v in primitive["min"])
        x1, y1 = (frac(v) for v in primitive["max"])
        if recipe.get("family") != "s3_face_pocket_v1":
            raise ValueError(f"unexpected fixture family in {path.name}")
        if op.get("support_feature") != "rect_x/zmin" or op.get("frame", {}).get("origin_mm") != [0, 0, 4]:
            raise ValueError(f"unexpected support identity in {path.name}")
        if op.get("direction") != [0, 0, -1]:
            raise ValueError(f"unexpected pocket normal in {path.name}")
        cases[name] = (x0, x1, y0, y1)
        inputs[name] = {"recipe": str(path.relative_to(ROOT)),
                        "recipe_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "family": recipe["family"], "support_feature": op["support_feature"],
                        "loop_id": loop["id"], "primitive": primitive,
                        "frame_origin_mm": op["frame"]["origin_mm"],
                        "direction": op["direction"], "depth_mm": op["depth_mm"],
                        "parent_block_count": len(recipe["parent_model"]["blocks"])}
    return cases, inputs


def de_boor(degree, knots, controls, t):
    xw, yw, w = homogeneous_de_boor(degree, knots, controls, t)
    return xw / w, yw / w


def homogeneous_de_boor(degree, knots, controls, t):
    n = len(controls) - 1
    if t >= knots[n + 1]:
        span = n
    else:
        span = next(i for i in range(degree, n + 1)
                    if knots[i] <= t < knots[i + 1])
    d = [list(map(float, controls[j])) for j in range(span-degree, span+1)]
    for r in range(1, degree + 1):
        for j in range(degree, r - 1, -1):
            i = span - degree + j
            den = knots[i + degree - r + 1] - knots[i]
            alpha = 0.0 if den == 0 else (t - knots[i]) / den
            d[j] = [(1-alpha)*d[j-1][k] + alpha*d[j][k]
                    for k in range(3)]
    return d[degree]


def source_jet(t):
    """Euclidean first derivative from the derivative homogeneous spline."""
    p = DEG
    derivative_controls = []
    for i in range(len(CONTROLS)-1):
        den = KNOTS[i+p+1]-KNOTS[i+1]
        derivative_controls.append([p*(CONTROLS[i+1][k]-CONTROLS[i][k])/den
                                    for k in range(3)])
    original = json.loads(RECIPE.read_text())["sketches"][0]["loops"][0]["curve"]
    n = len(original["controls"])-1
    span = n if t >= original["knots"][n+1] else next(
        i for i in range(p, n+1) if original["knots"][i] <= t < original["knots"][i+1])
    d = [list(map(float, original["controls"][j])) for j in range(span-p, span+1)]
    for r in range(1,p+1):
        for j in range(p,r-1,-1):
            i=span-p+j
            a=(t-original["knots"][i])/(original["knots"][i+p-r+1]-original["knots"][i])
            d[j]=[(1-a)*d[j-1][k]+a*d[j][k] for k in range(3)]
    X,Y,W=d[p]
    dknots=original["knots"][1:-1]
    Xp,Yp,Wp=homogeneous_de_boor(p-1,dknots,derivative_controls,t)
    return ((Xp*W-X*Wp)/(W*W), (Yp*W-Y*Wp)/(W*W))


def arc_length(a, b):
    value, error = quad(lambda t: math.hypot(*source_jet(t)), a, b,
                        epsabs=2e-13, epsrel=2e-13, limit=200)
    return value, error


DEG, KNOTS, CONTROLS, SOURCE_SKETCH_ID, SOURCE_LOOP_ID = curve_data()
N_CONTROL = len(CONTROLS)
SPANS = [(float(a), float(b)) for a, b in
         zip(KNOTS[DEG:N_CONTROL], KNOTS[DEG+1:N_CONTROL+1]) if b > a]
RECTS, FIXTURE_META = fixture_cases()


def crossings_for_x(x):
    roots = []
    for a, b in SPANS:
        xa, ya = de_boor(DEG, KNOTS, CONTROLS, a)
        xb, yb = de_boor(DEG, KNOTS, CONTROLS, b)
        fa, fb = xa-x, xb-x
        if abs(fa) < 2e-13:
            roots.append((a, ya))
        if fa * fb < 0:
            t = brentq(lambda q: de_boor(DEG, KNOTS, CONTROLS, q)[0]-x,
                       a, b, xtol=5e-15)
            roots.append((t, de_boor(DEG, KNOTS, CONTROLS, t)[1]))
        if abs(fb) < 2e-13:
            roots.append((b, yb))
    uniq = []
    for t, y in sorted(roots):
        if not uniq or abs(t-uniq[-1][0]) > 1e-10:
            uniq.append((t, y))
    ys = sorted(y for _, y in uniq)
    return ys


def vertical_hole_length(x, ylo, yhi):
    ys = crossings_for_x(x)
    return sum(max(0.0, min(yhi, b)-max(ylo, a))
               for a, b in zip(ys[::2], ys[1::2]))


def x_limits_at_y(y):
    return sorted(x for _, x in parameters_at_y(y))


def parameters_at_y(y):
    roots = []
    for a, b in SPANS:
        ya = de_boor(DEG, KNOTS, CONTROLS, a)[1]
        yb = de_boor(DEG, KNOTS, CONTROLS, b)[1]
        fa, fb = ya-y, yb-y
        if abs(fa) < 2e-12:
            roots.append((a, de_boor(DEG, KNOTS, CONTROLS, a)[0]))
        if fa*fb < 0:
            t = brentq(lambda q: de_boor(DEG, KNOTS, CONTROLS, q)[1]-y,
                       a, b, xtol=5e-15)
            roots.append((t, de_boor(DEG, KNOTS, CONTROLS, t)[0]))
        if abs(fb) < 2e-12:
            roots.append((b, de_boor(DEG, KNOTS, CONTROLS, b)[0]))
    vals = []
    for t, x in sorted(roots):
        if not vals or abs(t-vals[-1][0]) > 1e-10:
            vals.append((t, x))
    return vals


def integrate_hole(xlo, xhi, ylo, yhi):
    # Split quadrature at x-extrema and all curve-span endpoint x values.
    cuts = [xlo, xhi]
    for a, b in SPANS:
        for t in (a, (a+b)/2, b):
            x, _ = de_boor(DEG, KNOTS, CONTROLS, t)
            if xlo < x < xhi:
                cuts.append(x)
    cuts = sorted(set(cuts))
    total, errors = 0.0, 0.0
    for a, b in zip(cuts, cuts[1:]):
        value, err = quad(lambda x: vertical_hole_length(x, ylo, yhi),
                          a, b, epsabs=2e-11, epsrel=2e-12, limit=300)
        total += value
        errors += err
    return total, errors


def classify_source(x, y):
    # Vertical ray to +infinity, using source spline spans only.
    hits = crossings_for_x(x)
    return "inside_preexisting_nurbs_void" if any(v > y for v in hits) and \
        (sum(v > y for v in hits) % 2 == 1) else "material_before_pocket"


def main():
    depth = Z_SUPPORT-Z_FLOOR
    cases = {}
    for name, (x0, x1, y0, y1) in RECTS.items():
        hole_area, quad_error = integrate_hole(x0, x1, y0, y1)
        rectangle_area = (x1-x0)*(y1-y0)
        material_area = rectangle_area-hole_area
        lower_hits = x_limits_at_y(y0)
        upper_hits = x_limits_at_y(y1)
        def material_trace_length(hits):
            # For this single simple loop, the interval between the two roots is
            # the old hole; the complement within the footprint is material.
            if len(hits) < 2:
                return x1-x0
            return max(0.0, min(x1, hits[0])-x0) + max(0.0, x1-max(x0, hits[1]))
        lower_wall_len = material_trace_length(lower_hits)
        upper_wall_len = material_trace_length(upper_hits)
        x_min_is_stock_edge = x0 == 0.0
        x_min_wall_area = 0.0 if x_min_is_stock_edge else depth*(y1-y0)
        side_wall_area = depth*(lower_wall_len+upper_wall_len)+x_min_wall_area
        # Derive both source parameters from the fixture's exact rectangle
        # ordinates and select the left-side crossings inside its x interval.
        candidates_start = [(t, x) for t, x in parameters_at_y(y1) if x0 <= x <= x1]
        candidates_end = [(t, x) for t, x in parameters_at_y(y0) if x0 <= x <= x1]
        cut_start = min(candidates_start, key=lambda item: abs(item[0]-3.0))[0]
        cut_end = min((item[0] for item in candidates_end if item[0] > cut_start),
                      key=lambda t: t-cut_start)
        edge_start, edge_end = 2.5, 3.5
        removed_start, removed_end = max(edge_start, cut_start), min(edge_end, cut_end)
        removed_start = min(removed_start, removed_end)
        left_end = min(edge_end, removed_start)
        right_start = max(edge_start, removed_end)
        before_edge_len, before_edge_err = arc_length(edge_start, edge_end)
        retained_left_len, retained_left_err = arc_length(edge_start, left_end)
        truncated_len, truncated_err = arc_length(removed_start, removed_end)
        retained_right_len, retained_right_err = arc_length(right_start, edge_end)
        # Existing YZ cross-cut occupies this y band only at z>=4; below z=4
        # the source NURBS hole is the only pre-existing void.
        witnesses = []
        for p in ((x0+(x1-x0)*0.1, y0+(y1-y0)*0.1, 3),
                   (x0+(x1-x0)*0.5, y0+(y1-y0)*0.5, 3),
                   (x1-1e-4, y1-1e-4, 3), (x1, (y0+y1)/2, 3)):
            before = classify_source(p[0], p[1])
            covered = x0 < p[0] < x1 and y0 < p[1] < y1
            after = "pocket_void" if covered and before == "material_before_pocket" else before
            witnesses.append({"point_mm": p, "source_class_before_pocket": before,
                              "source_class_after_pocket": after})
        cases[name] = {
            "footprint_xy_mm": [[x0, x1], [y0, y1]],
            "recipe_identity_and_inputs": FIXTURE_META[name],
            "areas_mm2": {"footprint": rectangle_area, "preexisting_nurbs_void_inside_footprint": hole_area,
                          "material_cut_by_pocket": material_area,
                          "quadrature_error_estimate_void_area": quad_error},
            "volumes_mm3": {"incremental_removed": material_area*depth,
                            "pocket_envelope": rectangle_area*depth,
                            "envelope_already_void_before_pocket": hole_area*depth},
            "old_arc_at_z4": {"status": "the source wall arc is truncated at the footprint crossings; the portion adjacent to pre-existing hole becomes internal to the merged void",
                              "source_parameter_interval_removed": [cut_start, cut_end],
                              "old_floor_left_edge_arc_lengths_mm": {
                                  "edge_parameter_interval": [edge_start, edge_end],
                                  "old_edge_length_before_pocket": before_edge_len,
                                  "retained_left_portion_interval": [edge_start, left_end],
                                  "retained_left_portion_length": retained_left_len,
                                  "truncated_portion_interval": [removed_start, removed_end],
                                  "truncated_portion_length": truncated_len,
                                  "retained_right_portion_interval": [right_start, edge_end],
                                  "retained_right_portion_length": retained_right_len,
                                  "truncated_portion_relocated_on_floor_z2_length": truncated_len,
                                  "partition_sum_check": retained_left_len+truncated_len+retained_right_len,
                                  "partition_residual_mm": before_edge_len-(retained_left_len+truncated_len+retained_right_len),
                                  "quadrature_error_estimates_mm": {"old_edge": before_edge_err,
                                      "retained_left": retained_left_err, "truncated": truncated_err,
                                      "retained_right": retained_right_err}},
                              "closed_profile_perimeter_mm_t0_to_8": arc_length(0.0, 8.0)[0],
                              "source_curve_y_band_intersections_x_at_lower_upper": {
                                  "lower_y": lower_hits, "upper_y": upper_hits}},
            "new_surfaces": {"floor": {"z_mm": Z_FLOOR, "exposed_area_mm2": material_area},
                             "side_walls": {"y_min_wall_length_mm": lower_wall_len, "y_max_wall_length_mm": upper_wall_len,
                                            "x_min_wall_area_mm2": x_min_wall_area,
                                            "total_exposed_side_wall_area_mm2": side_wall_area,
                                            "x_max_wall": "absent; x=75 is inside the existing NURBS void"},
                             "side_wall_line_intervals_x": {"lower_y": lower_hits, "upper_y": upper_hits}},
            "classification_witnesses": witnesses,
        }
    report = {
        "family": "s3_face_pocket_v1",
        "method": "direct periodic rational B-spline source evaluation; Cox/de Boor; vertical-ray intersections; adaptive 1D quadrature; no BRep or constructor geometry",
        "source_recipe": "tests/references/g2_cross_sketches.recipe.json",
        "source_profile": {"degree": DEG, "knots": KNOTS, "controls_homogeneous": CONTROLS,
                            "active_domain": [0, 8], "source_sketch_id": SOURCE_SKETCH_ID,
                            "source_loop_id": SOURCE_LOOP_ID},
        "pocket": {"support_z_mm": Z_SUPPORT, "floor_z_mm": Z_FLOOR, "normal": [0, 0, -1],
                   "intermediate_y_exact_at_t_13_over_4": "5893/266"},
        "cases": cases,
        "limits": ["Boundary coincidences at y=y0 and x=75 are treated as source boundaries; quadrature is numerical and does not certify topology."],
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"report": str(REPORT), "cases": cases}, indent=2))


if __name__ == "__main__":
    main()
