#!/usr/bin/env python3
"""Independent 3D seam-locus oracle for perpendicular NURBS/rectangular cuts.

The source periodic curve is evaluated from its sketch controls. No constructor
intersection, seam, projection, or classification output is used here.
"""
from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import math
import re
import subprocess
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "g2_plate_jets_source", ROOT / "tests/references/g2_plate_jets.py")
jets = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(jets)


def exact_point(curve: dict, t: Fraction) -> tuple[Fraction, Fraction]:
    h = jets.de_boor(int(curve["degree"]), [Fraction(k) for k in curve["knots"]],
                     [[Fraction(v) for v in p] for p in curve["controls"]], t, "right")
    return h[0] / h[2], h[1] / h[2]


def _homogeneous_exact(curve: dict, t: Fraction, side: str = "right") -> list[Fraction]:
    return jets.de_boor(int(curve["degree"]), [Fraction(k) for k in curve["knots"]],
                        [[Fraction(v) for v in p] for p in curve["controls"]], t, side)


def _solve_fraction(matrix: list[list[Fraction]], rhs: list[Fraction]) -> list[Fraction]:
    a = [row[:] + [value] for row, value in zip(matrix, rhs)]
    for col in range(len(a)):
        pivot = next(r for r in range(col, len(a)) if a[r][col])
        a[col], a[pivot] = a[pivot], a[col]
        q = a[col][col]
        a[col] = [v / q for v in a[col]]
        for r in range(len(a)):
            if r != col and a[r][col]:
                q = a[r][col]
                a[r] = [x - q * y for x, y in zip(a[r], a[col])]
    return [a[i][-1] for i in range(len(a))]


def span_level_polynomial(curve: dict, span: tuple[int, int], level: Fraction) -> list[int]:
    a, b = map(Fraction, span)
    nodes = [Fraction(0), Fraction(1, 3), Fraction(2, 3), Fraction(1)]
    matrix = [[x**i for i in range(4)] for x in nodes]
    values = []
    for s in nodes:
        t = a + (b - a) * s
        h = _homogeneous_exact(curve, t, "left" if s == 1 else "right")
        values.append(level * h[2] - h[1])
    rational = _solve_fraction(matrix, values)
    den = 1
    for c in rational:
        den = math.lcm(den, c.denominator)
    integers = [int(c * den) for c in rational]
    common = math.gcd(*[abs(c) for c in integers if c]) or 1
    integers = [c // common for c in integers]
    if integers[-1] < 0:
        integers = [-c for c in integers]
    return integers  # c0,c1,c2,c3 in the span-local parameter s


def rational_root_candidates(coefficients: list[int]) -> dict:
    c0, c1, c2, c3 = coefficients
    def divisors(n: int) -> list[int]:
        n = abs(n)
        return sorted({i for i in range(1, n + 1) if n % i == 0}) if n else [0]
    nums = divisors(c0)
    dens = divisors(c3)
    candidates = sorted({Fraction(sign * p, q) for p in nums for q in dens if q
                         for sign in (-1, 1)})
    roots = [r for r in candidates if c0 + c1 * r + c2 * r**2 + c3 * r**3 == 0]
    return {"integer_coefficients_c0_to_c3": coefficients,
            "rational_root_theorem_candidates_checked": len(candidates),
            "rational_roots": [str(r) for r in roots], "has_rational_root": bool(roots)}


def find_level_roots(curve: dict, level: float, subdivisions: int = 128, axis: int = 1) -> list[float]:
    degree, knots, controls = int(curve["degree"]), curve["knots"], curve["controls"]
    lo, hi = float(knots[degree]), float(knots[len(controls)])
    spans = sorted({(float(knots[i]), float(knots[i + 1]))
                    for i in range(degree, len(controls)) if knots[i] < knots[i + 1]})
    roots: list[float] = []
    def f(t: float) -> float:
        return jets.rational_jets(curve, t, "right")["position_mm"][axis] - level
    for a, b in spans:
        xs = [a + (b - a) * i / subdivisions for i in range(subdivisions + 1)]
        fs = [f(x) for x in xs]
        for i in range(subdivisions):
            left, right, fl, fr = xs[i], xs[i + 1], fs[i], fs[i + 1]
            if abs(fl) < 1e-12:
                roots.append(left)
            if fl * fr < 0:
                for _ in range(80):
                    mid = (left + right) * 0.5
                    fm = f(mid)
                    if fl * fm <= 0:
                        right, fr = mid, fm
                    else:
                        left, fl = mid, fm
                roots.append((left + right) * 0.5)
        if abs(fs[-1]) < 1e-12:
            roots.append(b)
    roots.sort()
    unique = []
    for t in roots:
        if lo <= t <= hi and (not unique or abs(t - unique[-1]) > 1e-8):
            unique.append(t)
    return unique


def line_locus(curve: dict, t: float, zrange: list[float], branch: str) -> dict:
    j = jets.rational_jets(curve, t, "right")
    x, y = j["position_mm"]
    a, b = zrange
    endpoints = [[x, y, a], [x, y, b]]
    tx, ty = j["unit_tangent_oriented"]
    profile_normal = [ty, -tx, 0.0]
    rect_side_normal = [0.0, 1.0, 0.0]
    transversality = math.sqrt(max(0.0, 1.0 - sum(profile_normal[i] * rect_side_normal[i]
                                                 for i in range(3)) ** 2))
    return {"branch": branch, "parameter": t, "source_position_xy_mm": [x, y],
            "unit_tangent_oriented_xy": j["unit_tangent_oriented"],
            "curvature_vector_xy_per_mm": j["curvature_vector_per_mm"],
            "sketch_provenance": ["free:boundary", "cross:boundary"],
            "intersection_classification": "transverse" if transversality > 1e-10 else "indeterminate_tangent",
            "surface_normal_sine": transversality,
            "3d_endpoints_mm": endpoints,
            "generatrix_tangent_xyz": [0.0, 0.0, 1.0],
            "endpoint_sketch_junction_classification": "G0_sharp_rectangle_corner; no G2 claim",
            "projections": {"xy": [[x, y], [x, y]], "xz": [[x, a], [x, b]],
                            "yz": [[y, a], [y, b]]},
            "length_mm": b - a}


def arc_locus(curve: dict, ta: float, tb: float, z: float, branch: str, count: int = 9) -> dict:
    points = []
    for i in range(count):
        t = ta + (tb - ta) * i / (count - 1)
        x, y = jets.rational_jets(curve, t, "right")["position_mm"]
        points.append([x, y, z])
    import numpy as np
    q, w = np.polynomial.legendre.leggauss(32)
    length = sum(weight * (tb - ta) * 0.5
                 * jets.rational_jets(curve, (ta + tb) * 0.5 + (tb - ta) * x * 0.5,
                                      "right")["speed_mm_per_parameter"]
                 for x, weight in zip(q, w))
    middle = jets.rational_jets(curve, (ta + tb) * 0.5, "right")
    tx, ty = middle["unit_tangent_oriented"]
    return {"branch": branch, "parameter_interval": [ta, tb], "z_mm": z,
            "sampled_curve_3d_mm": points,
            "sketch_provenance": ["free:boundary", "cross:boundary"],
            "intersection_classification": "transverse_surface_intersection",
            "surface_normal_sine": 1.0,
            "profile_surface_normal_midpoint": [ty, -tx, 0.0],
            "rectangle_face_normal": [0.0, 0.0, 1.0],
            "profile_arc_length_mm_q32": length,
            "endpoint_sketch_junction_classification": "G0_sharp_rectangle_corner; no G2 claim",
            "projections": {"xy": [[p[0], p[1]] for p in points],
                            "xz": [[p[0], p[2]] for p in points],
                            "yz": [[p[1], p[2]] for p in points]}}


def insert_knot(curve: dict, knot: Fraction, integerize: bool = True) -> dict:
    p = int(curve["degree"])
    U = [Fraction(k) for k in curve["knots"]]
    P = [[Fraction(v) for v in row] for row in curve["controls"]]
    n = len(P) - 1
    k = next(i for i in range(p, n + 1) if U[i] <= knot < U[i + 1])
    s = sum(1 for u in U if u == knot)
    Q = [None] * (n + 2)
    for i in range(0, k - p + 1):
        Q[i] = P[i]
    for i in range(k - s, n + 1):
        Q[i + 1] = P[i]
    for i in range(k - p + 1, k - s + 1):
        den = U[i + p] - U[i]
        alpha = (knot - U[i]) / den if den else Fraction(0)
        Q[i] = [alpha * P[i][c] + (1 - alpha) * P[i - 1][c]
                for c in range(len(P[i]))]
    V = U[:k + 1] + [knot] + U[k + 1:]
    knot_scale = 1
    for v in V:
        knot_scale = math.lcm(knot_scale, v.denominator)
    control_scale = 1
    for row in Q:
        for v in row:
            control_scale = math.lcm(control_scale, v.denominator)
    if integerize:
        out_knots = [int(v * knot_scale) for v in V]
        out_controls = [[int(v * control_scale) for v in row] for row in Q]
    else:
        out_knots = [float(v) for v in V]
        out_controls = [[float(v) for v in row] for row in Q]
    return {"degree": p, "knots": out_knots, "controls": out_controls}


def reindexed_curve(curve: dict, shift: int = 1) -> dict:
    p = int(curve["degree"])
    core = curve["controls"][:-p]
    core = core[shift:] + core[:shift]
    return {**curve, "controls": core + core[:p]}


def reversed_curve(curve: dict) -> dict:
    U = curve["knots"]
    a, b = U[int(curve["degree"])], U[len(curve["controls"])]
    return {"degree": curve["degree"],
            "knots": [a + b - u for u in reversed(U)],
            "controls": list(reversed(curve["controls"]))}


def transformation_checks(curve: dict, levels: list[float]) -> dict:
    degree = int(curve["degree"])
    lo, hi = float(curve["knots"][degree]), float(curve["knots"][len(curve["controls"])])
    sample = [lo + (hi - lo) * i / 64 for i in range(65)]
    base = [jets.rational_jets(curve, t, "right")["position_mm"] for t in sample]
    inserted = insert_knot(curve, Fraction(9, 2), integerize=False)
    inserted_error = max(math.dist(a, jets.rational_jets(inserted, t, "right")["position_mm"])
                         for t, a in zip(sample, base))
    inserted_integer = insert_knot(curve, Fraction(9, 2), integerize=True)
    parameter_scale = (float(inserted_integer["knots"][len(inserted_integer["controls"])])
                       - float(inserted_integer["knots"][int(inserted_integer["degree"])])
                       ) / (hi - lo)
    inserted_integer_error = max(
        math.dist(a, jets.rational_jets(inserted_integer, t * parameter_scale, "right")["position_mm"])
        for t, a in zip(sample, base))
    reverse = reversed_curve(curve)
    reverse_error = max(math.dist(a, jets.rational_jets(reverse, lo + hi - t, "right")["position_mm"])
                        for t, a in zip(sample, base))
    reindexed = reindexed_curve(curve, 1)
    period = hi - lo
    reindex_error = max(math.dist(a, jets.rational_jets(reindexed, lo + ((t - lo - 1) % period), "right")["position_mm"])
                        for t, a in zip(sample, base))
    roots = {str(level): find_level_roots(curve, level) for level in levels}
    roots_after_insert = {str(level): find_level_roots(inserted, level) for level in levels}
    roots_after_integer_insert = {
        str(level): [t / parameter_scale for t in find_level_roots(inserted_integer, level)]
        for level in levels}
    roots_after_reverse = {str(level): find_level_roots(reverse, level) for level in levels}
    roots_after_reindex = {str(level): find_level_roots(reindexed, level) for level in levels}
    return {"sample_count": len(sample), "inserted_knot": 4.5,
            "max_position_error_after_knot_insertion_mm": inserted_error,
            "max_position_error_after_integer_recipe_knot_insertion_mm": inserted_integer_error,
            "integer_recipe_knot_parameter_scale": parameter_scale,
            "max_position_error_after_parameter_reversal_mm": reverse_error,
            "max_position_error_after_periodic_control_reindex_mm": reindex_error,
            "intersection_parameters_by_y_after_insertion": roots_after_insert,
            "intersection_parameters_by_y_after_integer_recipe_insertion": roots_after_integer_insert,
            "intersection_parameters_by_y_after_reversal": roots_after_reverse,
            "intersection_parameters_by_y_after_periodic_reindex": roots_after_reindex,
            "invariance_threshold_mm": 1e-10,
            "invariant": max(inserted_error, inserted_integer_error, reverse_error, reindex_error) <= 1e-10,
            "roots": roots}


def sketch_curve(sketch: dict) -> dict:
    loops = sketch.get("loops", [])
    if len(loops) != 1:
        raise ValueError(f"expected one closed loop for sketch {sketch.get('id')}")
    return loops[0]["curve"]


def world_point(sketch: dict, local: tuple[float, float]) -> list[float]:
    frame = sketch["frame"]
    o, u, v = frame["origin_mm"], frame["u"], frame["v"]
    return [o[i] + local[0] * u[i] + local[1] * v[i] for i in range(3)]


def _weighted_integer_controls(points: list[tuple[Fraction, Fraction]]) -> list[list[int]]:
    import math as _math
    den = 1
    for p in points:
        for x in p:
            den = _math.lcm(den, x.denominator)
    return [[int(x * den), int(y * den), den] for x, y in points]


def elevated_rectangle(rectangle_curve: dict, split_first_side: bool = False) -> dict:
    degree, knots, controls = int(rectangle_curve["degree"]), rectangle_curve["knots"], rectangle_curve["controls"]
    if degree != 1:
        raise ValueError("base rectangle fixture must be degree one")
    vertices = [(Fraction(p[0], p[2]), Fraction(p[1], p[2])) for p in controls]
    segments = []
    for side, (a, b) in enumerate(zip(vertices, vertices[1:])):
        endpoints = [a, b]
        if split_first_side and side == 0:
            endpoints.insert(1, tuple((a[k] + b[k]) / 2 for k in range(2)))
        segments.extend(zip(endpoints, endpoints[1:]))
    bspline_controls = [segments[0][0]]
    for a, b in segments:
        bspline_controls.extend([
            tuple((2 * a[k] + b[k]) / 3 for k in range(2)),
            tuple((a[k] + 2 * b[k]) / 3 for k in range(2)),
            b,
        ])
    internal = []
    for side in range(len(vertices) - 1):
        if split_first_side and side == 0:
            internal.append(Fraction(1, 2))
        internal.append(Fraction(side + 1))
    # The last endpoint is the closure and is clamped with multiplicity four.
    knot_values = [Fraction(0)] * 4
    for k in internal[:-1]:
        knot_values.extend([k] * 3)
    knot_values.extend([Fraction(4)] * 4)
    # The rational side split is at 1/2. Scale the whole parameter domain by
    # two so every knot remains an exact integer in the versioned recipe.
    knot_scale = 2 if split_first_side else 1
    return {"degree": 3, "knots": [int(k * knot_scale) for k in knot_values],
            "controls": _weighted_integer_controls(bspline_controls)}


def max_reparameterized_trace_error(reference: dict, candidate: dict,
                                    samples: int = 128) -> float:
    ref_a, ref_b = reference["knots"][int(reference["degree"])], reference["knots"][len(reference["controls"])]
    cand_a = candidate["knots"][int(candidate["degree"]) ]
    cand_b = candidate["knots"][len(candidate["controls"])]
    errors = []
    for i in range(samples + 1):
        fraction = i / samples
        ref_t = ref_a + (ref_b - ref_a) * fraction
        cand_t = cand_a + (cand_b - cand_a) * fraction
        ref_side = "right" if i < samples else "left"
        cand_side = "right" if i < samples else "left"
        p_ref = jets.rational_jets(reference, ref_t, ref_side)["position_mm"]
        p_cand = jets.rational_jets(candidate, cand_t, cand_side)["position_mm"]
        errors.append(math.dist(p_ref, p_cand))
    return max(errors, default=0.0)


def rectangle_variant_diagnostics(curve: dict) -> dict:
    degree, knots = int(curve["degree"]), curve["knots"]
    internal = sorted({k for k in knots if knots[degree] < k < knots[len(curve["controls"])]})
    reports = [jets.compare_jets(curve, t, "left", "right", 1e-10, 1e-10, 1e-10)
               for t in internal]
    for row in reports:
        left = row["left"]["unit_tangent_oriented"]
        right = row["right"]["unit_tangent_oriented"]
        left_normal, right_normal = [0.0, left[1], -left[0]], [0.0, right[1], -right[0]]
        row["extruded_surface"] = {
            "extrusion_axis_world": [1.0, 0.0, 0.0],
            "oriented_unit_normal_left": left_normal,
            "oriented_unit_normal_right": right_normal,
            "oriented_unit_normal_jump": math.dist(left_normal, right_normal),
            "section_curvature_jump_per_mm": row["curvature_vector_jump_per_mm"],
            "generator_curvature_left_per_mm": 0.0,
            "generator_curvature_right_per_mm": 0.0,
            "generator_curvature_jump_per_mm": 0.0,
        }
    closure = jets.rational_jets(curve, knots[len(curve["controls"])], "left")
    start = jets.rational_jets(curve, knots[degree], "right")
    closure_position = math.dist(closure["position_mm"], start["position_mm"])
    closure_tangent = math.dist(closure["unit_tangent_oriented"], start["unit_tangent_oriented"])
    closure_normal_jump = math.dist([0.0, closure["unit_tangent_oriented"][1],
                                     -closure["unit_tangent_oriented"][0]],
                                    [0.0, start["unit_tangent_oriented"][1],
                                     -start["unit_tangent_oriented"][0]])
    corners = sum(row["oriented_unit_tangent_jump"] > 1e-10 for row in reports)
    closure_corner = closure_tangent > 1e-10
    return {"degree": degree,
            "knot_multiplicities": {str(k): knots.count(k) for k in internal},
            "internal_junctions": [{**row, "classification": "G0_sharp"
                                     if row["oriented_unit_tangent_jump"] > 1e-10
                                     else row["classification"]} for row in reports],
            "closure": {"position_jump_mm": closure_position,
                        "oriented_unit_tangent_jump": closure_tangent,
                        "extruded_surface_normal_jump": closure_normal_jump,
                        "generator_curvature_jump_per_mm": 0.0,
                        "classification": "G0_sharp" if closure_corner else "smooth"},
            "true_sharp_corners": corners + int(closure_corner),
            "smooth_internal_knots": sum(row["oriented_unit_tangent_jump"] <= 1e-10 for row in reports)}


def sketch_geometry(recipe: dict):
    if "sketches" not in recipe:
        profile = recipe["profile"]
        ybounds = [Fraction(*p) if isinstance(p, list) else Fraction(p)
                   for p in recipe["rect_y_mm"]]
        zrange = [Fraction(z) for z in recipe["rect_z_mm"]]
        return profile, [float(y) for y in ybounds], [float(z) for z in zrange], None
    free = next(s for s in recipe["sketches"] if s["direction"] == [0, 0, 1])
    cross = next(s for s in recipe["sketches"] if s["direction"] == [1, 0, 0])
    profile, rect = sketch_curve(free), sketch_curve(cross)
    frame = cross["frame"] if "frame" in cross else cross
    # Sketch frames are orthonormal coordinate maps; the rectangle's local u
    # and v axes are world y and z for this prospect.
    cf = cross["frame"]
    local = [(Fraction(p[0], p[2]), Fraction(p[1], p[2])) for p in rect["controls"]]
    world = [tuple(Fraction(cf["origin_mm"][k])
                   + uv[0] * Fraction(cf["u"][k]) + uv[1] * Fraction(cf["v"][k])
                   for k in range(3)) for uv in local]
    ybounds = [min(p[1] for p in world), max(p[1] for p in world)]
    zrange = [min(p[2] for p in world), max(p[2] for p in world)]
    return profile, [float(y) for y in ybounds], [float(z) for z in zrange], rect


def exact_exploratory_levels(profile: dict, levels: list[int]) -> dict:
    p, knots, controls = int(profile["degree"]), profile["knots"], profile["controls"]
    spans = [(knots[i], knots[i + 1]) for i in range(p, len(controls))
             if knots[i] < knots[i + 1]]
    report = {}
    for level_value in levels:
        roots = find_level_roots(profile, float(level_value))
        rows = []
        all_rational_rootless = True
        for t in roots:
            span = next((ab for ab in spans if ab[0] - 1e-9 <= t <= ab[1] + 1e-9), None)
            if span is None:
                continue
            poly = rational_root_candidates(span_level_polynomial(profile, span, Fraction(level_value)))
            rows.append({"span": span, "root_parameter_numeric": t,
                         "local_parameter_numeric": t - span[0], **poly})
            all_rational_rootless &= not poly["has_rational_root"]
        transversality = [abs(jets.rational_jets(profile, t, "right")["first_derivative"][1])
                          for t in roots]
        exists = len(roots) >= 2 and all(v > 1e-8 for v in transversality)
        report[str(level_value)] = {
            "roots_numeric": roots,
            "polynomial_by_root_span": rows,
            "transverse_level_crossings": exists,
            "no_rational_root_in_each_candidate_span": all_rational_rootless and bool(rows),
            "classification": "refused_exact_rational_parameter_domain_intersection_exists"
                if exists and all_rational_rootless else "admitted"
                if exists else "indeterminate",
            "interpretation": "The level set crosses the source profile transversely in real parameters; refusal denotes the exact rational-parameter representation limit, not absence of a geometric intersection. A B-rep approximation can represent this geometry approximately.",
            "rational_root_theorem": "all signed divisors of the constant term checked against reduced integer span polynomial",
        }
    return report


def clipped_profile_area(profile: dict, ylo: float, yhi: float, order: int) -> float:
    import numpy as np
    q, w = np.polynomial.legendre.leggauss(order)
    area = 0.0
    for x, weight in zip(q, w):
        y = (ylo + yhi) * 0.5 + (yhi - ylo) * x * 0.5
        roots = find_level_roots(profile, y, subdivisions=48)
        if len(roots) < 2:
            raise ArithmeticError(f"unresolved profile intersections at y={y}")
        xs = [jets.rational_jets(profile, t, "right")["position_mm"][0] for t in roots]
        area += weight * (yhi - ylo) * 0.5 * (max(xs) - min(xs))
    return area


def run_cli(cli: Path, command: list[str]) -> dict:
    process = subprocess.run([str(cli), *command], cwd=ROOT, capture_output=True, text=True)
    if process.returncode:
        raise RuntimeError(f"CLI {' '.join(command)} failed: {process.stderr or process.stdout}")
    return json.loads(process.stdout)


def model_volume_and_crossing_checks(model: dict, source_report: dict) -> dict:
    """Check model volume by layer and its exact two-sketch crossing plan."""
    spec = importlib.util.spec_from_file_location(
        "g2_cross_volume_oracle", ROOT / "tests/references/g2_cross_model_volume.py")
    volume_oracle = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(volume_oracle)
    by_layer = volume_oracle.integrate(model, 32)
    bounds = model["source_recipe"]["bounds_xyz_mm"]
    stock_area = ((bounds[0][1] - bounds[0][0]) * (bounds[1][1] - bounds[1][0]))
    profile_area = source_report["periodic_profile_removed_volume_mm3"] / (bounds[2][1] - bounds[2][0])
    lo, hi = [a / b for a, b in model["sketch_plan"]["z_cut_levels_mm"]]
    cross_area = ((source_report["rectangle_removed_volume_mm3"]
                   - source_report["intersection_volume_q64_mm3"]) / (hi - lo))
    layers = []
    for (z0, z1), value in sorted(by_layer.items()):
        removed_depth = max(0.0, min(hi, z1) - max(lo, z0))
        expected = (stock_area - profile_area) * (z1 - z0) - removed_depth * cross_area
        layers.append({"z_mm": [z0, z1], "volume_mm3": value,
                       "expected_volume_mm3": expected,
                       "relative_error": abs(value - expected) / max(abs(expected), 1.0)})
    total = sum(by_layer.values())
    expected_total = source_report["remaining_stock_volume_mm3"]
    volume = {"order": 32, "volume_mm3": total, "expected_volume_mm3": expected_total,
              "relative_total_error": abs(total - expected_total) / expected_total,
              "max_relative_layer_error": max((row["relative_error"] for row in layers), default=math.inf),
              "layers": layers}
    crossings = model.get("sketch_plan", {}).get("crossings", [])
    expected_y = {"719/36", "1081/36"}
    y_keys = [f"{row['y_mm'][0]}/{row['y_mm'][1]}" for row in crossings]
    crossing_types_ok = all(row.get("profile") == "free"
                            and row.get("profile_loop") == "boundary"
                            and row.get("transverse_sketch") == "cross"
                            and row.get("transverse_loop") == "boundary"
                            and len(row.get("parameter", [])) == 2
                            and len(row.get("xy_mm", [])) == 2 for row in crossings)
    crossing_check = {"count": len(crossings), "expected_count": 4,
                      "y_level_multiplicities": {key: y_keys.count(key) for key in sorted(expected_y)},
                      "all_are_free_profile_vs_cross_rectangle_boundary_events": crossing_types_ok,
                      "events": crossings,
                      "status": "verified" if len(crossings) == 4 and crossing_types_ok
                          and all(y_keys.count(key) == 2 for key in expected_y) else "failed"}
    volume["status"] = "verified" if volume["relative_total_error"] < 1e-8 \
        and volume["max_relative_layer_error"] < 1e-8 else "failed"
    return {"volume_q32": volume, "sketch_plan_crossings": crossing_check}


def c0_iga_check(cli: Path, recipe_path: Path, work: Path,
                 volume_reference: dict) -> dict:
    import hashlib
    spec = importlib.util.spec_from_file_location(
        "g2_cross_space_oracle", ROOT / "tests/references/nurbs_space_oracle.py")
    space = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(space)
    spec = importlib.util.spec_from_file_location(
        "g2_cross_patch_oracle", ROOT / "tests/references/nurbs_patch_equilibrium.py")
    patch = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(patch)
    model_path, package_path = work / "cross.model.json", work / "cross.space.json"
    build = run_cli(cli, ["build", str(recipe_path), str(model_path)])
    model = json.loads(model_path.read_text())
    if model.get("source_recipe") != json.loads(recipe_path.read_text()):
        raise AssertionError("C0 model does not preserve the exact two-sketch recipe")
    export = run_cli(cli, ["export-nurbs-space", str(model_path), str(package_path)])
    package = json.loads(package_path.read_text())
    trace = space.package_checks(model, package)
    equilibrium = [patch.integrate(package, q) for q in (4, 8, 16)]
    if equilibrium[-1]["relative_residual_linf"] > 1e-8:
        raise AssertionError(f"q16 constant-stress equilibrium exceeds threshold: {equilibrium[-1]}")
    features = Counter(face["feature"] for face in model["design_faces"])
    model_diagnostics = model_volume_and_crossing_checks(model, volume_reference)
    if model_diagnostics["volume_q32"]["status"] != "verified" \
            or model_diagnostics["sketch_plan_crossings"]["status"] != "verified":
        raise AssertionError(f"base model volume or crossing plan failed: {model_diagnostics}")
    wall_samples, wall_max_error = 0, 0.0
    source_recipe = json.loads(recipe_path.read_text())
    free = next(s for s in source_recipe["sketches"] if s["direction"] == [0, 0, 1])
    source_curve = sketch_curve(free)
    package_blocks = {b["id"]: b for b in package["blocks"]}
    for face in model["design_faces"]:
        if face["feature"] != "nurbs_z_hole/wall":
            continue
        axis, upper = "uvw".index(face["side"][0]), face["side"][1] == "1"
        free_axes = [a for a in range(3) if a != axis]
        block = package_blocks[face["block"]]
        for i in range(4):
            for j in range(4):
                uvw = [Fraction(0)] * 3
                uvw[axis] = Fraction(int(upper))
                uvw[free_axes[0]], uvw[free_axes[1]] = Fraction(i, 3), Fraction(j, 3)
                point = space.map_from_json(block, tuple(uvw))
                roots_x = find_level_roots(source_curve, float(point[0]), axis=0)
                if not roots_x:
                    wall_max_error = math.inf
                else:
                    error = min(abs(jets.rational_jets(source_curve, t, "right")["position_mm"][1]
                                    - float(point[1])) for t in roots_x)
                    wall_max_error = max(wall_max_error, error)
                wall_samples += 1
    wall_check = {"source_parameterization": "periodic B-spline from raw recipe controls",
                  "sample_count": wall_samples,
                  "max_xy_source_curve_distance_mm": wall_max_error,
                  "threshold_mm": 1e-10,
                  "status": "verified" if wall_samples and wall_max_error <= 1e-10
                      else "unavailable" if not wall_samples else "failed"}
    constructed_jets = constructed_wall_junctions(model, source_curve)
    return {"status": "verified", "cli": str(cli),
            "cli_sha256": hashlib.sha256(cli.read_bytes()).hexdigest(),
            "build_audit": build, "space_export": export,
            "design_face_feature_counts": dict(sorted(features.items())),
            "periodic_wall_vs_source_curve": wall_check,
            "constructed_wall_one_sided_jets": constructed_jets,
            **model_diagnostics,
            "exact_c0_trace_and_affine_checks": trace,
            "constant_stress_equilibrium_orders": equilibrium,
            "q16_relative_linf_threshold": 1e-8,
            "source_recipe_preserved": True,
            "artifacts": {"model": str(model_path.relative_to(ROOT)),
                          "q3_package": str(package_path.relative_to(ROOT))}}


def nearest_source_parameter(curve: dict, target: tuple[float, float], interval: tuple[float, float]) -> float:
    lo, hi = interval
    count = 32
    samples = [lo + (hi - lo) * i / count for i in range(count + 1)]
    def distance2(t: float) -> float:
        point = jets.rational_jets(curve, t, "right")["position_mm"]
        return (point[0] - target[0]) ** 2 + (point[1] - target[1]) ** 2
    index = min(range(len(samples)), key=lambda i: distance2(samples[i]))
    a, b = samples[max(0, index - 1)], samples[min(count, index + 1)]
    ratio = (math.sqrt(5.0) - 1.0) * 0.5
    x1, x2 = b - ratio * (b - a), a + ratio * (b - a)
    f1, f2 = distance2(x1), distance2(x2)
    for _ in range(60):
        if f1 > f2:
            a, x1, f1 = x1, x2, f2
            x2 = a + ratio * (b - a)
            f2 = distance2(x2)
        else:
            b, x2, f2 = x2, x1, f1
            x1 = b - ratio * (b - a)
            f1 = distance2(x1)
    return (a + b) * 0.5


def constructed_wall_junctions(model: dict, source_curve: dict) -> dict:
    blocks = {block["id"]: block for block in model["blocks"]}
    selected = {}
    for face in model["design_faces"]:
        if face["feature"] != "nurbs_z_hole/wall":
            continue
        block_id = face["block"]
        match = re.fullmatch(r"g2_s(\d+)_p(\d+)_z(\d+)_src(\d+)", block_id)
        if not match:
            continue
        sector, piece, layer, source_span = map(int, match.groups())
        patch_key = (sector, piece, source_span)
        if patch_key in selected:
            continue
        block = blocks[block_id]
        # The wall is the v=1 face. At w=0 the homogeneous control row is an
        # independent cubic rational curve in the source sketch plane.
        controls = [block["controls"][(i * 2 + 1) * 2] for i in range(4)]
        curve = {"degree": 3, "knots": [0, 0, 0, 0, 1, 1, 1, 1],
                 "controls": [[p[0], p[1], p[3]] for p in controls]}
        xy0 = jets.rational_jets(curve, 0.0, "right")["position_mm"]
        xy1 = jets.rational_jets(curve, 1.0, "left")["position_mm"]
        t0 = nearest_source_parameter(source_curve, tuple(xy0), (float(source_span), float(source_span + 1)))
        t1 = nearest_source_parameter(source_curve, tuple(xy1), (float(source_span), float(source_span + 1)))
        if t1 < t0:
            curve["controls"] = list(reversed(curve["controls"]))
            xy0, xy1 = xy1, xy0
            t0, t1 = t1, t0
        sample_error = 0.0
        for u in (0.0, 1 / 3, 2 / 3, 1.0):
            point = jets.rational_jets(curve, u, "right")["position_mm"]
            t = nearest_source_parameter(source_curve, tuple(point),
                                         (float(source_span), float(source_span + 1)))
            source_point = jets.rational_jets(source_curve, t, "right")["position_mm"]
            sample_error = max(sample_error, math.dist(point, source_point))
        selected[patch_key] = {"sector": sector, "piece": piece, "source_span": source_span,
                              "curve": curve, "parameter_interval": [t0, t1],
                              "max_sample_error_mm": sample_error}

    ordered = sorted(selected.values(), key=lambda row: (row["sector"], row["piece"]))
    junctions = []
    if len(ordered) >= 2:
        for i, left_row in enumerate(ordered):
            right_row = ordered[(i + 1) % len(ordered)]
            left = jets.rational_jets(left_row["curve"], 1.0, "left")
            right = jets.rational_jets(right_row["curve"], 0.0, "right")
            pos = math.dist(left["position_mm"], right["position_mm"])
            tangent = math.dist(left["unit_tangent_oriented"], right["unit_tangent_oriented"])
            curvature = math.dist(left["curvature_vector_per_mm"], right["curvature_vector_per_mm"])
            t = left_row["parameter_interval"][1]
            parameter = t % 8.0
            if abs(parameter - round(parameter)) < 1e-7:
                kind = "source_periodic_knot"
            elif any(abs(parameter - v) < 1e-7 for v in (2.5, 3.5, 6.5, 7.5)):
                kind = "rectangle_intersection_cut"
            else:
                kind = "sector_partition"
            n_left = [left["unit_tangent_oriented"][1], -left["unit_tangent_oriented"][0], 0.0]
            n_right = [right["unit_tangent_oriented"][1], -right["unit_tangent_oriented"][0], 0.0]
            junctions.append({"after_sector": left_row["sector"], "before_sector": right_row["sector"],
                              "kind": kind, "source_parameter": t,
                              "position_jump_mm": pos,
                              "oriented_unit_tangent_jump": tangent,
                              "extruded_surface_normal_jump": math.dist(n_left, n_right),
                              "curvature_vector_jump_per_mm": curvature,
                              "classification": "G2" if max(pos, tangent, curvature) <= 1e-9
                                  else "G1_only" if max(pos, tangent) <= 1e-9
                                  else "C0_only" if pos <= 1e-9 else "discontinuous"})
    return {"face_patch_count": len(ordered),
            "pieces": [{k: row[k] for k in ("sector", "piece", "source_span", "parameter_interval", "max_sample_error_mm")}
                       for row in ordered],
            "junction_count": len(junctions), "junctions": junctions,
            "source_position_tangent_curvature_thresholds": {"position_mm": 1e-9,
                                                               "unit_tangent": 1e-9,
                                                               "curvature_mm_inverse": 1e-9},
            "all_constructed_profile_joins_g2": bool(junctions) and all(j["classification"] == "G2" for j in junctions),
            "maximum_span_sample_error_mm": max((row["max_sample_error_mm"] for row in ordered), default=None)}


def recipe_variants(recipe: dict, rectangle_curve: dict) -> dict[str, dict]:
    variants = {}
    free_sketch = next(s for s in recipe["sketches"] if s["direction"] == [0, 0, 1])
    base_profile = sketch_curve(free_sketch)
    elevated = copy.deepcopy(recipe)
    next(s for s in elevated["sketches"] if s["direction"] == [1, 0, 0])["loops"][0]["curve"] = elevated_rectangle(rectangle_curve)
    variants["rectangle_degree3_elevated"] = elevated
    subdivided = copy.deepcopy(recipe)
    next(s for s in subdivided["sketches"] if s["direction"] == [1, 0, 0])["loops"][0]["curve"] = elevated_rectangle(rectangle_curve, True)
    variants["rectangle_degree3_side_knot"] = subdivided
    inserted = copy.deepcopy(recipe)
    next(s for s in inserted["sketches"] if s["direction"] == [0, 0, 1])["loops"][0]["curve"] = insert_knot(base_profile, Fraction(9, 2))
    variants["periodic_profile_knot_inserted"] = inserted
    reversed_input = copy.deepcopy(recipe)
    next(s for s in reversed_input["sketches"] if s["direction"] == [0, 0, 1])["loops"][0]["curve"] = reversed_curve(base_profile)
    variants["periodic_profile_reversed"] = reversed_input
    reindexed = copy.deepcopy(recipe)
    next(s for s in reindexed["sketches"] if s["direction"] == [0, 0, 1])["loops"][0]["curve"] = reindexed_curve(base_profile, 1)
    variants["periodic_profile_reindexed"] = reindexed
    return variants


def cli_variant_checks(cli: Path, recipe: dict, rectangle_curve: dict, work: Path,
                       volume_reference: dict) -> dict:
    results = {}
    for name, variant in recipe_variants(recipe, rectangle_curve).items():
        rp, mp = work / f"{name}.recipe.json", work / f"{name}.model.json"
        rp.write_text(json.dumps(variant, separators=(",", ":")) + "\n")
        try:
            report = run_cli(cli, ["build", str(rp), str(mp)])
            model = json.loads(mp.read_text())
            diagnostics = model_volume_and_crossing_checks(model, volume_reference)
            checks_verified = all(row["status"] == "verified" for row in diagnostics.values())
            results[name] = {"status": "accepted" if checks_verified else "accepted_validation_failed",
                             "build_audit": report,
                             "source_recipe_preserved": model.get("source_recipe") == variant,
                             "blocks": len(model["blocks"]),
                             "design_faces": len(model["design_faces"]),
                             **diagnostics}
        except Exception as error:
            results[name] = {"status": "rejected_or_unavailable", "detail": str(error)}
    return results


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("recipe", type=Path)
    ap.add_argument("--cli", type=Path,
                    help="optional production CLI for exact C0 Q3 traces, affine field, equilibrium, and variants")
    ap.add_argument("--work", type=Path, default=ROOT / "build/g2-cross-junctions")
    ap.add_argument("--output", type=Path, default=ROOT / "docs/results/g2-cross-junctions.json")
    args = ap.parse_args()
    recipe = json.loads(args.recipe.read_text())
    if "sketches" in recipe:
        sketches = recipe["sketches"]
        free = next(s for s in sketches if s["direction"] == [0, 0, 1])
        cross = next(s for s in sketches if s["direction"] == [1, 0, 0])
        profile = sketch_curve(free)
        rectangle_curve = sketch_curve(cross)
        rect_points = [world_point(cross, (p[0] / p[2], p[1] / p[2]))
                       for p in rectangle_curve["controls"]]
        profile, ylevels, zrange, rectangle_curve = sketch_geometry(recipe)
        rectangle_diagnostics = rectangle_variant_diagnostics(rectangle_curve)
        elevated = elevated_rectangle(rectangle_curve)
        subdivided = elevated_rectangle(rectangle_curve, split_first_side=True)
        elev_error = max_reparameterized_trace_error(rectangle_curve, elevated)
        subdiv_error = max_reparameterized_trace_error(rectangle_curve, subdivided)
        rectangle_variants = {
            "degree_one": {**rectangle_diagnostics, "max_position_error_to_degree_one_mm": 0.0},
            "degree_three_elevated": {**rectangle_variant_diagnostics(elevated),
                                      "max_position_error_to_degree_one_mm": elev_error},
            "degree_three_with_side_knot": {**rectangle_variant_diagnostics(subdivided),
                                             "max_position_error_to_degree_one_mm": subdiv_error,
                                             "parameter_domain_scale_from_degree_one": 2},
        }
    else:
        profile = recipe["profile"]
        ybounds = recipe["rect_y_mm"]
        ylevels = [float(Fraction(*v)) if isinstance(v, list) else float(v) for v in ybounds]
        zrange = recipe["rect_z_mm"]
        if isinstance(zrange[0], list):
            zrange = [float(Fraction(*v)) for v in zrange]
        else:
            zrange = list(map(float, zrange))
        rectangle_variants = None
    roots = {str(y): find_level_roots(profile, y) for y in ylevels}
    exact_intersections = []
    if "sketches" in recipe:
        cross_sketch = next(s for s in recipe["sketches"] if s["direction"] == [1, 0, 0])
        rectangle_points_exact = [tuple(Fraction(p[k], p[2]) for k in (0, 1))
                                  for p in sketch_curve(cross_sketch)["controls"]]
        exact_y_levels = [min(p[0] for p in rectangle_points_exact),
                          max(p[0] for p in rectangle_points_exact)]
    else:
        exact_y_levels = [Fraction(*p) if isinstance(p, list) else Fraction(p)
                          for p in recipe["rect_y_mm"]]
    for level_exact in exact_y_levels:
        for t in roots[str(float(level_exact))]:
            candidate = Fraction(round(t * 2), 2)
            point_exact = exact_point(profile, candidate)
            exact_intersections.append({"y_mm": str(level_exact), "parameter": str(candidate),
                                        "exact_source_point_xy_mm": [str(v) for v in point_exact],
                                        "y_matches_rectangle_level_exactly": point_exact[1] == level_exact})
    seam_lines = []
    for y in ylevels:
        for t in roots[str(y)]:
            seam_lines.append(line_locus(profile, t, zrange, f"y={y:g}"))
    ordered = sorted((t, y) for y in ylevels for t in roots[str(y)])
    seam_arcs = []
    for z in zrange:
        for (ta, ya), (tb, yb) in zip(ordered, ordered[1:]):
            if ya != yb:
                mid = (ta + tb) * 0.5
                ym = jets.rational_jets(profile, mid, "right")["position_mm"][1]
                if min(ylevels) < ym < max(ylevels):
                    seam_arcs.append(arc_locus(profile, ta, tb, z, f"arc@z={z:g}"))
    exploratory = exact_exploratory_levels(profile, [20, 30])
    transformed = transformation_checks(profile, ylevels)
    profile_jets = jets.analyze_curve(
        {"id": "periodic-profile", "curve": profile, "closed": True,
         "expect": "G2", "test_nonconic": True},
        {"position_tol": 1e-10, "tangent_tol": 1e-10, "curvature_tol": 1e-10})
    volume_report = None
    if "bounds_xyz_mm" in recipe:
        (xmin, xmax), (ymin, ymax), (zmin, zmax) = recipe["bounds_xyz_mm"]
        (ylo, yhi) = ylevels
        (zlo, zhi) = zrange
        profile_area = jets.periodic_area(profile, 64)
        intersection_area_32 = clipped_profile_area(profile, ylo, yhi, 32)
        intersection_area_64 = clipped_profile_area(profile, ylo, yhi, 64)
        vz = profile_area * (zmax - zmin)
        vx = (xmax - xmin) * (yhi - ylo) * (zhi - zlo)
        overlap = intersection_area_64 * (zhi - zlo)
        volume_report = {"periodic_profile_removed_volume_mm3": vz,
                         "rectangle_removed_volume_mm3": vx,
                         "intersection_volume_q32_mm3": intersection_area_32 * (zhi - zlo),
                         "intersection_volume_q64_mm3": overlap,
                         "intersection_relative_order_delta": abs(intersection_area_64 - intersection_area_32) / max(intersection_area_64, 1.0),
                         "union_removed_volume_mm3": vz + vx - overlap,
                         "remaining_stock_volume_mm3": (xmax - xmin) * (ymax - ymin) * (zmax - zmin) - vz - vx + overlap,
                         "method": "periodic B-spline Green integral and independent horizontal-slice integration of source-profile intersections"}
        coordinator_reference = {
            "periodic_profile_removed_volume_mm3": 3489.731852993521,
            "rectangle_removed_volume_mm3": 4022.222222222222,
            "intersection_volume_mm3": 712.3493494134791,
            "union_removed_volume_mm3": 6799.604725802264,
            "remaining_stock_volume_mm3": 53200.395274197734,
            "provenance": "independent SciPy coordinate reference supplied by the parent coordinator",
        }
        volume_report["coordinator_reference"] = coordinator_reference
        volume_report["relative_difference_q64"] = {
            key: abs(volume_report[field] - coordinator_reference[key]) / max(abs(coordinator_reference[key]), 1.0)
            for field, key in (("periodic_profile_removed_volume_mm3", "periodic_profile_removed_volume_mm3"),
                               ("rectangle_removed_volume_mm3", "rectangle_removed_volume_mm3"),
                               ("intersection_volume_q64_mm3", "intersection_volume_mm3"),
                               ("union_removed_volume_mm3", "union_removed_volume_mm3"),
                               ("remaining_stock_volume_mm3", "remaining_stock_volume_mm3"))}
    checks_passed = (
        transformed["invariant"]
        and bool(exact_intersections)
        and all(row["y_matches_rectangle_level_exactly"] for row in exact_intersections)
        and len(seam_lines) == 4
        and all(row["intersection_classification"] == "transverse" for row in seam_lines)
        and len(seam_arcs) == 4
        and all(row["profile_arc_length_mm_q32"] > 0.0 for row in seam_arcs)
        and (rectangle_variants is None or all(row["true_sharp_corners"] == 4
                                                for row in rectangle_variants.values()))
        and (volume_report is None or volume_report["intersection_relative_order_delta"] < 1e-8
             and all(v < 1e-8 for v in volume_report["relative_difference_q64"].values()))
    )
    report = {"schema": "g2_cross_seam_junctions_v1", "status": "verified" if checks_passed else "failed",
              "method": "independent periodic NURBS source evaluation; spanwise level-set bisection; exact homogeneous knot insertion; reverse and cyclic reindex metamorphic checks",
              "recipe": str(args.recipe), "rectangle_y_mm": ylevels, "rectangle_z_mm": zrange,
              "source_junctions_and_closure": profile_jets,
              "intersection_roots_by_y": roots,
              "exact_intersection_witnesses": exact_intersections,
              "rectangle_loop_variants": rectangle_variants,
              "vertical_generatrix_seam_loci": seam_lines,
              "horizontal_profile_arc_seam_loci": seam_arcs,
              "exploratory_levels": exploratory,
              "independent_volume_reference": volume_report,
              "metamorphic_curve_checks": transformed,
              "limits": ["Only transverse intersections of the supplied source periodic curve with the two rectangle y-levels are certified by this oracle.",
                         "Seam loci are source-curve constructions; they do not prove that the production Boolean or stitched surface complex is valid.",
                         "Exploratory level classifications are root-isolation diagnostics, not a general intersection certificate."]}
    if args.cli:
        try:
            work = args.work.resolve()
            work.mkdir(parents=True, exist_ok=True)
            report["c0_iga"] = c0_iga_check(args.cli.resolve(), args.recipe.resolve(), work,
                                             volume_report)
            report["production_variant_builds"] = cli_variant_checks(
                args.cli.resolve(), recipe, rectangle_curve, work, volume_report)
            if (report["c0_iga"]["status"] != "verified"
                    or any(row["status"] != "accepted"
                           for row in report["production_variant_builds"].values())):
                report["status"] = "failed"
        except Exception as error:
            report["c0_iga"] = {"status": "failed", "error": str(error)}
            report["status"] = "failed"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "seam_generatrices": len(seam_lines),
                      "seam_arcs": len(seam_arcs), "output": str(args.output)}))
    return 0 if report["status"] == "verified" else 1


if __name__ == "__main__":
    raise SystemExit(main())
