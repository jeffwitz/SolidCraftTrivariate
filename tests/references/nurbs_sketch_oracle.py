#!/usr/bin/env python3
"""Independent De Boor, quadrature and geometry oracle for cubic sketch strips."""
from __future__ import annotations

try:
    import FreeCAD as App
    import Part
except ImportError:
    App = None
    Part = None

import argparse
import hashlib
import json
import math
import platform
import random
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CORPUS = Path(__file__).with_name("nurbs_sketch_corpus.json")
SEED = 20261125


def basis(i, p, u, knots, memo):
    key = (i, p)
    if key in memo:
        return memo[key]
    if p == 0:
        result = 1.0 if knots[i] <= u < knots[i + 1] else 0.0
    else:
        left_den = knots[i + p] - knots[i]
        right_den = knots[i + p + 1] - knots[i + 1]
        left = 0.0 if left_den == 0 else (u - knots[i]) / left_den * basis(i, p - 1, u, knots, memo)
        right = 0.0 if right_den == 0 else (knots[i + p + 1] - u) / right_den * basis(i + 1, p - 1, u, knots, memo)
        result = left + right
    memo[key] = result
    return result


def curve(curve_spec, u, derivative=False):
    p = curve_spec["degree"]
    knots = curve_spec["knots"]
    controls = np.asarray(curve_spec["controls"], dtype=np.float64)
    end = knots[len(controls)]
    if u == end:
        u = np.nextafter(float(end), -math.inf)
    memo = {}
    n = len(controls)
    N = np.array([basis(i, p, u, knots, memo) for i in range(n)])
    dN = np.zeros(n)
    if p:
        low = [basis(i, p - 1, u, knots, memo) for i in range(n + 1)]
        for i in range(n):
            a = knots[i + p] - knots[i]
            b = knots[i + p + 1] - knots[i + 1]
            dN[i] = (0.0 if a == 0 else p * low[i] / a) - (0.0 if b == 0 else p * low[i + 1] / b)
    H = N @ controls
    Hd = dN @ controls
    point = H[:2] / H[2]
    if not derivative:
        return point
    return point, (Hd[:2] * H[2] - H[:2] * Hd[2]) / H[2]**2


def gauss_integral(function, a, b, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    x = (a + b) * 0.5 + (b - a) * 0.5 * nodes
    return (b - a) * 0.5 * sum(float(w) * function(float(t)) for t, w in zip(x, weights))


def active_spans(spec):
    p, U = spec["degree"], spec["knots"]
    return [(U[i], U[i + 1]) for i in range(p, len(spec["controls"])) if U[i] < U[i + 1]]


def exact_volume(recipe, order):
    lo, hi = recipe["lower"], recipe["upper"]
    total = 0.0
    for a, b in active_spans(lo):
        def integrand(t):
            lower, d_lower = curve(lo, t, True)
            upper, d_upper = curve(hi, t, True)
            if abs(lower[0] - upper[0]) > 1e-10:
                raise ValueError("x-boundaries do not coincide")
            if abs(d_lower[0] - d_upper[0]) > 1e-9:
                raise ValueError("x-derivatives do not coincide")
            gap = upper[1] - lower[1]
            return gap * d_lower[0] * recipe["height_mm"]
        total += gauss_integral(integrand, a, b, order)
    return total


def physical_state(recipe, xyz):
    x, y, z = xyz
    if z < 0 or z > recipe["height_mm"]:
        return "OUT"
    lo, hi = recipe["lower"], recipe["upper"]
    spans = active_spans(lo)
    x0 = curve(lo, spans[0][0])[0]
    x1 = curve(lo, spans[-1][1])[0]
    if x < x0 - 1e-9 or x > x1 + 1e-9:
        return "OUT"
    a, b = spans[0][0], spans[-1][1]
    for _ in range(70):
        m = (a + b) / 2
        if curve(lo, m)[0] < x:
            a = m
        else:
            b = m
    t = (a + b) / 2
    yl, yu = curve(lo, t)[1], curve(hi, t)[1]
    if y < yl - 1e-8 or y > yu + 1e-8:
        return "OUT"
    if min(abs(y - yl), abs(y - yu), abs(x - x0), abs(x - x1), abs(z), abs(z - recipe["height_mm"])) < 1e-8:
        return "ON"
    return "IN"


def bernstein(n, t):
    return [math.comb(n, i) * t**i * (1 - t)**(n-i) for i in range(n + 1)]


def evaluate_block(block, uvw):
    degree = block["degrees"]
    net = np.asarray(block["controls"], dtype=float).reshape(*(d + 1 for d in degree), 4)
    bs = [bernstein(n, t) for n, t in zip(degree, uvw)]
    h = np.zeros(4)
    for i in range(degree[0] + 1):
        for j in range(degree[1] + 1):
            for k in range(degree[2] + 1):
                h += bs[0][i] * bs[1][j] * bs[2][k] * net[i, j, k]
    return h[:3] / h[3]


def model_state(model, xyz):
    """Invert the serialized Bernstein strip maps without using the builder."""
    x, y, z = xyz
    for block in model["blocks"]:
        p0 = evaluate_block(block, (0.0, 0.5, 0.5))
        p1 = evaluate_block(block, (1.0, 0.5, 0.5))
        if not min(p0[0], p1[0]) - 1e-9 <= x <= max(p0[0], p1[0]) + 1e-9:
            continue
        if not -1e-9 <= z <= max(evaluate_block(block, (0.5, 0.5, 0.0))[2],
                                 evaluate_block(block, (0.5, 0.5, 1.0))[2]) + 1e-9:
            continue
        lo, hi = 0.0, 1.0
        for _ in range(64):
            mid = 0.5 * (lo + hi)
            if evaluate_block(block, (mid, 0.5, 0.5))[0] < x:
                lo = mid
            else:
                hi = mid
        u = 0.5 * (lo + hi)
        lower = evaluate_block(block, (u, 0.0, 0.5))[1]
        upper = evaluate_block(block, (u, 1.0, 0.5))[1]
        if min(lower, upper) - 1e-8 <= y <= max(lower, upper) + 1e-8:
            if min(abs(y-lower), abs(y-upper), abs(z), abs(z-model["height_mm"])) < 1e-8:
                return "ON"
            return "IN"
    return "OUT"


def occ_reference_shape(recipe):
    """Build the input boundary directly as an OCC rational B-spline wire."""
    def edge(spec):
        unique = []
        mults = []
        for knot in spec["knots"]:
            if unique and unique[-1] == knot:
                mults[-1] += 1
            else:
                unique.append(knot); mults.append(1)
        poles = [App.Vector(x/w, y/w, 0) for x, y, w in spec["controls"]]
        weights = [w for _, _, w in spec["controls"]]
        curve_obj = Part.BSplineCurve()
        curve_obj.buildFromPolesMultsKnots(poles, mults, unique, False, spec["degree"], weights)
        return curve_obj.toShape()
    lower = edge(recipe["lower"])
    upper_source = edge(recipe["upper"])
    upper = Part.Edge(upper_source.Curve, upper_source.LastParameter, upper_source.FirstParameter)
    start_l, end_l = lower.Vertexes[0].Point, lower.Vertexes[-1].Point
    start_u, end_u = upper_source.Vertexes[-1].Point, upper_source.Vertexes[0].Point
    right = Part.makeLine(end_l, start_u)
    left = Part.makeLine(end_u, start_l)
    face = Part.Face(Part.Wire([lower, right, upper, left]))
    return face.extrude(App.Vector(0, 0, recipe["height_mm"]))


def audit_case(cli, case, model_path, recipe_path, rng):
    recipe = case["recipe"]
    recipe_path.write_text(json.dumps(recipe, separators=(",", ":")))
    previous = model_path.read_bytes() if model_path.exists() else None
    start = time.monotonic()
    build = subprocess.run([str(cli), "build", str(recipe_path), str(model_path)], capture_output=True, text=True)
    elapsed = time.monotonic() - start
    try:
        diag = json.loads(build.stdout)
    except json.JSONDecodeError:
        diag = {"stdout": build.stdout, "stderr": build.stderr}
    row = {"id": case["id"], "expected_accepted": case["accepted"], "build_audit": diag,
           "build_s": elapsed, "qualified": (build.returncode == 0) == case["accepted"]}
    if build.returncode:
        row["previous_model_preserved"] = (model_path.read_bytes() if model_path.exists() else None) == previous
        row["qualified"] = not case["accepted"] and row["previous_model_preserved"]
        return row
    if not case["accepted"]:
        row["qualified"] = False
        return row
    model = json.loads(model_path.read_text())
    volume16, volume32, volume64 = (exact_volume(recipe, q) for q in (16, 32, 64))
    expected = case.get("expected_volume_mm3")
    volume_error = abs(volume64 - expected) if expected is not None else None
    direct_eval_error = 0.0
    trace_error = 0.0
    intervals = [tuple(i) for i in model.get("span_intervals", [])]
    blocks = model["blocks"]
    if len(intervals) != len(blocks):
        row["qualified"] = False
    for block, (a, b) in zip(blocks, intervals):
        for u in (0.137, 0.381, 0.719, 0.913):
            t = a + u * (b - a)
            for v in (0.0, 0.23, 0.67, 1.0):
                got = evaluate_block(block, (u, v, 0.37))
                low, high = curve(recipe["lower"], t), curve(recipe["upper"], t)
                expected_xy = low * (1-v) + high * v
                trace_error = max(trace_error, float(np.linalg.norm(got[:2] - expected_xy)))
                direct_eval_error = max(direct_eval_error, abs(got[2] - 0.37 * recipe["height_mm"]))
    points = []
    spans = active_spans(recipe["lower"])
    for a, b in spans:
        for q in (0.11, 0.37, 0.68, 0.91):
            t = a + q*(b-a)
            l, h = curve(recipe["lower"], t), curve(recipe["upper"], t)
            points.extend([(l[0], (l[1]+h[1])/2, recipe["height_mm"]/2),
                           (l[0], l[1]-0.01, recipe["height_mm"]/2),
                           (l[0], h[1]+0.01, recipe["height_mm"]/2)])
    states = {"IN": 0, "OUT": 0, "ON": 0}; mismatches = []
    for point in points:
        expected_state = physical_state(recipe, point)
        observed_state = model_state(model, point)
        states[expected_state] += 1
        if observed_state != expected_state:
            mismatches.append({"xyz": list(point), "expected": expected_state, "observed": observed_state})
    row.update(block_count=len(blocks), span_intervals=intervals,
               volume_quadrature_mm3={"q16": volume16, "q32": volume32, "q64": volume64},
               volume_quadrature_delta_q32_q64=abs(volume32-volume64),
               expected_volume_mm3=expected, absolute_volume_error_mm3=volume_error,
               max_control_map_error_mm=trace_error, max_extrusion_z_error_mm=direct_eval_error,
               classification_sample_count=len(points), classification_states=states,
               classification_mismatches=mismatches)
    row["qualified"] = bool(row["qualified"] and len(blocks) == len(spans) and len(blocks) == 2
                             and trace_error < 1e-9 and direct_eval_error < 1e-10
                             and abs(volume64-volume32) < 1e-9
                             and (volume_error is None or volume_error < 1e-9)
                             and diag.get("method") == "exact_monotone_nurbs_strip_extrusion")
    if Part is not None:
        step_path = model_path.with_name(case["id"] + ".step")
        exported = subprocess.run([str(cli), "export-step", str(model_path), str(step_path)],
                                  capture_output=True, text=True)
        try:
            export_diag = json.loads(exported.stdout)
        except json.JSONDecodeError:
            export_diag = {"stdout": exported.stdout, "stderr": exported.stderr}
        shape = Part.read(str(step_path)) if exported.returncode == 0 else None
        step_volume = float(shape.Volume) if shape else None
        adaptive_measure = None
        if shape is not None:
            brep = step_path.with_suffix(".brep")
            shape.exportBrep(str(brep))
            helper = model_path.parent / "occ-volume"
            if not helper.exists():
                prefix = Path(sys.prefix)
                subprocess.run(["g++", "-std=c++17", "-O2", "-I" + str(prefix / "include/opencascade"),
                                str(Path(__file__).with_name("occ_volume.cpp")), "-L" + str(prefix / "lib"),
                                "-Wl,-rpath," + str(prefix / "lib"), "-lTKTopAlgo", "-lTKBRep", "-lTKG3d",
                                "-lTKMath", "-lTKernel", "-o", str(helper)], check=True)
            raw_measure = json.loads(subprocess.check_output([str(helper), str(brep), "--centroid"], text=True))
            adaptive_measure = raw_measure[-1] if isinstance(raw_measure, list) else raw_measure
        reference_measure = None
        reference_shape_valid = None
        if shape is not None:
            reference = occ_reference_shape(recipe)
            reference_shape_valid = bool(reference.isValid() and len(reference.Solids) == 1
                                         and all(s.isClosed() for s in reference.Shells))
            reference_brep = step_path.with_name(step_path.stem + "-reference.brep")
            reference.exportBrep(str(reference_brep))
            raw_reference = json.loads(subprocess.check_output([str(helper), str(reference_brep), "--centroid"], text=True))
            reference_measure = raw_reference[-1] if isinstance(raw_reference, list) else raw_reference
        step_ok = bool(shape is not None and shape.isValid() and len(shape.Solids) == 1
                       and all(s.isClosed() for s in shape.Shells))
        row["step_occ"] = {"export": export_diag, "valid": bool(shape and shape.isValid()),
                            "solid_count": len(shape.Solids) if shape else 0,
                            "shell_count": len(shape.Shells) if shape else 0,
                            "all_shells_closed": bool(shape and all(s.isClosed() for s in shape.Shells)),
                            "volume_mm3": step_volume,
                            "adaptive_volume_mm3": adaptive_measure.get("volume_mm3") if adaptive_measure else None,
                            "adaptive_volume_error_vs_independent_mm3": abs(adaptive_measure["volume_mm3"]-volume64) if adaptive_measure else None,
                            "direct_occ_reference_volume_mm3": reference_measure.get("volume_mm3") if reference_measure else None,
                            "direct_occ_reference_error_vs_independent_mm3": abs(reference_measure["volume_mm3"]-volume64) if reference_measure else None,
                            "direct_occ_reference_valid_closed": reference_shape_valid,
                            "step_minus_direct_occ_reference_mm3": (adaptive_measure["volume_mm3"]-reference_measure["volume_mm3"]) if adaptive_measure and reference_measure else None,
                            "volume_error_vs_independent_mm3": abs(step_volume-volume64) if step_volume is not None else None,
                            "volume_comparison": "diagnostic_only_no_predeclared_tolerance",
                            "qualified": step_ok}
        row["qualified"] &= step_ok
    else:
        row["step_occ"] = {"status": "not_run", "reason": "FreeCAD/OCC Python bindings are unavailable in this interpreter"}
    return row


def run(cli, work):
    started = time.monotonic(); sha = hashlib.sha256(cli.read_bytes()).hexdigest()
    work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    model = work / "model.json"; recipe_path = work / "recipe.json"
    rows = [audit_case(cli, case, model, recipe_path, rng) for case in json.loads(CORPUS.read_text())]
    if hashlib.sha256(cli.read_bytes()).hexdigest() != sha:
        raise RuntimeError("CLI changed during oracle execution")
    return {"schema": "nurbs_sketch_strip_independent_oracle_v1",
            "provenance": {"cli": str(cli), "cli_sha256": sha,
                           "revision": subprocess.check_output(["rtk", "git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                           "python": platform.python_version(), "numpy": np.__version__, "seed": SEED,
                           "occ": getattr(Part, "OCC_VERSION", None),
                           "elapsed_s": time.monotonic()-started},
            "protocol": {"curve_oracle": "independent Cox-de Boor basis and rational quotient derivatives",
                         "volume_oracle": "Gauss-Legendre integration of height*(upper_y-lower_y)*lower_x' over original knot spans",
                         "geometry_oracle": "Bernstein evaluation of each extruded block against direct original B-spline evaluation",
                         "coverage": "three positive degree-3 two-span cases and three malformed or nonconforming cases"},
            "case_count": len(rows), "cases": rows, "qualified": all(row["qualified"] for row in rows)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.cli.resolve(), args.work)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    raise SystemExit(int(not report["qualified"]))
