#!/usr/bin/env python3
"""Independent one-sided geometric jet checks for rational B-spline curves.

The input curves are read from their source sketch representation. Derivatives
are obtained from derivative control polygons and the rational quotient rule;
no constructor-generated certificates are consumed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return f"<work>/{path.name}"


def _span(degree: int, knots: list, controls: list, t, side: str) -> int:
    n = len(controls) - 1
    if t == knots[n + 1] or side == "left" and any(t == u for u in knots):
        candidates = [i for i in range(degree, n + 1) if knots[i] < knots[i + 1] and knots[i + 1] == t]
        if candidates:
            return candidates[-1]
    candidates = [i for i in range(degree, n + 1) if knots[i] <= t < knots[i + 1]]
    if not candidates and t == knots[n + 1]:
        return n
    if not candidates:
        raise ValueError(f"parameter {t} lies outside B-spline domain")
    return candidates[0] if side == "right" else candidates[-1]


def de_boor(degree: int, knots: list, controls: list, t, side: str):
    k = _span(degree, knots, controls, t, side)
    work = [list(controls[j]) for j in range(k - degree, k + 1)]
    for r in range(1, degree + 1):
        for j in range(degree, r - 1, -1):
            i = k - degree + j
            den = knots[i + degree - r + 1] - knots[i]
            a = (t - knots[i]) / den if den else 0
            work[j] = [(1 - a) * work[j - 1][c] + a * work[j][c]
                       for c in range(len(work[j]))]
    return work[degree]


def derivative_spline(degree: int, knots: list, controls: list):
    if degree <= 0:
        return 0, [knots[0], knots[-1]], [[0 for _ in controls[0]]]
    out = []
    for i in range(len(controls) - 1):
        den = knots[i + degree + 1] - knots[i + 1]
        out.append([degree * (controls[i + 1][c] - controls[i][c]) / den if den else 0
                    for c in range(len(controls[i]))])
    return degree - 1, knots[1:-1], out


def homogeneous_jets(curve: dict, t, side: str):
    degree = int(curve["degree"])
    knots = list(curve["knots"])
    controls = [list(map(float, p)) for p in curve["controls"]]
    h0 = de_boor(degree, knots, controls, float(t), side)
    p1, k1, c1 = derivative_spline(degree, knots, controls)
    h1 = de_boor(p1, k1, c1, float(t), side)
    p2, k2, c2 = derivative_spline(p1, k1, c1)
    h2 = de_boor(p2, k2, c2, float(t), side)
    return h0, h1, h2


def rational_jets(curve: dict, t, side: str):
    h, dh, ddh = homogeneous_jets(curve, t, side)
    weight, dweight, ddweight = h[-1], dh[-1], ddh[-1]
    if weight <= 0 or not math.isfinite(weight):
        raise ArithmeticError("nonpositive homogeneous weight at jet sample")
    point = [h[c] / weight for c in range(2)]
    first = [(dh[c] - point[c] * dweight) / weight for c in range(2)]
    second = [(ddh[c] - 2 * first[c] * dweight - point[c] * ddweight) / weight
              for c in range(2)]
    speed = math.hypot(*first)
    if speed <= 1e-14 or not math.isfinite(speed):
        raise ArithmeticError("zero or nonfinite tangent at knot jet")
    tangent = [v / speed for v in first]
    tangential_accel = sum(tangent[c] * second[c] for c in range(2))
    curvature = [(second[c] - tangent[c] * tangential_accel) / (speed * speed)
                 for c in range(2)]
    return {"position_mm": point, "first_derivative": first,
            "unit_tangent_oriented": tangent, "second_derivative": second,
            "curvature_vector_per_mm": curvature, "speed_mm_per_parameter": speed}


def unique_internal_knots(curve: dict) -> list:
    p, knots, controls = int(curve["degree"]), curve["knots"], curve["controls"]
    lo, hi = knots[p], knots[len(controls)]
    return sorted({u for u in knots if lo < u < hi})


def compare_jets(curve: dict, parameter, left_side: str, right_side: str,
                 position_tol: float, tangent_tol: float, curvature_tol: float) -> dict:
    left = rational_jets(curve, parameter, left_side)
    right = rational_jets(curve, parameter, right_side)
    position_jump = math.dist(left["position_mm"], right["position_mm"])
    tangent_jump = math.dist(left["unit_tangent_oriented"], right["unit_tangent_oriented"])
    curvature_jump = math.dist(left["curvature_vector_per_mm"], right["curvature_vector_per_mm"])
    c0 = position_jump <= position_tol
    g1 = c0 and tangent_jump <= tangent_tol
    g2 = g1 and curvature_jump <= curvature_tol
    left_normal = [left["unit_tangent_oriented"][1], -left["unit_tangent_oriented"][0], 0.0]
    right_normal = [right["unit_tangent_oriented"][1], -right["unit_tangent_oriented"][0], 0.0]
    return {"parameter": parameter, "left": left, "right": right,
            "position_jump_mm": position_jump, "oriented_unit_tangent_jump": tangent_jump,
            "curvature_vector_jump_per_mm": curvature_jump,
            "extruded_surface": {
                "oriented_unit_normal_jump": math.dist(left_normal, right_normal),
                "section_curvature_jump_per_mm": curvature_jump,
                "generator_curvature_left_per_mm": 0.0,
                "generator_curvature_right_per_mm": 0.0,
                "generator_curvature_jump_per_mm": 0.0,
            },
            "classification": "G2" if g2 else "G1_only" if g1 else "C0_only" if c0 else "discontinuous",
            "c0": c0, "g1": g1, "g2": g2}


def _fraction_deboor(curve: dict, t: Fraction) -> tuple[Fraction, Fraction]:
    degree = int(curve["degree"])
    knots = [Fraction(k) for k in curve["knots"]]
    controls = [[Fraction(v) for v in p] for p in curve["controls"]]
    h = de_boor(degree, knots, controls, t, "right")
    return h[0] / h[-1], h[1] / h[-1]


def fraction_rank(rows: list[list[Fraction]]) -> int:
    a = [row[:] for row in rows]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((r for r in range(rank, len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        q = a[rank][col]
        a[rank] = [v / q for v in a[rank]]
        for r in range(len(a)):
            if r != rank and a[r][col]:
                q = a[r][col]
                a[r] = [x - q * y for x, y in zip(a[r], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def conic_design_rank(curve: dict) -> dict:
    degree, knots, controls = int(curve["degree"]), curve["knots"], curve["controls"]
    lo, hi = Fraction(knots[degree]), Fraction(knots[len(controls)])
    parameters = [lo + (hi - lo) * Fraction(i, 8) for i in range(1, 8)]
    points = [_fraction_deboor(curve, t) for t in parameters]
    matrix = [[x*x, x*y, y*y, x, y, Fraction(1)] for x, y in points]
    rank = fraction_rank(matrix)
    return {"exact_sample_count": len(points), "conic_design_matrix_rank": rank,
            "certifies_not_one_conic": rank == 6}


def analyze_curve(item: dict, tolerances: dict) -> dict:
    curve = item.get("curve", item)
    knots, degree, controls = curve["knots"], int(curve["degree"]), curve["controls"]
    domain = [knots[degree], knots[len(controls)]]
    knots_report = [compare_jets(curve, t, "left", "right", **tolerances)
                    for t in unique_internal_knots(curve)]
    closure = None
    if item.get("closed", curve.get("closed", False)):
        closure = compare_jets(curve, domain[0], "left", "right", **tolerances)
        # Explicitly evaluate the left limit at the final parameter against the
        # right limit at the initial parameter, as periodicity requires.
        end = rational_jets(curve, domain[1], "left")
        start = rational_jets(curve, domain[0], "right")
        end_normal = [end["unit_tangent_oriented"][1], -end["unit_tangent_oriented"][0], 0.0]
        start_normal = [start["unit_tangent_oriented"][1], -start["unit_tangent_oriented"][0], 0.0]
        closure["periodic_end_to_start"] = {
            "position_jump_mm": math.dist(end["position_mm"], start["position_mm"]),
            "oriented_unit_tangent_jump": math.dist(end["unit_tangent_oriented"], start["unit_tangent_oriented"]),
            "curvature_vector_jump_per_mm": math.dist(end["curvature_vector_per_mm"], start["curvature_vector_per_mm"]),
            "extruded_surface_normal_jump": math.dist(end_normal, start_normal),
            "generator_curvature_jump_per_mm": 0.0,
        }
        p = closure["periodic_end_to_start"]
        closure["periodic_end_to_start"]["classification"] = (
            "G2" if p["position_jump_mm"] <= tolerances["position_tol"]
            and p["oriented_unit_tangent_jump"] <= tolerances["tangent_tol"]
            and p["curvature_vector_jump_per_mm"] <= tolerances["curvature_tol"]
            else "G1_only" if p["position_jump_mm"] <= tolerances["position_tol"]
            and p["oriented_unit_tangent_jump"] <= tolerances["tangent_tol"]
            else "C0_only" if p["position_jump_mm"] <= tolerances["position_tol"]
            else "discontinuous")
    all_junctions = knots_report + ([closure] if closure else [])
    expected = item.get("expect", curve.get("expect"))
    rank = conic_design_rank(curve) if item.get("test_nonconic", curve.get("test_nonconic", False)) else None
    return {"id": item.get("id", "curve"), "domain": domain,
            "internal_junctions": knots_report, "periodic_closure": closure,
            "all_internal_and_closure_g2": bool(all_junctions) and all(j["g2"] for j in knots_report)
                and (closure is None or closure.get("periodic_end_to_start", {}).get("classification") == "G2"),
            "expected_classification": expected,
            "expectation_met": None if expected is None else (
                all(j["g2"] for j in knots_report) and (closure is None or closure["periodic_end_to_start"]["classification"] == "G2")
                if expected == "G2" else (
                    any(j["g1"] and not j["g2"] for j in knots_report)
                    or closure is not None and closure["periodic_end_to_start"]["classification"] == "G1_only")
                if expected == "G1_only" else True),
            "nonconic_diagnostic": rank}


def algorithm_self_check(tolerances: dict) -> dict:
    parabola = {"degree": 2, "knots": [0, 0, 0, 1, 1, 1],
                "controls": [[0, 0, 1], [0.5, 0, 1], [1, 1, 1]]}
    jet = rational_jets(parabola, 0.5, "right")
    position_error = math.dist(jet["position_mm"], [0.5, 0.25])
    curvature_error = math.dist(jet["curvature_vector_per_mm"], [-0.5, 0.5])
    if position_error > 1e-14 or curvature_error > 1e-14:
        raise AssertionError("analytic parabola jet self-check failed")
    g1_witness = {"degree": 2, "knots": [0, 0, 0, 1, 2, 2, 2],
                  "controls": [[0, 0, 1], [0.5, 0, 1], [1, 1, 1], [2, 5, 1]]}
    witness = compare_jets(g1_witness, 1, "left", "right", **tolerances)
    if witness["classification"] != "G1_only":
        raise AssertionError("G1 non-G2 witness was not distinguished")
    return {"analytic_parabola_position_error_mm": position_error,
            "analytic_parabola_curvature_error_per_mm": curvature_error,
            "g1_non_g2_witness": witness}


def rectangle_corner_jets(rectangle: list) -> list[dict]:
    (xmin, ymin), (xmax, ymax) = rectangle
    points = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
    result = []
    for i, point in enumerate(points):
        previous, following = points[i - 1], points[(i + 1) % len(points)]
        incoming = [point[k] - previous[k] for k in range(2)]
        outgoing = [following[k] - point[k] for k in range(2)]
        incoming = [x / math.hypot(*incoming) for x in incoming]
        outgoing = [x / math.hypot(*outgoing) for x in outgoing]
        result.append({"position_mm": list(point), "incoming_unit_tangent": incoming,
                       "outgoing_unit_tangent": outgoing,
                       "oriented_tangent_jump": math.dist(incoming, outgoing),
                       "classification": "G0_sharp"})
    return result


def periodic_area(curve: dict, order: int = 32) -> float:
    import numpy as np
    points, weights = np.polynomial.legendre.leggauss(order)
    spans = sorted({(curve["knots"][i], curve["knots"][i + 1])
                    for i in range(int(curve["degree"]), len(curve["controls"]))
                    if curve["knots"][i] < curve["knots"][i + 1]})
    area = 0.0
    for lo, hi in spans:
        for x, w in zip(points, weights):
            t = (lo + hi) * 0.5 + (hi - lo) * x * 0.5
            jet = rational_jets(curve, t, "right")
            px, py = jet["position_mm"]
            dx, dy = jet["first_derivative"]
            area += 0.25 * (hi - lo) * w * (px * dy - py * dx)
    return area


def _bezier(points: list[list[float]]) -> dict:
    return {"degree": 3, "knots": [0, 0, 0, 0, 1, 1, 1, 1],
            "controls": [[x, y, 1.0] for x, y in points]}


def _curve_join(left_curve: dict, right_curve: dict, tolerances: dict) -> dict:
    left = rational_jets(left_curve, left_curve["knots"][len(left_curve["controls"])], "left")
    right = rational_jets(right_curve, right_curve["knots"][right_curve["degree"]], "right")
    position = math.dist(left["position_mm"], right["position_mm"])
    tangent = math.dist(left["unit_tangent_oriented"], right["unit_tangent_oriented"])
    curvature = math.dist(left["curvature_vector_per_mm"], right["curvature_vector_per_mm"])
    return {"position_jump_mm": position, "oriented_unit_tangent_jump": tangent,
            "curvature_vector_jump_per_mm": curvature,
            "classification": "G2" if position <= tolerances["position_tol"]
                and tangent <= tolerances["tangent_tol"]
                and curvature <= tolerances["curvature_tol"] else "G1_only"
                if position <= tolerances["position_tol"] and tangent <= tolerances["tangent_tol"]
                else "C0_only" if position <= tolerances["position_tol"] else "discontinuous"}


def rounded_rectangle_jets(rectangle: list, radius: float, tolerances: dict) -> dict:
    (xmin, ymin), (xmax, ymax) = rectangle
    if not 0 < radius < min(xmax - xmin, ymax - ymin) / 2:
        raise ValueError("rounded-rectangle radius must fit inside the rectangle")
    bl_out, br_in = [xmin + radius, ymin], [xmax - radius, ymin]
    br_out, tr_in = [xmax, ymin + radius], [xmax, ymax - radius]
    tr_out, tl_in = [xmax - radius, ymax], [xmin + radius, ymax]
    tl_out, bl_in = [xmin, ymax - radius], [xmin, ymin + radius]
    corners = [[xmax, ymin], [xmax, ymax], [xmin, ymax], [xmin, ymin]]
    starts = [br_in, tr_in, tl_in, bl_in]
    ends = [br_out, tr_out, tl_out, bl_out]
    # Order is bottom edge, BR arc, right edge, TR arc, top edge, TL arc,
    # left edge, BL arc. The [A,C,C,B] arcs have zero normal curvature at
    # each endpoint, matching their straight neighbours through G2.
    segments = [
        _bezier([bl_out, [(2 * bl_out[0] + br_in[0]) / 3, ymin],
                 [(bl_out[0] + 2 * br_in[0]) / 3, ymin], br_in]),
        _bezier([br_in, corners[0], corners[0], br_out]),
        _bezier([br_out, [xmax, (2 * br_out[1] + tr_in[1]) / 3],
                 [xmax, (br_out[1] + 2 * tr_in[1]) / 3], tr_in]),
        _bezier([tr_in, corners[1], corners[1], tr_out]),
        _bezier([tr_out, [(2 * tr_out[0] + tl_in[0]) / 3, ymax],
                 [(tr_out[0] + 2 * tl_in[0]) / 3, ymax], tl_in]),
        _bezier([tl_in, corners[2], corners[2], tl_out]),
        _bezier([tl_out, [xmin, (2 * tl_out[1] + bl_in[1]) / 3],
                 [xmin, (tl_out[1] + 2 * bl_in[1]) / 3], bl_in]),
        _bezier([bl_in, corners[3], corners[3], bl_out]),
    ]
    joins = [_curve_join(segments[i], segments[(i + 1) % len(segments)], tolerances)
             for i in range(len(segments))]
    return {"radius_mm": radius, "segments": len(segments), "joins_including_closure": joins,
            "all_joins_g2": all(j["classification"] == "G2" for j in joins),
            "construction": "four straight cubic Bezier segments and four cubic [A,C,C,B] corner segments"}


def run_cli(cli: Path, command: list[str]) -> dict:
    process = subprocess.run([str(cli), *command], cwd=ROOT, capture_output=True, text=True)
    if process.returncode:
        raise RuntimeError(f"CLI {' '.join(command)} failed: {process.stderr or process.stdout}")
    return json.loads(process.stdout)


def c0_iga_check(cli: Path, recipe_path: Path, work: Path) -> dict:
    import hashlib
    import importlib.util
    def load(name, filename):
        spec = importlib.util.spec_from_file_location(name, ROOT / "tests/references" / filename)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    space = load("g2_plate_space_oracle", "nurbs_space_oracle.py")
    patch = load("g2_plate_equilibrium_oracle", "nurbs_patch_equilibrium.py")
    model_path, package_path = work / "plate.model.json", work / "plate.space.json"
    build = run_cli(cli, ["build", str(recipe_path), str(model_path)])
    model = json.loads(model_path.read_text())
    if model.get("source_recipe") != json.loads(recipe_path.read_text()):
        raise AssertionError("built plate does not preserve the exact source recipe")
    export = run_cli(cli, ["export-nurbs-space", str(model_path), str(package_path)])
    package = json.loads(package_path.read_text())
    trace = space.package_checks(model, package)
    orders = [patch.integrate(package, q) for q in (4, 8, 16)]
    if orders[-1]["relative_residual_linf"] > 1e-8:
        raise AssertionError(f"q16 C0 constant-stress patch residual too high: {orders[-1]}")
    return {"status": "verified", "cli": display_path(cli),
            "cli_sha256": hashlib.sha256(cli.read_bytes()).hexdigest(),
            "build_audit": build, "space_export": export,
            "exact_c0_trace_and_affine_checks": trace,
            "constant_stress_equilibrium_orders": orders,
            "q16_relative_linf_threshold": 1e-8,
            "source_recipe_preserved": True,
            "artifacts": {"model": display_path(model_path),
                          "q3_package": display_path(package_path)}}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("recipe", type=Path)
    ap.add_argument("--cli", type=Path,
                    help="optional production CLI for C0 Q3 trace, affine field, and patch checks")
    ap.add_argument("--work", type=Path, default=ROOT / "build/g2-plate-jets")
    ap.add_argument("--output", type=Path, default=ROOT / "docs/results/g2-plate-jets.json")
    ap.add_argument("--position-tol", type=float, default=1e-10)
    ap.add_argument("--tangent-tol", type=float, default=1e-10)
    ap.add_argument("--curvature-tol", type=float, default=1e-10)
    args = ap.parse_args()
    recipe = json.loads(args.recipe.read_text())
    tolerances = {"position_tol": args.position_tol, "tangent_tol": args.tangent_tol,
                  "curvature_tol": args.curvature_tol}
    self_check = algorithm_self_check(tolerances)
    raw_curves = recipe.get("curves")
    if raw_curves is None and recipe.get("profile"):
        raw_curves = [{"id": "periodic-free-profile", "curve": recipe["profile"],
                       "closed": True, "expect": "G2", "test_nonconic": True}]
        if recipe.get("g1_witness"):
            raw_curves.append({"id": "g1-not-g2-witness", "curve": recipe["g1_witness"],
                               "expect": "G1_only"})
    if raw_curves is None and recipe.get("edges"):
        raw_curves = [{"id": edge["id"], "curve": edge["curve"]} for edge in recipe["edges"]]
    if raw_curves is None:
        raise SystemExit("recipe must provide curves:[{id,curve,closed,expect}]")
    rows = [analyze_curve(item, tolerances) for item in raw_curves]
    expected_rows = [row for row in rows if row["expected_classification"]]
    rounded = rounded_rectangle_jets(recipe["rectangle_mm"],
                                     float(recipe.get("rounded_witness_radius_mm", 5.0)), tolerances)
    passed = all(row["expectation_met"] for row in expected_rows) and rounded["all_joins_g2"]
    profile_area = periodic_area(recipe["profile"]) if "profile" in recipe else None
    volume_checks = None
    if profile_area is not None and "bounds_xyz_mm" in recipe and "rectangle_mm" in recipe:
        (xmin, xmax), (ymin, ymax), (zmin, zmax) = recipe["bounds_xyz_mm"]
        (rx0, ry0), (rx1, ry1) = recipe["rectangle_mm"]
        stock = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        rectangle_hole = (rx1 - rx0) * (ry1 - ry0)
        volume_checks = {"stock_volume_mm3": stock,
                         "rectangle_hole_area_mm2": rectangle_hole,
                         "periodic_profile_area_mm2": profile_area,
                         "independent_remaining_volume_mm3": stock - (rectangle_hole + profile_area) * (zmax - zmin),
                         "method": "exact rectangular stock/rectangle volume plus independent periodic B-spline Green integral"}
    recipe_path = args.recipe.resolve()
    report = {"schema": "g2_curve_jets_report_v1", "status": "verified" if passed else "failed",
              "method": "independent homogeneous B-spline derivative control polygons and rational quotient-rule jets; oriented unit tangent and geometric curvature vector",
              "tolerances": tolerances,
              "source_recipe": display_path(recipe_path),
              "source_recipe_sha256": hashlib.sha256(recipe_path.read_bytes()).hexdigest(),
              "source_recipe_renamed_from": "tests/references/g2_plate.recipe.json",
              "curves": rows,
              "algorithm_self_check": self_check,
              "rectangle_sharp_junctions": rectangle_corner_jets(recipe["rectangle_mm"])
                  if "rectangle_mm" in recipe else None,
              "rounded_rectangle_g2_witness": rounded if "rectangle_mm" in recipe else None,
              "independent_periodic_profile_area_mm2": profile_area,
              "independent_volume_check": volume_checks,
              "limits": ["These one-sided differential checks qualify only the supplied source curves and junctions.",
                         "A finite exact conic design rank of six proves that the sampled curve is not one conic; it does not establish topology or injectivity."]}
    if args.cli:
        try:
            report["c0_iga"] = c0_iga_check(args.cli.resolve(), args.recipe.resolve(), args.work.resolve())
        except Exception as error:
            report["c0_iga"] = {"status": "failed", "error": str(error)}
            report["status"] = "failed"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "curves": len(rows), "output": str(args.output)}))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
