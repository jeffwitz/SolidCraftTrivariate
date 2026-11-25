#!/usr/bin/env python3
"""Independent de Boor oracle for variable-profile intersecting X/Z bores.

The oracle reads only the design recipe and the CLI-published block model. It
evaluates the four weighted cubic profiles itself and integrates section area
with composite Simpson quadrature. No production geometry routines are used.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
from pathlib import Path
from fractions import Fraction

HERE = Path(__file__).resolve().parent
RECIPE = HERE / "nurbs_cross_bores.recipe.json"


def de_boor(curve: dict, t: float) -> tuple[float, float]:
    p, knots = curve["degree"], curve["knots"]
    controls = curve["controls"]
    n = len(controls) - 1
    if not knots[p] <= t <= knots[n + 1]:
        raise ValueError("parameter outside NURBS domain")
    span = n if t == knots[n + 1] else next(
        i for i in range(p, n + 1) if knots[i] <= t < knots[i + 1]
    )
    d = [list(map(float, controls[j])) for j in range(span - p, span + 1)]
    for r in range(1, p + 1):
        for j in range(p, r - 1, -1):
            i = span - p + j
            den = knots[i + p - r + 1] - knots[i]
            a = 0.0 if den == 0 else (t - knots[i]) / den
            d[j] = [(1 - a) * d[j - 1][k] + a * d[j][k] for k in range(3)]
    xw, yw, w = d[p]
    return xw / w, yw / w


def at_y(curve: dict, y: float) -> float:
    lo, hi = curve["knots"][curve["degree"]], curve["knots"][-curve["degree"] - 1]
    for _ in range(60):
        mid = (lo + hi) / 2
        _, ym = de_boor(curve, mid)
        if ym < y:
            lo = mid
        else:
            hi = mid
    return de_boor(curve, (lo + hi) / 2)[0]


def derivative_xy(curve: dict, y: float, h: float = 1.0e-4) -> float:
    """Independent centered d(x)/d(y), using the inverse parameter map."""
    return (at_y(curve, y + h) - at_y(curve, y - h)) / (2 * h)


def exact_rank(rows: list[list[Fraction]]) -> int:
    """Small exact Gauss-Jordan rank, used for the conic regression witness."""
    a = [row[:] for row in rows]
    rank, col = 0, 0
    while rank < len(a) and col < len(a[0]):
        pivot = next((r for r in range(rank, len(a)) if a[r][col]), None)
        if pivot is None:
            col += 1
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [v / scale for v in a[rank]]
        for r in range(len(a)):
            if r != rank and a[r][col]:
                factor = a[r][col]
                a[r] = [u - factor * v for u, v in zip(a[r], a[rank])]
        rank += 1
        col += 1
    return rank


def conic_witness(curve: dict) -> dict:
    """Show exact rank six for [x²,xy,y²,x,y,1] on the first knot span."""
    from bisect import bisect_right
    p, knots, controls = curve["degree"], curve["knots"], curve["controls"]
    lo, hi = knots[p], next(k for k in knots[p + 1:] if k > knots[p])
    # Seven rational parameters in one span overdetermine the six-dimensional
    # conic polynomial space. Homogeneous curve values stay exact.
    samples = [Fraction(lo) + Fraction(hi - lo) * Fraction(i, 6) for i in range(7)]
    rows = []
    for t in samples:
        tt = Fraction(t)
        span = len(controls) - 1 if tt == knots[-p - 1] else bisect_right(knots, tt) - 1
        d = [[Fraction(v) for v in controls[j]] for j in range(span - p, span + 1)]
        for r in range(1, p + 1):
            for j in range(p, r - 1, -1):
                i = span - p + j
                den = knots[i + p - r + 1] - knots[i]
                alpha = Fraction(0) if den == 0 else (tt - knots[i]) / den
                d[j] = [(1-alpha)*d[j-1][k] + alpha*d[j][k] for k in range(3)]
        xw, yw, w = d[p]
        x, y = xw / w, yw / w
        rows.append([x*x, x*y, y*y, x, y, Fraction(1)])
    return {"span": [lo, hi], "rational_samples": len(samples),
            "conic_design_matrix_exact_rank": exact_rank(rows),
            "columns": ["x^2", "x*y", "y^2", "x", "y", "1"]}


def bernstein(n: int, t: float) -> list[float]:
    return [math.comb(n, i) * t**i * (1 - t) ** (n - i) for i in range(n + 1)]


def chart_value(block: dict, uvw: tuple[float, float, float]) -> tuple[list[float], list[list[float]]]:
    degrees, controls = block["degrees"], block["controls"]
    basis = [bernstein(d, t) for d, t in zip(degrees, uvw)]
    deriv = []
    for axis, (degree, t) in enumerate(zip(degrees, uvw)):
        low = bernstein(degree - 1, t)
        deriv.append([degree * ((low[i - 1] if i else 0.0) - (low[i] if i < degree else 0.0))
                      for i in range(degree + 1)])
    h = [0.0] * 4
    dh = [[0.0] * 4 for _ in range(3)]
    for i, bi in enumerate(basis[0]):
        for j, bj in enumerate(basis[1]):
            for k, bk in enumerate(basis[2]):
                p = controls[(i * (degrees[1] + 1) + j) * (degrees[2] + 1) + k]
                b = bi * bj * bk
                for c in range(4):
                    h[c] += b * p[c]
                    dh[0][c] += deriv[0][i] * bj * bk * p[c]
                    dh[1][c] += bi * deriv[1][j] * bk * p[c]
                    dh[2][c] += bi * bj * deriv[2][k] * p[c]
    xyz = [h[c] / h[3] for c in range(3)]
    jac = [[(dh[a][c] * h[3] - h[c] * dh[a][3]) / h[3] ** 2 for a in range(3)] for c in range(3)]
    return xyz, jac


def gauss_legendre(n: int) -> tuple[list[float], list[float]]:
    nodes, weights = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(20):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            pn = x if n == 1 else p1
            pnm1 = 1.0 if n == 1 else p0
            dp = n * (x * pn - pnm1) / (x * x - 1)
            step = pn / dp
            x -= step
            if abs(step) < 2e-16:
                break
        nodes.append(x)
        weights.append(2 / ((1 - x * x) * dp * dp))
    pairs = sorted(zip(nodes, weights))
    return [x for x, _ in pairs], [w for _, w in pairs]


def det3(j: list[list[float]]) -> float:
    return (j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
            - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
            + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]))


def block_volume(blocks: list[dict], order: int) -> float:
    qx, qw = gauss_legendre(order)
    points = [(0.5 * (x + 1), 0.5 * w) for x, w in zip(qx, qw)]
    total = 0.0
    for block in blocks:
        for u, wu in points:
            for v, wv in points:
                for w, ww in points:
                    _, jac = chart_value(block, (u, v, w))
                    det = det3(jac)
                    if det <= 0:
                        raise ValueError(f"nonpositive chart determinant in {block.get('id')}: {det}")
                    total += wu * wv * ww * det
    return total


def independent_classification(recipe: dict, model: dict, subdivisions: int = 3) -> dict:
    profiles, bounds = recipe["profiles"], recipe["bounds_mm"]
    checked, disagreements = 0, []
    for block in model.get("blocks", []):
        for i in range(1, subdivisions + 1):
            for j in range(1, subdivisions + 1):
                for k in range(1, subdivisions + 1):
                    p, _ = chart_value(block, (i / (subdivisions + 1), j / (subdivisions + 1), k / (subdivisions + 1)))
                    x, y, z = p
                    removed = False
                    if 8 < y < 30:
                        removed = (at_y(profiles["z_min"], y) < z < at_y(profiles["z_max"], y)
                                   or at_y(profiles["x_min"], y) < x < at_y(profiles["x_max"], y))
                    checked += 1
                    if removed:
                        disagreements.append({"block": block.get("id"), "point_mm": p})
                        if len(disagreements) == 5:
                            break
                if len(disagreements) == 5:
                    break
            if len(disagreements) == 5:
                break
    return {"strict_interior_points_checked": checked,
            "source_removed_region_disagreements": disagreements}


def model_profiles(model: dict, bounds: list[list[int]]) -> dict[str, list[tuple[float, float]]]:
    """Recover four section walls from actual homogeneous Rust block nets."""
    rows: dict[float, tuple[set[float], set[float]]] = {}
    for block in model.get("blocks", []):
        du, dv, dw = block["degrees"]
        controls = block["controls"]
        for j in range(dv + 1):
            xs, zs = [], []
            for i in (0, du):
                for k in (0, dw):
                    p = controls[(i * (dv + 1) + j) * (dw + 1) + k]
                    xs.append(p[0] / p[3])
                    zs.append(p[2] / p[3])
                    y = p[1] / p[3]
            key = round(y, 12)
            gx, gz = rows.setdefault(key, (set(), set()))
            gx.update(xs)
            gz.update(zs)
    result = {name: [] for name in ("x_min", "x_max", "z_min", "z_max")}
    for y, (xs, zs) in sorted(rows.items()):
        ix = sorted(x for x in xs if bounds[0][0] + 1e-9 < x < bounds[0][1] - 1e-9)
        iz = sorted(z for z in zs if bounds[2][0] + 1e-9 < z < bounds[2][1] - 1e-9)
        if 8 - 1e-9 <= y <= 30 + 1e-9:
            if len(ix) == 2 and len(iz) == 2:
                for name, value in zip(result, (*ix, *iz)):
                    result[name].append((y, value))
    return result


def polyline_value(points: list[tuple[float, float]], y: float) -> float:
    for (y0, x0), (y1, x1) in zip(points, points[1:]):
        if y0 - 1e-11 <= y <= y1 + 1e-11:
            t = min(1.0, max(0.0, (y - y0) / (y1 - y0)))
            return x0 + t * (x1 - x0)
    return points[0][1] if y < points[0][0] else points[-1][1]


def profile_deviation(recipe: dict, model: dict) -> dict:
    actual = model_profiles(model, recipe["bounds_mm"])
    output = {}
    for name, curve in recipe["profiles"].items():
        points = actual[name]
        max_error = 0.0
        worst_y = None
        if len(points) < 2:
            output[name] = {"model_profile_rows": len(points), "max_lateral_error_mm": None}
            continue
        for (y0, _), (y1, _) in zip(points, points[1:]):
            for i in range(33):
                y = y0 + (y1 - y0) * i / 32
                if 8 <= y <= 30:
                    error = abs(at_y(curve, y) - polyline_value(points, y))
                    if error > max_error:
                        max_error, worst_y = error, y
        output[name] = {"model_profile_rows": len(points), "max_lateral_error_mm": max_error,
                        "worst_y_mm": worst_y}
    return output


def volume_from_model_profiles(model: dict, recipe: dict, panels: int = 4096) -> float | None:
    profiles = model_profiles(model, recipe["bounds_mm"])
    if any(len(p) < 2 for p in profiles.values()):
        return None
    dy = 22.0 / panels
    values = []
    for i in range(panels + 1):
        y = 8 + i * dy
        xw = polyline_value(profiles["x_max"], y) - polyline_value(profiles["x_min"], y)
        zw = polyline_value(profiles["z_max"], y) - polyline_value(profiles["z_min"], y)
        values.append((60 - xw) * (30 - zw))
    active_material = dy / 3 * (values[0] + values[-1] + 4 * sum(values[1:-1:2]) + 2 * sum(values[2:-1:2]))
    return 60 * (40 - 22) * 30 + active_material


def area_removed(recipe: dict, y: float) -> float:
    profiles = recipe["profiles"]
    x_width = at_y(profiles["x_max"], y) - at_y(profiles["x_min"], y)
    z_width = at_y(profiles["z_max"], y) - at_y(profiles["z_min"], y)
    # X-axis bore removes a z strip across the stock; Z-axis bore removes
    # an x strip. Inclusion-exclusion handles their true overlap.
    return 60 * z_width + 30 * x_width - x_width * z_width


def expected_volume(recipe: dict, panels: int = 8192) -> float:
    # Profiles are active only between their cap planes y=8 and y=30.
    # Outside that interval the stock is uncut.
    dy = 22.0 / panels
    values = [area_removed(recipe, 8 + i * dy) for i in range(panels + 1)]
    removed = dy / 3 * (values[0] + values[-1] + 4 * sum(values[1:-1:2]) + 2 * sum(values[2:-1:2]))
    return 60 * 40 * 30 - removed


def run(cli: Path, work: Path) -> dict:
    recipe = json.loads(RECIPE.read_text())
    work.mkdir(parents=True, exist_ok=True)
    model_path = work / "model.json"
    built = subprocess.run([str(cli), "build", str(RECIPE), str(model_path)], text=True, capture_output=True)
    if built.returncode:
        raise RuntimeError(f"CLI build failed ({built.returncode}): {built.stderr or built.stdout}")
    report = json.loads(built.stdout)
    model = json.loads(model_path.read_text())
    expected = expected_volume(recipe)
    volume_q12 = block_volume(model.get("blocks", []), 12)
    volume_q16 = block_volume(model.get("blocks", []), 16)
    volume_error = volume_q16 - expected

    # End-point, interior, and derivative observations from the source curve.
    profile_observations = {}
    for name, curve in recipe["profiles"].items():
        profile_observations[name] = {
            "y8": at_y(curve, 8), "y18": at_y(curve, 18), "y30": at_y(curve, 30),
            "dx_dy_at18": derivative_xy(curve, 18),
            "strictly_increasing_y": all(de_boor(curve, curve["knots"][3] + (curve["knots"][-4] - curve["knots"][3]) * i / 100)[1]
                                          < de_boor(curve, curve["knots"][3] + (curve["knots"][-4] - curve["knots"][3]) * (i + 1) / 100)[1]
                                          for i in range(100)),
            "nonconic_first_span": conic_witness(curve),
        }

    # Rejection must leave an existing publication byte-for-byte intact.
    sentinel = work / "preserved.json"
    sentinel.write_text('{"sentinel":"keep"}\n')
    invalid = json.loads(json.dumps(recipe))
    invalid["profiles"]["x_min"]["controls"][2][1] = 36
    invalid_recipe = work / "desynchronized.json"
    invalid_recipe.write_text(json.dumps(invalid))
    refused = subprocess.run([str(cli), "build", str(invalid_recipe), str(sentinel)], text=True, capture_output=True)
    preserved_desynchronized = sentinel.read_text() == '{"sentinel":"keep"}\n'
    invalid_weight = json.loads(json.dumps(recipe))
    invalid_weight["profiles"]["x_min"]["controls"][0][2] = 0
    invalid_weight_recipe = work / "invalid-weight.json"
    invalid_weight_recipe.write_text(json.dumps(invalid_weight))
    refused_weight = subprocess.run([str(cli), "build", str(invalid_weight_recipe), str(sentinel)], text=True, capture_output=True)
    preserved_weight = sentinel.read_text() == '{"sentinel":"keep"}\n'

    # Distinct interior ordinate controls exercise the explicitly approximate,
    # certified graph mode while preserving monotonicity and both cap levels.
    bounded = json.loads(json.dumps(recipe))
    bounded["family"] = "nurbs_cross_bores_bounded_v1"
    bounded["profiles"]["x_min"]["controls"][2][1] = 19
    bounded["lateral_tolerance_mm"] = [1, 4]
    bounded["max_depth"] = 8
    bounded_recipe = work / "bounded.json"
    bounded_path = work / "bounded-model.json"
    bounded_recipe.write_text(json.dumps(bounded))
    bounded_build = subprocess.run([str(cli), "build", str(bounded_recipe), str(bounded_path)], text=True, capture_output=True)
    bounded_report = json.loads(bounded_build.stdout) if bounded_build.stdout.strip().startswith("{") else {}
    bounded_model = json.loads(bounded_path.read_text()) if bounded_path.exists() else {}
    bounded_profile_errors = profile_deviation(bounded, bounded_model) if bounded_model else {}
    bounded_profile_volume = volume_from_model_profiles(bounded_model, bounded) if bounded_model else None
    bounded_block_volume = block_volume(bounded_model.get("blocks", []), 10) if bounded_model else None
    bounded_source_volume = expected_volume(bounded, panels=2048) if bounded_model else None

    zero_budget = json.loads(json.dumps(bounded))
    zero_budget["max_depth"] = 0
    zero_budget_path = work / "zero-budget.json"
    zero_budget_path.write_text(json.dumps(zero_budget))
    refused_budget = subprocess.run([str(cli), "build", str(zero_budget_path), str(sentinel)], text=True, capture_output=True)
    preserved_budget = sentinel.read_text() == '{"sentinel":"keep"}\n'
    tiny_tolerance = json.loads(json.dumps(bounded))
    tiny_tolerance["lateral_tolerance_mm"] = [1, 1000000]
    tiny_recipe_path = work / "tiny-tolerance.json"
    tiny_recipe_path.write_text(json.dumps(tiny_tolerance))
    refused_tiny = subprocess.run([str(cli), "build", str(tiny_recipe_path), str(sentinel)], text=True, capture_output=True)
    preserved_tiny = sentinel.read_text() == '{"sentinel":"keep"}\n'
    return {
        "oracle": "nurbs_cross_bores_independent_de_boor_simpson_v1",
        "expected_volume_mm3": expected,
        "rust_export_block_volume_q12_mm3": volume_q12,
        "rust_export_block_volume_q16_mm3": volume_q16,
        "rust_volume_quadrature_delta_mm3": volume_q16 - volume_q12,
        "volume_error_mm3": volume_error,
        "volume_relative_error": abs(volume_error) / expected,
        "simpson_panels": 8192,
        "profile_observations_mm": profile_observations,
        "valid_model_schema": model.get("schema"),
        "valid_model_blocks": len(model.get("blocks", [])),
        "independent_block_point_classification": independent_classification(recipe, model),
        "desynchronized_ordinate_refusal_exit": refused.returncode,
        "desynchronized_ordinate_refusal_output": (refused.stderr or refused.stdout).strip(),
        "desynchronized_atomic_output_preserved": preserved_desynchronized,
        "invalid_weight_refusal_exit": refused_weight.returncode,
        "invalid_weight_refusal_output": (refused_weight.stderr or refused_weight.stdout).strip(),
        "invalid_weight_atomic_output_preserved": preserved_weight,
        "bounded_family_exit": bounded_build.returncode,
        "bounded_schema": bounded_model.get("schema"),
        "bounded_cli_audit": bounded_report,
        "bounded_profile_errors_from_rust_blocks_mm": bounded_profile_errors,
        "bounded_profile_tolerance_mm": 0.25,
        "bounded_profile_volume_from_rust_blocks_mm3": bounded_profile_volume,
        "bounded_block_jacobian_volume_q10_mm3": bounded_block_volume,
        "bounded_volume_crosscheck_error_mm3": None if bounded_profile_volume is None or bounded_block_volume is None else bounded_block_volume - bounded_profile_volume,
        "bounded_independent_source_volume_mm3": bounded_source_volume,
        "bounded_block_vs_source_volume_error_mm3": None if bounded_source_volume is None or bounded_block_volume is None else bounded_block_volume - bounded_source_volume,
        "zero_depth_refusal_exit": refused_budget.returncode,
        "zero_depth_refusal_output": (refused_budget.stderr or refused_budget.stdout).strip(),
        "zero_depth_atomic_output_preserved": preserved_budget,
        "tiny_tolerance_refusal_exit": refused_tiny.returncode,
        "tiny_tolerance_refusal_output": (refused_tiny.stderr or refused_tiny.stdout).strip(),
        "tiny_tolerance_atomic_output_preserved": preserved_tiny,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cli", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    result = run(args.cli.resolve(), args.work.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
