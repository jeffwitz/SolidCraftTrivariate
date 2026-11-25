#!/usr/bin/env python3
"""Independent exact oracle for NURBS geometry-degree space packages.

This script rebuilds one two-span strip and one crossed-bore model from their
versioned recipes. It checks degree elevation in Fraction arithmetic, physical
maps, D4 trace-node identities, and affine value/gradient reproduction. Volumes
are integrated separately from the source curves with Gauss-Legendre rules.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
from fractions import Fraction
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "tests/references/nurbs_sketch_corpus.json"
CROSS_RECIPE = ROOT / "tests/references/nurbs_cross_bores.recipe.json"


def choose(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    out = 1
    for i in range(1, k + 1):
        out = out * (n - i + 1) // i
    return out


def bernstein(degree: int, t: Fraction) -> list[Fraction]:
    return [Fraction(choose(degree, i)) * t**i * (1 - t)**(degree - i) for i in range(degree + 1)]


def bernstein_derivative(degree: int, t: Fraction) -> list[Fraction]:
    lower = bernstein(degree - 1, t)
    return [degree * ((lower[i - 1] if i else 0) - (lower[i] if i < degree else 0)) for i in range(degree + 1)]


def flat_index(index: tuple[int, int, int], degrees: tuple[int, int, int]) -> int:
    i, j, k = index
    return (i * (degrees[1] + 1) + j) * (degrees[2] + 1) + k


def eval_homogeneous(degrees: tuple[int, int, int], controls: list[list[int | Fraction]], uvw: tuple[Fraction, Fraction, Fraction]) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    bases = [bernstein(degrees[a], uvw[a]) for a in range(3)]
    out = [Fraction(0) for _ in range(4)]
    for i in range(degrees[0] + 1):
        for j in range(degrees[1] + 1):
            for k in range(degrees[2] + 1):
                factor = bases[0][i] * bases[1][j] * bases[2][k]
                p = controls[flat_index((i, j, k), degrees)]
                for c in range(4):
                    out[c] += factor * p[c]
    return tuple(out)  # type: ignore[return-value]


def physical(h: tuple[Fraction, ...]) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(h[c] / h[3] for c in range(3))  # type: ignore[return-value]


def elevate_axis(controls: list[list[Fraction]], degrees: tuple[int, int, int], axis: int) -> tuple[list[list[Fraction]], tuple[int, int, int]]:
    target = list(degrees)
    target[axis] += 1
    target_tuple = tuple(target)
    out = [[Fraction(0) for _ in range(4)] for _ in range(math.prod(d + 1 for d in target_tuple))]
    n = target_tuple[axis]
    for i in range(target_tuple[0] + 1):
        for j in range(target_tuple[1] + 1):
            for k in range(target_tuple[2] + 1):
                q = [i, j, k]
                m = q[axis]
                left = q.copy()
                right = q.copy()
                left[axis] = max(0, m - 1)
                right[axis] = min(m, n - 1)
                a = controls[flat_index(tuple(left), degrees)]
                b = controls[flat_index(tuple(right), degrees)]
                out[flat_index((i, j, k), target_tuple)] = [(m * a[c] + (n - m) * b[c]) / n for c in range(4)]
    return out, target_tuple


def elevated_expected(block: dict, target: int) -> tuple[list[list[Fraction]], int]:
    degrees = tuple(block["degrees"])
    controls = [[Fraction(v) for v in p] for p in block["controls"]]
    scale = 1
    for axis in range(3):
        while degrees[axis] < target:
            controls, degrees = elevate_axis(controls, degrees, axis)
            scale *= degrees[axis]
    return controls, scale


def gauss_legendre(order: int) -> list[tuple[float, float]]:
    result = []
    for i in range(order):
        x = math.cos(math.pi * (i + 0.75) / (order + 0.5))
        for _ in range(20):
            p0, p1 = 1.0, x
            for n in range(2, order + 1):
                p0, p1 = p1, ((2 * n - 1) * x * p1 - (n - 1) * p0) / n
            derivative = order * (x * p1 - p0) / (x * x - 1.0)
            delta = p1 / derivative
            x -= delta
            if abs(delta) < 2e-16:
                break
        p0, p1 = 1.0, x
        for n in range(2, order + 1):
            p0, p1 = p1, ((2 * n - 1) * x * p1 - (n - 1) * p0) / n
        derivative = order * (x * p1 - p0) / (x * x - 1.0)
        weight = 2.0 / ((1.0 - x * x) * derivative * derivative)
        result.append(((x + 1.0) / 2.0, weight / 2.0))
    return result


def de_boor(degree: int, knots: list[int], controls: list[list[int | float]], t: float) -> list[float]:
    n = len(controls) - 1
    span = n if t == knots[n + 1] else next(i for i in range(degree, n + 1) if knots[i] <= t < knots[i + 1])
    work = [list(map(float, controls[j])) for j in range(span - degree, span + 1)]
    for level in range(1, degree + 1):
        for j in range(degree, level - 1, -1):
            i = span - degree + j
            denominator = knots[i + degree + 1 - level] - knots[i]
            alpha = (t - knots[i]) / denominator if denominator else 0.0
            work[j] = [(1.0 - alpha) * work[j - 1][c] + alpha * work[j][c] for c in range(3)]
    return work[degree]


def de_boor_derivative(degree: int, knots: list[int], controls: list[list[int | float]], t: float) -> list[float]:
    derivative_controls = []
    for i in range(len(controls) - 1):
        denominator = knots[i + degree + 1] - knots[i + 1]
        derivative_controls.append([degree * (controls[i + 1][c] - controls[i][c]) / denominator if denominator else 0.0 for c in range(3)])
    return de_boor(degree - 1, knots[1:-1], derivative_controls, t)


def curve_xy(curve: dict, t: float) -> tuple[float, float, float, float]:
    p, knots, controls = curve["degree"], curve["knots"], curve["controls"]
    h = de_boor(p, knots, controls, t)
    dh = de_boor_derivative(p, knots, controls, t)
    x, y, w = (h[c] / h[2] for c in range(3))
    dx = (dh[0] * h[2] - h[0] * dh[2]) / h[2] ** 2
    dy = (dh[1] * h[2] - h[1] * dh[2]) / h[2] ** 2
    return x, y, dx, dy


def knot_intervals(curve: dict) -> list[tuple[int, int]]:
    p, u = curve["degree"], curve["knots"]
    n = len(curve["controls"]) - 1
    return [(u[i], u[i + 1]) for i in range(p, n + 1) if u[i] < u[i + 1]]


def expected_strip_volume(recipe: dict, order: int = 64) -> float:
    lower, upper, height = recipe["lower"], recipe["upper"], recipe["height_mm"]
    assert knot_intervals(lower) == knot_intervals(upper)
    total = 0.0
    for lo, hi in knot_intervals(lower):
        for q, weight in gauss_legendre(order):
            t = lo + (hi - lo) * q
            x0, y0, dx0, _ = curve_xy(lower, t)
            x1, y1, dx1, _ = curve_xy(upper, t)
            assert abs(x0 - x1) < 1e-11 and abs(dx0 - dx1) < 1e-9
            total += weight * (hi - lo) * (y1 - y0) * dx0
    return height * total


def expected_cross_volume(recipe: dict, order: int = 64) -> float:
    profiles = recipe["profiles"]
    base = math.prod(hi - lo for lo, hi in recipe["bounds_mm"])
    representative = profiles["x_min"]
    total_removed = 0.0
    for lo, hi in knot_intervals(representative):
        for q, weight in gauss_legendre(order):
            t = lo + (hi - lo) * q
            _, y, _, dy = curve_xy(representative, t)
            xmin, _, _, _ = curve_xy(profiles["x_min"], t)
            xmax, _, _, _ = curve_xy(profiles["x_max"], t)
            zmin, _, _, _ = curve_xy(profiles["z_min"], t)
            zmax, _, _, _ = curve_xy(profiles["z_max"], t)
            dx, dz = xmax - xmin, zmax - zmin
            cross_section_removed = (recipe["bounds_mm"][0][1] - recipe["bounds_mm"][0][0]) * dz
            cross_section_removed += (recipe["bounds_mm"][2][1] - recipe["bounds_mm"][2][0]) * dx
            cross_section_removed -= dx * dz
            total_removed += weight * (hi - lo) * cross_section_removed * dy
    return base - total_removed


def side(side_json: list[str]) -> tuple[str, int, int]:
    block, name = side_json
    return block, "uvw".index(name[0]), int(name[1])


def face_axes(axis: int) -> tuple[int, int]:
    return ((1, 2), (0, 2), (0, 1))[axis]


def map_from_json(block: dict, uvw: tuple[Fraction, Fraction, Fraction]) -> tuple[Fraction, Fraction, Fraction]:
    h = eval_homogeneous(tuple(block["degrees"]), block["controls"], uvw)
    return physical(h)


def package_checks(model: dict, package: dict) -> dict:
    assert package["source_model"] == model
    source_blocks = {b["id"]: b for b in model["blocks"]}
    package_blocks = {b["id"]: b for b in package["blocks"]}
    assert source_blocks.keys() == package_blocks.keys()
    target = package["tensor_degree"]
    assert all(b["degrees"] == [target] * 3 for b in package_blocks.values())
    scalar_ids = [d for b in package_blocks.values() for d in b["dofs"]]
    assert all(len(b["controls"]) == (target + 1) ** 3 and len(b["dofs"]) == (target + 1) ** 3 for b in package_blocks.values())
    unique_ids = set(scalar_ids)
    assert package["scalar_ndof"] == len(unique_ids)
    assert unique_ids == set(range(package["scalar_ndof"]))

    elevation_checks = 0
    source_package_map_samples = 0
    affine_checks = 0
    dof_positions: dict[int, tuple[Fraction, Fraction, Fraction]] = {}
    affine_matrix = ((3, -2, 1), (1, 4, -3), (-2, 1, 2))
    affine_offset = (2, -1, 5)
    affine_values: dict[int, tuple[Fraction, Fraction, Fraction]] = {}
    for block_id, source in source_blocks.items():
        package_block = package_blocks[block_id]
        elevated, gauge = elevated_expected(source, target)
        for expected, actual in zip(elevated, package_block["controls"]):
            if any(Fraction(actual[c]) != expected[c] * gauge for c in range(4)):
                raise AssertionError(f"{block_id}: exact degree elevation differs")
            elevation_checks += 1
        for dof, h in zip(package_block["dofs"], package_block["controls"]):
            point = tuple(Fraction(h[c], h[3]) for c in range(3))
            if dof in dof_positions and dof_positions[dof] != point:
                raise AssertionError(f"shared scalar DOF {dof} has inconsistent physical nodes")
            dof_positions[dof] = point
            values = tuple(Fraction(affine_offset[r]) + sum(Fraction(affine_matrix[r][c]) * point[c] for c in range(3)) for r in range(3))
            if dof in affine_values and affine_values[dof] != values:
                raise AssertionError(f"shared scalar DOF {dof} has inconsistent affine data")
            affine_values[dof] = values

        source_degrees = tuple(source["degrees"])
        target_degrees = tuple(package_block["degrees"])
        for t0 in (Fraction(0), Fraction(1, 3), Fraction(2, 3), Fraction(1)):
            for t1 in (Fraction(0), Fraction(1, 3), Fraction(2, 3), Fraction(1)):
                for t2 in (Fraction(0), Fraction(1, 3), Fraction(2, 3), Fraction(1)):
                    uvw = (t0, t1, t2)
                    x_source = physical(eval_homogeneous(source_degrees, source["controls"], uvw))
                    x_package = physical(eval_homogeneous(target_degrees, package_block["controls"], uvw))
                    if x_source != x_package:
                        raise AssertionError(f"{block_id}: source/elevated map mismatch at {uvw}")
                    source_package_map_samples += 1

        # A nontrivial spatial-affine Qp field is interpolated with rational
        # geometry weights. Its value and parametric gradient must agree exactly.
        controls = package_block["controls"]
        degrees = target_degrees
        dofs = package_block["dofs"]
        for uvw in product((Fraction(0), Fraction(1, 3), Fraction(2, 3), Fraction(1)), repeat=3):
            bases = [bernstein(target, uvw[a]) for a in range(3)]
            derivs = [bernstein_derivative(target, uvw[a]) for a in range(3)]
            H = [Fraction(0)] * 4
            dH = [[Fraction(0)] * 4 for _ in range(3)]
            S = [Fraction(0)] * 3
            dS = [[Fraction(0)] * 3 for _ in range(3)]
            for i in range(target + 1):
                for j in range(target + 1):
                    for k in range(target + 1):
                        idx = flat_index((i, j, k), degrees)
                        b = (bases[0][i], bases[1][j], bases[2][k])
                        h = controls[idx]
                        dof_value = affine_values[dofs[idx]]
                        for c in range(4):
                            H[c] += b[0] * b[1] * b[2] * h[c]
                            for axis in range(3):
                                factors = list(b)
                                factors[axis] = derivs[axis][(i, j, k)[axis]]
                                dH[axis][c] += factors[0] * factors[1] * factors[2] * h[c]
                        for c in range(3):
                            wh = Fraction(h[3])
                            S[c] += b[0] * b[1] * b[2] * wh * dof_value[c]
                            for axis in range(3):
                                factors = list(b)
                                factors[axis] = derivs[axis][(i, j, k)[axis]]
                                dS[axis][c] += factors[0] * factors[1] * factors[2] * wh * dof_value[c]
            point = tuple(H[c] / H[3] for c in range(3))
            field = tuple(S[c] / H[3] for c in range(3))
            expected_field = tuple(Fraction(affine_offset[r]) + sum(Fraction(affine_matrix[r][c]) * point[c] for c in range(3)) for r in range(3))
            if field != expected_field:
                raise AssertionError(f"{block_id}: affine field value is not reproduced")
            for axis in range(3):
                geom_derivative = tuple((dH[axis][c] - point[c] * dH[axis][3]) / H[3] for c in range(3))
                field_derivative = tuple((dS[axis][r] - field[r] * dH[axis][3]) / H[3] for r in range(3))
                expected_derivative = tuple(sum(Fraction(affine_matrix[r][c]) * geom_derivative[c] for c in range(3)) for r in range(3))
                if field_derivative != expected_derivative:
                    raise AssertionError(f"{block_id}: affine field gradient is not reproduced")
            affine_checks += 1

    interfaces = model["topology"]["interfaces"]
    trace_nodes, trace_samples = 0, 0
    for interface in interfaces:
        first_id, first_axis, first_upper = side(interface["first"])
        second_id, second_axis, second_upper = side(interface["second"])
        first, second = package_blocks[first_id], package_blocks[second_id]
        orientation = interface["orientation"]
        first_free = face_axes(first_axis)
        second_free = face_axes(second_axis)
        for a in range(target + 1):
            for b in range(target + 1):
                af, bf = (target - a if orientation["flip_first"] else a), (target - b if orientation["flip_second"] else b)
                sa, sb = (bf, af) if orientation["swap"] else (af, bf)
                qa, qb = [0, 0, 0], [0, 0, 0]
                qa[first_axis], qb[second_axis] = first_upper * target, second_upper * target
                qa[first_free[0]], qa[first_free[1]] = a, b
                qb[second_free[0]], qb[second_free[1]] = sa, sb
                ia, ib = flat_index(tuple(qa), (target,) * 3), flat_index(tuple(qb), (target,) * 3)
                if first["dofs"][ia] != second["dofs"][ib]:
                    raise AssertionError(f"D4 trace DOF mismatch on {first_id}/{second_id}")
                pa = tuple(Fraction(first["controls"][ia][c], first["controls"][ia][3]) for c in range(3))
                pb = tuple(Fraction(second["controls"][ib][c], second["controls"][ib][3]) for c in range(3))
                if pa != pb:
                    raise AssertionError("shared trace DOF has mismatched physical control point")
                trace_nodes += 1
        for s, t in product((Fraction(0), Fraction(1, 3), Fraction(2, 3), Fraction(1)), repeat=2):
            first_uvw, second_uvw = [Fraction(0)] * 3, [Fraction(0)] * 3
            first_uvw[first_axis] = Fraction(first_upper)
            second_uvw[second_axis] = Fraction(second_upper)
            first_uvw[first_free[0]], first_uvw[first_free[1]] = s, t
            a, b = (1 - s if orientation["flip_first"] else s), (1 - t if orientation["flip_second"] else t)
            second_uvw[second_free[0]], second_uvw[second_free[1]] = (b, a) if orientation["swap"] else (a, b)
            if map_from_json(source_blocks[first_id], tuple(first_uvw)) != map_from_json(source_blocks[second_id], tuple(second_uvw)):
                raise AssertionError("source geometry trace mismatch under recorded D4 map")
            trace_samples += 1

    return {"source_blocks": len(source_blocks), "elevated_blocks": len(package_blocks),
            "target_tensor_degree": target, "scalar_ndof": package["scalar_ndof"],
            "shared_interfaces": len(interfaces), "exact_elevated_control_checks": elevation_checks,
            "source_vs_package_map_samples": source_package_map_samples,
            "exact_trace_dof_nodes": trace_nodes, "source_trace_map_samples": trace_samples,
            "affine_value_and_gradient_checks": affine_checks}


def run_cli(cli: Path, command: list[str]) -> dict:
    process = subprocess.run([str(cli), *command], cwd=ROOT, capture_output=True, text=True)
    if process.returncode:
        raise RuntimeError(f"CLI {' '.join(command)} failed: {process.stderr or process.stdout}")
    try:
        return json.loads(process.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError(f"invalid CLI JSON for {' '.join(command)}: {process.stdout}") from error


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cli", type=Path, required=True)
    ap.add_argument("--work", type=Path, default=ROOT / "build/nurbs/nurbs-space-oracle")
    ap.add_argument("--output", type=Path, default=ROOT / "docs/results/nurbs-space-20261125.json")
    args = ap.parse_args()
    cli, work = args.cli.resolve(), args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    strip_case = next(x for x in json.loads(CORPUS.read_text()) if x["id"] == "two_spans_variable_gap")
    recipes = {"strip": strip_case["recipe"], "cross": json.loads(CROSS_RECIPE.read_text())}
    rows = {}
    for name, recipe in recipes.items():
        recipe_path, model_path = work / f"{name}.recipe.json", work / f"{name}.model.json"
        package_path, analysis_path = work / f"{name}.space.json", work / f"{name}.analysis.json"
        recipe_path.write_text(json.dumps(recipe, separators=(",", ":")) + "\n")
        build_report = run_cli(cli, ["build", str(recipe_path), str(model_path)])
        model = json.loads(model_path.read_text())
        package_report = run_cli(cli, ["export-nurbs-space", str(model_path), str(package_path)])
        package = json.loads(package_path.read_text())
        integrated_by_order = {}
        analysis = None
        for order in (8, 16):
            target_path = analysis_path if order == 8 else work / f"{name}.analysis{order}.json"
            run_cli(cli, ["export-analysis", str(model_path), str(order), str(target_path)])
            current = json.loads(target_path.read_text())
            integrated_by_order[str(order)] = current["integrated_volume_mm3"]
            if order == 8:
                analysis = current
        assert analysis is not None
        space_metrics = package_checks(model, package)
        expected_volume = expected_strip_volume(recipe) if name == "strip" else expected_cross_volume(recipe)
        integrated = analysis["integrated_volume_mm3"]
        relative_error = abs(integrated - expected_volume) / max(abs(expected_volume), 1.0)
        converged_error = abs(integrated_by_order["16"] - expected_volume) / max(abs(expected_volume), 1.0)
        if converged_error > 2e-9:
            raise AssertionError(f"{name}: order-16 integrated volume differs from independent source-curve quadrature by {converged_error}")
        if analysis["source_model"] != model or analysis["blocks"] != package["blocks"]:
            raise AssertionError(f"{name}: analysis package changed the audited source or space");
        rows[name] = {"schema": model["schema"], "build_audit": build_report,
                      "export_scalar_ndof": package_report["scalar_ndof"],
                      "space": space_metrics,
                      "independent_volume_mm3": expected_volume,
                      "package_gauss_order": analysis["quadrature_order"],
                      "package_integrated_volume_by_order_mm3": integrated_by_order,
                      "volume_relative_error_by_order": {
                          order: abs(value - expected_volume) / max(abs(expected_volume), 1.0)
                          for order, value in integrated_by_order.items()},
                      "order16_relative_error": converged_error,
                      "files": {"recipe": str(recipe_path.relative_to(ROOT)),
                                "model": str(model_path.relative_to(ROOT)),
                                "space": str(package_path.relative_to(ROOT)),
                                "analysis": str(analysis_path.relative_to(ROOT))}}
    report = {"schema": "nurbs_space_oracle_report_v1",
              "protocol": "Exact Fraction tensor Bernstein evaluation/elevation and D4 trace identity; rational affine field/gradient checks; independent Gauss-Legendre integration of source De Boor curves.",
              "environment": {"cli": str(cli), "exact_arithmetic": "fractions.Fraction", "quadrature_order_independent": 64, "package_quadrature_order": 8},
              "models": rows,
              "companion_patch_test": None,
              "limits": ["This qualifies the two frozen strip/cross recipes only; it is not a global approximation theorem for arbitrary NURBS partitions.",
                            "Patch-test equilibrium residual was not assembled in this oracle."]}
    patch_path = ROOT / "docs/results/nurbs-patch-equilibrium-20261125.json"
    if patch_path.exists():
        patch = json.loads(patch_path.read_text())
        report["companion_patch_test"] = {
            "report": str(patch_path.relative_to(ROOT)),
            "status": patch.get("status"),
            "method": patch.get("method"),
            "orders": patch.get("orders"),
            "relative_residual_threshold": patch.get("q16_threshold_relative_linf"),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "passed", "models": list(rows), "output": str(args.output.relative_to(ROOT))}))


if __name__ == "__main__":
    main()
