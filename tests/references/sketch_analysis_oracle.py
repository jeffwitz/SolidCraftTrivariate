#!/usr/bin/env python3
"""Independent PA7 oracle for sketch Q2 packages and the fixed bracket baseline.

The Rust CLI constructs, audits, numbers, and assembles each space. This script
replays the published package using exact Python fractions: it checks every
whole-face orientation, shared Q2 control, partition-of-unity patch, affine
field and physical gradient. Conditioning is a sampled diagnostic on the same
3x3x3 tensor grid in every block; it is not a continuous quality certificate.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import subprocess
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


GRID = (F(0), F(1, 2), F(1))
GRADIENT = (
    (F(13, 1000), F(-4, 1000), F(2, 1000)),
    (F(1, 1000), F(9, 1000), F(-3, 1000)),
    (F(-2, 1000), F(5, 1000), F(6, 1000)),
)
OFFSET = (F(7, 10), F(-1, 5), F(3, 10))


def run(cli: Path, args: list[str]) -> tuple[str, float]:
    start = time.perf_counter()
    completed = subprocess.run(
        [str(cli), *args], check=True, text=True, capture_output=True
    )
    return completed.stdout.strip(), time.perf_counter() - start


def bernstein(t: F) -> tuple[tuple[F, ...], tuple[F, ...]]:
    return (
        ((1 - t) ** 2, 2 * t * (1 - t), t**2),
        (-2 * (1 - t), 2 - 4 * t, 2 * t),
    )


def bernstein_values(degree: int, t: F) -> tuple[F, ...]:
    """Exact Bernstein basis values for arbitrary small nonnegative degree."""
    return tuple(
        F(math.comb(degree, i)) * t**i * (1 - t) ** (degree - i)
        for i in range(degree + 1)
    )


def position(control: list[int]) -> tuple[F, F, F]:
    return tuple(F(control[a], control[3]) for a in range(3))


def tensor_eval(block: dict, uvw: tuple[F, F, F]) -> tuple[tuple[F, ...], list[list[F]], tuple[F, ...], list[list[F]]]:
    basis = [bernstein(t) for t in uvw]
    homogeneous = [F(0) for _ in range(4)]
    derivatives = [[F(0) for _ in range(4)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                n = (i * 3 + j) * 3 + k
                c = block["controls"][n]
                indices = (i, j, k)
                value = math.prod(basis[a][0][indices[a]] for a in range(3))
                for component in range(4):
                    homogeneous[component] += value * c[component]
                for axis in range(3):
                    deriv = basis[axis][1][indices[axis]]
                    for a in range(3):
                        if a != axis:
                            deriv *= basis[a][0][indices[a]]
                    for component in range(4):
                        derivatives[axis][component] += deriv * c[component]
    weight = homogeneous[3]
    physical = tuple(homogeneous[a] / weight for a in range(3))
    jac = [
        [
            (derivatives[axis][a] * weight - homogeneous[a] * derivatives[axis][3])
            / (weight * weight)
            for a in range(3)
        ]
        for axis in range(3)
    ]
    return physical, jac, tuple(homogeneous), derivatives


def affine_reproduction(block: dict) -> tuple[float, float]:
    max_position = F(0)
    max_gradient = F(0)
    for uvw in __import__("itertools").product(GRID, repeat=3):
        p, jac, homogeneous, derivatives = tensor_eval(block, uvw)
        bases = [bernstein(t) for t in uvw]
        for a in range(3):
            numerator = F(0)
            d_numerator = [F(0) for _ in range(3)]
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        n = (i * 3 + j) * 3 + k
                        control = block["controls"][n]
                        ix = (i, j, k)
                        bval = math.prod(bases[q][0][ix[q]] for q in range(3))
                        field_value = OFFSET[a] + sum(
                            GRADIENT[a][q] * position(control)[q] for q in range(3)
                        )
                        numerator += bval * control[3] * field_value
                        for axis in range(3):
                            dval = bases[axis][1][ix[axis]]
                            for q in range(3):
                                if q != axis:
                                    dval *= bases[q][0][ix[q]]
                            d_numerator[axis] += dval * control[3] * field_value
            field_value = numerator / homogeneous[3]
            expected_value = OFFSET[a] + sum(GRADIENT[a][q] * p[q] for q in range(3))
            max_position = max(max_position, abs(field_value - expected_value))
            for axis in range(3):
                field_derivative = (
                    d_numerator[axis] * homogeneous[3]
                    - numerator * derivatives[axis][3]
                ) / homogeneous[3] ** 2
                expected = sum(GRADIENT[a][b] * jac[axis][b] for b in range(3))
                max_gradient = max(max_gradient, abs(field_derivative - expected))
    return float(max_position), float(max_gradient)


def side_indices(side: str, a: int, b: int) -> tuple[int, int, int]:
    axes = {"u": 0, "v": 1, "w": 2}
    axis = axes[side[0]]
    fixed = 2 if side[1] == "1" else 0
    free = {0: (1, 2), 1: (0, 2), 2: (0, 1)}[axis]
    index = [0, 0, 0]
    index[axis] = fixed
    index[free[0]] = a
    index[free[1]] = b
    return tuple(index)


def oriented_index(orientation: dict, a: int, b: int) -> tuple[int, int]:
    x = 2 - a if orientation["flip_first"] else a
    y = 2 - b if orientation["flip_second"] else b
    return (y, x) if orientation["swap"] else (x, y)


def flat(index: tuple[int, int, int]) -> int:
    return (index[0] * 3 + index[1]) * 3 + index[2]


def check_interfaces(package: dict) -> tuple[int, Counter]:
    blocks = {b["id"]: b for b in package["analysis_space"]["blocks"]}
    orientations: Counter = Counter()
    checked = 0
    for face in package["topology"]["interfaces"]:
        (first_id, first_side) = face["first"]
        (second_id, second_side) = face["second"]
        orient = face["orientation"]
        orientations[tuple(orient[k] for k in ("swap", "flip_first", "flip_second"))] += 1
        first, second = blocks[first_id], blocks[second_id]
        for a in range(3):
            for b in range(3):
                ix = side_indices(first_side, a, b)
                sa, sb = oriented_index(orient, a, b)
                iy = side_indices(second_side, sa, sb)
                fi, si = flat(ix), flat(iy)
                if first["dofs"][fi] != second["dofs"][si]:
                    raise AssertionError((first_id, first_side, second_id, second_side, ix, iy, "DOF mismatch"))
                if position(first["controls"][fi]) != position(second["controls"][si]):
                    raise AssertionError((first_id, first_side, second_id, second_side, ix, iy, "physical trace mismatch"))
                checked += 1
    return checked, orientations


def jacobian_quality(block: dict, uvw: tuple[F, F, F]) -> tuple[float, float]:
    _, jac, _, _ = tensor_eval(block, uvw)
    a = [[sum(float(jac[i][k]) * float(jac[j][k]) for k in range(3)) for j in range(3)] for i in range(3)]
    for _ in range(40):
        changed = False
        for p, q in ((0, 1), (0, 2), (1, 2)):
            if abs(float(a[p][q])) < 1e-28:
                continue
            angle = 0.5 * math.atan2(2 * float(a[p][q]), float(a[q][q] - a[p][p]))
            c, s = math.cos(angle), math.sin(angle)
            app, aqq, apq = float(a[p][p]), float(a[q][q]), float(a[p][q])
            a[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq
            a[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq
            a[p][q] = a[q][p] = 0.0
            for k in range(3):
                if k not in (p, q):
                    apk, aqk = float(a[p][k]), float(a[q][k])
                    a[p][k] = a[k][p] = c * apk - s * aqk
                    a[q][k] = a[k][q] = s * apk + c * aqk
            changed = True
        if not changed:
            break
    eig = sorted(max(0.0, float(a[i][i])) for i in range(3))
    cond = math.sqrt(eig[2] / eig[0]) if eig[0] > 0 else math.inf
    columns = [[float(jac[i][k]) for i in range(3)] for k in range(3)]
    lengths = [math.sqrt(sum(x * x for x in col)) for col in columns]
    det = float(
        jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1])
        - jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0])
        + jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0])
    )
    scaled = det / math.prod(lengths)
    if det <= 0 or not math.isfinite(cond):
        raise AssertionError("nonpositive or singular sampled Q2 geometry")
    return cond, scaled


def inspect_package(package: dict, expected_blocks: int, expected_interfaces: int) -> dict:
    if package["schema"] not in (
        "solidcraft_sketch_analysis_package_v1",
        "solidcraft_solver_package_v1",
    ):
        raise AssertionError("unexpected sketch analysis package schema")
    analysis = package["analysis_space"]
    blocks = analysis["blocks"]
    if len(blocks) != expected_blocks:
        raise AssertionError((len(blocks), expected_blocks))
    if len(package["topology"]["interfaces"]) != expected_interfaces:
        raise AssertionError("interface count mismatch")
    dofs = set()
    max_position = max_gradient = 0.0
    max_condition = 0.0
    min_scaled_det = math.inf
    for block in blocks:
        if block["degrees"] != [2, 2, 2] or len(block["controls"]) != 27 or len(block["dofs"]) != 27:
            raise AssertionError((block["id"], "not a 27-node Q2 block"))
        if any(weight <= 0 for _, _, _, weight in block["controls"]):
            raise AssertionError((block["id"], "nonpositive Q2 weight"))
        dofs.update(block["dofs"])
        p_err, g_err = affine_reproduction(block)
        max_position = max(max_position, p_err)
        max_gradient = max(max_gradient, g_err)
        for uvw in __import__("itertools").product(GRID, repeat=3):
            cond, scaled_det = jacobian_quality(block, uvw)
            max_condition = max(max_condition, cond)
            min_scaled_det = min(min_scaled_det, scaled_det)
    if max(dofs) + 1 != analysis["scalar_ndof"] or len(dofs) != analysis["scalar_ndof"]:
        raise AssertionError("global DOF numbering has holes")
    interface_nodes, orientations = check_interfaces(package)
    if interface_nodes != 9 * expected_interfaces:
        raise AssertionError("not every oriented face-node pair was checked")
    boundary = {name: set(values) for name, values in analysis["boundary_dofs"].items()}
    expected_boundary: dict[str, set[int]] = defaultdict(set)
    for face in package["design_faces"]:
        block = next(b for b in blocks if b["id"] == face["block"])
        for a in range(3):
            for b in range(3):
                ix = side_indices(face["side"], a, b)
                expected_boundary[face["feature"]].add(block["dofs"][flat(ix)])
    if {k: v for k, v in expected_boundary.items()} != boundary:
        raise AssertionError("design face selection differs from Q2 boundary DOFs")
    diagnostics = package.get("analysis_diagnostics")
    rigid = None
    production_energy = None
    if diagnostics is not None:
        rigid = diagnostics["mounted_q2"]["rigid_mode_relative_residuals"]
        if len(rigid) != 6 or not all(math.isfinite(x) for x in rigid):
            raise AssertionError("six rigid-mode diagnostics are missing or nonfinite")
        production_energy = diagnostics["mounted_q2"]["affine_field_energy"]
        if production_energy["relative_error"] > 1e-10:
            raise AssertionError("production Q2 affine-strain energy mismatch")
    if max_position > 1e-12 or max_gradient > 1e-11:
        raise AssertionError((max_position, max_gradient, "affine field reproduction failed"))
    return {
        "blocks": len(blocks),
        "interfaces": expected_interfaces,
        "oriented_face_node_pairs_checked": interface_nodes,
        "orientation_counts": {str(k): v for k, v in orientations.items()},
        "scalar_dofs": analysis["scalar_ndof"],
        "vector_dofs": analysis["vector_ndof"],
        "affine_max_position_error_mm": max_position,
        "affine_max_parametric_gradient_error": max_gradient,
        "sampled_grid_per_block": [3, 3, 3],
        "sampled_evaluations": len(blocks) * 27,
        "sampled_max_condition_number": max_condition,
        "sampled_min_scaled_determinant": min_scaled_det,
        "production_affine_energy": production_energy,
        "production_six_rigid_mode_relative_residuals": rigid,
        "production_assembly_diagnostics": diagnostics is not None,
        "selection_feature_count": len(boundary),
    }


def check_steinmetz_poles(package: dict) -> dict:
    model = package["source_model"]
    if model.get("schema") != "orthogonal_sketch_features_v1" or len(model.get("bores", [])) != 2:
        return {"status": "not_applicable"}
    axes = {"x": 0, "y": 1, "z": 2}
    first, second = model["bores"]
    a0, a1 = axes[first["axis"]], axes[second["axis"]]
    if a0 == a1:
        return {"status": "not_applicable"}
    normal = next(axis for axis in range(3) if axis not in (a0, a1))
    center = first["center_mm"]
    if first["center_mm"] != second["center_mm"] or first["radius_mm"] != second["radius_mm"]:
        return {"status": "not_applicable"}
    radius = first["radius_mm"]
    poles = []
    for sign in (-1, 1):
        point = list(center)
        point[normal] += sign * radius
        poles.append(tuple(F(v) for v in point))
    witnesses: dict[tuple[F, ...], list[dict]] = {pole: [] for pole in poles}
    for block in package["analysis_space"]["blocks"]:
        for uvw in __import__("itertools").product(GRID, repeat=3):
            point = tensor_eval(block, uvw)[0]
            if point in witnesses:
                witnesses[point].append({"block": block["id"], "uvw": [str(v) for v in uvw]})
    if any(not matches for matches in witnesses.values()):
        raise AssertionError({"Steinmetz pole not represented on the common Q2 sample grid": witnesses})
    return {
        "status": "both_axis_intersection_poles_represented",
        "physical_points_mm": [[str(x) for x in pole] for pole in poles],
        "witnesses": {str(pole): matches[:4] for pole, matches in witnesses.items()},
    }


def exact_wall_audit(package: dict) -> dict:
    model = package["source_model"]
    axis_index = {"x": 0, "y": 1, "z": 2}
    bore_data = {
        bore["id"]: (axis_index[bore["axis"]], bore["center_mm"], bore["radius_mm"])
        for bore in model["bores"]
    }
    blocks = {block["id"]: block for block in model["blocks"]}
    counts: Counter = Counter()
    checked_nodes = 0
    for face in package["design_faces"]:
        feature = face["feature"]
        if not feature.startswith("bore/") or not feature.endswith("/wall"):
            continue
        bore_id = feature[len("bore/") : -len("/wall")]
        axis, center, radius = bore_data[bore_id]
        block = blocks[face["block"]]
        side = face["side"]
        param_axis = {"u": 0, "v": 1, "w": 2}[side[0]]
        fixed = block["degrees"][param_axis] if side[1] == "1" else 0
        free_axes = {0: (1, 2), 1: (0, 2), 2: (0, 1)}[param_axis]
        degrees = [block["degrees"][free_axes[0]], block["degrees"][free_axes[1]]]
        samples = [2 * degrees[0] + 1, 2 * degrees[1] + 1]
        transverse = [coordinate for coordinate in range(3) if coordinate != axis]
        for ia in range(samples[0]):
            a = F(ia, samples[0] - 1)
            ba = bernstein_values(degrees[0], a)
            for ib in range(samples[1]):
                b = F(ib, samples[1] - 1)
                bb = bernstein_values(degrees[1], b)
                homogeneous = [F(0) for _ in range(4)]
                for i, wa in enumerate(ba):
                    for j, wb in enumerate(bb):
                        index = [0, 0, 0]
                        index[param_axis] = fixed
                        index[free_axes[0]] = i
                        index[free_axes[1]] = j
                        local = (
                            index[0] * (block["degrees"][1] + 1) + index[1]
                        ) * (block["degrees"][2] + 1) + index[2]
                        control = block["controls"][local]
                        for component in range(4):
                            homogeneous[component] += wa * wb * control[component]
                if homogeneous[3] <= 0:
                    raise AssertionError((face, "nonpositive boundary weight"))
                delta = [homogeneous[q] - center[q] * homogeneous[3] for q in transverse]
                residual = sum(d * d for d in delta) - radius**2 * homogeneous[3] ** 2
                if residual != 0:
                    raise AssertionError((face, "stored boundary is not the declared cylinder", residual))
                checked_nodes += 1
        counts[bore_id] += 1
    if set(counts) != set(bore_data) or not all(counts.values()):
        raise AssertionError(("cylinder wall selection coverage mismatch", counts, bore_data))
    return {
        "method": "exact Fraction tensor-unisolvent checks of the squared cylinder polynomial",
        "faces_by_bore": dict(counts),
        "unisolvent_nodes_checked": checked_nodes,
        "continuous_identity": "proved: degree (2p,2q) polynomial vanishes on (2p+1)x(2q+1) tensor grid",
    }


def analytic_target(build_report: dict, family: str) -> tuple[float, str]:
    volume = build_report.get("volume_exact")
    if volume is None:
        rational, pi_coefficient = F(34080), F(-616)
        source = "frozen filleted-bracket analytical target"
    else:
        if "twice_rational_part_mm3" in volume:
            rational = F(volume["twice_rational_part_mm3"], 2)
        else:
            rational = F(*volume["rational_part_mm3"])
        pi_coefficient = F(*volume["pi_coefficient_mm3"])
        source = "exact rational and pi coefficients from the cold audited build"
    expected = {
        "orthogonal_pa5": (F(34080), F(-616)),
        "fixed_bracket": (F(34080), F(-616)),
        "steinmetz_cross": (F(72144), F(-810)),
        "steinmetz_permuted": (F(72144), F(-630)),
    }[family]
    if (rational, pi_coefficient) != expected:
        raise AssertionError((family, rational, pi_coefficient, expected))
    # Fraction.from_float records the actual binary64 value used only for the
    # comparison; the exact rational/pi identity above remains in the report.
    target = float(rational + pi_coefficient * F.from_float(math.pi))
    return target, source


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, default=Path("docs/results/sketch-analysis-20261125.json"))
    args = parser.parse_args()
    root = args.root.resolve()
    cli = args.cli.resolve()
    cases = [
        ("orthogonal_pa5", root / "build/pa5/bracket.json", "orthogonal", 48, 74),
        ("fixed_bracket", root / "tests/references/filleted_bracket.json", "bracket", 36, 55),
        ("steinmetz_cross", root / "tests/references/steinmetz_box.recipe.json", "steinmetz", 40, 68),
        ("steinmetz_permuted", root / "tests/references/steinmetz_box_permuted.recipe.json", "steinmetz-permuted", 40, 68),
    ]
    results = {}
    for name, recipe, slug, nblocks, ninterfaces in cases:
        model_path = root / f"build/pa7/{slug}.model.json"
        package_path = root / f"build/pa7/{slug}.analysis.json"
        package16_path = root / f"build/pa7/{slug}.analysis16.json"
        build_stdout, build_seconds = run(cli, ["build", str(recipe), str(model_path)])
        build_report = json.loads(build_stdout)
        audit_stdout, audit_seconds = run(cli, ["audit", str(model_path)])
        audit_report = json.loads(audit_stdout)
        export_stdout, export_seconds = run(
            cli, ["export-sketch-analysis", str(model_path), "8", str(package_path)]
        )
        export_report = json.loads(export_stdout)
        export16_stdout, export16_seconds = run(
            cli, ["export-sketch-analysis", str(model_path), "16", str(package16_path)]
        )
        export16_report = json.loads(export16_stdout)
        if build_report.get("status") != "valid" or audit_report.get("status") != "valid":
            raise AssertionError((name, build_report, audit_report))
        if export_report.get("status") != "exported" or export16_report.get("status") != "exported":
            raise AssertionError((name, export_report, export16_report))
        package8 = json.loads(package_path.read_text())
        package = json.loads(package16_path.read_text())
        check = inspect_package(package, nblocks, ninterfaces)
        ideal_volume, target_source = analytic_target(build_report, name)
        volume8 = package8["measures"]["integrated_volume_mm3"]
        volume16 = package["measures"]["integrated_volume_mm3"]
        check["volume"] = {
            "target_formula": "72144-810*pi mm^3" if name == "steinmetz_cross" else "34080-616*pi mm^3",
            "target_source": target_source,
            "target_mm3": float(ideal_volume),
            "q8_mm3": volume8,
            "q16_mm3": volume16,
            "q8_to_q16_delta_mm3": volume16 - volume8,
            "q16_absolute_error_mm3": volume16 - float(ideal_volume),
            "q16_relative_error": (volume16 - float(ideal_volume)) / float(ideal_volume),
            "scope": "quadrature comparison only; it does not prove exact volume of the represented maps",
        }
        if name == "steinmetz_cross":
            check["steinmetz_poles"] = check_steinmetz_poles(package)
            check["cylindrical_wall_identity"] = exact_wall_audit(package)
        elif name == "steinmetz_permuted":
            check["steinmetz_poles"] = check_steinmetz_poles(package)
            check["cylindrical_wall_identity"] = exact_wall_audit(package)
        results[name] = {
            "recipe": str(recipe.relative_to(root)),
            "model": str(model_path.relative_to(root)),
            "package": str(package_path.relative_to(root)),
            "package_order16": str(package16_path.relative_to(root)),
            "build_seconds": build_seconds,
            "cold_audit_seconds": audit_seconds,
            "production_export_seconds": export_seconds,
            "production_export_order16_seconds": export16_seconds,
            "build_report": build_report,
            "cold_audit_report": audit_report,
            "export_report": export_report,
            "check": check,
        }
    payload = {
        "schema": "solidcraft_sketch_analysis_oracle_v1",
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "cli": str(cli.relative_to(root) if cli.is_relative_to(root) else cli),
            "commands": [
                "pa5-cli build RECIPE MODEL",
                "pa5-cli audit MODEL",
                "pa5-cli export-sketch-analysis MODEL 8 PACKAGE",
                "pa5-cli export-sketch-analysis MODEL 16 PACKAGE",
            ],
        },
        "sampling_scope": "exact Fraction evaluation on a common 3x3x3 Q2 parameter grid; sampled quality is diagnostic only",
        "cases": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "verified", "output": str(args.output), "cases": list(results)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
