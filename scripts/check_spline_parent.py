#!/usr/bin/env python3
"""Independent exact oracle for a NURBS strip's recovered tensor parent."""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RECIPE = ROOT / "tests/references/nurbs_strip.recipe.json"


def basis_values(degree: int, knots: list[int], parameter: F) -> list[F]:
    count = len(knots) - degree - 1
    if parameter == knots[count]:
        return [F(int(i == count - 1)) for i in range(count)]
    values = [
        F(int(F(knots[i]) <= parameter < F(knots[i + 1])))
        for i in range(len(knots) - 1)
    ]
    for order in range(1, degree + 1):
        next_values = []
        for i in range(len(knots) - order - 1):
            left_den = knots[i + order] - knots[i]
            right_den = knots[i + order + 1] - knots[i + 1]
            left = (
                F(0) if left_den == 0 else (parameter - knots[i]) * values[i] / left_den
            )
            right = (
                F(0)
                if right_den == 0
                else (knots[i + order + 1] - parameter) * values[i + 1] / right_den
            )
            next_values.append(left + right)
        values = next_values
    return values


def source_curve(curve: dict, parameter: F) -> tuple[F, F]:
    degree = curve["degree"]
    basis = basis_values(degree, curve["knots"], parameter)
    homogeneous = [
        sum(basis[i] * curve["controls"][i][axis] for i in range(len(basis)))
        for axis in range(3)
    ]
    if homogeneous[2] <= 0:
        raise AssertionError("the source denominator must remain positive")
    return homogeneous[0] / homogeneous[2], homogeneous[1] / homogeneous[2]


def source_curve_jet(curve: dict, parameter: F) -> tuple[tuple[F, F], tuple[F, F]]:
    degree = curve["degree"]
    knots = curve["knots"]
    values = basis_values(degree, knots, parameter)
    lower = basis_values(degree - 1, knots, parameter)
    derivatives = []
    for i in range(len(values)):
        left_den = knots[i + degree] - knots[i]
        right_den = knots[i + degree + 1] - knots[i + 1]
        left = F(0) if left_den == 0 else F(degree, left_den) * lower[i]
        right = F(0) if right_den == 0 else F(degree, right_den) * lower[i + 1]
        derivatives.append(left - right)
    h = [
        sum(values[i] * curve["controls"][i][axis] for i in range(len(values)))
        for axis in range(3)
    ]
    dh = [
        sum(derivatives[i] * curve["controls"][i][axis] for i in range(len(values)))
        for axis in range(3)
    ]
    point = (h[0] / h[2], h[1] / h[2])
    derivative = tuple(
        (dh[axis] * h[2] - h[axis] * dh[2]) / h[2] ** 2 for axis in range(2)
    )
    return point, derivative


def bernstein(degree: int, i: int, t: F) -> F:
    from math import comb

    return F(comb(degree, i)) * t**i * (1 - t) ** (degree - i)


def block_point(block: dict, u: F, v: F, w: F) -> tuple[F, F, F]:
    degree = block["degrees"][0]
    homogeneous = [F(0) for _ in range(4)]
    for i in range(degree + 1):
        bu = bernstein(degree, i, u)
        for j in range(2):
            bv = (1 - v) if j == 0 else v
            for k in range(2):
                bw = (1 - w) if k == 0 else w
                control = block["controls"][(i * 2 + j) * 2 + k]
                factor = bu * bv * bw
                for axis in range(4):
                    homogeneous[axis] += factor * control[axis]
    if homogeneous[3] <= 0:
        raise AssertionError("the source block denominator must remain positive")
    return tuple(homogeneous[axis] / homogeneous[3] for axis in range(3))


def block_u_jet(block: dict, u: F, v: F, w: F) -> tuple[tuple[F, F, F], tuple[F, F, F]]:
    degree = block["degrees"][0]
    h = [F(0) for _ in range(4)]
    dh = [F(0) for _ in range(4)]
    for i in range(degree + 1):
        bu = bernstein(degree, i, u)
        dbu = F(0)
        if i:
            dbu += degree * bernstein(degree - 1, i - 1, u)
        if i < degree:
            dbu -= degree * bernstein(degree - 1, i, u)
        for j in range(2):
            bv = (1 - v) if j == 0 else v
            for k in range(2):
                bw = (1 - w) if k == 0 else w
                control = block["controls"][(i * 2 + j) * 2 + k]
                for axis in range(4):
                    h[axis] += bu * bv * bw * control[axis]
                    dh[axis] += dbu * bv * bw * control[axis]
    point = tuple(h[axis] / h[3] for axis in range(3))
    derivative = tuple(
        (dh[axis] * h[3] - h[axis] * dh[3]) / h[3] ** 2 for axis in range(3)
    )
    return point, derivative


def direct_parent(recipe: dict, t: F, v: F, w: F) -> tuple[F, F, F]:
    lower = source_curve(recipe["lower"], t)
    upper = source_curve(recipe["upper"], t)
    return (
        lower[0] + v * (upper[0] - lower[0]),
        lower[1] + v * (upper[1] - lower[1]),
        F(recipe["height_mm"]) * w,
    )


def check(recipe: dict, model: dict) -> dict:
    assert model["schema"] == "nurbs_sketch_strip_v1"
    assert model["lower"] == recipe["lower"] and model["upper"] == recipe["upper"]
    assert model["height_mm"] == recipe["height_mm"]
    intervals = model["span_intervals"]
    blocks = model["blocks"]
    degree = recipe["lower"]["degree"]
    source_knots = recipe["lower"]["knots"]
    source_count = len(recipe["lower"]["controls"])
    expected_intervals = [
        [source_knots[i], source_knots[i + 1]]
        for i in range(degree, source_count)
        if source_knots[i] < source_knots[i + 1]
    ]
    first_x = recipe["lower"]["controls"][0][0] * recipe["lower"]["controls"][-1][2]
    last_x = recipe["lower"]["controls"][-1][0] * recipe["lower"]["controls"][0][2]
    source_reversed = first_x > last_x
    if source_reversed:
        expected_intervals = [[b, a] for a, b in reversed(expected_intervals)]
    assert len(intervals) == len(blocks) == len(expected_intervals)
    assert [list(map(int, interval)) for interval in intervals] == expected_intervals
    assert [block["degrees"] for block in blocks] == [[degree, 1, 1]] * len(blocks)

    # Each coordinate of a rational cubic map is a quotient of cubics. Equality
    # at seven distinct exact parameters proves the cross-multiplied degree-six
    # polynomial identity on the entire source span.
    location_count = 2 * degree + 1
    exact_locations = [F(i, location_count - 1) for i in range(location_count)]
    comparisons = 0
    derivative_comparisons = 0
    for interval, block in zip(intervals, blocks):
        a, b = map(F, interval)
        for local in exact_locations:
            t = a + (b - a) * local
            for v in (F(0), F(1)):
                for w in (F(0), F(1)):
                    expected = direct_parent(recipe, t, v, w)
                    actual = block_point(block, local, v, w)
                    assert actual == expected, (block["id"], t, v, w, actual, expected)
                    comparisons += 1
        for local in (F(1, 5), F(2, 5), F(3, 5), F(4, 5)):
            t = a + (b - a) * local
            for v in (F(0), F(1)):
                direct, direct_du = (
                    source_curve_jet(recipe["lower"], t)
                    if v == 0
                    else source_curve_jet(recipe["upper"], t)
                )
                actual, actual_du = block_u_jet(block, local, v, F(0))
                assert actual[:2] == direct
                assert actual_du[:2] == tuple(
                    component * (b - a) for component in direct_du
                )
                derivative_comparisons += 1

    # The knot multiplicity theorem gives C^(p-m)=C2 for homogeneous curve
    # coordinates. Shared positive weights and x controls make it the parent
    # volume's exact u-direction continuity as well.
    interior = sorted(
        set(source_knots[degree + 1 : source_count])
        & set(source_knots[degree:source_count])
    )
    interior = [
        knot
        for knot in interior
        if source_knots[degree] < knot < source_knots[source_count]
    ]
    multiplicities = [source_knots.count(knot) for knot in interior]
    continuity_orders = [degree - multiplicity for multiplicity in multiplicities]
    assert all(order >= 0 for order in continuity_orders)
    assert recipe["upper"]["knots"] == recipe["lower"]["knots"]
    return {
        "status": "exact_parent_identity",
        "spans": len(blocks),
        "control_shape": [source_count, 2, 2],
        "parent_controls": source_count * 4,
        "source_span_occurrences": len(blocks) * (degree + 1) * 4,
        "c0_unique_controls": len(blocks) * (degree + 1) * 4 - (len(blocks) - 1) * 4,
        "exact_location_comparisons": comparisons,
        "exact_derivative_comparisons": derivative_comparisons,
        "interior_knots": interior,
        "interior_knot_multiplicities": multiplicities,
        "inherited_continuity_orders": continuity_orders,
        "oracle": "Fraction Cox-de Boor direct parent vs homogeneous Bernstein block maps",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cli", required=True, help="path to solidcraft-trivariate CLI"
    )
    parser.add_argument("--recipe", type=Path, default=DEFAULT_RECIPE)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--nonuniform-source", action="store_true")
    parser.add_argument("--reverse-source", action="store_true")
    args = parser.parse_args()
    recipe = json.loads(args.recipe.read_text())
    if args.nonuniform_source or args.reverse_source:
        recipe = json.loads(json.dumps(recipe))
        for curve_name in ("lower", "upper"):
            if args.nonuniform_source:
                recipe[curve_name]["knots"] = [-5, -5, -5, -5, 0, 4, 4, 4, 4]
            if args.reverse_source:
                knots = recipe[curve_name]["knots"]
                lo = knots[recipe[curve_name]["degree"]]
                hi = knots[len(recipe[curve_name]["controls"])]
                recipe[curve_name]["knots"] = [
                    lo + hi - knot for knot in reversed(knots)
                ]
                recipe[curve_name]["controls"].reverse()
    with tempfile.TemporaryDirectory(prefix="spline-parent-") as folder:
        temp = Path(folder)
        recipe_path = temp / "source.recipe.json"
        recipe_path.write_text(json.dumps(recipe))
        output = args.output or temp / "strip.model.json"
        output.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([args.cli, "build", recipe_path, output], check=True)
        model = json.loads(output.read_text())
        report = check(recipe, model)
        parent_file = temp / "parent.json"
        subprocess.run(
            [args.cli, "recover-spline-parent", output, parent_file], check=True
        )
        parent = json.loads(parent_file.read_text())
        actual = [
            [F(int(a["numerator"]), int(a["denominator"])) for a in h]
            for h in parent["homogeneous_controls"]
        ]
        expected = [
            [x, y, z * w, w]
            for curves in zip(
                recipe["lower"]["controls"], recipe["upper"]["controls"], strict=True
            )
            for x, y, w in curves
            for z in (0, recipe["height_mm"])
        ]
        assert actual == expected and parent["control_shape"] == report["control_shape"]
        assert parent["knots"] == [recipe["lower"]["knots"], [0, 0, 1, 1], [0, 0, 1, 1]]
        for span, interval, block in zip(
            parent["source_spans"],
            model["span_intervals"],
            model["blocks"],
            strict=True,
        ):
            restriction = span["restriction"]
            assert restriction["axes"] == [0, 1, 2] and restriction["intervals"] == [
                sorted(interval),
                [0, 1],
                [0, 1],
            ]
            assert (
                restriction["reversals"] == [interval[0] > interval[1], False, False]
                and span["source_block_id"] == block["id"]
            )
            # Gauge at the first corner is determined by the parent weight
            # evaluated at the (possibly reversed) source parameter.
            N = basis_values(
                recipe["lower"]["degree"], recipe["lower"]["knots"], F(interval[0])
            )
            W = sum(
                value * p[2]
                for value, p in zip(N, recipe["lower"]["controls"], strict=True)
            )
            assert F(*map(int, span["gauge"])) == F(block["controls"][0][3]) / W
        galerkin = json.loads(
            subprocess.check_output([args.cli, "check-parent-galerkin", output, "8"])
        )
        assert galerkin["exact_transfer_rank"] == len(expected)
        report["exported_homogeneous_parent_components"] = len(expected) * 4
        report["galerkin"] = galerkin
        print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
