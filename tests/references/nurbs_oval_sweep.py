"""Reproduce the 16-case synchronized NURBS meridian sweep."""

import argparse
import itertools
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import BSpline


def meridian(lateral, heights, weights, knots):
    points = np.array([lateral, heights, np.ones(len(weights), dtype=int)]).T
    controls = points * np.array(weights)[:, None]
    recipe = {"degree": 3, "knots": knots, "controls": controls.tolist()}
    return recipe, BSpline(knots, controls, 3)


def volume(curves, knots):
    # Outside the common 10 mm height the full 60 by 30 section is material.
    def area_times_dy(u):
        a, y, w = curves[0](u)
        b = curves[1](u)[0]
        _, dy, dw = curves[0](u, nu=1)
        return (60 - 2 * a / w) * (30 - 2 * b / w) * (dy * w - y * dw) / w**2

    area_integral, _ = quad(area_times_dy, 0, 4, points=sorted(set(knots[4:-4])),
                            epsabs=1e-9, epsrel=1e-12)
    return 54000 + 2 * area_integral


def sweep(cli, work):
    fixture = json.loads(Path(__file__).with_name("nurbs_cross_oval.recipe.json").read_text())
    rows = []
    for width, weights, knot, spans in itertools.product(
        (5, 8), ([1, 2, 1, 2, 1], [2, 1, 3, 1, 2]), (1, 2), (2, 3)
    ):
        heights = [0, 1, 3, 5, 5]
        knots = [0] * 4 + ([knot] if spans == 2 else [knot, 3]) + [4] * 4
        if spans == 3:
            heights.insert(2, 2)
            weights = weights[:3] + [2] + weights[3:]
        curves = []
        for axis, lateral in zip("xz", ([width, width, 3, 1, 0], [9, 9, 5, 3, 0])):
            if spans == 3:
                lateral.insert(3, lateral[2] - 1)
            entry, curve = meridian(lateral, heights, weights, knots)
            fixture["meridians"][axis] = entry
            curves.append(curve)
        recipe, model, space = (work / name for name in ("recipe", "model", "space"))
        recipe.write_text(json.dumps(fixture))
        for command in (["build", recipe, model], ["audit", model],
                        ["export-analysis", model, "16", space]):
            subprocess.run([cli, *command], check=True, stdout=subprocess.DEVNULL)
        native = json.loads(space.read_text())["integrated_volume_mm3"]
        reference = volume(curves, knots)
        count = len(json.loads(model.read_text())["blocks"])
        error = abs(native - reference) / reference
        if count != 8 * (2 * spans + 1) or error > 1e-9:
            raise ValueError(f"Unexpected block count or volume: {count}, {error}")
        rows.append(dict(width=width, weights=weights, knots=knots, spans=spans,
                         blocks=count, independent_volume_mm3=reference,
                         q16_volume_mm3=native, relative_error=error))
    return {"cases": rows, "max_relative_error": max(r["relative_error"] for r in rows)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cli", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as directory:
        result = sweep(args.cli.resolve(), Path(directory))
    args.output.write_text(json.dumps(result, indent=2) + "\n")
