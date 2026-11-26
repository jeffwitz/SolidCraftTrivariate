"""Compare persisted Boolean charts with de Boor evaluation of the source tools."""

import argparse
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import brentq


def ordinate_roots(curve, knots, y):
    roots = []
    for a, b in zip(knots, knots[1:]):

        def residual(t):
            _, wy, w = curve(t)
            return wy / w - y

        if residual(a) * residual(b) < 0:
            wx, _, w = curve(brentq(residual, a, b))
            roots.append(wx / w)
    return sorted(roots)


def section(blocks, y):
    """Intersect each serialized chart with the horizontal section at y."""

    def bernstein(t):
        return np.array([(1 - t) ** 3, 3 * t * (1 - t) ** 2, 3 * t * t * (1 - t), t**3])

    slabs = []
    for block in blocks:
        net = np.array(block["controls"], dtype=float).reshape(4, 2, 2, 4)
        row = net[:, 0, 0]

        def residual(t):
            return bernstein(t) @ (row[:, 1] - y * row[:, 3])

        if residual(0) * residual(1) >= 0:
            continue
        t = brentq(residual, 0, 1)
        homogeneous = np.einsum("i,ivwk->vwk", bernstein(t), net)
        physical = homogeneous[:, :, :3] / homogeneous[:, :, 3:]
        slabs.append(physical)
    return slabs


def check_points(blocks, curve, knots, rise, operation):
    checked = 0
    for y in [
        *np.linspace(14.037, 35.937, 29),
        *(level + d for level in (719 / 36, 1081 / 36) for d in (-1e-4, 1e-4)),
    ]:
        roots = ordinate_roots(curve, knots, y)
        xs = [
            *np.linspace(0.37, 99.37, 19),
            *(x + d for x in roots for d in [-1e-4, 1e-4]),
        ]
        slabs = section(blocks, y)
        for x in xs:
            inside_profile = sum(x > root for root in roots) % 2 == 1
            shift = rise * x / 100
            for z in [
                2.37,
                6.12,
                10.62,
                *(level + shift + d for level in (4, 8) for d in [-1e-4, 1e-4]),
            ]:
                a = inside_profile and 0 < z < 12
                b = 719 / 36 < y < 1081 / 36 and 4 < z - shift < 8
                expected = {
                    "intersection": a and b,
                    "union": a or b,
                    "a-minus-b": a and not b,
                    "b-minus-a": b and not a,
                }[operation]
                count = 0
                for slab in slabs:
                    right, left = slab[:, 0, 0]
                    if left < x < right:
                        v = (right - x) / (right - left)
                        low, high = (1 - v) * slab[0, :, 2] + v * slab[1, :, 2]
                        count += int(low < z < high)
                if count != int(expected):
                    raise AssertionError((operation, rise, x, y, z, expected, count))
                checked += 1
    return checked


def integrated_volume(blocks, order):
    # Cross-sectional area at a fixed u is width times mean height.
    # Integrating area * dy/du tests the persisted homogeneous controls.
    nodes, weights = np.polynomial.legendre.leggauss(order)
    total = 0.0
    for block in blocks:
        net = np.array(block["controls"], dtype=float).reshape(4, 2, 2, 4)
        for u, weight in zip((nodes + 1) / 2, weights / 2):
            v = 1 - u
            basis = np.array([v**3, 3 * v * v * u, 3 * v * u * u, u**3])
            derivative = np.array(
                [-3 * v * v, 3 * v * v - 6 * v * u, 6 * v * u - 3 * u * u, 3 * u * u]
            )
            h = np.einsum("i,ivwk->vwk", basis, net)
            d = derivative @ net[:, 0, 0, :]
            xyz = h[:, :, :3] / h[:, :, 3:]
            dy = (d[1] - xyz[0, 0, 1] * d[3]) / h[0, 0, 3]
            width = xyz[0, 0, 0] - xyz[1, 0, 0]
            height = np.mean(xyz[:, 1, 2] - xyz[:, 0, 2])
            total += weight * dy * width * height
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    args = parser.parse_args()
    path = Path(__file__).with_name("g2_cross_sketches.recipe.json")
    pair = json.loads(path.read_text())
    source = next(s for s in pair["sketches"] if s["direction"] == [0, 0, 1])["loops"][
        0
    ]["curve"]
    curve = BSpline(source["knots"], source["controls"], source["degree"])
    p = source["degree"]
    knots = np.unique(source["knots"][p:-p])
    with tempfile.TemporaryDirectory() as directory:
        recipe, output = (
            Path(directory) / name for name in ("input.json", "output.json")
        )
        for rise in [-3, 0, 3]:
            recipe.write_text(
                json.dumps(
                    dict(
                        family="inclined_crossed_tools_v1",
                        sketch_pair=pair,
                        transverse_rise_mm=rise,
                    )
                )
            )
            for op in ["intersection", "a-minus-b", "b-minus-a", "union"]:
                subprocess.run(
                    [args.binary.resolve(), "crossed-boolean", recipe, op, output],
                    check=True,
                    capture_output=True,
                )
                blocks = json.loads(output.read_text())["blocks"]
                probes = check_points(blocks, curve, knots, rise, op)
                volume = integrated_volume(blocks, 20)
                common, a, b = 712.3493494134794, 3489.731852993524, 4022.2222222222226
                target = {
                    "intersection": common,
                    "union": a + b - common,
                    "a-minus-b": a - common,
                    "b-minus-a": b - common,
                }[op]
                if (
                    abs(volume - target) > 1e-6
                    or abs(integrated_volume(blocks, 12) - volume) > 1e-6
                ):
                    raise AssertionError((op, rise, volume, target))
                print(rise, op, len(blocks), probes, f"{volume:.10f}")


if __name__ == "__main__":
    main()
