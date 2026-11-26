"""Verify stock-rim homogeneous charts using integer Bernstein coefficients and quadrature."""

import itertools
import json
from fractions import Fraction as F
import math
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss


def normalized_coefficients(net):
    """Coefficients in the unnormalized tensor basis t^i (1-t)^(n-i)."""
    result = net.copy()
    for axis, count in enumerate(net.shape):
        shape = [1, 1, 1]
        shape[axis] = count
        result *= np.array(
            [math.comb(count - 1, i) for i in range(count)], dtype=object
        ).reshape(shape)
    return result


def product(a, b):
    result = np.zeros(np.add(a.shape, b.shape) - 1, dtype=object)
    for i in np.ndindex(a.shape):
        if a[i]:
            result[tuple(slice(k, k + n) for k, n in zip(i, b.shape))] += a[i] * b
    return result


def minors(a, b):
    return {
        (i, j): product(a[i], b[j]) - product(a[j], b[i])
        for i, j in itertools.combinations(range(4), 2)
    }


def certify(net):
    """Expand det[Hu,Hv,Hw,H]. Equal control scaling cancels against W^4."""
    rows = []
    for axis in range(3):
        d = (net.shape[axis] - 1) * np.diff(net, axis=axis)
        rows.append([normalized_coefficients(d[..., i]) for i in range(4)])
    h = [normalized_coefficients(net[..., i]) for i in range(4)]
    first, second = minors(*rows[:2]), minors(rows[2], h)
    total = 0
    for (i, j), minor in first.items():
        other = tuple(sorted(set(range(4)) - {i, j}))
        total += (-1) ** (i + j - 1) * product(minor, second[other])
    # On an open parameter stratum, every Bernstein basis function in its
    # restriction is positive. One positive coefficient suffices there.
    if total.min() < 0:
        raise ValueError("negative Bernstein coefficient; further subdivision required")
    zero_strata = []
    for face in itertools.product((-1, 0, 1), repeat=3):
        sl = tuple(slice(None) if k == 0 else (0 if k < 0 else -1) for k in face)
        if np.max(total[sl]) == 0:
            zero_strata.append(face)
    return zero_strata


def basis(degree, t):
    return np.stack(
        [
            math.comb(degree, i) * t**i * (1 - t) ** (degree - i)
            for i in range(degree + 1)
        ],
        axis=1,
    )


def integrate(net, order):
    t, wt = leggauss(order)
    ts = (t + 1) / 2
    deg = np.array(net.shape[:3]) - 1
    bases = [basis(d, ts) for d in deg]

    def evaluate(p, bs):
        return np.einsum("ia,jb,kc,abch->ijkh", *bs, p, optimize=True)

    p = np.array(net, dtype=float)
    h = evaluate(p, bases)
    derivatives = []
    for a, d in enumerate(deg):
        bs = bases.copy()
        bs[a] = basis(d - 1, ts)
        q = evaluate(d * np.diff(p, axis=a), bs)
        derivatives.append(
            (q[..., :3] - h[..., :3] * q[..., 3:] / h[..., 3:]) / h[..., 3:]
        )
    jac = np.linalg.det(np.stack(derivatives, axis=-1))
    return np.einsum("i,j,k,ijk", wt / 2, wt / 2, wt / 2, jac)


def verify_plane_section(net, feature, state):
    # Every u row must satisfy the same cross-section law. Five distinct
    # t values establish the circle identity after clearing its quartic denominator.
    normal = 0 if "/x" in feature else 1
    limit = (100, 50)[normal] if feature.endswith("max") else 0
    inward = -1 if limit else 1
    for row in net[:, 0]:
        bottom = row[0, 2] == 0
        cap = 0 if bottom else 12
        for t in (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)):
            h = [
                sum(
                    k * p[c]
                    for k, p in zip(((1 - t) ** 2, 2 * t * (1 - t), t * t), row)
                )
                for c in range(4)
            ]
            inset = inward * (h[normal] / h[3] - limit)
            depth = abs(h[2] / h[3] - cap)
            r = F(1, 4)
            assert 0 <= inset <= r and 0 <= depth <= r
            if state == "sharp":
                assert inset == 0
            elif state == "chamfer":
                assert inset + depth == r
            else:
                assert (inset - r) ** 2 + (depth - r) ** 2 == r * r


def main(parent_file, directory):
    parent = json.loads(Path(parent_file).read_text())
    folder = Path(directory)
    volumes = {}
    outer = {
        f["block"]: f["feature"]
        for f in parent["design_faces"]
        if f["feature"] in {"stock/xmin", "stock/xmax", "stock/ymin", "stock/ymax"}
    }
    for state in ("sharp", "chamfer", "fillet"):
        model = json.loads((folder / (state + ".json")).read_text())
        assert model["topology"] == parent["topology"]
        assert [v["id"] for v in model["volumes"]] == [
            b["id"] for b in parent["blocks"]
        ]
        positive = singular = 0
        measures = []
        for volume in model["volumes"]:
            for s in volume["spans"]:
                net = np.array(s["controls"], dtype=object).reshape(
                    *(d + 1 for d in s["degrees"]), 4
                )
                assert net[..., 3].min() > 0
                if volume["id"] in outer and net.shape[2] == 3:
                    verify_plane_section(net, outer[volume["id"]], state)
                zeros = certify(net)
                if zeros:
                    assert (
                        state == "fillet"
                        and volume["id"] in outer
                        and net.shape[2] == 3
                    )
                    cap = -1 if net[0, 0, 0, 2] == 0 else 1
                    assert set(zeros) == {(u, -1, cap) for u in (-1, 0, 1)}, (
                        s["id"],
                        zeros,
                    )
                    singular += 1
                else:
                    positive += 1
                measures.append([integrate(net, q) for q in (12, 16)])
        assert (singular > 0) == (state == "fillet")
        volumes[state] = np.sum(measures, axis=0)
        print(
            state,
            "closed positive spans",
            positive,
            "spans singular only at cap-wall rail",
            singular,
            "volumes q12/q16",
            volumes[state],
            flush=True,
        )
    # Each perimeter shaves off integrand 2(L+B)d - 4d². The
    # through-hole area cancels when taking differences from the sharp state.
    r = 0.25
    for state, a, b in [
        ("chamfer", 0.5, 1 / 3),
        ("fillet", 1 - math.pi / 4, 5 / 3 - math.pi / 2),
    ]:
        expected = 2 * (300 * r * r * a - 4 * r**3 * b)
        removed = volumes["sharp"] - volumes[state]
        assert np.max(abs(removed - expected)) < 1e-7, (state, removed, expected)
        print(state, "removed mm^3", expected)
    assert max(abs(v[0] - v[1]) for v in volumes.values()) < 2e-6


if __name__ == "__main__":
    main(*sys.argv[1:])
