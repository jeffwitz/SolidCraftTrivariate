"""Independent dyadic-knot insertion, THB truncation, trace and Ritz controls.

Uses exact Fraction knot insertion for the basis, then a fine-grid tensor
Galerkin pullback with SciPy for the independent harmonic solution.
"""

import argparse
import bisect
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import linear_sum_assignment
from scipy.sparse import coo_matrix, diags, kron
from scipy.sparse.linalg import spsolve

ROOT = Path(__file__).resolve().parents[1]


def knots(level):
    n = 8 << level
    return [F(0)] * 3 + [F(i, n) for i in range(n + 1)] + [F(1)] * 3


@lru_cache(None)
def dyadic(level):
    """Old coefficient -> fine coefficients, using single-knot Boehm insertions."""
    k = knots(level)
    count = len(k) - 4
    rows = [{i: F(1)} for i in range(count)]
    for t in sorted(set(knots(level + 1)) - set(k)):
        s = bisect.bisect_right(k, t) - 1
        nxt = rows[: s - 2]
        for i in range(s - 2, s + 1):
            a = (t - k[i]) / (k[i + 3] - k[i])
            r = defaultdict(F)
            for j, v in rows[i - 1].items():
                r[j] += (1 - a) * v
            for j, v in rows[i].items():
                r[j] += a * v
            nxt.append(dict(r))
        rows = nxt + rows[s:]
        k.insert(s + 1, t)
    assert k == knots(level + 1) and all(sum(r.values()) == 1 for r in rows)
    columns = [[] for _ in range(count)]
    for i, r in enumerate(rows):
        for j, v in r.items():
            if v:
                columns[j].append((i, v))
    return columns


@lru_cache(None)
def patch(level, uniform, width=4):
    """Fine uniform tensor rows of a local chart in its active THB columns."""
    coeff = {}
    for depth in range(level + 1):
        k = knots(depth)
        all_ids = [(i, j) for i in range(len(k) - 4) for j in range(len(k) - 4)]
        if uniform:
            if depth == level:
                coeff = {ij: {(depth, *ij): F(1)} for ij in all_ids}
            continue
        end = F(width, 8 << (depth - 1)) if depth else F(1)
        next_end = F(width, 8 << depth)
        if depth:
            child = defaultdict(lambda: defaultdict(F))
            mapping = dyadic(depth - 1)
            for (i, j), r in coeff.items():
                for ii, a in mapping[i]:
                    for jj, b in mapping[j]:
                        if max(k[ii + 4], k[jj + 4]) > end:
                            for key, v in r.items():
                                child[ii, jj][key] += a * b * v
            coeff = dict(child)
        for i, j in all_ids:
            reach = max(k[i + 4], k[j + 4])
            if reach <= end and (depth == level or reach > next_end):
                coeff[i, j] = {(depth, i, j): F(1)}
    assert all(sum(r.values()) == 1 for r in coeff.values())
    keys = sorted({key for r in coeff.values() for key in r})
    dim = (8 << level) + 3
    assert len(coeff) == dim**2
    # The two chart axes use identical knots, so complete seam polynomials
    # have identical ordered coefficient rows after reflection of chart labels.
    for j in range(dim):
        assert coeff[0, j] == {
            (depth, b, a): c for (depth, a, b), c in coeff[j, 0].items()
        }
    return coeff, keys


def chart(p, u, v):
    return [(-u, -v), (v, -u), (-v, u)][p]


def target(p, u, v):
    x, y = chart(p, u, v)
    theta = np.mod(np.arctan2(-y, x), 2 * np.pi)
    r = np.hypot(x, y)
    a = 2 / 3
    s, c = np.sin(a * theta), np.cos(a * theta)
    scale = a * r ** (a - 2)
    return r**a * s, scale * (s * x + c * y), scale * (s * y - c * x)


def gluing_id(p, key):
    depth, i, j = key
    if p == 1 and j == 0:
        return 0, depth, 0, i
    if p == 2 and i == 0:
        return 0, depth, j, 0
    return p, depth, i, j


def matrices(level, uniform, width=4):
    coeff, keys = patch(level, uniform, width)
    labels = sorted({gluing_id(p, k) for p in range(3) for k in keys})
    ids = {label: i for i, label in enumerate(labels)}
    d = (8 << level) + 3
    matrices = []
    for p in range(3):
        r, c, v = [], [], []
        for (i, j), terms in coeff.items():
            for k, a in terms.items():
                r.append(d * i + j)
                c.append(ids[gluing_id(p, k)])
                v.append(float(a))
        matrices.append(coo_matrix((v, (r, c)), shape=(d * d, len(ids))).tocsr())
    return matrices


def axis(level, order):
    n = 8 << level
    k = np.array(knots(level), dtype=float)
    q, w = np.polynomial.legendre.leggauss(order)
    q = ((np.arange(n)[:, None] + (q + 1) / 2) / n).ravel()
    w = np.tile(w / (2 * n), n)
    B = BSpline.design_matrix(q, k, 3).tocsr()
    # Derivative from independently constructed quadratic B-splines.
    D = BSpline.design_matrix(q, k[1:-1], 2).tocsr()
    delta = 3 / (k[4:-1] - k[1:-4])
    D = D @ diags((-delta, delta), (0, 1), shape=(len(delta), len(delta) + 1))
    return q, w, B, D


def replay(row, level, uniform, width=4):
    A = matrices(level, uniform, width)
    n = 8 << level
    d = n + 3
    q, w, B, D = axis(level, 16)
    M = B.T @ diags(w) @ B
    S = D.T @ diags(w) @ D
    Kone = kron(M, S) + kron(S, M)
    K = sum(a.T @ Kone @ a for a in A).tocsr()
    boundary = coo_matrix(K.shape).tocsr()
    load = np.zeros(K.shape[0])
    traces = []
    for p, a in enumerate(A):
        faces = [(0, 1), (1, 1)] + ([(0, 0)] if p == 1 else [(1, 0)] if p == 2 else [])
        for axis_id, hi in faces:
            idx = (
                hi * (d - 1) * d + np.arange(d)
                if axis_id == 0
                else hi * (d - 1) + np.arange(d) * d
            )
            trace = B @ a[idx, :]
            boundary += trace.T @ diags(w) @ trace
            t = np.full_like(q, hi)
            values = target(p, t, q)[0] if axis_id == 0 else target(p, q, t)[0]
            load += trace.T @ (w * values)
            traces.append((trace, values))
    fixed = np.flatnonzero(boundary.diagonal() > 0)
    free = np.setdiff1d(np.arange(K.shape[0]), fixed)
    solution = np.zeros(K.shape[0])
    solution[fixed] = spsolve(boundary[fixed][:, fixed].tocsc(), load[fixed])
    solution[free] = spsolve(K[free][:, free].tocsc(), -(K @ solution)[free])
    assert (len(solution), len(free)) == (
        row["solve"]["active_scalar_dofs"],
        row["solve"]["free_scalar_dofs"],
    )
    edge_error = np.sqrt(
        sum(w @ (values - trace @ solution) ** 2 for trace, values in traces)
        / sum(w @ values**2 for _, values in traces)
    )
    if not np.isclose(
        edge_error, row["solve"]["relative_boundary_L2"], rtol=1e-4, atol=1e-12
    ):
        raise ValueError("projected boundary error differs from the native measurement")
    # Corner square integrated separately on two triangles with a radial
    # cubic substitution; all other squares use ordinary tensor quadrature.
    accum = np.zeros(4)
    k = np.array(knots(level), dtype=float)
    spline = BSpline(k, np.eye(d), 3)

    def integrate(values, gradx, grady, truth, weights):
        val, gx, gy = truth
        return np.array(
            [
                np.sum(weights * ((gradx - gx) ** 2 + (grady - gy) ** 2)),
                np.sum(weights * (gx * gx + gy * gy)),
                np.sum(weights * (values - val) ** 2),
                np.sum(weights * val**2),
            ]
        )

    weights = w[:, None] * w[None, :]
    weights[:16, :16] = 0
    points, ww = np.polynomial.legendre.leggauss(24)
    r = (points + 1) / 2
    corner_weights = (ww[:, None] * ww[None, :] / 4) * 3 * r[:, None] ** 5 / n**2
    for p, a in enumerate(A):
        ctrl = (a @ solution).reshape(d, d)
        ub = (B @ ctrl) @ B.T
        gx, gy = chart(p, (D @ ctrl) @ B.T, (B @ ctrl) @ D.T)
        accum += integrate(ub, gx, gy, target(p, q[:, None], q[None, :]), weights)
        for swap in [False, True]:
            u, v = np.broadcast_arrays(
                r[:, None] ** 3 / n, r[:, None] ** 3 * r[None, :] / n
            )
            if swap:
                u, v = v, u
            bu, bv = spline(u.ravel()), spline(v.ravel())
            du, dv = spline(u.ravel(), nu=1), spline(v.ravel(), nu=1)

            def sample(x, y):
                return np.einsum("ij,ij->i", x @ ctrl, y).reshape(u.shape)

            gx, gy = chart(p, sample(du, bv), sample(bu, dv))
            accum += integrate(sample(bu, bv), gx, gy, target(p, u, v), corner_weights)
    norms = np.sqrt(accum[[0, 2]] / accum[[1, 3]])
    official = [row["solve"]["relative_H1_seminorm"], row["solve"]["relative_L2"]]
    if not np.allclose(norms, official, rtol=3e-6, atol=2e-10):
        raise ValueError(f"corner Ritz norms {norms} vs native {official}")
    for probe in row["basis_probes"]:
        p = probe["patch"]
        u, v = probe["uv"]

        def fn(bu, bv):
            return np.asarray(np.kron(bu, bv) @ A[p]).ravel()

        val = fn(spline(u), spline(v))
        gx, gy = chart(
            p, fn(spline(u, nu=1), spline(v)), fn(spline(u), spline(v, nu=1))
        )
        found = np.stack([val, gx, gy], axis=1)[val != 0]
        reference = np.asarray(probe["values_and_xy_gradients"])
        if reference.shape != found.shape:
            raise ValueError("native and oracle active counts disagree")
        distance = np.abs(found[:, None, :] - reference[None, :, :]).max(axis=2)
        pairs = linear_sum_assignment(distance)
        if distance[pairs].max() > 1e-9:
            raise ValueError(f"basis or physical gradients differ at {probe['uv']}")
    return dict(
        relative_H1_seminorm=float(norms[0]),
        relative_L2=float(norms[1]),
        basis_probes=len(row["basis_probes"]),
        free_scalar_dofs=len(free),
        relative_boundary_L2=float(edge_error),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=ROOT / "build/reentrant-ritz")
    args = parser.parse_args()
    results = []
    for depth in range(4):
        for mode in ["thb", "uniform"]:
            row = json.loads((args.work / f"{mode}-{depth}.json").read_text())
            results.append(
                dict(curve=mode, level=depth, **replay(row, depth, mode == "uniform"))
            )
            print(mode, depth, results[-1], flush=True)
    wide = []
    for level in range(1, 4):
        row = json.loads((args.work / f"wide-{level}.json").read_text())
        evidence = replay({"solve": row, "basis_probes": []}, level, False, width=6)
        wide.append(dict(curve="wide", level=level, **evidence))
    (args.work / "independent-wide.json").write_text(json.dumps(wide, indent=2) + "\n")
    published = json.loads(
        (ROOT / "docs/results/continuity-followup/reentrant-ritz.json").read_text()
    )
    if len(published["states"]) != len(results):
        raise ValueError("missing singular convergence states")
    records = results + wide
    references = published["states"] + published["supplementary_width_six"]
    if len(records) != len(references):
        raise ValueError("missing supplementary corner results")
    for current, frozen in zip(records, references):
        previous = frozen["independent"]
        if (current["curve"], current["level"], current["free_scalar_dofs"]) != (
            previous["curve"],
            previous["level"],
            previous["free_scalar_dofs"],
        ):
            raise ValueError("changed corner space dimension or level")
        for k in ("relative_H1_seminorm", "relative_L2", "relative_boundary_L2"):
            if not np.isclose(current[k], previous[k], rtol=5e-7, atol=2e-12):
                raise ValueError(f"singular result changed: {current} vs {previous}")
    if (
        json.loads((args.work / "efficiency.json").read_text())
        != published["efficiency_gate"]
    ):
        raise ValueError("the prospective error-per-DOF verdict changed")
    (args.work / "independent.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
