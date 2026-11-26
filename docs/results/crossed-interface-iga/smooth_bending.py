"""Compare C0 and C1 quadratic fields on the same conforming S2 integration cells."""

import argparse
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

import conforming_bending_study as conforming


def prolongation(elements, lookup, labels, free_map):
    """Map independent coefficients into the rational C0 field at fixed weights."""
    relations = {}
    checks = []
    for axis in (1, 2):
        for e, (chart, iv, iw, net) in enumerate(elements):
            neighbor = lookup.get((chart, iv + (axis == 1), iw + (axis == 2)))
            if neighbor is None:
                continue
            # Uniform spans and equal degrees give equal derivative factors.
            # Check continuity of the homogeneous geometry and its denominator.
            a = np.take(net, -2, axis=axis)
            m = np.take(net, -1, axis=axis)
            b = np.take(elements[neighbor][3], 1, axis=axis)
            np.testing.assert_allclose((a + b) / 2, m, rtol=3e-13, atol=3e-12)
            left = np.take(labels[e], -2, axis=axis)
            center = np.take(labels[e], -1, axis=axis)
            right = np.take(labels[neighbor], 1, axis=axis)
            for i in np.ndindex(center.shape):
                left_id, c, r = (int(free_map[v[i]]) for v in (left, center, right))
                alpha, beta = a[i][-1] / (2 * m[i][-1]), b[i][-1] / (2 * m[i][-1])
                if left_id == c == r:
                    continue
                if min(left_id, c, r) < 0:
                    raise ValueError("A knot removal crosses the Dirichlet boundary")
                if c in (left_id, r):
                    raise ValueError("Unexpected partially degenerate trace constraint")
                eq = ((left_id, float(alpha)), (r, float(beta)))
                relations.setdefault(c, eq)
                checks.append((c, eq))

    count = np.count_nonzero(free_map >= 0)
    independent = {
        i: j for j, i in enumerate(i for i in range(count) if i not in relations)
    }

    @lru_cache(None)
    def expand(i):
        if i in independent:
            return {independent[i]: 1.0}
        row = {}
        for child, weight in relations[i]:
            for j, value in expand(child).items():
                row[j] = row.get(j, 0.0) + weight * value
        return row

    rows, cols, entries = [], [], []
    for i in range(count):
        for j, value in expand(i).items():
            rows.append(i)
            cols.append(j)
            entries.append(value)
    transfer = sparse.csr_matrix(
        (entries, (rows, cols)), shape=(count, len(independent))
    )
    # All equations, including duplicates at intersecting planes and chart
    # interfaces, must be satisfied by every column of the prolongation.
    error = 0.0
    for i, eq in checks:
        diff = expand(i).copy()
        for child, weight in eq:
            for j, value in expand(child).items():
                diff[j] = diff.get(j, 0.0) - weight * value
        error = max(error, max(map(abs, diff.values()), default=0.0))
    if error > 1e-12 or not np.allclose(transfer.sum(axis=1), 1, rtol=0, atol=1e-12):
        raise ValueError("Inconsistent C1 prolongation")
    return sparse.kron(transfer, sparse.eye(3), format="csr"), error


def solve(binary, rise, nv, nw, q):
    base, locations = conforming.load_geometry(binary, rise)
    k, f, free_map, elements, lookup, labels, volume = conforming.solve_case(
        base, locations, (3, 2, 2), nv, nw, q, return_system=True
    )
    p, err = prolongation(elements, lookup, labels, free_map)
    rows = []
    for name, matrix, rhs in [("C0", k, f), ("C1", (p.T @ k @ p).tocsr(), p.T @ f)]:
        u = spsolve(matrix, rhs)
        residual = float(np.linalg.norm(matrix @ u - rhs) / np.linalg.norm(rhs))
        full = u if name == "C0" else p @ u
        work = float(f @ full)
        if residual > 2e-8 or not np.isclose(
            full @ (k @ full), work, rtol=2e-9, atol=0
        ):
            raise ValueError("Equilibrium or energy check failed")
        rows.append(
            dict(
                continuity=name,
                free_vector_dofs=len(u),
                mean_ux_mm=work / f[0::3].sum(),
                energy_N_mm=work / 2,
                relative_residual=residual,
            )
        )
    if rows[1]["mean_ux_mm"] > rows[0]["mean_ux_mm"] * (1 + 1e-9):
        raise ValueError("Nested Galerkin compliance ordering failed")
    return dict(
        rise_mm=rise,
        degree=[3, 2, 2],
        spans_vw=[nv, nw],
        gauss_order=q,
        geometry_charts=len(base),
        volume_mm3=volume,
        force_N=f[0::3].sum(),
        max_constraint_coefficient_residual=err,
        fields=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--order", type=int, default=8)
    args = parser.parse_args()
    runs = [solve(args.binary, 3, n, n, args.order) for n in (2, 3)]
    reference = runs[0]["fields"][0]
    smooth = runs[1]["fields"][1]
    if reference["free_vector_dofs"] != smooth["free_vector_dofs"]:
        raise ValueError("Expected equal dimensions for C0/2 spans and C1/3 spans")
    result = {
        "runs": runs,
        "equal_budget": {
            "free_vector_dofs": reference["free_vector_dofs"],
            "c0_two_spans_mean_ux_mm": reference["mean_ux_mm"],
            "c1_three_spans_mean_ux_mm": smooth["mean_ux_mm"],
            "c0_integration_cells": 22 * 2**2,
            "c1_integration_cells": 22 * 3**2,
        },
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["equal_budget"]))
