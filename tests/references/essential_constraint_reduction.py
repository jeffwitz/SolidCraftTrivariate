#!/usr/bin/env python3
"""Independent SciPy oracle for exact elimination of dependent supports.

This checks the frozen v2 refusal/v3 acceptance regression. It rebuilds the
elastic stiffness and traction from exported geometry, checks the affine lift
against the physical support, and independently solves the reduced equations.
Small rational systems exercise signed coefficients, offsets and contradiction
handling without calling the Rust elimination implementation.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.sparse import block_diag, coo_matrix
from scipy.sparse.linalg import spsolve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mechanical_job as oracle


def exact_reduce(matrix, rhs):
    """Return RREF equations and pivot columns over exact rationals."""
    a = [[Fraction(x) for x in row] + [Fraction(y)] for row, y in zip(matrix, rhs)]
    n = len(matrix[0])
    pivot_cols = []
    row = 0
    for col in range(n):
        pivot = next((r for r in range(row, len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        scale = a[row][col]
        a[row] = [v / scale for v in a[row]]
        for r in range(len(a)):
            if r != row and a[r][col]:
                scale = a[r][col]
                a[r] = [x - scale * y for x, y in zip(a[r], a[row])]
        pivot_cols.append(col)
        row += 1
        if row == len(a):
            break
    if any(not any(line[:n]) and line[n] for line in a):
        raise ValueError("inconsistent affine constraints")
    return a, pivot_cols


def synthetic_checks():
    # Signed and dependent equations with a nonzero affine offset.
    matrix = [[1, -2, 0, 3], [-2, 4, 0, -6], [0, 1, -1, 2]]
    rhs = [7, -14, -5]
    reduced, pivots = exact_reduce(matrix, rhs)
    assert len(pivots) == 2
    x = [Fraction(1), Fraction(-3), Fraction(2), Fraction(0)]
    assert all(sum(Fraction(a) * b for a, b in zip(row, x)) == Fraction(y)
               for row, y in zip(matrix, rhs))
    assert all(sum(line[j] * x[j] for j in range(4)) == line[4]
               for line in reduced if any(line[:4]))
    try:
        exact_reduce([[1, -2], [2, -4]], [1, 3])
    except ValueError:
        contradiction_rejected = True
    else:
        contradiction_rejected = False
    assert contradiction_rejected
    return {"rank": len(pivots), "free_columns": 4 - len(pivots),
            "signed_coefficients": True, "nonzero_affine_rhs": True,
            "contradiction_rejected": contradiction_rejected}


def lift_matrix(result, ndof):
    maps, offsets = [], []
    for lift in result["component_lifts"]:
        rr, cc, vv = [], [], []
        g = np.zeros(ndof)
        for line in lift["lines"]:
            dof = int(line["dof"])
            for root, num, den in line["terms"]:
                rr.append(dof); cc.append(int(root)); vv.append(int(num) / int(den))
            g[dof] = int(line["offset"][0]) / int(line["offset"][1])
        maps.append(coo_matrix((vv, (rr, cc)),
                               shape=(ndof, len(lift["free_roots"]))).tocsr())
        offsets.append(g)
    p_component = block_diag(maps, format="csr")
    physical_rows = np.asarray([c * ndof + d for d in range(ndof) for c in range(3)])
    return p_component[physical_rows].tocsr(), np.concatenate(offsets)


def run(binary, root):
    fixture = json.loads(Path("tests/regressions/adaptive_dependent_support.json").read_text())
    recipe = root / "recipe.json"
    model_file = root / "model.json"
    job2_file = root / "job-v2.json"
    job3_file = root / "job-v3.json"
    result_file = root / "result-v3.json"
    recipe.write_text(json.dumps(fixture["recipe"], separators=(",", ":")))
    job2 = dict(fixture["job"]); job2["schema"] = "solidcraft_mechanical_job_v2"
    job3 = dict(fixture["job"]); job3["schema"] = "solidcraft_mechanical_job_v3"
    job2_file.write_text(json.dumps(job2, separators=(",", ":")))
    job3_file.write_text(json.dumps(job3, separators=(",", ":")))
    subprocess.run([str(binary), "build", str(recipe), str(model_file)],
                   check=True, capture_output=True, text=True, timeout=180)
    model = json.loads(model_file.read_text())
    old = subprocess.run([str(binary), "solve-job", str(model_file), str(job2_file), str(root / "old.json")],
                         capture_output=True, text=True, timeout=180)
    old_text = old.stdout + old.stderr
    expected_refusal = old.returncode != 0 and "prescribed dependent DOF requires elimination or basis adaptation" in old_text
    if not expected_refusal:
        raise AssertionError(f"v2 did not preserve the frozen refusal: rc={old.returncode} {old_text}")
    start = time.monotonic()
    subprocess.run([str(binary), "solve-job", str(model_file), str(job3_file), str(result_file)],
                   check=True, capture_output=True, text=True, timeout=1800)
    elapsed = time.monotonic() - start
    result = json.loads(result_file.read_text())
    job = result["job"]
    k, ndof = oracle.assemble(result["analysis_blocks"], job["material"]["young_mpa"],
                              job["material"]["poisson"], job["order"])
    # Integrate the exported selected exterior faces directly. The generic M5
    # selector assumes a single-segment source ID; this recipe's triangle IDs
    # themselves contain slashes, so selection is certified here from the
    # exact support plane on every serialized face control.
    by_id = {b["id"]: b for b in result["analysis_blocks"]}
    f = np.zeros(3 * ndof); applied = np.zeros(3); moment = np.zeros(3)
    area = 0.0; load_scale = 0.0; max_radius = 0.0
    for block_id, side in result["loads"][0]["faces"]:
        block = by_id[block_id]
        axis, upper = oracle.parse_side({"side": side})
        for local, h in enumerate(block["controls"]):
            ijk = (local // 9, (local // 3) % 3, local % 3)
            if ijk[axis] == 2 * upper and int(h[0]) != 16 * int(h[3]):
                raise AssertionError(f"selected load face is not exactly on x=16: {block_id}/{side}")
        ids = np.asarray(block["dofs"], dtype=int)
        for uvw, point, _grad, _normal, measure in oracle.face_rule(block, {"side": side}, job["order"]):
            traction = np.asarray(job["loads"][0]["vector_mpa"], dtype=float)
            area += measure; load_scale += float(np.linalg.norm(traction)) * measure
            applied += traction * measure; moment += np.cross(point, traction) * measure
            max_radius = max(max_radius, float(np.linalg.norm(point)))
            values = oracle.shape(block["controls"], uvw)[0]
            for i, node in enumerate(ids):
                f[3 * node:3 * node + 3] += values[i] * traction * measure
    p, g_component = lift_matrix(result, ndof)
    g = g_component.reshape(3, ndof).T.reshape(3 * ndof)
    reduced_k = p.T @ k @ p
    reduced_rhs = np.asarray(p.T @ (f - k @ g)).ravel()
    reduced = spsolve(reduced_k.tocsc(), reduced_rhs)
    independent_u = np.asarray(p @ reduced).ravel() + g
    reported_u = np.asarray(result["displacement_mm_xyz"], dtype=float)
    solution_error = float(np.linalg.norm(reported_u - independent_u) /
                           max(np.linalg.norm(independent_u), 1e-30))
    residual = k @ reported_u - f
    stationarity = float(np.linalg.norm(p.T @ residual) / max(load_scale, 1.0))
    support_max = 0.0
    u_xyz = reported_u.reshape(ndof, 3)
    for support in result["supports"]:
        for block_id, side in support["faces"]:
            block = by_id[block_id]
            axis, upper = oracle.parse_side({"side": side})
            for local, node in enumerate(block["dofs"]):
                ijk = (local // 9, (local // 3) % 3, local % 3)
                if ijk[axis] == 2 * upper:
                    for component in support["components"]:
                        support_max = max(support_max, abs(float(u_xyz[int(node), component])))
    reaction = residual.reshape(ndof, 3).sum(axis=0)
    resultant_error = float(np.linalg.norm(reaction + applied) / max(load_scale, 1.0))
    record = {
        "fixture": fixture.get("id", "adaptive_dependent_support"), "schema": job["schema"], "old_schema_refused": expected_refusal,
        "old_refusal": old_text.strip(), "blocks": len(result["analysis_blocks"]),
        "scalar_dofs": ndof, "free_roots_per_component": len(result["component_lifts"][0]["free_roots"]),
        "prescribed_support_faces": sum(len(s["faces"]) for s in result["supports"]),
        "support_max_displacement_mm": support_max,
        "scipy_reduced_solution_relative_error": solution_error,
        "free_space_virtual_work_relative": stationarity,
        "global_force_balance_relative": resultant_error,
        "relative_residual_rust": result["relative_residual"],
        "relative_force_balance_rust": result["relative_force_balance"],
        "volume_mm3": result["volume_mm3"], "load_area_mm2": area,
        "solve_wall_s": elapsed, "synthetic_exact_checks": synthetic_checks(),
    }
    if support_max > 1e-12 or solution_error > 2e-8 or stationarity > 2e-8 or resultant_error > 2e-8:
        raise AssertionError(f"independent qualification failed: {record}")
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=Path("target/release/solidcraft-trivariate"))
    parser.add_argument("--out", type=Path, default=Path("docs/results/essential-constraint-reduction-20261125.json"))
    parser.add_argument("--work", type=Path, default=Path("build/essential-oracle"))
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    record = run(args.binary.resolve(), args.work.resolve())
    report = {"protocol": "M8c independent exact-support reduction oracle v1",
              "environment": {"python": sys.version.split()[0], "numpy": np.__version__},
              "command": "python3 tests/references/essential_constraint_reduction.py",
              "result": record,
              "scope": ["frozen v2 refusal", "same geometry and supports in v3", "independent SciPy stiffness/load solve",
                        "exact affine lift support satisfaction", "free-space virtual work", "synthetic signed rational constraints"],
              "limits": ["synthetic exact elimination is an oracle check, not a proof of all integer-capacity cases",
                         "the FE solve qualification uses the complete frozen fixture plan and its 2968-cell final partition"]}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
