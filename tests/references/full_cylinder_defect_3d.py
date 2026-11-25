#!/usr/bin/env python3
"""Independent SciPy reference for the archived h4 full-K3 control.

This oracle expands the exported Rust block-CSR K3 into SciPy, factors the
free 3D operator, and evaluates each available f64 field with ordinary,
compensated, and extended-precision row sums. It deliberately does not infer
that the h4 result qualifies h5/h6; those larger fields are checked against
their Rust-reassembled K3 operators and remain separate observations.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import resource
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.sparse import bsr_matrix
from scipy.sparse.linalg import splu

LIMIT_BYTES = 6 * 1024**3
resource.setrlimit(resource.RLIMIT_AS, (LIMIT_BYTES, LIMIT_BYTES))
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def relative_norm(vector: np.ndarray, rhs: np.ndarray) -> float:
    numerator = np.sqrt(np.sum(vector.astype(np.longdouble) ** 2, dtype=np.longdouble))
    denominator = max(float(np.sqrt(np.sum(rhs.astype(np.longdouble) ** 2, dtype=np.longdouble))), 1.0)
    return float(numerator / denominator)


def residual_metrics(matrix, field: np.ndarray, rhs: np.ndarray, free: np.ndarray) -> dict:
    plain = matrix @ field - rhs
    fsum_rows = np.empty(rhs.size, dtype=np.float64)
    longdouble_rows = np.empty(rhs.size, dtype=np.longdouble)
    data = matrix.data
    indices = matrix.indices
    indptr = matrix.indptr
    xlong = field.astype(np.longdouble)
    datalong = data.astype(np.longdouble)
    rhslong = rhs.astype(np.longdouble)
    for row in range(rhs.size):
        lo, hi = indptr[row : row + 2]
        ids = indices[lo:hi]
        fsum_rows[row] = math.fsum(float(a) * float(b) for a, b in zip(data[lo:hi], field[ids])) - rhs[row]
        longdouble_rows[row] = np.sum(datalong[lo:hi] * xlong[ids], dtype=np.longdouble) - rhslong[row]
    return {
        "relative_plain_f64": relative_norm(plain[free], rhs[free]),
        "relative_fsum_products_f64": relative_norm(fsum_rows[free], rhs[free]),
        "relative_longdouble_products_and_sum": relative_norm(longdouble_rows[free], rhs[free]),
        "max_abs_plain_vs_longdouble_n": float(np.max(np.abs(plain[free].astype(np.longdouble) - longdouble_rows[free]))),
    }


def audit_h4(csr_path: Path, direct_path: Path, corrected_path: Path) -> dict:
    exported = read_json(csr_path)
    if exported.get("schema") != "full_cylinder_sparse_system_v1" or exported["level"] != 4 or exported["order"] != 8:
        raise ValueError("unexpected h4 K3 CSR schema or discretization")
    block = exported["matrix_block_csr"]
    if block["block_size"] != 3:
        raise ValueError("expected interleaved xyz 3x3 block CSR")
    data = np.asarray(block["values_3x3"], dtype=np.float64)
    matrix = bsr_matrix(
        (data, np.asarray(block["column_indices"], dtype=np.int64), np.asarray(block["row_ptr"], dtype=np.int64)),
        shape=(3 * block["block_rows"], 3 * block["block_rows"]),
        blocksize=(3, 3),
    ).tocsr()
    rhs = np.asarray(exported["rhs_xyz"], dtype=np.float64)
    fixed = np.asarray(exported["fixed_vector_dofs"], dtype=np.int64)
    all_dofs = np.arange(rhs.size, dtype=np.int64)
    free = np.setdiff1d(all_dofs, fixed, assume_unique=True)
    if matrix.shape[0] != rhs.size or matrix.nnz == 0 or free.size == 0:
        raise ValueError("inconsistent h4 K3 export")

    direct = read_json(direct_path)
    corrected = read_json(corrected_path)
    for solution in (direct, corrected):
        if solution.get("schema") != "full_cylinder_solution_v1" or solution["level"] != 4 or solution["order"] != 8:
            raise ValueError("unexpected h4 solution schema")
    x_direct = np.asarray(direct["displacement_xyz"], dtype=np.float64)
    x_corrected = np.asarray(corrected["displacement_xyz"], dtype=np.float64)
    factor_started = time.monotonic()
    reduced = matrix[free][:, free].tocsc()
    factor = splu(reduced)
    factor_s = time.monotonic() - factor_started
    oracle = np.zeros_like(rhs)
    oracle[free] = factor.solve(rhs[free])
    return {
        "schema": "full_cylinder_defect_3d_reference_v1",
        "case": {"level": 4, "order": 8, "nu": exported["nu"], "formulation": exported["formulation"]},
        "method": "SciPy expands archived 3x3 block CSR, sparse-LU solves the free K3 equations, and evaluates f-K3*u with naive, math.fsum product, and longdouble row sums.",
        "factorization_s": factor_s,
        "dimension": int(rhs.size),
        "nnz_scalar": int(matrix.nnz),
        "free_dofs": int(free.size),
        "fields": {
            "p4_direct": residual_metrics(matrix, x_direct, rhs, free),
            "p5b_k3_corrected": residual_metrics(matrix, x_corrected, rhs, free),
            "scipy_sparse_lu": residual_metrics(matrix, oracle, rhs, free),
        },
        "field_agreement": {
            "relative_l2_p5b_vs_scipy_lu": float(
                np.linalg.norm((x_corrected[free] - oracle[free]).astype(np.longdouble))
                / max(float(np.linalg.norm(oracle[free].astype(np.longdouble))), np.finfo(float).tiny)
            ),
            "max_abs_p5b_minus_scipy_lu_mm": float(np.max(np.abs(x_corrected[free] - oracle[free]))),
        },
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "limits": {"memory_address_space_bytes": LIMIT_BYTES, "residual_acceptance": 1e-10},
        "scope_note": "This independent sparse-LU reference covers the h4 control only. Larger h5/h6 K3 operators are reassembled and residual-checked by the Rust route; no h4 result is extrapolated to those grids.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csr", type=Path)
    parser.add_argument("p4_solution", type=Path)
    parser.add_argument("p5b_solution", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    report = audit_h4(args.csr, args.p4_solution, args.p5b_solution)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False, sort_keys=True) + "\n")
    print(json.dumps({"case": report["case"], "fields": report["fields"], "factorization_s": report["factorization_s"]}, allow_nan=False))


if __name__ == "__main__":
    main()
