#!/usr/bin/env python3
"""Independent LAPACK spectrum audit for an exported bracket stiffness.

Usage: python3 tests/references/rigid_spectrum.py build/order4

Reads Matrix Market stiffness plus zero-based scalar `y50` constraints, computes
the unconstrained and reduced symmetric spectra with SciPy/LAPACK, and writes a
deterministic JSON report under build/ unless an output path is supplied.
"""

import argparse
import json
import os
from pathlib import Path
import time

# Fix common threaded BLAS implementations before NumPy/SciPy are imported.
for _name in (
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import numpy as np
from scipy.io import mmread
from scipy.linalg import eigvalsh

try:
    from threadpoolctl import threadpool_limits
except ImportError:  # Environment variables above still constrain standard BLAS builds.
    threadpool_limits = None


RELATIVE_ZERO_TOLERANCE = 1e-10


def _fixed_dofs(path, scalar_count):
    values = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        dof = int(line)
        if not 0 <= dof < scalar_count:
            raise ValueError(f"y50 scalar dof outside matrix: {dof}")
        values.append(dof)
    if len(values) != len(set(values)):
        raise ValueError("duplicate y50 scalar dof")
    return np.asarray(sorted(values), dtype=np.int64)


def _spectrum_report(matrix, threads_context):
    started = time.perf_counter()
    if threads_context is None:
        eigenvalues = eigvalsh(matrix.copy(), check_finite=False, driver="evd")
    else:
        with threads_context:
            eigenvalues = eigvalsh(matrix.copy(), check_finite=False, driver="evd")
    elapsed = time.perf_counter() - started
    largest = float(eigenvalues[-1])
    if not np.isfinite(largest) or largest <= 0.0:
        raise ValueError("spectrum has no positive finite largest eigenvalue")
    threshold = RELATIVE_ZERO_TOLERANCE * largest
    near_zero = np.abs(eigenvalues) <= threshold
    true_negative = eigenvalues < -threshold
    tiny_negative = (eigenvalues < 0.0) & near_zero
    return {
        "dimension": int(eigenvalues.size),
        "lambda_max": largest,
        "lambda_min": float(eigenvalues[0]),
        "relative_zero_tolerance": RELATIVE_ZERO_TOLERANCE,
        "absolute_zero_threshold": threshold,
        "near_zero_count": int(np.count_nonzero(near_zero)),
        "tiny_negative_within_zero_tolerance_count": int(np.count_nonzero(tiny_negative)),
        "true_negative_count": int(np.count_nonzero(true_negative)),
        "most_negative_eigenvalue": float(eigenvalues[0]),
        "lambda_6_one_based": float(eigenvalues[5]) if eigenvalues.size >= 6 else None,
        "lambda_7_one_based": float(eigenvalues[6]) if eigenvalues.size >= 7 else None,
        "eigensolve_seconds": elapsed,
    }


def analyze(directory, output_path=None):
    directory = Path(directory)
    matrix_path = directory / "stiffness.mtx"
    fixed_path = directory / "y50_dofs.txt"
    sparse = mmread(matrix_path)
    matrix = sparse.toarray() if hasattr(sparse, "toarray") else np.asarray(sparse)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("stiffness matrix must be square")
    if not np.isfinite(matrix).all():
        raise ValueError("stiffness matrix contains nonfinite values")
    scale = max(float(np.max(np.abs(matrix))), np.finfo(float).tiny)
    symmetry_error = float(np.max(np.abs(matrix - matrix.T)) / scale)
    if symmetry_error > 1e-12:
        raise ValueError(f"exported stiffness is not symmetric: relative error {symmetry_error:.3e}")
    # Remove only serialization roundoff before invoking the symmetric solver.
    matrix = 0.5 * (matrix + matrix.T)
    scalar_count = matrix.shape[0] // 3
    if 3 * scalar_count != matrix.shape[0]:
        raise ValueError("matrix dimension is not three times the scalar DOF count")
    fixed_scalar = _fixed_dofs(fixed_path, scalar_count)
    fixed_vector = np.concatenate(
        [3 * fixed_scalar + component for component in range(3)]
    )
    fixed_vector.sort()
    free_mask = np.ones(matrix.shape[0], dtype=bool)
    free_mask[fixed_vector] = False
    free_vector = np.flatnonzero(free_mask)
    if free_vector.size == 0:
        raise ValueError("y50 constraint removes every vector DOF")

    thread_context = threadpool_limits(limits=1) if threadpool_limits else None
    full = _spectrum_report(matrix, thread_context)
    reduced = matrix[np.ix_(free_vector, free_vector)]
    restrained = _spectrum_report(reduced, thread_context)
    report = {
        "schema": "solidcraft_rigid_spectrum_v1",
        "input_directory": str(directory),
        "matrix_market": str(matrix_path),
        "fixed_dof_file": str(fixed_path),
        "material": {"young_mpa": 210000.0, "poisson": 0.3},
        "quadrature_order": int(directory.name.removeprefix("order"))
        if directory.name.startswith("order")
        and directory.name.removeprefix("order").isdigit()
        else None,
        "scalar_dof_count": scalar_count,
        "vector_dof_count": int(matrix.shape[0]),
        "encastres_scalar_dof_count": int(fixed_scalar.size),
        "encastres_vector_dof_count": int(fixed_vector.size),
        "free_vector_dof_count": int(free_vector.size),
        "matrix_nnz": int(sparse.nnz) if hasattr(sparse, "nnz") else int(np.count_nonzero(matrix)),
        "symmetry_relative_error": symmetry_error,
        "blas_threads": 1,
        "free_body_spectrum": full,
        "y50_restrained_spectrum": restrained,
    }
    if output_path is None:
        repo = Path(__file__).resolve().parents[2]
        output_path = repo / "build" / f"rigid_spectrum_{directory.name}.json"
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report, output_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="directory containing stiffness.mtx and y50_dofs.txt")
    parser.add_argument("--output", type=Path, help="JSON report path (default: build/rigid_spectrum_<dir>.json)")
    args = parser.parse_args()
    report, output = analyze(args.directory, args.output)
    print(json.dumps({"report": str(output), **report}, sort_keys=True))


if __name__ == "__main__":
    main()
