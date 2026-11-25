#!/usr/bin/env python3
"""Independent SciPy/LAPACK spectrum oracle for the three-cell T-junction."""

import argparse
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.io import mmread
from scipy.linalg import eigvalsh


def inspect(path, threshold_ratio, expected_modes):
    loaded = mmread(path)
    matrix = np.asarray(loaded.toarray() if hasattr(loaded, "toarray") else loaded, dtype=np.float64)
    eigenvalues = eigvalsh(matrix, check_finite=True, driver="evd")
    scale = float(np.max(np.abs(eigenvalues)))
    threshold = threshold_ratio * scale
    near_zero = eigenvalues[np.abs(eigenvalues) <= threshold]
    positive = eigenvalues[eigenvalues > threshold]
    negative = eigenvalues[eigenvalues < -threshold]
    return {
        "path": str(path),
        "dimension": int(matrix.shape[0]),
        "lambda_min": float(eigenvalues[0]),
        "lambda_max": float(eigenvalues[-1]),
        "threshold": threshold,
        "near_zero_count": int(near_zero.size),
        "near_zero_extrema": [float(near_zero[0]), float(near_zero[-1])] if near_zero.size else [],
        "near_zero_values": [float(value) for value in near_zero],
        "eigenvalue_after_expected_null_modes": float(eigenvalues[expected_modes]) if expected_modes < eigenvalues.size else None,
        "first_positive": float(positive[0]) if positive.size else None,
        "negative_below_threshold": int(negative.size),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", default="build/nested_spectrum")
    parser.add_argument("--output", default="build/nested_spectrum/spectrum.json")
    args = parser.parse_args()
    directory = Path(args.directory)
    metadata = json.loads((directory / "metadata.json").read_text())
    assert metadata["family"] == "three_cell_q2_t_junction_v1"
    threshold_ratio = metadata["threshold_ratio"]
    results = {
        "family": metadata["family"],
        "threshold_ratio": threshold_ratio,
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "lapack_driver": "evd",
        "orders": {},
    }
    for order in (3, 5):
        sub = directory / f"q{order}"
        free = inspect(sub / "free.mtx", threshold_ratio, 6)
        pinned = inspect(sub / "pinned.mtx", threshold_ratio, 0)
        assert free["near_zero_count"] == 6, (order, free)
        assert free["negative_below_threshold"] == 0, (order, free)
        assert pinned["near_zero_count"] == 0, (order, pinned)
        assert pinned["lambda_min"] > pinned["threshold"], (order, pinned)
        assert pinned["negative_below_threshold"] == 0, (order, pinned)
        results["orders"][str(order)] = {"free": free, "pinned": pinned}
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
