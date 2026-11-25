#!/usr/bin/env python3
"""Independent SciPy full-3D correction runs for the pressure-cylinder P3 study."""

from __future__ import annotations

import argparse
import json
import math
import os
import resource
import subprocess
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


def command(args: list[str], timeout: int = 1800) -> float:
    start = time.monotonic()
    result = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(args)}\n{result.stdout}\n{result.stderr}")
    return time.monotonic() - start


def compensated_residual(matrix, x: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Compute each CSR row dot product with math.fsum compensation."""
    result = np.empty(rhs.size)
    for row in range(rhs.size):
        start, end = matrix.indptr[row : row + 2]
        result[row] = math.fsum(float(a) * float(b) for a, b in zip(matrix.data[start:end], x[matrix.indices[start:end]])) - rhs[row]
    return result


def solve_one(binary: Path, work: Path, level: int, order: int, nu: float, projected: bool, plane_index: dict):
    formulation = "projected" if projected else "full"
    name = f"h{level}-q{order}-nu{nu}-{formulation}"
    export_path = work / f"{name}.csr.json"
    solution_path = work / f"{name}.solution.json"
    readback_path = work / f"{name}.readback.json"
    export_wall = command([str(binary), "full-export", str(level), str(order), str(nu), "100", "0", formulation, str(export_path)])
    exported = json.loads(export_path.read_text())
    block = exported["matrix_block_csr"]
    block_data = np.asarray(block["values_3x3"], dtype=np.float64)
    matrix = bsr_matrix((block_data, np.asarray(block["column_indices"], dtype=np.int32),
                         np.asarray(block["row_ptr"], dtype=np.int32)),
                        shape=(3 * int(block["block_rows"]),) * 2).tocsr()
    rhs = np.asarray(exported["rhs_xyz"], dtype=np.float64)
    fixed = np.asarray(exported["fixed_vector_dofs"], dtype=np.int64)
    free_mask = np.ones(rhs.size, dtype=bool)
    free_mask[fixed] = False
    free = np.flatnonzero(free_mask)
    reduced = matrix[free][:, free].tocsc()
    reduced_rhs = rhs[free]
    factor_start = time.monotonic()
    factor = splu(reduced)
    factor_s = time.monotonic() - factor_start
    solve_start = time.monotonic()
    displacement = np.zeros(rhs.size)
    displacement[free] = factor.solve(reduced_rhs)
    refinement = 0
    residual_vector = compensated_residual(matrix, displacement, rhs)
    residual = float(np.linalg.norm(residual_vector[free]) / max(np.linalg.norm(rhs[free]), 1.0))
    while residual > 1e-12 and refinement < 8:
        correction_rhs = residual_vector[free]
        displacement[free] -= factor.solve(correction_rhs)
        refinement += 1
        residual_vector = compensated_residual(matrix, displacement, rhs)
        residual = float(np.linalg.norm(residual_vector[free]) / max(np.linalg.norm(rhs[free]), 1.0))
    solve_s = time.monotonic() - solve_start
    residual_vector = compensated_residual(matrix, displacement, rhs)
    residual = float(np.linalg.norm(residual_vector[free]) / max(np.linalg.norm(rhs[free]), 1.0))
    solution_path.write_text(json.dumps({
        "schema": "full_cylinder_solution_v1", "level": level, "order": order,
        "nu": nu, "inner_mpa": 100.0, "outer_mpa": 0.0,
        "formulation": exported["formulation"], "displacement_xyz": displacement.tolist(),
        "relative_residual_compensated": residual,
    }, separators=(",", ":")) + "\n")
    readback_wall = command([str(binary), "full-read", str(level), str(order), str(nu), "100", "0", formulation, str(solution_path), str(readback_path)])
    readback = json.loads(readback_path.read_text())
    plane_case = plane_index[(level, order, nu, exported["formulation"])]
    plane = json.loads(Path(plane_case["artifacts"]["readback"]).read_text())
    plane_u = np.asarray(plane["displacement_mm_xyz"], dtype=np.float64)
    difference = displacement - plane_u
    record = {
        "status": "valid" if residual <= 1e-10 and readback["relative_residual_3d_free"] <= 1e-10 else "indeterminate",
        "level": level, "order": order, "nu": nu, "formulation": exported["formulation"],
        "vector_dofs_3d": int(rhs.size), "nnz_3d": int(matrix.nnz), "free_vector_dofs_3d": int(free.size),
        "relative_residual_compensated_original_3d": residual,
        "relative_residual_rust_readback": readback["relative_residual_3d_free"],
        "plane_strain_comparison": {
            "relative_displacement_coefficient_error": float(np.linalg.norm(difference) / max(np.linalg.norm(displacement), 1e-30)),
            "max_abs_displacement_coefficient_error_mm": float(np.max(np.abs(difference))),
            "relative_physical_energy_difference": abs(readback["physical_strain_energy_n_mm"] - plane["physical_strain_energy_n_mm"]) / max(abs(readback["physical_strain_energy_n_mm"]), 1.0),
        },
        "physical_metrics": {key: readback[key] for key in ("relative_l2", "relative_energy_norm", "physical_strain_energy_n_mm", "work_n_mm", "strain_energy_matrix_n_mm", "energy_work_relative", "pressure_resultant_n", "inner", "outer", "volume_mm3")},
        "timing_s": {"rust_preparation": exported["preparation_s"], "rust_assembly": exported["assembly_s"], "rust_export_wall": export_wall, "scipy_factorization": factor_s, "scipy_solve_and_refinement": solve_s, "rust_readback_wall": readback_wall},
        "iterative_refinement_steps": refinement,
        "peak_rss_mib_scipy_process": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
        "peak_rss_mib_child_process": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0,
        "artifacts": {"csr": str(export_path), "solution": str(solution_path), "readback": str(readback_path)},
    }
    return record


def add_equilibrium_checks(record: dict) -> None:
    """Recompute compensated free residual and resultant balance from saved CSR."""
    artifacts = record["artifacts"]
    exported = json.loads(Path(artifacts["csr"]).read_text())
    solved = json.loads(Path(artifacts["solution"]).read_text())
    b = exported["matrix_block_csr"]
    matrix = bsr_matrix((np.asarray(b["values_3x3"], dtype=np.float64),
                         np.asarray(b["column_indices"], dtype=np.int32),
                         np.asarray(b["row_ptr"], dtype=np.int32)),
                        shape=(3 * int(b["block_rows"]),) * 2).tocsr()
    rhs = np.asarray(exported["rhs_xyz"], dtype=np.float64)
    u = np.asarray(solved["displacement_xyz"], dtype=np.float64)
    fixed = np.asarray(exported["fixed_vector_dofs"], dtype=np.int64)
    residual = compensated_residual(matrix, u, rhs)
    fixed_mask = np.zeros(rhs.size, dtype=bool)
    fixed_mask[fixed] = True
    free = ~fixed_mask
    residual_relative = float(np.linalg.norm(residual[free]) / max(np.linalg.norm(rhs[free]), 1.0))
    external = rhs.reshape(-1, 3).sum(axis=0)
    reaction = np.zeros(3)
    for c in range(3):
        reaction[c] = residual[fixed[fixed % 3 == c]].sum()
    balance = float(np.linalg.norm(external + reaction) / max(np.linalg.norm(external), 1.0))
    record["relative_residual_compensated_original_3d"] = residual_relative
    record["equilibrium"] = {"external_resultant_n": external.tolist(),
        "support_reaction_resultant_n": reaction.tolist(), "force_balance_relative": balance}
    record["status"] = "valid" if residual_relative <= 1e-10 and balance <= 1e-8 and record["physical_metrics"]["energy_work_relative"] <= 1e-8 else "indeterminate"


def plane_only_extensions(binary: Path, work: Path) -> list[dict]:
    """Run h5 plane reductions; h6 is gated by the recorded fill estimate."""
    from plane_strain_cylinder import vector_solve

    progress = work / "plane-extension-progress.json"
    runs = json.loads(progress.read_text()) if progress.is_file() else []
    done = {(r.get("level"), r.get("order"), r.get("nu"), r.get("formulation")) for r in runs}
    for order in (8, 12):
        for nu in (0.49, 0.4999):
            for projected in (False, True):
                formulation = "projected_bulk_q1" if projected else "displacement"
                if (5, order, nu, formulation) in done:
                    continue
                try:
                    result = vector_solve(binary, work, 5, order, nu, 100.0, 0.0, projected)
                    result.pop("_displacement_mm_xyz", None)
                except (RuntimeError, MemoryError, subprocess.TimeoutExpired) as error:
                    result = {"status": "indeterminate", "level": 5, "order": order, "nu": nu,
                              "formulation": formulation, "reason": str(error),
                              "memory_budget_bytes": LIMIT_BYTES, "timeout_seconds": 1800}
                runs.append(result)
                progress.write_text(json.dumps(runs, indent=2, allow_nan=False) + "\n")
    return runs


def level_six_pilot(binary: Path, work: Path) -> dict:
    """Run the requested first h6 plane-only solve under the child 6 GiB cap."""
    from plane_strain_cylinder import vector_solve

    checkpoint = work / "h6-pilot.json"
    if checkpoint.is_file():
        return json.loads(checkpoint.read_text())
    try:
        run = vector_solve(binary, work, 6, 8, 0.49, 100.0, 0.0, False)
        run.pop("_displacement_mm_xyz", None)
    except (RuntimeError, MemoryError, subprocess.TimeoutExpired) as error:
        run = {"status": "indeterminate", "level": 6, "order": 8, "nu": 0.49,
               "formulation": "displacement", "reason": str(error),
               "memory_budget_bytes": LIMIT_BYTES, "timeout_seconds": 1800}
    checkpoint.write_text(json.dumps(run, indent=2, allow_nan=False) + "\n")
    return run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("plane_report", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    plane = json.loads(args.plane_report.read_text())
    plane_index = {(r["level"], r["order"], r["nu"], r["formulation"]): r for r in plane["h4_runs"]}
    progress = args.work / "progress.json"
    runs = json.loads(progress.read_text()) if progress.is_file() else []
    done = {(r.get("level"), r.get("order"), r.get("nu"), r.get("formulation")) for r in runs}
    for nu, order in ((0.4999, 8), (0.4999, 12), (0.49, 8)):
        for projected in (False, True):
            formulation = "projected_bulk_q1" if projected else "displacement"
            if (4, order, nu, formulation) in done:
                continue
            try:
                runs.append(solve_one(args.binary.resolve(), args.work, 4, order, nu, projected, plane_index))
            except (RuntimeError, MemoryError, subprocess.TimeoutExpired) as error:
                runs.append({"status": "indeterminate", "level": 4, "order": order, "nu": nu,
                             "formulation": formulation,
                             "reason": str(error), "budget_bytes": LIMIT_BYTES, "timeout_seconds": 1800})
            progress.write_text(json.dumps(runs, indent=2, allow_nan=False) + "\n")
    for run in runs:
        if run.get("artifacts") and "equilibrium" not in run:
            add_equilibrium_checks(run)
    extensions = plane_only_extensions(args.binary.resolve(), args.work)
    pilot = level_six_pilot(args.binary.resolve(), args.work)
    payload = {
        "schema": "plane_strain_full_3d_correction_v1",
        "scope": "full 3D Q2 pressure-cylinder solve at h4; comparison only against the exact axial-invariant subspace",
        "method": "independent SciPy sparse LU on the original assembled 3D CSR, compensated row residual, iterative refinement",
        "constraints": "same axial end fixation and in-plane symmetry gauges as the P3 cylinder example",
        "load_mpa": {"inner": 100.0, "outer": 0.0},
        "material_young_mpa": 210000.0,
        "environment": {"python": __import__("platform").python_version(), "numpy": np.version.version, "scipy": scipy.__version__, "memory_budget_mib": LIMIT_BYTES // 1024**2, "timeout_per_case_seconds": 1800},
        "runs": runs,
        "plane_only_extensions": extensions,
        "plane_only_extension_preflight": {
            "method": "measured h4 q8 displacement LU nonzeros scaled with 2D tensor-grid N^1.5; 16 bytes per predicted LU nonzero, then 2x workspace allowance",
            "h4_scalar_vector_dofs": 8448,
            "h4_lu_nonzeros": 5091073,
            "h4_measured_peak_rss_mib": 194.4,
            "h5_predicted_lu_storage_mib": 621.5,
            "h5_predicted_with_workspace_mib": 1243.0,
            "h6_predicted_lu_storage_mib": 4971.8,
            "h6_predicted_with_workspace_mib": 9943.5,
            "h5_included": True,
            "h6_skipped_over_6_gib_estimate": True,
        },
        "updated_h6_preflight": {
            "method": "measured h4 and h5 sparse LU fill; extrapolate the observed level ratio, 16 bytes per L/U nonzero and 2x workspace allowance",
            "h4_vector_dofs": 8448,
            "h4_lu_nonzeros": 5091073,
            "h4_peak_rss_mib": 194.5,
            "h5_vector_dofs": 33280,
            "h5_lu_nonzeros": 25426109,
            "h5_peak_rss_mib": 598.4,
            "observed_lu_growth_h5_over_h4": 4.9943,
            "h6_predicted_lu_nonzeros": 126984433,
            "h6_predicted_factor_storage_mib": 1937.6,
            "h6_predicted_with_workspace_mib": 3875.3,
            "memory_budget_mib": LIMIT_BYTES // 1024**2,
            "h6_pilot_launched": True,
        },
        "h6_plane_only_pilot_h4_q8_nu049_full": pilot,
        "qualification": {"residual_threshold": 1e-10, "status": "valid" if all(r["status"] == "valid" for r in runs) else "indeterminate"},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps(payload, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
