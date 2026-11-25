#!/usr/bin/env python3
"""P5 residual diagnostic for archived cylinder fields and the original K3.

SciPy only produces diagnostic correction candidates in the already mounted
plane-strain space. The Rust ``plane-diagnose`` route independently rebuilds
the original three-component elasticity operator and evaluates every field;
only its residual and the normal Rust readback can establish qualification.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import resource
import subprocess
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu

LIMIT_BYTES = 6 * 1024**3
resource.setrlimit(resource.RLIMIT_AS, (LIMIT_BYTES, LIMIT_BYTES))
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False, sort_keys=True) + "\n")


def run(command: list[str], timeout: int = 1800) -> float:
    start = time.monotonic()
    result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}\n{result.stdout}\n{result.stderr}")
    return time.monotonic() - start


def matrix_from_export(value: dict) -> csr_matrix:
    c = value["matrix_csr"]
    n = int(c["dimension"])
    matrix = csr_matrix((np.asarray(c["values"], dtype=np.float64),
                         np.asarray(c["column_indices"], dtype=np.int64),
                         np.asarray(c["row_ptr"], dtype=np.int64)), shape=(n, n))
    if matrix.shape != (n, n) or matrix.nnz == 0:
        raise ValueError("invalid archived plane-strain CSR")
    return matrix


def residuals(matrix: csr_matrix, x: np.ndarray, rhs: np.ndarray) -> dict[str, float]:
    plain = matrix @ x - rhs
    compensated = np.empty(rhs.size, dtype=np.float64)
    extended = np.empty(rhs.size, dtype=np.longdouble)
    lx = x.astype(np.longdouble)
    ld = matrix.data.astype(np.longdouble)
    lrhs = rhs.astype(np.longdouble)
    for row in range(rhs.size):
        lo, hi = matrix.indptr[row : row + 2]
        ids = matrix.indices[lo:hi]
        compensated[row] = math.fsum(float(a) * float(b) for a, b in zip(matrix.data[lo:hi], x[ids])) - rhs[row]
        extended[row] = np.sum(ld[lo:hi] * lx[ids], dtype=np.longdouble) - lrhs[row]

    def relnorm(value: np.ndarray) -> float:
        return float(np.sqrt(np.sum(value.astype(np.longdouble) ** 2, dtype=np.longdouble)) /
                     max(float(np.sqrt(np.sum(rhs.astype(np.longdouble) ** 2, dtype=np.longdouble))), 1.0))

    return {"relative_plain_f64": relnorm(plain),
            "relative_fsum_products_f64": relnorm(compensated),
            "relative_longdouble_products_and_sum": relnorm(extended),
            "max_abs_plain_vs_longdouble": float(np.max(np.abs(plain.astype(np.longdouble) - extended)))}


def candidate_sequence(work: Path, level: int, order: int, nu: float, formulation: str) -> tuple[dict, dict, list[dict], list[dict], dict]:
    stem = f"h{level}-q{order}-nu{nu}-{formulation}"
    exported = json.loads((work / f"{stem}.csr.json").read_text())
    solution = json.loads((work / f"{stem}.solution.json").read_text())
    if exported["schema"] != "plane_strain_sparse_system_v1" or solution["schema"] != "plane_strain_solution_v1":
        raise ValueError(f"unexpected schemas for {stem}")
    matrix = matrix_from_export(exported)
    n = matrix.shape[0]
    free = np.asarray(exported["free_vector_dofs"], dtype=np.int64)
    rhs = np.asarray(exported["rhs_xy"], dtype=np.float64)
    if rhs.size != n or free.size == 0 or solution["formulation"] != exported["formulation"]:
        raise ValueError(f"inconsistent archived system {stem}")
    free_matrix = matrix[free][:, free].tocsc()
    factor_start = time.monotonic()
    factor = splu(free_matrix)
    factor_s = time.monotonic() - factor_start
    current = np.asarray(solution["displacement_xy"], dtype=np.float64)
    if current.size != n or not np.isfinite(current).all():
        raise ValueError(f"invalid archived displacement {stem}")
    reduced_rhs = rhs[free]
    free_csr = free_matrix.tocsr()
    rows = []
    candidates = []
    # Correction zero is the published P4 field. Each later vector applies
    # one correction from the same sparse factorization, using long-double
    # reduced residual evaluation before converting the correction RHS to f64.
    for step in range(7):
        measures = residuals(free_csr, current[free], reduced_rhs)
        rows.append({"corrections_applied": step, **measures})
        candidates.append({"refinement": step, "displacement_xy": current.tolist()})
        if step < 6:
            r_ext = np.empty(reduced_rhs.size, dtype=np.longdouble)
            x_ext = current[free].astype(np.longdouble)
            data_ext = free_csr.data.astype(np.longdouble)
            rhs_ext = reduced_rhs.astype(np.longdouble)
            for i in range(free_csr.shape[0]):
                lo, hi = free_csr.indptr[i : i + 2]
                ids = free_csr.indices[lo:hi]
                r_ext[i] = rhs_ext[i] - np.sum(data_ext[lo:hi] * x_ext[ids], dtype=np.longdouble)
            delta = factor.solve(np.asarray(r_ext, dtype=np.float64))
            current = current.copy()
            current[free] += delta
    timing = {"scipy_factorization_s": factor_s,
              "scipy_candidate_refinement_s": time.monotonic() - factor_start - factor_s}
    return exported, solution, rows, candidates, timing


def diagnose_case(binary: Path, work: Path, output_dir: Path, level: int, order: int, nu: float, formulation: str) -> dict:
    stem = f"h{level}-q{order}-nu{nu}-{formulation}"
    exported, base_solution, reduced_rows, candidates, timing = candidate_sequence(work, level, order, nu, formulation)
    candidates_path = output_dir / f"{stem}-candidates.json"
    diagnostic_path = output_dir / f"{stem}-k3-diagnostic.json"
    write_json(candidates_path, {"schema": "plane_residual_candidates_v1", "level": level,
                                "order": order, "nu": nu, "inner_mpa": 100.0, "outer_mpa": 0.0,
                                "formulation": exported["formulation"], "candidates": candidates})
    rust_diagnostic_s = 0.0
    if not diagnostic_path.is_file():
        rust_diagnostic_s = run([str(binary), "plane-diagnose", str(level), str(order), str(nu), "100", "0",
                                 "projected" if formulation == "projected" else "full",
                                 str(candidates_path), str(diagnostic_path)])
    diagnostic = json.loads(diagnostic_path.read_text())
    if diagnostic["schema"] != "full_cylinder_residual_diagnostic_v1" or len(diagnostic["rows"]) != 7:
        raise ValueError(f"invalid Rust K3 diagnostic for {stem}")
    previous = json.loads((output_dir / f"{stem}-previous-readback.json").read_text()) if (output_dir / f"{stem}-previous-readback.json").is_file() else None
    p4_path = work / f"{stem}.readback.json"
    if p4_path.is_file():
        previous = json.loads(p4_path.read_text())
    rows = []
    for reduced, k3 in zip(reduced_rows, diagnostic["rows"]):
        rows.append({**reduced,
                     "k3_relative_residual_naive_full_xyz": k3["relative_residual_naive_full_xyz"],
                     "k3_relative_residual_compensated_full_xyz": k3["relative_residual_compensated_full_xyz"],
                     "k3_relative_residual_naive_by_component": k3["relative_residual_naive_by_component"],
                     "k3_relative_residual_compensated_by_component": k3["relative_residual_compensated_by_component"],
                     "k3_max_abs_naive_compensated_difference_n": k3["max_abs_naive_compensated_difference_n"],
                     "k3_external_force_n": k3["external_force_n"],
                     "k3_force_balance_n": k3["force_balance_n"],
                     "k3_external_moment_n_mm": k3["external_moment_n_mm"],
                     "k3_moment_balance_n_mm": k3["moment_balance_n_mm"]})
    accepted = [row["corrections_applied"] for row in rows if row["k3_relative_residual_compensated_full_xyz"] <= 1e-10]
    selected = min(accepted) if accepted else min(rows, key=lambda row: row["k3_relative_residual_compensated_full_xyz"])["corrections_applied"]
    final_solution = {"schema": "plane_strain_solution_v1", "level": level, "order": order, "nu": nu,
                      "inner_mpa": 100.0, "outer_mpa": 0.0, "formulation": exported["formulation"],
                      "displacement_xy": candidates[selected]["displacement_xy"],
                      "relative_residual_reduced": rows[selected]["relative_longdouble_products_and_sum"]}
    solution_path = output_dir / f"{stem}-selected.solution.json"
    readback_path = output_dir / f"{stem}-selected.readback.json"
    write_json(solution_path, final_solution)
    readback_s = run([str(binary), "plane-read", str(level), str(order), str(nu), "100", "0",
                      "projected" if formulation == "projected" else "full", str(solution_path), str(readback_path)])
    readback = json.loads(readback_path.read_text())
    force_balance_relative = float(np.linalg.norm(rows[selected]["k3_force_balance_n"]) /
                                   max(np.linalg.norm(rows[selected]["k3_external_force_n"]), 1.0))
    moment_balance_relative = float(np.linalg.norm(rows[selected]["k3_moment_balance_n_mm"]) /
                                     max(np.linalg.norm(rows[selected]["k3_external_moment_n_mm"]), 1.0))
    maximum_equilibrium = max(force_balance_relative, moment_balance_relative)
    k3_residual = rows[selected]["k3_relative_residual_compensated_full_xyz"]
    status = "valid" if k3_residual <= 1e-10 and maximum_equilibrium <= 1e-8 and readback["energy_work_relative"] <= 1e-8 else "indeterminate"
    return {"id": stem, "level": level, "order": order, "nu": nu, "formulation": exported["formulation"],
            "status": status, "selected_corrections": selected, "candidate_metrics": rows,
            "previous_p4_status": previous.get("status") if previous else None,
            "previous_p4_residual_3d_free": previous.get("relative_residual_3d_free") if previous else None,
            "selected_rust_k3_residual_all_free_xyz": k3_residual,
            "selected_rust_plane_readback_residual_xy": readback["relative_residual_3d_free"],
            "selected_force_balance_relative": force_balance_relative,
            "selected_moment_balance_relative": moment_balance_relative,
            "selected_physical_metrics": {k: readback[k] for k in ("relative_l2", "relative_energy_norm", "physical_strain_energy_n_mm", "strain_energy_matrix_n_mm", "work_n_mm", "energy_work_relative", "pressure_resultant_n", "inner", "outer", "volume_mm3")},
            "timing_s": {**timing, "rust_k3_diagnostic_s": rust_diagnostic_s, "rust_selected_readback_s": readback_s},
            "candidate_file": str(candidates_path), "diagnostic_file": str(diagnostic_path),
            "selected_solution_file": str(solution_path), "selected_readback_file": str(readback_path)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    progress = args.output_dir / "progress.json"
    cases = json.loads(progress.read_text())["cases"] if progress.is_file() else []
    completed = {case["id"] for case in cases}
    plan = [
        (5, 8, "projected"),
        (5, 8, "full"),
        (5, 12, "full"),
        (5, 12, "projected"),
        *[(6, order, formulation)
          for order in (8, 12)
          for formulation in ("full", "projected")],
    ]
    for level, order, formulation in plan:
        stem = f"h{level}-q{order}-nu0.4999-{formulation}"
        if stem in completed:
            continue
        case = diagnose_case(args.binary.resolve(), args.work, args.output_dir,
                             level, order, 0.4999, formulation)
        cases.append(case)
        write_json(progress, {"cases": cases})
    plan_order = {f"h{level}-q{order}-nu0.4999-{formulation}": index
                  for index, (level, order, formulation) in enumerate(plan)}
    cases.sort(key=lambda case: plan_order[case["id"]])
    payload = {
        "schema": "pressure_cylinder_residual_diagnostic_v1",
        "protocol": "P5, frozen before numerical runs",
        "status": "valid" if all(case["status"] == "valid" for case in cases) else "indeterminate",
        "method": "SciPy corrections in archived Kxy factors; Rust rebuilds original K3 once per case and compares plain sparse f64 residuals with its FMA-compensated product/addition residual; Rust plane-read independently checks the published field and physical objectives.",
        "material": {"young_mpa": 210000.0, "nu": 0.4999, "pressure_inner_mpa": 100.0, "pressure_outer_mpa": 0.0},
        "limits": {"residual_threshold": 1e-10, "equilibrium_threshold": 1e-8, "energy_work_threshold": 1e-8, "memory_budget_bytes": LIMIT_BYTES, "timeout_per_case_seconds": 1800},
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "correction_method": "single archived SciPy sparse LU factor reused for corrections 0..6; long-double residual of Kxy is converted to f64 for each solve; no candidate is adopted unless Rust K3 and readback satisfy unchanged limits.",
        "cases": cases,
        "limitations": ["extended precision is evaluated independently for the reduced Kxy matrix only; K3 is checked in Rust using plain f64 versus compensated FMA summation", "the diagnostic candidate generator is not a production Rust solver", "high-nu fields that fail the original K3 residual remain indeterminate"],
    }
    write_json(args.report, payload)
    print(json.dumps({"status": payload["status"], "cases": [(x["id"], x["status"], x["selected_rust_k3_residual_all_free_xyz"]) for x in cases]}, allow_nan=False))


if __name__ == "__main__":
    main()
