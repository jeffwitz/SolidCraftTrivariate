#!/usr/bin/env python3
"""Solve the P1/P2 cylinder in an exact axial-invariant subspace with SciPy.

Rust assembles the original 3D P1/P2 CSR operator and right hand side. This
oracle only applies the exact map u(x,y,z)=(ux(x,y),uy(x,y),0), factorizes the
resulting sparse 2D operator, and asks Rust to evaluate the reconstructed 3D
field and its residual in the original 3D free test space.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import resource
import subprocess
import time
from pathlib import Path

LIMIT_BYTES = 6 * 1024**3
resource.setrlimit(resource.RLIMIT_AS, (LIMIT_BYTES, LIMIT_BYTES))
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import scipy
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu


def run_command(command, timeout=1800):
    started = time.monotonic()
    env = os.environ.copy()
    result = subprocess.run(command, capture_output=True, text=True, timeout=timeout, env=env)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(map(str, command))}\n{result.stdout}\n{result.stderr}")
    return result.stdout, time.monotonic() - started


def vector_solve(binary: Path, work: Path, level: int, order: int, nu: float,
                 inner: float, outer: float, projected: bool):
    stem = f"h{level}-q{order}-nu{nu}-{'projected' if projected else 'full'}"
    export_path = work / f"{stem}.csr.json"
    solution_path = work / f"{stem}.solution.json"
    report_path = work / f"{stem}.readback.json"
    formulation = "projected" if projected else "full"
    base = [str(binary)]
    _, export_wall = run_command(base + ["plane-export", str(level), str(order), str(nu), str(inner), str(outer), formulation, str(export_path)])
    exported = json.loads(export_path.read_text())
    csr = exported["matrix_csr"]
    n = int(csr["dimension"])
    matrix = csr_matrix((np.asarray(csr["values"], dtype=np.float64),
                         np.asarray(csr["column_indices"], dtype=np.int64),
                         np.asarray(csr["row_ptr"], dtype=np.int64)), shape=(n, n))
    rhs = np.asarray(exported["rhs_xy"], dtype=np.float64)
    free = np.asarray(exported["free_vector_dofs"], dtype=np.int64)
    if matrix.shape != (n, n) or rhs.size != n or len(free) == 0:
        raise ValueError("invalid plane-strain CSR export")
    started = time.monotonic()
    factor = splu(matrix[free][:, free].tocsc())
    factorization_s = time.monotonic() - started
    started = time.monotonic()
    displacement = np.zeros(n)
    displacement[free] = factor.solve(rhs[free])
    free_matrix = matrix[free][:, free]
    rhs_free = rhs[free]
    refinement_steps = 0
    while refinement_steps < 8:
        reduced_residual = free_matrix @ displacement[free] - rhs_free
        residual = float(np.linalg.norm(reduced_residual) / max(np.linalg.norm(rhs_free), 1.0))
        if residual <= 1.0e-12:
            break
        displacement[free] -= factor.solve(reduced_residual)
        refinement_steps += 1
    # Report the residual of the exact vector written below, including a
    # correction applied by the final permitted refinement step.
    reduced_residual = free_matrix @ displacement[free] - rhs_free
    residual = float(np.linalg.norm(reduced_residual) / max(np.linalg.norm(rhs_free), 1.0))
    solve_s = time.monotonic() - started
    reasons=[]
    if not math.isfinite(residual) or residual > 1.0e-10:
        reasons.append("reduced SciPy residual exceeds 1e-10")
    solution_path.write_text(json.dumps({
        "schema": "plane_strain_solution_v1", "level": level, "order": order,
        "nu": nu, "inner_mpa": inner, "outer_mpa": outer,
        "formulation": exported["formulation"], "displacement_xy": displacement.tolist(),
        "relative_residual_reduced": residual,
    }, separators=(",", ":")) + "\n")
    _, read_wall = run_command(base + ["plane-read", str(level), str(order), str(nu), str(inner), str(outer), formulation, str(solution_path), str(report_path)])
    report = json.loads(report_path.read_text())
    if report["relative_residual_3d_free"] > 1.0e-10:
        reasons.append("reconstructed initial-3D free residual exceeds 1e-10")
    if report["energy_work_relative"] > 1.0e-8:
        reasons.append("energy/work identity exceeds 1e-8")
    mem_mib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    return {"status":"valid" if not reasons else "indeterminate", "indeterminate_reasons":reasons,
            "level": level, "order": order, "nu": nu, "pressures_mpa": [inner, outer],
            "formulation": exported["formulation"], "vector_dofs_2d": n,
            "scalar_dofs_2d": exported["scalar_dofs_2d"], "nnz_2d": matrix.nnz,
            "free_vector_dofs_2d": len(free), "reduced_residual": residual,
            "iterative_refinement_steps": refinement_steps,
            "relative_residual_3d_free": report["relative_residual_3d_free"],
            "physical_metrics": {k: v for k, v in report.items() if k in (
                "relative_l2", "relative_energy_norm", "physical_strain_energy_n_mm",
                "strain_energy_matrix_n_mm", "work_n_mm", "energy_work_relative",
                "pressure_resultant_n", "inner", "outer", "volume_mm3")},
            "timing_s": {"rust_preparation": exported["preparation_s"],
                         "rust_assembly": exported["assembly_s"],
                         "rust_reduction": exported["reduction_s"],
                         "export_process_wall": export_wall,
                         "scipy_factorization": factorization_s,
                         "scipy_backsolve": solve_s,
                         "rust_field_readback_wall": read_wall},
            "peak_rss_mib_campaign_process": mem_mib,
            "artifacts": {"csr": str(export_path), "solution": str(solution_path), "readback": str(report_path),
                          "readback_sha256": hashlib.sha256(report_path.read_bytes()).hexdigest()},
            "_displacement_mm_xyz": report["displacement_mm_xyz"]}


def compare_with_3d(binary, work, result, inner, outer):
    stem = f"3d-h{result['level']}-q{result['order']}-nu{result['nu']}-{'projected' if result['formulation']=='projected_bulk_q1' else 'full'}"
    path = work / f"{stem}.json"
    projected = result["formulation"] == "projected_bulk_q1"
    command = [str(binary), "field", str(result["level"]), str(result["order"]), str(result["nu"]), str(inner), str(outer), "false", "projected" if projected else "full", str(path)]
    _, elapsed = run_command(command)
    full = json.loads(path.read_text())
    u2 = np.asarray(result["_displacement_mm_xyz"], dtype=np.float64)
    u3 = np.asarray(full["displacement_mm_xyz"], dtype=np.float64)
    if u2.shape != u3.shape:
        raise ValueError("3D and axial-invariant coefficient vector dimensions differ")
    relative = float(np.linalg.norm(u2-u3)/max(np.linalg.norm(u3),1.0e-30))
    max_error = float(np.max(np.abs(u2-u3)))
    energy_relative = abs(result["physical_metrics"]["physical_strain_energy_n_mm"] - full["physical_strain_energy_n_mm"]) / max(abs(full["physical_strain_energy_n_mm"]),1.0)
    result["comparison_3d"] = {"relative_displacement_coefficient_error": relative,
        "max_abs_displacement_coefficient_error_mm": max_error,
        "relative_physical_energy_error": energy_relative,
        "three_dimensional_solve_wall_s": elapsed,
        "three_dimensional_iterations": full["pcg_iterations"],
        "three_dimensional_relative_residual": full["relative_residual"]}
    if relative > 1e-8 or energy_relative > 1e-8:
        result["status"]="indeterminate"
        result.setdefault("indeterminate_reasons",[]).append("h1/h2 direct 3D comparison exceeds 1e-8")
    del result["_displacement_mm_xyz"]


def failures_from_prior(path):
    data = json.loads(path.read_text())
    failures = []
    for case in data.get("cases", []):
        for run in case.get("runs", []):
            if run.get("status") != "completed" or not run.get("checks", {}).get("threshold_satisfied", False):
                failures.append({"nu": case.get("nu"), "pressures_mpa": case.get("pressures_mpa"),
                    "level": run.get("level"), "order": run.get("order"),
                    "projected_bulk": data.get("projected_bulk", False),
                    "status": run.get("status"), "checks": run.get("checks"),
                    "diagnostic": run.get("diagnostic")})
    return failures


def campaign(binary: Path, work: Path, output: Path, p1_path: Path, p2_path: Path):
    work.mkdir(parents=True, exist_ok=True)
    runs = []
    for nu in (0.49, 0.4999):
        for level in (1, 2):
            for projected in (False, True):
                item = vector_solve(binary, work, level, 8, nu, 100.0, 0.0, projected)
                compare_with_3d(binary, work, item, 100.0, 0.0)
                runs.append(item)
    comparisons = [{"level": r["level"], "nu": r["nu"], "formulation": r["formulation"], **r["comparison_3d"]} for r in runs]
    for nu in (0.49, 0.4999):
        for order in (8, 12):
            for projected in (False, True):
                item = vector_solve(binary, work, 4, order, nu, 100.0, 0.0, projected)
                item.pop("_displacement_mm_xyz", None)
                runs.append(item)
    failures = {"P1": failures_from_prior(p1_path), "P2": failures_from_prior(p2_path)}
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {"schema": "plane_strain_lame_cylinder_campaign_v1",
        "scope": "constant-height annulus under z-invariant pressure; not a general 3D reduction",
        "method": "exact axial-invariant displacement map, 3D CSR reassembly, sparse SciPy LU",
        "material_young_mpa": 210000.0, "pressures_mpa": [100.0,0.0],
        "environment": {"python": __import__("platform").python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
            "memory_budget_mib": 6144, "timeout_seconds_per_child": 1800},
        "h1_h2_3d_equivalence": comparisons, "h4_runs": runs[8:],
        "prior_p1_p2_failures_preserved": failures,
        "qualification": {"h1_h2_displacement_and_energy_relative": "<=1e-8",
            "reconstructed_3d_free_residual": "<=1e-10",
            "status": "valid" if all(r["status"]=="valid" for r in runs) else "indeterminate",
            "indeterminate_cases":[{"level":r["level"],"order":r["order"],"nu":r["nu"],"formulation":r["formulation"],"reasons":r.get("indeterminate_reasons",[])} for r in runs if r["status"]!="valid"]}}
    output.write_text(json.dumps(payload, indent=2, allow_nan=False)+"\n")
    print(json.dumps(payload, indent=2, allow_nan=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary",type=Path)
    parser.add_argument("work",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--p1",type=Path,default=Path("docs/results/pressure-cylinder-20261125.json"))
    parser.add_argument("--p2",type=Path,default=Path("docs/results/pressure-cylinder-projected-20261125.json"))
    args=parser.parse_args()
    started=time.monotonic()
    campaign(args.binary.resolve(),args.work,args.output,args.p1,args.p2)
    print(f"elapsed_seconds={time.monotonic()-started:.3f}")


if __name__ == "__main__":
    main()
