#!/usr/bin/env python3
"""Independent volume/traction equilibrium check for a cubic NURBS strip.

The CLI supplies an audited model and its production C0 Q3 numbering. This
oracle independently evaluates rational Q3 geometry and basis derivatives in
NumPy, then compares the assembled constant-stress volume term with traction
integrals over every selected exterior face. No linear system is solved.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np


SIGMA = np.array(
    [[2.0, 0.4, -0.3], [0.4, 3.0, 0.2], [-0.3, 0.2, 4.0]],
    dtype=np.float64,
)


def basis(t: float) -> tuple[np.ndarray, np.ndarray]:
    s = 1.0 - t
    return (
        np.array([s**3, 3 * t * s**2, 3 * t**2 * s, t**3]),
        np.array([-3 * s**2, 3 * s * (1 - 3 * t), 3 * t * (2 - 3 * t), 3 * t**2]),
    )


def quadrature(order: int) -> tuple[np.ndarray, np.ndarray]:
    points, weights = np.polynomial.legendre.leggauss(order)
    return (points + 1.0) / 2.0, weights / 2.0


def local_data(block: dict, q: tuple[float, float, float]):
    controls = np.asarray(block["controls"], dtype=np.float64).reshape((4, 4, 4, 4))
    values, derivatives = zip(*(basis(t) for t in q))
    b = np.einsum("i,j,k->ijk", *values)
    db = np.stack(
        (
            np.einsum("i,j,k->ijk", derivatives[0], values[1], values[2]),
            np.einsum("i,j,k->ijk", values[0], derivatives[1], values[2]),
            np.einsum("i,j,k->ijk", values[0], values[1], derivatives[2]),
        ),
        axis=0,
    )
    weights = controls[..., 3]
    weighted_basis = b * weights
    denominator = float(weighted_basis.sum())
    if denominator <= 0 or not math.isfinite(denominator):
        raise ArithmeticError("nonpositive rational basis denominator")
    ddenominator = np.einsum("aijk,ijk->a", db, weights)
    # Controls are already homogeneous [wx, wy, wz, w]. Apply the polynomial
    # basis once here; rational weights belong in the denominator and basis.
    homogeneous = np.einsum("ijk,ijkc->c", b, controls)
    dhomogeneous = np.einsum("aijk,ijkc->ac", db, controls)
    position = homogeneous[:3] / denominator
    jac = (
        dhomogeneous[:, :3] * denominator
        - ddenominator[:, None] * homogeneous[None, :3]
    ) / denominator**2
    jac = jac.T  # physical-coordinate rows, parametric-coordinate columns
    det = float(np.linalg.det(jac))
    if not math.isfinite(det) or det <= 0:
        raise ArithmeticError("nonpositive quadrature Jacobian")
    rational_basis = weighted_basis / denominator
    drational_param = (
        db * weights[None, ...] * denominator
        - weighted_basis[None, ...] * ddenominator[:, None, None, None]
    ) / denominator**2
    gradients = np.linalg.solve(jac.T, drational_param.reshape(3, -1)).T
    dofs = np.asarray(block["dofs"], dtype=np.int64)
    return position, jac, det, rational_basis.reshape(-1), gradients, dofs


def volume_contribution(block: dict, order: int, ndof: int) -> np.ndarray:
    result = np.zeros((ndof, 3), dtype=np.float64)
    points, weights = quadrature(order)
    for i, u in enumerate(points):
        for j, v in enumerate(points):
            for k, w in enumerate(points):
                _, _, det, _, gradients, dofs = local_data(block, (u, v, w))
                measure = weights[i] * weights[j] * weights[k] * det
                local = gradients @ SIGMA.T * measure
                np.add.at(result, dofs, local)
    return result


def face_contribution(block: dict, side: str, order: int, ndof: int) -> np.ndarray:
    result = np.zeros((ndof, 3), dtype=np.float64)
    axis = {"u": 0, "v": 1, "w": 2}[side[0]]
    upper = side[1] == "1"
    face_axes = {0: (1, 2), 1: (0, 2), 2: (0, 1)}[axis]
    points, weights = quadrature(order)
    for i, a in enumerate(points):
        for j, b in enumerate(points):
            q = [0.0, 0.0, 0.0]
            q[axis] = 1.0 if upper else 0.0
            q[face_axes[0]], q[face_axes[1]] = a, b
            _, jac, _, shape, _, dofs = local_data(block, tuple(q))
            tangent = [jac[:, a] for a in range(3)]
            if axis == 0:
                area = np.cross(tangent[1], tangent[2])
            elif axis == 1:
                area = np.cross(tangent[2], tangent[0])
            else:
                area = np.cross(tangent[0], tangent[1])
            if not upper:
                area = -area
            traction_times_area = SIGMA @ area
            measure = weights[i] * weights[j]
            local = shape[:, None] * traction_times_area[None, :] * measure
            np.add.at(result, dofs, local)
    return result


def integrate(package: dict, order: int) -> dict:
    ndof = package["scalar_ndof"]
    blocks = {block["id"]: block for block in package["blocks"]}
    source = package["source_model"]
    volume = np.zeros((ndof, 3), dtype=np.float64)
    surface = np.zeros((ndof, 3), dtype=np.float64)
    for block in blocks.values():
        volume += volume_contribution(block, order, ndof)
    seen = set()
    for face in source["design_faces"]:
        key = (face["block"], face["side"])
        if key in seen:
            continue
        seen.add(key)
        surface += face_contribution(blocks[face["block"]], face["side"], order, ndof)
    residual = volume - surface
    scale = max(float(np.max(np.abs(volume))), float(np.max(np.abs(surface))), 1e-300)
    return {
        "order": order,
        "exterior_faces": len(seen),
        "volume_term_linf": float(np.max(np.abs(volume))),
        "traction_term_linf": float(np.max(np.abs(surface))),
        "residual_linf": float(np.max(np.abs(residual))),
        "relative_residual_linf": float(np.max(np.abs(residual))) / scale,
        "volume_resultant": volume.sum(axis=0).tolist(),
        "traction_resultant": surface.sum(axis=0).tolist(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    start = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="nurbs-equilibrium-") as temp_dir:
        model_path = Path(temp_dir) / "strip.model.json"
        package_path = Path(temp_dir) / "strip.analysis.json"
        subprocess.run([str(args.cli), "build", str(args.recipe), str(model_path)], check=True)
        subprocess.run(
            [str(args.cli), "export-nurbs-space", str(model_path), str(package_path)],
            check=True,
        )
        package = json.loads(package_path.read_text())
    cases = [integrate(package, order) for order in (4, 8, 16)]
    if cases[-1]["relative_residual_linf"] > 1e-8:
        raise SystemExit(f"q16 equilibrium residual exceeds threshold: {cases[-1]}")
    report = {
        "schema": "nurbs_constant_stress_patch_equilibrium_v1",
        "status": "verified",
        "method": "independent NumPy rational-Q3 volume and exterior-face quadrature; no solve",
        "recipe": str(args.recipe),
        "model_blocks": len(package["blocks"]),
        "scalar_dofs": package["scalar_ndof"],
        "tensor_degree": package["tensor_degree"],
        "stress_mpa": SIGMA.tolist(),
        "traction": "sigma*n on every distinct selected exterior face",
        "volume_term": "integral grad(N_i) dot sigma over the volume",
        "free_surface_volume_force": "-div(sigma)=0 for constant sigma",
        "orders": cases,
        "q16_threshold_relative_linf": 1e-8,
        "elapsed_seconds": time.perf_counter() - start,
        "limitations": [
            "floating-point quadrature convergence is a diagnostic, not an exact integral proof",
            "the production model, audit, Q3 elevation, and global DOF numbering are input artifacts; integral evaluation is independent",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "output": str(args.output), "q16": cases[-1]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
