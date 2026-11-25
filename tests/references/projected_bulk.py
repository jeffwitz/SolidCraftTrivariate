#!/usr/bin/env python3
"""Independent SciPy oracle for projected-bulk Q2 elasticity fixtures.

The fixture contains serialized rational Q2 geometry and shared scalar DOF maps.
This script rebuilds geometry, deviatoric and full stiffness, Q1 pressure mass,
and divergence coupling without using the Rust matrices. It also solves the
uncondensed mixed saddle system for a manufactured deterministic load.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh, solve
from scipy.sparse import bmat, block_diag, coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve


def bernstein2(t: float) -> tuple[np.ndarray, np.ndarray]:
    """Degree-two Bernstein basis and derivative on [0, 1]."""
    s = 1.0 - t
    return np.array([s * s, 2.0 * s * t, t * t]), np.array([-2.0 * s, 2.0 - 4.0 * t, 2.0 * t])


def geometry(controls: list[list[int]], uvw: tuple[float, float, float]):
    """Evaluate a rational tensor-product Q2 map and its Jacobian."""
    b = [bernstein2(t) for t in uvw]
    poly = np.empty(27)
    dpoly = np.empty((27, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                n = (i * 3 + j) * 3 + k
                poly[n] = b[0][0][i] * b[1][0][j] * b[2][0][k]
                dpoly[n] = (
                    b[0][1][i] * b[1][0][j] * b[2][0][k],
                    b[0][0][i] * b[1][1][j] * b[2][0][k],
                    b[0][0][i] * b[1][0][j] * b[2][1][k],
                )

    h = np.asarray(controls, dtype=np.float64)
    weights = h[:, 3]
    scale = float(np.max(weights))
    if not math.isfinite(scale) or scale <= 0.0:
        raise ValueError("invalid homogeneous gauge")
    weights = weights / scale
    points = h[:, :3] / h[:, 3, None]
    weighted = poly * weights
    denom = float(np.sum(weighted))
    ddenom = dpoly.T @ weights
    if not math.isfinite(denom) or denom <= 0.0:
        raise ValueError("nonpositive rational denominator")
    basis = weighted / denom
    grad_param = (dpoly * weights[:, None] - basis[:, None] * ddenom[None, :]) / denom
    position = basis @ points
    jac = points.T @ grad_param
    det = float(np.linalg.det(jac))
    if not math.isfinite(det) or det <= 0.0:
        raise ValueError("nonpositive geometry Jacobian")
    grad_phys = grad_param @ np.linalg.inv(jac)
    return position, basis, grad_phys, det


def q1_basis(uvw: tuple[float, float, float]) -> np.ndarray:
    u, v, w = uvw
    coordinates = ((0,0,0),(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1))
    return np.asarray([
        (u if i else 1-u) * (v if j else 1-v) * (w if k else 1-w)
        for i,j,k in coordinates
    ])


def gauss01(order: int) -> list[tuple[float, float]]:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return list(zip((nodes + 1.0) * 0.5, weights * 0.5))


def constitutive(young: float, poisson: float) -> tuple[float, float, float]:
    mu = young / (2.0 * (1.0 + poisson))
    bulk = young / (3.0 * (1.0 - 2.0 * poisson))
    lam = bulk - 2.0 * mu / 3.0
    return mu, bulk, lam


def element_matrices(controls: list[list[int]], young: float, poisson: float, order: int):
    """Rebuild K_dev, K_full, Q1 mass M, and divergence coupling B locally."""
    mu, _bulk, lam = constitutive(young, poisson)
    kdev = np.zeros((81, 81))
    kfull = np.zeros((81, 81))
    mass = np.zeros((8, 8))
    coupling = np.zeros((8, 81))
    rule = gauss01(order)
    for u, wu in rule:
        for v, wv in rule:
            for w, ww in rule:
                uvw = (u, v, w)
                _x, _values, grad, det = geometry(controls, uvw)
                dv = wu * wv * ww * det
                q = q1_basis(uvw)
                mass += np.outer(q, q) * dv
                strains = np.zeros((81, 3, 3))
                divs = np.zeros(81)
                for node in range(27):
                    for component in range(3):
                        dof = 3 * node + component
                        strains[dof, component, :] += 0.5 * grad[node]
                        strains[dof, :, component] += 0.5 * grad[node]
                        divs[dof] = grad[node, component]
                traces = np.trace(strains, axis1=1, axis2=2)
                deviatoric = strains - traces[:, None, None] * np.eye(3)[None, :, :] / 3.0
                kdev += 2.0 * mu * np.einsum("iab,jab->ij", deviatoric, deviatoric) * dv
                kfull += (
                    2.0 * mu * np.einsum("iab,jab->ij", strains, strains)
                    + lam * np.outer(divs, divs)
                ) * dv
                coupling += np.outer(q, divs) * dv
    return kdev, kfull, mass, coupling


def dense_relative_error(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.linalg.norm(actual - expected) / max(np.linalg.norm(expected), 1.0))


def physical_nodes(case: dict) -> np.ndarray:
    """Recover physical Q2 geometry nodes, checking shared DOF identities."""
    points: dict[int, np.ndarray] = {}
    for block in case["blocks"]:
        controls = np.asarray(block["controls"], dtype=np.float64)
        dofs = block["dofs"]
        for local, (control, node) in enumerate(zip(controls, dofs)):
            point = control[:3] / control[3]
            prior = points.get(node)
            if prior is not None and np.linalg.norm(prior - point) > 1.0e-10:
                raise ValueError(f"incompatible physical coordinates at scalar DOF {node}")
            points[node] = point
    if len(points) != case["ndof"]:
        raise ValueError("fixture does not cover all scalar geometry DOFs")
    result = np.empty((case["ndof"], 3))
    for node, point in points.items():
        result[node] = point
    return result


def rigid_mode_count(matrix: np.ndarray) -> tuple[int, float]:
    symmetric = 0.5 * (matrix + matrix.T)
    values = eigvalsh(symmetric, check_finite=True, driver="evd")
    scale = max(float(np.max(np.abs(values))), 1.0)
    near_zero = int(np.count_nonzero(np.abs(values) <= 1.0e-9 * scale))
    return near_zero, float(values[6] / scale) if len(values) > 6 else math.nan


def clamp_nodes(points: np.ndarray) -> np.ndarray:
    """Select three non-collinear geometry nodes to eliminate rigid motion."""
    origin = points[0]
    distances = np.linalg.norm(points - origin, axis=1)
    second = int(np.argmax(distances))
    axis = points[second] - origin
    axis_norm = np.linalg.norm(axis)
    if axis_norm <= 1.0e-12:
        raise ValueError("degenerate geometry nodes")
    cross = np.cross(points - origin, axis)
    third = int(np.argmax(np.linalg.norm(cross, axis=1)))
    if third == 0 or np.linalg.norm(cross[third]) <= 1.0e-10:
        raise ValueError("no non-collinear clamp nodes")
    return np.asarray([0, second, third], dtype=int)


def run_case(case: dict) -> dict:
    ndof = int(case["ndof"])
    young = float(case["young_mpa"])
    poisson = float(case["poisson"])
    order = int(case["order"])
    if not 0.0 < young or not -1.0 < poisson < 0.5 or order < 2:
        raise ValueError(f"invalid material/quadrature in {case['name']}")
    size = 3 * ndof
    mu, bulk, lam = constitutive(young, poisson)
    kdev_rows: list[int] = []
    kdev_cols: list[int] = []
    kdev_data: list[float] = []
    full_rows: list[int] = []
    full_cols: list[int] = []
    full_data: list[float] = []
    pressure_masses = []
    local_couplings = []
    pressure_ids = []
    for block in case["blocks"]:
        dofs = np.asarray(block["dofs"], dtype=int)
        if len(dofs) != 27 or np.any(dofs < 0) or np.any(dofs >= ndof):
            raise ValueError(f"invalid Q2 DOF map in {block['id']}")
        if len(block["controls"]) != 27:
            raise ValueError(f"invalid Q2 control net in {block['id']}")
        ke_dev, ke_full, mass, coupling = element_matrices(block["controls"], young, poisson, order)
        global_dofs = (3 * dofs[:, None] + np.arange(3)[None, :]).reshape(-1)
        rr, cc = np.meshgrid(global_dofs, global_dofs, indexing="ij")
        kdev_rows.extend(rr.ravel()); kdev_cols.extend(cc.ravel()); kdev_data.extend(ke_dev.ravel())
        full_rows.extend(rr.ravel()); full_cols.extend(cc.ravel()); full_data.extend(ke_full.ravel())
        pressure_masses.append(mass)
        local_couplings.append((global_dofs, coupling))
        pressure_ids.append(block["id"])

    kdev = coo_matrix((kdev_data, (kdev_rows, kdev_cols)), shape=(size, size)).tocsr()
    kfull = coo_matrix((full_data, (full_rows, full_cols)), shape=(size, size)).tocsr()
    mass = block_diag(pressure_masses, format="csr")
    b_rows: list[int] = []
    b_cols: list[int] = []
    b_data: list[float] = []
    for element, (dofs, coupling) in enumerate(local_couplings):
        for p in range(8):
            for local, global_dof in enumerate(dofs):
                value = coupling[p, local]
                if value != 0.0:
                    b_rows.append(8 * element + p)
                    b_cols.append(int(global_dof))
                    b_data.append(float(value))
    bmat_u = coo_matrix((b_data, (b_rows, b_cols)), shape=(8 * len(pressure_ids), size)).tocsr()
    m_solve_b = solve(mass.toarray(), bmat_u.toarray(), assume_a="pos", check_finite=True)
    projected = (kdev + bulk * csr_matrix(bmat_u.T @ m_solve_b)).tocsr()
    projected = 0.5 * (projected + projected.T)
    kfull = 0.5 * (kfull + kfull.T)

    rust_projected = np.asarray(case["projected_stiffness"], dtype=float)
    rust_full = np.asarray(case["full_stiffness"], dtype=float)
    if rust_projected.shape != (size, size) or rust_full.shape != (size, size):
        raise ValueError(f"matrix dimensions mismatch in {case['name']}")
    project_error = dense_relative_error(rust_projected, projected.toarray())
    full_error = dense_relative_error(rust_full, kfull.toarray())

    # The projected volumetric term is an L2 best approximation and cannot
    # increase the volumetric energy relative to the primitive full form.
    u_test = np.asarray(case["displacement_test"], dtype=float)
    if u_test.shape != (size,) or not np.all(np.isfinite(u_test)):
        raise ValueError(f"invalid deterministic displacement in {case['name']}")
    projected_energy = float(u_test @ (projected @ u_test))
    full_energy = float(u_test @ (kfull @ u_test))
    energy_ordering = projected_energy <= full_energy + 1.0e-9 * max(abs(full_energy), 1.0)

    recovered_pressure = bulk * m_solve_b @ u_test
    rust_pressure = case["mean_stress_coefficients"]
    pressure_errors = {}
    for elem, block_id in enumerate(pressure_ids):
        values = np.asarray(rust_pressure[block_id], dtype=float)
        expected = recovered_pressure[8 * elem:8 * (elem + 1)]
        pressure_errors[block_id] = dense_relative_error(values, expected)

    # An affine displacement is represented exactly by the rational Q2 basis
    # when its coefficients are evaluated at the homogeneous control points.
    points = physical_nodes(case)
    affine = np.array([[0.013, -0.004, 0.002], [0.003, -0.006, 0.005], [-0.001, 0.004, 0.009]])
    offset = np.array([0.4, -0.2, 0.1])
    patch = (points @ affine.T + offset).reshape(-1)
    patch_projected_energy = float(patch @ (projected @ patch))
    strain = 0.5 * (affine + affine.T)
    trace = float(np.trace(strain))
    dev = strain - trace * np.eye(3) / 3.0
    volume = float(sum(
        sum(wu * wv * ww * geometry(block["controls"], (u,v,w))[3]
            for u,wu in gauss01(order)
            for v,wv in gauss01(order)
            for w,ww in gauss01(order))
        for block in case["blocks"]
    ))
    patch_exact_energy = volume * (2.0 * mu * float(np.sum(dev * dev)) + bulk * trace * trace)
    patch_error = abs(patch_projected_energy - patch_exact_energy) / max(abs(patch_exact_energy), 1.0)
    patch_pressure = bulk * m_solve_b @ patch
    patch_pressure_spread = max(
        float(np.ptp(patch_pressure[8*i:8*(i+1)])) for i in range(len(pressure_ids))
    )

    # Verify all six rigid fields and count the null eigenvalues on the small
    # fixture matrices. A reduced constrained saddle solve follows below.
    node_points = points
    center = np.mean(node_points, axis=0)
    relative = node_points - center
    rigid = []
    for axis in range(3):
        value = np.zeros((ndof, 3)); value[:,axis] = 1.0; rigid.append(value.reshape(-1))
    for axis in np.eye(3):
        value = np.cross(axis[None,:], relative)
        rigid.append(value.reshape(-1))
    rigid_residuals = [float(np.linalg.norm(projected @ mode) / max(np.linalg.norm(projected.data) * np.linalg.norm(mode), 1.0)) for mode in rigid]
    projected_zero_modes, projected_first_positive = rigid_mode_count(projected.toarray())
    full_zero_modes, full_first_positive = rigid_mode_count(kfull.toarray())

    # Manufacture a load from the deterministic field and solve the explicit
    # discontinuous-Q1 mixed system after fixing three non-collinear nodes.
    clamp = clamp_nodes(node_points)
    fixed = np.sort(np.concatenate([3 * clamp[:,None] + np.arange(3)[None,:]]).ravel())
    free = np.setdiff1d(np.arange(size), fixed, assume_unique=True)
    force = projected @ u_test
    u_fixed = u_test[fixed]
    kff = kdev[free][:,free]
    kfd = kdev[free][:,fixed]
    bf = bmat_u[:,free]
    bd = bmat_u[:,fixed]
    saddle = bmat([
        [kff, bf.T],
        [bf, -mass / bulk],
    ], format="csc")
    rhs = np.concatenate([force[free] - kfd @ u_fixed, -(bd @ u_fixed)])
    solved = spsolve(saddle, rhs)
    recovered_u = u_test.copy()
    recovered_u[free] = solved[:len(free)]
    saddle_pressure = solved[len(free):]
    mixed_u_error = dense_relative_error(recovered_u, u_test)
    mixed_p_error = dense_relative_error(saddle_pressure, recovered_pressure)
    saddle_residual = float(np.linalg.norm(saddle @ solved-rhs) / max(np.linalg.norm(rhs),1.0))

    result = {
        "name": case["name"],
        "ndof": ndof,
        "blocks": len(pressure_ids),
        "projected_matrix_relative_error": project_error,
        "full_matrix_relative_error": full_error,
        "projected_energy_n_mm": projected_energy,
        "full_energy_n_mm": full_energy,
        "projected_not_greater_than_full": bool(energy_ordering),
        "mean_stress_relative_errors": pressure_errors,
        "patch_energy_relative_error": patch_error,
        "patch_pressure_coefficient_spread_mpa": patch_pressure_spread,
        "projected_rigid_zero_modes": projected_zero_modes,
        "full_rigid_zero_modes": full_zero_modes,
        "projected_first_positive_eigenvalue_ratio": projected_first_positive,
        "full_first_positive_eigenvalue_ratio": full_first_positive,
        "rigid_mode_action_relative_residuals": rigid_residuals,
        "explicit_saddle_displacement_relative_error": mixed_u_error,
        "explicit_saddle_pressure_relative_error": mixed_p_error,
        "explicit_saddle_relative_residual": saddle_residual,
        "volume_mm3": volume,
        "patch_energy_mm3_times_mpa": patch_exact_energy,
    }
    checks = {
        "projected_matrix": project_error <= 1.0e-10,
        "full_matrix": full_error <= 1.0e-10,
        "pressure_recovery": max(pressure_errors.values(), default=0.0) <= 1.0e-9,
        "energy_ordering": energy_ordering,
        "patch_energy": patch_error <= 1.0e-9,
        "patch_constant_pressure": patch_pressure_spread <= 1.0e-5,
        "six_projected_rigid_modes": projected_zero_modes == 6 and projected_first_positive > 0.0,
        "six_full_rigid_modes": full_zero_modes == 6 and full_first_positive > 0.0,
        "rigid_mode_action": max(rigid_residuals) <= 1.0e-10,
        "explicit_saddle_displacement": mixed_u_error <= 1.0e-8,
        "explicit_saddle_pressure": mixed_p_error <= 1.0e-8,
        "explicit_saddle_residual": saddle_residual <= 1.0e-10,
    }
    result["checks"] = checks
    result["qualified"] = all(checks.values())
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixture", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.loads(args.fixture.read_text())
    if data.get("schema") != "projected_bulk_verification_v1":
        raise SystemExit("unsupported projected-bulk fixture schema")
    cases = [run_case(case) for case in data["cases"]]
    output = {"schema": "projected_bulk_reference_v1", "cases": cases,
              "qualified_cases": sum(bool(case["qualified"]) for case in cases)}
    payload = json.dumps(output, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(payload)
    else:
        print(payload, end="")
    failures = [case["name"] for case in cases if not case["qualified"]]
    if failures:
        raise SystemExit("projected-bulk checks failed: " + ", ".join(failures))


if __name__ == "__main__":
    main()
