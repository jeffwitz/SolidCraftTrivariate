#!/usr/bin/env python3
"""Independent SciPy oracle and campaign runner for versioned mechanical jobs.

For each frozen case this driver asks the CLI to build its recipe and solve its
job, then rebuilds the Q2 stiffness, surface traction, and constrained solution
from the serialized effective mesh without calling Rust analysis routines.
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
import scipy
from scipy.sparse import coo_matrix, eye, kron
from scipy.sparse.linalg import spsolve


def basis1(t: float):
    s = 1.0 - t
    return np.array([s * s, 2 * s * t, t * t]), np.array([-2 * s, 2 - 4 * t, 2 * t])


def shape(controls, uvw):
    b = [basis1(float(t)) for t in uvw]
    p = np.empty(27)
    dp = np.empty((27, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                n = (i * 3 + j) * 3 + k
                p[n] = b[0][0][i] * b[1][0][j] * b[2][0][k]
                dp[n] = (b[0][1][i] * b[1][0][j] * b[2][0][k],
                         b[0][0][i] * b[1][1][j] * b[2][0][k],
                         b[0][0][i] * b[1][0][j] * b[2][1][k])
    h = np.asarray(controls, dtype=float)
    w = h[:, 3]
    points = h[:, :3] / w[:, None]
    den = p @ w
    weighted = p * w / den
    dden = dp.T @ w
    dweighted = (dp * w[:, None] - weighted[:, None] * dden[None, :]) / den
    x = weighted @ points
    jac = points.T @ dweighted
    det = float(np.linalg.det(jac))
    if not math.isfinite(det) or det <= 0:
        raise ValueError("nonpositive geometry Jacobian in independent oracle")
    grad = dweighted @ np.linalg.inv(jac)
    return weighted, x, jac, grad, det


def gauss01(order):
    x, w = np.polynomial.legendre.leggauss(order)
    return (x + 1) * 0.5, w * 0.5


def assemble(blocks, young, poisson, order):
    ndof = max(int(d) for block in blocks for d in block["dofs"]) + 1
    rows, cols, vals = [], [], []
    u, wu = gauss01(order)
    mu = young / (2 * (1 + poisson))
    lam = young * poisson / ((1 + poisson) * (1 - 2 * poisson))
    for block in blocks:
        ids = np.asarray(block["dofs"], dtype=int)
        local = np.zeros((81, 81))
        for a, wa in zip(u, wu):
            for b, wb in zip(u, wu):
                for c, wc in zip(u, wu):
                    _n, _x, _j, grad, det = shape(block["controls"], (a, b, c))
                    eps = np.zeros((81, 3, 3))
                    div = np.zeros(81)
                    for i in range(27):
                        for q in range(3):
                            r = 3 * i + q
                            eps[r, q, :] += 0.5 * grad[i]
                            eps[r, :, q] += 0.5 * grad[i]
                            div[r] = grad[i, q]
                    local += (2 * mu * np.einsum("iab,jab->ij", eps, eps)
                              + lam * np.outer(div, div)) * (wa * wb * wc * det)
        gdof = (3 * ids[:, None] + np.arange(3)[None, :]).ravel()
        rr, cc = np.meshgrid(gdof, gdof, indexing="ij")
        rows.extend(rr.ravel()); cols.extend(cc.ravel()); vals.extend(local.ravel())
    return coo_matrix((vals, (rows, cols)), shape=(3 * ndof, 3 * ndof)).tocsr(), ndof


def parse_side(face):
    name = face["side"]
    return "uvw".index(name[0]), int(name[1])


def expected_root_faces(model, selector):
    """Independently apply a named feature and exact homogeneous plane filter."""
    blocks = {b["id"]: b for b in model["blocks"]}
    result = set()
    for face in model["design_faces"]:
        if selector["feature"] not in ("*", face["feature"]):
            continue
        if "plane" in selector:
            plane = selector["plane"]
            axis, coordinate = int(plane["axis"]), int(plane["coordinate_mm"])
            block = blocks[face["block"]]
            ax, upper = parse_side(face)
            dims = [int(d) + 1 for d in block["degrees"]]
            points = block["controls"]
            matches = True
            for i in range(dims[0]):
                for j in range(dims[1]):
                    for k in range(dims[2]):
                        index = [i, j, k]
                        if index[ax] != block["degrees"][ax] * upper:
                            continue
                        h = points[(i * dims[1] + j) * dims[2] + k]
                        if int(h[axis]) != int(h[3]) * coordinate:
                            matches = False
            if not matches:
                continue
        result.add((face["block"], face["side"]))
    return result


def face_rule(block, side, order):
    axis, upper = parse_side(side)
    axes = [i for i in range(3) if i != axis]
    nodes, weights = gauss01(order)
    samples = []
    for i, a in enumerate(nodes):
        for j, b in enumerate(nodes):
            p = [0.0, 0.0, 0.0]
            p[axis] = float(upper)
            p[axes[0]], p[axes[1]] = float(a), float(b)
            _n, x, jac, grad, _det = shape(block["controls"], p)
            normal = np.cross(jac[:, axes[0]], jac[:, axes[1]])
            # Outward signs for the two increasing free parameter axes.
            sign = {(0, 0): -1, (0, 1): 1, (1, 0): 1,
                    (1, 1): -1, (2, 0): -1, (2, 1): 1}[(axis, upper)]
            normal *= sign
            area_jac = float(np.linalg.norm(normal))
            samples.append((p, x, grad, normal / area_jac, area_jac * weights[i] * weights[j]))
    return samples


def law_value(spec, normal):
    if "vector_mpa" in spec:
        return np.asarray(spec["vector_mpa"], dtype=float)
    if "pressure_mpa" in spec:
        return -float(spec["pressure_mpa"]) * normal
    return np.asarray(spec["stress_mpa"], dtype=float) @ normal


def independent_load(model, result, job, order, ndof):
    block_by_id = {b["id"]: b for b in result["analysis_blocks"]}
    face_features = {(f["block"], f["side"]): f["feature"] for f in model["design_faces"]}
    f = np.zeros(3 * ndof)
    total = np.zeros(3)
    moment = np.zeros(3)
    area_total = 0.0
    absolute_force_scale = 0.0
    max_radius = 0.0
    load_reports = []
    for load, report in zip(job["loads"], result["loads"]):
        expected_faces = {(x[0], x[1]) for x in report["faces"]}
        expected_roots = expected_root_faces(model, load["select"])
        actual_roots = {(block.split("/")[0], side) for block, side in expected_faces}
        if actual_roots != expected_roots:
            raise ValueError(f"independent whole-face load selection mismatch: {actual_roots ^ expected_roots}")
        selected = set()
        selector = load["select"]
        for key in expected_faces:
            feature = face_features.get((key[0].split("/")[0], key[1]))
            if feature is None:
                raise ValueError(f"analysis face lost source provenance: {key}")
            if selector["feature"] != "*" and selector["feature"] != feature:
                continue
            if "plane" in selector:
                plane = selector["plane"]
                axis = int(plane["axis"])
                coordinate = int(plane["coordinate_mm"])
                block = block_by_id.get(key[0])
                if block is None:
                    continue
                ax, upper = parse_side({"side": key[1]})
                free = [i for i in range(3) if i != ax]
                on_plane = True
                for i in range(27):
                    ijk = (i // 9, (i // 3) % 3, i % 3)
                    if ijk[ax] != 2 * upper:
                        continue
                    h = block["controls"][i]
                    if h[axis] != h[3] * coordinate:
                        on_plane = False; break
                if not on_plane:
                    continue
            selected.add(key)
        if selected != expected_faces:
            raise ValueError(f"independent selected-face mismatch: {selected ^ expected_faces}")
        area = 0.0
        case_total = np.zeros(3)
        case_moment = np.zeros(3)
        for bid, side in sorted(expected_faces):
            block = block_by_id[bid]
            ids = np.asarray(block["dofs"], dtype=int)
            for uvw, x, _grad, normal, measure in face_rule(block, {"side": side}, order):
                traction = law_value(load, normal)
                area += measure
                absolute_force_scale += float(np.linalg.norm(traction)) * measure
                max_radius = max(max_radius, float(np.linalg.norm(x)))
                case_total += traction * measure
                case_moment += np.cross(x, traction) * measure
                values = shape(block["controls"], uvw)[0]
                for i, node in enumerate(ids):
                    f[3 * node:3 * node + 3] += values[i] * traction * measure
        total += case_total; moment += case_moment; area_total += area
        load_reports.append({"area_mm2": area, "resultant_n": case_total.tolist(),
                             "moment_n_mm": case_moment.tolist()})
    return f, total, moment, area_total, load_reports, absolute_force_scale, max_radius


def projection_from_lifts(result, ndof):
    columns = []
    for lift in result["component_lifts"]:
        rr, cc, vv = [], [], []
        for line in lift["lines"]:
            d = int(line["dof"])
            for root, num, den in line["terms"]:
                rr.append(d); cc.append(int(root)); vv.append(int(num) / int(den))
            if int(line["offset"][0]) != 0:
                raise ValueError("M5 homogeneous lift unexpectedly has a nonzero offset")
        columns.append(coo_matrix((vv, (rr, cc)), shape=(ndof, len(lift["free_roots"]))).tocsr())
    return columns


def independently_solve(k, f, result, ndof):
    pcomp = projection_from_lifts(result, ndof)
    # Component-major reduced columns, scalar-major physical vector rows.
    from scipy.sparse import hstack
    p = hstack([kron(pcomp[c], np.eye(3)[:, c].reshape(3, 1), format="csr") for c in range(3)], format="csr")
    reduced = (p.T @ k @ p).tocsr()
    rhs = p.T @ f
    x = spsolve(reduced, rhs)
    return np.asarray(p @ x).ravel(), reduced


def run_one(binary: Path, fixture: Path, outdir: Path):
    data = json.loads(fixture.read_text())
    if "expected_error" in data:
        with tempfile.TemporaryDirectory() as temp:
            model_path = Path(temp) / "model.json"
            recipe_path = Path(temp) / "recipe.json"
            recipe_path.write_text(json.dumps(data["recipe"]))
            built = subprocess.run([str(binary), "build", str(recipe_path), str(model_path)], capture_output=True, text=True)
            if built.returncode != 0:
                raise RuntimeError(f"model build failed: {built.stdout} {built.stderr}")
            job_path = Path(temp) / "job.json"; job_path.write_text(json.dumps(data["job"]))
            refused_path = Path(temp) / "refused.json"
            refused_path.write_text("previous-good-output\n")
            out = subprocess.run([str(binary), "solve-job", str(model_path), str(job_path), str(refused_path)], capture_output=True, text=True)
            report = json.loads(out.stdout)
            if out.returncode == 0 or report.get("reason") != data["expected_error"]:
                raise AssertionError(f"refusal mismatch {data['id']}: {report}")
            preserved = refused_path.read_text() == "previous-good-output\n"
            if not preserved:
                raise AssertionError("refusal mutated the existing output")
            return {"id": data["id"], "status": "refused", "reason": report["reason"], "output_unchanged": preserved}
    outdir.mkdir(parents=True, exist_ok=True)
    model_path = outdir / f"{data['id']}.model.json"
    job_path = outdir / f"{data['id']}.job.json"
    result_path = outdir / f"{data['id']}.result.json"
    recipe_path = outdir / f"{data['id']}.recipe.json"
    recipe_path.write_text(json.dumps(data["recipe"]))
    job_path.write_text(json.dumps(data["job"], indent=2) + "\n")
    built = subprocess.run([str(binary), "build", str(recipe_path), str(model_path)], capture_output=True, text=True)
    if built.returncode != 0:
        raise RuntimeError(f"build failed for {data['id']}: {built.stdout} {built.stderr}")
    solved = subprocess.run([str(binary), "solve-job", str(model_path), str(job_path), str(result_path)], capture_output=True, text=True)
    response = json.loads(solved.stdout)
    if solved.returncode != 0:
        raise RuntimeError(f"solve failed for {data['id']}: {response}")
    model = json.loads(model_path.read_text()); result = json.loads(result_path.read_text()); job = data["job"]
    for support, support_result in zip(job["supports"], result["supports"]):
        actual = {(face[0].split("/")[0], face[1]) for face in support_result["faces"]}
        expected = expected_root_faces(model, support["select"])
        if actual != expected:
            raise ValueError(f"independent support selection mismatch: {actual ^ expected}")
    k, ndof = assemble(result["analysis_blocks"], job["material"]["young_mpa"], job["material"]["poisson"], job["order"])
    force, resultant, moment, area, load_reports, load_scale, max_radius = independent_load(model, result, job, job["order"], ndof)
    rust_force = np.asarray(result["force_n_xyz"], dtype=float)
    rhs_error = float(np.max(np.abs(force - rust_force)))
    u_rust = np.asarray(result["displacement_mm_xyz"], dtype=float)
    u_scipy, reduced = independently_solve(k, force, result, ndof)
    displacement_error = float(np.linalg.norm(u_rust - u_scipy) / max(np.linalg.norm(u_scipy), 1.0e-30))
    coefficient_max_error = float(np.max(np.abs(u_rust - u_scipy)))
    energy = float(u_rust @ (k @ u_rust) * 0.5)
    external_work = float(force @ u_rust)
    reaction = np.asarray(result["support_reaction_n"], dtype=float)
    force_balance = float(np.linalg.norm(reaction + resultant) / max(np.linalg.norm(resultant), 1.0))
    moment_balance = float(np.linalg.norm(np.asarray(result["support_moment_n_mm"]) + moment) / max(load_scale * max_radius, 1.0))
    scipy_residual = k @ u_rust - force
    points = np.full((ndof, 3), np.nan)
    for block in result["analysis_blocks"]:
        for h, node in zip(block["controls"], block["dofs"]):
            points[node] = np.asarray(h[:3], dtype=float) / h[3]
    reaction_oracle = scipy_residual.reshape(-1, 3).sum(axis=0)
    reaction_moment_oracle = np.cross(points, scipy_residual.reshape(-1, 3)).sum(axis=0)
    reaction_error = float(np.linalg.norm(np.asarray(result["support_reaction_n"]) - reaction_oracle) / max(load_scale, 1.0))
    reaction_moment_error = float(np.linalg.norm(np.asarray(result["support_moment_n_mm"]) - reaction_moment_oracle) / max(load_scale * max_radius, 1.0))
    k_symmetry = float(np.linalg.norm((k - k.T).data) / max(np.linalg.norm(k.data), 1.0))
    affine_error = None
    if data["oracle"]["kind"] == "affine_displacement":
        grad = np.asarray(data["oracle"]["gradient"], dtype=float)
        offset = np.asarray(data["oracle"]["offset"], dtype=float)
        points = np.full((ndof, 3), np.nan)
        for block in result["analysis_blocks"]:
            for h, node in zip(block["controls"], block["dofs"]):
                x = np.asarray(h[:3], dtype=float) / h[3]
                if np.all(np.isfinite(points[node])) and np.linalg.norm(points[node] - x) > 1e-9:
                    raise ValueError("shared analysis DOF has inconsistent physical coordinates")
                points[node] = x
        expected = points @ grad.T + offset
        affine_error = float(np.max(np.abs(u_rust.reshape(-1, 3) - expected)))
    report = {"id": data["id"], "status": "qualified", "ndof": ndof,
              "effective_blocks": len(result["analysis_blocks"]), "load_area_mm2": area,
              "reintegrated_resultant_n": resultant.tolist(), "reintegrated_moment_n_mm": moment.tolist(),
              "load_rhs_max_abs_error_n": rhs_error, "load_groups": load_reports,
              "scipy_displacement_relative_error": displacement_error,
              "scipy_displacement_max_abs_error_mm": coefficient_max_error,
              "displacement_coefficient_layout": "scalar_dof_major_xyz_interleaved",
              "displacement_coefficient_count": int(u_rust.size),
              "affine_displacement_max_error_mm": affine_error,
              "energy_from_scipy_stiffness_n_mm": energy,
              "external_work_n_mm": external_work,
              "energy_work_relative_difference": abs(2 * energy - external_work) / max(abs(external_work), 1.0),
              "support_reaction_n": reaction.tolist(), "force_balance_relative": force_balance,
              "moment_balance_relative": moment_balance, "scipy_reaction_relative_error": reaction_error,
              "scipy_reaction_moment_relative_error": reaction_moment_error,
              "stiffness_symmetry_relative": k_symmetry,
              "rust_relative_residual": result["relative_residual"], "quadrature_order": job["order"],
              "matrix_order": int(k.shape[0]), "reduced_order": int(reduced.shape[0])}
    if rhs_error > 1e-8 or displacement_error > 1e-8 or force_balance > 1e-8 or moment_balance > 1e-8 or reaction_error > 1e-8 or reaction_moment_error > 1e-8 or report["energy_work_relative_difference"] > 1e-8 or k_symmetry > 1e-12 or (affine_error is not None and affine_error > 1e-8):
        raise AssertionError(f"M5 independent check failed for {data['id']}: {report}")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--jobs", type=Path, default=Path(__file__).parent / "mechanical_jobs")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, default=Path("build/mechanical-jobs"))
    args = parser.parse_args()
    start = time.monotonic()
    results = [run_one(args.binary.resolve(), p, args.build_dir) for p in sorted(args.jobs.glob("*.json"))]
    expected = 6
    if len(results) != expected or sum(x["status"] == "qualified" for x in results) != 5:
        raise AssertionError(f"expected five admitted jobs and one refusal, received {len(results)} records")
    doc = {"schema": "solidcraft_mechanical_job_oracle_v1", "oracle": "independent SciPy Q2 assembly, boundary quadrature, and constrained sparse solve",
           "environment": {"python": __import__("platform").python_version(), "numpy": np.__version__, "scipy": scipy.__version__}, "cases": results,
           "elapsed_seconds": time.monotonic() - start,
           "thresholds": {"rhs_max_abs_error_n": 1e-8, "displacement_relative_error": 1e-8,
                          "force_and_moment_balance_relative": 1e-8, "energy_work_relative": 1e-8,
                          "affine_displacement_max_error_mm": 1e-9}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(doc, indent=2, allow_nan=False) + "\n")
    print(json.dumps(doc, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
