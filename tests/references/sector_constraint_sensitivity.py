"""Independent S7 oracle for angle-dependent adaptive field constraints."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mechanical_job as fem
import sector_job_sensitivity as s6


def write(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def call(cli: Path, *args: object, success=True):
    proc = subprocess.run([str(cli), *map(str, args)], capture_output=True, text=True)
    if success and proc.returncode:
        raise RuntimeError(f"{args[0]} failed ({proc.returncode}): {proc.stdout} {proc.stderr}")
    return proc


def plan(t: Fraction) -> dict:
    n, d = t.numerator, t.denominator
    split = n > d and (n - d) ** 2 > 2 * d * d
    if not split:
        # One chart: refine a parent, then one angular/radial/axial child.
        root = "whole"
        steps = [
            {"refine": [{"root": root, "level": 0, "index": [0, 0, 0]}]},
            {"refine": [{"root": root, "level": 1, "index": [1, 0, 0]}]},
        ]
    else:
        # Two charts: first create the chart seam on both sides, then refine
        # just the first chart's seam child to expose a nonconforming raccord.
        steps = [
            {"refine": [
                {"root": "first", "level": 0, "index": [0, 0, 0]},
                {"root": "second", "level": 0, "index": [0, 0, 0]},
            ]},
            {"refine": [{"root": "first", "level": 1, "index": [1, 0, 0]}]},
        ]
    return {"schema": "adaptive_refinement_plan_v1", "max_depth": 2,
            "max_cells": 256, "steps": steps}


def job(law: dict, t: Fraction) -> dict:
    out = s6.job(law, order=16)
    out["plan"] = plan(t)
    return out


def solve(cli: Path, folder: Path, name: str, t: Fraction, j: dict):
    folder.mkdir(parents=True, exist_ok=True)
    rp, mp, jp, out = (folder / f"{name}.{ext}" for ext in
                       ("recipe.json", "model.json", "job.json", "result.json"))
    write(rp, s6.recipe(t)); write(jp, j)
    call(cli, "build", rp, mp)
    call(cli, "solve-job", mp, jp, out)
    return json.loads(mp.read_text()), json.loads(out.read_text()), mp, jp, rp


def root_bounds(block_id: str, t: Fraction) -> tuple[tuple[float, float], ...]:
    parts = block_id.split("/")
    if not parts:
        raise ValueError(f"invalid adaptive sector ID {block_id}")
    n, d = t.numerator, t.denominator
    split = n > d and (n - d) ** 2 > 2 * d * d
    base = {"whole": [(0., 1.), (0., 1.), (0., 1.)],
            "first": [(0., .5), (0., 1.), (0., 1.)],
            "second": [(.5, 1.), (0., 1.), (0., 1.)]}
    if parts[0] not in base or (split and parts[0] == "whole") or (not split and parts[0] != "whole"):
        raise ValueError(f"chart mismatch for {block_id} at t={t}")
    bounds = list(base[parts[0]])
    for path in parts[1:]:
        if len(path) != 3 or any(bit not in "01" for bit in path):
            raise ValueError(f"invalid dyadic path {block_id}")
        for axis, bit in enumerate(path):
            lo, hi = bounds[axis]
            mid = (lo + hi) * 0.5
            bounds[axis] = (lo, mid) if bit == "0" else (mid, hi)
    return tuple(bounds)


def node_keys(result: dict) -> tuple[dict[int, tuple], dict[tuple, int]]:
    occurrences: dict[int, list[tuple[str, int]]] = {}
    for block in result["analysis_blocks"]:
        for local, dof in enumerate(block["dofs"]):
            occurrences.setdefault(int(dof), []).append((block["id"], local))
    by_dof = {dof: tuple(sorted(locations)) for dof, locations in occurrences.items()}
    by_key = {key: dof for dof, key in by_dof.items()}
    return by_dof, by_key


def lift_maps(result: dict, keys: dict[int, tuple]) -> list[dict[tuple, dict]]:
    out = []
    for lift in result["component_lifts"]:
        free = [int(x) for x in lift["free_roots"]]
        rows = {}
        for line in lift["lines"]:
            row = {}
            for j, num, den in line["terms"]:
                root = free[int(j)]
                key = keys[root]
                row[key] = row.get(key, 0.0) + int(num) / int(den)
            rows[keys[int(line["dof"])]] = {"terms": row,
                                           "offset": int(line["offset"][0]) / int(line["offset"][1])}
        out.append(rows)
    return out


def rate_maps(derivative: dict, keys: dict[int, tuple]) -> list[dict[tuple, dict]]:
    result = []
    for component in derivative["component_lift_rates"]:
        rows = {}
        for line in component:
            terms = {}
            for root, rate in line["terms"]:
                key = keys[int(root)]
                terms[key] = terms.get(key, 0.0) + float(rate)
            rows[keys[int(line["dof"])]] = {"terms": terms,
                                           "offset": float(line["offset_rate"])}
        result.append(rows)
    return result


def weight_log_rates(result: dict, t: Fraction) -> tuple[dict[int, float], dict[int, tuple]]:
    by_dof, _ = node_keys(result)
    rates: dict[int, float] = {}
    for block in result["analysis_blocks"]:
        bounds = root_bounds(block["id"], t)
        a, b = bounds[0]
        for local, dof0 in enumerate(block["dofs"]):
            dof = int(dof0)
            angular = local // 9
            left, right = ((a, a), (a, b), (b, b))[angular]
            log_rate = 2.0 * float(t) * left * right / (1.0 + float(t) ** 2 * left * right)
            prior = rates.get(dof)
            if prior is not None and abs(prior - log_rate) > 2e-12:
                raise ValueError(f"inconsistent analytic weight rate at DOF {dof}: {prior}, {log_rate}")
            rates[dof] = log_rate
    if set(rates) != set(by_dof):
        raise ValueError("homogeneous weights do not cover every material DOF")
    return rates, by_dof


def lift_rate_formula(result: dict, t: Fraction, lifts: list[dict[tuple, dict]]) -> list[dict[tuple, dict]]:
    log_rates, by_dof = weight_log_rates(result, t)
    by_key = {key: dof for dof, key in by_dof.items()}
    predicted = []
    for component_rows in lifts:
        rows = {}
        for row_key, row in component_rows.items():
            wi = log_rates[by_key[row_key]]
            terms = {root_key: value * (log_rates[by_key[root_key]] - wi)
                     for root_key, value in row["terms"].items()}
            rows[row_key] = {"terms": terms, "offset": 0.0}
        predicted.append(rows)
    return predicted


def rate_fd(base: dict, plus: dict, minus: dict, t: Fraction, h: float) -> list[dict[tuple, dict]]:
    kb, _ = node_keys(base); kp, _ = node_keys(plus); km, _ = node_keys(minus)
    lb, lp, lm = lift_maps(base, kb), lift_maps(plus, kp), lift_maps(minus, km)
    out = []
    for c in range(3):
        rows = {}
        row_keys = set(lb[c]) & set(lp[c]) & set(lm[c])
        for row_key in row_keys:
            pb, pp, pm = lb[c][row_key], lp[c][row_key], lm[c][row_key]
            masters = set(pb["terms"]) | set(pp["terms"]) | set(pm["terms"])
            terms = {key: (pp["terms"].get(key, 0.) - pm["terms"].get(key, 0.)) / (2*h)
                     for key in masters}
            rows[row_key] = {"terms": terms,
                             "offset": (pp["offset"] - pm["offset"]) / (2*h)}
        out.append(rows)
    return out


def combine_rates(fine: list[dict], coarse: list[dict], fine_scale: float,
                  coarse_scale: float) -> list[dict]:
    out = []
    for fcomp, ccomp in zip(fine, coarse):
        rows = {}
        for row_key in set(fcomp) | set(ccomp):
            fr, cr = fcomp.get(row_key, {"terms": {}, "offset": 0.}), ccomp.get(row_key, {"terms": {}, "offset": 0.})
            terms = {key: fine_scale * fr["terms"].get(key, 0.) + coarse_scale * cr["terms"].get(key, 0.)
                     for key in set(fr["terms"]) | set(cr["terms"])}
            rows[row_key] = {"terms": terms,
                             "offset": fine_scale * fr["offset"] + coarse_scale * cr["offset"]}
        out.append(rows)
    return out


def error_summary(expected: list[dict], actual: list[dict]) -> dict:
    max_abs, max_fraction, count = 0.0, 0.0, 0
    worst = None
    for c, (ecomp, acomp) in enumerate(zip(expected, actual)):
        for row in set(ecomp) | set(acomp):
            erow = ecomp.get(row, {"terms": {}, "offset": 0.})
            arow = acomp.get(row, {"terms": {}, "offset": 0.})
            for master in set(erow["terms"]) | set(arow["terms"]):
                e, a = erow["terms"].get(master, 0.), arow["terms"].get(master, 0.)
                err = abs(e - a)
                frac = err / max(1e-8, 1e-5 * abs(e))
                count += 1
                max_fraction = max(max_fraction, frac)
                if err > max_abs:
                    max_abs, worst = err, (c, row, master, e, a)
            err = abs(erow["offset"] - arow["offset"])
            frac = err / max(1e-8, 1e-5 * abs(erow["offset"]))
            max_fraction = max(max_fraction, frac)
            if err > max_abs:
                max_abs, worst = err, (c, row, "offset", erow["offset"], arow["offset"])
    return {"compared_coefficients": count, "maximum_absolute_error": max_abs,
            "maximum_tolerance_fraction": max_fraction, "worst": repr(worst)}


def field_map(result: dict, key_name: str, keys: dict[int, tuple] | None = None) -> dict[tuple, np.ndarray]:
    if keys is None:
        keys, _ = node_keys(result)
    value = np.asarray(result[key_name], dtype=float).reshape(-1, 3)
    return {key: value[dof] for dof, key in keys.items()}


def fd_map(base: dict, plus: dict, minus: dict, field: str, h: float) -> dict[tuple, np.ndarray]:
    a, b, c = field_map(base, field), field_map(plus, field), field_map(minus, field)
    common = set(a) & set(b) & set(c)
    if len(common) != len(a) or len(common) != len(b) or len(common) != len(c):
        raise ValueError(f"material DOF mismatch across cold perturbations for {field}")
    return {key: (b[key] - c[key]) / (2*h) for key in common}


def field_combine(fine: dict, coarse: dict, fine_scale: float, coarse_scale: float) -> dict:
    return {key: fine_scale * fine[key] + coarse_scale * coarse[key] for key in fine}


def field_error(expected: dict, actual: dict) -> dict:
    diffs = [np.asarray(expected[k]) - np.asarray(actual[k]) for k in expected]
    norm = math.sqrt(sum(float(np.dot(v, v)) for v in diffs))
    ref = math.sqrt(sum(float(np.dot(expected[k], expected[k])) for k in expected))
    return {"relative_error": norm / max(ref, 1e-30),
            "tolerance_fraction": norm / max(1e-8, 1e-5 * ref),
            "absolute_norm": norm, "reference_norm": ref,
            "matched_material_dofs": len(diffs)}


def scalar_fd(plus: dict, minus: dict, key: str, h: float) -> float:
    return (float(plus[key]) - float(minus[key])) / (2*h)


def displacement_lift_identity(base: dict, derivative: dict, lifts: list[dict[tuple, dict]],
                               rates: list[dict[tuple, dict]]) -> dict:
    keys, by_key = node_keys(base)
    u = np.asarray(base["displacement_mm_xyz"], dtype=float).reshape(-1, 3)
    du = np.asarray(derivative["physical_displacement_rate_mm_xyz"], dtype=float).reshape(-1, 3)
    errors = []
    for c in range(3):
        for row_key, row in lifts[c].items():
            dof = by_key[row_key]
            rate = rates[c].get(row_key, {"terms": {}, "offset": 0.})
            expected = rate["offset"]
            for master_key, coefficient_rate in rate["terms"].items():
                expected += coefficient_rate * u[by_key[master_key], c]
            for master_key, coefficient in row["terms"].items():
                expected += coefficient * du[by_key[master_key], c]
            errors.append(abs(expected - du[dof, c]))
    return {"rows": len(errors), "maximum_absolute_error_mm_per_t": max(errors, default=0.),
            "mean_absolute_error_mm_per_t": float(np.mean(errors)) if errors else 0.}


def interfaces(result: dict, t: Fraction) -> list[tuple]:
    blocks = result["analysis_blocks"]
    bounds = {b["id"]: root_bounds(b["id"], t) for b in blocks}
    pairs = []
    for i, first in enumerate(blocks):
        a = bounds[first["id"]]
        for second in blocks[i+1:]:
            b = bounds[second["id"]]
            for axis in range(3):
                if abs(a[axis][1] - b[axis][0]) < 1e-13:
                    upper_a, upper_b, plane = True, False, a[axis][1]
                elif abs(b[axis][1] - a[axis][0]) < 1e-13:
                    upper_a, upper_b, plane = False, True, a[axis][0]
                else:
                    continue
                free = [j for j in range(3) if j != axis]
                overlap = [(max(a[j][0], b[j][0]), min(a[j][1], b[j][1])) for j in free]
                if all(hi - lo > 1e-13 for lo, hi in overlap):
                    pairs.append((first, second, axis, upper_a, upper_b, plane, free, overlap))
    return pairs


def trace_continuity(result: dict, t: Fraction, field_name: str,
                     field_result: dict | None = None) -> dict:
    blocks = result["analysis_blocks"]
    by_id = {b["id"]: b for b in blocks}
    keys, _ = node_keys(result)
    field_source = result if field_result is None else field_result
    field = np.asarray(field_source[field_name], dtype=float).reshape(-1, 3)
    primal_field = np.asarray(result["displacement_mm_xyz"], dtype=float).reshape(-1, 3)
    log_weights, _ = weight_log_rates(result, t)
    worst_geom, worst_field, count = 0.0, 0.0, 0
    for first, second, axis, upper_a, upper_b, plane, free, overlap in interfaces(result, t):
        for s in (0., .5, 1.):
            for q in (0., .5, 1.):
                global_param = [0., 0., 0.]
                global_param[axis] = plane
                global_param[free[0]] = overlap[0][0] + s * (overlap[0][1] - overlap[0][0])
                global_param[free[1]] = overlap[1][0] + q * (overlap[1][1] - overlap[1][0])
                points, values = [], []
                for block, upper in ((first, upper_a), (second, upper_b)):
                    bb = root_bounds(block["id"], t)
                    uvw = [(global_param[j] - bb[j][0]) / (bb[j][1] - bb[j][0]) for j in range(3)]
                    uvw[axis] = float(upper)
                    basis, x, *_ = fem.shape(block["controls"], uvw)
                    local_u = np.asarray(block["dofs"], dtype=int)
                    if field_result is None:
                        values.append(basis @ field[local_u])
                    else:
                        alpha = np.asarray([log_weights[int(d)] for d in local_u])
                        basis_rate = basis * (alpha - float(basis @ alpha))
                        values.append(basis @ field[local_u] + basis_rate @ primal_field[local_u])
                    points.append(x)
                worst_geom = max(worst_geom, float(np.linalg.norm(points[0] - points[1])))
                worst_field = max(worst_field, float(np.linalg.norm(values[0] - values[1])))
                count += 1
    return {"interface_pairs": len(interfaces(result, t)), "sampled_trace_pairs": count,
            "maximum_geometry_gap_mm": worst_geom, "maximum_field_gap": worst_field}


def interface_rate_witness(result: dict, t: Fraction,
                           rates: list[dict[tuple, dict]]) -> dict:
    keys, _ = node_keys(result)
    witnesses = []
    for first, second, axis, upper_a, upper_b, _plane, _free, _overlap in interfaces(result, t):
        local_dofs = []
        for block, upper in ((first, upper_a), (second, upper_b)):
            for local, dof in enumerate(block["dofs"]):
                ijk = (local // 9, (local // 3) % 3, local % 3)
                if ijk[axis] == 2 * int(upper):
                    local_dofs.append(int(dof))
        nonzero = []
        for dof in set(local_dofs):
            key = keys[dof]
            for component, rows in enumerate(rates):
                row = rows.get(key)
                if row and any(abs(value) > 1e-12 for value in row["terms"].values()):
                    nonzero.append((dof, component))
        if nonzero:
            witnesses.append({"first": first["id"], "second": second["id"],
                              "axis": "uvw"[axis], "nonzero_rate_rows": len(nonzero),
                              "sample": nonzero[0]})
    return {"interface_pairs": len(interfaces(result, t)),
            "interfaces_with_nonzero_constraint_rate": len(witnesses),
            "witnesses": witnesses[:12]}


def independent_stiffness(blocks: list[dict], young: float, poisson: float, order: int):
    """Vectorized independent Q2 isotropic stiffness quadrature on exported maps."""
    ndof = max(int(d) for block in blocks for d in block["dofs"]) + 1
    mu = young / (2 * (1 + poisson))
    lam = young * poisson / ((1 + poisson) * (1 - 2 * poisson))
    material = np.zeros((6, 6))
    material[:3, :3] = lam
    for i in range(3):
        material[i, i] += 2 * mu
    material[3, 3] = material[4, 4] = material[5, 5] = mu
    nodes, weights = fem.gauss01(order)
    grid = [(u, v, w, wu * wv * ww) for u, wu in zip(nodes, weights)
            for v, wv in zip(nodes, weights) for w, ww in zip(nodes, weights)]
    rows, cols, vals = [], [], []
    for block in blocks:
        strain = np.zeros((len(grid), 6, 81))
        measures = np.zeros(len(grid))
        for q, (u, v, w, quadrature_weight) in enumerate(grid):
            _basis, _position, _jac, gradient, determinant = fem.shape(block["controls"], (u, v, w))
            measures[q] = quadrature_weight * determinant
            for node in range(27):
                gx, gy, gz = gradient[node]
                col = 3 * node
                strain[q, 0, col] = gx
                strain[q, 3, col] = gy
                strain[q, 5, col] = gz
                strain[q, 3, col + 1] = gx
                strain[q, 1, col + 1] = gy
                strain[q, 4, col + 1] = gz
                strain[q, 5, col + 2] = gx
                strain[q, 4, col + 2] = gy
                strain[q, 2, col + 2] = gz
        material_strain = np.einsum("ab,qbj->qaj", material, strain, optimize=True)
        local = np.einsum("qai,qaj,q->ij", strain, material_strain, measures, optimize=True)
        gdof = (3 * np.asarray(block["dofs"], dtype=int)[:, None]
                + np.arange(3)[None, :]).ravel()
        rr, cc = np.meshgrid(gdof, gdof, indexing="ij")
        rows.extend(rr.ravel()); cols.extend(cc.ravel()); vals.extend(local.ravel())
    return coo_matrix((vals, (rows, cols)), shape=(3 * ndof, 3 * ndof)).tocsr(), ndof


def cprime_ablation(model: dict, base: dict, j: dict, derivative: dict,
                    rates: list[dict[tuple, dict]]) -> dict:
    blocks = base["analysis_blocks"]
    stiffnesses = {}
    for order in (8, 12):
        stiffnesses[order], ndof = independent_stiffness(
            blocks, j["material"]["young_mpa"], j["material"]["poisson"], order)
    force, _resultant, _moment, _area, _loads, _scale, _radius = fem.independent_load(
        model, base, j, j["order"], ndof)
    u = np.asarray(base["displacement_mm_xyz"], dtype=float)
    keys, by_key = node_keys(base)
    cp_x = np.zeros_like(u).reshape(-1, 3)
    lifts = lift_maps(base, keys)
    for c in range(3):
        for row_key, row in rates[c].items():
            dof = by_key[row_key]
            value = 0.0
            for root_key, rate in row["terms"].items():
                value += rate * u[3 * by_key[root_key] + c]
            cp_x[dof, c] = value
    corrections = {}
    residual_norms = {}
    for order, k in stiffnesses.items():
        residual = force - k @ u
        corrections[order] = 2.0 * float(np.dot(residual, cp_x.ravel()))
        residual_norms[order] = float(np.linalg.norm(residual))
    difference = stiffnesses[12] - stiffnesses[8]
    stiffness_difference = float(np.sqrt(difference.multiply(difference).sum()) /
                                  max(1e-30, np.sqrt(stiffnesses[12].multiply(stiffnesses[12]).sum())))
    tangent = float(derivative["compliance_rate_tangent_n_mm"])
    correction = corrections[12]
    return {"independent_stiffness_shape": list(stiffnesses[12].shape),
            "independent_stiffness_orders": [8, 12],
            "stiffness_relative_difference_q8_q12": stiffness_difference,
            "independent_force_norm": float(np.linalg.norm(force)),
            "native_force_difference_norm": float(np.linalg.norm(force - np.asarray(base["force_n_xyz"]))),
            "residual_norm_by_order": {str(k): v for k, v in residual_norms.items()},
            "cprime_x_norm": float(np.linalg.norm(cp_x)),
            "2_residual_dot_cprime_x_n_mm": correction,
            "correction_by_order_n_mm": {str(k): v for k, v in corrections.items()},
            "correction_q8_q12_relative_difference": abs(corrections[12]-corrections[8]) /
                max(abs(corrections[12]), 1e-30),
            "compliance_rate_with_cprime_n_mm": tangent,
            "ablated_rate_without_cprime_n_mm": tangent - correction,
            "relative_ablation_effect": abs(correction) / max(abs(tangent), 1e-30)}


def case_specs():
    return [
        (Fraction(1, 2), {"select": {"feature": "inner"}, "pressure_mpa": 0.3}, "pressure_inner"),
        (Fraction(3, 2), {"select": {"feature": "end"}, "vector_mpa": [0.13, -0.2, 0.09]}, "vector_end"),
        (Fraction(3, 1), {"select": {"feature": "angle_end"},
                          "stress_mpa": [[0.3, 0.15, 0.0], [0.15, -0.2, 0.3], [0.0, 0.3, 0.4]]}, "stress_angle_end"),
    ]


def run_case(cli: Path, work: Path, t: Fraction, law: dict, name: str) -> dict:
    folder = work / name
    j = job(law, t)
    base_model, base, mp, jp, rp = solve(cli, folder, "base", t, j)
    direction = {"schema": "solidcraft_sector_angle_velocity_v2"}
    dp, derivative_path = folder / "direction.json", folder / "derivative.json"
    write(dp, direction)
    call(cli, "differentiate-job", mp, jp, dp, derivative_path)
    derivative = json.loads(derivative_path.read_text())
    if derivative.get("schema") != "solidcraft_sector_angle_derivative_v1":
        raise ValueError(f"unexpected S7 derivative schema: {derivative.get('schema')}")

    base_keys, _ = node_keys(base)
    base_lifts = lift_maps(base, base_keys)
    native_rates = rate_maps(derivative, base_keys)
    analytic_rates = lift_rate_formula(base, t, base_lifts)
    formula_error = error_summary(analytic_rates, native_rates)

    differences = {"force_n_xyz": [], "displacement_mm_xyz": [], "external_work_n_mm": [],
                   "lift_rates": []}
    fd_records = []
    for denominator in (16, 32, 64):
        h = Fraction(1, denominator)
        results = []
        for sign, side in ((1, "plus"), (-1, "minus")):
            perturb = t + sign * h
            sub = folder / f"{side}_{denominator}"
            _mod, result, _mp, _jp, _rp = solve(cli, sub, "solve", perturb, j)
            results.append(result)
        plus, minus = results
        force_fd = fd_map(base, plus, minus, "force_n_xyz", float(h))
        disp_fd = fd_map(base, plus, minus, "displacement_mm_xyz", float(h))
        lift_fd = rate_fd(base, plus, minus, t, float(h))
        differences["force_n_xyz"].append(force_fd)
        differences["displacement_mm_xyz"].append(disp_fd)
        differences["external_work_n_mm"].append(scalar_fd(plus, minus, "external_work_n_mm", float(h)))
        differences["lift_rates"].append(lift_fd)
        fd_records.append({"step": str(h), "step_float": float(h),
                           "plus_blocks": len(plus["analysis_blocks"]),
                           "minus_blocks": len(minus["analysis_blocks"]),
                           "plus_scalar_dofs": plus["scalar_dofs"],
                           "minus_scalar_dofs": minus["scalar_dofs"]})

    # Sixth-order centered estimates from the three frozen fixed-regime steps.
    force_r4_2 = field_combine(differences["force_n_xyz"][1], differences["force_n_xyz"][0], 4/3, -1/3)
    force_r4_1 = field_combine(differences["force_n_xyz"][2], differences["force_n_xyz"][1], 4/3, -1/3)
    force_r6 = field_combine(force_r4_1, force_r4_2, 16/15, -1/15)
    disp_r4_2 = field_combine(differences["displacement_mm_xyz"][1], differences["displacement_mm_xyz"][0], 4/3, -1/3)
    disp_r4_1 = field_combine(differences["displacement_mm_xyz"][2], differences["displacement_mm_xyz"][1], 4/3, -1/3)
    disp_r6 = field_combine(disp_r4_1, disp_r4_2, 16/15, -1/15)
    force_report = field_map(derivative, "nodal_force_rate_n_xyz", base_keys)
    disp_report = field_map(derivative, "physical_displacement_rate_mm_xyz", base_keys)
    force_error = field_error(force_r6, force_report)
    disp_error = field_error(disp_r6, disp_report)
    work_r4_2 = (4*differences["external_work_n_mm"][1] - differences["external_work_n_mm"][0]) / 3
    work_r4_1 = (4*differences["external_work_n_mm"][2] - differences["external_work_n_mm"][1]) / 3
    work_r6 = (16*work_r4_1 - work_r4_2) / 15
    tangent = float(derivative["compliance_rate_tangent_n_mm"])
    adjoint = float(derivative["compliance_rate_adjoint_n_mm"])
    work_error = abs(work_r6 - tangent)
    volume_expected = 7 * (8**2 - 3**2) / (1 + float(t)**2)
    exact_resultant = s6.exact_rate(float(t), law)
    resultant = np.asarray(derivative["applied_resultant_rate_n"], dtype=float)

    lift_fd_r4_2 = combine_rates(differences["lift_rates"][1], differences["lift_rates"][0], 4/3, -1/3)
    lift_fd_r4_1 = combine_rates(differences["lift_rates"][2], differences["lift_rates"][1], 4/3, -1/3)
    lift_fd_r6 = combine_rates(lift_fd_r4_1, lift_fd_r4_2, 16/15, -1/15)
    lift_fd_error = error_summary(lift_fd_r6, native_rates)
    lift_formula_error = formula_error
    continuity = trace_continuity(base, t, "displacement_mm_xyz")
    tangent_continuity = trace_continuity(base, t, "physical_displacement_rate_mm_xyz", derivative)
    identity = displacement_lift_identity(base, derivative, base_lifts, native_rates)
    rate_witness = interface_rate_witness(base, t, native_rates)
    weight_rates, _ = weight_log_rates(base, t)
    ablation = (cprime_ablation(base_model, base, j, derivative, native_rates)
                if name == "pressure_inner" else {"not_run": "one independent case is sufficient for C-prime ablation"})

    checks = {
        "volume_relative_error": abs(float(derivative["volume_rate_mm3"]) / volume_expected - 1),
        "resultant_relative_error": float(np.linalg.norm(resultant - exact_resultant) / max(np.linalg.norm(exact_resultant), 1e-30)),
        "tangent_adjoint_absolute_error": abs(tangent - adjoint),
        "tangent_adjoint_tolerance_fraction": abs(tangent-adjoint) / max(1e-9, 1e-8*abs(tangent)),
        "force_fd_relative_error": force_error["relative_error"],
        "force_fd_tolerance_fraction": force_error["tolerance_fraction"],
        "displacement_fd_relative_error": disp_error["relative_error"],
        "displacement_fd_tolerance_fraction": disp_error["tolerance_fraction"],
        "work_fd_relative_error": abs(work_r6 - tangent) / max(1e-8, 1e-5*abs(tangent)),
        "lift_fd_max_tolerance_fraction": lift_fd_error["maximum_tolerance_fraction"],
        "lift_formula_max_tolerance_fraction": lift_formula_error["maximum_tolerance_fraction"],
        "lift_identity_max_abs_mm_per_t": identity["maximum_absolute_error_mm_per_t"],
        "primal_interface_gap_mm": continuity["maximum_field_gap"],
        "tangent_interface_gap_mm_per_t": tangent_continuity["maximum_field_gap"],
        "nonzero_constraint_rate_interface_count": rate_witness["interfaces_with_nonzero_constraint_rate"],
    }
    passed = (checks["volume_relative_error"] <= 1e-8
              and checks["resultant_relative_error"] <= 1e-8
              and checks["tangent_adjoint_tolerance_fraction"] <= 1
              and checks["force_fd_tolerance_fraction"] <= 1
              and checks["displacement_fd_tolerance_fraction"] <= 1
              and checks["work_fd_relative_error"] <= 1
              and checks["lift_fd_max_tolerance_fraction"] <= 1
              and checks["lift_formula_max_tolerance_fraction"] <= 1
              and checks["lift_identity_max_abs_mm_per_t"] <= 1e-8
              and continuity["maximum_geometry_gap_mm"] <= 1e-8
              and continuity["maximum_field_gap"] <= 1e-8
              and tangent_continuity["maximum_geometry_gap_mm"] <= 1e-8
              and tangent_continuity["maximum_field_gap"] <= 1e-8
              and rate_witness["interfaces_with_nonzero_constraint_rate"] > 0)
    return {"id": name, "t": [t.numerator, t.denominator], "recipe": s6.recipe(t), "job": j,
            "direction": direction, "chart_count": base["source_audit"].get("blocks"),
            "adaptive_blocks": len(base["analysis_blocks"]), "scalar_dofs": base["scalar_dofs"],
            "relative_residual": base["relative_residual"],
            "native_derivative": derivative,
            "reported_lift_rate_nonzero_lines": sum(len(rows) for rows in derivative["component_lift_rates"]),
            "lift_formula_check": lift_formula_error, "lift_fd_check": lift_fd_error,
            "lift_identity": identity, "cold_finite_differences": fd_records,
            "weight_log_rate_range": [min(weight_rates.values()), max(weight_rates.values())],
            "interface_constraint_rate_witness": rate_witness,
            "sixth_order_force_fd": force_error, "sixth_order_displacement_fd": disp_error,
            "sixth_order_work_fd_n_mm": work_r6,
            "analytic_resultant_rate_n": exact_resultant.tolist(),
            "analytic_volume_rate_mm3": volume_expected,
            "primal_interface_continuity": continuity,
            "tangent_interface_continuity": tangent_continuity,
            "cprime_ablation": ablation, "checks": checks, "passed": bool(passed)}


def refusals(cli: Path, work: Path) -> list[dict]:
    base_t = Fraction(1, 2)
    cases = []
    spec = {"schema": "solidcraft_sector_angle_velocity_v2"}
    rec = s6.recipe(base_t)
    j = job({"select": {"feature": "end"}, "vector_mpa": [0.13, -0.2, 0.09]}, base_t)
    bads = [
        ("incompatible_family", {"family": "reference_rectangle_v1"}, j, spec),
        ("spatial_plane_support", rec, dict(j, supports=[dict(select={"plane": {"axis": 2, "coordinate_mm": 0}}, components=[0, 1, 2])]), spec),
        ("nonhomogeneous_support", rec, dict(j, supports=[dict(select={"feature": "angle_start"}, components=[0, 1, 2],
                          affine_mm={"denominator": 100, "numerators": [[1,0,0,0],[0,0,0,0],[0,0,0,0]]})]), spec),
    ]
    for name, recipe_value, jj, direction in bads:
        folder = work / "refusals" / name; folder.mkdir(parents=True, exist_ok=True)
        rp, mp, jp, dp, out = [folder / x for x in
            ("recipe.json", "model.json", "job.json", "direction.json", "preserved.json")]
        write(rp, recipe_value); write(jp, jj); write(dp, direction)
        call(cli, "build", rp, mp)
        sentinel = b"prior accepted derivative target\n"; out.write_bytes(sentinel)
        proc = call(cli, "differentiate-job", mp, jp, dp, out, success=False)
        cases.append({"name": name, "refused": proc.returncode != 0,
                      "preserved": out.read_bytes() == sentinel,
                      "diagnostic": (proc.stdout + proc.stderr).strip()})
    return cases


def run(cli: Path, work: Path) -> dict:
    started = time.monotonic(); work.mkdir(parents=True, exist_ok=True)
    sha = hashlib.sha256(cli.read_bytes()).hexdigest()
    rows = [run_case(cli, work, t, law, name) for t, law, name in case_specs()]
    refusal_rows = refusals(cli, work)
    nonzero_ablation = [x["cprime_ablation"]["relative_ablation_effect"] for x in rows
                        if "relative_ablation_effect" in x["cprime_ablation"]]
    return {"schema": "sector_constraint_sensitivity_oracle_v1",
            "provenance": {"cli": str(cli), "cli_sha256": sha,
                           "checkout_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE.parents[1], text=True).strip(),
                           "python": sys.version.split()[0], "elapsed_s": time.monotonic()-started},
            "protocol": {"angles_t": [[1, 2], [3, 2], [3, 1]],
                         "finite_difference_steps": ["1/16", "1/32", "1/64"],
                         "quadrature_order": 16, "adaptive_depth": 2,
                         "lift_relative_tolerance": 1e-5, "lift_absolute_tolerance": 1e-8,
                         "same_chart_regime_required": True,
                         "ablation": "independent mechanical_job.assemble and independent_load; 2(f-Ku)^T C'x"},
            "cases": rows, "refusals": refusal_rows,
            "maximum_ablation_relative_effect": max(nonzero_ablation, default=0.),
            "qualified": all(x["passed"] for x in rows)
                         and all(x["refused"] and x["preserved"] for x in refusal_rows)
                         and max(nonzero_ablation, default=0.) > 1e-10}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.cli.resolve(), args.work.resolve())
    write(args.output, report)
    print("S7 sector constraint sensitivity", report["qualified"])
    raise SystemExit(int(not report["qualified"]))
