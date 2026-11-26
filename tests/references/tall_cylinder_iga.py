#!/usr/bin/env python3
"""Independent scalar H1 solve on the tall-cylinder Bezier/NURBS covers.

The geometry controls and knot vectors are read as input data. SciPy evaluates
the rational tensor-product bases; SciPy sparse matrices assemble and solve a
manufactured Poisson problem. Coincident control points identify C0 traces.
The closed annular patch can be run with its periodic self-seam either tied or
deliberately left open.
"""

from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_patch(path: Path):
    raw = json.loads(path.read_text())
    if isinstance(raw, dict):
        seam_patch_ids = {
            end["patch_id"]
            for seam in raw.get("self_seams", [])
            for end in (seam["first"], seam["second"])
        }
        raw = raw["patches"]
    else:
        seam_patch_ids = set()
    patches = []
    for item in raw:
        shape = tuple(item["control_shape"])
        controls = item.get("controls_homogeneous", item.get("controls"))
        ctrl = np.asarray(controls, dtype=float).reshape((*shape, 4))
        closed = (
            bool(item.get("closed_u", False))
            or bool(item.get("self_seams"))
            or item.get("id") in seam_patch_ids
        )
        patches.append(
            {
                "degrees": item["degrees"],
                "knots": item["knots"],
                "ctrl": ctrl,
                "closed_u": closed,
                "sectors": item.get(
                    "sectors", item.get("cells", [item.get("id", "patch")])
                ),
            }
        )
    return patches


def load_bezier_blocks(path: Path):
    model = json.loads(path.read_text())
    result = []
    for block in model["blocks"]:
        degrees = block["degrees"]
        shape = tuple(p + 1 for p in degrees)
        ctrl = np.asarray(block["controls"], dtype=float).reshape((*shape, 4))
        knots = [[0.0] * (p + 1) + [1.0] * (p + 1) for p in degrees]
        result.append(
            {
                "degrees": degrees,
                "knots": knots,
                "ctrl": ctrl,
                "closed_u": False,
                "sectors": [block["id"]],
            }
        )
    return result


def basis_at(knots, degree, nctrl, x):
    eye = np.eye(nctrl)
    spline = BSpline(knots, eye, degree, extrapolate=False, axis=0)
    return np.asarray(spline(x)), np.asarray(spline.derivative()(x))


def knot_spans(knots, degree):
    return [(float(a), float(b)) for a, b in zip(knots[:-1], knots[1:]) if b > a]


def control_keys(patches):
    keys, counts = [], []
    for p in patches:
        euclid = p["ctrl"][..., :3] / p["ctrl"][..., 3, None]
        # Geometry coordinates are rational and exact in these fixtures. The
        # decimal key only unifies equal points after an overall homogeneous scale.
        arr = [tuple(np.round(q, 9)) for q in euclid.reshape(-1, 3)]
        keys.extend(arr)
        counts.append(len(arr))
    uniq = sorted(set(keys))
    ids = {key: i for i, key in enumerate(uniq)}
    merged = [
        np.asarray(
            [ids[k] for k in keys[sum(counts[:i]) : sum(counts[: i + 1])]], dtype=int
        )
        for i in range(len(patches))
    ]
    return merged, len(uniq), uniq


def quadrature(order):
    x, w = np.polynomial.legendre.leggauss(order)
    return (x + 1) * 0.5, w * 0.5


def assemble(patches, tie_closed_seam=True, order=5):
    maps, ndof, key_list = control_keys(patches)
    extra_coords = []
    if not tie_closed_seam:
        # Split only duplicated first/last u rows of the one closed annulus patch.
        occurrences = {}
        for patch in patches:
            euclid = patch["ctrl"][..., :3] / patch["ctrl"][..., 3, None]
            for q in euclid.reshape(-1, 3):
                key = tuple(np.round(q, 9))
                occurrences[key] = occurrences.get(key, 0) + 1
        for pi, p in enumerate(patches):
            if not p["closed_u"]:
                continue
            shape = p["ctrl"].shape[:3]
            ids = maps[pi].reshape(shape).copy()
            fresh = ndof
            for j in range(shape[1]):
                for k in range(shape[2]):
                    point = p["ctrl"][-1, j, k, :3] / p["ctrl"][-1, j, k, 3]
                    key = tuple(np.round(point, 9))
                    # Keep ties to other patches (notably the core) intact;
                    # split only a seam pair that is unique to this patch.
                    first = p["ctrl"][0, j, k, :3] / p["ctrl"][0, j, k, 3]
                    first_key = tuple(np.round(first, 9))
                    if key == first_key and occurrences.get(key, 0) == 2:
                        ids[-1, j, k] = fresh
                        fresh += 1
                        extra_coords.append(point.tolist())
            maps[pi] = ids.reshape(-1)
            ndof = fresh
    rows, cols, vals, f = [], [], [], np.zeros(ndof)
    qx, qw = quadrature(order)
    zall = np.concatenate(
        [(p["ctrl"][..., 2] / p["ctrl"][..., 3]).ravel() for p in patches]
    )
    zmin, zmax = float(zall.min()), float(zall.max())
    radius = 5.0
    height = zmax - zmin

    for pi, p in enumerate(patches):
        deg = p["degrees"]
        shape = p["ctrl"].shape[:3]
        ctrl = p["ctrl"].reshape((-1, 4))
        gids = maps[pi]
        local_matrix = np.zeros((len(gids), len(gids)))
        spans = [knot_spans(p["knots"][a], deg[a]) for a in range(3)]
        for span in product(*spans):
            for iq in product(range(order), repeat=3):
                point = [
                    span[a][0] + qx[iq[a]] * (span[a][1] - span[a][0]) for a in range(3)
                ]
                wt = [qw[iq[a]] * (span[a][1] - span[a][0]) for a in range(3)]
                b, db = [], []
                for a in range(3):
                    ba, da = basis_at(p["knots"][a], deg[a], shape[a], point[a])
                    b.append(ba)
                    db.append(da)
                poly = np.einsum("i,j,k->ijk", *b)
                dpoly = np.stack(
                    (
                        np.einsum("i,j,k->ijk", db[0], b[1], b[2]),
                        np.einsum("i,j,k->ijk", b[0], db[1], b[2]),
                        np.einsum("i,j,k->ijk", b[0], b[1], db[2]),
                    )
                )
                wctrl = ctrl[:, 3].reshape(shape)
                den = float(np.sum(poly * wctrl))
                dden = np.asarray([np.sum(dpoly[a] * wctrl) for a in range(3)])
                hom = np.einsum("ijk,ijkc->c", poly, p["ctrl"])
                dhom = np.einsum("aijk,ijkc->ac", dpoly, p["ctrl"])
                xyz = hom[:3] / den
                jac = ((dhom[:, :3] * den - dden[:, None] * hom[None, :3]) / den**2).T
                det = float(np.linalg.det(jac))
                if det <= 1e-12:
                    raise ValueError(
                        f"nonpositive Jacobian {det} in patch {p['sectors']}"
                    )
                rshape = (poly * wctrl / den).ravel()
                dr = (
                    (
                        dpoly * wctrl[None, ...] * den
                        - poly[None, ...] * wctrl[None, ...] * dden[:, None, None, None]
                    )
                    / den**2
                ).reshape(3, -1)
                grad = np.linalg.solve(jac.T, dr).T
                measure = math.prod(wt) * det
                ke = (grad @ grad.T) * measure
                # Manufactured nonuniform source excites angular modes, exposing a broken seam.
                source = (
                    1.0
                    + 0.25 * xyz[0] / radius
                    + 0.15 * xyz[1] / radius
                    + 0.1 * (xyz[2] - zmin) / height
                )
                fe = rshape * source * measure
                local_matrix += ke
                np.add.at(f, gids, fe)
        rows.extend(np.repeat(gids, len(gids)))
        cols.extend(np.tile(gids, len(gids)))
        vals.extend(local_matrix.ravel())
    K = coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsr()
    K.eliminate_zeros()
    # Homogeneous Dirichlet on the two end disks; natural Neumann elsewhere.
    coords = np.asarray(
        key_list
        + extra_coords
        + [(math.nan,) * 3] * max(0, ndof - len(key_list) - len(extra_coords))
    )
    fixed = np.flatnonzero(
        np.isclose(coords[:, 2], zmin, atol=1e-8)
        | np.isclose(coords[:, 2], zmax, atol=1e-8)
    )
    free = np.setdiff1d(np.arange(ndof), fixed)
    u = np.zeros(ndof)
    u[free] = spsolve(K[free][:, free], f[free])
    residual = K @ u - f
    rinf = float(np.linalg.norm(residual[free], ord=np.inf))
    energy = float(u @ (K @ u))
    work = float(u @ f)
    # Per-patch coefficient values retain duplicates for the seam diagnostic.
    patch_coeff = [u[m] for m in maps]
    return {
        "K": K,
        "u": u,
        "f": f,
        "energy": energy,
        "work": work,
        "residual_linf_free": rinf,
        "free_dofs": len(free),
        "dofs": ndof,
        "maps": maps,
        "patch_coeff": patch_coeff,
        "keys": key_list,
        "zmin": zmin,
        "zmax": zmax,
        "patches": patches,
    }


def seam_jump(solved):
    for i, p in enumerate(solved["patches"]):
        if not p["closed_u"]:
            continue
        c = solved["patch_coeff"][i].reshape(p["ctrl"].shape[:3])
        # Evaluate the scalar field itself on both sides, not just control rows.
        d, knots = p["degrees"], p["knots"]
        shape = p["ctrl"].shape[:3]
        vals = []
        for v in np.linspace(float(knots[1][d[1]]), float(knots[1][-d[1] - 1]), 9):
            bv, _ = basis_at(knots[1], d[1], shape[1], v)
            for w in np.linspace(float(knots[2][d[2]]), float(knots[2][-d[2] - 1]), 9):
                bw, _ = basis_at(knots[2], d[2], shape[2], w)
                traces = []
                for upar in (float(knots[0][d[0]]), float(knots[0][-d[0] - 1])):
                    bu, _ = basis_at(knots[0], d[0], shape[0], upar)
                    b = np.einsum("i,j,k->ijk", bu, bv, bw)
                    den = float(np.sum(b * p["ctrl"][..., 3]))
                    traces.append(float(np.sum(b * p["ctrl"][..., 3] * c) / den))
                vals.append(traces[0] - traces[1])
        arr = np.asarray(vals)
        return float(np.max(np.abs(arr))), float(np.linalg.norm(arr)), int(arr.size)
    return 0.0, 0.0, 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--model",
        type=Path,
        default=Path("docs/site/box-cylinder/tall/cylinder.model.json"),
    )
    ap.add_argument(
        "--two",
        type=Path,
        default=Path("docs/site/box-cylinder/tall/cylinder-2.patches.json"),
    )
    ap.add_argument(
        "--three",
        type=Path,
        default=Path("docs/site/box-cylinder/tall/cylinder-3.patches.json"),
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=Path("docs/site/box-cylinder/tall/cylinder-iga.json"),
    )
    ap.add_argument(
        "--check",
        action="store_true",
        help="verify the report at --output without rewriting it",
    )
    ap.add_argument(
        "--write",
        action="store_true",
        help="write the deterministic report (default behavior)",
    )
    args = ap.parse_args()
    cases = {
        "15_bezier_cells": load_bezier_blocks(args.model),
        "3_nurbs_patches": load_patch(args.three),
        "2_nurbs_patches": load_patch(args.two),
    }
    solved = {name: assemble(patches) for name, patches in cases.items()}
    bad_seam = assemble(cases["2_nurbs_patches"], tie_closed_seam=False)
    base = solved["15_bezier_cells"]
    base_f = base["f"]
    results = {}
    for name, case in solved.items():
        common = sorted(set(base["keys"]) & set(case["keys"]))
        bmap = {k: base["u"][i] for i, k in enumerate(base["keys"])}
        cmap = {k: case["u"][i] for i, k in enumerate(case["keys"])}
        diff = np.asarray([bmap[k] - cmap[k] for k in common])
        matrix_delta = case["K"] - base["K"]
        matrix_relative = float(
            np.sqrt(matrix_delta.multiply(matrix_delta).sum())
            / max(np.sqrt(base["K"].multiply(base["K"]).sum()), 1e-300)
        )
        energy_relative = abs(case["energy"] - base["energy"]) / max(
            abs(base["energy"]), 1e-300
        )
        rhs_relative = float(
            np.linalg.norm(case["f"] - base_f) / max(np.linalg.norm(base_f), 1e-300)
        )
        results[name] = {
            "patch_count": len(cases[name]),
            "dofs": case["dofs"],
            "free_dofs": case["free_dofs"],
            "stiffness_nnz": int(case["K"].nnz),
            "quadratic_form_uTKu": case["energy"],
            "strain_energy_half_uTKu": 0.5 * case["energy"],
            "load_work_uTf": case["work"],
            "relative_stiffness_frobenius_difference_vs_15_cell": matrix_relative,
            "relative_rhs_l2_difference_vs_15_cell": rhs_relative,
            "relative_energy_difference_vs_15_cell": energy_relative,
            "relative_energy_work_gap": abs(case["energy"] - case["work"])
            / max(abs(case["energy"]), 1e-300),
            "free_residual_linf": case["residual_linf_free"],
            "common_control_value_max_abs_vs_15_cell": float(np.max(np.abs(diff))),
            "common_control_values_compared": len(common),
        }
    jmax, jnorm, jcount = seam_jump(bad_seam)
    tied_jump, tied_jump_l2, tied_samples = seam_jump(tied := solved["2_nurbs_patches"])
    report = {
        "schema": "tall_cylinder_independent_iga_poisson_v1",
        "status": "verified",
        "method": "SciPy rational tensor-product NURBS/Bézier scalar H1 Galerkin assembly and sparse direct solve",
        "problem": {
            "operator": "integral grad(N_i) dot grad(u) dV",
            "source": "1 + 0.25*x/5 + 0.15*y/5 + 0.1*(z-zmin)/height",
            "boundary": "homogeneous Dirichlet on both end disks; natural Neumann on cylindrical wall",
            "units": "dimensionless",
        },
        "inputs": {
            "source_model": str(args.model),
            "source_model_sha256": sha(args.model),
            "two_patch_controls": str(args.two),
            "two_patch_sha256": sha(args.two),
            "three_patch_controls": str(args.three),
            "three_patch_sha256": sha(args.three),
        },
        "quadrature": {
            "gauss_points_per_nonzero_knot_span_axis": 5,
            "same tensor rule used for all representations": True,
        },
        "cases": results,
        "seam_open_negative_control": {
            "periodic_seam_tied": False,
            "physical_trace_sample_count": jcount,
            "max_abs_physical_trace_jump": jmax,
            "l2_physical_trace_jump": jnorm,
            "quadratic_form_uTKu": bad_seam["energy"],
            "strain_energy_half_uTKu": 0.5 * bad_seam["energy"],
            "relative_energy_change_vs_tied": abs(bad_seam["energy"] - tied["energy"])
            / abs(tied["energy"]),
            "free_residual_linf": bad_seam["residual_linf_free"],
        },
        "closed_seam_positive_control": {
            "max_abs_physical_trace_jump": tied_jump,
            "l2_physical_trace_jump": tied_jump_l2,
            "sample_count": tied_samples,
        },
        "criteria": {
            "max_relative_energy_difference_vs_15_cells": 1e-8,
            "max_common_control_value_difference": 2e-8,
            "max_free_residual_linf": 2e-9,
            "negative_control_seam_jump_min": 1e-5,
        },
        "limitations": [
            "floating-point quadrature and solve are numerical evidence, not exact integration or a proof of global injectivity",
            "coordinate-key trace identification is independent of production interface numbering but uses a 1e-9 mm geometric key",
            "the manufactured scalar Poisson field tests H1 trace assembly; it is not a calibrated elasticity or stress analysis",
        ],
    }
    errs = []
    for name, c in results.items():
        if not all(math.isfinite(v) for v in c.values()):
            errs.append(name + ": nonfinite metric")
        if c["relative_rhs_l2_difference_vs_15_cell"] > 1e-8:
            errs.append(name + ": load mismatch")
        if c["dofs"] != 68 or c["common_control_values_compared"] != 68:
            errs.append(name + ": missing coefficients")
        if c["relative_energy_difference_vs_15_cell"] > 1e-8:
            errs.append(name + ": energy mismatch")
        if c["relative_stiffness_frobenius_difference_vs_15_cell"] > 1e-8:
            errs.append(name + ": stiffness mismatch")
        if c["common_control_value_max_abs_vs_15_cell"] > 2e-8:
            errs.append(name + ": solution coefficient mismatch")
        if c["free_residual_linf"] > 2e-9:
            errs.append(name + ": free residual")
    if jmax < 1e-5:
        errs.append("open seam negative control did not produce detectable jump")
    if tied_jump > 2e-8:
        errs.append("tied seam has a physical trace jump")
    report["errors"] = errs
    report["status"] = "verified" if not errs else "failed"
    rendered = json.dumps(report, indent=2) + "\n"
    if args.check:
        if not args.output.exists():
            raise SystemExit(f"report is missing: {args.output}")
        expected = json.loads(args.output.read_text())

        def compare(a, b, trail="report"):
            if isinstance(a, dict) and isinstance(b, dict):
                if set(a) != set(b):
                    raise SystemExit(f"{trail}: field names differ")
                for key in a:
                    compare(a[key], b[key], f"{trail}.{key}")
            elif isinstance(a, list) and isinstance(b, list):
                if len(a) != len(b):
                    raise SystemExit(f"{trail}: list lengths differ")
                for i, (x, y) in enumerate(zip(a, b)):
                    compare(x, y, f"{trail}[{i}]")
            elif (
                isinstance(a, (int, float))
                and not isinstance(a, bool)
                and isinstance(b, (int, float))
                and not isinstance(b, bool)
            ):
                if not math.isclose(float(a), float(b), rel_tol=2e-12, abs_tol=2e-10):
                    raise SystemExit(f"{trail}: {a!r} differs from stored {b!r}")
            elif a != b:
                raise SystemExit(f"{trail}: {a!r} differs from stored {b!r}")

        compare(report, expected)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(
        json.dumps(
            {
                "status": report["status"],
                "output": str(args.output),
                "cases": results,
                "seam_open": report["seam_open_negative_control"],
            }
        )
    )
    return 0 if not errs else 2


if __name__ == "__main__":
    raise SystemExit(main())
