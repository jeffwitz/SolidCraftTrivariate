"""Compare the original cylinder analysis space on compact geometry."""

import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import BSpline

import cylinder_torsion as reference

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "docs/site/box-cylinder/tall/cylinder-2.patches.json"
DEFAULT_COMPACT = ROOT / "build/axial-knot-reduction/cylinder-2.json"
REPORT = ROOT / "docs/results/axial-knot-reduction/torsion.json"


def _rational_basis(knots, degree, bounds, q, flip=False):
    """Return active indices, basis and derivative with respect to source q."""
    a, b = bounds
    if not a < b:
        raise ValueError("Invalid geometry parameter interval")
    lo, hi = (b, a) if flip else (a, b)
    midpoint = (lo + hi) / 2
    span = int(np.searchsorted(knots, midpoint, side="right") - 1)
    span = min(max(span, degree), len(knots) - degree - 2)
    ids = np.arange(span - degree, span + 1)
    spline = BSpline(knots, np.eye(len(knots) - degree - 1)[:, ids], degree)
    t = b - q * (b - a) if flip else a + q * (b - a)
    sign = -1 if flip else 1
    scale = sign * (b - a)
    return ids, spline(t), spline(t, nu=1) * scale


def geometry_basis(patch, span, q):
    """Evaluate compact NURBS geometry and its Jacobian in source-cell coordinates."""
    axis_map = span["axis_map"]
    flips = span["flips"]
    rows, derivs, ids_by_axis = [], [], []
    for patch_axis, (degree, knots, bounds) in enumerate(
        zip(patch["degrees"], patch["knots"], span["parameter_intervals"])
    ):
        source_axis = axis_map[patch_axis]
        axis_q = q[source_axis] if isinstance(q, (tuple, list)) else q
        ids, values, derivatives = _rational_basis(
            knots, degree, bounds, axis_q, flips[patch_axis]
        )
        rows.append(values)
        derivs.append(derivatives)
        ids_by_axis.append(ids)

    def tensor(factors):
        return np.einsum("ai,bj,ck->abcijk", *factors).reshape(
            -1, np.prod([len(ids) for ids in ids_by_axis])
        )

    n = tensor(rows)
    dn_patch = [
        tensor([derivs[i] if i == axis else rows[i] for i in range(3)])
        for axis in range(3)
    ]
    dn = [None, None, None]
    for patch_axis, source_axis in enumerate(axis_map):
        dn[source_axis] = dn_patch[patch_axis]
    dn = np.stack(dn, axis=-1)
    indices = np.ravel_multi_index(np.ix_(*ids_by_axis), patch["control_shape"]).ravel()
    ctrl = np.asarray(patch["controls_homogeneous"], dtype=float)[indices]
    w = n @ ctrl[:, 3]
    r = n * ctrl[:, 3] / w[:, None]
    dr = (
        dn * ctrl[None, :, 3, None]
        - r[:, :, None] * np.einsum("qia,i->qa", dn, ctrl[:, 3])[:, None, :]
    ) / w[:, None, None]
    xyz = ctrl[:, :3] / ctrl[:, 3, None]
    x = r @ xyz
    jac = np.einsum("ic,qia->qca", xyz, dr)
    return x, jac


def _topology(compact, original):
    """Use the original field controls and source traces without renumbering."""
    fields = []
    knots_by_id = {
        item["patch_id"]: item["knots"] for item in compact["analysis_knots"]
    }
    provenance_by_id = {item["id"]: item for item in compact["provenance"]["patches"]}
    geometries = []
    for original_patch, reduced_patch in zip(
        original["patches"], compact["patches"], strict=True
    ):
        if original_patch["id"] != reduced_patch["id"]:
            raise ValueError("Compact geometry and field patch identities differ")
        provenance = provenance_by_id[reduced_patch["id"]]
        field = dict(original_patch)
        field["knots"] = knots_by_id[field["id"]]
        field["spans"] = provenance["spans"]
        fields.append(field)
        geometry = dict(reduced_patch)
        geometry["spans"] = provenance["spans"]
        geometries.append(geometry)
    return fields, geometries


def _couplings(records):
    return [
        {
            side: [item[side]["source_block"], item[side]["side"]]
            for side in ("first", "second")
        }
        for item in records
    ]


def solve_mixed(fields, geometries, interfaces, order):
    offsets, labels, xyz = reference.topology(fields, interfaces)
    n = len(xyz)
    stiffness = np.zeros((3 * n, 3 * n))
    force = np.zeros((n, 3))
    qg, wg = leggauss(order)
    q, weights = (qg + 1) / 2, wg / 2
    torque = 0.0
    geo_by_id = {patch["id"]: patch for patch in geometries}
    field_by_id = {patch["id"]: (i, patch) for i, patch in enumerate(fields)}
    for patch_id, (patch_index, field_patch) in field_by_id.items():
        geom_patch = geo_by_id[patch_id]
        geo_spans = {span["source_block"]: span for span in geom_patch["spans"]}
        for span in field_patch["spans"]:
            block = span["source_block"]
            geo_span = geo_spans[block]
            indices, _, d_r, _, _, _ = reference.basis(field_patch, span, q, weights)
            x, jac = geometry_basis(geom_patch, geo_span, q)
            local = labels[offsets[patch_index] + indices]
            det = np.linalg.det(jac)
            if det.min() <= 0:
                raise ValueError(f"Nonpositive compact-geometry Jacobian on {block}")
            deriv = np.einsum("qia,qac->qic", d_r, np.linalg.inv(jac))
            w3 = np.einsum("i,j,k->ijk", weights, weights, weights).ravel()
            bmat = np.zeros((len(w3), 6, len(local), 3))
            for c in range(3):
                bmat[:, c, :, c] = deriv[:, :, c]
            for row, (a, b) in enumerate(((0, 1), (1, 2), (2, 0)), 3):
                bmat[:, row, :, a], bmat[:, row, :, b] = deriv[:, :, b], deriv[:, :, a]
            bmat = bmat.reshape(len(w3), 6, -1)
            dofs = (3 * local[:, None] + np.arange(3)).ravel()
            weighted = (reference.ELASTIC @ bmat) * (w3 * det)[:, None, None]
            stiffness[np.ix_(dofs, dofs)] += bmat.reshape(
                -1, len(dofs)
            ).T @ weighted.reshape(-1, len(dofs))

            if block.startswith("disk_upper_"):
                # Integrate the unchanged traction on the upper physical cap.
                _, r, _, _, _, _ = reference.basis(
                    field_patch, span, q, weights, top=True
                )
                # Re-evaluate geometry at source-cell top (q_w=1).
                top_q = (q, q, np.array([1.0]))
                top_x, top_jac = geometry_basis(geom_patch, geo_span, top_q)
                top_area_w = (
                    np.linalg.norm(np.cross(top_jac[:, :, 0], top_jac[:, :, 1]), axis=1)
                    * np.einsum("i,j,k->ijk", weights, weights, [1.0]).ravel()
                )
                radius = top_x[:, :2] - 20.0
                traction = reference.TWIST * np.stack(
                    (-radius[:, 1], radius[:, 0], np.zeros(len(radius))), axis=1
                )
                force[local] += r.T @ (top_area_w[:, None] * traction)
                torque += np.dot(
                    top_area_w,
                    radius[:, 0] * traction[:, 1] - radius[:, 1] * traction[:, 0],
                )

    free = np.repeat(xyz[:, 2] != -12, 3)
    rhs = force.ravel()
    displacement = np.zeros(3 * n)
    displacement[free] = np.linalg.solve(stiffness[np.ix_(free, free)], rhs[free])
    energy = float(displacement @ stiffness @ displacement / 2)
    torque_ref = reference.TWIST * np.pi * reference.RADIUS**4 / 2
    energy_ref = torque_ref * reference.TWIST * reference.HEIGHT / 2
    exact = (
        reference.TWIST
        * (xyz[:, 2] + 12)[:, None]
        * np.stack((20 - xyz[:, 1], xyz[:, 0] - 20, np.zeros(n)), axis=1)
    ).ravel()
    return {
        "vector_dofs": int(3 * n),
        "free_vector_dofs": int(free.sum()),
        "stiffness": stiffness,
        "rhs": rhs,
        "displacement": displacement,
        "energy_N_mm": energy,
        "torque_N_mm": float(torque),
        "torque_relative_error": float(abs(torque / torque_ref - 1)),
        "energy_relative_error": float(abs(energy / energy_ref - 1)),
        "displacement_relative_error": float(
            np.linalg.norm(displacement - exact) / np.linalg.norm(exact)
        ),
        "free_residual": float(
            np.linalg.norm((stiffness @ displacement - rhs)[free])
            / np.linalg.norm(rhs[free])
        ),
    }


def relative_norm(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(a), np.linalg.norm(b)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compact", type=Path, default=DEFAULT_COMPACT)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    original = json.loads(SOURCE.read_text())
    compact = json.loads(args.compact.read_text())
    if compact.get("schema") != "gallery_axial_knot_reduction_v1":
        raise ValueError("Unexpected compact-geometry schema")
    fields, geometries = _topology(compact, original)
    seams = original.get("self_seams", [])
    interfaces = _couplings(original["interfaces"] + seams)
    baseline = solve_mixed(
        fields,
        original["patches"],
        interfaces,
        16,
    )
    reduced = solve_mixed(fields, geometries, interfaces, 16)
    if (baseline["vector_dofs"], baseline["free_vector_dofs"]) != (
        reduced["vector_dofs"],
        reduced["free_vector_dofs"],
    ):
        raise ValueError("Compact geometry changed the analysis DOF space")
    comparison = {
        "stiffness_relative_frobenius": relative_norm(
            baseline["stiffness"], reduced["stiffness"]
        ),
        "rhs_relative_l2": relative_norm(baseline["rhs"], reduced["rhs"]),
        "displacement_relative_l2": relative_norm(
            baseline["displacement"], reduced["displacement"]
        ),
        "energy_relative": abs(reduced["energy_N_mm"] / baseline["energy_N_mm"] - 1),
    }
    thresholds = {
        "stiffness_relative_frobenius": 2e-10,
        "rhs_relative_l2": 2e-10,
        "displacement_relative_l2": 2e-10,
        "energy_relative": 2e-10,
    }
    if any(
        not np.isfinite(comparison[k]) or comparison[k] > limit
        for k, limit in thresholds.items()
    ):
        raise ValueError(f"Compact geometry changes the analysis result: {comparison}")
    report = {
        "schema": "gallery_axial_knot_reduction_torsion_v1",
        "analysis": "3D linear elasticity Saint-Venant torsion on the unchanged cylinder-2 field basis and trace DOFs",
        "compact_geometry": str(args.compact.name),
        "compact_geometry_sha256": hashlib.sha256(
            args.compact.read_bytes()
        ).hexdigest(),
        "source_cover_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "gauss_order": 16,
        "shear_modulus_MPa": 1,
        "twist_per_mm": reference.TWIST,
        "field_patch_count": len(fields),
        "field_cell_count": sum(len(p["spans"]) for p in fields),
        "vector_dofs": reduced["vector_dofs"],
        "free_vector_dofs": reduced["free_vector_dofs"],
        "baseline": {
            k: baseline[k]
            for k in (
                "torque_N_mm",
                "energy_N_mm",
                "torque_relative_error",
                "energy_relative_error",
                "displacement_relative_error",
                "free_residual",
            )
        },
        "compact": {
            k: reduced[k]
            for k in (
                "torque_N_mm",
                "energy_N_mm",
                "torque_relative_error",
                "energy_relative_error",
                "displacement_relative_error",
                "free_residual",
            )
        },
        "comparison": comparison,
        "acceptance_thresholds": thresholds,
        "limitations": [
            "Quadrature and floating-point assembly are numerical checks, not a proof of equality of the geometric maps.",
            "Geometry reduction is evaluated only for this cylinder-2 torsion load case.",
        ],
    }
    if args.write:
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        REPORT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    else:
        saved = json.loads(REPORT.read_text())
        for key in (
            "schema",
            "analysis",
            "compact_geometry_sha256",
            "source_cover_sha256",
            "gauss_order",
            "shear_modulus_MPa",
            "field_patch_count",
            "field_cell_count",
            "vector_dofs",
            "free_vector_dofs",
            "acceptance_thresholds",
        ):
            if saved.get(key) != report.get(key):
                raise ValueError(f"Report metadata changed for {key}")
        for group_name in ("baseline", "compact", "comparison"):
            for key, value in report[group_name].items():
                if not np.isclose(saved[group_name][key], value, rtol=2e-9, atol=1e-11):
                    raise ValueError(f"Report value changed: {group_name}.{key}")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
