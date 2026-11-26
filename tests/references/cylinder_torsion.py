"""Solve three-dimensional Saint-Venant torsion on each published cylinder cover."""

import argparse
import json
import os
from fractions import Fraction
from pathlib import Path

# This test uses small dense element matrices.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import BSpline

ROOT = Path(__file__).resolve().parents[2]
SITE = ROOT / "docs/site/box-cylinder/tall"
TWIST = 0.001
RADIUS, HEIGHT = 5.0, 36.0
# E=2.6, nu=0.3 gives shear modulus G=1 and first Lamé coefficient 1.5.
ELASTIC = np.diag([2.0, 2.0, 2.0, 1.0, 1.0, 1.0])
ELASTIC[:3, :3] += 1.5


def identity_patches(model):
    patches = []
    for b in model["blocks"]:
        p = b["degrees"]
        patches.append(
            dict(
                degrees=p,
                knots=[[0] * (d + 1) + [1] * (d + 1) for d in p],
                control_shape=[d + 1 for d in p],
                controls_homogeneous=b["controls"],
                spans=[
                    dict(
                        source_block=b["id"],
                        parameter_intervals=[[0, 1]] * 3,
                        axis_map=[0, 1, 2],
                        flips=[False] * 3,
                    )
                ],
            )
        )
    return patches


def topology(patches, interfaces):
    """Number local basis functions and identify the oriented source trace DOFs."""
    offsets, points, cells = [], [], {}
    for patch in patches:
        offsets.append(len(points))
        for point in patch["controls_homogeneous"]:
            values = [
                Fraction(*x) if isinstance(x, list) else Fraction(x) for x in point
            ]
            points.append(tuple(v / values[-1] for v in values[:3]))
        degree = np.array(patch["degrees"])
        for span in patch["spans"]:
            origin = np.array([a for a, _ in span["parameter_intervals"]]) * degree
            indices = np.indices(tuple(degree + 1)).reshape(3, -1).T
            ids = offsets[-1] + np.ravel_multi_index(
                (indices + origin).T, patch["control_shape"]
            )
            grid = ids.reshape(tuple(degree + 1))
            for axis, flip in enumerate(span["flips"]):
                if flip:
                    grid = np.flip(grid, axis=axis)
            cells[span["source_block"]] = grid.transpose(np.argsort(span["axis_map"]))
    parents = list(range(len(points)))

    def root(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i

    def trace(face):
        block, side = face
        ids = np.take(
            cells[block], -1 if side[1] == "1" else 0, axis="uvw".index(side[0])
        )
        order = {points[i]: i for i in ids.flat}
        if len(order) != ids.size:
            raise ValueError("Repeated trace control within one face")
        return order

    for face in interfaces:
        a, b = (trace(face[key]) for key in ("first", "second"))
        if a.keys() != b.keys():
            raise ValueError("Nonconforming control traces")
        for point in a:
            parents[root(a[point])] = root(b[point])
    _, labels = np.unique([root(i) for i in range(len(points))], return_inverse=True)
    coords = np.zeros((labels.max() + 1, 3))
    coords[labels] = np.asarray(points, dtype=float)
    return offsets, labels, coords


def basis(patch, span, quadrature, weights, top=False):
    rows, derivs = [], []
    local_ids = []
    for axis, (degree, knots, bounds) in enumerate(
        zip(patch["degrees"], patch["knots"], span["parameter_intervals"])
    ):
        a, b = bounds
        if b != a + 1:
            raise ValueError("Expected uniform parameter span")
        ids = np.arange(a * degree, (a + 1) * degree + 1)
        spline = BSpline(knots, np.eye(len(knots) - degree - 1)[:, ids], degree)
        t = np.array([float(b)]) if top and axis == 2 else a + quadrature
        rows.append(spline(t))
        derivs.append(spline(t, nu=1))
        local_ids.append(ids)

    def tensor(rows):
        return np.einsum("ai,bj,ck->abcijk", *rows).reshape(
            -1, np.prod([len(i) for i in local_ids])
        )

    N = tensor(rows)
    DN = np.stack(
        [
            tensor([derivs[i] if i == a else rows[i] for i in range(3)])
            for a in range(3)
        ],
        axis=-1,
    )
    indices = np.ravel_multi_index(np.ix_(*local_ids), patch["control_shape"]).ravel()
    ctrl = np.array(patch["controls_homogeneous"], dtype=float)[indices]
    W = N @ ctrl[:, 3]
    R = N * ctrl[:, 3] / W[:, None]
    dR = (
        DN * ctrl[None, :, 3, None]
        - R[:, :, None] * np.einsum("qia,i->qa", DN, ctrl[:, 3])[:, None, :]
    ) / W[:, None, None]
    xyz = ctrl[:, :3] / ctrl[:, 3, None]
    J = np.einsum("ic,qia->qca", xyz, dR)
    w = np.einsum("i,j,k->ijk", weights, weights, [1.0] if top else weights).ravel()
    return indices, R, dR, R @ xyz, J, w


def solve(patches, interfaces, order):
    offsets, labels, xyz = topology(patches, interfaces)
    n = len(xyz)
    K, F = np.zeros((3 * n, 3 * n)), np.zeros((n, 3))
    q, g = leggauss(order)
    torque = 0.0
    for pi, patch in enumerate(patches):
        for span in patch["spans"]:
            idx, _, dR, _, J, w = basis(patch, span, (q + 1) / 2, g / 2)
            local = labels[offsets[pi] + idx]
            der = np.einsum("qia,qac->qic", dR, np.linalg.inv(J))
            det = np.linalg.det(J)
            if det.min() <= 0:
                raise ValueError("Nonpositive sampled Jacobian")
            # Engineering strain ordering: xx, yy, zz, xy, yz, zx.
            B = np.zeros((len(w), 6, len(local), 3))
            for c in range(3):
                B[:, c, :, c] = der[:, :, c]
            for row, (a, b) in enumerate(((0, 1), (1, 2), (2, 0)), 3):
                B[:, row, :, a], B[:, row, :, b] = der[:, :, b], der[:, :, a]
            B = B.reshape(len(w), 6, -1)
            size = 3 * len(local)
            D = (ELASTIC @ B) * (w * det)[:, None, None]
            dofs = (3 * local[:, None] + np.arange(3)).ravel()
            K[np.ix_(dofs, dofs)] += B.reshape(-1, size).T @ D.reshape(-1, size)
            if span["source_block"].startswith("disk_upper_"):
                _, R, _, x, J, w = basis(patch, span, (q + 1) / 2, g / 2, top=True)
                area = np.linalg.norm(np.cross(J[:, :, 0], J[:, :, 1]), axis=1) * w
                r = x[:, :2] - 20.0
                traction = TWIST * np.stack(
                    (-r[:, 1], r[:, 0], np.zeros(len(r))), axis=1
                )
                F[local] += R.T @ (area[:, None] * traction)
                torque += np.dot(
                    area, r[:, 0] * traction[:, 1] - r[:, 1] * traction[:, 0]
                )
    free = np.repeat(xyz[:, 2] != -12, 3)
    force = F.ravel()
    u = np.zeros(3 * n)
    u[free] = np.linalg.solve(K[np.ix_(free, free)], force[free])
    phi = TWIST * (xyz[:, 2] + 12)
    exact = (
        phi[:, None] * np.stack((20 - xyz[:, 1], xyz[:, 0] - 20, np.zeros(n)), axis=1)
    ).ravel()
    internal = (K @ u).reshape(-1, 3)
    bottom = xyz[:, 2] == -12
    reaction = internal[bottom] - F[bottom]
    arm = xyz[bottom, :2] - 20
    clamp_torque = np.sum(arm[:, 0] * reaction[:, 1] - arm[:, 1] * reaction[:, 0])
    energy = float(u @ K @ u / 2)
    torque_ref = TWIST * np.pi * RADIUS**4 / 2
    energy_ref = torque_ref * TWIST * HEIGHT / 2
    return dict(
        vector_dofs=3 * n,
        free_vector_dofs=int(free.sum()),
        torque_N_mm=float(torque),
        support_torque_N_mm=float(clamp_torque),
        torque_balance_relative=float(abs(clamp_torque + torque) / torque),
        energy_N_mm=energy,
        torque_relative_error=abs(torque / torque_ref - 1),
        energy_relative_error=abs(energy / energy_ref - 1),
        displacement_relative_error=float(
            np.linalg.norm(u - exact) / np.linalg.norm(exact)
        ),
        free_residual=float(
            np.linalg.norm((K @ u - force)[free]) / np.linalg.norm(force[free])
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    model = json.loads((SITE / "cylinder.model.json").read_text())
    covers = {"bezier_15": (identity_patches(model), model["topology"]["interfaces"])}

    def couplings(records):
        return [
            {
                side: [iface[side]["source_block"], iface[side]["side"]]
                for side in ("first", "second")
            }
            for iface in records
        ]

    for n in (3, 2):
        cover = json.loads((SITE / f"cylinder-{n}.patches.json").read_text())
        seams = cover.get("self_seams", [])
        if n == 2 and not seams:
            raise ValueError("Closed annulus must declare its self seam")
        covers[f"nurbs_{n}"] = (
            cover["patches"],
            couplings(cover["interfaces"] + seams),
        )
    results = {name: solve(*data, 16) for name, data in covers.items()}
    fine = solve(*covers["nurbs_2"], 20)
    opened = solve(cover["patches"], couplings(cover["interfaces"]), 16)
    if (
        opened["vector_dofs"] <= fine["vector_dofs"]
        or opened["displacement_relative_error"] < 1e-5
    ):
        raise ValueError(
            "Removing the closed seam failed to expose the expected independent traces"
        )
    for r in [*results.values(), fine]:
        if (
            r["vector_dofs"] != results["bezier_15"]["vector_dofs"]
            or max(
                r[key]
                for key in (
                    "free_residual",
                    "torque_balance_relative",
                    "displacement_relative_error",
                    "torque_relative_error",
                    "energy_relative_error",
                )
            )
            > 1e-8
        ):
            raise ValueError(f"Torsion verification failed: {r}")
    report = dict(
        test="3D linear elasticity, exact circular-cylinder torsion",
        gauss_order=16,
        shear_modulus_MPa=1,
        twist_per_mm=TWIST,
        analytic_torque_N_mm=TWIST * np.pi * RADIUS**4 / 2,
        analytic_energy_N_mm=TWIST**2 * np.pi * RADIUS**4 * HEIGHT / 4,
        covers=results,
        closed_ring_order20=fine,
        negative_control_opened_seam=opened,
    )
    if args.write:
        (SITE / "cylinder-torsion.json").write_text(json.dumps(report, indent=2) + "\n")
    else:
        saved = json.loads((SITE / "cylinder-torsion.json").read_text())

        def compare(expected, actual):
            if isinstance(actual, dict):
                if expected.keys() != actual.keys():
                    raise ValueError("Torsion report fields differ")
                for key in actual:
                    compare(expected[key], actual[key])
            elif isinstance(actual, float):
                if not np.isclose(expected, actual, rtol=1e-8, atol=1e-11):
                    raise ValueError(f"Torsion report differs: {expected}, {actual}")
            elif expected != actual:
                raise ValueError("Torsion report metadata differs")

        compare(saved, report)
    print(json.dumps(report))


if __name__ == "__main__":
    main()
