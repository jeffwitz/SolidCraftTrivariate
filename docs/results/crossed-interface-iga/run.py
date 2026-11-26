"""3-D linear-elasticity test on the native crossed-tool union, in mm and N."""

import argparse
import json
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import block_diag, eigvalsh, null_space, solve

ROOT = Path(__file__).resolve().parents[3]
YOUNG, POISSON = 1000.0, 0.3
MU = YOUNG / (2 * (1 + POISSON))
LAMBDA = YOUNG * POISSON / ((1 + POISSON) * (1 - 2 * POISSON))
ELASTIC = np.diag([2 * MU] * 3 + [MU] * 3)
ELASTIC[:3, :3] += LAMBDA


def points(q, axis=None, side=None):
    t, w = leggauss(q)
    axes = [(t + 1) / 2] * 3
    weights = [w / 2] * 3
    if axis is not None:
        axes[axis] = [side]
        weights[axis] = [1.0]
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3), np.prod(
        np.meshgrid(*weights, indexing="ij"), axis=0
    ).ravel()


def basis(net, uvw):
    u, v, w = uvw.T
    a = 1 - u
    U = np.array([a**3, 3 * u * a * a, 3 * u * u * a, u**3]).T
    dU = np.array(
        [-3 * a * a, 3 * a * a - 6 * u * a, 6 * u * a - 3 * u * u, 3 * u * u]
    ).T
    V, W = (np.stack([1 - t, t], axis=1) for t in [v, w])
    d = np.broadcast_to([-1.0, 1.0], V.shape)

    def tensor(a, b, c):
        return np.einsum("qi,qj,qk->qijk", a, b, c).reshape(-1, 16)

    weights = net[:, 3]
    p = tensor(U, V, W) * weights
    D = p.sum(axis=1)
    N = p / D[:, None]
    dp = (
        np.stack([tensor(dU, V, W), tensor(U, d, W), tensor(U, V, d)], axis=-1)
        * weights[None, :, None]
    )
    dN = (dp - N[:, :, None] * dp.sum(axis=1)[:, None, :]) / D[:, None, None]
    xyz = net[:, :3] / weights[:, None]
    return N, dN, N @ xyz, np.einsum("ia,qij->qaj", xyz, dN)


def partial_pair(nets, slave, master, col, uv):
    # u follows the shared cubic. Native z is affine in master w.
    slave_v = float(col == 2)
    uvw = np.column_stack([uv[:, 0], np.full(len(uv), slave_v), uv[:, 1]])
    S, _, x, J = basis(nets[slave], uvw)
    uvw[:, 1] = 1 - slave_v
    uvw[:, 2] = x[:, 2] / 12
    M, _, xm, _ = basis(nets[master], uvw)
    assert np.max(np.abs(x - xm)) < 5e-12
    return S, M, J


def constraints(nets, ids, nodes):
    size = 16 * len(nets)
    full = []
    partial = []

    def stitch(rows, block, coeff, other, opposite):
        matrix = np.zeros((len(coeff), size))
        matrix[:, 16 * block : 16 * (block + 1)] = coeff
        matrix[:, 16 * other : 16 * (other + 1)] -= opposite
        rows.append(matrix)

    e = np.eye(16).reshape(4, 2, 2, 16)
    for (band, col), block in ids.items():
        if (band + 1, col) in ids:
            stitch(
                full,
                block,
                e[-1].reshape(4, 16),
                ids[band + 1, col],
                e[0].reshape(4, 16),
            )
        if col != 1:
            uv = np.array([(u, w) for u in nodes for w in (0.0, 1.0)])
            master = ids[band, 1]
            S, M, _ = partial_pair(nets, block, master, col, uv)
            stitch(partial, block, S, master, M)
    # Single physical displacement on each collapsed polar line.
    for band, end in [(min(b for b, c in ids), 0), (max(b for b, c in ids), -1)]:
        block = ids[band, 1]
        stitch(full, block, e[end, 0], block, e[end, 1])
    return np.vstack(full), np.vstack(partial)


def exact_trace_space(nets, ids, full):
    # Form constraints with the integer geometry gauges, then orthonormalize
    # only after their exact null space has been determined.
    import sympy as sp

    rows = [list(map(int, row)) for row in full]
    for (band, col), slave in ids.items():
        if col == 1:
            continue
        master = ids[band, 1]
        for k in range(7):
            u = sp.Rational(k, 6)
            U = [(1 - u) ** 3, 3 * u * (1 - u) ** 2, 3 * u * u * (1 - u), u**3]
            for side in [0, 1]:
                v = int(col == 2)
                indices = [4 * i + 2 * v + side for i in range(4)]
                S = [a * int(nets[slave][i, 3]) for a, i in zip(U, indices)]
                den = sum(S)
                w = sum(a * int(nets[slave][i, 2]) for a, i in zip(U, indices)) / (
                    12 * den
                )
                row = [0] * full.shape[1]
                for i, idx in enumerate(indices):
                    row[16 * slave + idx] = S[i] / den
                mden = sum(
                    a * int(nets[master][4 * i + 2 * (1 - v), 3])
                    for i, a in enumerate(U)
                )
                for i, a in enumerate(U):
                    for j, h in enumerate([1 - w, w]):
                        idx = 4 * i + 2 * (1 - v) + j
                        row[16 * master + idx] = (
                            -a * h * int(nets[master][idx, 3]) / mden
                        )
                rows.append(row)
    matrix = sp.polys.matrices.DomainMatrix.from_Matrix(
        sp.SparseMatrix(rows)
    ).to_field()
    vectors = matrix.nullspace().to_Matrix().T
    # Scale columns before conversion to avoid amplifying the chart gauge.
    A = np.array(vectors.tolist(), dtype=float)
    return np.linalg.qr(A / np.linalg.norm(A, axis=0))[0]


def assemble(nets, ids, q, fitted=False):
    # Manufactured u=(x,-nu*y,-nu*z)/E gives exact uniaxial stress of 1 MPa.
    sigma = np.diag([1.0, 0.0, 0.0])
    volume = 0.0
    matrices = []
    boundary_load = []
    traction = np.zeros((len(nets), 16, 3))
    uvw, weights = points(q)
    for block, net in enumerate(nets):
        _, dN, _, J = basis(net, uvw)
        det = np.linalg.det(J)
        assert np.all(det > 0)
        grad = dN @ np.linalg.inv(J)
        B = np.zeros((len(J), 6, 16, 3))
        for i in range(3):
            B[:, i, :, i] = grad[:, :, i]
        for k, (i, j) in enumerate([(0, 1), (1, 2), (2, 0)], 3):
            B[:, k, :, i] = grad[:, :, j]
            B[:, k, :, j] = grad[:, :, i]
        B = B.reshape(-1, 6, 48)
        volume += weights @ det
        matrices.append(
            np.einsum("qai,ab,qbj,q->ij", B, ELASTIC, B, weights * det, optimize=True)
        )
        load = np.zeros((16, 3))
        band, col, *_ = next(key for key, value in ids.items() if value == block)
        for axis in range(3):
            for side in [0, 1]:
                if fitted and (col, axis, side) != (2, 1, 0):
                    continue
                if axis == 0 and (band + 2 * side - 1, col) in ids:
                    continue
                if axis == 1 and col != 1 and side == int(col == 2):
                    continue
                face, fw = points(q, axis, side)
                N, _, _, Jf = basis(net, face)
                area = (2 * side - 1) * np.cross(
                    Jf[:, :, (axis + 1) % 3], Jf[:, :, (axis + 2) % 3]
                )
                load += N.T @ ((fw[:, None] * area) @ sigma)
                if col == 2 and axis == 1 and side == 0:
                    traction[block, :, 0] = N.T @ (fw * np.linalg.norm(area, axis=1))
        boundary_load.append(load)
    # The exposed parts of a master face carry traction. Remove only the
    # footprint covered by the slave, using the opposite slave normal.
    if fitted:
        return block_diag(*matrices), traction.ravel(), volume
    face, fw = points(q, 1, 0)
    for (band, col), slave in ids.items():
        if col == 1:
            continue
        master = ids[band, 1]
        _, M, J = partial_pair(nets, slave, master, col, face[:, [0, 2]])
        area = (2 * int(col == 2) - 1) * np.cross(J[:, :, 2], J[:, :, 0])
        boundary_load[master] += M.T @ ((fw[:, None] * area) @ sigma)
    return block_diag(*matrices), np.ravel(boundary_load), traction.ravel(), volume


def matching_reference(nets, ids, rise, q):
    # Extend both inclined planes throughout the central profile, including
    # the bands outside the transverse tool. All shared faces now match.
    fitted = []
    locations = {}
    for (band, col), index in ids.items():
        net = nets[index]
        for layer in range(3) if col == 1 else [1]:
            patch = net.copy().reshape(4, 2, 2, 4)
            if col == 1:
                for side in [0, 1]:
                    level = layer + side
                    row = patch[:, :, side, :]
                    wx, weight = row[:, :, 0], row[:, :, 3]
                    row[:, :, 2] = [0, 4, 8, 12][level] * weight
                    if level in (1, 2):
                        row[:, :, 2] += rise / 100 * wx
            locations[band, col, layer] = len(fitted)
            fitted.append(patch.reshape(16, 4))

    # Merge matching physical control rows along declared neighbours only.
    local = np.arange(16 * len(fitted)).reshape(-1, 4, 2, 2)
    parents = list(range(local.size))

    def root(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i

    def join(one, two):
        for i, j in zip(one.ravel(), two.ravel()):
            parents[root(i)] = root(j)

    for (b, c, l), i in locations.items():
        for axis, other in enumerate([(b + 1, c, l), (b, c - 1, l), (b, c, l + 1)]):
            if other not in locations:
                continue
            j = locations[other]
            hi = local[i].take(-1, axis=axis)
            lo = local[j].take(0, axis=axis)
            neti = fitted[i].reshape(4, 2, 2, 4).take(-1, axis=axis)
            netj = fitted[j].reshape(4, 2, 2, 4).take(0, axis=axis)
            np.testing.assert_allclose(
                neti / neti.flat[3], netj / netj.flat[3], rtol=1e-13, atol=1e-12
            )
            join(hi, lo)
        if b in (0, 5) and c == 1:
            pole = local[i, 0 if b == 0 else -1]
            join(pole[0], pole[1])
    fixed = {
        root(k)
        for (_, c, _), i in locations.items()
        if c == 0
        for k in local[i, :, 1, :].ravel()
    }
    groups = np.array([root(i) for i in range(local.size)])
    unique, compact = np.unique(groups, return_inverse=True)
    T = np.kron(
        np.eye(len(unique))[compact][:, ~np.isin(unique, list(fixed))], np.eye(3)
    )
    K, t, V = assemble(fitted, locations, q, fitted=True)
    kr = T.T @ K @ T
    u = T @ solve(kr, T.T @ t, assume_a="pos")
    return dict(
        fitted_charts=len(fitted),
        fitted_scalar_dofs=len(unique),
        fitted_volume_mm3=float(V),
        fitted_mean_ux_mm=float(t @ u / t[::3].sum()),
    )


def run(model, q, rise):
    blocks = model["blocks"]
    assert all(abs(c) < 2**53 for b in blocks for h in b["controls"] for c in h)
    nets = [np.asarray(b["controls"], float) for b in blocks]
    ids = {
        tuple(map(int, re.fullmatch(r"band(\d+)_col(\d+)_z.*", b["id"]).groups())): i
        for i, b in enumerate(blocks)
    }
    nd = 16 * len(nets)
    full, _ = constraints(nets, ids, np.linspace(0, 1, 7))
    # Off-collocation evaluations test continuity for the whole constrained space.
    _, holdout = constraints(nets, ids, (leggauss(17)[0] + 1) / 2)
    _, coarse = constraints(nets, ids, np.linspace(0, 1, 4))
    Z = exact_trace_space(nets, ids, full)
    loose = null_space(np.vstack([full, coarse]), rcond=1e-11)
    leakage = lambda z: float(np.max(np.linalg.norm(holdout @ z, axis=1)))
    K, f, t, V = assemble(nets, ids, q)
    exact = (
        np.vstack([n[:, :3] / n[:, 3:] for n in nets]) * [1, -POISSON, -POISSON] / YOUNG
    ).ravel()
    Z3 = np.kron(Z, np.eye(3))

    def modes(z):
        z3 = np.kron(z, np.eye(3))
        ev = eigvalsh(z3.T @ K @ z3)
        return int(np.sum(abs(ev) < 1e-10 * ev[-1]))

    disconnected = null_space(full, rcond=1e-11)
    patch_res = np.linalg.norm(Z3.T @ (K @ exact - f)) / np.linalg.norm(Z3.T @ f)
    # Clamp the two left end patches (x=0); pull the opposite end at 1 MPa.
    clamp = []
    for (band, col), block in ids.items():
        if col == 0:
            clamp.extend(16 * block + np.arange(16).reshape(4, 2, 2)[:, 1, :].ravel())
    Zbc = np.kron(Z @ null_space(Z[clamp, :], rcond=1e-11), np.eye(3))
    Kr = Zbc.T @ K @ Zbc
    fr = Zbc.T @ t
    u = Zbc @ solve(Kr, fr, assume_a="pos")
    force = t[::3].sum()
    # Rigid virtual translation of the entire right wing measures interface
    # resultant via its internal elastic forces minus external end traction.
    right = np.concatenate(
        [np.arange(48 * b, 48 * (b + 1), 3) for (_, c), b in ids.items() if c == 2]
    )
    internal = K @ u
    exchange = (internal - t)[right].sum()
    metrics = dict(
        charts=len(nets),
        scalar_dofs=Z.shape[1],
        rigid_modes=modes(Z),
        disconnected_modes=modes(disconnected),
        trace_operator_defect=leakage(Z),
        four_node_defect=leakage(loose),
        affine_residual=patch_res,
        affine_space_error=float(
            np.linalg.norm(exact - Z3 @ (Z3.T @ exact)) / np.linalg.norm(exact)
        ),
        relative_affine_energy_error=abs(exact @ K @ exact / (V / YOUNG) - 1),
        force_N=force,
        right_interface_force_N=exchange,
        loaded_end_mean_ux_mm=float(t @ u / force),
        strain_energy_Nmm=float(0.5 * u @ K @ u),
        max_solved_jump_mm=float(np.max(np.abs(holdout @ u.reshape(nd, 3)))),
        reduced_residual=float(
            np.linalg.norm(Zbc.T @ (internal - t)) / np.linalg.norm(fr)
        ),
    )
    assert metrics["rigid_modes"] == 6 and metrics["disconnected_modes"] == 18
    assert (
        max(
            patch_res,
            metrics["relative_affine_energy_error"],
            metrics["trace_operator_defect"],
            metrics["reduced_residual"],
            metrics["max_solved_jump_mm"],
        )
        < 1e-9
    ), metrics
    assert metrics["affine_space_error"] < 1e-11
    if rise:
        assert metrics["four_node_defect"] > 1e-5
    assert abs(exchange + force) < 1e-7
    metrics.update(matching_reference(nets, ids, rise, q))
    assert abs(metrics["fitted_volume_mm3"] - V) < 1e-7
    metrics["relative_compliance_difference"] = (
        metrics["loaded_end_mean_ux_mm"] / metrics["fitted_mean_ux_mm"] - 1
    )
    return metrics


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("binary", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    pair = json.loads(
        (ROOT / "tests/references/g2_cross_sketches.recipe.json").read_text()
    )
    rows = []
    with tempfile.TemporaryDirectory() as directory:
        inp, out = [Path(directory) / n for n in ["in.json", "charts.json"]]
        for rise in [-3, 0, 3]:
            inp.write_text(
                json.dumps(
                    dict(
                        family="inclined_crossed_tools_v1",
                        sketch_pair=pair,
                        transverse_rise_mm=rise,
                    )
                )
            )
            subprocess.run(
                [args.binary.resolve(), "crossed-boolean", inp, "union", out],
                check=True,
                capture_output=True,
            )
            for q in [8, 14]:
                rows.append(
                    dict(
                        rise_mm=rise,
                        order=q,
                        **run(json.loads(out.read_text()), q, rise),
                    )
                )
    args.output.write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
