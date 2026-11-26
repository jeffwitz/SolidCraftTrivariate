"""Diagnostics for the S2 elasticity study. All force components use world XYZ."""

import numpy as np
from scipy.linalg import qr


def residual_forces(ids, residual, Zfree):
    # Use actual wing trace coefficients as independent virtual coordinates.
    # Shared band endpoints appear once. Interior variations have zero work.
    labels, candidates = [], []
    wing = np.zeros(Zfree.shape[0], dtype=bool)
    for (band, col), i in ids.items():
        if col == 1:
            continue
        wing[16 * i : 16 * i + 16] = True
        for a in range(4):
            for w in range(2):
                candidates.append(16 * i + 4 * a + 2 * int(col == 2) + w)
                labels.append([band, col, a, w])
    trace = Zfree[candidates]
    _, R, perm = qr(trace.T, pivoting=True, mode="economic")
    rank = np.sum(abs(R.diagonal()) > 1e-10)
    idx = perm[:rank]
    lift = Zfree @ np.linalg.pinv(trace[idx], rcond=1e-12)
    np.testing.assert_allclose(
        lift[np.asarray(candidates)[idx]], np.eye(rank), atol=1e-12
    )
    force = residual.reshape(-1, 3)
    W = lift[wing].T @ force[wing]
    M = lift[~wing].T @ force[~wing]
    gap = W + M
    return dict(
        independent_trace_coordinates=int(rank),
        max_gap_N=float(np.linalg.norm(gap, axis=1).max()),
        max_wing_force_N=float(np.linalg.norm(W, axis=1).max()),
        nodes=[
            dict(
                wing_coefficient=labels[k],
                wing_N=w.tolist(),
                core_N=m.tolist(),
                gap_N=g.tolist(),
            )
            for k, w, m, g in zip(idx, W, M, gap)
        ],
    )


def traction_forces(nets, ids, disp, q, fitted, basis, points, mu, lam):
    """Use eight wing face rational basis functions as common physical test weights."""

    def traction(i, uvw, normal):
        N, dN, x, J = basis(nets[i], uvw)
        grad_u = np.einsum("ia,qib->qab", disp[i], dN @ np.linalg.inv(J))
        stress = mu * (grad_u + grad_u.transpose(0, 2, 1))
        stress += lam * np.trace(grad_u, axis1=1, axis2=2)[:, None, None] * np.eye(3)
        return N, x, np.einsum("qab,qb->qa", stress, normal)

    rows = []
    for key, slave in ids.items():
        band, col = key[:2]
        if col == 1:
            continue
        master = ids[(band, 1, 1) if fitted else (band, 1)]
        side = int(col == 2)
        face, fw = points(q, 1, side)
        _, _, x, J = basis(nets[slave], face)
        area = (2 * side - 1) * np.cross(J[:, :, 2], J[:, :, 0])
        dA = np.linalg.norm(area, axis=1)
        normal = area / dA[:, None]
        N, _, ts = traction(slave, face, normal)
        master_face = face.copy()
        master_face[:, 1] = 1 - side
        if not fitted:
            master_face[:, 2] = x[:, 2] / 12
        _, xm, tm = traction(master, master_face, -normal)
        np.testing.assert_allclose(x, xm, atol=5e-12, rtol=0)
        weights = (
            N.reshape(-1, 4, 2, 2)[:, :, side, :].reshape(-1, 8) * (fw * dA)[:, None]
        )
        S = weights.T @ ts
        M = weights.T @ tm
        G = S + M
        rows.append(
            dict(
                band=band,
                wing_col=col,
                wing_N=S.tolist(),
                core_N=M.tolist(),
                gap_N=G.tolist(),
                max_gap_N=float(np.linalg.norm(G, axis=1).max()),
                resultant_gap_N=G.sum(axis=0).tolist(),
                relative_L2_traction_jump=float(
                    np.sqrt(
                        (fw * dA)
                        @ np.sum((ts + tm) ** 2, axis=1)
                        / ((fw * dA) @ (0.5 * np.sum(ts**2 + tm**2, axis=1)))
                    )
                ),
            )
        )
    return rows
