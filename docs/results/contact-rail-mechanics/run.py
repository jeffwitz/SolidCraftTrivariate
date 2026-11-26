"""Anti-plane shear on one triangle, with P2 FEM or a singular Q2 spline chart."""

import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import BSpline
from scipy.linalg import solve, eigvalsh


def field(x, y, affine):
    a, b = (0.0, 0.0) if affine else (0.6, 0.4)
    e = np.exp(a * x + b * y)
    l = 1 - x - y
    return (
        l * e,
        np.stack((e * (a * l - 1), e * (b * l - 1)), axis=-1),
        e * (2 * (a + b) - (a * a + b * b) * l),
    )


def gauss(q):
    t, w = leggauss(q)
    return (t + 1) / 2, w / 2


def traction(x, y, normal, affine):
    return field(x, y, affine)[1] @ normal


def assemble(kind, n, q, affine, solution=None):
    t, w = gauss(q)
    a, b = np.meshgrid(t, t, indexing="ij")
    a, b = a.ravel(), b.ravel()
    prod = np.outer(w, w).ravel()
    if kind == "IGA":
        knots = np.r_[0, 0, np.linspace(0, 1, n + 1), 1, 1]
        m = n + 2
        spline = BSpline(knots, np.eye(m), 2)
        deriv = spline.derivative()
        dofs = m * m
        free = np.arange((m - 1) * m)
    else:
        points = [(i, j) for i in range(2 * n + 1) for j in range(2 * n + 1 - i)]
        index = {p: i for i, p in enumerate(points)}
        dofs = len(points)
        free = np.array([k for k, (i, j) in enumerate(points) if i + j != 2 * n])

    K = np.zeros((dofs, dofs))
    f = np.zeros(dofs)
    norms = np.zeros(4)

    def volume(ids, N, grad, x, y, weight):
        u, g, load = field(x, y, affine)
        if solution is None:
            K[np.ix_(ids, ids)] += np.einsum(
                "qia,qja,q->ij", grad, grad, weight, optimize=True
            )
            f[ids] += N.T @ (weight * load)
        else:
            uh = N @ solution[ids]
            gh = np.einsum("qia,i->qa", grad, solution[ids])
            norms[:] += np.array(
                [
                    weight @ ((uh - u) ** 2),
                    weight @ (u * u),
                    weight @ np.sum((gh - g) ** 2, axis=1),
                    weight @ np.sum(g * g, axis=1),
                ]
            )

    if kind == "IGA":
        for i in range(n):
            for j in range(n):
                ids = np.array(
                    [(i + k) * m + j + l for k in range(3) for l in range(3)]
                )

                def basis(s, h):
                    U = spline(s)[:, i : i + 3]
                    V = spline(h)[:, j : j + 3]

                    def tensor(L, R):
                        return (L[:, :, None] * R[:, None, :]).reshape(-1, 9)

                    return (
                        tensor(U, V),
                        tensor(deriv(s)[:, i : i + 3], V),
                        tensor(U, deriv(h)[:, j : j + 3]),
                    )

                # Blow up the corner into two integration triangles: s+t
                # then contains a radial factor cancelled by the quadrature Jacobian.
                quadratures = (
                    [(a, a * b, prod * a), (a * b, a, prod * a)]
                    if i == j == 0
                    else [(a, b, prod)]
                )
                for s, h, weight in quadratures:
                    s, h = (i + s) / n, (j + h) / n
                    N, ds, dt = basis(s, h)
                    det = s + h
                    grad = np.stack(
                        (
                            (s[:, None] * ds - h[:, None] * dt) / det[:, None],
                            (ds + dt) / det[:, None],
                        ),
                        axis=-1,
                    )
                    volume(ids, N, grad, s - h, s * h, weight * det / n**2)
                if solution is not None:
                    continue
                for active, s, h, normal, length in [
                    (i == 0, 0 * t, (j + t) / n, [0, -1], 1),
                    (j == 0, (i + t) / n, 0 * t, [0, -1], 1),
                    (
                        j == n - 1,
                        (i + t) / n,
                        0 * t + 1,
                        np.array([-1, 1]) / math.sqrt(2),
                        math.sqrt(2),
                    ),
                ]:
                    if active:
                        N = basis(s, h)[0]
                        f[ids] += N.T @ (
                            w * length / n * traction(s - h, s * h, normal, affine)
                        )
    else:
        # Reference simplex coordinates map affinely from (-1,0),(1,0),(0,1).
        def xy(p):
            i, j = p
            return np.array([-1 + (2 * i + j) / (2 * n), j / (2 * n)])

        for i in range(n):
            for j in range(n - i):
                A = np.array([2 * i, 2 * j])
                B = A + [2, 0]
                C = A + [0, 2]
                tris = [(A, B, C)]
                if i + j < n - 1:
                    tris.append((B, A + [2, 2], C))
                for tri in tris:
                    pairs = [(0, 1), (1, 2), (2, 0)]
                    ids = [
                        index[tuple(p)]
                        for p in [*tri, *[(tri[k] + tri[l]) // 2 for k, l in pairs]]
                    ]
                    P = np.array([xy(p) for p in tri])
                    J = (P[1:] - P[0]).T
                    dL = np.array([[-1, -1], [1, 0], [0, 1]]) @ np.linalg.inv(J)

                    def shape(L):
                        return np.column_stack(
                            [
                                *(v * (2 * v - 1) for v in L.T),
                                *(4 * L[:, k] * L[:, l] for k, l in pairs),
                            ]
                        )

                    L = np.stack((1 - a - (1 - a) * b, a, (1 - a) * b), axis=1)
                    g = np.stack(
                        [
                            *((4 * L[:, k, None] - 1) * dL[k] for k in range(3)),
                            *(
                                4 * (L[:, k, None] * dL[l] + L[:, l, None] * dL[k])
                                for k, l in pairs
                            ),
                        ],
                        axis=1,
                    )
                    x, y = (L @ P).T
                    volume(ids, shape(L), g, x, y, prod * (1 - a) * np.linalg.det(J))
                    if solution is not None:
                        continue
                    for k, l in pairs:
                        if tri[k][1] == tri[l][1] == 0:
                            normal = [0, -1]
                        elif tri[k][0] == tri[l][0] == 0:
                            normal = np.array([-1, 1]) / math.sqrt(2)
                        else:
                            continue
                        L = np.zeros((q, 3))
                        L[:, k] = 1 - t
                        L[:, l] = t
                        x, y = (L @ P).T
                        f[ids] += shape(L).T @ (
                            w
                            * np.linalg.norm(P[l] - P[k])
                            * traction(x, y, normal, affine)
                        )
    if solution is not None:
        return norms
    A = K[np.ix_(free, free)]
    eig = eigvalsh(A)
    assert eig[0] > 0
    u = np.zeros(dofs)
    u[free] = solve(A, f[free], assume_a="pos")
    return u, dict(
        free_dofs=len(free),
        condition=eig[-1] / eig[0],
        energy=0.5 * u @ K @ u,
        residual=np.linalg.norm((K @ u - f)[free]) / np.linalg.norm(f[free]),
    )


def run(kind, n, q=10, affine=False):
    u, row = assemble(kind, n, q, affine)
    l, L, e, E = assemble(kind, n, 20, affine, solution=u)
    if kind == "IGA" and n == 8 and not affine:
        spline = BSpline(np.r_[0, 0, np.linspace(0, 1, n + 1), 1, 1], np.eye(n + 2), 2)
        coefficients = u.reshape(n + 2, n + 2)
        probes = []
        for r in [1e-2, 1e-3, 1e-4]:
            v, d = spline(r), spline.derivative()(r)
            shear = (d @ coefficients @ v + v @ coefficients @ d) / (2 * r)
            probes.append(
                dict(x=0.0, y=r * r, shear=shear, exact=field(0.0, r * r, False)[1][1])
            )
        row["shear_near_rail"] = probes
    return dict(
        method=kind,
        n=n,
        quadrature=q,
        affine=affine,
        **row,
        relative_L2=np.sqrt(l / L),
        relative_energy_norm=np.sqrt(e / E),
        reference_energy=E / 2,
    )


if __name__ == "__main__":
    rows = [run(k, n) for k in ["FEM", "IGA"] for n in (2, 4, 8, 16)]
    patch = [run(k, 4, affine=True) for k in ["FEM", "IGA"]]
    check = run("IGA", 8, q=18)
    for r in rows:
        print(
            r["method"],
            r["n"],
            r["free_dofs"],
            r["relative_L2"],
            r["relative_energy_norm"],
            r["condition"],
        )
    print("Patch", patch, "quadrature", check)
    assert all(r["relative_energy_norm"] < 1e-11 for r in patch)
    assert all(r["residual"] < 1e-10 for r in rows)
    assert (
        abs(check["relative_energy_norm"] / rows[6]["relative_energy_norm"] - 1) < 1e-6
    )
    Path(__file__).with_name("results.json").write_text(
        json.dumps(
            dict(refinement=rows, affine_patch=patch, quadrature_check=check), indent=2
        )
        + "\n"
    )
