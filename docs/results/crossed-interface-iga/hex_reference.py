"""Conforming, fully integrated HEX8 and HEX20 elasticity on S2 sampled from its ten charts."""

from functools import lru_cache
import argparse
import json
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import run as study

# Counterclockwise quadrilateral followed by the corresponding top four nodes.
BITS = np.array(
    [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ]
)


def shape(p):
    a = np.where(BITS, p, 1 - p)
    N = a.prod(axis=1)
    d = np.stack(
        [
            (2 * BITS[:, k] - 1) * a[:, (k + 1) % 3] * a[:, (k + 2) % 3]
            for k in range(3)
        ],
        axis=1,
    )
    return N, d


# VTK quadratic hexahedron edge order.
EDGES = np.array(
    [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ]
)


@lru_cache(maxsize=256)
def element_basis(p, quadratic):
    p = np.array(p)
    if not quadratic:
        return shape(p)
    t = 2 * p - 1
    signs = 2 * BITS - 1
    factors = 1 + signs * t
    base = factors.prod(axis=1) / 8
    multiplier = signs @ t - 2
    N = list(base * multiplier)
    d = list(
        np.column_stack(
            [
                signs[:, k]
                * (
                    factors[:, (k + 1) % 3] * factors[:, (k + 2) % 3] * multiplier / 8
                    + base
                )
                for k in range(3)
            ]
        )
    )
    for a, b in EDGES:
        axis = np.flatnonzero(BITS[a] != BITS[b])[0]
        j, k = (axis + 1) % 3, (axis + 2) % 3
        f = factors[a]
        N.append((1 - t[axis] ** 2) * f[j] * f[k] / 4)
        grad = np.zeros(3)
        grad[axis] = -t[axis] * f[j] * f[k] / 2
        grad[j] = (1 - t[axis] ** 2) * signs[a, j] * f[k] / 4
        grad[k] = (1 - t[axis] ** 2) * f[j] * signs[a, k] / 4
        d.append(grad)
    return np.array(N), 2 * np.asarray(d)


def elevate_hexes(X, C):
    pairs = np.sort(C[:, EDGES], axis=2).reshape(-1, 2)
    edges, indices = np.unique(pairs, axis=0, return_inverse=True)
    return np.vstack([X, X[edges].mean(axis=1)]), np.column_stack(
        [C, len(X) + indices.reshape(-1, 12)]
    )


def mesh(model, n, rise, *, circumferential=None, nv_by_col=None, height=None):
    xyz = []
    ids = {}
    conn = []
    loads = []
    traces = {}
    nets = {}

    def node(p):
        key = tuple(np.round(p, 10))
        if key not in ids:
            ids[key] = len(xyz)
            xyz.append(p)
        return ids[key]

    for block in model["blocks"]:
        b, c = map(int, re.fullmatch(r"band(\d+)_col(\d+)_z.*", block["id"]).groups())
        net = np.array(block["controls"], float).reshape(4, 2, 2, 4)
        nets[b, c] = net

        def xy(u, v):
            U = np.array(
                [(1 - u) ** 3, 3 * u * (1 - u) ** 2, 3 * u * u * (1 - u), u**3]
            )
            H = U @ ((1 - v) * net[:, 0, 0, :] + v * net[:, 1, 0, :])
            return H[:2] / H[3]

        mesh_height = n if height is None else height
        nv = [12, 3, 4][c] * n if nv_by_col is None else nv_by_col[c]
        nu = 3 * n if circumferential is None else circumferential
        for i in range(nu):
            for j in range(nv):
                # u grows with y and v grows from right to left.
                poly = np.array(
                    [
                        (i / nu, j / nv),
                        ((i + 1) / nu, j / nv),
                        ((i + 1) / nu, (j + 1) / nv),
                        (i / nu, (j + 1) / nv),
                    ]
                )
                p = np.array([xy(*uv) for uv in poly])
                # Coincident endpoints form a triangular fan at either pole.
                # Quadrangulate each triangle/quad using its edge midpoints.
                valid = np.linalg.norm(p - np.roll(p, -1, axis=0), axis=1) > 1e-10
                verts = p[valid]
                miduv = ((poly + np.roll(poly, -1, axis=0)) / 2)[valid]
                mids = (verts + np.roll(verts, -1, axis=0)) / 2
                center = np.mean(verts, axis=0)
                valid_indices = np.flatnonzero(valid)
                for h in range(len(verts)):
                    section = np.array([verts[h], mids[h], center, mids[h - 1]])
                    layers = range(3) if c == 1 else [1]
                    caps = np.array(
                        [
                            np.zeros(4),
                            4 + rise * section[:, 0] / 100,
                            8 + rise * section[:, 0] / 100,
                            np.full(4, 12),
                        ]
                    )
                    for l in layers:
                        for k in range(mesh_height):
                            heights = [
                                caps[l] + (caps[l + 1] - caps[l]) * (k + s) / mesh_height
                                for s in (0, 1)
                            ]
                            cell = [
                                node([x, y, z])
                                for zz in heights
                                for (x, y), z in zip(section, zz)
                            ]
                            assert len(set(cell)) == 8
                            eid = len(conn)
                            conn.append(cell)
                            # Edge h starts at vertex h. Edge h-1 ends there.
                            for ed, local, start in [
                                (h, 0, True),
                                ((h - 1) % len(verts), 3, False),
                            ]:
                                orig = valid_indices[ed]
                                va, vb = poly[orig], poly[(orig + 1) % 4]
                                if va[1] != vb[1]:
                                    continue
                                if c == 2 and va[1] == 0:
                                    loads.append((eid, local))
                                if b not in (2, 3) or l != 1:
                                    continue
                                side = None
                                if (c == 0 and va[1] == 0) or (c == 1 and va[1] == 1):
                                    side = 0
                                if (c == 2 and va[1] == 1) or (c == 1 and va[1] == 0):
                                    side = 2
                                if side is None:
                                    continue
                                # local edge is directed counterclockwise
                                endpoints = (
                                    (va[0], miduv[ed, 0])
                                    if start
                                    else (miduv[ed, 0], vb[0])
                                )
                                key = tuple(
                                    sorted(
                                        cell[t]
                                        for t in [
                                            local,
                                            (local + 1) % 4,
                                            local + 4,
                                            (local + 1) % 4 + 4,
                                        ]
                                    )
                                )
                                traces.setdefault(key, []).append(
                                    (eid, local, endpoints, k, b, side, c)
                                )
    assert all(len(p) == 2 for p in traces.values()) and len(traces) == 8 * nu * mesh_height, (
        len(traces),
        [v for v in traces.values() if len(v) != 2][:4],
    )
    return np.array(xyz), np.array(conn), loads, traces, nets


def face_point(edge, t, w):
    return [(t, 0, w), (1, t, w), (1 - t, 1, w), (0, 1 - t, w)][edge]


def solve_mesh(model, n, rise, vtk_directory=None, quadratic=False, mesh_options=None):
    options = {} if mesh_options is None else mesh_options
    X, C, loads, traces, nets = mesh(model, n, rise, **options)
    mesh_height = n if options.get("height") is None else options["height"]
    if quadratic:
        X, C = elevate_hexes(X, C)
    nodes_per_cell = C.shape[1]
    nd = 3 * nodes_per_cell
    fixed = np.flatnonzero(X[:, 0] < 1e-9)
    q, g = np.polynomial.legendre.leggauss(3 if quadratic else 2)
    q = (1 + q) / 2
    g = g / 2
    # Assemble element matrices in batches to bound peak memory.
    S = sparse.csr_matrix((X.size, X.size))
    volume = 0.0
    min_det = np.inf
    min_corner = np.inf
    for chunk in np.array_split(C, max(1, len(C) // 400)):
        P = X[chunk]
        ke = np.zeros((len(chunk), nd, nd))
        # For vertically extruded planar quads the determinant is the
        # bilinear height times the affine planar Jacobian. Positivity of
        # each factor at the four planar corners extends over the cube.
        assert np.allclose(P[:, :4, :2], P[:, 4:8, :2], atol=1e-12, rtol=0)
        assert np.all(P[:, 4:8, 2] > P[:, :4, 2])
        for p in BITS[:4]:
            _, dc = shape(p)
            planar_det = np.linalg.det((P[:, :8].transpose(0, 2, 1) @ dc)[:, :2, :2])
            assert np.all(planar_det > 0)
            min_corner = min(min_corner, planar_det.min())
        for a, wa in zip(q, g):
            for b, wb in zip(q, g):
                for c, wc in zip(q, g):
                    _, d = element_basis((a, b, c), quadratic)
                    J = P.transpose(0, 2, 1) @ d
                    det = np.linalg.det(J)
                    min_det = min(min_det, det.min())
                    assert det.min() > 0
                    grad = d @ np.linalg.inv(J)
                    B = np.zeros((len(chunk), 6, nodes_per_cell, 3))
                    for i in range(3):
                        B[:, i, :, i] = grad[:, :, i]
                    for idx, (i, j) in enumerate([(0, 1), (1, 2), (2, 0)], 3):
                        B[:, idx, :, i] = grad[:, :, j]
                        B[:, idx, :, j] = grad[:, :, i]
                    B = B.reshape(-1, 6, nd)
                    ke += (B.transpose(0, 2, 1) @ study.ELASTIC @ B) * (
                        det * wa * wb * wc
                    )[:, None, None]
                    volume += det.sum() * wa * wb * wc
        dofs = (3 * chunk[:, :, None] + np.arange(3)).reshape(-1, nd)
        r = np.broadcast_to(dofs[:, :, None], ke.shape).ravel()
        c = np.broadcast_to(dofs[:, None, :], ke.shape).ravel()
        S += sparse.coo_matrix((ke.ravel(), (r, c)), shape=S.shape).tocsr()
    F = np.zeros_like(X)
    for eid, edge in loads:
        nodes = C[eid]
        for t, wt in zip(q, g):
            for w, ww in zip(q, g):
                N, d = element_basis(tuple(face_point(edge, t, w)), quadratic)
                J = X[nodes].T @ d
                F[nodes, 0] += (
                    N * np.linalg.norm(np.cross(J[:, edge % 2], J[:, 2])) * wt * ww
                )
    affine = X * np.array([1, -study.POISSON, -study.POISSON]) / study.YOUNG
    patch_error = abs(
        affine.ravel() @ (S @ affine.ravel()) / (volume / study.YOUNG) - 1
    )
    assert patch_error < 1e-8
    free = np.repeat(X[:, 0] >= 1e-9, 3)
    f = F.ravel()
    u = np.zeros_like(f)
    u[free] = spsolve(S[free, :][:, free], f[free], permc_spec="MMD_AT_PLUS_A")
    residual = (S @ u - f).reshape(-1, 3)
    assert np.linalg.norm(residual[X[:, 0] >= 1e-9]) / np.linalg.norm(f) < 1e-7
    force = F[:, 0].sum()
    assert abs(force - 40.2222222222222) < 1e-8
    assert np.max(np.abs(residual[fixed].sum(axis=0) + [force, 0, 0])) < 1e-6
    # Consistent traction moments on the original eight face weights, evaluated
    # at corresponding face parameters on the planarized hex boundary.
    gaps = {}
    t4, w4 = np.polynomial.legendre.leggauss(4)
    for pair in traces.values():
        for eid, edge, endpoints, k, b, side, col in pair:
            nodes = C[eid]
            P = X[nodes]
            U = u.reshape(-1, 3)[nodes]
            G = gaps.setdefault((b, side), np.zeros((2, 8, 3)))
            cw = nets[b, side][:, 0, 0, 3]
            for t, wt in zip((1 + t4) / 2, w4 / 2):
                uc = (1 - t) * endpoints[0] + t * endpoints[1]
                B = (
                    np.array(
                        [
                            (1 - uc) ** 3,
                            3 * (1 - uc) ** 2 * uc,
                            3 * uc**2 * (1 - uc),
                            uc**3,
                        ]
                    )
                    * cw
                )
                B /= B.sum()
                for w, ww in zip((1 + t4) / 2, w4 / 2):
                    _, d = element_basis(tuple(face_point(edge, t, w)), quadratic)
                    J = P.T @ d
                    grad = U.T @ d @ np.linalg.inv(J)
                    sig = study.MU * (grad + grad.T) + study.LAMBDA * np.trace(
                        grad
                    ) * np.eye(3)
                    area = np.cross(J[:, edge % 2], J[:, 2]) * (1 if edge < 2 else -1)
                    height = (k + w) / mesh_height
                    weights = np.outer(B, [1 - height, height]).ravel()
                    G[int(col == 1)] += weights[:, None] * (sig @ area) * wt * ww
    if vtk_directory:
        vtk_directory.mkdir(parents=True, exist_ok=True)
        with (
            vtk_directory / f"crossed-union-hex{nodes_per_cell}-rise{rise}-n{n}.vtk"
        ).open("w") as file:
            file.write(
                "# vtk DataFile Version 3.0\nS2 hexahedral elasticity, mm\nASCII\nDATASET UNSTRUCTURED_GRID\n"
            )
            file.write(f"POINTS {len(X)} double\n")
            np.savetxt(file, X, fmt="%.12g")
            file.write(f"CELLS {len(C)} {(nodes_per_cell + 1) * len(C)}\n")
            np.savetxt(
                file, np.column_stack([np.full(len(C), nodes_per_cell), C]), fmt="%d"
            )
            file.write(
                f"CELL_TYPES {len(C)}\n" + ("25\n" if quadratic else "12\n") * len(C)
            )
            file.write(f"POINT_DATA {len(X)}\nVECTORS displacement_mm double\n")
            np.savetxt(file, u.reshape(-1, 3), fmt="%.12g")
    faces = [
        dict(
            band=b,
            wing=side,
            wing_N=g[0].tolist(),
            core_N=g[1].tolist(),
            gap_N=g.sum(axis=0).tolist(),
            max_gap_N=float(np.linalg.norm(g.sum(axis=0), axis=1).max()),
        )
        for (b, side), g in sorted(gaps.items())
    ]
    return dict(
        n=n,
        rise_mm=rise,
        element=f"HEX{nodes_per_cell}",
        hexes=len(C),
        free_dofs=int(free.sum()),
        volume_mm3=volume,
        force_N=force,
        affine_energy_relative_error=patch_error,
        min_corner_planar_det=min_corner,
        min_gauss_det=min_det,
        ux_mm=f @ u / force,
        energy_Nmm=0.5 * f @ u,
        max_gap_N=max(r["max_gap_N"] for r in faces),
        faces=faces,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("binary", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--levels", type=int, nargs="+", default=[1, 2, 3, 4])
    ap.add_argument("--rises", type=int, nargs="+", default=[0, 3])
    ap.add_argument("--vtk-directory", type=Path)
    ap.add_argument("--quadratic", action="store_true")
    args = ap.parse_args()
    rows = []
    pair = json.loads(
        (study.ROOT / "tests/references/g2_cross_sketches.recipe.json").read_text()
    )
    with tempfile.TemporaryDirectory() as d:
        inp, out = [Path(d) / f for f in ["in", "out"]]
        for rise in args.rises:
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
            for n in args.levels:
                row = solve_mesh(
                    json.loads(out.read_text()),
                    n,
                    rise,
                    args.vtk_directory,
                    args.quadratic,
                )
                rows.append(row)
                print(n, rise, row["hexes"], row["ux_mm"], row["max_gap_N"], flush=True)
                args.output.write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
