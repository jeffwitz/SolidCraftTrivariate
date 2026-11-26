"""Refine the transverse v analysis direction on the ten original S2 charts."""

import argparse
import json
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.linalg import null_space
from scipy import sparse
from scipy.sparse.linalg import spsolve

import run as study
from nodal_forces import traction_forces


def refine(nets, ids, trace, n, q):
    # Each inserted plane has its own coefficients. Across neighbouring bands
    # identify the same full-face endpoint. A collapsed polar line has only
    # its original coefficient, regardless of the number of transverse knots.
    columns = {}
    local_maps, stiffness, loads, children = [], [], [], {}
    rank = trace.shape[1]
    for (band, col), b in ids.items():
        parent = nets[b].reshape(4, 2, 2, 4)
        assert np.all(parent[:, 0, :, 3] == parent[:, 1, :, 3])
        for j in range(n):
            planes = []
            for s in [j, j + 1]:
                v = s / n
                planes.append((1 - v) * parent[:, 0] + v * parent[:, 1])
            net = np.stack(planes, axis=1).reshape(16, 4)
            loaded = col == 2 and j == 0
            K, t, _ = study.assemble(
                [net], {(band, 2 if loaded else 1, j): 0}, q, fitted=True
            )
            stiffness.append(K)
            loads.append(t)
            children[band, col, j] = net
            rows = []
            for a in range(4):
                for s in [j, j + 1]:
                    for w in range(2):
                        pole = col == 1 and (band, a) in [
                            (0, 0),
                            (max(k[0] for k in ids), 3),
                        ]
                        if pole or s in (0, n):
                            side = int(s == n and not pole)
                            z = trace[16 * b + 4 * a + 2 * side + w]
                            rows.append(dict(zip(np.flatnonzero(z), z[z != 0])))
                        else:
                            # Successive bands share only their endpoint a=3/0.
                            key = (col, 3 * band + a, s, w)
                            rows.append(
                                {columns.setdefault(key, rank + len(columns)): 1.0}
                            )
            local_maps.append(rows)
    rr = []
    cc = []
    vv = []
    for i, row in enumerate(r for m in local_maps for r in m):
        for k, value in row.items():
            rr.append(i)
            cc.append(k)
            vv.append(value)
    T = sparse.kron(
        sparse.csr_matrix(
            (vv, (rr, cc)), shape=(len(local_maps) * 16, rank + len(columns))
        ),
        sparse.eye(3),
        format="csr",
    )
    K = sparse.block_diag(stiffness, format="csr")
    f = np.concatenate(loads)
    Kr = (T.T @ K @ T).tocsc()
    rhs = T.T @ f
    reduced = spsolve(Kr, rhs)
    u = T @ reduced
    relres = np.linalg.norm(T.T @ (K @ u - f)) / np.linalg.norm(rhs)
    assert relres < 1e-8
    fields = dict(zip(children, u.reshape(-1, 16, 3)))
    gaps = []
    for band, col in ids:
        if col == 1:
            continue
        right = col == 2
        keys = [(band, col, n - 1 if right else 0), (band, 1, 0 if right else n - 1)]
        rows = traction_forces(
            [children[k] for k in keys],
            {(band, col): 0, (band, 1): 1},
            [fields[k] for k in keys],
            q,
            False,
            study.basis,
            study.points,
            study.MU,
            study.LAMBDA,
        )
        gaps.extend(rows)
    return dict(
        intervals_v=n,
        gauss_order=q,
        free_vector_dofs=int(T.shape[1]),
        reduced_residual=float(relres),
        loaded_end_mean_ux_mm=float(f @ u / f[::3].sum()),
        max_gap_N=max(row["max_gap_N"] for row in gaps),
        faces=gaps,
    )


def tangential_split(nets, ids, count):
    # de Casteljau subdivision of each cubic in u. Keep the original basis
    # pullback for measuring all levels with exactly the same test functions.
    matrices = [np.eye(4)]
    for _ in range(count.bit_length() - 1):
        children = []
        for P in matrices:
            levels = [P]
            for _ in range(3):
                A = levels[-1]
                levels.append((A[:-1] + A[1:]) / 2)
            children.extend(
                [
                    np.stack([A[0] for A in levels]),
                    np.stack([A[-1] for A in levels[::-1]]),
                ]
            )
        matrices = children
    assert len(matrices) == count
    split = []
    new_ids = {}
    transfer = {}
    for (b, c), i in ids.items():
        parent = nets[i].reshape(4, 2, 2, 4)
        for j, D in enumerate(matrices):
            key = (b * count + j, c)
            net = np.einsum("ji,ivwk->jvwk", D, parent) * count**3
            assert np.max(abs(net)) < 2**53 and np.all(net == np.round(net))
            new_ids[key] = len(split)
            split.append(net.reshape(-1, 4))
            w = parent[:, 0, 0, 3]
            transfer[key] = (b, c, D * w[None, :] / (D @ w)[:, None])
    return split, new_ids, transfer


def original_face_moments(faces, transfer):
    accumulated = {}
    for f in faces:
        b, c, P = transfer[f["band"], f["wing_col"]]
        pair = accumulated.setdefault((b, c), np.zeros((2, 4, 2, 3)))
        for i, name in enumerate(["wing_N", "core_N"]):
            pair[i] += np.einsum("ji,jwk->iwk", P, np.reshape(f[name], (4, 2, 3)))
    rows = []
    for (b, c), (wing, core) in accumulated.items():
        gap = wing + core
        rows.append(
            dict(
                band=b,
                wing_col=c,
                max_gap_N=float(np.linalg.norm(gap, axis=2).max()),
                gap_N=gap.reshape(8, 3).tolist(),
            )
        )
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("binary", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--combined", action="store_true")
    args = ap.parse_args()
    pair = json.loads(
        (study.ROOT / "tests/references/g2_cross_sketches.recipe.json").read_text()
    )
    baseline = json.loads(Path(__file__).with_name("nodal-forces.json").read_text())
    records = []
    with tempfile.TemporaryDirectory() as tmp:
        inp, out = [Path(tmp) / name for name in ("in", "out")]
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
            blocks = json.loads(out.read_text())["blocks"]
            nets = [np.array(b["controls"], float) for b in blocks]
            ids = {
                tuple(
                    map(int, re.fullmatch(r"band(\d+)_col(\d+)_z.*", b["id"]).groups())
                ): i
                for i, b in enumerate(blocks)
            }
            original_nets, original_ids = nets, ids
            for n in [1, 2, 4] if args.combined else [1, 2, 4, 8]:
                nu = n if args.combined else 1
                nets, ids, transfer = tangential_split(original_nets, original_ids, nu)
                full, _ = study.constraints(nets, ids, np.linspace(0, 1, 7))
                Z = study.exact_trace_space(nets, ids, full)
                clamp = [
                    16 * b + i
                    for (_, c), b in ids.items()
                    if c == 0
                    for i in np.arange(16).reshape(4, 2, 2)[:, 1, :].ravel()
                ]
                trace = Z @ null_space(Z[clamp, :], rcond=1e-11)
                old = next(
                    r for r in baseline if r["rise_mm"] == rise and r["order"] == 14
                )
                r = refine(nets, ids, trace, n, 10)
                r["intervals_u"] = nu
                if nu > 1:
                    r["faces"] = original_face_moments(r["faces"], transfer)
                    r["max_gap_N"] = max(f["max_gap_N"] for f in r["faces"])
                if n == 1:
                    np.testing.assert_allclose(
                        r["max_gap_N"],
                        max(v["max_gap_N"] for v in old["partial_faces"]),
                        rtol=1e-8,
                    )
                if n == (4 if args.combined else 8):
                    check = refine(nets, ids, trace, n, 14)
                    if nu > 1:
                        check["max_gap_N"] = max(
                            f["max_gap_N"]
                            for f in original_face_moments(check["faces"], transfer)
                        )
                    r["quadrature_max_gap_change_N"] = abs(
                        r["max_gap_N"] - check["max_gap_N"]
                    )
                records.append(dict(rise_mm=rise, **r))
                print(
                    rise,
                    n,
                    r["free_vector_dofs"],
                    r["max_gap_N"],
                    r["loaded_end_mean_ux_mm"],
                    flush=True,
                )
    args.output.write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    main()
