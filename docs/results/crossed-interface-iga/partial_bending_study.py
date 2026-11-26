"""P322 h-refined elasticity study on the native ten-chart S2 union.

The partial interfaces are imposed by overdetermined, preimage-aware trace
collocation. Independent off-grid trace checks are reported; this is a
numerical trace coupling study, not an exact mortar certificate.
"""

import argparse
import json
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.linalg import null_space, solve

ROOT = Path(__file__).resolve().parents[3]
YOUNG, POISSON = 1000.0, 0.3
MU = YOUNG / (2 * (1 + POISSON))
LAMBDA = YOUNG * POISSON / ((1 + POISSON) * (1 - 2 * POISSON))
ELASTIC = np.diag([2 * MU] * 3 + [MU] * 3)
ELASTIC[:3, :3] += LAMBDA
PU, PV, PW = 3, 2, 2
NLOCAL = (PU + 1) * (PV + 1) * (PW + 1)


def elevate_linear_axis(net, axis):
    """Elevate a homogeneous degree-one axis to degree two."""
    moved = np.moveaxis(net, axis, 0)
    elevated = np.stack([moved[0], (moved[0] + moved[1]) / 2, moved[1]])
    return np.moveaxis(elevated, 0, axis)


def subcell_net(net, iv, iw, nv, nw):
    source = net.reshape(4, 2, 2, 4)
    v0, v1 = iv / nv, (iv + 1) / nv
    w0, w1 = iw / nw, (iw + 1) / nw
    planes = []
    for w in (w0, w1):
        row = []
        for v in (v0, v1):
            at_v = (1 - v) * source[:, 0] + v * source[:, 1]
            row.append((1 - w) * at_v[:, 0, :] + w * at_v[:, 1, :])
        planes.append(np.stack(row, axis=1))
    child = np.stack(planes, axis=2)
    child = elevate_linear_axis(child, 1)
    child = elevate_linear_axis(child, 2)
    return child.reshape(NLOCAL, 4)


def make_elements(model, nv, nw):
    elements, by_key, parents = [], {}, {}
    for block in model["blocks"]:
        match = re.fullmatch(r"band(\d+)_col(\d+)_z.*", block["id"])
        band, col = map(int, match.groups())
        parent = np.asarray(block["controls"], float)
        elevated_parent = elevate_linear_axis(
            elevate_linear_axis(parent.reshape(4, 2, 2, 4), 1), 2
        ).reshape(NLOCAL, 4)
        parents[band, col] = elevated_parent
        for iv in range(nv):
            for iw in range(nw):
                by_key[band, col, iv, iw] = len(elements)
                elements.append(
                    dict(
                        band=band,
                        col=col,
                        iv=iv,
                        iw=iw,
                        net=subcell_net(parent, iv, iw, nv, nw),
                        source_id=block["id"],
                    )
                )
    return elements, by_key, parents


def bernstein(degree, x):
    if degree == 3:
        t = x
        a = 1 - t
        b = np.stack([a**3, 3 * t * a**2, 3 * t**2 * a, t**3], axis=1)
        d = np.stack(
            [-3 * a**2, 3 * a**2 - 6 * t * a, 6 * t * a - 3 * t**2, 3 * t**2], axis=1
        )
    else:
        t = x
        a = 1 - t
        b = np.stack([a**2, 2 * t * a, t**2], axis=1)
        d = np.stack([-2 * a, 2 - 4 * t, 2 * t], axis=1)
    return b, d


def basis(net, uvw):
    u, v, w = uvw.T
    U, dU = bernstein(PU, u)
    V, dV = bernstein(PV, v)
    W, dW = bernstein(PW, w)

    def tensor(a, b, c):
        return np.einsum("qi,qj,qk->qijk", a, b, c).reshape(-1, NLOCAL)

    weights = net[:, 3]
    p = tensor(U, V, W) * weights
    denominator = p.sum(axis=1)
    N = p / denominator[:, None]
    dp = (
        np.stack([tensor(dU, V, W), tensor(U, dV, W), tensor(U, V, dW)], axis=-1)
        * weights[None, :, None]
    )
    dN = (dp - N[:, :, None] * dp.sum(axis=1)[:, None, :]) / denominator[:, None, None]
    xyz = net[:, :3] / weights[:, None]
    x = N @ xyz
    J = np.einsum("ia,qij->qaj", xyz, dN)
    return N, dN, x, J


def tensor_points(q, dimensions=3):
    x, w = leggauss(q)
    x, w = (x + 1) / 2, w / 2
    grid = np.meshgrid(*([x] * dimensions), indexing="ij")
    weights = np.meshgrid(*([w] * dimensions), indexing="ij")
    return np.stack(grid, axis=-1).reshape(-1, dimensions), np.prod(
        weights, axis=0
    ).ravel()


class DisjointSet:
    def __init__(self, n):
        self.parent = np.arange(n)

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def join(self, a, b):
        self.parent[self.find(a)] = self.find(b)


def local_index(element, a, b, c):
    return element * NLOCAL + (a * (PV + 1) + b) * (PW + 1) + c


def identify_full_traces(elements, ids, nv, nw):
    ds = DisjointSet(len(elements) * NLOCAL)
    # True C0 h-refinement in v and w within each chart.
    for (band, col, iv, iw), e in ids.items():
        if iv + 1 < nv:
            other = ids[band, col, iv + 1, iw]
            for a in range(PU + 1):
                for c in range(PW + 1):
                    ds.join(local_index(e, a, PV, c), local_index(other, a, 0, c))
        if iw + 1 < nw:
            other = ids[band, col, iv, iw + 1]
            for a in range(PU + 1):
                for b in range(PV + 1):
                    ds.join(local_index(e, a, b, PW), local_index(other, a, b, 0))
        if (band + 1, col, iv, iw) in ids:
            other = ids[band + 1, col, iv, iw]
            for b in range(PV + 1):
                for c in range(PW + 1):
                    ds.join(local_index(e, PU, b, c), local_index(other, 0, b, c))
    # The two central polar lines collapse all v coefficients at fixed w.
    for band, a in ((0, 0), (5, PU)):
        if (band, 1, 0, 0) not in ids:
            continue
        for global_w in range(nw * PW + 1):
            owners = []
            iw = min(global_w // PW, nw - 1)
            c = global_w - iw * PW
            for iv in range(nv):
                e = ids[band, 1, iv, iw]
                for b in range(PV + 1):
                    owners.append(local_index(e, a, b, c))
            for x in owners[1:]:
                ds.join(owners[0], x)
    roots = np.array([ds.find(i) for i in range(len(elements) * NLOCAL)])
    unique, compact = np.unique(roots, return_inverse=True)
    return compact, len(unique)


def element_point(elem, u, v, w):
    _, _, x, _ = basis(elem["net"], np.array([[u, v, w]], float))
    return x[0]


def find_u_breaks(elem, side, nw):
    """Partition u where either slave edge crosses a master w knot."""
    from scipy.optimize import brentq

    targets = 12 * np.arange(1, nw) / nw
    samples = np.linspace(0, 1, 65)
    roots = [0.0, 1.0]
    for w in np.linspace(0.0, 1.0, nw + 1):
        z = np.array([element_point(elem, u, side, w)[2] for u in samples])
        for target in targets:
            values = z - target
            if np.max(np.abs(values)) < 1e-12:
                continue
            for i in range(len(samples) - 1):
                if abs(values[i]) < 1e-12:
                    roots.append(samples[i])
                elif values[i] * values[i + 1] < 0:
                    roots.append(
                        brentq(
                            lambda u: element_point(elem, u, side, w)[2] - target,
                            samples[i],
                            samples[i + 1],
                            xtol=1e-14,
                        )
                    )
    roots = sorted(roots)
    return [r for i, r in enumerate(roots) if i == 0 or r - roots[i - 1] > 1e-10]


def trace_samples(slave_parent, side, nv, nw, u_count, w_count):
    """Over-sample on each segment split by slave and master w knots."""
    u_breaks = find_u_breaks(slave_parent, side, nw)
    us = []
    for lo, hi in zip(u_breaks[:-1], u_breaks[1:]):
        us.extend(np.linspace(lo, hi, u_count))
    us = sorted(set(round(float(x), 14) for x in us))
    samples = []
    for u in us:
        z0 = element_point(slave_parent, u, side, 0.0)[2]
        z1 = element_point(slave_parent, u, side, 1.0)[2]
        cuts = list(np.linspace(0.0, 1.0, nw + 1))
        if abs(z1 - z0) > 1e-13:
            cuts.extend((12 * k / nw - z0) / (z1 - z0) for k in range(1, nw))
        cuts = sorted(c for c in cuts if -1e-12 <= c <= 1 + 1e-12)
        cuts = [min(1.0, max(0.0, c)) for c in cuts]
        cuts = [c for i, c in enumerate(cuts) if i == 0 or c - cuts[i - 1] > 1e-11]
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            if hi - lo < 1e-12:
                continue
            for w in np.linspace(lo, hi, w_count):
                samples.append((u, float(w)))
    return samples, u_breaks


def partial_trace_matrix(
    elements, ids, parents, compact, nbase, nv, nw, u_count=21, w_count=4
):
    rows, max_segments = [], 0
    interface_meta = []
    max_geometry_gap = 0.0
    for band in (2, 3):
        for col in (0, 2):
            side = float(col == 2)
            slave_parent = dict(net=parents[band, col])
            samples, ubreaks = trace_samples(
                slave_parent, side, nv, nw, u_count, w_count
            )
            max_segments = max(max_segments, len(ubreaks) - 1)
            siv = 0 if col == 0 else nv - 1
            sv = 0.0 if col == 0 else 1.0
            miv = nv - 1 if col == 0 else 0
            mv = 1.0 if col == 0 else 0.0
            master_cells = set()
            for u, global_w in samples:
                siw = min(int(np.floor(global_w * nw)), nw - 1)
                sw = global_w * nw - siw
                slave_id = ids[band, col, siv, siw]
                xyz = element_point(elements[slave_id], u, sv, sw)
                global_mw = xyz[2] / 12.0
                if not -1e-10 <= global_mw <= 1 + 1e-10:
                    raise ValueError("partial interface left the central w range")
                global_mw = min(1.0, max(0.0, global_mw))
                miw = min(int(np.floor(global_mw * nw)), nw - 1)
                master_cells.add(miw)
                mw = global_mw * nw - miw
                master_id = ids[band, 1, miv, miw]
                master_x = element_point(elements[master_id], u, mv, mw)
                max_geometry_gap = max(
                    max_geometry_gap, float(np.linalg.norm(xyz - master_x))
                )
                Ns = basis(elements[slave_id]["net"], np.array([[u, sv, sw]]))[0][0]
                Nm = basis(elements[master_id]["net"], np.array([[u, mv, mw]]))[0][0]
                row = np.zeros(nbase)
                np.add.at(row, compact[slave_id * NLOCAL : (slave_id + 1) * NLOCAL], Ns)
                np.add.at(
                    row, compact[master_id * NLOCAL : (master_id + 1) * NLOCAL], -Nm
                )
                rows.append(row)
            interface_meta.append(
                dict(
                    band=band,
                    wing_col=col,
                    u_breaks=ubreaks,
                    samples=len(samples),
                    master_w_cells=sorted(master_cells),
                )
            )
    return np.asarray(rows), interface_meta, max_geometry_gap


def heldout_trace(elements, ids, parents, compact, Z, nv, nw, rise, count=37):
    max_jump_space = 0.0
    records = []
    for band in (2, 3):
        for col in (0, 2):
            side = float(col == 2)
            parent = dict(net=parents[band, col])
            samples, ubreaks = trace_samples(parent, side, nv, nw, count, 7)
            max_row = 0.0
            for u, global_w in samples:
                siv = 0 if col == 0 else nv - 1
                sv = side
                siw = min(int(np.floor(global_w * nw)), nw - 1)
                sw = global_w * nw - siw
                s = ids[band, col, siv, siw]
                x = element_point(elements[s], u, sv, sw)
                gm = min(1.0, max(0.0, x[2] / 12.0))
                miv, mv = (nv - 1, 1.0) if col == 0 else (0, 0.0)
                miw = min(int(np.floor(gm * nw)), nw - 1)
                mw = gm * nw - miw
                m = ids[band, 1, miv, miw]
                Ns = basis(elements[s]["net"], np.array([[u, sv, sw]]))[0][0]
                Nm = basis(elements[m]["net"], np.array([[u, mv, mw]]))[0][0]
                row = np.zeros(Z.shape[0])
                np.add.at(row, compact[s * NLOCAL : (s + 1) * NLOCAL], Ns)
                np.add.at(row, compact[m * NLOCAL : (m + 1) * NLOCAL], -Nm)
                max_row = max(max_row, float(np.linalg.norm(row @ Z)))
            max_jump_space = max(max_jump_space, max_row)
            records.append(
                dict(
                    band=band,
                    wing_col=col,
                    u_intervals=len(ubreaks) - 1,
                    heldout_points=len(samples),
                    max_unit_space_jump=max_row,
                )
            )
    return max_jump_space, records


def assemble(elements, ids, nv, nw, q):
    qp, qw = tensor_points(q, 3)
    qf, qfw = tensor_points(q, 2)
    rr, cc, vv = [], [], []
    force = np.zeros(len(elements) * NLOCAL * 3)
    volume = 0.0
    for e, elem in enumerate(elements):
        N, dN, _, J = basis(elem["net"], qp)
        det = np.linalg.det(J)
        if np.any(det <= 0):
            raise ValueError(f"nonpositive Jacobian on {elem['source_id']}")
        grad = dN @ np.linalg.inv(J)
        B = np.zeros((len(qp), 6, NLOCAL, 3))
        B[:, 0, :, 0] = grad[:, :, 0]
        B[:, 1, :, 1] = grad[:, :, 1]
        B[:, 2, :, 2] = grad[:, :, 2]
        for row, (i, j) in zip((3, 4, 5), ((0, 1), (1, 2), (2, 0))):
            B[:, row, :, i] = grad[:, :, j]
            B[:, row, :, j] = grad[:, :, i]
        B = B.reshape(-1, 6, NLOCAL * 3)
        Ke = np.einsum("qai,ab,qbj,q->ij", B, ELASTIC, B, qw * det, optimize=True)
        volume += float(qw @ det)
        dofs = (e * NLOCAL * 3 + np.arange(NLOCAL * 3)).reshape(NLOCAL, 3)
        rr.extend(np.repeat(dofs.ravel(), NLOCAL * 3))
        cc.extend(np.tile(dofs.ravel(), NLOCAL * 3))
        vv.extend(Ke.ravel())
        # The right end is the exposed v=0 face of the two col=2 wing charts.
        if elem["col"] == 2 and elem["iv"] == 0:
            face = np.column_stack([qf[:, 0], np.zeros(len(qf)), qf[:, 1]])
            Nf, _, _, Jf = basis(elem["net"], face)
            area = -np.cross(Jf[:, :, 2], Jf[:, :, 0])
            traction_x = area[:, 0]  # sigma=diag(1,0,0) MPa
            local = np.zeros((NLOCAL, 3))
            local[:, 0] = Nf.T @ (qfw * traction_x)
            force[dofs.ravel()] += local.ravel()
    K = sparse.coo_matrix((vv, (rr, cc)), shape=(len(force), len(force))).tocsr()
    return K, force, volume


def solve_case(model, rise, nv, nw, q, trace_u_count=21):
    elements, ids, parents = make_elements(model, nv, nw)
    compact, nbase = identify_full_traces(elements, ids, nv, nw)
    C, interfaces, interface_geometry_gap = partial_trace_matrix(
        elements, ids, parents, compact, nbase, nv, nw, u_count=trace_u_count, w_count=4
    )
    # SVD of the overdetermined partial-trace matrix; the rank is numerical and reported.
    _, s, vh = np.linalg.svd(C, full_matrices=True)
    cutoff = s[0] * 1e-12 if len(s) else 0.0
    rank = int(np.count_nonzero(s > cutoff))
    rank_loose = int(np.count_nonzero(s > (s[0] * 1e-10 if len(s) else 0.0)))
    Z = vh[rank:].T
    heldout_defect, heldout = heldout_trace(
        elements, ids, parents, compact, Z, nv, nw, rise
    )

    clamp = set()
    for (band, col, iv, iw), e in ids.items():
        if col == 0 and iv == nv - 1:
            for a in range(PU + 1):
                for c in range(PW + 1):
                    clamp.add(compact[local_index(e, a, PV, c)])
    clamp = np.array(sorted(clamp), dtype=int)
    Zbc = Z @ null_space(Z[clamp], rcond=1e-12)
    E = sparse.coo_matrix(
        (np.ones(len(compact)), (np.arange(len(compact)), compact)),
        shape=(len(compact), nbase),
    ).tocsr()
    Tscalar_unclamped = np.asarray(E @ Z)
    T_unclamped = np.kron(Tscalar_unclamped, np.eye(3))
    Tscalar = np.asarray(E @ Zbc)
    T = np.kron(Tscalar, np.eye(3))
    K, f, volume = assemble(elements, ids, nv, nw, q)
    symmetry_error = float(np.max(np.abs((K - K.T).data))) if (K - K.T).nnz else 0.0
    eig = np.linalg.eigvalsh(T_unclamped.T @ (K @ T_unclamped))
    rigid_modes = int(np.count_nonzero(np.abs(eig) < 1e-10 * eig[-1]))
    affine = np.concatenate(
        [
            (
                elem["net"][:, :3]
                / elem["net"][:, 3, None]
                * np.array([1.0, -POISSON, -POISSON])
                / YOUNG
            ).reshape(-1)
            for elem in elements
        ]
    )
    affine_energy_error = abs(float(affine @ (K @ affine)) / (volume / YOUNG) - 1)
    reduced_K = T.T @ (K @ T)
    reduced_f = T.T @ f
    reduced_u = solve(reduced_K, reduced_f, assume_a="pos")
    u = T @ reduced_u
    residual = K @ u - f
    relative_residual = float(
        np.linalg.norm(T.T @ residual) / np.linalg.norm(reduced_f)
    )
    local_u = u.reshape(len(elements), NLOCAL, 3)

    force_N = float(f[::3].sum())
    mean_tip_ux = float(f @ u / force_N)
    energy = float(0.5 * u @ (K @ u))
    solved_jump = 0.0
    for band in (2, 3):
        for col in (0, 2):
            side = float(col == 2)
            parent = dict(net=parents[band, col])
            samples, _ = trace_samples(parent, side, nv, nw, 37, 7)
            for gu, gw in samples:
                siw = min(int(np.floor(gw * nw)), nw - 1)
                sw = gw * nw - siw
                siv, sv = (0, 0.0) if col == 0 else (nv - 1, 1.0)
                s = ids[band, col, siv, siw]
                xyz = element_point(elements[s], gu, sv, sw)
                gm = min(1.0, max(0.0, xyz[2] / 12.0))
                miv, mv = (nv - 1, 1.0) if col == 0 else (0, 0.0)
                miw = min(int(np.floor(gm * nw)), nw - 1)
                mw = gm * nw - miw
                m = ids[band, 1, miv, miw]
                Ns = basis(elements[s]["net"], np.array([[gu, sv, sw]]))[0][0]
                Nm = basis(elements[m]["net"], np.array([[gu, mv, mw]]))[0][0]
                left = Ns @ local_u[s]
                right = Nm @ local_u[m]
                solved_jump = max(solved_jump, float(np.linalg.norm(left - right)))

    return dict(
        rise_mm=rise,
        degree=[PU, PV, PW],
        refine_v=nv,
        refine_w=nw,
        quadrature_order=q,
        charts=len(model["blocks"]),
        elements=len(elements),
        scalar_dofs_after_full_identifications=nbase,
        partial_trace_constraint_rows=int(C.shape[0]),
        partial_trace_rank=rank,
        partial_trace_rank_rcond_1e_10=rank_loose,
        scalar_dofs_after_partial_trace=int(Z.shape[1]),
        clamped_scalar_dofs=int(np.linalg.matrix_rank(Z[clamp], tol=1e-12)),
        free_scalar_dofs=int(Zbc.shape[1]),
        free_vector_dofs=int(3 * Zbc.shape[1]),
        rigid_modes=rigid_modes,
        affine_energy_error=affine_energy_error,
        trace_svd_rcond=1e-12,
        trace_svd_cutoff=float(cutoff),
        trace_operator_heldout_defect=heldout_defect,
        max_partial_face_geometry_mismatch_mm=interface_geometry_gap,
        stiffness_symmetry_error=symmetry_error,
        max_solved_partial_interface_jump_mm=solved_jump,
        reduced_equilibrium_residual=relative_residual,
        volume_mm3=float(volume),
        force_N=force_N,
        loaded_end_mean_ux_mm=mean_tip_ux,
        strain_energy_Nmm=energy,
        interfaces=interfaces,
        heldout_interfaces=heldout,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("binary", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--rises", type=int, nargs="+", default=[3, 0, -3])
    ap.add_argument("--v-refinements", type=int, nargs="+", default=[1, 2, 4])
    ap.add_argument("--w-refinements", type=int, nargs="+", default=[1])
    ap.add_argument("--quadrature", type=int, default=8)
    ap.add_argument(
        "--plan", type=Path, help="JSON list of explicit (rise,v,w,q) cases"
    )
    args = ap.parse_args()
    recipe = json.loads(
        (ROOT / "tests/references/g2_cross_sketches.recipe.json").read_text()
    )
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        inp, out = [Path(tmp) / n for n in ("request.json", "model.json")]
        if args.plan:
            cases = json.loads(args.plan.read_text())["cases"]
        else:
            cases = [
                dict(
                    rise_mm=rise,
                    refine_v=nv,
                    refine_w=nw,
                    quadrature_order=args.quadrature,
                )
                for rise in args.rises
                for nv in args.v_refinements
                for nw in args.w_refinements
            ]
        for case in cases:
            rise = int(case["rise_mm"])
            nv = int(case["refine_v"])
            nw = int(case["refine_w"])
            q = int(case["quadrature_order"])
            inp.write_text(
                json.dumps(
                    dict(
                        family="inclined_crossed_tools_v1",
                        sketch_pair=recipe,
                        transverse_rise_mm=rise,
                    )
                )
            )
            subprocess.run(
                [args.binary.resolve(), "crossed-boolean", inp, "union", out],
                check=True,
                capture_output=True,
            )
            model = json.loads(out.read_text())
            if len(model["blocks"]) != 10:
                raise ValueError(
                    f"expected native ten-map S2 union, received {len(model['blocks'])} maps"
                )
            results.append(solve_case(model, rise, nv, nw, q))
    args.output.write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
