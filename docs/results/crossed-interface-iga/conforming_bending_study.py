"""Degree and mesh refinement study on the conforming 22-chart S2 union."""

import argparse
import importlib.util
import json
import math
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import spsolve

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("crossed_study", HERE / "run.py")
study = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(study)
RECIPE = study.ROOT / "tests/references/g2_cross_sketches.recipe.json"


def fitted_geometry(model, rise):
    blocks = model["blocks"]
    nets = [np.asarray(block["controls"], dtype=float) for block in blocks]
    ids = {
        tuple(map(int, re.fullmatch(r"band(\d+)_col(\d+)_z.*", b["id"]).groups())): i
        for i, b in enumerate(blocks)
    }
    fitted, locations = [], {}
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
    return [n.reshape(4, 2, 2, 4) for n in fitted], locations


def elevate_axis(net, axis, target_degree):
    degree = net.shape[axis] - 1
    if target_degree < degree:
        raise ValueError("This study only elevates the geometry degrees")
    while degree < target_degree:
        moved = np.moveaxis(net, axis, 0)
        elevated = np.empty((degree + 2, *moved.shape[1:]), dtype=float)
        elevated[0], elevated[-1] = moved[0], moved[-1]
        for i in range(1, degree + 1):
            elevated[i] = (i / (degree + 1)) * moved[i - 1] + (
                1 - i / (degree + 1)
            ) * moved[i]
        net = np.moveaxis(elevated, 0, axis)
        degree += 1
    return net


def split_curve(curve, t):
    levels = [curve]
    while len(levels[-1]) > 1:
        levels.append((1 - t) * levels[-1][:-1] + t * levels[-1][1:])
    degree = len(curve) - 1
    left = np.stack([levels[i][0] for i in range(degree + 1)])
    right = np.stack([levels[degree - i][i] for i in range(degree + 1)])
    return left, right


def restrict_curve(curve, a, b):
    if a == 0 and b == 1:
        return curve.copy()
    if a == 0:
        return split_curve(curve, b)[0]
    tail = split_curve(curve, a)[1]
    return split_curve(tail, (b - a) / (1 - a))[0]


def restrict_axis(net, axis, a, b):
    moved = np.moveaxis(net, axis, 0)
    shape = moved.shape
    curves = moved.reshape(shape[0], -1)
    result = np.stack(
        [restrict_curve(curves[:, i], a, b) for i in range(curves.shape[1])], axis=1
    )
    return np.moveaxis(result.reshape(shape), 0, axis)


def element_nets(parent, degrees, nv, nw):
    net = elevate_axis(parent, 1, degrees[1])
    net = elevate_axis(net, 2, degrees[2])
    elements = {}
    for iv in range(nv):
        vnet = restrict_axis(net, 1, iv / nv, (iv + 1) / nv)
        for iw in range(nw):
            elements[iv, iw] = restrict_axis(vnet, 2, iw / nw, (iw + 1) / nw)
    return elements


def bernstein(degree, q):
    q = np.asarray(q)
    values = np.stack(
        [
            math.comb(degree, i) * q**i * (1 - q) ** (degree - i)
            for i in range(degree + 1)
        ],
        axis=1,
    )
    lower = (
        np.stack(
            [
                math.comb(degree - 1, i) * q**i * (1 - q) ** (degree - 1 - i)
                for i in range(degree)
            ],
            axis=1,
        )
        if degree
        else np.ones((len(q), 1))
    )
    derivative = np.zeros_like(values)
    if degree:
        derivative[:, 0] = -degree * lower[:, 0]
        derivative[:, -1] = degree * lower[:, -1]
        for i in range(1, degree):
            derivative[:, i] = degree * (lower[:, i - 1] - lower[:, i])
    return values, derivative


def tensor_basis(degrees, uvw):
    one_d = [bernstein(p, uvw[:, axis]) for axis, p in enumerate(degrees)]
    n = np.einsum("qi,qj,qk->qijk", *(v[0] for v in one_d)).reshape(len(uvw), -1)
    dn = []
    for axis in range(3):
        factors = [v[0] for v in one_d]
        factors[axis] = one_d[axis][1]
        dn.append(np.einsum("qi,qj,qk->qijk", *factors).reshape(len(uvw), -1))
    return n, np.stack(dn, axis=-1)


def evaluate(net, degrees, uvw):
    flat = net.reshape(-1, 4)
    n, dn = tensor_basis(degrees, uvw)
    weights = flat[:, 3]
    den = n @ weights
    r = n * weights / den[:, None]
    dr = (
        dn * weights[None, :, None]
        - r[:, :, None] * np.einsum("qia,i->qa", dn, weights)[:, None, :]
    ) / den[:, None, None]
    xyz = flat[:, :3] / weights[:, None]
    x = r @ xyz
    jac = np.einsum("ic,qia->qca", xyz, dr)
    return r, dr, x, jac


class DisjointSet:
    def __init__(self, size):
        self.parent = np.arange(size)

    def find(self, i):
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return int(i)

    def join(self, a, b):
        self.parent[self.find(int(a))] = self.find(int(b))


def build_mesh(base, locations, degrees, nv, nw):
    elements = []
    lookup = {}
    for chart, net in enumerate(base):
        key = next(k for k, value in locations.items() if value == chart)
        restricted = element_nets(net, degrees, nv, nw)
        for (iv, iw), geom in restricted.items():
            ei = len(elements)
            elements.append((key, iv, iw, geom))
            lookup[key, iv, iw] = ei
    shape = tuple(p + 1 for p in degrees)
    local_count = math.prod(shape)
    ids = np.arange(len(elements) * local_count).reshape((len(elements), *shape))
    sets = DisjointSet(ids.size)

    def face_equal(net_a, net_b, axis, side_a=1, side_b=0):
        aa = np.take(net_a, -1 if side_a else 0, axis=axis)
        bb = np.take(net_b, -1 if side_b else 0, axis=axis)
        # Homogeneous nets on corresponding traces must agree after one chart gauge.
        aa = aa / aa.reshape(-1, 4)[0, 3]
        bb = bb / bb.reshape(-1, 4)[0, 3]
        if not np.allclose(aa, bb, rtol=2e-12, atol=2e-12):
            raise ValueError("Elevated/refined geometry traces do not match")

    def join_face(ea, eb, axis):
        na, nb = elements[ea][3], elements[eb][3]
        face_equal(na, nb, axis)
        for idx in np.ndindex(tuple(shape[d] for d in range(3) if d != axis)):
            one = list(idx)
            two = list(idx)
            one.insert(axis, shape[axis] - 1)
            two.insert(axis, 0)
            sets.join(ids[(ea, *one)], ids[(eb, *two)])

    # C0 joins between children within each chart.
    for key, iv, iw, _ in elements:
        ei = lookup[key, iv, iw]
        if iv + 1 < nv:
            join_face(ei, lookup[key, iv + 1, iw], 1)
        if iw + 1 < nw:
            join_face(ei, lookup[key, iv, iw + 1], 2)

    # Full-face joins on the conforming 22-chart cover.
    for (band, col, layer), ei0 in locations.items():
        neighbors = [
            ((band + 1, col, layer), 0),
            ((band, col - 1, layer), 1),
            ((band, col, layer + 1), 2),
        ]
        for other, axis in neighbors:
            if other not in locations:
                continue
            if axis == 1:
                pairs = [((nv - 1, iw), (0, iw)) for iw in range(nw)]
            elif axis == 2:
                pairs = [((iv, nw - 1), (iv, 0)) for iv in range(nv)]
            else:
                pairs = [((iv, iw), (iv, iw)) for iv in range(nv) for iw in range(nw)]
            for (iv_a, iw_a), (iv_b, iw_b) in pairs:
                a = lookup[(band, col, layer), iv_a, iw_a]
                b = lookup[other, iv_b, iw_b]
                join_face(a, b, axis)

    # Collapse the cubic chart's two polar lines to single-valued traces.
    max_band = max(b for b, _, _ in locations)
    for band, end in ((0, 0), (max_band, degrees[0])):
        for layer in range(3):
            key = (band, 1, layer)
            if key not in locations:
                continue
            for iw in range(nw):
                for kw in range(shape[2]):
                    anchor = ids[(lookup[key, 0, iw], end, 0, kw)]
                    for iv in range(nv):
                        ei = lookup[key, iv, iw]
                        for jv in range(shape[1]):
                            sets.join(anchor, ids[(ei, end, jv, kw)])

    roots = np.array([sets.find(i) for i in range(ids.size)])
    unique, labels = np.unique(roots, return_inverse=True)
    labels = labels.reshape(ids.shape)
    fixed = set()
    for ei, (key, iv, iw, _) in enumerate(elements):
        band, col, layer = key
        if col == 0 and iv == nv - 1:
            fixed.update(labels[ei, :, -1, :].ravel().tolist())
    return elements, lookup, labels, len(unique), fixed


def gauss_points(order, face_axis=None, side=0):
    x, w = leggauss(order)
    q, gw = (x + 1) / 2, w / 2
    if face_axis is None:
        xyz = np.stack(np.meshgrid(q, q, q, indexing="ij"), axis=-1).reshape(-1, 3)
        weights = np.einsum("i,j,k->ijk", gw, gw, gw).ravel()
    else:
        axes = [q, q, q]
        axis_w = [gw, gw, gw]
        axes[face_axis] = np.array([float(side)])
        axis_w[face_axis] = np.array([1.0])
        xyz = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
        weights = np.einsum("i,j,k->ijk", *axis_w).ravel()
    return xyz, weights


def element_stiffness(net, degrees, order):
    uvw, weights = gauss_points(order)
    _, dn, _, jac = evaluate(net, degrees, uvw)
    det = np.linalg.det(jac)
    if np.min(det) <= 0:
        raise ValueError("Nonpositive Jacobian in refined geometry")
    grad = np.einsum("qia,qac->qic", dn, np.linalg.inv(jac))
    b = np.zeros((len(weights), 6, grad.shape[1], 3))
    for axis in range(3):
        b[:, axis, :, axis] = grad[:, :, axis]
    for row, (a, c) in enumerate(((0, 1), (1, 2), (2, 0)), 3):
        b[:, row, :, a], b[:, row, :, c] = grad[:, :, c], grad[:, :, a]
    b = b.reshape(len(weights), 6, -1)
    local_k = np.einsum(
        "qai,ab,qbj,q->ij", b, study.ELASTIC, b, weights * det, optimize=True
    )
    volume = float(weights @ det)
    return local_k, volume


def solve_case(base, locations, degrees, nv, nw, order, return_system=False):
    elements, lookup, labels, scalar_dofs, fixed = build_mesh(
        base, locations, degrees, nv, nw
    )
    free_scalar = np.array([i for i in range(scalar_dofs) if i not in fixed], dtype=int)
    free_map = np.full(scalar_dofs, -1, dtype=int)
    free_map[free_scalar] = np.arange(len(free_scalar))
    total = 3 * len(free_scalar)
    row_parts, col_parts, data_parts = [], [], []
    rhs = np.zeros(total)
    total_volume = 0.0
    uvw, face_weights = gauss_points(order, face_axis=1, side=0)
    for ei, (key, iv, iw, net) in enumerate(elements):
        local_k, volume = element_stiffness(net, degrees, order)
        total_volume += volume
        local_scalar = labels[ei].ravel()
        local_free = free_map[local_scalar]
        local_dofs = (3 * local_free[:, None] + np.arange(3)).ravel()
        keep = local_dofs >= 0
        kept_dofs = local_dofs[keep]
        selected = local_k[np.ix_(keep, keep)]
        row_parts.append(np.repeat(kept_dofs, len(kept_dofs)))
        col_parts.append(np.tile(kept_dofs, len(kept_dofs)))
        data_parts.append(selected.ravel())

        band, col, layer = key
        if col == 2 and iv == 0:
            n, _, _, jac = evaluate(net, degrees, uvw)
            area = -np.cross(jac[:, :, 2], jac[:, :, 0])
            local_f = np.zeros((len(local_scalar), 3))
            local_f[:, 0] = n.T @ (face_weights * np.linalg.norm(area, axis=1))
            for i, root in enumerate(local_scalar):
                mapped = free_map[root]
                if mapped >= 0:
                    rhs[3 * mapped : 3 * mapped + 3] += local_f[i]

    stiffness = sparse.coo_matrix(
        (
            np.concatenate(data_parts),
            (np.concatenate(row_parts), np.concatenate(col_parts)),
        ),
        shape=(total, total),
    ).tocsr()
    stiffness.sum_duplicates()
    if return_system:
        return stiffness, rhs, free_map, elements, lookup, labels, total_volume
    displacement = spsolve(stiffness, rhs)
    residual = float(
        np.linalg.norm(stiffness @ displacement - rhs) / np.linalg.norm(rhs)
    )
    energy = float(0.5 * displacement @ (stiffness @ displacement))
    force = float(rhs[::3].sum())
    mean_ux = float(rhs @ displacement / force)
    return {
        "degrees_uvw": list(degrees),
        "subdivisions_vw": [nv, nw],
        "gauss_order": order,
        "element_count": len(elements),
        "scalar_dofs_free": len(free_scalar),
        "vector_dofs_free": total,
        "loaded_force_x_N": force,
        "loaded_end_mean_ux_mm": mean_ux,
        "strain_energy_N_mm": energy,
        "reduced_residual": residual,
        "geometry_volume_mm3": total_volume,
    }


def load_geometry(binary, rise):
    recipe = json.loads(RECIPE.read_text())
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        inp, out = root / "input.json", root / "charts.json"
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
            [str(binary.resolve()), "crossed-boolean", str(inp), "union", str(out)],
            check=True,
            capture_output=True,
        )
        model = json.loads(out.read_text())
    return fitted_geometry(model, rise)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--rise", type=int, default=3)
    args = parser.parse_args()
    base, locations = load_geometry(args.binary, args.rise)
    baseline = json.loads((HERE / "results.json").read_text())
    standard = next(
        row for row in baseline if row["rise_mm"] == args.rise and row["order"] == 14
    )
    cases = [
        ("p311", (3, 1, 1), 1, 1),
        ("p321", (3, 2, 1), 1, 1),
        ("p312", (3, 1, 2), 1, 1),
        ("p322", (3, 2, 2), 1, 1),
        ("p322_v2", (3, 2, 2), 2, 1),
        ("p322_w2", (3, 2, 2), 1, 2),
        ("p322_vw2", (3, 2, 2), 2, 2),
        ("p322_v4", (3, 2, 2), 4, 1),
        ("p322_w4", (3, 2, 2), 1, 4),
        ("p322_v4w2", (3, 2, 2), 4, 2),
        ("p322_v2w4", (3, 2, 2), 2, 4),
    ]
    result = {
        "schema": "crossed_interface_conforming_bending_study_v1",
        "case": {
            "rise_mm": args.rise,
            "geometry": "original S2 tilted crossed-tool union, extended to its conforming 22-chart matching-face control",
            "geometry_source": "tests/references/g2_cross_sketches.recipe.json via crossed-boolean union",
            "material": {"young_MPa": study.YOUNG, "poisson": study.POISSON},
            "boundary_conditions": "x=0 clamped; 1 MPa traction on x=100",
            "observable": "traction-weighted mean loaded-end ux",
            "parameter_axes": {
                "u": "shared cubic ordinate parameter along the crossed profile",
                "v": "longitudinal center-to-wing direction; p321 enriches this direction",
                "w": "physical z thickness/height direction; p312 enriches this direction",
            },
            "geometry_map": "unchanged rational (3,1,1) map; homogeneous controls are exactly degree-elevated for analysis quadrature",
            "analysis_basis": "rational Bernstein tensor basis with p=(3,pv,pw)",
            "solver": "SciPy sparse direct solve, with reduced residual reported",
            "trace": "C0 coefficient identification on every matching face, all mesh seams, and collapsed polar lines",
            "quadrature": "tensor Gauss-Legendre, recorded per row",
            "interior_continuity": "single-span elements joined C0; internal knot multiplicity equals degree",
        },
        "baseline_reference": {
            "degrees_uvw": [3, 1, 1],
            "fitted_charts": standard["fitted_charts"],
            "fitted_scalar_dofs": standard["fitted_scalar_dofs"],
            "loaded_end_mean_ux_mm": standard["fitted_mean_ux_mm"],
        },
        "cases": {},
    }
    for name, degrees, nv, nw in cases:
        # q=8 preserves the established crossed-interface study integration order.
        row = solve_case(base, locations, degrees, nv, nw, 8)
        row["relative_volume_difference"] = abs(
            row["geometry_volume_mm3"] / standard["fitted_volume_mm3"] - 1
        )
        result["cases"][name] = row
        print(
            name,
            row["vector_dofs_free"],
            row["loaded_end_mean_ux_mm"],
            row["strain_energy_N_mm"],
            row["reduced_residual"],
            flush=True,
        )
    last = solve_case(base, locations, (3, 2, 2), 2, 2, 12)
    last["relative_volume_difference"] = abs(
        last["geometry_volume_mm3"] / standard["fitted_volume_mm3"] - 1
    )
    result["quadrature_convergence"] = {
        "case": "p322_vw2",
        "order_8": result["cases"]["p322_vw2"],
        "order_12": last,
        "loaded_end_mean_ux_relative_change": abs(
            last["loaded_end_mean_ux_mm"]
            / result["cases"]["p322_vw2"]["loaded_end_mean_ux_mm"]
            - 1
        ),
        "strain_energy_relative_change": abs(
            last["strain_energy_N_mm"]
            / result["cases"]["p322_vw2"]["strain_energy_N_mm"]
            - 1
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
