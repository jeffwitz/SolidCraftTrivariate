#!/usr/bin/env python3
"""Independent mapped anti-plane elasticity sweep on the bilinear blind pocket."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import qr
from scipy.sparse import coo_matrix, bmat, csc_matrix
from scipy.sparse.linalg import eigsh, spsolve

ROOT = Path(__file__).resolve().parents[2]
RECIPE = ROOT / "tests/references/blind_bilinear_pocket.recipe.json"
P_DEGREES = (3, 4)
REFINEMENTS = (1, 2, 4, 8)
MU = 1.0
LX, LY = 60.0, 40.0
POCKET_X, POCKET_Y = (20.0, 40.0), (15.0, 25.0)
POCKET_WALL_HEIGHT = 6.0


@dataclass
class Patch:
    name: str
    control: np.ndarray  # P00, P10, P11, P01 in physical x,y
    thickness: float
    degree: int
    n: int
    continuity: int
    knots_u: np.ndarray
    knots_v: np.ndarray
    nb_u: int
    nb_v: int
    offset: int = 0


def knots_for(p: int, n: int, continuity: int) -> np.ndarray:
    interior_mult = p if continuity == 0 else p - 1
    return np.asarray([0.0] * (p + 1)
                      + [k / n for k in range(1, n) for _ in range(interior_mult)]
                      + [1.0] * (p + 1))


def spline_pair(knots: np.ndarray, p: int):
    count = len(knots) - p - 1
    base = BSpline(knots, np.eye(count), p, axis=0, extrapolate=False)
    return base, base.derivative(1)


def geometry(patch: Patch, u: float, v: float):
    p00, p10, p11, p01 = patch.control
    point = ((1-u)*(1-v)*p00 + u*(1-v)*p10 + u*v*p11 + (1-u)*v*p01)
    du = (1-v)*(p10-p00) + v*(p11-p01)
    dv = (1-u)*(p01-p00) + u*(p11-p10)
    jac = np.column_stack((du, dv))
    return point, jac


def edge_uv(side: str, t: float) -> tuple[float, float]:
    return ((0.0, t) if side == "u0" else
            (1.0, t) if side == "u1" else
            (t, 0.0) if side == "v0" else
            (t, 1.0))


def edge_axis(side: str) -> int:
    return 0 if side[0] == "u" else 1


def edge_basis_indices(patch: Patch, side: str) -> list[int]:
    if side[0] == "u":
        fixed = 0 if side[1] == "0" else patch.nb_u - 1
        return [patch.offset + fixed * patch.nb_v + j for j in range(patch.nb_v)]
    fixed = 0 if side[1] == "0" else patch.nb_v - 1
    return [patch.offset + i * patch.nb_v + fixed for i in range(patch.nb_u)]


def local_dof(patch: Patch, i: int, j: int) -> int:
    return patch.offset + i * patch.nb_v + j


class UnionFind:
    def __init__(self, count: int):
        self.parent = list(range(count))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def join(self, a: int, b: int):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def decode_patch_model(model: dict, p: int, n: int, continuity: int):
    blocks = {b["id"]: b for b in model["blocks"]}
    # The section complex is the five lower charts. The four upper charts repeat
    # the four outer footprints and contribute the additional 6 mm thickness.
    names = [b["id"] for b in model["blocks"] if b["id"].endswith("_below")]
    if len(names) != 5:
        raise ValueError(f"expected five lower quadrilateral charts, got {len(names)}")
    patches: dict[str, Patch] = {}
    cursor = 0
    for name in names:
        block = blocks[name]
        if block["degrees"] != [1, 1, 1] or any(c[3] != 1 for c in block["controls"]):
            raise ValueError(f"{name}: expected an exact unit-weight bilinear extrusion")
        controls = block["controls"]
        pts = np.asarray([controls[((i*2+j)*2)][:2]
                          for i, j in ((0, 0), (1, 0), (1, 1), (0, 1))], dtype=float)
        thickness = 14.0 if "footprint" in name else 20.0
        ku = knots_for(p, n, continuity)
        kv = knots_for(p, n, continuity)
        nb_u, nb_v = len(ku)-p-1, len(kv)-p-1
        patch = Patch(name, pts, thickness, p, n, continuity, ku, kv, nb_u, nb_v, cursor)
        # det(dX/du,dX/dv) is affine on a bilinear quadrilateral, so corner
        # positivity is an exact closed-domain check up to input arithmetic.
        dets = [np.linalg.det(geometry(patch, u, v)[1])
                for u, v in ((0, 0), (1, 0), (1, 1), (0, 1))]
        if min(dets) <= 0:
            raise ValueError(f"{name}: nonpositive bilinear Jacobian {dets}")
        patches[name] = patch
        cursor += nb_u*nb_v

    # Keep only the eight horizontal interfaces in the lower five-chart complex.
    interfaces = []
    for face in model["topology"]["interfaces"]:
        first, second = face["first"], face["second"]
        if (first[0] in patches and second[0] in patches
                and first[1][0] in "uv" and second[1][0] in "uv"):
            interfaces.append((first[0], first[1], second[0], second[1]))
    if len(interfaces) != 8:
        raise ValueError(f"expected 8 in-plane chart interfaces, got {len(interfaces)}")
    return patches, interfaces, cursor


def edge_orientation(pa: Patch, sa: str, pb: Patch, sb: str) -> int:
    a0 = geometry(pa, *edge_uv(sa, 0.0))[0]
    a1 = geometry(pa, *edge_uv(sa, 1.0))[0]
    b0 = geometry(pb, *edge_uv(sb, 0.0))[0]
    b1 = geometry(pb, *edge_uv(sb, 1.0))[0]
    tol = 1e-10
    if np.linalg.norm(a0-b0) < tol and np.linalg.norm(a1-b1) < tol:
        return 1
    if np.linalg.norm(a0-b1) < tol and np.linalg.norm(a1-b0) < tol:
        return -1
    raise ValueError(f"unmatched physical interface endpoints: {pa.name}/{sa}, {pb.name}/{sb}")


def build_global_map(model, patches, interfaces, raw_count):
    uf = UnionFind(raw_count)
    for an, ass, bn, bs in interfaces:
        a, b = patches[an], patches[bn]
        if a.degree != b.degree or a.n != b.n:
            raise ValueError("interface spaces are not aligned")
        orient = edge_orientation(a, ass, b, bs)
        ia, ib = edge_basis_indices(a, ass), edge_basis_indices(b, bs)
        if len(ia) != len(ib):
            raise ValueError("interface trace basis counts differ")
        for k, da in enumerate(ia):
            kb = k if orient > 0 else len(ib)-1-k
            uf.join(da, ib[kb])
    roots = sorted({uf.find(i) for i in range(raw_count)})
    root_index = {r: k for k, r in enumerate(roots)}
    mapping = np.asarray([root_index[uf.find(i)] for i in range(raw_count)], dtype=int)

    fixed = set()
    # The physical outer contour, not chart seams or the pocket wall, is clamped.
    outer = []
    for face in model["design_faces"]:
        name = face["block"]
        side = face["side"]
        if name in patches and side[0] in "uv" and face["feature"].startswith("outer/"):
            outer.append((name, side))
    for name, side in outer:
        fixed.update(mapping[d] for d in edge_basis_indices(patches[name], side))
    if not fixed:
        raise ValueError("no outer stock Dirichlet trace was identified")
    return mapping, len(roots), fixed


def active_span(knots, p, midpoint):
    k = int(np.searchsorted(knots, midpoint, side="right") - 1)
    return list(range(k-p, k+1))


def basis_at(spline, derivative, x, ids):
    vals = np.asarray(spline(float(x)), dtype=float)[ids]
    ders = np.asarray(derivative(float(x)), dtype=float)[ids]
    return vals, ders


def chart_point_and_gradient(patch: Patch, u: float, v: float,
                             basis_u, deriv_u, basis_v, deriv_v):
    point, jac = geometry(patch, u, v)
    det = float(np.linalg.det(jac))
    if det <= 0:
        raise ValueError(f"{patch.name}: nonpositive sampled Jacobian")
    Nu = np.asarray(basis_u(u), dtype=float)
    Nv = np.asarray(basis_v(v), dtype=float)
    dNu = np.asarray(deriv_u(u), dtype=float)
    dNv = np.asarray(deriv_v(v), dtype=float)
    param = np.stack((np.outer(dNu, Nv).ravel(), np.outer(Nu, dNv).ravel()), axis=1)
    physical = np.linalg.solve(jac.T, param.T).T
    shape = np.outer(Nu, Nv).ravel()
    return point, jac, det, shape, physical


def assemble(model, patches, interfaces, raw_count, qorder):
    mapping, nglobal, fixed = build_global_map(model, patches, interfaces, raw_count)
    free = np.asarray([i for i in range(nglobal) if i not in fixed], dtype=int)
    free_index = {g: i for i, g in enumerate(free)}
    qx, qw = np.polynomial.legendre.leggauss(qorder)
    qx, qw = (qx+1)/2, qw/2
    rows, cols, vals = [], [], []
    rhs = np.zeros(nglobal)
    kx, ky = math.pi/LX, math.pi/LY
    for patch in patches.values():
        bu, dbu = spline_pair(patch.knots_u, patch.degree)
        bv, dbv = spline_pair(patch.knots_v, patch.degree)
        unique_u = np.unique(patch.knots_u)
        unique_v = np.unique(patch.knots_v)
        for ua, ub in zip(unique_u, unique_u[1:]):
            if ub <= ua:
                continue
            iu = active_span(patch.knots_u, patch.degree, (ua+ub)/2)
            gx = ua + (ub-ua)*qx
            U = np.asarray(bu(gx))[:, iu]
            dU = np.asarray(dbu(gx))[:, iu]
            for va, vb in zip(unique_v, unique_v[1:]):
                if vb <= va:
                    continue
                iv = active_span(patch.knots_v, patch.degree, (va+vb)/2)
                gy = va + (vb-va)*qx
                V = np.asarray(bv(gy))[:, iv]
                dV = np.asarray(dbv(gy))[:, iv]
                local_ids = np.asarray([mapping[local_dof(patch, i, j)]
                                        for i in iu for j in iv])
                Ke = np.zeros((len(local_ids), len(local_ids)))
                fe = np.zeros(len(local_ids))
                for a in range(qorder):
                    u = gx[a]
                    for b in range(qorder):
                        v = gy[b]
                        point, jac = geometry(patch, u, v)
                        det = float(np.linalg.det(jac))
                        if det <= 0:
                            raise ValueError(f"{patch.name}: invalid quadrature Jacobian")
                        par = np.stack((np.outer(dU[a], V[b]).ravel(),
                                        np.outer(U[a], dV[b]).ravel()), axis=1)
                        grad = np.linalg.solve(jac.T, par.T).T
                        shape = np.outer(U[a], V[b]).ravel()
                        x, y = point
                        w = math.sin(kx*x)*math.sin(ky*y)
                        body = (kx*kx+ky*ky)*w
                        weight = qw[a]*qw[b]*(ub-ua)*(vb-va)*det*patch.thickness*MU
                        Ke += weight*(grad @ grad.T)
                        fe += weight*shape*body
                rows.extend(np.repeat(local_ids, len(local_ids)))
                cols.extend(np.tile(local_ids, len(local_ids)))
                vals.extend(Ke.ravel())
                np.add.at(rhs, local_ids, fe)

    # Manufactured natural shear traction on the four real pocket walls.
    footprint = next(p for p in patches.values() if "footprint" in p.name)
    wall_links = [link for link in interfaces if link[0] == footprint.name or link[2] == footprint.name]
    if len(wall_links) != 4:
        raise ValueError(f"expected four footprint wall interfaces, got {len(wall_links)}")
    line_q, line_w = np.polynomial.legendre.leggauss(qorder)
    line_q, line_w = (line_q+1)/2, line_w/2
    for an, ass, bn, bs in wall_links:
        # The ring chart's outward normal points from material into the pocket.
        ring, side = (patches[an], ass) if an != footprint.name else (patches[bn], bs)
        bu, dbu = spline_pair(ring.knots_u, ring.degree)
        bv, dbv = spline_pair(ring.knots_v, ring.degree)
        is_u = side[0] == "u"
        trace_ids = edge_basis_indices(ring, side)
        tangential_knots = ring.knots_v if is_u else ring.knots_u
        line_breaks = np.unique(tangential_knots)
        line_points = np.concatenate([ta+(tb-ta)*line_q
                                      for ta, tb in zip(line_breaks, line_breaks[1:])])
        line_weights = np.concatenate([(tb-ta)*line_w
                                       for ta, tb in zip(line_breaks, line_breaks[1:])])
        for q, qwgt in zip(line_points, line_weights):
            u, v = edge_uv(side, float(q))
            point, jac = geometry(ring, u, v)
            tangent = jac[:, 1] if is_u else jac[:, 0]
            ds = float(np.linalg.norm(tangent))
            normal = np.array([tangent[1], -tangent[0]])/ds
            mid = geometry(ring, *edge_uv(side, 0.5))[0]
            center = geometry(ring, 0.5, 0.5)[0]
            if np.dot(normal, mid-center) < 0:
                normal = -normal
            gx, gy = kx*math.cos(kx*point[0])*math.sin(ky*point[1]), \
                     ky*math.sin(kx*point[0])*math.cos(ky*point[1])
            traction = POCKET_WALL_HEIGHT*MU*np.dot([gx, gy], normal)
            Nu = np.asarray(bu(u), dtype=float)
            Nv = np.asarray(bv(v), dtype=float)
            # The trace coefficients are in the side's canonical edge order.
            if side[0] == "u":
                fixed_u = 0 if side[1] == "0" else ring.nb_u-1
                trace = Nu[fixed_u]*Nv
            else:
                fixed_v = 0 if side[1] == "0" else ring.nb_v-1
                trace = Nu*Nv[fixed_v]
            for j, value in enumerate(trace):
                raw = trace_ids[j]
                rhs[mapping[raw]] += qwgt*ds*traction*value

    K = coo_matrix((vals, (rows, cols)), shape=(nglobal, nglobal)).tocsc()
    K.sum_duplicates()
    keep = np.asarray([free_index[g] for g in range(nglobal) if g not in fixed], dtype=int)
    Kff = K[free, :][:, free].tocsc()
    bff = rhs[free]
    return K, rhs, free, fixed, Kff, bff, mapping, keep


def c1_constraints(patches, interfaces, mapping, nglobal, fixed, p, sample_count=None,
                   normalize_rows=True):
    # On each tangential spline span, p+3 samples determine the numerator
    # polynomial after clearing the two linear bilinear-map Jacobians.
    rows, cols, data = [], [], []
    row = 0
    for an, ass, bn, bs in interfaces:
        a, b = patches[an], patches[bn]
        orient = edge_orientation(a, ass, b, bs)
        tangential_knots = a.knots_v if ass[0] == "u" else a.knots_u
        breaks = np.unique(tangential_knots)
        per_span_count = sample_count or p+3
        sample = np.concatenate([
            np.linspace(left, right, per_span_count)
            for left, right in zip(breaks, breaks[1:])
        ])
        ba_u, dba_u = spline_pair(a.knots_u, p)
        ba_v, dba_v = spline_pair(a.knots_v, p)
        bb_u, dbb_u = spline_pair(b.knots_u, p)
        bb_v, dbb_v = spline_pair(b.knots_v, p)
        u0, v0 = edge_uv(ass, 0.5)
        _, ja = geometry(a, u0, v0)
        tangent = ja[:, 1] if ass[0] == "u" else ja[:, 0]
        normal = np.array([tangent[1], -tangent[0]])
        mid = geometry(a, *edge_uv(ass, 0.5))[0]
        center = geometry(a, 0.5, 0.5)[0]
        normal /= np.linalg.norm(normal)
        if np.dot(normal, mid-center) < 0:
            normal = -normal
        for t in sample:
            ta = float(t)
            tb = ta if orient > 0 else 1-ta
            ua, va = edge_uv(ass, ta)
            ub, vb = edge_uv(bs, tb)
            _, ja, _, _, grad_a = chart_point_and_gradient(a, ua, va, ba_u, dba_u, ba_v, dba_v)
            _, jb, _, _, grad_b = chart_point_and_gradient(b, ub, vb, bb_u, dbb_u, bb_v, dbb_v)
            wa = grad_a @ normal
            wb = grad_b @ normal
            ca = np.flatnonzero(np.abs(wa) > 1e-13)
            cb = np.flatnonzero(np.abs(wb) > 1e-13)
            coeffs = {}
            for flat in ca:
                i, j = divmod(int(flat), a.nb_v)
                g = mapping[local_dof(a, i, j)]
                coeffs[g] = coeffs.get(g, 0.0) + wa[flat]
            for flat in cb:
                i, j = divmod(int(flat), b.nb_v)
                g = mapping[local_dof(b, i, j)]
                coeffs[g] = coeffs.get(g, 0.0) - wb[flat]
            norm = math.sqrt(sum(v*v for v in coeffs.values()))
            if norm < 1e-12:
                continue
            scale = norm if normalize_rows else 1.0
            for g, value in coeffs.items():
                if g not in fixed and abs(value/scale) > 1e-14:
                    rows.append(row); cols.append(g); data.append(value/scale)
            row += 1
    C = coo_matrix((data, (rows, cols)), shape=(row, nglobal)).tocsr()
    free = np.asarray([g for g in range(nglobal) if g not in fixed], dtype=int)
    C = C[:, free]
    if C.shape[0] == 0:
        return C
    # Remove dependent collocation rows by rank-revealing pivoted QR.
    dense_t = C.toarray().T
    _, R, piv = qr(dense_t, mode="economic", pivoting=True)
    diag = np.abs(np.diag(R))
    rank = int(np.sum(diag > max(diag[0] if len(diag) else 0.0, 1.0)*1e-10))
    return C[np.sort(piv[:rank]), :]


def exact_integrals(patches, qorder):
    q, w = np.polynomial.legendre.leggauss(qorder)
    q, w = (q+1)/2, w/2
    kx, ky = math.pi/LX, math.pi/LY
    mass = energy = 0.0
    for patch in patches.values():
        for u, wu in zip(q, w):
            for v, wv in zip(q, w):
                point, jac = geometry(patch, float(u), float(v))
                det = np.linalg.det(jac)
                x, y = point
                val = math.sin(kx*x)*math.sin(ky*y)
                grad = np.array([kx*math.cos(kx*x)*math.sin(ky*y),
                                 ky*math.sin(kx*x)*math.cos(ky*y)])
                wt = wu*wv*det*patch.thickness
                mass += wt*val*val
                energy += wt*float(grad@grad)
    return mass, energy


def manufactured_weak_identity(patches, interfaces, qorder):
    """Check exact-solution work against body force and pocket-wall traction."""
    q, w = np.polynomial.legendre.leggauss(qorder)
    q, w = (q+1)/2, w/2
    kx, ky = math.pi/LX, math.pi/LY
    body_work = exact_energy = wall_work = 0.0
    for patch in patches.values():
        for u, wu in zip(q, w):
            for v, wv in zip(q, w):
                point, jac = geometry(patch, float(u), float(v))
                det = float(np.linalg.det(jac))
                x, y = point
                value = math.sin(kx*x)*math.sin(ky*y)
                grad = np.array([kx*math.cos(kx*x)*math.sin(ky*y),
                                 ky*math.sin(kx*x)*math.cos(ky*y)])
                weight = wu*wv*det*patch.thickness*MU
                body_work += weight*(kx*kx+ky*ky)*value*value
                exact_energy += weight*float(grad@grad)
    footprint = next(p for p in patches.values() if "footprint" in p.name)
    walls = [f for f in interfaces if f[0] == footprint.name or f[2] == footprint.name]
    for an, ass, bn, bs in walls:
        ring, side = (patches[an], ass) if an != footprint.name else (patches[bn], bs)
        is_u = side[0] == "u"
        tangential = np.unique(ring.knots_v if is_u else ring.knots_u)
        for ta, tb in zip(tangential, tangential[1:]):
            for qi, wi in zip(q, w):
                t = ta+(tb-ta)*qi
                point, jac = geometry(ring, *edge_uv(side, float(t)))
                tangent = jac[:, 1] if is_u else jac[:, 0]
                ds = float(np.linalg.norm(tangent))
                normal = np.array([tangent[1], -tangent[0]])/ds
                midpoint = geometry(ring, *edge_uv(side, 0.5))[0]
                center = geometry(ring, 0.5, 0.5)[0]
                if np.dot(normal, midpoint-center) < 0:
                    normal = -normal
                x, y = point
                value = math.sin(kx*x)*math.sin(ky*y)
                grad = np.array([kx*math.cos(kx*x)*math.sin(ky*y),
                                 ky*math.sin(kx*x)*math.cos(ky*y)])
                wall_work += (wi*(tb-ta)*ds*POCKET_WALL_HEIGHT*MU
                              * float(grad@normal)*value)
    return {"body_work": body_work, "pocket_wall_work": wall_work,
            "exact_strain_norm_squared": exact_energy,
            "relative_weak_identity_defect": abs(body_work+wall_work-exact_energy)/exact_energy}


def errors(patches, mapping, coeff, fixed, p, qorder):
    q, w = np.polynomial.legendre.leggauss(qorder)
    q, w = (q+1)/2, w/2
    kx, ky = math.pi/LX, math.pi/LY
    l2 = en = 0.0
    for patch in patches.values():
        bu, dbu = spline_pair(patch.knots_u, p)
        bv, dbv = spline_pair(patch.knots_v, p)
        spans_u, spans_v = np.unique(patch.knots_u), np.unique(patch.knots_v)
        for ua, ub in zip(spans_u, spans_u[1:]):
            iu = active_span(patch.knots_u, p, (ua+ub)/2)
            for va, vb in zip(spans_v, spans_v[1:]):
                iv = active_span(patch.knots_v, p, (va+vb)/2)
                ids = [mapping[local_dof(patch, i, j)] for i in iu for j in iv]
                c = coeff[ids]
                for a, qx in enumerate(q):
                    u = ua+(ub-ua)*qx
                    Nu, dNu = basis_at(bu, dbu, u, iu)
                    for b, qy in enumerate(q):
                        v = va+(vb-va)*qy
                        Nv, dNv = basis_at(bv, dbv, v, iv)
                        point, jac = geometry(patch, float(u), float(v))
                        det = np.linalg.det(jac)
                        shp = np.outer(Nu, Nv).ravel()
                        par = np.stack((np.outer(dNu, Nv).ravel(),
                                        np.outer(Nu, dNv).ravel()), axis=1)
                        grad_h = np.linalg.solve(jac.T, par.T).T
                        uh = float(shp@c)
                        gh = c@grad_h
                        x, y = point
                        ue = math.sin(kx*x)*math.sin(ky*y)
                        ge = np.array([kx*math.cos(kx*x)*math.sin(ky*y),
                                       ky*math.sin(kx*x)*math.cos(ky*y)])
                        wt = w[a]*w[b]*(ub-ua)*(vb-va)*det*patch.thickness
                        l2 += wt*(uh-ue)**2
                        en += wt*float((gh-ge)@(gh-ge))
    return math.sqrt(l2), math.sqrt(en)


def solve_case(model, p, n, continuity, qorder):
    patches, interfaces, raw = decode_patch_model(model, p, n, continuity)
    K, rhs, free, fixed, Kff, bff, mapping, _ = assemble(
        model, patches, interfaces, raw, qorder)
    multipliers = np.zeros(0)
    C = None
    if continuity == 0:
        coeff_free = spsolve(Kff, bff)
        c_rank = 0
        constraint_defect = 0.0
        independent_jump = 0.0
        multipliers = np.zeros(0)
        C = None
    else:
        C = c1_constraints(patches, interfaces, mapping, K.shape[0], fixed, p)
        saddle = bmat([[Kff, C.T], [C, None]], format="csc")
        solved = spsolve(saddle, np.r_[bff, np.zeros(C.shape[0])])
        coeff_free = solved[:len(free)]
        multipliers = solved[len(free):]
        c_rank = C.shape[0]
        constraint_defect = float(np.max(np.abs(C@coeff_free))) if C.shape[0] else 0.0
        Ccheck = c1_constraints(patches, interfaces, mapping, K.shape[0], fixed, p,
                                sample_count=2*p+7, normalize_rows=False)
        independent_jump = float(np.max(np.abs(Ccheck@coeff_free))) if Ccheck.shape[0] else 0.0
    coeff = np.zeros(K.shape[0])
    coeff[free] = coeff_free
    residual = Kff@coeff_free-bff
    if C is not None and C.shape[0]:
        residual += C.T@multipliers
    relres = float(np.linalg.norm(residual)/max(np.linalg.norm(bff), 1e-30))
    smallest_eigenvalue = float(eigsh(Kff, k=1, which="SA", return_eigenvectors=False,
                                      tol=1e-7, maxiter=max(1000, 20*Kff.shape[0]))[0])
    l2, en = errors(patches, mapping, coeff, fixed, p, qorder)
    mass_exact, energy_exact = exact_integrals(patches, qorder+8)
    # Also compare exact analytic norms at a distinct higher Gauss order.
    mass_hi, energy_hi = exact_integrals(patches, qorder+12)
    norm_change = max(abs(mass_exact-mass_hi)/mass_hi,
                      abs(energy_exact-energy_hi)/energy_hi)
    return {
        "p": p, "n_per_parent_direction": n,
        "continuity": "C1" if continuity else "C0",
        "scalar_free_dofs_before_C1": int(len(free)),
        "scalar_free_dofs": int(len(free)-c_rank),
        "constraint_rank": int(c_rank), "quadrature_order": qorder,
        "relative_energy_error": en/math.sqrt(energy_exact),
        "relative_l2_error": l2/math.sqrt(mass_exact),
        "strain_energy": 0.5*float(coeff@(K@coeff)),
        "relative_free_residual": relres,
        "max_c1_constraint_defect": constraint_defect,
        "independent_dense_c1_jump_defect": independent_jump,
        "stiffness_smallest_eigenvalue": smallest_eigenvalue,
        "exact_norm_q_change": norm_change,
        "matrix_nnz": int(Kff.nnz),
    }


def load_model(binary: Path, work: Path):
    work.mkdir(parents=True, exist_ok=True)
    path = work / "bilinear-pocket.model.json"
    completed = subprocess.run([str(binary), "build", str(RECIPE), str(path)], check=True,
                               cwd=ROOT, capture_output=True, text=True)
    model = json.loads(path.read_text())
    model["_build_audit"] = json.loads(completed.stdout)
    if model.get("schema") != "sketch_layer_partition_v1":
        raise ValueError(f"wrong source schema: {model.get('schema')}")
    if len(model["blocks"]) != 9 or any(b["degrees"] != [1, 1, 1] for b in model["blocks"]):
        raise ValueError("CLI did not produce the expected nine trilinear blocks")
    return model, hashlib.sha256(path.read_bytes()).hexdigest()


def qualifies(base, q_delta):
    return (base["relative_free_residual"] <= 1e-10
            and base["max_c1_constraint_defect"] <= 1e-9
            and base["independent_dense_c1_jump_defect"] <= 1e-9
            and base["stiffness_smallest_eigenvalue"] > 0
            and q_delta <= 1e-9
            and base["exact_norm_q_change"] <= 1e-12)


def run(binary: Path, output: Path, qorders=(16, 24), resume=False,
        initial_gate=None):
    model, model_sha = load_model(binary.resolve(), ROOT/"build/pocket-lateral-c1")
    existing = {}
    if resume and output.exists():
        previous = json.loads(output.read_text())
        if previous["source_model_sha256"] != model_sha:
            raise ValueError("resume output source-model SHA differs")
        existing = {(r["p"], r["n_per_parent_direction"], r["continuity"]): r
                    for r in previous["results"]}
    records = []
    for p in P_DEGREES:
        for n in REFINEMENTS:
            for c in (0, 1):
                key = (p, n, "C1" if c else "C0")
                old = existing.get(key)
                if old and old.get("status") == "valid":
                    records.append(old)
                    continue
                attempts = list(old.get("quadrature_attempts", [])) if old else []
                if old and not attempts:
                    attempts.append({"orders": old.get("quadrature_pair", list(qorders)),
                                     "result": old, "status": old.get("status")})
                pair = (tuple(qorders) if not old else
                        (old.get("quadrature_pair", list(qorders))[1],
                         old.get("quadrature_pair", list(qorders))[1]+8))
                final = old
                while True:
                    base, check = [solve_case(model, p, n, c, q) for q in pair]
                    q_delta = max(abs(base[k]-check[k])/max(abs(check[k]), 1e-30)
                                  for k in ("relative_energy_error", "relative_l2_error", "strain_energy"))
                    status = "valid" if qualifies(base, q_delta) else "failed"
                    attempts.append({"orders": list(pair), "base": base,
                                    "replay": check, "relative_observable_change": q_delta,
                                    "status": status})
                    base["quadrature_pair"] = list(pair)
                    base["quadrature_observable_relative_change"] = q_delta
                    base[f"q{pair[1]}_replay"] = check
                    base["quadrature_attempts"] = attempts
                    base["status"] = status
                    final = base
                    print(p, n, base["continuity"], base["scalar_free_dofs"],
                          base["relative_energy_error"], base["relative_l2_error"],
                          base["relative_free_residual"], q_delta, status, flush=True)
                    if status == "valid" or not resume or pair[1] >= 32:
                        break
                    pair = (pair[1], pair[1]+8)
                records.append(final)
                write_results(output, model, model_sha, records, existing, initial_gate)
    # Observed h-rates, separately for energy and L2, do not presuppose a rate.
    for c in ("C0", "C1"):
        for p in P_DEGREES:
            selected = sorted((r for r in records if r["continuity"] == c and r["p"] == p),
                              key=lambda r: r["n_per_parent_direction"])
            for a, b in zip(selected, selected[1:]):
                rate_e = math.log(a["relative_energy_error"]/b["relative_energy_error"], 2)
                rate_l = math.log(a["relative_l2_error"]/b["relative_l2_error"], 2)
                b["rate_from_previous_n"] = {"energy": rate_e, "l2": rate_l}
    write_results(output, model, model_sha, records, existing, initial_gate)
    failed = [r for r in records if r["status"] != "valid"]
    if failed:
        raise RuntimeError(f"{len(failed)} sweep cases failed; see {output}")
    return records


def write_results(output, model, model_sha, records, existing, initial_gate):
    # Include already-resolved cases when writing an incremental resume.
    by_key = {(r["p"], r["n_per_parent_direction"], r["continuity"]): r
              for r in [*existing.values(), *records]}
    ordered = [by_key[(p, n, c)] for p in P_DEGREES for n in REFINEMENTS
               for c in ("C0", "C1") if (p, n, c) in by_key]
    control_patches, control_interfaces, _ = decode_patch_model(model, 3, 1, 0)
    weak_identity = {f"q{q}": manufactured_weak_identity(control_patches, control_interfaces, q)
                     for q in (24, 48)}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({
        "schema": "blind_bilinear_pocket_lateral_sweep_v1",
        "source_recipe": str(RECIPE.relative_to(ROOT)),
        "source_model_sha256": model_sha,
        "source_family": model["schema"],
        "source_geometry_audit": model["_build_audit"],
        "section_chart_count": 5,
        "lower_footprint_chart_count": 5,
        "upper_outer_chart_count": 4,
        "in_plane_interface_count": 8,
        "field_scope": "restricted anti-plane Ritz subproblem of 3D linear elasticity; displacement constant through physical z",
        "postrun_manufactured_weak_identity": weak_identity,
        "initial_q12_q16_gate": initial_gate,
        "results": ordered,
    }, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT/"target/release/solidcraft-trivariate")
    parser.add_argument("--output", type=Path, default=ROOT/"build/pocket-lateral-c1/replay.json")
    parser.add_argument("--work", type=Path, default=ROOT/"build/pocket-lateral-c1")
    parser.add_argument("--resume", action="store_true",
                        help="continue failed cases in an existing result using the next quadrature pair")
    args = parser.parse_args()
    start = time.monotonic()
    model, sha = load_model(args.binary.resolve(), args.work)
    print(f"source={sha} blocks={len(model['blocks'])}", flush=True)
    # The first n=1 state for each p/continuity is the preregistered gate.
    gates = []
    for p in P_DEGREES:
        for c in (0, 1):
            gate_results = []
            for q in (12, 16):
                first = solve_case(model, p, 1, c, q)
                gate_results.append(first)
                print("gate", p, first["continuity"], q, first["scalar_free_dofs"],
                      first["relative_energy_error"], first["relative_l2_error"],
                      first["relative_free_residual"], flush=True)
            q_delta = max(abs(gate_results[0][k]-gate_results[1][k])
                          / max(abs(gate_results[1][k]), 1e-30)
                          for k in ("relative_energy_error", "relative_l2_error", "strain_energy"))
            gates.append({"p": p, "continuity": "C1" if c else "C0",
                          "q12": gate_results[0], "q16": gate_results[1],
                          "relative_observable_change": q_delta,
                          "status_at_original_1e-9_threshold": "valid" if qualifies(gate_results[1], q_delta) else "failed"})
    run(args.binary, args.output, resume=args.resume, initial_gate=gates)
    print(f"elapsed_s={time.monotonic()-start:.3f}", flush=True)
