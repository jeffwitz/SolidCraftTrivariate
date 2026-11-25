#!/usr/bin/env python3
"""Independent source-curve quadrature and optional OCC oracle for M1 recipes."""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path as FilePath

DEFAULT_RECIPE = FilePath(__file__).with_name("nurbs_multiloop.recipe.json")


def basis(i, p, t, knots):
    if p == 0:
        return float(knots[i] <= t < knots[i + 1])
    a, b = knots[i + p] - knots[i], knots[i + p + 1] - knots[i + 1]
    left = 0.0 if a == 0 else (t - knots[i]) / a * basis(i, p - 1, t, knots)
    right = 0.0 if b == 0 else (knots[i + p + 1] - t) / b * basis(i + 1, p - 1, t, knots)
    return left + right


def curve(spec, t):
    p, u, controls = spec["degree"], spec["knots"], spec["controls"]
    if t == u[len(controls)]:
        t = math.nextafter(float(t), -math.inf)
    ns = [basis(i, p, t, u) for i in range(len(controls))]
    d_ns = []
    for i in range(len(controls)):
        a, b = u[i + p] - u[i], u[i + p + 1] - u[i + 1]
        lo = basis(i, p - 1, t, u) if p else 0.0
        hi = basis(i + 1, p - 1, t, u) if p else 0.0
        d_ns.append(0.0 if p == 0 else p * ((lo / a if a else 0.0) - (hi / b if b else 0.0)))
    h = [sum(n * c[k] for n, c in zip(ns, controls)) for k in range(3)]
    dh = [sum(n * c[k] for n, c in zip(d_ns, controls)) for k in range(3)]
    pt = (h[0] / h[2], h[1] / h[2])
    der = tuple((dh[k] * h[2] - h[k] * dh[2]) / h[2]**2 for k in range(2))
    return pt, der


def gauss_legendre(n):
    """Compute Gauss-Legendre nodes and weights by Newton iteration."""
    xs, ws = [], []
    for i in range(1, n + 1):
        z = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, z
            for j in range(2, n + 1):
                p0, p1 = p1, ((2*j-1)*z*p1-(j-1)*p0)/j
            pn = z*p1 if n == 1 else p1
            pnm1 = 1.0 if n == 1 else p0
            dp = n * (z * pn - pnm1) / (z*z - 1.0)
            dz = pn / dp
            z -= dz
            if abs(dz) < 2e-16:
                break
        xs.append(z)
        ws.append(2.0 / ((1-z*z) * dp*dp))
    return list(zip(xs, ws))


def integrate(f, a, b, order):
    return (b-a)*0.5*sum(w*f((a+b)*0.5 + (b-a)*0.5*x) for x,w in gauss_legendre(order))


def loop_area(loop, order):
    low, high = loop["lower"], loop["upper"]
    knots, p = low["knots"], low["degree"]
    spans = [(knots[i], knots[i+1]) for i in range(p, len(low["controls"])) if knots[i] < knots[i+1]]
    total = 0.0
    for a, b in spans:
        def f(t):
            (xl, yl), (dx, _) = curve(low, t)
            (xh, yh), _ = curve(high, t)
            if abs(xl-xh) > 1e-11:
                raise ValueError(f"unsynchronized profile x at t={t}: {xl} vs {xh}")
            return (yh-yl)*dx
        total += integrate(f, a, b, order)
    return total


def numeric_rank(rows, tol=1e-10):
    a = [list(map(float, row)) for row in rows]
    m, n, rank = len(a), len(a[0]), 0
    for col in range(n):
        pivot = max(range(rank, m), key=lambda r: abs(a[r][col]), default=rank)
        if rank >= m or abs(a[pivot][col]) <= tol:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        q = a[rank][col]
        a[rank] = [x/q for x in a[rank]]
        for r in range(m):
            if r != rank:
                q = a[r][col]
                a[r] = [x-q*y for x,y in zip(a[r],a[rank])]
        rank += 1
        if rank == m:
            break
    return rank


def conic_rank_diagnostics(loop):
    rows = []
    for key in ("lower", "upper"):
        spec = loop[key]
        degree = spec["degree"]
        knots = [Fraction(k) for k in spec["knots"]]
        controls = [[Fraction(c) for c in row] for row in spec["controls"]]
        a, b = knots[degree], knots[degree+1]
        def exact_curve(t):
            ns = []
            for i in range(len(controls)):
                # Exact Cox/de Boor basis for this first non-empty knot span.
                vals = [Fraction(int(knots[j] <= t < knots[j+1])) for j in range(len(knots)-1)]
                for p in range(1,degree+1):
                    nxt=[]
                    for j in range(len(vals)-1):
                        dl=knots[j+p]-knots[j]; dr=knots[j+p+1]-knots[j+1]
                        nxt.append((Fraction(0) if not dl else (t-knots[j])/dl*vals[j]) +
                                    (Fraction(0) if not dr else (knots[j+p+1]-t)/dr*vals[j+1]))
                    vals=nxt
                ns.append(vals[i])
            h=[sum(n*c[k] for n,c in zip(ns,controls)) for k in range(3)]
            return h[0]/h[2], h[1]/h[2]
        points = [exact_curve(a+(b-a)*Fraction(q,7)) for q in range(7)]
        matrix = [[x*x, x*y, y*y, x, y, Fraction(1)] for x,y in points]
        rank = fraction_rank(matrix)
        rows.append({"curve": key, "samples": len(points), "first_span_exact_fraction_conic_design_rank": rank,
                     "interpretation": "exact Fraction rank 6 means the first span is not contained in a conic; therefore the full curve is not one conic"})
    return rows


def fraction_rank(rows):
    a=[list(row) for row in rows]
    m,n=len(a),len(a[0]); rank=0
    for col in range(n):
        pivot=next((r for r in range(rank,m) if a[r][col]),None)
        if pivot is None: continue
        a[rank],a[pivot]=a[pivot],a[rank]
        q=a[rank][col]; a[rank]=[x/q for x in a[rank]]
        for r in range(m):
            if r != rank:
                q=a[r][col]
                if q: a[r]=[x-q*y for x,y in zip(a[r],a[rank])]
        rank+=1
        if rank==m: break
    return rank


def occ_reference(recipe):
    """Construct exact source B-spline wire cutters and boolean them from a stock box."""
    saved_globals = {name: globals().get(name) for name in ("argparse", "json", "math", "FilePath")}
    try:
        import FreeCAD as App
        import Part
    except ImportError:
        return {"available": False, "reason": "FreeCAD Part modules are not importable in this Python runtime"}
    globals().update(saved_globals)

    def source_edge(spec, z):
        unique, mults = [], []
        for knot in spec["knots"]:
            if unique and unique[-1] == knot:
                mults[-1] += 1
            else:
                unique.append(knot); mults.append(1)
        poles = [App.Vector(x/w, y/w, z) for x,y,w in spec["controls"]]
        weights = [w for _,_,w in spec["controls"]]
        c = Part.BSplineCurve()
        c.buildFromPolesMultsKnots(poles, mults, unique, False, spec["degree"], weights)
        return c.toShape()

    xmin,xmax = recipe["bounds_xy_mm"][0]
    ymin,ymax = recipe["bounds_xy_mm"][1]
    h = recipe["height_mm"]
    solid = Part.makeBox(xmax-xmin, ymax-ymin, h, App.Vector(xmin,ymin,0))
    cutters = []
    for loop in recipe["loops"]:
        z0,z1 = loop["z_interval_mm"]
        low, high = source_edge(loop["lower"],z0), source_edge(loop["upper"],z0)
        rev_high = Part.Edge(high.Curve, high.LastParameter, high.FirstParameter)
        l0,l1=low.Vertexes[0].Point,low.Vertexes[-1].Point
        h0,h1=high.Vertexes[0].Point,high.Vertexes[-1].Point
        wire=Part.Wire([low, Part.makeLine(l1,h1), rev_high, Part.makeLine(h0,l0)])
        cutter=Part.Face(wire).extrude(App.Vector(0,0,z1-z0))
        cutters.append({"id":loop["id"], "volume_mm3":cutter.Volume, "area_q32_times_depth_mm3":loop_area(loop,32)*(z1-z0), "valid":cutter.isValid()})
        solid=solid.cut(cutter)
    return {"available":True, "valid":solid.isValid(), "solids":len(solid.Solids), "faces":len(solid.Faces), "volume_mm3":solid.Volume, "cutters":cutters}


def audit(recipe):
    rows = []
    removed = 0.0
    for loop in recipe["loops"]:
        area16, area32 = loop_area(loop, 16), loop_area(loop, 32)
        z0, z1 = loop["z_interval_mm"]
        vol = area32 * (z1-z0)
        removed += vol
        rows.append({"id": loop["id"], "area_q16_mm2": area16, "area_q32_mm2": area32,
                     "area_delta_mm2": area32-area16, "z_interval_mm": [z0,z1], "removed_volume_q32_mm3": vol,
                     "nonconic_samples":conic_rank_diagnostics(loop)})
    stock = ((recipe["bounds_xy_mm"][0][1]-recipe["bounds_xy_mm"][0][0]) *
             (recipe["bounds_xy_mm"][1][1]-recipe["bounds_xy_mm"][1][0]) * recipe["height_mm"])
    return {"oracle":"independent_cox_de_boor_gauss_legendre_v1", "source_revision":"81cb61d", "method":"homogeneous Cox/de Boor basis evaluation; rational quotient derivative; per-knot-span Gauss-Legendre integration", "order_comparison":[16,32], "stock_volume_mm3":stock, "removed_volume_mm3":removed,
            "remaining_volume_mm3":stock-removed, "loops":rows,
            "occ_reference_exploratory_not_qualifying":occ_reference(recipe),
            "counterexamples":[{"path":"tests/counterexamples/nurbs_multiloop_stationary.json", "case":"stationary-profile-counterexample", "reason":"second profile x'(2)=0 at an interior knot because physical x controls are [55,62,70,70,70,80]", "expected":"reject input before construction", "disposition":"preserved; regular mixed fixture changes the physical x controls to [55,62,67,70,74,80]"}],
            "qualification_limits":["OCC cutter volumes reported here are exploratory and not a qualification oracle; use adaptive OCC mass after STEP/BREP export for candidate comparisons", "No candidate CLI/STEP was available during this measurement"],
            "checks":{"positive_loop_areas":all(r["area_q32_mm2"]>0 for r in rows), "removed_volume_below_stock":removed<stock,
                       "profiles_synchronized_x":True, "all_first_span_exact_conic_ranks_six":all(c["first_span_exact_fraction_conic_design_rank"]==6 for r in rows for c in r["nonconic_samples"])}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recipe", nargs="?", type=FilePath, default=DEFAULT_RECIPE)
    ap.add_argument("--out", type=FilePath)
    args = ap.parse_args()
    recipe = json.loads(args.recipe.read_text())
    result = audit(recipe)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
