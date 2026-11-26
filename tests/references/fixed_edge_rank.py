#!/usr/bin/env python3
"""Audit exposed block edges that force a zero Jacobian when smoothed.

The exact conclusion is topological/differential: two exposed cube sides share
one block edge. If their boundary tangent planes are made equal along that
edge, the three-dimensional block Jacobian vanishes there. This script reads
the published homogeneous Bezier maps and independently reports the exposed
edge incidences and local physical crease/smooth diagnostics.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


AXES = "uvw"
SIDES = [(axis, upper) for axis in range(3) for upper in (False, True)]


def bernstein(n: int, x: float) -> list[float]:
    if n == 0:
        return [1.0]
    return [math.comb(n, i) * x**i * (1.0 - x) ** (n - i) for i in range(n + 1)]


def bernstein_derivative(n: int, x: float) -> list[float]:
    if n == 0:
        return [0.0]
    lower = bernstein(n - 1, x)
    return [n * ((lower[i - 1] if i else 0.0) - (lower[i] if i < n else 0.0))
            for i in range(n + 1)]


def eval_jacobian(block: dict, q: tuple[float, float, float]) -> list[tuple[float, float, float]]:
    p = block["degrees"]
    controls = block["controls"]
    bases = [bernstein(p[d], q[d]) for d in range(3)]
    derivs = [bernstein_derivative(p[d], q[d]) for d in range(3)]
    h = [0.0] * 4
    dh = [[0.0] * 4 for _ in range(3)]
    for i, j, k in itertools.product(range(p[0] + 1), range(p[1] + 1), range(p[2] + 1)):
        c = controls[(i * (p[1] + 1) + j) * (p[2] + 1) + k]
        indices = (i, j, k)
        factors = (bases[0][i], bases[1][j], bases[2][k])
        weight = factors[0] * factors[1] * factors[2]
        for a in range(4):
            h[a] += c[a] * weight
            for d in range(3):
                f = derivs[d][(i, j, k)[d]]
                for e in range(3):
                    if e != d:
                        f *= bases[e][indices[e]]
                dh[d][a] += c[a] * f
    w = h[3]
    return [tuple((dh[d][a] * w - h[a] * dh[d][3]) / (w * w) for a in range(3)) for d in range(3)]


def fbernstein(n: int, x: Fraction) -> list[Fraction]:
    return [Fraction(math.comb(n, i)) * x**i * (1 - x) ** (n - i) for i in range(n + 1)]


def fbernstein_derivative(n: int, x: Fraction) -> list[Fraction]:
    if n == 0:
        return [Fraction(0)]
    lower = fbernstein(n - 1, x)
    return [n * ((lower[i - 1] if i else 0) - (lower[i] if i < n else 0)) for i in range(n + 1)]


def eval_jacobian_exact(block: dict, q: tuple[Fraction, Fraction, Fraction]):
    """Return exact Cartesian derivatives for integer homogeneous controls."""
    p = block["degrees"]
    controls = block["controls"]
    bases = [fbernstein(p[d], q[d]) for d in range(3)]
    derivs = [fbernstein_derivative(p[d], q[d]) for d in range(3)]
    h = [Fraction(0) for _ in range(4)]
    dh = [[Fraction(0) for _ in range(4)] for _ in range(3)]
    for i, j, k in itertools.product(range(p[0] + 1), range(p[1] + 1), range(p[2] + 1)):
        indices = (i, j, k)
        c = controls[(i * (p[1] + 1) + j) * (p[2] + 1) + k]
        factors = (bases[0][i], bases[1][j], bases[2][k])
        weight = factors[0] * factors[1] * factors[2]
        for a in range(4):
            h[a] += c[a] * weight
            for d in range(3):
                f = derivs[d][indices[d]]
                for e in range(3):
                    if e != d:
                        f *= bases[e][indices[e]]
                dh[d][a] += c[a] * f
    w = h[3]
    return [tuple((dh[d][a] * w - h[a] * dh[d][3]) / (w * w) for a in range(3)) for d in range(3)]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def norm(a):
    return math.sqrt(sum(x * x for x in a))


def face_normal(jac, fixed_axis, free_axis, fixed_side):
    # Side normal is J_free1 x J_free2; use the block's two other axes.
    other = [a for a in range(3) if a != fixed_axis]
    return cross(jac[other[0]], jac[other[1]])


def scan(path: Path) -> dict:
    model = json.loads(path.read_text())
    blocks = {b["id"]: b for b in model["blocks"]}
    exposed = {(b, side): None for b, side in model["topology"]["exterior"]}
    for row in model.get("design_faces", []):
        exposed[(row["block"], row["side"])] = row["feature"]
    edges = []
    for bid, block in blocks.items():
        for free in range(3):
            fixed = [a for a in range(3) if a != free]
            for bits in itertools.product((False, True), repeat=2):
                sides = [(fixed[z], bits[z]) for z in range(2)]
                face_keys = [(bid, f"{AXES[a]}{int(upper)}") for a, upper in sides]
                if all(key in exposed for key in face_keys):
                    records = []
                    angles = []
                    exact_sharp = False
                    for t in (Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1)):
                        q = [Fraction(1, 2), Fraction(1, 2), Fraction(1, 2)]
                        q[free] = t
                        for (a, upper) in sides:
                            q[a] = Fraction(int(upper))
                        jac_exact = eval_jacobian_exact(block, tuple(q))
                        jac = [tuple(float(x) for x in row) for row in jac_exact]
                        normals = [face_normal(jac, a, free, upper) for a, upper in sides]
                        n0, n1 = (norm(n) for n in normals)
                        cosine = abs(sum(normals[0][i] * normals[1][i] for i in range(3)) / (n0 * n1))
                        angle = math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
                        angles.append(angle)
                        exact_normals = [cross(jac_exact[[x for x in range(3) if x != a][0]], jac_exact[[x for x in range(3) if x != a][1]]) for a, _ in sides]
                        exact_sharp = exact_sharp or any(x != 0 for x in cross(exact_normals[0], exact_normals[1]))
                        # J determinant in standard parameter orientation.
                        det = sum(jac[0][i] * cross(jac[1], jac[2])[i] for i in range(3))
                        records.append(abs(det))
                    sharp = exact_sharp
                    length_axis = AXES[free]
                    edges.append({
                        "block": bid,
                        "edge": f"{length_axis} with " + ",".join(f"{AXES[a]}{int(u)}" for a, u in sides),
                        "faces": [{"side": key[1], "feature": exposed[key]} for key in face_keys],
                        "classification": "physical_crease" if sharp else "smooth_seam_or_tangent_edge",
                        "normal_angle_deg_samples": [round(x, 8) for x in angles],
                        "abs_jacobian_samples": [round(x, 12) for x in records],
                    })
    sharp = [e for e in edges if e["classification"] == "physical_crease"]
    smooth = [e for e in edges if e["classification"] != "physical_crease"]
    return {
        "input": str(path),
        "schema": model.get("schema"),
        "block_count": len(blocks),
        "interface_count": len(model["topology"]["interfaces"]),
        "exterior_face_count": len(exposed),
        "two_exposed_face_edges": len(edges),
        "sampled_physical_creases": len(sharp),
        "sampled_smooth_or_tangent_edges": len(smooth),
        "sample_method": "exact Fraction rational Bezier first derivatives at edge parameters 0, 1/4, 1/2, 3/4, 1; angles and determinant magnitudes are rounded for display",
        "exact_argument": "At any regular point of a block edge, the edge tangent F_r is nonzero. Equal tangent planes of its two exposed coordinate faces force the remaining two transverse derivatives into one common plane containing F_r; hence det(DF)=0. Simultaneous vertex blends cannot change this edge-local implication.",
        "edges": edges,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("models", nargs="*", type=Path, default=[
        Path("docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json"),
        Path("docs/results/s3-face-pocket.model.json"),
    ])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"schema": "fixed_exposed_edge_rank_v1", "models": [scan(p) for p in args.models]}
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
