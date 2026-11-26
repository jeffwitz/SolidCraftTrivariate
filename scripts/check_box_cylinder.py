#!/usr/bin/env python3
"""Independent checks for the exact four-chart box/cylinder control case."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from fractions import Fraction
from pathlib import Path

import numpy as np


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bernstein(degree: int, x: float) -> np.ndarray:
    from math import comb
    return np.array([comb(degree, i) * x**i * (1 - x)**(degree-i)
                     for i in range(degree + 1)])


def bernstein_derivative(degree: int, x: float) -> np.ndarray:
    if degree == 0:
        return np.zeros(1)
    lower = bernstein(degree - 1, x)
    return np.array([degree * ((lower[i - 1] if i else 0.0)
                               - (lower[i] if i < degree else 0.0))
                     for i in range(degree + 1)])


def tensor_value(block: dict, uvw: tuple[float, float, float]) -> tuple[np.ndarray, np.ndarray]:
    degrees = block["degrees"]
    controls = np.asarray(block["controls"], dtype=float).reshape(
        (degrees[0] + 1, degrees[1] + 1, degrees[2] + 1, 4))
    b = [bernstein(d, x) for d, x in zip(degrees, uvw)]
    db = [bernstein_derivative(d, x) for d, x in zip(degrees, uvw)]
    h = np.einsum("i,j,k,ijkc->c", b[0], b[1], b[2], controls)
    dh = np.stack([
        np.einsum("i,j,k,ijkc->c", db[0], b[1], b[2], controls),
        np.einsum("i,j,k,ijkc->c", b[0], db[1], b[2], controls),
        np.einsum("i,j,k,ijkc->c", b[0], b[1], db[2], controls),
    ])
    xyz = h[:3] / h[3]
    jac = (dh[:, :3] - xyz[None, :] * dh[:, 3, None]) / h[3]
    return xyz, jac


def check_circle_identity(block: dict) -> None:
    """Exact degree-four Bernstein identity on the bore's v=1 trace."""
    p, q, r = block["degrees"]
    assert (p, q, r) == (2, 1, 1)
    net = [block["controls"][(i * 2 + 1) * 2] for i in range(3)]
    assert all(
        all(a[k] * b[3] == b[k] * a[3] for k in (0, 1))
        for i in range(3)
        for a, b in [(block["controls"][i * 4 + 2], block["controls"][i * 4 + 3])]
    ), "cylindrical trace varies through thickness"
    assert all(
        block["controls"][i * 4 + 2][2] == 0
        and block["controls"][i * 4 + 3][2] == 12 * block["controls"][i * 4 + 3][3]
        for i in range(3)
    ), "wall extrusion is not exactly z=0..12"
    # Power coefficients for (X-20W)^2+(Y-20W)^2-(5W)^2.
    # Product of degree-two Bernstein polynomials has degree four.
    from math import comb
    vectors = [[Fraction(c[0] - 20*c[3]), Fraction(c[1] - 20*c[3]),
                Fraction(c[3])] for c in net]

    def product_coeff(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(5)]
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i+j] += x*y*Fraction(comb(2, i)*comb(2, j), comb(4, i+j))
        return out

    x2 = product_coeff([v[0] for v in vectors], [v[0] for v in vectors])
    y2 = product_coeff([v[1] for v in vectors], [v[1] for v in vectors])
    w2 = product_coeff([v[2] for v in vectors], [v[2] for v in vectors])
    assert all(x+y-25*w == 0 for x, y, w in zip(x2, y2, w2)), \
        f"exact circle polynomial failed for {block['id']}"


def side_net(block: dict, axis: int, upper: bool) -> list[list[int]]:
    deg = block["degrees"]
    shape = tuple(x + 1 for x in deg)
    net = np.asarray(block["controls"], dtype=object).reshape((*shape, 4))
    return net.take(deg[axis] if upper else 0, axis=axis).reshape(-1, 4).tolist()


def proportional_net(a: list[list[int]], b: list[list[int]]) -> bool:
    if len(a) != len(b):
        return False
    global_scale: Fraction | None = None
    for x, y in zip(a, b):
        ratios = [Fraction(i, j) for i, j in zip(x, y) if i and j]
        if any((i == 0) != (j == 0) for i, j in zip(x, y)) or not ratios:
            return False
        if len(set(ratios)) != 1 or (global_scale is not None and ratios[0] != global_scale):
            return False
        scale = ratios[0]
        global_scale = scale
        if any(Fraction(i) != scale * j for i, j in zip(x, y)):
            return False
    return True


def verify(model: dict) -> dict:
    assert model["schema"] in ("orthogonal_sketch_features_v1", "box_cylinder_boolean_v1")
    assert model["partition_plan"]["method"] == "four_sector_local_sweeps"
    blocks = {b["id"]: b for b in model["blocks"]}
    assert len(blocks) == 4
    assert len(model["topology"]["interfaces"]) == 4
    assert len(model["topology"]["exterior"]) == 16
    walls = [f for f in model["design_faces"] if f["feature"].endswith("/wall")]
    assert len(walls) == 4 and all(f["side"] == "v1" for f in walls)
    for block in blocks.values():
        check_circle_identity(block)
    for interface in model["topology"]["interfaces"]:
        first, s1 = interface["first"]
        second, s2 = interface["second"]
        axis1, up1 = "uvw".index(s1[0]), s1[1] == "1"
        axis2, up2 = "uvw".index(s2[0]), s2[1] == "1"
        a, b = side_net(blocks[first], axis1, up1), side_net(blocks[second], axis2, up2)
        orient = interface["orientation"]
        assert not orient["swap"] and not orient["flip_first"] and not orient["flip_second"]
        assert proportional_net(a, b), f"interface trace mismatch: {first} {s1} / {second} {s2}"

    nodes, weights = np.polynomial.legendre.leggauss(20)
    nodes, weights = (nodes + 1) / 2, weights / 2
    volume = 0.0
    min_det = math.inf
    for block in blocks.values():
        for i, u in enumerate(nodes):
            for j, v in enumerate(nodes):
                for k, w in enumerate(nodes):
                    _, jac = tensor_value(block, (u, v, w))
                    det = float(np.linalg.det(jac))
                    assert det > 0
                    min_det = min(min_det, det)
                    volume += weights[i]*weights[j]*weights[k]*det
    exact = 19200 - 300*math.pi
    rel = abs(volume-exact)/exact
    assert rel < 2e-13, (volume, exact, rel)
    return {
        "status": "valid",
        "model_sha256": None,
        "independent_checks": {
            "chart_count": 4,
            "interfaces_exact_homogeneous": 4,
            "exterior_faces": 16,
            "cylindrical_wall_traces": 4,
            "circle_identity": "exact_degree_four_bernstein",
            "circle_center_mm": [20, 20],
            "circle_radius_mm": 5,
            "through_thickness_mm": [0, 12],
            "volume_q20_mm3": volume,
            "volume_exact_expression_mm3": "19200 - 300*pi",
            "volume_exact_numeric_mm3": exact,
            "volume_relative_error": rel,
            "minimum_sampled_jacobian_mm3": min_det,
        },
    }


def compare_map_sets(first: dict, second: dict) -> dict:
    """Compare the two constructor routes after discarding only chart labels."""
    def maps(model: dict) -> list[tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]]:
        return sorted((tuple(block["degrees"],), tuple(map(tuple, block["controls"])))
                      for block in model["blocks"])
    a, b = maps(first), maps(second)
    return {"exact_control_maps_equal_ignoring_ids": a == b,
            "classic_charts": len(a), "boolean_charts": len(b)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cli", type=Path, default=Path("target/release/solidcraft-trivariate"))
    ap.add_argument("--recipe", type=Path, default=Path("tests/references/box_cylinder_classic.recipe.json"))
    ap.add_argument("--model", type=Path, default=Path("build/box-cylinder-classic.model.json"))
    ap.add_argument("--output", type=Path, default=Path("docs/results/box-cylinder/classic-independent-audit.json"))
    ap.add_argument("--boolean-recipe", type=Path, default=Path("tests/references/box_minus_cylinder.recipe.json"))
    ap.add_argument("--boolean-model", type=Path, default=Path("build/box-minus-cylinder.model.json"))
    args = ap.parse_args()
    args.model.parent.mkdir(parents=True, exist_ok=True)
    args.boolean_model.parent.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([str(args.cli), "build", str(args.recipe), str(args.model)], check=True)
    subprocess.run([str(args.cli), "audit", str(args.model)], check=True)
    model = json.loads(args.model.read_text())
    report = verify(model)
    report["model_sha256"] = sha(args.model)
    report["recipe_sha256"] = sha(args.recipe)
    subprocess.run([str(args.cli), "build", str(args.boolean_recipe), str(args.boolean_model)], check=True)
    subprocess.run([str(args.cli), "audit", str(args.boolean_model)], check=True)
    boolean_model = json.loads(args.boolean_model.read_text())
    boolean_report = verify(boolean_model)
    route_comparison = compare_map_sets(model, boolean_model)
    assert route_comparison["classic_charts"] == route_comparison["boolean_charts"] == 4
    assert route_comparison["exact_control_maps_equal_ignoring_ids"], route_comparison
    report["two_operand_route"] = {
        "model_sha256": sha(args.boolean_model),
        "recipe_sha256": sha(args.boolean_recipe),
        "independent_checks": boolean_report["independent_checks"],
        "route_comparison": route_comparison,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
