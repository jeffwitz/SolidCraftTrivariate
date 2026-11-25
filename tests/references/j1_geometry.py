#!/usr/bin/env python3
"""Independent geometric measurements for the frozen J1 CLI corpus.

Models are consumed from the native CLI JSON format. ``--baseline`` asks the
CLI to build every frozen recipe and records refusals as data; it never fills
missing geometry from the Python constructor.
"""

import argparse
import itertools
import json
import math
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "tests/references/j1_geometry_corpus.json"
ANALYTIC_VOLUME = {
    "rectangle": 3200.0,
    "disk": 128.0 * math.pi,
    "annulus": 96.0 * math.pi,
    "pierced": 3200.0 - 128.0 * math.pi,
    "corner": 8.0 * (196.0 - 9.0 * math.pi),
    "pocket": 3200.0 - 64.0 * math.pi,
    "bracket": 34080.0 - 616.0 * math.pi,
}
# Each entry gives (interior interfaces, Euler characteristic of the complete boundary).
TOPOLOGY = {
    "rectangle": (0, 2), "disk": (8, 2), "annulus": (4, 0),
    "pierced": (4, 0), "corner": (1, 2), "pocket": (24, 2), "bracket": (55, -2),
}
GRID = (0.0, 0.21132486540518713, 0.5, 0.7886751345948129, 1.0)


def basis(degree, x):
    return tuple(math.comb(degree, i) * x**i * (1.0 - x) ** (degree - i)
                 for i in range(degree + 1))


def dbasis(degree, x):
    if degree == 1:
        return (-1.0, 1.0)
    low = basis(degree - 1, x)
    return tuple(degree * ((low[i - 1] if i else 0.0)
                           - (low[i] if i < degree else 0.0))
                 for i in range(degree + 1))


def evaluate(block, xyz):
    degrees, controls = block["degrees"], block["controls"]
    b = [basis(d, x) for d, x in zip(degrees, xyz)]
    d = [dbasis(degree, x) for degree, x in zip(degrees, xyz)]
    h = [0.0] * 4
    dh = [[0.0] * 4 for _ in range(3)]
    ix = 0
    for ijk in itertools.product(*(range(n + 1) for n in degrees)):
        p = controls[ix]
        ix += 1
        factor = math.prod(b[a][ijk[a]] for a in range(3))
        for c in range(4):
            h[c] += factor * p[c]
        for axis in range(3):
            term = d[axis][ijk[axis]]
            for other in range(3):
                if other != axis:
                    term *= b[other][ijk[other]]
            for c in range(4):
                dh[axis][c] += term * p[c]
    w = h[3]
    point = tuple(c / w for c in h[:3])
    jac = tuple(tuple((row[c] - point[c] * row[3]) / w for c in range(3))
                for row in dh)
    return point, jac


def determinant(j):
    a, b, c = j
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def gauss(order):
    nodes, weights = [], []
    for k in range(1, order + 1):
        x = math.cos(math.pi * (k - 0.25) / (order + 0.5))
        for _ in range(40):
            p0, p1 = 1.0, x
            for n in range(2, order + 1):
                p0, p1 = p1, ((2*n-1)*x*p1 - (n-1)*p0)/n
            if order == 1:
                p0, p1 = 1.0, x
            dp = order * (x*p1 - p0) / (x*x - 1.0)
            step = p1 / dp
            x -= step
            if abs(step) < 2e-16:
                break
        p0, p1 = 1.0, x
        for n in range(2, order + 1):
            p0, p1 = p1, ((2*n-1)*x*p1 - (n-1)*p0)/n
        if order == 1:
            p0, p1 = 1.0, x
        dp = order * (x*p1 - p0) / (x*x - 1.0)
        nodes.append((x + 1.0) / 2.0)
        weights.append(1.0 / ((1.0 - x*x) * dp*dp))
    return nodes, weights


def volume(blocks, order):
    nodes, weights = gauss(order)
    total = 0.0
    for block in blocks:
        for i, u in enumerate(nodes):
            for j, v in enumerate(nodes):
                for k, w in enumerate(nodes):
                    _, jac = evaluate(block, (u, v, w))
                    total += weights[i] * weights[j] * weights[k] * determinant(jac)
    return total


def face_samples(block, axis, side, subdivisions=2):
    varying = [i for i in range(3) if i != axis]
    samples = []
    grid = [i / subdivisions for i in range(subdivisions + 1)]
    for a, b in itertools.product(grid, repeat=2):
        q = [0.0, 0.0, 0.0]
        q[axis] = float(side)
        q[varying[0]], q[varying[1]] = a, b
        point, _ = evaluate(block, q)
        samples.append(tuple(round(x, 10) for x in point))
    return samples


def topology(blocks):
    faces = {}
    for bi, block in enumerate(blocks):
        for axis in range(3):
            for side in (0, 1):
                key = tuple(sorted(face_samples(block, axis, side)))
                faces.setdefault(key, []).append((bi, axis, side))
    shared = sum(len(uses) == 2 for uses in faces.values())
    exterior = sum(len(uses) == 1 for uses in faces.values())
    overused = sum(len(uses) > 2 for uses in faces.values())

    # Triangulate a shared 4x4 sample grid on each exterior face and compute
    # V-E+F after coordinate welding. This is a sampled boundary Euler check.
    vertices, edges, triangles = set(), set(), 0
    n = 4
    for uses in faces.values():
        if len(uses) != 1:
            continue
        bi, axis, side = uses[0]
        varying = [i for i in range(3) if i != axis]
        ids = {}
        for i in range(n + 1):
            for j in range(n + 1):
                q = [0.0, 0.0, 0.0]
                q[axis] = float(side)
                q[varying[0]], q[varying[1]] = i/n, j/n
                point, _ = evaluate(blocks[bi], q)
                ids[i, j] = tuple(round(x, 9) for x in point)
                vertices.add(ids[i, j])
        for i in range(n):
            for j in range(n):
                a, b, c, d = ids[i,j], ids[i+1,j], ids[i,j+1], ids[i+1,j+1]
                for x, y in ((a,b), (a,c), (b,d), (c,d), (a,d)):
                    edges.add(tuple(sorted((x,y))))
                triangles += 2
    return {"interfaces": shared, "exterior_faces": exterior,
            "nonmanifold_faces": overused,
            "boundary_euler_sampled": len(vertices) - len(edges) + triangles}


def circle_error(case, blocks):
    tests = []
    for block in blocks:
        name = block.get("id", "")
        center_radius = None
        edge = None
        if name.startswith(("disk_sector_", "pocket_disk_sector_")):
            center_radius, edge = ((0.0, 0.0, 4.0), (1, 0.0))
        elif name.startswith("annulus_sector_"):
            # Both radial boundaries of each annular sector are circular.
            for v, radius in ((0.0, 4.0), (1.0, 2.0)):
                for u in GRID:
                    p, _ = evaluate(block, (u, v, 0.5))
                    tests.append(abs(math.hypot(p[0], p[1]) - radius))
            continue
        elif name.startswith("pierced_sector_") or name.startswith("pocket_outer_sector_"):
            center_radius, edge = ((0.0, 0.0, 4.0), (1, 1.0))
        elif name.startswith("corner_") or name.startswith("elbow_"):
            center_radius, edge = ((14.0, 14.0, 6.0), (1, 1.0))
        elif name.startswith("hole_x_sector_"):
            center_radius, edge = ((40.0, 4.0, 4.0), (1, 1.0))
        elif name.startswith("hole_y_sector_"):
            center_radius, edge = ((4.0, 32.0, 4.0), (1, 1.0))
        if center_radius is not None:
            cx, cy, radius = center_radius
            axis, parameter = edge
            for u in GRID:
                q = [0.5, 0.5, 0.5]
                q[axis] = parameter
                q[0] = u
                p, _ = evaluate(block, q)
                # Hole axes differ by family; evaluate radius in the transverse plane.
                if name.startswith("hole_x_sector_"):
                    actual = math.hypot(p[0] - cx, p[2] - 20.0)
                elif name.startswith("hole_y_sector_"):
                    actual = math.hypot(p[1] - cy, p[2] - 20.0)
                else:
                    actual = math.hypot(p[0] - cx, p[1] - cy)
                tests.append(abs(actual - radius))
    return max(tests) if tests else None


def load_models(args, manifest):
    models, refusals = {}, {}
    if args.models:
        for case in manifest["cases"]:
            path = Path(args.models) / f"{case['id']}.json"
            if path.exists():
                models[case["id"]] = json.loads(path.read_text())
            else:
                refusals[case["id"]] = f"missing model JSON: {path}"
        return models, refusals
    if not args.baseline:
        raise SystemExit("provide --models DIR or --baseline")
    with tempfile.TemporaryDirectory(prefix="j1-baseline-") as tmp:
        tmp = Path(tmp)
        for case in manifest["cases"]:
            recipe = tmp / f"{case['id']}-recipe.json"
            target = tmp / f"{case['id']}.json"
            recipe.write_text(json.dumps(case["recipe"]))
            result = subprocess.run([args.cli, "build", str(recipe), str(target)],
                                    text=True, capture_output=True, check=False)
            try:
                status = json.loads(result.stdout)
            except json.JSONDecodeError:
                refusals[case["id"]] = f"CLI output was not JSON (exit {result.returncode})"
                continue
            if result.returncode == 0 and target.exists():
                models[case["id"]] = json.loads(target.read_text())
            else:
                refusals[case["id"]] = status.get("reason", f"exit {result.returncode}")
    return models, refusals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", default="target/debug/solidcraft-trivariate")
    parser.add_argument("--models", help="directory containing <case-id>.json CLI exports")
    parser.add_argument("--baseline", action="store_true",
                        help="ask the CLI to build each frozen manifest recipe")
    parser.add_argument("--require-qualified", action="store_true",
                        help="fail if any required model is missing or outside the frozen thresholds")
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text())
    models, refusals = load_models(args, manifest)
    report = {"manifest_version": manifest["version"], "frozen_base": manifest["frozen_base"],
              "volume_reference_mm3": ANALYTIC_VOLUME, "cases": {}}
    for case in manifest["cases"]:
        cid = case["id"]
        if cid not in models:
            report["cases"][cid] = {"status": "refused", "reason": refusals.get(cid, "no model")}
            continue
        model = models[cid]
        blocks = model.get("blocks", [])
        expected_count = case["expected_blocks"]
        orders = {str(order): volume(blocks, order) for order in (4, 8, 12, 16)}
        topo = topology(blocks)
        radial_error = circle_error(cid, blocks)
        target_volume = ANALYTIC_VOLUME[cid]
        interfaces, euler = TOPOLOGY[cid]
        thresholds = manifest["thresholds"]
        qualified = (model.get("schema") == case["recipe"]["family"]
                     and len(blocks) == expected_count
                     and topo["interfaces"] == interfaces
                     and topo["exterior_faces"] == 6 * expected_count - 2 * interfaces
                     and topo["nonmanifold_faces"] == 0 and topo["boundary_euler_sampled"] == euler
                     and abs(orders["16"] - target_volume) <= thresholds["volume_relative"] * target_volume
                     and (cid == "rectangle" or
                          radial_error is not None and radial_error <= thresholds["radius_error_mm"]))
        report["cases"][cid] = {
            "qualified_by_oracle": qualified,
            "status": "measured", "schema": model.get("schema"),
            "block_count": len(blocks), "expected_block_count": expected_count,
            "block_count_matches": len(blocks) == expected_count,
            "volume_by_gauss_order_mm3": orders,
            "volume_analytic_mm3": ANALYTIC_VOLUME[cid],
            "volume_error_order16_mm3": abs(orders["16"] - ANALYTIC_VOLUME[cid]),
            "quadrature_delta_12_16_mm3": abs(orders["16"] - orders["12"]),
            "max_circle_or_fillet_error_mm": radial_error,
            **topo,
        }
    report["coverage"] = {"measured": sum(c["status"] == "measured" for c in report["cases"].values()),
                          "required": len(manifest["cases"]),
                          "refused": sum(c["status"] == "refused" for c in report["cases"].values())}
    qualified = sum(c.get("qualified_by_oracle", False) for c in report["cases"].values())
    report["coverage"]["qualified_by_oracle"] = qualified
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.require_qualified and qualified != len(manifest["cases"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
