"""Exact independent coefficient oracle for lateral joins of audited extrusions.

Check coplanarity of the two tangent directions of the physical graph (x,y,q)
on one side and the transverse direction on the other. Work with polynomials
in the face parameter, independently of the Rust unisolvent gradient equations.
"""

import argparse
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def add(a, b, scale=1):
    c = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        c[i] += x
    for i, x in enumerate(b):
        c[i] += scale * x
    while c and not c[-1]:
        c.pop()
    return c


def mul(a, b):
    if not a or not b:
        return []
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return add(c, [])


def bernstein(row):
    # Coefficients in ascending monomials on [0,1].
    d = len(row) - 1
    out = [F(0)] * (d + 1)
    for i, value in enumerate(row):
        for j in range(d - i + 1):
            out[i + j] += value * comb(d, i) * comb(d - i, j) * (-1) ** j
    return add(out, [])


def deriv(p):
    return [i * p[i] for i in range(1, len(p))]


def cross(a, b):
    return add(mul(a[0], b[1]), mul(a[1], b[0]), -1)


def elevated_plane(block, degrees):
    du, dv, _ = block["degrees"]
    points = [
        [list(map(F, block["controls"][2 * (i * (dv + 1) + j)])) for j in range(dv + 1)]
        for i in range(du + 1)
    ]

    def elevate(rows):
        n = len(rows) - 1
        new = [rows[0]]
        for i in range(1, n + 1):
            alpha = F(i, n + 1)
            new.append(
                [
                    [alpha * a + (1 - alpha) * b for a, b in zip(x, y)]
                    for x, y in zip(rows[i - 1], rows[i])
                ]
            )
        return new + [rows[-1]]

    for _ in range(du, degrees[0]):
        points = elevate(points)
    for _ in range(dv, degrees[1]):
        points = list(zip(*elevate(list(zip(*points)))))
    return points


def plane_jets(net, degrees, side, reverse=False):
    d = degrees["uv".index(side[0])]
    rows = net if side[0] == "u" else list(zip(*net))
    edge, near = (rows[-1], rows[-2]) if side[1] == "1" else (rows[0], rows[1])
    # Reversing either transverse direction does not change coplanarity.
    normal = [d * (b - a) for a, b in zip(edge, near)]
    if reverse:
        edge, normal = edge[::-1], normal[::-1]
    value = bernstein(edge)
    return value, deriv(value), bernstein(normal)


def lifted_derivative(jets, weight):
    h, t, n = jets
    w, wt, wn = weight
    return add(mul(t, w), mul(h, wt), -1), add(mul(n, w), mul(h, wn), -1)


def traces(block, cell, side, reverse):
    du, dv, p = cell["degrees"]
    grid = elevated_plane(block, (du, dv))

    def component(c):
        return [[pt[c] for pt in row] for row in grid]

    W = plane_jets(component(3), [du, dv], side, reverse)
    derivatives = [
        lifted_derivative(plane_jets(component(c), [du, dv], side, reverse), W)
        for c in [0, 1]
    ]
    tangent, normal = map(tuple, zip(*derivatives))
    rows = [
        {i: F(int(r[0]), int(r[1])) for i, r in row}
        for row in cell["rational_extraction"]
    ]

    def field_at(layer):
        # Only the first two transverse rows can contribute to the boundary jet.
        boundary = 0 if side[1] == "0" else (du if side[0] == "u" else dv)
        ids = set()
        for i in range(du + 1):
            for j in range(dv + 1):
                if abs((i if side[0] == "u" else j) - boundary) <= 1:
                    ids.update(rows[(i * (dv + 1) + j) * (p + 1) + layer])
        result = {}
        for g in ids:
            net = [
                [
                    grid[i][j][3] * rows[(i * (dv + 1) + j) * (p + 1) + layer].get(g, 0)
                    for j in range(dv + 1)
                ]
                for i in range(du + 1)
            ]
            jets = plane_jets(net, [du, dv], side, reverse)
            dt, dn = lifted_derivative(jets, W)
            result[g] = (jets[0], dt, dn)
        return result

    return W[0], tangent, normal, [field_at(k) for k in range(p + 1)]


def verify(package, check_retained=False):
    model = package["geometry"]
    blocks = {b["id"]: b for b in model["blocks"]}
    cells = {c["id"]: c for c in package["field"]["cells"]}
    orientations = {
        (tuple(f["first"]), tuple(f["second"])): f["orientation"]
        for f in model["topology"]["interfaces"]
    }
    joins = identities = bad = 0
    for e in package["field"]["interfaces"]:
        if e["role"] != "distinct_support_join" and not check_retained:
            continue
        a, sa = e["first"]
        b, sb = e["second"]
        if "w" in (sa[0], sb[0]):
            continue
        o = orientations[(tuple(e["first"]), tuple(e["second"]))]
        assert not o["swap"] and not o["flip_second"]
        w1, t, n1, A = traces(blocks[a], cells[a], sa, False)
        w2, _, n2, B = traces(blocks[b], cells[b], sb, o["flip_first"])
        ab, bc, ca = cross(t, n1), cross(n1, n2), cross(n2, t)
        joins += 1
        for arow, brow in zip(A, B, strict=True):
            for g in arow.keys() | brow.keys():
                a0, at, an = arow.get(g, ([], [], []))
                b0, _, bn = brow.get(g, ([], [], []))
                # Determinant of three cleared graph tangents at the common
                # physical point. Factors w^2 may differ between charts.
                coplanar = add(add(mul(bn, ab), mul(at, bc)), mul(an, ca))
                c0 = add(mul(a0, w2), mul(b0, w1), -1)
                if c0:
                    raise ValueError("nonzero trace difference")
                bad += bool(coplanar)
                identities += 1
    return dict(
        joins=joins, column_layer_identities=identities, nonzero_polynomials=bad
    )


def run(binary, work, quick=False):
    cases = {
        "cylinder": "docs/site/box-cylinder/tall/cylinder.model.json",
        "pocket": "docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json",
    }
    work.mkdir(parents=True, exist_ok=True)
    summary = {}
    for name, source in cases.items():
        model = json.loads((ROOT / source).read_text())
        lateral = [
            [f["first"][0], f["second"][0]]
            for f in model["topology"]["interfaces"]
            if f["first"][1][0] != "w"
        ]
        bucket = {}
        states = [
            ("source", 0, [], 0),
            ("inherited", 1, [], 0),
            ("first", 1, lateral[:1], 0),
            ("all", 1, lateral, 0),
        ]
        if not quick:
            states.append(("all_cubic_section", 1, lateral, 3))
        for label, c, links, p in states:
            plan = {
                "schema": "inherited_field_plan_v1",
                "axial_degree": 2,
                "axial_continuity": c,
                "restore_strip": False,
                "transverse_interfaces": links,
                "plane_degree": p,
            }
            path = work / f"{name}-{label}.model.json"
            req = path.with_suffix(".plan.json")
            req.write_text(json.dumps(plan))
            subprocess.run(
                [binary, "export-inherited-space", ROOT / source, req, path],
                check=True,
                capture_output=True,
            )
            pack = json.loads(path.read_text())
            # Use the unchanged source as a negative control on the full face set.
            results = verify(pack, check_retained=(label == "source"))
            assert (results["nonzero_polynomials"] > 0) == (label == "source"), (
                name,
                label,
                results,
            )
            info = pack["field"]
            bucket[label] = {
                **results,
                **{
                    k: info[k]
                    for k in [
                        "scalar_dofs",
                        "inherited_join_count",
                        "transverse_c1_join_count",
                        "retained_c0_join_count",
                    ]
                },
            }
            print(name, label, bucket[label], flush=True)
        summary[name] = bucket
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/solidcraft-trivariate"
    )
    parser.add_argument("--work", type=Path, default=ROOT / "build/inherited-oracle")
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    result = run(args.binary.resolve(), args.work, args.quick)
    (args.work / "jets.json").write_text(json.dumps(result, indent=2) + "\n")
