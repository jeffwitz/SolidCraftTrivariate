"""Independent HEX8 witness on rational trivariate experiment models.

Run from the repository root with one or more model JSON paths. The fixed
levels 1, 2, 4 share nodes by exact Cartesian position. Test local validity
on the whole reference cube, both before and after VTK binary64 rounding.
The check of declared interface facets and exterior facets addresses mesh
conformity. This diagnostic does not prove global intersection freedom.
"""

import hashlib
import itertools
import json
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


def nodes(block, n):
    degrees = block["degrees"]
    bases = [
        [
            [math.comb(d, i) * k**i * (n - k) ** (d - i) for i in range(d + 1)]
            for k in range(n + 1)
        ]
        for d in degrees
    ]
    for ijk in itertools.product(range(n + 1), repeat=3):
        rows = [base[k] for k, base in zip(ijk, bases)]
        h = [0] * 4
        for raw, abw in zip(block["controls"], itertools.product(*rows)):
            factor = math.prod(abw)
            for i in range(4):
                h[i] += raw[i] * factor
        yield ijk, tuple(F(x, h[3]) for x in h[:3])


cube = list(itertools.product(range(2), repeat=3))


def det(a, b, c):
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        + a[1] * (b[2] * c[0] - b[0] * c[2])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def bound(vertices):
    den = math.lcm(*(v.denominator for p in vertices for v in p))
    p = [[v.numerator * (den // v.denominator) for v in xyz] for xyz in vertices]
    edges = [
        [
            [p[i + s][k] - p[i][k] for k in range(3)]
            for i, t in enumerate(cube)
            if not t[d]
        ]
        for d, s in enumerate([4, 2, 1])
    ]
    # Each derivative is a convex combination of four axis-aligned reference edges.
    return F(min(det(a, b, c) for a, b, c in itertools.product(*edges)), den**3)


def run(path, n):
    model = json.loads(path.read_text())
    uniq = {}
    xyz = []
    rounded = []
    conn = []
    sides = Counter()
    outer_faces = defaultdict(set)
    bounds = []
    vtk_bounds = []
    for block in model["blocks"]:
        grid = {}
        for pos, pt in nodes(block, n):
            if pt not in uniq:
                uniq[pt] = len(xyz)
                xyz.append(pt)
                rounded.append(tuple(F.from_float(float(v)) for v in pt))
            grid[pos] = uniq[pt]
        for cell in itertools.product(range(n), repeat=3):
            ids = [grid[tuple(i + j for i, j in zip(cell, delta))] for delta in cube]
            assert len(set(ids)) == 8
            conn.append([ids[i] for i in [0, 4, 6, 2, 1, 5, 7, 3]])
            bounds.append(bound([xyz[i] for i in ids]))
            vtk_bounds.append(bound([rounded[i] for i in ids]))
            for d in range(3):
                for u in range(2):
                    key = tuple(sorted(ids[i] for i, t in enumerate(cube) if t[d] == u))
                    sides[key] += 1
                    if cell[d] == (n - 1 if u else 0):
                        outer_faces[block["id"], "uvw"[d] + str(u)].add(key)
    topology = model["topology"]
    for interface in topology["interfaces"]:
        assert (
            outer_faces[tuple(interface["first"])]
            == outer_faces[tuple(interface["second"])]
        )
    exposed = set().union(*(outer_faces[tuple(s)] for s in topology["exterior"]))
    assert (
        exposed == {f for f, count in sides.items() if count == 1}
        and max(sides.values()) == 2
    )
    assert len(set(rounded)) == len(xyz), "serialization merged distinct nodes"
    out = Path("build/hex-pilot") / (path.name.split(".")[0] + "-n" + str(n) + ".vtk")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as stream:
        stream.write(
            "# vtk DataFile Version 3.0\nTrilinear hexahedra from trivariate blocks\nASCII\nDATASET UNSTRUCTURED_GRID\n"
        )
        print("POINTS", len(xyz), "double", file=stream)
        for p in xyz:
            print(*(repr(float(v)) for v in p), file=stream)
        print("CELLS", len(conn), 9 * len(conn), file=stream)
        for c in conn:
            print(8, *c, file=stream)
        print("CELL_TYPES", len(conn), file=stream)
        stream.write("12\n" * len(conn))
    return dict(
        model=path.name,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        n=n,
        vertices=len(xyz),
        hex8=len(conn),
        interfaces=len(topology["interfaces"]),
        exposed_quads=len(exposed),
        min_bound_mm3=float(min(bounds)),
        certified_positive=sum(b > 0 for b in bounds),
        vtk_min_bound_mm3=float(min(vtk_bounds)),
        vtk_certified_positive=sum(b > 0 for b in vtk_bounds),
        vtk_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
        vtk=str(out),
    )


if __name__ == "__main__":
    rows = [run(Path(p), n) for p in sys.argv[1:] for n in (1, 2, 4)]
    print(json.dumps(rows, indent=2))
