"""Recompute volumes, material membership and topology of the tall-cylinder gallery."""

import argparse
from bisect import bisect_right
import json
import math
from pathlib import Path

SITE = Path(__file__).resolve().parents[1] / "docs/site/box-cylinder/tall"
# Each tuple is the material truth table ordered by 2*A + B, the analytic
# volume as constant + coefficient*pi, and the number of face-connected parts.
STATES = {
    "box": ("0011", 19200, 0, 1),
    "cylinder": ("0101", 0, 900, 1),
    "union": ("0111", 19200, 600, 1),
    "intersection": ("0001", 0, 300, 1),
    "a-minus-b": ("0010", 19200, -300, 1),
    "b-minus-a": ("0100", 0, 600, 2),
    "xor": ("0110", 19200, 300, 3),
}


def value(net, t):
    degree = len(net) - 1
    basis = [
        math.comb(degree, i) * t**i * (1 - t) ** (degree - i) for i in range(degree + 1)
    ]
    return [sum(b * q[k] for b, q in zip(basis, net)) for k in range(4)]


def section(block):
    nu, nv, nz = [d + 1 for d in block["degrees"]]
    controls = block["controls"]
    if nz != 2 or len(controls) != nu * nv * nz:
        raise ValueError("Expected an axially linear extrusion")
    lo, hi = [z / w for _, _, z, w in controls[:2]]
    if hi <= lo:
        raise ValueError("Invalid axial extent")
    for p, q in zip(controls[::2], controls[1::2]):
        if p[0:2] != q[0:2] or p[3] != q[3] or p[2] != lo * p[3] or q[2] != hi * q[3]:
            raise ValueError("The exported volume is not a straight extrusion")
    plane = controls[::2]
    edges = [plane[::nv], plane[-nv:], plane[nv - 1 :: nv][::-1], plane[:nv][::-1]]
    outline = []
    signed_area = 0.0
    for edge in edges:
        degree = len(edge) - 1
        diff = [
            [degree * (b - a) for a, b in zip(p, q)] for p, q in zip(edge, edge[1:])
        ]
        # Composite Simpson integration of 1/2*(x dy - y dx).
        n = 256
        for i in range(n + 1):
            h = value(edge, i / n)
            dh = value(diff, i / n)
            signed_area += (
                (1 if i in (0, n) else 4 if i % 2 else 2)
                * (h[0] * dh[1] - h[1] * dh[0])
                / (6 * n * h[3] ** 2)
            )
            if i < n:
                outline.append((h[0] / h[3], h[1] / h[3]))
    if signed_area <= 0:
        raise ValueError("Invalid oriented cross-section")
    return lo, hi, outline, signed_area * (hi - lo)


def membership_check(state, cells):
    xy = [(x + 0.231, y + 0.419) for x in range(-2, 43, 3) for y in range(-2, 43, 3)]
    for i in range(32):
        theta = 2 * math.pi * (i + 0.173) / 32
        xy.extend(
            (20 + r * math.cos(theta), 20 + r * math.sin(theta)) for r in (4.99, 5.01)
        )
    zs = (-12.01, -11.99, -6, -0.01, 0.01, 6, 11.99, 12.01, 18, 23.99, 24.01)
    for x, y in xy:
        intervals = []
        for lo, hi, loop, _ in cells:
            cuts = sorted(
                a + (y - b) * (c - a) / (d - b)
                for (a, b), (c, d) in zip(loop, loop[1:] + loop[:1])
                if (b > y) != (d > y)
            )
            if bisect_right(cuts, x) % 2:
                intervals.append((lo, hi))
        for z in zs:
            a = 0 < x < 40 and 0 < y < 40 and 0 < z < 12
            b = (x - 20) ** 2 + (y - 20) ** 2 < 25 and -12 < z < 24
            if sum(lo < z < hi for lo, hi in intervals) != int(
                STATES[state][0][2 * a + b]
            ):
                raise ValueError(f"{state} material membership differs at {(x, y, z)}")
    return len(xy) * len(zs)


def components(model):
    links = {b["id"]: set() for b in model["blocks"]}
    for interface in model["topology"]["interfaces"]:
        a, b = [interface[key][0] for key in ("first", "second")]
        links[a].add(b)
        links[b].add(a)
    remaining = set(links)
    count = 0
    while remaining:
        stack = [remaining.pop()]
        count += 1
        while stack:
            neighbors = links[stack.pop()] & remaining
            remaining.difference_update(neighbors)
            stack.extend(neighbors)
    return count


def manifest():
    classic = json.loads((SITE.parent / "classic.model.json").read_bytes())

    def control_maps(model):
        return sorted(
            (block["degrees"], block["controls"]) for block in model["blocks"]
        )

    records = []
    for state, (_, constant, coeff, count) in STATES.items():
        model = json.loads((SITE / f"{state}.model.json").read_bytes())
        patches = json.loads((SITE / f"{state}.patches.json").read_bytes())["patches"]
        if state == "a-minus-b" and control_maps(model) != control_maps(classic):
            raise ValueError(
                "Changing cutter height changed the original four-sector difference"
            )
        cells = [section(block) for block in model["blocks"]]
        volume = constant + math.pi * coeff
        if components(model) != count or not math.isclose(
            sum(c[3] for c in cells), volume, rel_tol=1e-9
        ):
            raise ValueError(f"{state} volume or component count differs")
        probes = membership_check(state, cells)
        records.append(
            dict(
                id=state,
                bezier_cells=len(cells),
                nurbs_volumes=len(patches),
                volume_mm3=volume,
                volume_constant_mm3=constant,
                volume_pi_coefficient_mm3=coeff,
                components=count,
                material_probes=probes,
            )
        )
    return dict(
        stock_size_mm=[40, 40, 12],
        cylinder_radius_mm=5,
        cylinder_interval_mm=[-12, 24],
        states=records,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    expected = manifest()
    path = SITE / "manifest.json"
    if args.write:
        path.write_text(json.dumps(expected, indent=2) + "\n")
    elif json.loads(path.read_bytes()) != expected:
        raise SystemExit(
            "Tall-cylinder manifest differs from independent geometry checks"
        )
    print("Tall-cylinder geometry passed for all seven states.")
