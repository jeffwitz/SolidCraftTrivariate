"""Cross-check the published NURBS control lattice against the native cell maps."""

from fractions import Fraction
from itertools import product
from math import comb, prod


def spline_basis(knots, degree, t):
    """Cox–de Boor basis on the closed clamped interval."""
    n = len(knots) - degree - 1
    if t == knots[-1]:
        return [0] * (n - 1) + [1]
    row = [int(a <= t < b) for a, b in zip(knots, knots[1:])]
    for k in range(1, degree + 1):
        updated = []
        for i in range(len(row) - 1):
            left = knots[i + k] - knots[i]
            right = knots[i + k + 1] - knots[i + 1]
            updated.append(
                (row[i] * (t - knots[i]) / left if left else 0)
                + (row[i + 1] * (knots[i + k + 1] - t) / right if right else 0)
            )
        row = updated
    return row


def weighted_point(controls, bases):
    x = [0] * 4
    if len(controls) != prod(map(len, bases)):
        raise ValueError("Unexpected tensor control lattice dimensions")
    for h, factors in zip(controls, product(*bases)):
        weight = prod(factors)
        if weight:
            for axis in range(4):
                x[axis] += h[axis] * weight
    if x[3] <= 0:
        raise ValueError("Nonpositive homogeneous denominator in gallery patch")
    return x


def exact_control_net(patch, span, source, offsets, shape, controls):
    """Compare every coefficient after closed-form tensor Bernstein degree elevation."""
    degrees = patch["degrees"]
    axes = span["axis_map"]
    flips = span["flips"]
    inversions = sum(x > y for i, x in enumerate(axes) for y in axes[i + 1 :])
    if (inversions + sum(flips)) % 2:
        raise ValueError(f"Orientation-reversing parameter map in {patch['id']}")
    source_degrees = source["degrees"]
    scale = None
    for out in product(*(range(d + 1) for d in degrees)):
        rows = [None] * 3
        for p, axis in enumerate(axes):
            n, m = source_degrees[axis], degrees[p]
            if n > m:
                raise ValueError("Target degree is smaller than source degree")
            i = m - out[p] if flips[p] else out[p]
            rows[axis] = [
                Fraction(comb(n, j) * comb(m - n, i - j), comb(m, i))
                if 0 <= i - j <= m - n
                else 0
                for j in range(n + 1)
            ]
        expected = weighted_point(source["controls"], rows)
        i, j, k = [v + base * d for v, base, d in zip(out, offsets, degrees)]
        pos = (i * shape[1] + j) * shape[2] + k
        actual = controls[pos]
        if scale is None:
            scale = Fraction(actual[3]) / expected[3]
        if scale <= 0 or any(scale * x != y for x, y in zip(expected, actual)):
            raise ValueError(
                f"Control mismatch in {patch['id']} for {source['id']} at {out}"
            )


def verify_patch_samples(patches, model):
    """Evaluate independent rational interior and seam witnesses on every source span."""
    cells = {b["id"]: b for b in model["blocks"]}
    visited = set()
    witnesses = [
        (Fraction(a), Fraction(b), Fraction(c))
        for a, b, c in [
            ("1/5", "2/7", "3/8"),
            ("3/4", "1/3", "4/7"),
            ("0", "2/5", "3/7"),
            ("1", "2/5", "3/7"),
        ]
    ]
    for patch in patches:
        degrees = patch["degrees"]
        count = len(patch["spans"])
        expected_knots = []
        for axis, degree in enumerate(degrees):
            n = patch["knots"][axis][-1]
            if not isinstance(n, int) or not 1 <= n <= count:
                raise ValueError(f"Invalid span count in {patch['id']}")
            expected_knots.append(
                [0] * (degree + 1)
                + [k for k in range(1, n) for _ in range(degree)]
                + [n] * (degree + 1)
            )
        if patch["knots"] != expected_knots:
            raise ValueError(f"Unexpected clamped C0 knot vectors in {patch['id']}")
        shape = [len(knots) - d - 1 for knots, d in zip(expected_knots, degrees)]
        controls = [
            [Fraction(*v) if isinstance(v, list) else Fraction(v) for v in h]
            for h in patch["controls_homogeneous"]
        ]
        if patch["control_shape"] != shape or len(controls) != prod(shape):
            raise ValueError("Invalid assembled control lattice dimensions")
        available = set(product(*(range(k[-1]) for k in expected_knots)))
        for interval, span in enumerate(patch["spans"]):
            name = span["source_block"]
            if name in visited or name not in cells:
                raise ValueError(f"Duplicated or unknown source cell {name}")
            visited.add(name)
            src = cells[name]
            if "parameter_intervals" in span:
                bounds = span["parameter_intervals"]
            else:
                bounds = [span["parameter_interval"], [0, 1], [0, 1]]
                if bounds[0] != [interval, interval + 1]:
                    raise ValueError("Unexpected legacy cell order")
            axes = span["axis_map"]
            if sorted(axes) != [0, 1, 2] or len(span["flips"]) != 3 or len(bounds) != 3:
                raise ValueError(f"Invalid parameter map for {name}")
            offsets = tuple(pair[0] for pair in bounds)
            if offsets not in available or any(b != a + 1 for a, b in bounds):
                raise ValueError(f"Missing or duplicated knot cell for {name}")
            available.remove(offsets)
            exact_control_net(patch, span, src, offsets, shape, controls)
            for uvw in witnesses:
                t = tuple(start + u for start, u in zip(offsets, uvw))
                loc = [0] * 3
                for axis, source_axis in enumerate(axes):
                    loc[source_axis] = (
                        1 - uvw[axis] if span["flips"][axis] else uvw[axis]
                    )
                bezier = [
                    [comb(d, i) * u**i * (1 - u) ** (d - i) for i in range(d + 1)]
                    for u, d in zip(loc, src["degrees"])
                ]
                spline = [
                    spline_basis(k, d, u)
                    for k, d, u in zip(patch["knots"], patch["degrees"], t)
                ]
                a = weighted_point(src["controls"], bezier)
                b = weighted_point(controls, spline)
                if any(a[i] * b[3] != b[i] * a[3] for i in range(3)):
                    raise ValueError(
                        f"Patch {patch['id']} diverges from {name} at {uvw}"
                    )
        if available:
            raise ValueError(f"Uncovered knot cells in {patch['id']}")
    if visited != set(cells):
        raise ValueError("The NURBS patches do not cover every source cell")
    return len(visited) * len(witnesses)


def verify_patch_topology(representation, model):
    """Check each remapped boundary face and external interface in source coordinates."""
    from collections import Counter

    owners = {name: p["id"] for p in representation["patches"] for name in p["cells"]}

    def endpoint(face):
        block, side = face
        return owners[block], block, side

    def record(face):
        return face["patch_id"], face["source_block"], face["side"]

    def oriented(interface, accessor):
        return (
            accessor(interface["first"]),
            accessor(interface["second"]),
            tuple(sorted(interface["orientation"].items())),
        )

    src = model["topology"]
    expected = Counter(endpoint(face) for face in src["exterior"])
    if expected != Counter(map(record, representation["boundary_faces"])):
        raise ValueError("NURBS assembly changed the exterior cell faces")
    interfaces = [
        i for i in src["interfaces"] if owners[i["first"][0]] != owners[i["second"][0]]
    ]
    if Counter(oriented(i, endpoint) for i in interfaces) != Counter(
        oriented(i, record) for i in representation["interfaces"]
    ):
        raise ValueError(
            "NURBS assembly changed a surviving interface or its orientation"
        )
