"""Exact Fraction audit of the emitted meridional annulus boundary and joins."""

from fractions import Fraction as F
import json
import sys
from pathlib import Path


MID = ((768, 0), (1024, 1024), (768, 2048), (512, 1024))
CORNER = ((1024, 0), (1024, 2048), (512, 2048), (512, 0))


def xy(chart, i, j=0, k=0):
    nx, ny, nz = chart["shape"]
    assert i < nx and j < ny and k < nz
    return chart["xyzw"][(i * ny + j) * nz + k]


def point(hp):
    # The projective denominator includes the fixed 256 coordinate scale.
    return F(hp[0], hp[3]), F(hp[2], hp[3])


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def sub(a, b):
    return a[0] - b[0], a[1] - b[1]


def span_controls(chart, span):
    assert chart["degrees"][0] == 2
    return [xy(chart, 2 * span + j) for j in range(3)]


def bezier_eval(controls, t):
    b = ((1 - t) ** 2, 2 * t * (1 - t), t**2)
    den = sum(bi * h[3] for bi, h in zip(b, controls, strict=True))
    return tuple(sum(bi * h[k] for bi, h in zip(b, controls, strict=True)) / den for k in (0, 2))


def poly_add(a, b, scale=1):
    out = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += scale * x
    return out


def poly_mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def circle_identity(controls, center, radius):
    # Convert homogeneous quadratic Bezier coordinates to power coefficients.
    power = []
    for axis in (0, 2, 3):
        p0, p1, p2 = (F(h[axis]) for h in controls)
        power.append((p0, 2 * (p1 - p0), p0 - 2 * p1 + p2))
    x, z, w = power
    x = poly_add(x, w, -center[0])
    z = poly_add(z, w, -center[1])
    identity = poly_add(poly_mul(x, x), poly_mul(z, z))
    identity = poly_add(identity, poly_mul(w, w), -(radius**2))
    assert all(c == 0 for c in identity), (center, radius, identity)


def endpoint_tangent(controls, at_end):
    a, b = (controls[1], controls[2]) if at_end else (controls[0], controls[1])
    # The positive scalar factor is irrelevant; this is the exact rational
    # derivative direction of X/W, Z/W at the corresponding endpoint.
    return (
        F(b[0] * a[3] - a[0] * b[3]),
        F(b[2] * a[3] - a[2] * b[3]),
    )


def projectively_equal(a, b):
    return all(a[i] * b[3] == b[i] * a[3] for i in range(4))


def audit_file(path, expected_state):
    model = json.loads(Path(path).read_text())
    charts = {c["id"]: c for c in model["charts"]}
    bands = [charts[f"turn-0-section-{i}"] for i in range(4)]
    assert all(c["shape"] == [9, 5, 2] for c in bands)
    for chart in bands:
        assert all(h[3] > 0 for h in chart["xyzw"])

    # The four radial interfaces (including the periodic last-to-first join)
    # must have the same exact homogeneous traces, for all angular controls.
    for i in range(4):
        left, right = bands[i], bands[(i + 1) % 4]
        for j in range(left["shape"][1]):
            for k in range(left["shape"][2]):
                assert xy(left, 8, j, k) == xy(right, 0, j, k), (i, j, k)

    for side, chart in enumerate(bands):
        a = tuple(F(x, 256) for x in MID[side])
        b = tuple(F(x, 256) for x in MID[(side + 1) % 4])
        c = tuple(F(x, 256) for x in CORNER[side])
        toward = lambda q: tuple(ci + F(1, 4) * (qi - ci > 0 and 1 or qi - ci < 0 and -1 or 0)
                                 for qi, ci in zip(q, c, strict=True))
        t1, t2 = toward(a), toward(b)
        controls = [span_controls(chart, s) for s in range(4)]

        # Exact, rational control endpoints document all four wall spans.
        assert point(controls[0][0]) == a
        assert point(controls[0][2]) == t1
        assert point(controls[3][0]) == t2
        assert point(controls[3][2]) == b

        if expected_state == "sharp":
            assert point(controls[1][0]) == t1
            assert point(controls[1][2]) == c
            assert point(controls[2][0]) == c
            assert point(controls[2][2]) == t2
            for span in controls:
                p0, p1, p2 = (point(h) for h in span)
                assert cross(sub(p1, p0), sub(p2, p0)) == 0
        elif expected_state == "chamfer":
            bevel_mid = tuple((x + y) / 2 for x, y in zip(t1, t2, strict=True))
            assert point(controls[1][0]) == t1 and point(controls[1][2]) == bevel_mid
            assert point(controls[2][0]) == bevel_mid and point(controls[2][2]) == t2
            for span in controls[1:3]:
                p0, p1, p2 = (point(h) for h in span)
                assert cross(sub(p1, p0), sub(p2, p0)) == 0
        else:
            center = tuple(x + y - z for x, y, z in zip(t1, t2, c, strict=True))
            for span in controls[1:3]:
                circle_identity(span, center, F(1, 4))
            # Fillet endpoint tangents equal the adjoining wall directions;
            # the join between its two exact rational spans is also tangent.
            tangent_a = endpoint_tangent(controls[1], False)
            tangent_join_l = endpoint_tangent(controls[1], True)
            tangent_join_r = endpoint_tangent(controls[2], False)
            tangent_b = endpoint_tangent(controls[2], True)
            incoming = sub(t1, a)
            outgoing = sub(b, t2)
            assert cross(tangent_a, incoming) == 0
            assert cross(tangent_b, outgoing) == 0
            assert cross(tangent_join_l, tangent_join_r) == 0

        # All emitted span weights are strictly positive, and endpoints join
        # projectively even where a span uses a different homogeneous gauge.
        for span in controls:
            assert all(h[3] > 0 for h in span)
        for i in range(3):
            assert projectively_equal(controls[i][2], controls[i + 1][0])

    print(f"{expected_state}: 4 exact corner traces; radius/tangency/joins pass")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        root = Path(__file__).resolve().parents[2] / "docs/results/annulus-fixed-topology"
        args = [root / f"ten-{name}.json" for name in ("sharp", "chamfer", "fillet")]
    else:
        args = list(map(Path, sys.argv[1:]))
    assert len(args) == 3
    for state, path in zip(("sharp", "chamfer", "fillet"), args, strict=True):
        audit_file(path, state)
