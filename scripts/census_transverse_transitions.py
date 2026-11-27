"""First-order transverse coordinate transitions on exact extruded sections.

Determine whether the inward chart derivatives admit constant or affine
transition coefficients with a fixed normalization on the second chart.
A positive result is a first-jet identity, not a global reparametrization.
"""

from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import subprocess
from check_inherited_joins import (
    add,
    mul,
    cross,
    elevated_plane,
    plane_jets,
    lifted_derivative,
)

ROOT = Path(__file__).resolve().parents[1]


def quotient_affine(a, b):
    if not b:
        raise ValueError("singular section jet")
    if not b[0]:
        raise ValueError("singular endpoint")
    const = (a[0] if a else 0) / b[0]
    slope = ((a[1] if len(a) > 1 else 0) - const * (b[1] if len(b) > 1 else 0)) / b[0]
    return (const, slope) if not add(a, mul([const, slope], b), -1) else None


def jet(block, side, flip):
    du, dv, _ = block["degrees"]
    grid = elevated_plane(block, [du, dv])

    def jets(c):
        return plane_jets([[p[c] for p in row] for row in grid], [du, dv], side, flip)

    w = jets(3)
    tangent, normal = zip(*(lifted_derivative(jets(c), w) for c in (0, 1)))
    return tangent, normal, mul(w[0], w[0])


def _divmod_poly(a, b):
    if not b:
        raise ZeroDivisionError("zero polynomial")
    r, q = list(a), [F(0)] * max(1, len(a) - len(b) + 1)
    while r and len(r) >= len(b):
        k = len(r) - len(b)
        c = r[-1] / b[-1]
        q[k] = c
        r = add(r, [F(0)] * k + [c * x for x in b], -1)
    return add(q, []), r


def _gcd_poly(a, b):
    a, b = add(a, []), add(b, [])
    while b:
        _, r = _divmod_poly(a, b)
        a, b = b, r
    return [x / a[-1] for x in a] if a else []


def _split_beta(alpha_first, alpha_second, beta):
    """Find linear beta-side data with beta=alpha_first*b_second+alpha_second*b_first."""
    at = lambda p, i: p[i] if i < len(p) else F(0)
    m = [
        [at(alpha_second, 0), F(0), at(alpha_first, 0), F(0), at(beta, 0)],
        [at(alpha_second, 1), at(alpha_second, 0), at(alpha_first, 1), at(alpha_first, 0), at(beta, 1)],
        [F(0), at(alpha_second, 1), F(0), at(alpha_first, 1), at(beta, 2)],
    ]
    row = 0
    for col in range(4):
        pivot = next((r for r in range(row, 3) if m[r][col]), None)
        if pivot is None:
            continue
        m[row], m[pivot] = m[pivot], m[row]
        d = m[row][col]
        m[row] = [x / d for x in m[row]]
        for r in range(3):
            if r != row:
                f = m[r][col]
                m[r] = [x - f * y for x, y in zip(m[r], m[row])]
        row += 1
        if row == 3:
            break
    if any(not any(m[r][:4]) and m[r][4] for r in range(row, 3)):
        return None
    x = [F(0)] * 4
    for r in range(row - 1, -1, -1):
        col = next(c for c, v in enumerate(m[r][:4]) if v)
        x[col] = (m[r][4] - sum(m[r][c] * x[c] for c in range(col + 1, 4))) / m[r][col]
    result = ([x[0], x[1]], [x[2], x[3]])
    if add(add(mul(alpha_first, result[1]), mul(alpha_second, result[0])), beta, -1):
        return None
    return result


def _power_to_bernstein(p):
    degree = len(p) - 1
    return [sum(p[k] * F(comb(i, k), comb(degree, k)) for k in range(i + 1))
            for i in range(degree + 1)]


def _half_split(values):
    row = list(values)
    left, right = [row[0]], [row[-1]]
    while len(row) > 1:
        row = [(a + b) / 2 for a, b in zip(row, row[1:])]
        left.append(row[0])
        right.append(row[-1])
    return left, list(reversed(right))


def _bernstein_sign(values, depth=12):
    if all(x > 0 for x in values):
        return 1
    if all(x < 0 for x in values):
        return -1
    if values[0] == 0 or values[-1] == 0:
        return 0
    if depth == 0:
        return None
    left, right = _half_split(values)
    a, b = _bernstein_sign(left, depth - 1), _bernstein_sign(right, depth - 1)
    if a is None or b is None:
        return None
    return 0 if a == 0 or b == 0 or a != b else a


def _nonzero_sign(p):
    return 0 if not p else _bernstein_sign(_power_to_bernstein(p))


def as_g1(block_a, side_a, block_b, side_b, flip):
    """Construct and verify fixed linear AS-G1 data on degree 1..3 extrusions."""
    da, db = block_a.get("degrees"), block_b.get("degrees")
    if (not isinstance(da, list) or not isinstance(db, list) or len(da) != 3 or len(db) != 3
            or da[2] != 1 or db[2] != 1 or da[:2] != db[:2]
            or any(not isinstance(d, int) or isinstance(d, bool) or not 1 <= d <= 3 for d in da[:2])):
        return {"status": "unsupported_family"}
    controls_a = _constant_weight_points(block_a)
    controls_b = _constant_weight_points(block_b)
    if controls_a is None or controls_b is None:
        return {"status": "unsupported_projective_weights"}
    points_a = _planar_extrusion_points(block_a)
    points_b = _planar_extrusion_points(block_b)
    if points_a is None or points_b is None:
        return {"status": "unsupported_extrusion_support"}
    trace_a = _oriented_endpoints(block_a, side_a, False)
    trace_b = _oriented_endpoints(block_b, side_b, flip)
    if trace_a is None or trace_b is None or trace_a != trace_b:
        return {"status": "unsupported_trace", "face_association_verified": False}
    ta, na, wa = jet(block_a, side_a, False)
    tb, nb, wb = jet(block_b, side_b, flip)
    if len(wa) != 1 or len(wb) != 1:
        return {"status": "unsupported_projective_weights"}
    ta, na = [[x / wa[0] for x in p] for p in ta], [[x / wa[0] for x in p] for p in na]
    tb, nb = [[x / wb[0] for x in p] for p in tb], [[x / wb[0] for x in p] for p in nb]
    if tb != ta:
        return {"status": "unsupported_trace", "face_association_verified": False}
    # Each vector component is a polynomial in the shared edge parameter.
    aa = cross(nb, ta)
    ab = cross(ta, na)
    beta = cross(na, nb)
    common = _gcd_poly(_gcd_poly(aa, ab), beta)
    if not common:
        return {"status": "indeterminate"}
    factors = []
    for p in (aa, ab, beta):
        q, r = _divmod_poly(p, common)
        assert not r
        factors.append(q)
    aa, ab, beta = factors
    assert all(
        not add(add(mul(aa, na[i]), mul(ab, nb[i])), mul(beta, ta[i]))
        for i in range(2)
    )
    degree_ok = len(aa) <= 2 and len(ab) <= 2 and len(beta) <= 3
    split = _split_beta(aa, ab, beta)
    if split is None:
        return {"status": "not_as_g1", "reason": "no_polynomial_beta_split_within_degree",
                "face_association_verified": True}
    beta_first, beta_second = split
    assert not add(add(mul(aa, beta_second), mul(ab, beta_first)), beta, -1)
    relative = len(_gcd_poly(aa, ab)) <= 1
    sign_a, sign_b = _nonzero_sign(aa), _nonzero_sign(ab)
    # Both transverse vectors point inward, so compatible gluing has equal signs.
    orientable = sign_a not in (0, None) and sign_a == sign_b
    regular_a, regular_b = _nonzero_sign(cross(ta, na)), _nonzero_sign(cross(ta, nb))
    regular = regular_a not in (0, None) and regular_b not in (0, None)
    inconclusive = None in (sign_a, sign_b, regular_a, regular_b)
    contradicted = 0 in (sign_a, sign_b, regular_a, regular_b)
    def text(p):
        return [str(x) for x in p] or ["0"]
    return {
        "status": "as_g1" if degree_ok and split is not None and relative and orientable and regular
                  else "indeterminate" if inconclusive and not contradicted else "not_as_g1",
        "gluing_factors": {"alpha_first": text(aa), "alpha_second": text(ab), "beta_first": text(beta_first), "beta_second": text(beta_second), "beta": text(beta)},
        "common_factor": text(common),
        "relatively_prime": relative,
        "orientable": orientable,
        "regular_on_closed_edge": regular,
        "face_association_verified": trace_a == trace_b,
        "curved_edge_jacobians_nonzero": regular,
        "sign_certificate": "exact Bernstein subdivision, depth 12",
        "trace_tangent_match": True,
        "beta_split_verified": True,
        "syzygy_verified": True,
        "proof": "exact determinant syzygy; polynomial gcd division has zero remainder",
    }


def _constant_weight_points(block):
    degrees = block.get("degrees")
    controls = block.get("controls")
    if not isinstance(degrees, list) or len(degrees) != 3:
        return None
    if any(not isinstance(d, int) or isinstance(d, bool) for d in degrees):
        return None
    du, dv, dw = degrees
    if du <= 0 or dv <= 0 or dw != 1 or du > 64 or dv > 64:
        return None
    if not isinstance(controls, list) or len(controls) != 2 * (du + 1) * (dv + 1):
        return None
    if any(not isinstance(p, list) or len(p) != 4 or any(not isinstance(x, int) or isinstance(x, bool) for x in p) for p in controls):
        return None
    weight = controls[0][3]
    if weight <= 0 or any(p[3] != weight for p in controls):
        return None
    return [[F(x, weight) for x in p[:3]] for p in controls]


def _planar_extrusion_points(block):
    degrees = block.get("degrees")
    if (not isinstance(degrees, list) or len(degrees) != 3
            or any(not isinstance(d, int) or isinstance(d, bool) for d in degrees)
            or any(not 1 <= d <= 3 for d in degrees[:2]) or degrees[2] != 1):
        return None
    points = _constant_weight_points(block)
    du, dv, _ = degrees
    if points is None or len(points) != 2 * (du + 1) * (dv + 1):
        return None
    base_z = points[0][2]
    height = None
    for i in range((du + 1) * (dv + 1)):
        low, high = points[2 * i], points[2 * i + 1]
        if low[2] != base_z:
            return None
        if high[:2] != low[:2]:
            return None
        delta = high[2] - low[2]
        if delta <= 0:
            return None
        if height is not None and height != delta:
            return None
        height = delta
        if i > dv and low[2] != points[2 * (i - (dv + 1))][2]:
            return None
    return points


def _oriented_endpoints(block, side, reverse):
    points = _constant_weight_points(block)
    degrees = block.get("degrees")
    if points is None or not isinstance(degrees, list) or len(degrees) != 3 or degrees[2] != 1:
        return None
    if side not in ("u0", "u1", "v0", "v1"):
        return None
    du, dv, _ = degrees
    axis, end = (0, int(side[1])) if side[0] == "u" else (1, int(side[1]))
    normal_degree, tangent_degree = (du, dv) if axis == 0 else (dv, du)
    if end > 1 or normal_degree == 0:
        return None
    result = []
    for t0 in range(tangent_degree + 1):
        t = tangent_degree - t0 if reverse else t0
        i, j = (end * normal_degree, t) if axis == 0 else (t, end * normal_degree)
        index = 2 * (i * (dv + 1) + j)
        result.extend((points[index], points[index + 1]))
    return result


def classify(model):
    blocks = {b["id"]: b for b in model["blocks"]}

    def sweep(block):
        degrees = block.get("degrees")
        controls = block.get("controls")
        if not isinstance(degrees, list) or len(degrees) != 3:
            return None
        if any(not isinstance(d, int) or isinstance(d, bool) for d in degrees):
            return None
        du, dv, dw = degrees
        if du <= 0 or dv <= 0 or dw != 1 or du > 64 or dv > 64:
            return None
        expected = 2 * (du + 1) * (dv + 1)
        if not isinstance(controls, list) or len(controls) != expected:
            return None
        if any(not isinstance(p, list) or len(p) != 4 or any(not isinstance(x, int) or isinstance(x, bool) for x in p) for p in controls):
            return None
        shifts = set()
        for a, b in zip(controls[::2], controls[1::2]):
            if a[3] != b[3] or a[3] <= 0:
                return None
            shifts.add(tuple(F(b[i] - a[i], a[3]) for i in range(3)))
        return shifts.pop() if len(shifts) == 1 else None

    shifts = {k: sweep(b) for k, b in blocks.items()}
    for link in model["topology"]["interfaces"]:
        a, sa = link["first"]
        b, sb = link["second"]
        o = link["orientation"]
        record = {"faces": [link["first"], link["second"]]}
        if "w" in [sa[0], sb[0]]:
            status = "axial_or_nonlateral"
            asg1 = {"status": "nonlateral_interface"}
        elif (
            o["swap"] or o["flip_second"] or shifts[a] is None or shifts[a] != shifts[b]
        ):
            status = "unsupported_transition_frame"
            asg1 = {"status": "unsupported_frame"}
        else:
            asg1 = as_g1(blocks[a], sa, blocks[b], sb, o["flip_first"])

            def directional(left, left_side, right, right_side):
                t, na, wa = jet(left, left_side, False)
                _, nb, wb = jet(right, right_side, o["flip_first"])
                den = mul(cross(na, t), wb)
                alpha = quotient_affine(mul(cross(nb, t), wa), den)
                beta = quotient_affine(mul(cross(na, nb), wa), den)
                if alpha is None or beta is None:
                    return "outside_affine_transition_test", {}
                if alpha[0] * (alpha[0] + alpha[1]) <= 0:
                    return "transition_singular_on_face", {}
                kind = (
                    "affine_transition"
                    if alpha[1] or beta[1]
                    else "constant_transition"
                )
                return kind, dict(
                    alpha=list(map(str, alpha)), beta=list(map(str, beta))
                )

            status, coeffs = directional(blocks[a], sa, blocks[b], sb)
            reverse, reverse_coeffs = directional(blocks[b], sb, blocks[a], sa)
            record.update(
                coeffs, reversed_status=reverse, reversed_coefficients=reverse_coeffs
            )
        yield dict(**record, status=status, as_g1=asg1)


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--models", type=Path, default=ROOT / "build/support-census/models")
    p.add_argument("--work", type=Path, default=ROOT / "build/transverse-transitions")
    p.add_argument("--cli", type=Path, help="compare the full census with the native CLI")
    args = p.parse_args()
    args.work.mkdir(exist_ok=True, parents=True)
    models = list(args.models.glob("*.model.json"))
    if not models:
        raise RuntimeError("run the support census first")
    families = {
        "polygon": models,
        "named": [
            ROOT / "docs/site/box-cylinder/tall/cylinder.model.json",
            ROOT / "docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json",
        ],
    }
    report = {}
    for group, files in families.items():
        total = Counter()
        either = Counter()
        cases = {}
        asg1_total = Counter()
        for path in sorted(files):
            key = path.name.split(".")[0]
            records = list(classify(json.loads(path.read_bytes())))
            if args.cli and path.name.endswith(".model.json"):
                native_path = args.work / (key + ".native.json")
                proc = subprocess.run(
                    [str(args.cli), "classify-transverse-transitions", str(path), str(native_path)],
                    capture_output=True,
                    text=True,
                )
                if proc.returncode:
                    raise RuntimeError(f"native AS-G1 census failed for {path.name}: {proc.stdout}{proc.stderr}")
                if json.loads(native_path.read_bytes()) != records:
                    raise AssertionError(f"native/Python AS-G1 mismatch for {path.name}")
            count = Counter(r["status"] for r in records)
            total.update(count)
            asg1_total.update(r["as_g1"]["status"] for r in records)
            for r in records:
                options = (r["status"], r.get("reversed_status"))
                either.update(
                    [
                        next(
                            (
                                k
                                for k in ["constant_transition", "affine_transition"]
                                if k in options
                            ),
                            r["status"],
                        )
                    ]
                )
            cases[key] = dict(count)
            (args.work / (key + ".interfaces.json")).write_text(
                json.dumps(records) + "\n"
            )
        report[group] = dict(
            models=len(files),
            total=dict(total),
            either_direction=dict(either),
            cases=cases,
            as_g1=dict(asg1_total),
        )
        print(group, len(files), dict(total), flush=True)
    (args.work / "summary.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
