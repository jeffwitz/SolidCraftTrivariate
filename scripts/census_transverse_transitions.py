"""First-order transverse coordinate transitions on exact extruded sections.

Determine whether the inward chart derivatives admit constant or affine
transition coefficients with a fixed normalization on the second chart.
A positive result is a first-jet identity, not a global reparametrization.
"""

from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as F
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
    return ([x[0], x[1]], [x[2], x[3]])


def as_g1(block_a, side_a, block_b, side_b, flip):
    """Construct and exactly verify AS-G1 linear gluing data for bilinear extrusions."""
    if block_a["degrees"] != [1, 1, 1] or block_b["degrees"] != [1, 1, 1]:
        return {"status": "unsupported_family"}
    ta, na, wa = jet(block_a, side_a, False)
    tb, nb, wb = jet(block_b, side_b, flip)
    if len(wa) != 1 or len(wb) != 1:
        return {"status": "unsupported_projective_weights"}
    ta, na = [[x / wa[0] for x in p] for p in ta], [[x / wa[0] for x in p] for p in na]
    tb, nb = [[x / wb[0] for x in p] for p in tb], [[x / wb[0] for x in p] for p in nb]
    if tb != ta:
        return {"status": "unsupported_trace"}
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
    split = _split_beta(aa, ab, beta)
    if split is None:
        return {"status": "not_as_g1", "reason": "no_linear_beta_split"}
    beta_first, beta_second = split
    assert not add(add(mul(aa, beta_second), mul(ab, beta_first)), beta, -1)
    relative = len(_gcd_poly(aa, ab)) <= 1
    def endpoint_sign(p):
        values = (sum(p, F(0)), p[0] if p else F(0))
        if not all(values) or values[0] * values[1] < 0:
            return 0
        return 1 if values[0] > 0 else -1
    def no_zero(p):
        s = endpoint_sign(p)
        if not s:
            return False
        if len(p) <= 2:
            return True
        if len(p) != 3 or not p[2]:
            return False
        vertex = -p[1] / (2 * p[2])
        if not 0 < vertex < 1:
            return True
        value = p[0] - p[1] ** 2 / (4 * p[2])
        return value > 0 if s > 0 else value < 0
    # Both transverse vectors point inward, so compatible gluing has equal signs.
    orientable = endpoint_sign(aa) == endpoint_sign(ab) and endpoint_sign(aa) != 0
    regular = no_zero(cross(ta, na)) and no_zero(cross(ta, nb))
    # AS-G1 uses linear alpha factors and permits a quadratic beta factor.
    linear = len(aa) <= 2 and len(ab) <= 2 and len(beta) <= 3
    def text(p):
        return [str(x) for x in p] or ["0"]
    return {
        "status": "as_g1" if linear and relative and orientable and regular else "not_as_g1",
        "gluing_factors": {"alpha_first": text(aa), "alpha_second": text(ab), "beta_first": text(beta_first), "beta_second": text(beta_second), "beta": text(beta)},
        "common_factor": text(common),
        "relatively_prime": relative,
        "orientable": orientable,
        "regular_on_closed_edge": regular,
        "trace_tangent_match": True,
        "beta_split_verified": True,
        "syzygy_verified": True,
        "proof": "exact determinant syzygy; polynomial gcd division has zero remainder",
    }


def classify(model):
    blocks = {b["id"]: b for b in model["blocks"]}

    def sweep(block):
        if block["degrees"][2] != 1:
            return None
        shifts = set()
        for a, b in zip(block["controls"][::2], block["controls"][1::2]):
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
