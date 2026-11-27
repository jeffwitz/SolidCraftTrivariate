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
        elif (
            o["swap"] or o["flip_second"] or shifts[a] is None or shifts[a] != shifts[b]
        ):
            status = "unsupported_transition_frame"
        else:

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
        yield dict(**record, status=status)


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--models", type=Path, default=ROOT / "build/support-census/models")
    p.add_argument("--work", type=Path, default=ROOT / "build/transverse-transitions")
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
        for path in sorted(files):
            key = path.name.split(".")[0]
            records = list(classify(json.loads(path.read_bytes())))
            count = Counter(r["status"] for r in records)
            total.update(count)
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
        )
        print(group, len(files), dict(total), flush=True)
    (args.work / "summary.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
