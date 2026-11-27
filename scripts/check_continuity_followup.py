"""Compare native uniform-section and cube mechanics to the published observables."""

import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/results/continuity-followup"
ERRORS = ("relative_energy_norm", "relative_displacement_L2")


def load(path):
    return json.loads(path.read_text())


def close(a, b, label, atol=2e-12, rtol=5e-7):
    if not math.isclose(a, b, rel_tol=rtol, abs_tol=atol):
        raise ValueError(f"{label}: {a} differs from {b}")


def pocket(directory):
    rows = []
    for p in (3, 4):
        for c in (0, 1):
            last = None
            for n in (1, 2, 4, 8):

                def row(q):
                    return load(directory / f"p{p}-n{n}-c{c}-q{q}.json")

                record = row(8)
                result = dict(p=p, n=n, continuity=c)
                for key in (*ERRORS, "free_vector_dofs"):
                    result[key] = record[key]
                if n >= 4:
                    repeat = row(10)
                    for key in ERRORS:
                        close(
                            record[key],
                            repeat[key],
                            f"q8/q10 {p}/{n}/{c} {key}",
                            2e-12,
                            2e-4,
                        )
                if last:
                    for key in ERRORS:
                        result[key + "_rate"] = math.log2(last[key] / record[key])
                rows.append(result)
                last = record
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, default=ROOT / "build")
    parser.add_argument("--write-pocket", action="store_true")
    args = parser.parse_args()
    rows = pocket(args.build / "pocket-vector")
    path = REPORT / "pocket-vector.json"
    if args.write_pocket:
        path.write_text(json.dumps(rows, indent=2) + "\n")
    expected = load(path)
    if len(rows) != len(expected):
        raise ValueError("pocket observation count changed")
    for current, published in zip(rows, expected):
        for key in current:
            close(current[key], published[key], str((current["p"], current["n"], key)))
    for pair in load(REPORT / "cube-matched.json")["results"]:
        n = pair["C0"]["cells_per_axis"]
        q = pair["job"]["quadrature"]
        raw = load(args.build / "cube-convergence" / f"matched-gaussian-n{n}-q{q}.json")
        if raw["job"] != pair["job"]:
            raise ValueError("cube load or quadrature changed")
        for space in ("C0", "C2"):
            for key in (
                *ERRORS,
                "scalar_dofs",
                "free_vector_dofs",
                "quadrature_points",
            ):
                close(
                    raw[space]["result"][key],
                    pair[space][key],
                    f"cube {n}/{q}/{space}/{key}",
                )
    print(f"{len(rows)} vector pocket states and three paired cube solves reproduced")


if __name__ == "__main__":
    main()
