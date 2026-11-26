"""Replay the native rim edits across the five published FieldCraft brackets."""

import argparse
import json
from math import pi
from pathlib import Path
import subprocess

import numpy as np
from s3_stock_rim import certify, integrate

ROOT = Path(__file__).resolve().parents[2]


def check_state(source, result, state):
    assert source["topology"] == result["parent_topology"]
    assert {b["id"] for b in source["blocks"]} == {v["id"] for v in result["volumes"]}
    total = 0.0
    regular = singular = 0
    for v in result["volumes"]:
        assert v["parameter_intervals_32nds"] == [[0, 1], [1, 31], [31, 32]]
        assert v["span_axis"] == ("w" if v["id"].startswith("wing_") else "v")
        for element in v["spans"]:
            net = np.array(element["controls"], dtype=object).reshape(
                *(d + 1 for d in element["degrees"]), 4
            )
            zeros = certify(net)
            if zeros:
                assert state == "fillet" and v["id"].startswith("wing_")
                # Original axes: the bore wall is v=1 and each contact circle
                # lies on a cap, at w=0 or w=1 depending on the mouth.
                assert len(zeros) == 3
                assert {u for u, _, _ in zeros} == {-1, 0, 1}
                assert {v_ for _, v_, _ in zeros} == {1}
                cap = {w for _, _, w in zeros}
                assert cap in ({-1}, {1})
                singular += 1
            else:
                regular += 1
            total += integrate(net, 24)
    assert singular == (16 if state == "fillet" else 0)
    return dict(
        volume_mm3=total,
        logical_volumes=len(result["volumes"]),
        closed_regular_spans=regular,
        boundary_singular_spans=singular,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cli", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    args = p.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    outcomes = []
    inputs = sorted(
        (ROOT / "docs/site/data").glob("fieldcraft_*/fieldcraft_*.model.json")
    )
    assert len(inputs) == 5
    for recipe in inputs:
        source = json.loads(recipe.read_text())
        scales = source["transform"]["denominator"]
        radii = [b["radius_mm"] / scales for b in source["source_model"]["bores"]]
        f = 0.25
        ref = {
            "chamfer": sum(2 * pi * (r * f**2 + f**3 / 3) for r in radii),
            "fillet": sum(
                2 * pi * (2 * r * f**2 * (1 - pi / 4) + f**3 * (5 / 3 - pi / 2))
                for r in radii
            ),
        }
        row = {"id": recipe.stem.removesuffix(".model")}
        for state in ("sharp", "chamfer", "fillet"):
            dest = args.work / (row["id"] + "-" + state + ".json")
            subprocess.run(
                [args.cli, "edit-fieldcraft-bore-rims", recipe, state, dest],
                check=True,
                stdout=subprocess.DEVNULL,
            )
            row[state] = check_state(source, json.loads(dest.read_text()), state)
            if state != "sharp":
                measured = row["sharp"]["volume_mm3"] - row[state]["volume_mm3"]
                row[state]["removed_mm3"] = measured
                row[state]["expected_removed_mm3"] = ref[state]
                assert abs(measured - ref[state]) < 1e-7, (
                    row["id"],
                    state,
                    measured,
                    ref[state],
                )
        outcomes.append(row)
        print(row["id"], "passed", flush=True)
    (args.work / "summary.json").write_text(json.dumps(outcomes, indent=2) + "\n")


if __name__ == "__main__":
    main()
