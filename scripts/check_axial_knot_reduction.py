"""Rebuild C0 nets without reference controls and check compact geometry in exact arithmetic."""

import argparse
import copy
from fractions import Fraction
from itertools import product
import json
from math import prod
from pathlib import Path
import subprocess

from check_gallery_patch_maps import (
    spline_basis,
    verify_patch_samples,
    verify_patch_topology,
    weighted_point,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/site/box-cylinder/tall"
OUT = ROOT / "build/axial-knot-reduction"
REPORT = ROOT / "docs/results/axial-knot-reduction/geometry.json"


def read(path):
    return json.loads(path.read_bytes())


def knot_cells(patch):
    return prod(len(set(k)) - 1 for k in patch["knots"])


def values(patch, point):
    bases = [
        spline_basis(k, d, u)
        for k, d, u in zip(patch["knots"], patch["degrees"], point)
    ]
    return weighted_point(patch["controls_homogeneous"], bases)


def check_reinsertion(old, new):
    axes = {span["axis_map"].index(2) for span in old["spans"]}
    assert (
        len(axes) == 1 and old["id"] == new["id"] and new["degrees"] == old["degrees"]
    )
    z = axes.pop()
    assert old["degrees"][z] == 1
    k = new["knots"][z]
    for a in range(3):
        if a != z:
            assert (
                old["knots"][a] == new["knots"][a]
                and old["control_shape"][a] == new["control_shape"][a]
            )
    assert new["control_shape"][z] == len(k) - 2 and len(
        new["controls_homogeneous"]
    ) == prod(new["control_shape"])
    assert (
        k == sorted(k)
        and k[:2] == old["knots"][z][:2]
        and k[-2:] == old["knots"][z][-2:]
    )
    assert set(k).issubset(old["knots"][z])
    _, ny, nz = new["control_shape"]
    for idx, expected in zip(
        product(*(range(n) for n in old["control_shape"])),
        old["controls_homogeneous"],
        strict=True,
    ):
        t = Fraction(old["knots"][z][idx[z] + 1])
        b = spline_basis(k, 1, t)
        restored = [0] * 4
        coord = list(idx)
        for j, coefficient in enumerate(b):
            coord[z] = j
            point = new["controls_homogeneous"][
                (coord[0] * ny + coord[1]) * nz + coord[2]
            ]
            for h in range(4):
                restored[h] += coefficient * point[h]
        assert restored == expected, (old["id"], idx, t, expected, restored)
    # Independent Cox–de Boor witnesses on the *old* source integration cells.
    for intervals in product(
        *[[(a, b) for a, b in zip(ks, ks[1:]) if a < b] for ks in old["knots"]]
    ):
        u = [a + Fraction(i + 2, 7) * (b - a) for i, (a, b) in enumerate(intervals)]
        assert values(old, u) == values(new, u)


def check_case(name, cli):
    ref = read(SOURCE / f"{name}.patches.json")
    model_path = SOURCE / ref["source_model"]
    layout = copy.deepcopy(ref)
    for patch in layout["patches"]:
        for key in (
            "controls_homogeneous",
            "control_shape",
            "knots",
            "homogeneous_common_scale_denominator",
        ):
            del patch[key]
    plan, assembled, compact = [
        OUT / f"{name}{suffix}.json" for suffix in (".layout", ".assembled", "")
    ]
    plan.write_text(json.dumps(layout))
    subprocess.run([cli, "assemble-patches", model_path, plan, assembled], check=True)
    native = read(assembled)
    model = read(model_path)
    verify_patch_samples(native["patches"], model)
    verify_patch_topology(native, model)
    for a, b in zip(native["patches"], ref["patches"], strict=True):
        assert a["knots"] == b["knots"] and a["control_shape"] == b["control_shape"]
        # A common positive homogeneous factor leaves the rational patch unchanged.
        factor = Fraction(
            a["controls_homogeneous"][0][3], b["controls_homogeneous"][0][3]
        )
        assert factor > 0 and all(
            x == y * factor
            for p, q in zip(
                a["controls_homogeneous"], b["controls_homogeneous"], strict=True
            )
            for x, y in zip(p, q, strict=True)
        )
    subprocess.run([cli, "reduce-axial-knots", assembled, compact], check=True)
    reduced = read(compact)
    assert reduced["metadata"] == {k: v for k, v in native.items() if k != "patches"}
    patches = native["patches"]
    for before, after, provenance, analysis in zip(
        patches,
        reduced["patches"],
        reduced["provenance"]["patches"],
        reduced["analysis_knots"],
        strict=True,
    ):
        assert analysis == {"patch_id": before["id"], "knots": before["knots"]}
        for field in ("id", "spans", "source_blocks", "cells"):
            assert provenance[field] == before[field]
        check_reinsertion(before, after)
    return dict(
        design_volumes=len(patches),
        source_bezier_cells=len(model["blocks"]),
        geometry_spans_before=sum(map(knot_cells, patches)),
        geometry_spans_after=sum(map(knot_cells, reduced["patches"])),
        stored_controls_before=sum(len(p["controls_homogeneous"]) for p in patches),
        stored_controls_after=sum(
            len(p["controls_homogeneous"]) for p in reduced["patches"]
        ),
        retained_subface_joins=len(native["interfaces"]),
        retained_self_seam_records=len(native.get("self_seams", [])),
        retained_boundary_records=len(native["boundary_faces"]),
    )


def check_atomic_refusals(cli):
    import tempfile

    with tempfile.TemporaryDirectory(dir=OUT) as tmp:
        bad = Path(tmp) / "bad.json"
        target = Path(tmp) / "existing.json"
        target.write_text("preserve this file")
        layout = read(OUT / "cylinder-2.layout.json")
        layout["self_seams"].pop()
        model = SOURCE / layout["source_model"]
        bad.write_text(json.dumps(layout))
        commands = [
            [cli, "assemble-patches", model, bad, target],
            [cli, "reduce-axial-knots", bad, target],
        ]
        for command in commands:
            attempt = subprocess.run(command, capture_output=True)
            assert attempt.returncode != 0 and b"panicked" not in attempt.stderr, (
                command
            )
            assert target.read_text() == "preserve this file", command


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate"
    )
    parser.add_argument("--write-report", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        name: check_case(name, args.cli)
        for name in ("cylinder-2", "cylinder-3", "union")
    }
    check_atomic_refusals(args.cli)
    print(json.dumps(report, indent=2))
    if args.write_report:
        REPORT.write_text(json.dumps(report, indent=2) + "\n")
    elif read(REPORT) != report:
        raise ValueError("Geometry report differs")


if __name__ == "__main__":
    main()
