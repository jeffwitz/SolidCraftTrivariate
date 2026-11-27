"""Compare each canonical extraction column to the source univariate basis."""

import argparse
import json
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import subprocess

from check_spline_parent import basis_values, bernstein

ROOT = Path(__file__).resolve().parents[1]


def rows(records):
    return [{i: F(int(a), int(b)) for i, (a, b) in r} for r in records]


def invoke(cli, *args):
    subprocess.run([cli, *args], check=True, capture_output=True)


def check(cli, work):
    work.mkdir(parents=True, exist_ok=True)
    source = ROOT / "tests/references/nurbs_strip.recipe.json"
    recipe = json.loads(source.read_bytes())
    model = work / "model.json"
    invoke(cli, "build", source, model)
    data = json.loads(model.read_bytes())
    path = work / "field.json"
    invoke(cli, "export-parent-field", model, path)
    invoke(cli, "audit-parent-field", model, path)
    galerkin = json.loads(
        subprocess.check_output([cli, "check-parent-galerkin", model, "8"])
    )
    assert (galerkin["c0_scalar_dofs"], galerkin["exact_transfer_rank"]) == (28, 20)
    field = json.loads(path.read_bytes())
    dofs = field["scalar_dofs"]
    weights = [q[2] for q in recipe["lower"]["controls"] for _ in range(4)]
    assert dofs == 20
    verified = 0
    for cell, span, block in zip(
        field["elements"], data["span_intervals"], data["blocks"], strict=True
    ):
        E, C = map(rows, (cell["E"], cell["C"]))
        assert cell["source_block"] == block["id"]
        for r, c in zip(E, C, strict=True):
            weight = sum(x * weights[i] for i, x in r.items())
            assert c == {i: x * weights[i] / weight for i, x in r.items()}
        for u, v, w in product(
            [F(0), F(1, 5), F(4, 5), F(1)], [F(0), F(2, 7), F(1)], [F(1, 7), F(1)]
        ):
            t = span[0] + u * (span[1] - span[0])
            direct = basis_values(3, recipe["lower"]["knots"], t)
            poly = [x * y * z for x in direct for y in (1 - v, v) for z in (1 - w, w)]
            # B_e times E, with a separate rational path using the source
            # cell controls and C, including all columns outside local support.
            b = [
                bernstein(3, i, u) * y * z
                for i in range(4)
                for y in (1 - v, v)
                for z in (1 - w, w)
            ]
            w_e = [h[3] for h in block["controls"]]
            W = sum(x * y for x, y in zip(b, w_e))
            denominator = sum(a * b for a, b in zip(poly, weights))
            for i in range(dofs):
                assert sum(x * r.get(i, 0) for x, r in zip(b, E)) == poly[i]
                assert (
                    sum(x * we / W * r.get(i, 0) for x, we, r in zip(b, w_e, C))
                    == poly[i] * weights[i] / denominator
                )
                verified += 1
    plan = work / "plan.json"
    plan.write_text(
        json.dumps(
            dict(
                schema="inherited_field_plan_v1",
                axial_degree=1,
                axial_continuity=0,
                restore_strip=True,
                transverse_interfaces=[],
            )
        )
    )
    package = work / "jet_field.json"
    invoke(cli, "export-inherited-space", model, plan, package)
    assert json.loads(package.read_bytes())["field"]["scalar_dofs"] == dofs
    solves = {}
    for q in (8, 12, 16, 20, 24, 32):
        job = work / "job.json"
        job.write_text(
            json.dumps(
                dict(
                    schema="section_quadratic_mms_v1",
                    x_origin_mm=0,
                    x_scale_mm=4,
                    quadrature=q,
                )
            )
        )
        first = work / f"canonical-{q}.json"
        second = work / f"jet-{q}.json"
        invoke(cli, "solve-parent-field", model, job, first)
        invoke(cli, "solve-inherited-space", package, job, second)
        a, b = [json.loads(p.read_bytes()) for p in (first, second)]
        differences = {
            k: abs(a[k] - b[k]) / max(abs(a[k]), 1e-25)
            for k in ("energy_N_mm", "mean_displacement_z_mm", "relative_energy_norm")
        }
        assert max(differences.values()) < 1e-9
        previous = list(solves.values())[-1]["canonical"] if solves else None
        solves[q] = dict(
            canonical=a, independent_coordinates=b, relative_differences=differences
        )
        if previous:
            drift = max(
                abs(a[k] - previous[k]) / abs(a[k])
                for k in ("energy_N_mm", "mean_displacement_z_mm")
            )
            solves[q]["relative_quadrature_change"] = drift
            if drift < 1e-8:
                break
    assert solves[q].get("relative_quadrature_change", 1) < 1e-8
    output = dict(
        galerkin=galerkin,
        exact_column_point_checks=verified,
        scalar_dofs=dofs,
        runs=solves,
    )
    (work / "summary.json").write_text(json.dumps(output, indent=2) + "\n")
    print(
        verified,
        "coefficient-column evaluations per representation;",
        dofs,
        "scalar DOFs; quadrature",
        q,
        "relative drift",
        drift,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate"
    )
    parser.add_argument(
        "--work", type=Path, default=ROOT / "build/canonical-parent-field"
    )
    args = parser.parse_args()
    check(args.cli.resolve(), args.work)
