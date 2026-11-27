"""Replay a user-declared native mechanical job and check every field coefficient."""

import argparse
import json
import subprocess
from fractions import Fraction
from pathlib import Path


def verify(cli, work):
    package = work / "package.json"
    model = json.loads(package.read_text())["geometry"]
    extraction = json.loads((work / "extraction.json").read_text())
    target = work / "traction.json"
    request = Path("tests/references/thb_end_traction.job.json")
    subprocess.run([cli, "solve-job", package, request, target], check=True)
    result = json.loads(target.read_text())
    height = Fraction(model["height_mm"])
    # Exact stress = (0,0,1 MPa) in the axial component and E = 10 MPa.
    # w varies fastest in the recovered parent. Physical amplitudes transfer
    # with the rational, weight-scaled operators.
    uz = [
        k * height / 10
        for _ in model["lower"]["controls"]
        for _ in (0, 1)
        for k in (0, 1)
    ]
    for step in ("parent_to_initial", "initial_to_prepared", "prepared_to_thb"):
        uz = [
            sum(Fraction(int(p), int(q)) * uz[i] for i, (p, q) in row)
            for row in extraction["rational_transfer_" + step]
        ]
    values = result["coefficients_xyz_mm"]
    assert len(values) == len(uz)
    err = max(
        max(abs(x), abs(y), abs(z - float(target)))
        for (x, y, z), target in zip(values, uz)
    )
    energy = result["volume_mm3"] / 20
    load = result["volume_mm3"] / float(height)
    assert err < 1e-8 and abs(result["energy_n_mm"] / energy - 1) < 1e-9
    assert abs(result["force_n"][2] / load - 1) < 1e-9
    assert (
        max(abs(v + w) for v, w in zip(result["force_n"], result["reaction_n"]))
        < 1e-9 * load
    )
    assert result["relative_residual"] < 1e-9
    wrong = json.loads(request.read_text())
    wrong["tractions"][0]["select"]["feature"] = "missing-face"
    refused = work / "invalid-traction.job.json"
    refused.write_text(json.dumps(wrong))
    saved = target.read_bytes()
    failure = subprocess.run(
        [cli, "solve-job", package, refused, target], capture_output=True
    )
    assert (
        failure.returncode
        and b"empty spline feature selection" in failure.stdout
        and target.read_bytes() == saved
    )
    report = dict(
        scalar_coefficients=len(values),
        max_coefficient_error_mm=err,
        energy_relative_defect=abs(result["energy_n_mm"] / energy - 1),
        load_n=result["force_n"],
        reaction_n=result["reaction_n"],
        free_residual=result["relative_residual"],
        unknown_feature_refused=True,
    )
    (work / "traction-audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(report)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--cli", type=Path, default=Path("target/release/solidcraft-trivariate")
    )
    p.add_argument("--work", type=Path, default=Path("build/thb-mechanics"))
    args = p.parse_args()
    verify(args.cli.resolve(), args.work)
