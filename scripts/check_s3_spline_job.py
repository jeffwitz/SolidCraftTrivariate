"""Verify the physical displacement and reaction for a source-tagged S3 stress load."""

import argparse
import json
import subprocess
from pathlib import Path


def run(cli, work):
    field = work / "s3_face_pocket.field.json"
    pkg = json.loads(field.read_text())
    out = work / "s3-stress.json"
    subprocess.run(
        [cli, "solve-job", field, "tests/references/s3_inherited_stress.job.json", out],
        check=True,
    )
    data = json.loads(out.read_text())
    displacement = data["coefficients_xyz_mm"]
    blocks = {b["id"]: b for b in pkg["geometry"]["blocks"]}
    max_error = 0.0
    count = 0
    for cell in pkg["field"]["cells"]:
        p = cell["degrees"][2]
        net = blocks[cell["id"]]["controls"]
        for i, row in enumerate(cell["rational_extraction"]):
            section, k = divmod(i, p + 1)
            low, high = net[2 * section : 2 * section + 2]
            z = ((p - k) * low[2] / low[3] + k * high[2] / high[3]) / p
            u = [
                sum(int(a) / int(b) * displacement[j][c] for j, (a, b) in row)
                for c in range(3)
            ]
            max_error = max(max_error, abs(u[0]), abs(u[1]), abs(u[2] - z / 10))
            count += 1
    reference = json.loads(Path("docs/results/g2-face-pocket-iga.json").read_text())
    volume = reference["volume_map_audit"]["expected_s3_volume_mm3"]
    expected = volume / 20
    rel = abs(data["energy_n_mm"] / expected - 1)
    balance = max(
        abs(a + b) for a, b in zip(data["force_n"], data["reaction_n"])
    ) / abs(data["force_n"][2])
    assert (
        max_error < 1e-7
        and rel < 1e-7
        and balance < 1e-9
        and data["relative_residual"] < 1e-9
    )
    report = dict(
        local_coefficients_checked=count,
        source_reference_volume_mm3=volume,
        native_strain_energy_n_mm=data["energy_n_mm"],
        max_coefficient_error_mm=max_error,
        relative_energy_defect=rel,
        relative_force_balance=balance,
        force_n=data["force_n"],
        reaction_n=data["reaction_n"],
        free_residual=data["relative_residual"],
    )
    (work / "s3-stress-audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(report)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--cli", type=Path, default=Path("target/release/solidcraft-trivariate")
    )
    p.add_argument("--work", type=Path, default=Path("build/curved-supports"))
    args = p.parse_args()
    run(args.cli.resolve(), args.work)
