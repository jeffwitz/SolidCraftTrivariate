"""Independent small-system spectrum, rigid-vector and affine-strain checks."""

import argparse
import json
from pathlib import Path
import subprocess
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def run(cli, work):
    work.mkdir(exist_ok=True, parents=True)
    model = work / "model.json"
    package = work / "stiffness.json"
    subprocess.run(
        [cli, "build", ROOT / "tests/references/nurbs_strip.recipe.json", model],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [cli, "export-parent-stiffness", model, "12", package],
        check=True,
        capture_output=True,
    )
    data = json.loads(package.read_bytes())
    x = np.asarray(data["control_positions_mm"])
    n = len(x)
    K = np.zeros((3 * n, 3 * n))
    for i, (a, b) in enumerate(zip(data["offsets"], data["offsets"][1:])):
        for k in range(a, b):
            j = data["columns"][k]
            K[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] = data["values"][k]
    assert np.max(np.abs(K - K.T)) < 1e-12
    positive = np.linalg.eigvalsh(K)
    tol = positive[-1] * 1e-10
    assert min(positive) > -tol and sum(abs(positive) <= tol) == 6
    rigid = [np.tile(e, (n, 1)) for e in np.eye(3)] + [
        np.cross(e, x) for e in np.eye(3)
    ]
    residuals = [
        float(np.linalg.norm(K @ v.ravel()) / (positive[-1] * np.linalg.norm(v)))
        for v in rigid
    ]
    assert max(residuals) < 1e-10
    dofs = [3 * i + j for i in data["free_scalar_roots"] for j in range(3)]
    reduced = K[np.ix_(dofs, dofs)]
    np.linalg.cholesky(reduced)
    G = np.array(
        [[0.013, -0.004, 0.002], [0.001, 0.009, -0.003], [-0.002, 0.005, 0.006]]
    )
    u = (x @ G.T + np.array([0.7, -0.2, 0.3])).ravel()
    strain = (G + G.T) / 2
    expected = (np.sum(strain**2) + 0.75 * np.trace(strain) ** 2) * data["volume_mm3"]
    energy = 0.5 * u @ K @ u
    assert abs(energy - expected) < 1e-10 * expected
    report = dict(
        matrix_dimension=3 * n,
        rigid_relative_residuals=residuals,
        scaled_zero_threshold=1e-10,
        zero_modes=int(sum(abs(positive) <= tol)),
        first_positive_over_max=float(positive[6] / positive[-1]),
        clamped_min_eigenvalue=float(np.linalg.eigvalsh(reduced)[0]),
        affine_energy_relative_defect=abs(energy - expected) / expected,
    )
    (work / "spectrum.json").write_text(json.dumps(report, indent=2) + "\n")
    print(report)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate"
    )
    p.add_argument("--work", type=Path, default=ROOT / "build/parent-spectrum")
    args = p.parse_args()
    run(args.cli.resolve(), args.work)
