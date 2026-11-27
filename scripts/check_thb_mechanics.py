"""Replay the prepared and locally refined spaces at common physical quadrature."""

import argparse
import json
import subprocess
from pathlib import Path


def main(cli, work):
    work.mkdir(parents=True, exist_ok=True)
    model, package = (work / f"{name}.json" for name in ("model", "package"))

    def run(*args):
        return subprocess.run([cli, *args], check=True, capture_output=True, text=True)

    run("build", "tests/references/nurbs_strip.recipe.json", model)
    run("refine-thb-strip", model, "tests/references/thb_interior.plan.json", package)
    # A forged successful package and a bad refinement request must both leave
    # an existing destination untouched.
    forged = json.loads(package.read_text())
    forged["field"]["boundary_sources"][0]["upper"] = "forged"
    invalid = work / "forged.json"
    invalid.write_text(json.dumps(forged))
    sentinel = work / "refusal-target.json"
    sentinel.write_text("preserve existing output\n")
    bad_plan = work / "boundary.plan.json"
    bad_plan.write_text(
        json.dumps(
            dict(schema="interior_thb_refinement_v1", start=[0, 1, 1], end=[6, 3, 3])
        )
    )
    for request in [
        ("export-thb-extraction", invalid, sentinel),
        ("refine-thb-strip", model, bad_plan, sentinel),
    ]:
        result = subprocess.run([cli, *request], capture_output=True, text=True)
        assert (
            result.returncode and sentinel.read_text() == "preserve existing output\n"
        )
    galerkin = json.loads(run("check-thb-galerkin", package, "4").stdout)
    results = {}
    for name, center in [
        ("interior", [13 / 6, 23 / 12, 5 / 2]),
        ("off_center", [0.42, 0.86, 1.25]),
    ]:
        levels = []
        for order in (4, 6, 8, 10, 12, 16):
            job = work / f"{name}-q{order}.job.json"
            job.write_text(
                json.dumps(
                    dict(
                        schema="local_gaussian_mms_v1",
                        center_mm=center,
                        lengths_mm=[0.4, 0.35, 0.55],
                        quadrature=order,
                    )
                )
            )
            current = {}
            for mode in ("prepared", "local"):
                target = work / f"{name}-{mode}-q{order}.json"
                run("solve-thb-strip", package, mode, job, target)
                current[mode] = json.loads(target.read_text())
            stable = bool(levels) and all(
                abs(current[m][key] - levels[-1][m][key]) <= 1e-8 * abs(current[m][key])
                for m in current
                for key in ("energy_N_mm", "mean_displacement_z_mm")
            )
            levels.append(current)
            print(
                name,
                order,
                {m: v["relative_energy_norm"] for m, v in current.items()},
                flush=True,
            )
            if stable and all(
                v["virtual_work_relative_defect"] < 1e-8
                and v["relative_residual"] < 1e-9
                for v in current.values()
            ):
                break
        else:
            raise RuntimeError(f"{name}: quadrature or equilibrium budget exhausted")
        assert (
            current["local"]["relative_energy_norm"]
            <= current["prepared"]["relative_energy_norm"] + 1e-9
        )
        results[name] = levels
    field = json.loads(package.read_text())["field"]
    report = {
        k: v
        for k, v in field.items()
        if k not in ("integration_sources", "boundary_sources")
    }
    report.update(galerkin=galerkin, loads=results)
    (work / "results.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cli", type=Path, default=Path("target/release/solidcraft-trivariate")
    )
    parser.add_argument("--work", type=Path, default=Path("build/thb-mechanics"))
    args = parser.parse_args()
    main(args.cli.resolve(), args.work)
