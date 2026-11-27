"""Replay the fixed-geometry section-loading comparison at common quadrature."""

import argparse
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "cylinder": ("docs/site/box-cylinder/tall/cylinder", 20, 10),
    "pocket": ("docs/site/data/plate_blind_pocket/plate_blind_pocket", 30, 60),
}


def run(cli, work, quick=False):
    results = {}
    for name, (prefix, center, width) in CASES.items():
        directory = work / name
        directory.mkdir(parents=True, exist_ok=True)
        model = directory / "source.json"
        subprocess.run(
            [cli, "build", ROOT / (prefix + ".recipe.json"), model],
            check=True,
            capture_output=True,
        )
        assert json.loads(model.read_bytes()) == json.loads(
            (ROOT / (prefix + ".model.json")).read_bytes()
        )
        geometry = json.loads(model.read_bytes())
        pairs = [
            [f["first"][0], f["second"][0]]
            for f in geometry["topology"]["interfaces"]
            if f["first"][1][0] != "w"
        ]
        # A single physical rail includes its restrictions on every axial layer.
        rail = [
            p
            for p in pairs
            if (
                name == "cylinder"
                and all(i.endswith(("_core", "_sector_west")) for i in p)
            )
            or (
                name == "pocket"
                and all("east" in i or "south" in i for i in p)
                and any("east" in i for i in p)
                and any("south" in i for i in p)
            )
        ]
        assert rail
        modes = [
            ("C0", 0, [], 0),
            ("axial", 1, [], 0),
            ("one_rail", 1, rail, 0),
            ("all_lateral", 1, pairs, 0),
        ]
        if not quick:
            modes.extend(
                [
                    ("cubic_C0", 0, [], 3),
                    ("cubic_axial", 1, [], 3),
                    ("cubic_rail", 1, rail, 3),
                    ("cubic_all", 1, pairs, 3),
                ]
            )
        records = {}
        for mode, continuity, joins, plane in modes:
            plan = dict(
                schema="inherited_field_plan_v1",
                axial_degree=2,
                axial_continuity=continuity,
                restore_strip=False,
                transverse_interfaces=joins,
                plane_degree=plane,
            )
            request = directory / (mode + ".plan.json")
            request.write_text(json.dumps(plan))
            package = directory / (mode + ".field.json")
            result = subprocess.run(
                [cli, "export-inherited-space", model, request, package],
                text=True,
                capture_output=True,
            )
            if result.returncode:
                records[mode] = dict(
                    status="compile_refused",
                    reason=(result.stdout + result.stderr).strip(),
                )
                print(name, mode, records[mode], flush=True)
                continue
            runs = []
            for q in [12] if quick else [12, 16, 20, 24, 32]:
                job = directory / (mode + ".job.json")
                answer = directory / f"{mode}.q{q}.json"
                job.write_text(
                    json.dumps(
                        dict(
                            schema="section_quadratic_mms_v1",
                            quadrature=q,
                            x_origin_mm=center,
                            x_scale_mm=width,
                        )
                    )
                )
                completed = subprocess.run(
                    [cli, "solve-inherited-space", package, job, answer],
                    text=True,
                    capture_output=True,
                )
                if completed.returncode:
                    runs.append(
                        dict(
                            status="solve_refused",
                            reason=(completed.stdout + completed.stderr).strip(),
                        )
                    )
                else:
                    runs.append(json.loads(answer.read_bytes()))
                print(
                    name,
                    mode,
                    q,
                    runs[-1].get("free_vector_dofs"),
                    runs[-1].get("relative_energy_norm"),
                    runs[-1].get("reason"),
                    flush=True,
                )
                if "energy_N_mm" not in runs[-1]:
                    raise RuntimeError(runs[-1])
                if len(runs) >= 2:
                    previous, current = runs[-2:]
                    delta = {
                        key: abs(current[key] - previous[key])
                        / max(abs(current[key]), 1e-30)
                        for key in ("energy_N_mm", "mean_displacement_z_mm")
                    }
                    current["relative_quadrature_change"] = delta
                    if max(delta.values()) <= 1e-8:
                        break
            final = runs[-1]
            if final["relative_residual"] > 1e-9:
                raise RuntimeError(("native equilibrium", name, mode, final))
            final["quadrature_status"] = (
                "single_order"
                if quick
                else (
                    "stabilized"
                    if max(final["relative_quadrature_change"].values()) <= 1e-8
                    else "budget_exhausted"
                )
            )
            records[mode] = runs
        results[name] = records
        (work / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    for name, modes in results.items():
        for mode, runs in modes.items():
            if not isinstance(runs, list):
                raise RuntimeError((name, mode, runs))
            last = runs[-1]
            if not quick and last["quadrature_status"] != "stabilized":
                raise RuntimeError(("quadrature budget", name, mode, last))
            if last["virtual_work_relative_defect"] > 1e-8:
                raise RuntimeError(("external load closure", name, mode, last))
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/solidcraft-trivariate"
    )
    p.add_argument("--work", type=Path, default=ROOT / "build/inherited-mechanics")
    p.add_argument("--quick", action="store_true")
    args = p.parse_args()
    run(args.binary.resolve(), args.work, args.quick)
