"""Compare complete homogeneous chart jets and physical field gradients on curved solids."""

import argparse
import json
import subprocess
from pathlib import Path

from audit_fieldcraft_parent_support import audit_model
from check_inherited_joins import verify


def mechanics(cli, work, name, model, matched):
    results = {}
    geometry = json.loads(model.read_text())
    xs = [p[0] / p[3] for b in geometry["blocks"] for p in b["controls"]]
    for mode in ("c0", "axial", "matched"):
        path = matched
        if mode != "matched":
            plan = work / f"{name}.{mode}.plan.json"
            plan.write_text(
                json.dumps(
                    dict(
                        schema="inherited_field_plan_v1",
                        axial_degree=2,
                        axial_continuity=int(mode == "axial"),
                        restore_strip=False,
                        transverse_interfaces=[],
                    )
                )
            )
            path = work / f"{name}.{mode}.json"
            subprocess.run(
                [cli, "export-inherited-space", model, plan, path],
                check=True,
                capture_output=True,
            )
        history = []
        for q in (8, 12, 16, 20, 24):
            job = work / f"{name}.q{q}.job.json"
            job.write_text(
                json.dumps(
                    dict(
                        schema="section_quadratic_mms_v1",
                        quadrature=q,
                        x_origin_mm=min(xs),
                        x_scale_mm=max(xs) - min(xs),
                    )
                )
            )
            out = work / f"{name}.{mode}.q{q}.json"
            subprocess.run(
                [cli, "solve-inherited-space", path, job, out],
                check=True,
                capture_output=True,
            )
            r = json.loads(out.read_text())
            stable = history and all(
                abs(r[k] - history[-1][k]) <= abs(r[k]) * 1e-8
                for k in ("energy_N_mm", "mean_displacement_z_mm")
            )
            history.append(r)
            if (
                stable
                and r["virtual_work_relative_defect"] < 1e-8
                and r["relative_residual"] < 1e-9
            ):
                break
        else:
            raise RuntimeError(
                f"{name} {mode}: prescribed quadrature or equilibrium budget exhausted"
            )
        print(
            name, mode, q, r["free_vector_dofs"], r["relative_energy_norm"], flush=True
        )
        results[mode] = history
    return results


def run(cli, work):
    work.mkdir(parents=True, exist_ok=True)
    results = {}
    cases = {"g2_cross_sketches": 4, "s3_face_pocket": 16, "nurbs_multiloop": 26}
    for name, expected in cases.items():
        model = work / (name + ".model.json")
        subprocess.run(
            [cli, "build", f"tests/references/{name}.recipe.json", model],
            check=True,
            capture_output=True,
        )
        jets = audit_model(model)
        oracle = {}
        for row in jets["interfaces"]:
            if any(row[k][1].startswith("w") for k in ("first", "second")):
                continue
            witness = row["jet_evidence"]
            if witness["jet_scaled_ratios"][0]["match"]:
                oracle[(tuple(row["first"]), tuple(row["second"]))] = witness[
                    "candidate_interval_length_ratio_first_to_second"
                ]
        assert len(oracle) == expected
        path = work / (name + ".field.json")
        subprocess.run(
            [
                cli,
                "export-inherited-space",
                model,
                "tests/references/matched_chart_jets.plan.json",
                path,
            ],
            check=True,
            capture_output=True,
        )
        package = json.loads(path.read_text())
        f = package["field"]
        matched = {
            (tuple(r["first"]), tuple(r["second"])): r["normal_interval_ratio"]
            for r in f["interfaces"]
            if r["role"] == "distinct_support_join"
        }
        assert matched == oracle
        physical = verify(package)
        assert physical["nonzero_polynomials"] == 0 and physical["joins"] == expected
        results[name] = {
            k: f[k]
            for k in (
                "scalar_dofs",
                "inherited_join_count",
                "transverse_c1_join_count",
                "retained_c0_join_count",
            )
        }
        results[name].update(
            complete_faces=jets["interface_count"], physical_gradient_oracle=physical
        )
        print(name, results[name], flush=True)
        results[name]["mechanics"] = mechanics(cli, work, name, model, path)
    (work / "chart-jets.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--cli", type=Path, default=Path("target/release/solidcraft-trivariate")
    )
    p.add_argument("--work", type=Path, default=Path("build/curved-supports"))
    args = p.parse_args()
    run(args.cli.resolve(), args.work)
