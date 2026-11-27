"""Compare physical jet constraints on two unchanged polygon benchmark solids."""

import argparse
import json
from pathlib import Path
import subprocess

from census_transverse_transitions import classify
from check_inherited_joins import verify

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "seeded_convex_03": (
        "8af707b273a37a6081cf1e34aa90ee15256c86fdfc0abee7418cfa28ac5b4cd7",
        2000,
    ),
    "blind_triangle_shallow": (
        "537f431f5c393ff5af15e97d90bc67d0d16fedf112489ce05292017c12bc25a9",
        64,
    ),
}


def run(cli, sources, work):
    work.mkdir(parents=True, exist_ok=True)
    results = {}
    for name, (fingerprint, width) in CASES.items():
        model_path = sources / (fingerprint + ".model.json")
        model = json.loads(model_path.read_bytes())
        oracle_records = list(classify(model))
        native_report = work / name / "transverse-transitions.json"
        native_report.parent.mkdir(parents=True, exist_ok=True)
        classified = subprocess.run(
            [cli, "classify-transverse-transitions", model_path, native_report],
            capture_output=True,
            text=True,
        )
        if classified.returncode:
            raise RuntimeError(
                f"native transverse classification failed for {name}: "
                f"{classified.stdout}{classified.stderr}"
            )
        records = json.loads(native_report.read_bytes())
        if records != oracle_records:
            raise AssertionError(f"native/Fraction transition mismatch for {name}")
        qualifying = {"constant_transition", "affine_transition"}
        eligible = [
            r for r in records if {r["status"], r.get("reversed_status")} & qualifying
        ]
        candidate_pairs = [[face[0] for face in r["faces"]] for r in eligible]
        outcomes = {}
        for mode, continuity, pairs in (
            ("c0", 0, []),
            ("axial", 1, []),
            ("affine_jet_faces", 1, candidate_pairs),
        ):
            folder = work / name / mode
            folder.mkdir(parents=True, exist_ok=True)
            plan = folder / "plan.json"
            plan.write_text(
                json.dumps(
                    dict(
                        schema="inherited_field_plan_v1",
                        axial_degree=2,
                        axial_continuity=continuity,
                        restore_strip=False,
                        transverse_interfaces=pairs,
                        plane_degree=2,
                    )
                )
            )
            package = folder / "field.json"
            cmd = subprocess.run(
                [cli, "export-inherited-space", model_path, plan, package],
                capture_output=True,
                text=True,
            )
            if cmd.returncode:
                outcomes[mode] = dict(
                    status="compile_refused", reason=cmd.stdout + cmd.stderr
                )
                continue
            data = json.loads(package.read_bytes())
            proof = verify(data, check_retained=False) if pairs else None
            if proof:
                assert proof["nonzero_polynomials"] == 0 and proof["joins"] == len(
                    pairs
                )
            runs = []
            for order in (5, 8, 12, 16, 20, 24, 32):
                job, answer = folder / "job.json", folder / f"q{order}.json"
                job.write_text(
                    json.dumps(
                        dict(
                            schema="section_cubic_mms_v1",
                            x_origin_mm=-width,
                            x_scale_mm=width,
                            quadrature=order,
                        )
                    )
                )
                solved = subprocess.run(
                    [cli, "solve-inherited-space", package, job, answer],
                    capture_output=True,
                    text=True,
                )
                if solved.returncode:
                    runs.append(
                        dict(
                            status="solve_refused", reason=solved.stdout + solved.stderr
                        )
                    )
                    break
                runs.append(json.loads(answer.read_bytes()))
                if len(runs) >= 2:
                    prev, current = runs[-2:]
                    delta = max(
                        abs(current[k] - prev[k]) / max(abs(current[k]), 1e-30)
                        for k in ("energy_N_mm", "mean_displacement_z_mm")
                    )
                    current["maximum_relative_quadrature_change"] = delta
                    if delta <= 1e-8:
                        break
            outcomes[mode] = dict(
                requested_joins=len(pairs), exact_jet_check=proof, runs=runs
            )
            print(
                name,
                mode,
                runs[-1].get("free_vector_dofs"),
                runs[-1].get("relative_energy_norm"),
                runs[-1].get("reason"),
                flush=True,
            )
        results[name] = dict(
            recipe_fingerprint=fingerprint,
            eligible_lateral_count=len(eligible),
            transition_oracle_parity="exact",
            outcomes=outcomes,
        )
    (work / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    for name, case in results.items():
        for mode, outcome in case["outcomes"].items():
            last = outcome.get("runs", [{}])[-1]
            if (
                last.get("maximum_relative_quadrature_change", 1) > 1e-8
                or last.get("virtual_work_relative_defect", 1) > 1e-8
                or last.get("relative_residual", 1) > 1e-9
            ):
                raise RuntimeError((name, mode, outcome))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate"
    )
    parser.add_argument(
        "--sources", type=Path, default=ROOT / "build/support-census/models"
    )
    parser.add_argument(
        "--work", type=Path, default=ROOT / "build/polygon-join-mechanics"
    )
    args = parser.parse_args()
    run(args.cli.resolve(), args.sources, args.work)
