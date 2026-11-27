"""Classify complete-parent admission and mounted axial fields on named native states."""

from pathlib import Path
import argparse, json, subprocess
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]


def attempt(cli, *args):
    result = subprocess.run([cli, *args], capture_output=True, text=True)
    return result.returncode, json.loads(result.stdout)


def run(cli, work, census):
    work.mkdir(exist_ok=True, parents=True)
    named = json.loads(census.read_bytes())["named_states"]
    rows = []
    for state in named["states"]:
        model = ROOT / state["source"]
        code, report = attempt(cli, "spline-eligibility", model)
        assert code == 0 and report["source_validity"] == "valid"
        item = dict(
            name=state["state_id"],
            source=state["source"],
            source_validity="valid",
            profile=report["conversion_status"],
            explanation=report.get("common_parent_obstruction", {}).get(
                "kind", report["reason"]
            ),
        )
        for continuity in (0, 1):
            plan = work / "plan.json"
            plan.write_text(
                json.dumps(
                    dict(
                        schema="inherited_field_plan_v1",
                        axial_degree=2,
                        axial_continuity=continuity,
                        restore_strip=False,
                        transverse_interfaces=[],
                    )
                )
            )
            output = work / (state["state_id"] + f".c{continuity}.json")
            code, observation = attempt(
                cli, "export-inherited-space", model, plan, output
            )
            if code == 0:
                field = json.loads(output.read_bytes())["field"]
                item[f"c{continuity}"] = {
                    k: field[k]
                    for k in (
                        "scalar_dofs",
                        "inherited_join_count",
                        "retained_c0_join_count",
                    )
                }
            else:
                item[f"c{continuity}"] = observation
        rows.append(item)
    for state in named["s2_crossed_union"]["states"]:
        rows.append(
            dict(
                name=state["state_id"],
                source=state["source"],
                source_validity=state["status"],
                profile="outside_whole_face_profile",
                partial_faces=state["partial_trace_count"],
                complete_faces=state["full_exact_trace_count"],
                polar_line_constraints=state["non_face_collapsed_polar_line_constraints"],
            )
        )
    assert len(rows) == 17 and rows[0]["profile"] == "eligible_single_complete_parent"
    report = dict(
        states=rows, profile_observations=dict(Counter(s["profile"] for s in rows))
    )
    (work / "named-supports.json").write_text(json.dumps(report, indent=2) + "\n")
    print(report["profile_observations"])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate"
    )
    p.add_argument("--work", type=Path, default=ROOT / "build/named-supports")
    p.add_argument(
        "--census", type=Path, default=ROOT / "build/support-census/census.json"
    )
    args = p.parse_args()
    run(args.cli.resolve(), args.work, args.census)
