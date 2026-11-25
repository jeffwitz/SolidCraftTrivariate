"""Solve the h2 fillet ladders, preserving every time limit and field checkpoint."""
import argparse
import json
import os
from pathlib import Path
import time

from bracket_cost_fields import read, save
from local_fillet_cost import assess
from refine_tetra_fields import measure


def run(executable, work, report):
    root = Path(__file__).resolve().parents[2]
    manifest = read(root/"tests/references/bracket_cost_corpus.json")
    baseline = read(root/"docs/results/bracket-cost-20261125.json")
    cases = manifest["mandatory_radii_sixteenths_mm"]
    refs = {tuple(row["radii_ticks"]): row for row in baseline["cases"]}
    work.mkdir(exist_ok=True, parents=True)
    summary = dict(schema="bracket_h2_field_tiers_v1", physical_residual_limit=1e-10,
                   internal_pcg_target=5e-11, execution_order="depth, radius, then q8 and q6",
                   states=[], mandatory_radii=cases)
    by_case = {tuple(case): [] for case in cases}
    for step in range(4):
        for rx, ry in cases:
            folder = work/f"x{rx}-y{ry}"
            folder.mkdir(exist_ok=True)
            states = by_case[rx, ry]
            row = dict(radii_ticks=[rx, ry], step=step, orders=[])
            if states and states[-1]["status"] != "completed":
                row["status"] = "not_started_after_incomplete_state"
            else:
                values = {}
                spent = 0.
                for q in (8, 6):
                    out = folder/f"s{step}-q{q}.json"
                    checkpoint = folder/f"s{step}.bin"
                    cmd = list(map(str, [executable,"--field-tier",rx,ry,step,q,out,checkpoint]))
                    measurement = measure(cmd, dict(os.environ), folder/f"s{step}-q{q}.log", seconds=max(0.,1800-spent))
                    spent += measurement["elapsed_s"]
                    row["orders"].append(dict(order=q, **measurement))
                    if measurement["status"] != "completed":
                        row["status"] = measurement["status"]
                        break
                    trial = read(out)
                    # Earlier diagnostic builds labelled both fields q8; the
                    # stored quadrature order is independently checked below.
                    value = trial.get(f"q{q}", trial.get("q8"))
                    if value["q"] != q: raise ValueError("incorrect quadrature recorded")
                    values[q] = value
                    row["leaves"] = trial["leaves"]
                    row["partition_cumulative_s"] = trial["uniform_partition_s"] + trial["local_partition_s"]
                    row["constraints_s"] = trial["constraints_s"]
                    row["geometry_build_and_audit_s"] = trial["geometry_build_and_audit_s"]
                if len(values) == 2:
                    row["status"] = "completed"
                    row["q6_q8"] = [values[6], values[8]]
            states.append(row)
            summary["states"].append(row)
            completed = [s for s in states if s["status"]=="completed"]
            if completed:
                ref = refs[rx, ry]
                tetra = next(t for t in ref["runs"] if t["method"]=="tetra" and t["resolution"]==1.25)
                checks = dict(status=ref["reference_status"], point_tip_y_mm=tetra["result"]["point_tip_y_mm"],
                              von_mises_mpa=tetra["probe"]["von_mises_mpa"])
                assess(dict(paths=[dict(states=completed)], geometry_build_and_audit_s=completed[0]["geometry_build_and_audit_s"]),
                       [rx,ry], checks, manifest)
            save(report, summary)
            print(rx,ry,step,row["status"],row.get("qualified"), row["orders"],flush=True)
    summary["finished"] = True
    save(report, summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "work", "report"): p.add_argument(name,type=Path)
    args = p.parse_args()
    run(args.executable.resolve(), args.work, args.report)
