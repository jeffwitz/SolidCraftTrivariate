"""Collect EF4 refinement attempts, calibrated C3D10 solves and field comparisons."""
import argparse
import json
from pathlib import Path


def read(path):
    return json.loads(path.read_text())


def summarize(root):
    attempts = {name: read(root/name/"metrics.json") for name in
                ("tet-fine", "tet-iterative", "tet-corrected", "tet-corrected-auto")}
    fields = {}
    fine_root = root/"tet-corrected-auto"
    fine_cases = {(r["radius_x_ticks"],r["size_mm"]):r for r in attempts["tet-corrected-auto"] if r["stage"]=="fine"}
    for ticks in (64,56):
        for size in (1.,.9):
            run=fine_root/f"x{ticks}-h{size:g}"
            field=read(run/"field.json")
            timing=run/"field.resources"
            raw=timing.read_text().split() if timing.exists() else []
            status=field["qualification"]
            selected=field["levels"][-1]
            case=fine_cases[ticks,size]
            fields[f"x{ticks}-h{size:g}"]={
                "solve_qualified":case["qualified"],
                "field_qualified":all(status[k] for k in ("stable_distances","all_points_located_once","energy_recovery")),
                "comparison_wall_s":float(raw[0]) if raw else field["elapsed_s"],
                "comparison_hwm_kib":int(raw[1]) if raw else None,
                "finest_points":selected["total_points"],
                "field":field,
            }
    return {
        "schema":"bracket_energy_field_reference_refinement_v1",
        "units":{"energy":"N mm","distance":"dimensionless","length":"mm","memory":"KiB"},
        "scope":"Distances between discrete fields, divided by the reference field norm. Randomized integration dispersion is distinct from discretization uncertainty.",
        "sampling":"Eight scrambled Sobol replicates on each of 36 original rational charts. The same physical points evaluate all interpolants; each curved tetrahedron uses its actual quadratic inverse.",
        "resource_scope":"External child-process-tree peaks sampled during solver and mesher calls; comparison high-water marks from GNU time. Calibrations and failed attempts are retained. Wall times include concurrent campaign activity.",
        "budget_per_solver":{"seconds":1800,"rss_kib":6291456},
        "attempts":attempts,"comparisons":fields,
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path);parser.add_argument("output",type=Path)
    a=parser.parse_args();result=summarize(a.root)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    if not all(x["field_qualified"] and x["solve_qualified"] for x in result["comparisons"].values()):
        raise SystemExit(1)
