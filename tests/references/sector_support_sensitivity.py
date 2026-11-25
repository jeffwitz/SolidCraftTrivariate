"""Compare variable-weight sector support and reaction jets to cold, mounted models."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

from sector_job_sensitivity import call, job, recipe, solve, write


def offsets(solution):
    return np.array([n / d for component in solution["component_lifts"]
                     for row in component["lines"] for n, d in [row["offset"]]])


def rates_of_offsets(rate, count):
    result = np.zeros((3, count))
    for c, lines in enumerate(rate["component_lift_rates"]):
        for row in lines:
            result[c, row["dof"]] = row["offset_rate"]
    return result.ravel()


def quantities(solution):
    return [np.asarray(solution[key]) for key in (
        "traction_work_n_mm", "support_work_n_mm", "strain_energy_n_mm",
        "displacement_mm_xyz", "force_n_xyz")] + [offsets(solution)]


def compare(values, reference):
    return [float(np.linalg.norm(a-b) / max(1e-8, 1e-5*np.linalg.norm(b)))
            for a, b in zip(values, reference)]


def mounted_job(t, version):
    load = dict(select=dict(feature="inner"), pressure_mpa=0.3)
    specification = job(load, 12)
    specification["schema"] = f"solidcraft_mechanical_job_v{version}"
    root = "first" if t > 2.5 else "whole"
    specification["plan"].update(max_depth=2, max_cells=64, steps=[
        dict(refine=[dict(root=root, level=0, index=[0, 0, 0])]),
        dict(refine=[dict(root=root, level=1, index=[0, 0, 0])])])
    specification["supports"] = [dict(select=dict(feature="start"), components=[0, 1, 2],
        affine_mm=dict(denominator=4096, numerators=[[3, 4, 6, 0], [1, -2, 3, 1], [-2, 3, -4, 5]]))]
    return specification


def run(cli, work):
    work.mkdir(exist_ok=True, parents=True)
    cases = []
    for t, version in [(Fraction(1, 2), 2), (Fraction(3), 3)]:
        folder = work / f"v{version}"
        folder.mkdir(exist_ok=True)
        j = mounted_job(t, version)
        base = solve(cli, folder, t, j)
        direction = folder / "direction.json"
        write(direction, dict(schema="solidcraft_sector_angle_velocity_v3"))
        call(cli, "differentiate-job", folder/"model.json", folder/"job.json", direction, folder/"jet.json")
        jet = json.loads((folder/"jet.json").read_text())
        target = [np.asarray(jet[key]) for key in ("applied_work_rate_tangent_n_mm",
            "support_work_rate_n_mm", "strain_energy_rate_tangent_n_mm",
            "physical_displacement_rate_mm_xyz", "nodal_force_rate_n_xyz")]
        target.append(rates_of_offsets(jet, base["scalar_dofs"]))
        previous = previous_fourth = None
        rows = []
        for denominator in [16, 32, 64]:
            h = Fraction(1, denominator)
            states = [solve(cli, folder, t+s*h, j) for s in [-1, 1]]
            for state in states:
                assert [l["free_roots"] for l in state["component_lifts"]] == [
                    l["free_roots"] for l in base["component_lifts"]]
                assert [(b["id"], b["dofs"]) for b in state["analysis_blocks"]] == [
                    (b["id"], b["dofs"]) for b in base["analysis_blocks"]]
            minus, plus = map(quantities, states)
            current = [(b-a)*denominator/2 for a,b in zip(minus, plus)]
            row = dict(step=1/denominator, centered=compare(current, target))
            if previous is not None:
                fourth = [(4*a-b)/3 for a,b in zip(current, previous)]
                row["fourth"] = compare(fourth, target)
                if previous_fourth is not None:
                    row["sixth"] = compare([(16*a-b)/15 for a,b in zip(fourth, previous_fourth)], target)
                previous_fourth = fourth
            rows.append(row)
            previous = current
        checks = dict(cold_models=max(rows[-1]["sixth"]) <= 1,
                      nonzero_offsets=bool(np.max(np.abs(target[-1])) > 1e-5),
                      nonzero_lift_terms=any(abs(v) > 1e-6 for comp in jet["component_lift_rates"]
                                             for row in comp for _,v in row["terms"]),
                      work_balance=abs(jet["work_rate_closure_n_mm"]) <= 1e-9)
        for label in ["applied_work", "strain_energy"]:
            a, b = [jet[f"{label}_rate_{method}_n_mm"] for method in ["tangent", "adjoint"]]
            checks[label+"_adjoint"] = abs(a-b) <= max(1e-9, 1e-8*abs(a))
        cases.append(dict(recipe=recipe(t), job=j, checks=checks, history=rows,
            work_rates_n_mm=[float(x) for x in target[:3]],
            max_offset_rate_mm=float(np.max(np.abs(target[-1]))), passed=all(checks.values())))
        print(t, checks, flush=True)
    refusals = []
    for label in ["incompatible_laws", "fixed_spatial_plane"]:
        j = mounted_job(Fraction(1,2),3)
        if label == "incompatible_laws":
            j["supports"].append(dict(select=dict(feature="end"), components=[1]))
        else:
            j["supports"][0]["select"] = dict(plane=dict(axis=2,coordinate_mm=0))
        write(work/"recipe.json", recipe(Fraction(1,2)))
        call(cli, "build", work/"recipe.json", work/"model.json")
        write(work/"job.json", j)
        out = work/"protected.json"
        out.write_bytes(b"previous")
        p = subprocess.run([cli,"differentiate-job",work/"model.json",work/"job.json",direction,out],
                           capture_output=True,text=True)
        refusals.append(dict(name=label, diagnostic=p.stdout.strip(),
                             passed=p.returncode!=0 and out.read_bytes()==b"previous"))
    return dict(schema="sector_support_sensitivity_oracle_v1", cases=cases, refusals=refusals,
                observables=["applied_work","support_work","strain_energy","field","nodal_force","lift_offsets"],
                executable_sha256=hashlib.sha256(cli.read_bytes()).hexdigest(),
                qualified=all(v["passed"] for v in cases+refusals))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("cli","work","output"):
        p.add_argument(name,type=Path)
    a = p.parse_args()
    data = run(a.cli.resolve(), a.work)
    write(a.output, data)
    raise SystemExit(int(not data["qualified"]))
