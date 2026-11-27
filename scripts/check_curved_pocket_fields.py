"""Replay the historical physical load on the published pocket and its local section fields."""

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/site/data/plate_blind_pocket"
OUT = ROOT / "build/curved-pocket-fields"
EVIDENCE = ROOT / "docs/results/curved-pocket-fields"


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def qualify(row, twice_energy):
    # Galerkin orthogonality on a homogeneous bottom clamp.
    defect = abs(2 * row["energy_N_mm"] / twice_energy + row["relative_energy_norm"]**2 - 1)
    checks = {"volume": abs(row["volume_mm3"] / 48976 - 1) < 1e-10,
              "equilibrium": row["relative_residual"] < 1e-9,
              "load_work": row["virtual_work_relative_defect"] < 1e-7,
              "galerkin_energy": defect < 1e-8}
    return checks, defect


def replay_archived_space(model, twice_energy):
    """Repeat the earlier 180-DOF calculation with its original axial degree."""
    cli = ROOT / "target/release/solidcraft-trivariate"
    topology = json.loads(model.read_text())["topology"]
    faces = [[x["first"][0],x["second"][0]] for x in topology["interfaces"] if x["first"][1][0] != "w"]
    plan,job,field,answer = [OUT/name for name in ("archive-plan.json","archive-job.json","archive-field.json","archive-result.json")]
    write(plan,{"schema":"inherited_field_plan_v1","axial_degree":2,"axial_continuity":1,"restore_strip":False,
                "transverse_interfaces":faces,"plane_degree":3})
    write(job,{"schema":"section_quadratic_mms_v1","quadrature":16,"x_origin_mm":30,"x_scale_mm":60})
    subprocess.run([cli,"export-inherited-space",model,plan,field],check=True,capture_output=True)
    subprocess.run([cli,"solve-inherited-space",field,job,answer],check=True,capture_output=True)
    row=json.loads(answer.read_text())
    old=json.loads((ROOT/"docs/results/inherited-continuity/mechanical-evidence.json").read_text())["section_mms"]["pocket"]["cubic_all"][-1]
    assert all(qualify(row,twice_energy)[0].values()) and row["free_vector_dofs"]==old["free_vector_dofs"]
    for k in ("relative_energy_norm","relative_displacement_L2","energy_N_mm"):
        assert math.isclose(row[k],old[k],rel_tol=1e-7), (k,row,old)
    print("Replayed axial Q2 historical field:",row["free_vector_dofs"],row["relative_energy_norm"])


def run(binary, expected):
    model = SOURCE / "plate_blind_pocket.model.json"
    twice_energy = float(Fraction(json.loads((EVIDENCE / "continuum.json").read_text())
                                  ["historical_quadratic"]["integrals"]["twice_strain_energy"]))
    replay_archived_space(model,twice_energy)
    rows = []
    for p in (3, 4):
        for n in (1, 2, 4, 8):
            for continuity in ("C0", "C1"):
                for q in ((8, 10) if n == 8 else (8,)):
                    job = OUT / f"quadratic-q{q}.job.json"
                    write(job, {"schema": "section_quadratic_mms_v1", "quadrature": q,
                                "x_origin_mm": 30, "x_scale_mm": 60})
                    command = [binary, model, str(p), str(n), continuity, job]
                    start = time.perf_counter()
                    proc = subprocess.run(command, capture_output=True, text=True, check=False)
                    item = {"p": p, "n": n, "continuity": continuity, "quadrature": q}
                    if proc.returncode != 0:
                        item.update(status="failed", reason=proc.stderr.strip()[-400:])
                    else:
                        row = json.loads(proc.stdout)
                        checks, defect = qualify(row, twice_energy)
                        if p == 4:
                            checks["polynomial_reproduction"] = row["relative_energy_norm"] < 1e-8
                        item.update(row, checks=checks, galerkin_energy_defect=defect,
                                    status="valid" if all(checks.values()) else "failed")
                    item["wall_seconds"] = time.perf_counter() - start
                    rows.append(item)
                    print(p, n, continuity, q, item.get("relative_energy_norm"), item["status"], flush=True)
                    write(OUT / "historical.json", rows)
    for row in rows:
        if row["quadrature"] != 10 or row["status"] != "valid":
            continue
        old = next(r for r in rows if all(r[k] == row[k] for k in ("p", "n", "continuity"))
                   and r["quadrature"] == 8)
        if old["status"] != "valid":
            continue
        # Absolute check for the exactly representable Q4 control near roundoff.
        ok = math.isclose(old["relative_energy_norm"], row["relative_energy_norm"],
                          rel_tol=1e-4, abs_tol=1e-8 if row["p"] == 4 else 0.0)
        row["checks"]["replicate"] = ok
        if not ok:
            row["status"] = "failed"
    write(OUT / "historical.json", rows)
    if expected:
        target = json.loads(expected.read_text())
        assert len(target) == len(rows)
        for a,b in zip(target,rows,strict=True):
            for k in ("p","n","continuity","quadrature","free_vector_dofs","status"):
                assert a.get(k) == b.get(k), (k,a,b)
            if a["status"] == "valid":
                for k in ("relative_energy_norm","relative_displacement_L2","energy_N_mm"):
                    assert math.isclose(a[k], b[k], rel_tol=1e-6, abs_tol=5e-10), (k,a,b)
    return all(r["status"] == "valid" for r in rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "target/release/examples/curved_pocket_field")
    parser.add_argument("--compare", type=Path)
    args = parser.parse_args()
    raise SystemExit(0 if run(args.binary, args.compare) else 1)
