"""Verify the native h-sweep against exact continuum integrals and exact trace ranks."""

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from check_curved_pocket_fields import EVIDENCE, ROOT, qualify, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, help="save observations after qualification")
    args = parser.parse_args()
    continuum = json.loads((EVIDENCE/"continuum.json").read_text())["quintic"]
    traces = json.loads((ROOT/"tests/references/curved_pocket_trace_richness.json").read_text())
    dimensions = {(r["p"],r["spans_per_chart_axis"]):
                  (r["whole_five_patch_dofs_after_value_tangent_gluing"],
                   r["whole_five_patch_scalar_dofs_after_physical_C1"]) for r in traces["results"]}
    energy = float(Fraction(continuum["integrals"]["twice_strain_energy"]))
    data = ROOT/"build/curved-pocket-vector"
    rows = []
    for p in (3,4):
        for n in (1,2,4,8):
            for mode in ("none","all"):
                last = None
                for q in (8,10):
                    r = json.loads((data/f"p{p}-n{n}-{mode}-q{q}.json").read_text())
                    checks, defect = qualify(r, energy)
                    assert all(checks.values()), (p,n,mode,q,checks,defect)
                    assert r["free_vector_dofs"]==3*dimensions[p,n][int(mode=="all")]
                    assert r["interpolation_jet_residual"] <= 2e-10
                    for group in ("bottom","top","floor","wall"):
                        exact = continuum["force_and_body_moment"][group]
                        for a,b in zip(r[group+"_traction_resultant_N"],exact,strict=True):
                            assert math.isclose(a,float(Fraction(b)),rel_tol=1e-8,abs_tol=1e-8)
                    if last is not None and n==8:
                        assert math.isclose(last,r["relative_energy_norm"],rel_tol=1e-4)
                    last = r["relative_energy_norm"]
                    keys = ("free_vector_dofs","c1_constraint_rank","relative_energy_norm","relative_displacement_L2",
                            "relative_stress_L2","relative_residual","virtual_work_relative_defect",
                            "energy_N_mm","interpolation_jet_residual","assembly_ms","solve_ms")
                    rows.append({"p":p,"n":n,"joins":mode,"quadrature":q,
                                 **{k:r[k] for k in keys},"galerkin_energy_defect":defect})
    controls = json.loads((EVIDENCE/"continuum.json").read_text())
    for field, ref in (("affine-p3","affine_shifted"),("quadratic-p4","quadratic_mixed")):
        truth = controls[ref]
        for join in ("c0","c1"):
            r = json.loads((data/"controls"/f"{field}-{join}-q10.json").read_text())
            assert r["relative_energy_norm"] < 1e-8
            assert math.isclose(2*r["energy_N_mm"],float(Fraction(truth["integrals"]["twice_strain_energy"])),rel_tol=1e-8)
            for key, item in (("body","body_force_resultant_N"),("body_moment","body_force_moment_N_mm"),
                              *((g,g+"_traction_resultant_N") for g in ("floor","bottom","wall","top"))):
                for a,b in zip(r[item],truth["force_and_body_moment"][key],strict=True):
                    assert math.isclose(a,float(Fraction(b)),rel_tol=1e-8,abs_tol=1e-8),(field,item,a,b)
    for mode in ("all","none"):
        a = json.loads((data/f"p4-n1-{mode}-q10.json").read_text())
        b = json.loads((data/f"p4-n1-{mode}-primary-q12.json").read_text())
        assert math.isclose(a["relative_energy_norm"],b["relative_energy_norm"],rel_tol=1e-7)
    if args.record:
        write(args.record,rows)
    else:
        expected = json.loads((EVIDENCE/"primary.json").read_text())
        assert len(rows)==len(expected)
        for r,e in zip(rows,expected,strict=True):
            for k in ("p","n","joins","quadrature","free_vector_dofs","c1_constraint_rank"):
                assert r[k]==e[k]
            for k in ("relative_energy_norm","relative_displacement_L2","energy_N_mm"):
                assert math.isclose(r[k],e[k],rel_tol=2e-6,abs_tol=5e-10),(k,r,e)
    print(f"{len(rows)} native solves agree with the physical energy, boundary resultants and section dimensions.")


if __name__=="__main__":
    main()
