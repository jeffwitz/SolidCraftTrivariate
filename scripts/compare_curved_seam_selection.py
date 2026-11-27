"""Measure the mechanical effect of the two geometrically distinct seam families."""

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
import subprocess
import time

from check_curved_pocket_fields import EVIDENCE, ROOT, OUT, SOURCE, qualify, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "target/release/examples/curved_pocket_field")
    parser.add_argument("--compare", type=Path)
    args = parser.parse_args()
    model = SOURCE / "plate_blind_pocket.model.json"
    job = OUT / "quintic-q10.job.json"
    write(job, {"schema":"section_affine_power_mms_v1", "coefficients":[1,2,0], "offset":10,
                "scale":150, "degree":5, "quadrature":10})
    exact = float(Fraction(json.loads((EVIDENCE / "continuum.json").read_text())
                           ["quintic"]["integrals"]["twice_strain_energy"]))
    results = []
    for p in (3,4):
        for n in (4,8):
            for family in ("C0","radial","footprint","C1"):
                start = time.perf_counter()
                proc = subprocess.run([args.binary,model,str(p),str(n),family,job],
                                      check=True,capture_output=True,text=True)
                r = json.loads(proc.stdout)
                checks, defect = qualify(r,exact)
                assert all(checks.values()), (p,n,family,checks,r)
                results.append({**r,"p":p,"n":n,"selection":family,"wall_seconds":time.perf_counter()-start,
                                "exact_energy_identity_defect":defect})
                print(p,n,family,r["free_vector_dofs"],r["relative_energy_norm"],flush=True)
                write(OUT/"seam-selection.json",results)
    if args.compare:
        expected = json.loads(args.compare.read_text())
        assert len(results)==len(expected)
        for r,e in zip(results,expected,strict=True):
            for k in ("p","n","selection","free_vector_dofs","c1_constraint_rank"):
                assert r[k]==e[k],(k,r,e)
            assert math.isclose(r["relative_energy_norm"],e["relative_energy_norm"],rel_tol=1e-6,abs_tol=2e-10)


if __name__=="__main__":
    main()
