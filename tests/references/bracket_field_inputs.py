"""Rebuild J4e comparison fields under the EF4 memory and wall budgets."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
from iterative_tetra_fields import defect_corrected, displacements, physics
from refine_tetra_fields import HERE, measure

ROOT=HERE.parents[1]


def save(path, data):
    tmp=path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    tmp.replace(path)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(work, prefix):
    work.mkdir(exist_ok=True,parents=True)
    corpus=json.loads((HERE/"bracket_cost_corpus.json").read_text())["mandatory_radii_sixteenths_mm"]
    ccx=prefix/"bin/ccx"
    base=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="2")
    env=dict(base,LD_LIBRARY_PATH=str(prefix/"lib"),PYTHONPATH=str(prefix/"lib"))
    cases=[]
    for rx,ry in corpus:
        folder=work/f"x{rx}-y{ry}";folder.mkdir(exist_ok=True)
        baseline=ROOT/f"build/bracket-cost/x{rx}-y{ry}/tetra/h1.25-qNone"
        # Exact decks, displacements and solver physics are retained outside the timing study.
        row={"radii_ticks":[rx,ry],"reference_h125_sha256":digest(baseline/"bracket.inp"),
             "reference_h125_output_sha256":digest(baseline/"bracket.dat")}
        if rx==corpus[0][0] and ry==corpus[0][1]:
            cal=defect_corrected(baseline,folder/"calibration",ccx,env,1800,5)
            row["calibration"]=cal
            if cal["status"]=="completed":
                direct_u=displacements(baseline);other=displacements(folder/"calibration")
                if not np.array_equal(direct_u[:,0],other[:,0]):raise ValueError("C3D10 node order changed")
                row["calibrated_displacement_error"]=float(np.linalg.norm(other[:,1:]-direct_u[:,1:])/np.linalg.norm(direct_u[:,1:]))
                direct=physics(baseline)
                row["calibrated_energy_error"]=abs(cal["result"]["strain_energy_n_mm"]/direct["strain_energy_n_mm"]-1)
            if cal["status"]!="completed" or max(row["calibrated_displacement_error"],row["calibrated_energy_error"],cal["result"]["reaction_relative_error"])>1e-5:
                cases.append(row);save(work/"preparation.json",cases)
                raise ValueError("Fine C3D10 solver calibration failed")
        coarse=folder/"immersed.npz"
        row["immersed"]=measure([sys.executable,str(HERE/"immersed_bracket.py"),"--solver","multigrid",
            "--h","1.25","--order","8","--radius-x",str(rx/16),"--radius-y",str(ry/16),"--field-output",str(coarse)],
            base,folder/"immersed.json")
        if row["immersed"]["status"]=="completed":row["immersed_sha256"]=digest(coarse)
        fine=folder/"h1";mesh=folder/"deck";start=time.monotonic()
        row["mesh"]=measure([str(prefix/"bin/python"),str(HERE/"calculix_bracket.py"),"--ccx",str(ccx),
            "--deck-only","--store-displacements","--radius-x",str(rx/16),"--radius-y",str(ry/16),"--size","1",str(mesh)],
            env,folder/"mesh.log")
        if row["mesh"]["status"]=="completed":
            row["fine"]=defect_corrected(mesh,fine,ccx,env,max(0,1800-(time.monotonic()-start)),5)
            if row["fine"]["status"]=="completed":
                row["fine_qualified"]=row["fine"]["result"]["reaction_relative_error"]<=1e-5
        cases.append(row);save(work/"preparation.json",cases)
        print(rx,ry,row.get("fine_qualified"),row["immersed"]["status"],flush=True)
    return cases


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path);p.add_argument("freecad",type=Path)
    args=p.parse_args();run(args.work.resolve(),args.freecad.resolve())
