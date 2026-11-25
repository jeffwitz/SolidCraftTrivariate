"""Recover a completed first CalculiX pass under a separate cumulative reference budget."""
import argparse
import os
from pathlib import Path
import sys
import time
from bracket_cost_fields import read,save,ROOT
from finer_matched_fields import measured_assessment
from iterative_tetra_fields import defect_corrected
from refine_tetra_fields import HERE,measure


def run(root,rx,ry,prefix):
    case=f"x{rx}-y{ry}";prior=ROOT/"build/j4j"/case
    old=read(prior/"preparation.json");passes=old["solve"]["passes"]
    if old["qualified"] or len(passes)!=2 or [p["status"] for p in passes]!=["completed","timeout"]:
        raise ValueError("reference is ineligible for first-pass reuse")
    case_dir=root/case;case_dir.mkdir(parents=True,exist_ok=True)
    budget=3600-old["mesh_and_solve_s"]
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="2",LD_LIBRARY_PATH=str(prefix/"lib"))
    args=[sys.executable,str(Path(__file__).resolve()),"worker",str(root),str(prefix),str(rx),str(ry)]
    process=measure(args,env,case_dir/"continuation.log",seconds=budget)
    prep=dict(radii_ticks=[rx,ry],size_mm=.9,prior_attempt=old,cumulative_limit_s=3600,
              continued_process=process,cumulative_wall_s=old["mesh_and_solve_s"]+process["elapsed_s"],qualified=False)
    if process["status"]=="completed":
        result=read(case_dir/"solve.json");prep["solve"]=result
        if result["status"]=="completed":
            prep["probe_run"]=measure([sys.executable,str(HERE/"fillet_probe.py"),str(case_dir/"fine"),
                str(case_dir/"probe"),"--ccx",str(prefix/"bin/ccx")],dict(os.environ,OPENBLAS_NUM_THREADS="1"),case_dir/"probe.log")
            if prep["probe_run"]["status"]=="completed":
                prep["probe"]=read(case_dir/"probe/result.json")
                prep["qualified"]=(prep["cumulative_wall_s"]<=3600 and result["free_residual_relative"][-1]<=1e-6
                    and result["result"]["reaction_relative_error"]<=1e-5 and prep["probe"]["distance_mm"]<.001)
    save(case_dir/"preparation.json",prep)
    measured_assessment(root,rx,ry,prefix)


def worker(root,rx,ry,prefix):
    case=f"x{rx}-y{ry}";prior=ROOT/"build/j4j"/case;dest=root/case
    old=read(prior/"preparation.json")
    result=defect_corrected(prior/"mesh",dest/"fine",prefix/"bin/ccx",dict(os.environ),
            3600-old["mesh_and_solve_s"],5,(prior/"fine/pass-1",old["solve"]["passes"][0]))
    save(dest/"solve.json",result)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=["run","worker"]);p.add_argument("root",type=Path)
    p.add_argument("prefix",type=Path);p.add_argument("rx",type=int);p.add_argument("ry",type=int)
    a=p.parse_args();(run if a.mode=="run" else worker)(a.root.resolve(),a.rx,a.ry,a.prefix.resolve())
