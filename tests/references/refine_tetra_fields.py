"""Run finer OCC/Netgen/CalculiX references under a shared process-tree budget."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
PAGE_KIB=os.sysconf("SC_PAGE_SIZE")//1024


def session_rss(sid):
    total=0
    for p in Path("/proc").iterdir():
        if not p.name.isdigit():continue
        try:
            # After comm, stat fields are state, parent, pgrp, session, ...
            fields=(p/"stat").read_text().rsplit(")",1)[1].split()
            if int(fields[3])==sid: total+=int(fields[21])*PAGE_KIB
        except (FileNotFoundError,PermissionError,ProcessLookupError):
            pass
    return total


def measure(cmd,env,log):
    start=time.monotonic();peak=0;status="completed"
    with log.open("w") as stream:
        child=subprocess.Popen(cmd,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        while child.poll() is None:
            peak=max(peak,session_rss(child.pid))
            if time.monotonic()-start>=1800 or peak>=6*1024**2:
                status="timeout" if time.monotonic()-start>=1800 else "rss_budget_exceeded"
                os.killpg(child.pid,signal.SIGKILL)
                break
            time.sleep(.15)
        rc=child.wait()
    if rc and status=="completed":status="execution_error"
    return {"status":status,"exit_code":rc,"elapsed_s":time.monotonic()-start,"tree_peak_rss_kib":peak}


def run(root,freecad):
    root.mkdir(exist_ok=True,parents=True)
    env=dict(os.environ,LD_LIBRARY_PATH=str(freecad/"lib"),PYTHONPATH=str(freecad/"lib"),
             OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="2")
    ccx=str(freecad/"bin/ccx")
    rows=[]
    for ticks in (64,56):
        for size in (1.,.9):
            work=root/f"x{ticks}-h{size:g}"
            work.mkdir(exist_ok=True)
            options=[str(freecad/"bin/python"),str(HERE/"calculix_bracket.py"),"--ccx",ccx,"--store-displacements",
                     "--radius-x",str(ticks/16.),"--size",str(size),str(work)]
            record=measure(options,env,work/"geometry-and-solve.log")
            row={"radius_x_ticks":ticks,"size_mm":size,"geometry_and_solve":record}
            if record["status"]=="completed":
                row["result"]=json.loads((work/"result.json").read_text())
                post=measure([sys.executable,str(HERE/"fillet_probe.py"),str(work),str(work/"probe"),"--ccx",ccx],
                             dict(os.environ,OPENBLAS_NUM_THREADS="1"),work/"probe.log")
                row["probe_run"]=post
                if post["status"]=="completed":
                    row["probe"]=json.loads((work/"probe/result.json").read_text())
            rows.append(row)
            temp=root/"metrics.tmp"
            temp.write_text(json.dumps(rows,indent=2)+"\n");temp.replace(root/"metrics.json")
            print(ticks,size,record,flush=True)
            if record["status"] in ("rss_budget_exceeded","timeout"):break


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path);p.add_argument("--freecad-prefix",type=Path,required=True)
    a=p.parse_args();run(a.directory,a.freecad_prefix)
