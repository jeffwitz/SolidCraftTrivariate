"""Monitor adaptive-bracket rounds and preserve accepted checkpoints at a cutoff."""
import argparse
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time

from refine_tetra_fields import session_rss


def run(binary,work,ticks,rounds,warm=False):
    work.mkdir(parents=True,exist_ok=True)
    log=work/"run.log"
    started=time.monotonic()
    records=[]
    last_round=None
    result="completed"
    with log.open("w") as output:
        task=subprocess.Popen([str(binary),str(ticks),str(work),str(rounds),*(["--warm"] if warm else [])],stdout=output,
                              stderr=subprocess.STDOUT,start_new_session=True)
        while task.poll() is None:
            matches=re.findall(r"adaptive round (\d+):",log.read_text())
            current=int(matches[-1]) if matches else 0
            now=time.monotonic()
            if current!=last_round:
                last_round=current;at=now
                records.append({"round":current,"peak_rss_kib":0,"elapsed_s":0.})
            row=records[-1]
            row["elapsed_s"]=now-at
            row["peak_rss_kib"]=max(row["peak_rss_kib"],session_rss(task.pid))
            if row["peak_rss_kib"]>=6*1024**2 or row["elapsed_s"]>=1800:
                result="rss_limit" if row["peak_rss_kib"]>=6*1024**2 else "time_limit"
                os.killpg(task.pid,signal.SIGKILL);break
            time.sleep(.2)
        rc=task.wait()
    if rc!=0 and result=="completed":result="execution_error"
    return {"radius_x_ticks":ticks,"process_status":result,"exit_code":rc,
            "elapsed_s":time.monotonic()-started,"round_resources":records,"log_tail":log.read_text().splitlines()[-2:]}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("binary",type=Path);p.add_argument("directory",type=Path)
    p.add_argument("--rounds",type=int,default=6)
    p.add_argument("--warm",action="store_true")
    args=p.parse_args()
    for ticks in (64,56):
        work=args.directory/f"x{ticks}"
        metrics=run(args.binary.resolve(),work,ticks,args.rounds,args.warm)
        (work/"resources.json").write_text(json.dumps(metrics,indent=2)+"\n")
        print(metrics,flush=True)
