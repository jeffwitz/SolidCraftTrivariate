"""Run the independent bracket series with a per-run RSS and time budget."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
WALL_LIMIT=1800.
RSS_LIMIT=6*1024*1024


def run(directory,h,order,solver,radius_x=4.,radius_y=4.):
    out=directory/f"h{h}-q{order}.json"
    err=out.with_suffix(".log")
    started=time.monotonic()
    peak=0
    status="completed"
    with out.open("w") as stdout,err.open("w") as stderr:
        task=subprocess.Popen([sys.executable,str(HERE/"immersed_bracket.py"),"--h",str(h),"--order",str(order),"--solver",solver,"--radius-x",str(radius_x),"--radius-y",str(radius_y)],
                              stdout=stdout,stderr=stderr,env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1"))
        while task.poll() is None:
            try:
                stat=Path(f"/proc/{task.pid}/status").read_text()
                rss=next((int(s.split()[1]) for s in stat.splitlines() if s.startswith("VmHWM:")),0)
                peak=max(rss,peak)
            except (OSError, StopIteration):
                pass
            if time.monotonic()-started>WALL_LIMIT or peak>RSS_LIMIT:
                status="timeout" if time.monotonic()-started>WALL_LIMIT else "rss_budget_exceeded"
                task.kill()
                break
            time.sleep(.1)
        code=task.wait()
    if code and status=="completed":
        status="execution_error"
    info={"h_mm":h,"order":order,"wall_s":time.monotonic()-started,"max_rss_kib":peak,"status":status,"exit_code":code}
    if code==0:
        info["result"]=json.loads(out.read_text())
    else:
        info["diagnostic"]=err.read_text()[-1200:]
    return info


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path)
    p.add_argument("--h",type=float,nargs="+",default=[5.,2.5,1.25,1.])
    p.add_argument("--order",type=int,nargs="+",default=[4,6,8])
    p.add_argument("--solver",choices=("lu","multigrid"),default="lu")
    p.add_argument("--radius-x",type=float,default=4.)
    p.add_argument("--radius-y",type=float,default=4.)
    args=p.parse_args()
    args.directory.mkdir(parents=True,exist_ok=True)
    report={"solver":args.solver,"resource_budget":{"wall_s":WALL_LIMIT,"rss_kib":RSS_LIMIT},"runs":[]}
    target=args.directory/"result.json"
    for h in args.h:
        for q in args.order:
            case=run(args.directory,h,q,args.solver,args.radius_x,args.radius_y)
            report["runs"].append(case)
            temp=target.with_suffix(".tmp")
            temp.write_text(json.dumps(report,indent=2)+"\n")
            temp.replace(target)
            print(h,q,case["status"],case["wall_s"],flush=True)
            if case["status"] in ("timeout","rss_budget_exceeded"):
                # At a finer quadrature the same matrix graph and solver are required.
                # Record the budget failure; leave other resolutions independently testable.
                break
