"""Run the frozen Kirsch pilot with process-level elapsed time and peak memory."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "tests/references/j4_pilot_corpus.json"
REF = ROOT / "tests/references"


def trials(method, radius, env_dir):
    if method == "native":
        for level in range(1,6):
            for order in (4,8):
                yield f"l{level}-q{order}", [
                    str(ROOT/"target/release/examples/kirsch"),str(level),str(order),str(radius),"quarter"]
    elif method == "immersed":
        for n in (4,8,16,32,64):
            for order in (6,12):
                yield f"n{n}-q{order}", [sys.executable,str(REF/"immersed_kirsch.py"),str(n),str(order),str(radius)]
    else:
        for h in (2,1,0.7,0.5,0.4):
            yield f"h{h}-q6", [str(env_dir/"bin/python"),str(REF/"calculix_kirsch.py"),
                               "--ccx",str(env_dir/"bin/ccx"),"--radius",str(radius),"--size",str(h)]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("method",choices=("native","immersed","tetra"))
    p.add_argument("directory",type=Path)
    p.add_argument("--freecad-env",type=Path,default=Path("/home/jeff/Code/FreeCAD/.pixi/envs/default"))
    p.add_argument("--verify-tetra-quadrature",action="store_true",
                   help="repeat the finest tetrahedral mesh using surface quadrature 12")
    args=p.parse_args()
    if args.verify_tetra_quadrature and args.method!="tetra":
        p.error("quadrature verification applies to the tetrahedral reference")
    manifest=json.loads(MANIFEST.read_text())
    protocol=manifest["analysis"]
    env=dict(os.environ,OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="1")
    if args.method=="tetra":
        env.update(PYTHONPATH=str(args.freecad_env/"lib"),LD_LIBRARY_PATH=str(args.freecad_env/"lib"))
    for case in manifest["mandatory"]:
        radius=case["bore_radius"]
        for name,cmd in trials(args.method,radius,args.freecad_env):
            if args.verify_tetra_quadrature:
                if name != "h0.4-q6":
                    continue
                name="h0.4-q12"
                cmd.extend(["--quadrature","12"])
            dest=(args.directory/args.method/case["id"]/name).resolve()
            dest.mkdir(parents=True,exist_ok=True)
            if args.method=="tetra":
                cmd.append(str(dest))
            start=time.monotonic()
            with (dest/"stdout.log").open("w") as out, (dest/"stderr.log").open("w") as err:
                try:
                    result=subprocess.run(["/usr/bin/time","-f","%e %M","-o",str(dest/"resources"),*cmd],
                                          cwd=ROOT,env=env,stdout=out,stderr=err,check=False,
                                          timeout=protocol["wall_clock_budget_seconds_per_method_case"])
                    status="measured" if result.returncode==0 else "execution_error"
                except subprocess.TimeoutExpired:
                    status="timeout"
            report={"method":args.method,"case":case["id"],"configuration":name,"status":status,
                    "elapsed_s":time.monotonic()-start}
            if (dest/"resources").is_file():
                tokens=(dest/"resources").read_text().strip().split()
                if len(tokens)>=2 and tokens[-1].isdigit():
                    report["process_peak_rss_kib"]=int(tokens[-1])
            if status=="measured":
                text=(dest/"result.json" if args.method=="tetra" else dest/"stdout.log").read_text()
                report["metrics"]=json.loads(text)
                if report["process_peak_rss_kib"]>1024*protocol["memory_budget_mib_per_run"]:
                    report["status"]="resource_exceeded"
            (dest/"trial.json").write_text(json.dumps(report,indent=2)+"\n")
            print(case["id"],name,report["status"],round(report["elapsed_s"],2),flush=True)


if __name__=="__main__":
    main()
