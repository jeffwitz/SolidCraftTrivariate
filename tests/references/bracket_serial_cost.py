"""Sequential, CPU-pinned J4r cost comparison; never reuses trial results."""
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
PROTO=ROOT/"tests/references/bracket_serial_protocol.json"
MANIFEST=ROOT/"tests/references/bracket_cost_corpus.json"
BASELINE=ROOT/"docs/results/bracket-cost-20261125.json"
LOCAL=ROOT/"docs/results/local-fillet-cost-20261125.json"
OUTPUT=ROOT/"docs/results/bracket-serial-cost-20261125.json"
FREECAD=Path("/home/jeff/Code/FreeCAD/.pixi/envs/default")

def load(path):return json.loads(path.read_text())
def save(data):
    tmp=OUTPUT.with_suffix(".tmp");tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n");tmp.replace(OUTPUT)
def loadavg():
    try:return [float(x) for x in Path("/proc/loadavg").read_text().split()[:3]]
    except OSError:return None

def active_campaign_processes():
    needles=("local_fillet_cost.py","local_fillet_cost --path","bracket_field_inputs.py",
             "full_cylinder_correction.py","equivalent_prism_mechanics.py","immersed_bracket.py")
    found=[]
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():continue
        try:
            cmd=(entry/"cmdline").read_bytes().replace(b"\0",b" ").decode(errors="replace").strip()
            if cmd and any(n in cmd for n in needles):found.append({"pid":int(entry.name),"command":cmd[:360]})
        except (FileNotFoundError,PermissionError,ProcessLookupError):pass
    return found

def first_local_tiers(local_report):
    candidates_by_radius={}
    for row in local_report["radius_runs"]:
        if row.get("status")!="completed" or "measurements" not in row:continue
        for path in row["measurements"].get("paths",[]):
            tier=path.get("first_qualified_state")
            if tier is not None:
                q8=next(s["q6_q8"][1] for s in path["states"] if s["step"]==tier["step"])
                candidates_by_radius.setdefault(tuple(row["radii_ticks"]),[]).append(
                    (tier.get("search_wall_s",math.inf),path["path"],tier,q8))
    return {radii:dict(path=path,step=tier["step"],leaves=tier["leaves"],search_wall_s=wall,j4d_q8=q8)
            for radii,choices in candidates_by_radius.items()
            for wall,path,tier,q8 in [min(choices,key=lambda x:x[0])]}

def metrics(method,result):
    if method=="tetra":
        r=result["result"];return {"tip_mm":r["point_tip_y_mm"],"stress_mpa":result["probe"]["von_mises_mpa"],
            "volume_mm3":r["occ_volume_mm3"],"reaction_n":r["reaction_n"],"residual":None}
    if method=="immersed":
        r=result["result"];return {"tip_mm":r["point_tip_y_mm"],"stress_mpa":r["von_mises_mpa"],
            "volume_mm3":r["volume_mm3"],"reaction_n":r["reaction_n"],"residual":r.get("relative_residual")}
    r=result["result"]["q8"]
    stress=max(abs(float(x)) for x in r["fillet_limits_mpa"].values())
    return {"tip_mm":r["tip_mm"][1],"stress_mpa":stress,"volume_mm3":r["volume_mm3"],
        "reaction_n":r["reaction_n"],"residual":r["true_projected_residual"]}

def physical_qualification(method,result,radii,baseline_case,manifest,local_tier=None,trial_dir=None):
    if result.get("status")!="completed" or "result" not in result:
        return {"status":"not_evaluated","reason":"trial did not complete"}
    if baseline_case.get("reference_status")!="stable":
        return {"status":"indeterminate","reason":"frozen finest tetra reference is not stable"}
    m=metrics(method,result)
    ref=next(r for r in baseline_case["runs"] if r["method"]=="tetra" and r["resolution"]==1.25 and r["status"]=="completed")
    expected_volume=34080-math.pi*(360+8*sum((r/16)**2 for r in radii))
    q=manifest["qoi"]
    checks={"tip_relative_to_finest_tetra":abs(m["tip_mm"]/ref["result"]["point_tip_y_mm"]-1),
        "fillet_relative_to_finest_tetra":abs(m["stress_mpa"]/ref["probe"]["von_mises_mpa"]-1),
        "volume_relative":abs(m["volume_mm3"]/expected_volume-1),
        "reaction_relative":math.dist(m["reaction_n"],[0.,1000.,0.])/1000.,
        "tree_rss_within_cap":result.get("tree_peak_rss_kib",0)<=manifest["cost"]["max_tree_rss_kib"]}
    checks["qoi_qualified"]=(checks["tip_relative_to_finest_tetra"]<=q["tip_relative_to_finest_tetra"]
        and checks["fillet_relative_to_finest_tetra"]<=q["fillet_von_mises_relative_to_finest_tetra"]
        and checks["volume_relative"]<=q["volume_relative"] and checks["reaction_relative"]<=q["reaction_relative"]
        and checks["tree_rss_within_cap"]
        and (m["residual"] is None or m["residual"]<=1e-10))
    if method=="tetra":
        deck=(trial_dir/"mesh"/"bracket.inp").read_text()
        checks["clamp_boundary_prescribed"]="*BOUNDARY\nCLAMP,1,3" in deck
        drift_ref=next(r for r in baseline_case["runs"] if r["method"]=="tetra" and r["resolution"]==1.5 and r["status"]=="completed")
        drift=max(abs(m["tip_mm"]/drift_ref["result"]["point_tip_y_mm"]-1),
            abs(m["stress_mpa"]/drift_ref["probe"]["von_mises_mpa"]-1),
            abs(m["volume_mm3"]/drift_ref["result"]["occ_volume_mm3"]-1),
            math.dist(m["reaction_n"],drift_ref["result"]["reaction_n"])/1000.)
    elif method=="immersed":
        drift_ref=next(r for r in baseline_case["runs"] if r["method"]=="immersed" and r["resolution"]==1.25 and r["q"]==8 and r["status"]=="completed")
        d=drift_ref["result"]
        drift=max(abs(m["tip_mm"]/d["point_tip_y_mm"]-1),abs(m["stress_mpa"]/d["von_mises_mpa"]-1),
            abs(m["volume_mm3"]/d["volume_mm3"]-1),math.dist(m["reaction_n"],d["reaction_n"])/1000.)
    else:
        d=local_tier["j4d_q8"]
        drift=max(abs(m["tip_mm"]/d["tip_mm"][1]-1),
            math.dist(m["reaction_n"],d["reaction_n"])/1000.,
            abs(m["volume_mm3"]/d["volume_mm3"]-1),abs(m["stress_mpa"]/max(abs(float(x)) for x in d["fillet_limits_mpa"].values())-1))
    checks["drift_from_frozen_replay_relative"]=drift
    checks["drift_qualified"]=drift<=1e-8
    checks["status"]="qualified" if checks["qoi_qualified"] and checks["drift_qualified"] and checks.get("clamp_boundary_prescribed",True) else "failed"
    return checks

def timed(cmd,env,folder,timeout):
    wallfile=folder/"time.txt";log=folder/"stdout.log";err=folder/"stderr.log"
    full=["/usr/bin/taskset","-c","2,3","/usr/bin/time","-f","%e %M","-o",str(wallfile),*map(str,cmd)]
    start=time.monotonic();load_before=loadavg();interference_before=active_campaign_processes()
    try:
        proc=subprocess.run(full,cwd=ROOT,env=env,capture_output=True,text=True,timeout=timeout)
        log.write_text(proc.stdout);err.write_text(proc.stderr)
        status="completed" if proc.returncode==0 else "failed"
        diagnostic=(proc.stderr+proc.stdout)[-2000:]
    except subprocess.TimeoutExpired as e:
        status="timeout_1800s";diagnostic=str(e)[-2000:]
        log.write_text((e.stdout or "") if isinstance(e.stdout,str) else "")
        err.write_text((e.stderr or "") if isinstance(e.stderr,str) else "")
    wall=time.monotonic()-start
    usage={}
    if wallfile.exists():
        vals=wallfile.read_text().split()
        if len(vals)==2:usage={"process_wall_s":float(vals[0]),"tree_peak_rss_kib":int(vals[1])}
    return dict(status=status,wall_s=wall,command=full,diagnostic=diagnostic,
                interference={"loadavg_before":load_before,"loadavg_after":loadavg(),
                              "campaign_processes_before":interference_before,
                              "campaign_processes_after":active_campaign_processes()},**usage)

def tetra_trial(radii,folder,env,remaining,runner=timed):
    folder.mkdir(parents=True,exist_ok=False)
    ccx=FREECAD/"bin/ccx";py=FREECAD/"bin/python"
    env=dict(env,LD_LIBRARY_PATH=str(FREECAD/"lib"),PYTHONPATH=str(FREECAD/"lib"))
    hx=folder/"mesh"
    args=[py,ROOT/"tests/references/calculix_bracket.py","--ccx",ccx,"--size","1.5","--store-displacements",
          "--radius-x",str(radii[0]/16),"--radius-y",str(radii[1]/16),hx]
    run=runner(args,env,folder,remaining)
    if run["status"]!="completed":return run
    remaining=max(0.,remaining-run["wall_s"])
    if remaining<=0:
        run["status"]="timeout_1800s";run["physical_audit_s"]=0.
        return run
    probe=folder/"probe";probe.mkdir()
    probe_env=dict(env,OPENBLAS_NUM_THREADS="1")
    extra=runner([sys.executable,ROOT/"tests/references/fillet_probe.py",hx,probe,"--ccx",ccx],
                probe_env,folder,remaining)
    run["physical_audit_s"]=extra["wall_s"];run["status"] = extra["status"]
    run["diagnostic"] += extra["diagnostic"]
    run["tree_peak_rss_kib"]=max(run.get("tree_peak_rss_kib",0),extra.get("tree_peak_rss_kib",0))
    if run["status"]=="completed":
        run["result"]=load(hx/"result.json");run["probe"]=load(probe/"result.json")
    return run

def immersed_trial(radii,folder,env,remaining):
    folder.mkdir(parents=True,exist_ok=False)
    result=timed([sys.executable,ROOT/"tests/references/immersed_bracket.py","--solver","multigrid","--h","1.25","--order","8",
                  "--radius-x",str(radii[0]/16),"--radius-y",str(radii[1]/16)],env,folder,remaining)
    if result["status"]=="completed":result["result"]=json.loads((folder/"stdout.log").read_text())
    return result

def local_trial(radii,tier,folder,env,remaining):
    folder.mkdir(parents=True,exist_ok=False)
    binary=ROOT/"target/release/examples/local_fillet_cost"
    result_path=folder/"physical.json"
    result=timed([binary,"--tier",radii[0],radii[1],tier["path"],tier["step"],result_path],env,folder,remaining)
    if result["status"]=="completed":result["result"]=load(result_path)
    return result

def main():
    protocol=load(PROTO);manifest=load(MANIFEST);baseline=load(BASELINE)
    if not LOCAL.exists():
        save({"schema":protocol["schema"],"status":"waiting_for_j4d","reason":"local J4d qualification report is not present; no substituted tier will be used"})
        raise SystemExit("J4d report absent; no serial timing run started")
    local_report=load(LOCAL);tiers=first_local_tiers(local_report)
    baseline_cases={tuple(c["radii_ticks"]):c for c in baseline["cases"]}
    env=dict(os.environ,OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="2",NUMEXPR_NUM_THREADS="1")
    binary=ROOT/"target/release/examples/local_fillet_cost"
    subprocess.run(["rtk","cargo","build","--release","--locked","--offline","--example","local_fillet_cost"],cwd=ROOT,check=True)
    original=sum(len(c["runs"]) for c in baseline["cases"])
    rows=[]
    for order in protocol["balanced_order_by_radius_and_repetition"]:
        radii=tuple(order["radii"]);case_dir=ROOT/"build/bracket-cost/serial"/("x%d-y%d-rep%d"%(radii[0],radii[1],order["repeat"]))
        case_dir.mkdir(parents=True,exist_ok=True)
        before=loadavg();processes_before=active_campaign_processes();start=time.monotonic();trials=[]
        for method in order["methods"]:
            if method=="local" and radii not in tiers:
                trials.append({"method":"local","status":"not_admissible_no_qualified_j4d_tier","substitution":False});continue
            trial_dir=case_dir/method
            if method=="tetra":result=tetra_trial(radii,trial_dir,env,1800)
            elif method=="immersed":result=immersed_trial(radii,trial_dir,env,1800)
            else:
                result=local_trial(radii,tiers[radii],trial_dir,env,1800)
                result["j4d_tier"]={k:v for k,v in tiers[radii].items() if k!="j4d_q8"}
            if result.get("status")=="completed":
                result["physical_qualification"]=physical_qualification(method,result,radii,baseline_cases[radii],manifest,
                    tiers.get(radii) if method=="local" else None,trial_dir)
            trials.append({"method":method,**result})
            save({"schema":protocol["schema"],"status":"running","completed_schedule_rows":rows,
                  "active_row":{"radii":radii,"repeat":order["repeat"],"trials":trials},
                  "interference":{"loadavg_before":before,"loadavg_current":loadavg(),"cpu_affinity":[2,3]}})
        rows.append({"radii_ticks":list(radii),"repeat":order["repeat"],"planned_method_order":order["methods"],
                     "elapsed_wall_s":time.monotonic()-start,"loadavg_before":before,"loadavg_after":loadavg(),
                     "campaign_processes_before":processes_before,"campaign_processes_after":active_campaign_processes(),"trials":trials})
        save({"schema":protocol["schema"],"status":"running","completed_schedule_rows":rows})
        print(radii,order["repeat"],[(x["method"],x["status"],round(x.get("wall_s",0),1)) for x in trials],flush=True)
    final={"schema":protocol["schema"],"status":"completed","protocol":protocol,
           "baseline_reference_revision":baseline["summary_revision"],"initial_j4b_trials_preserved":original,
           "python":platform.python_version(),"cpu_affinity":[2,3],"thread_environment":{k:env[k] for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS")},
           "local_tiers":{str(k):{key:value for key,value in v.items() if key!="j4d_q8"} for k,v in tiers.items()},"schedule_rows":rows}
    save(final)

if __name__=="__main__":main()
