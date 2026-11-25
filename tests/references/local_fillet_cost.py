"""Run and assess the preregistered J4d local fillet cost corpus."""
import json
import math
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
BASELINE = ROOT / "docs/results/bracket-cost-20261125.json"
MANIFEST = ROOT / "tests/references/bracket_cost_corpus.json"
OUTPUT = ROOT / "docs/results/local-fillet-cost-20261125.json"
WORK = ROOT / "build/bracket-cost/local-fillet"


def rel(a,b):
    return abs(a/b-1.0)


def stress_values(state):
    grouped={}
    for ident,value in state["fillet_limits_mpa"].items():
        root=ident.split("/")[0]
        grouped.setdefault(root,[]).append(float(value))
    return grouped


def pair_change(a,b):
    if not a or not b:
        return math.inf
    if a.keys()!=b.keys():
        return math.inf
    return max(rel(x,y) for key in a for x in a[key] for y in b[key])


def assess(measured, case, reference, manifest):
    expected_volume = 34080-math.pi*(360+8*sum((r/16)**2 for r in case))
    limits=manifest["qoi"]
    for path in measured["paths"]:
        prior=None
        first=None
        for state in path["states"]:
            if state.get("status")!="completed":
                state["checks"]={"completed":False}
                continue
            q6,q8=state["q6_q8"]
            stress8=stress_values(q8)
            q_stress=pair_change(stress_values(q6),stress8)
            q_tip=rel(q6["tip_mm"][1],q8["tip_mm"][1])
            q_volume=rel(q6["volume_mm3"],q8["volume_mm3"])
            q_area=rel(q6["loaded_area_mm2"],q8["loaded_area_mm2"])
            q_error=max(q_stress,q_tip,q_volume,q_area)
            reaction=math.dist(q8["reaction_n"],[0.,1000.,0.])/1000.
            checks={
                "tip_relative_to_finest_tetra":rel(q8["tip_mm"][1],reference["point_tip_y_mm"]),
                "fillet_relative_to_finest_tetra":max(rel(v,reference["von_mises_mpa"]) for values in stress8.values() for v in values),
                "volume_relative":rel(q8["volume_mm3"],expected_volume),
                "reaction_relative":reaction,
                "quadrature_relative":q_error,
                "true_residual":q8["true_projected_residual"],
                "mesh_stability":False,
            }
            if prior is not None:
                checks["tip_change_from_prior"] = rel(q8["tip_mm"][1],prior[0])
                checks["fillet_change_from_prior"] = pair_change(stress8,prior[1])
                checks["mesh_stability"]=(checks["tip_change_from_prior"]<=limits["successive_tip_relative_change"]
                    and checks["fillet_change_from_prior"]<=limits["successive_fillet_relative_change"])
            qualified=(reference["status"]=="stable" and checks["mesh_stability"]
                and checks["tip_relative_to_finest_tetra"]<=limits["tip_relative_to_finest_tetra"]
                and checks["fillet_relative_to_finest_tetra"]<=limits["fillet_von_mises_relative_to_finest_tetra"]
                and checks["volume_relative"]<=limits["volume_relative"]
                and checks["reaction_relative"]<=limits["reaction_relative"]
                and checks["quadrature_relative"]<=limits["quadrature_qoi_relative"]
                and checks["true_residual"]<=1e-10)
            state["checks"]=checks
            state["qualified"]=bool(qualified)
            if qualified and first is None:
                first={"step":state["step"],"leaves":state["leaves"],
                       "q8_wall_s":sum(q8[k] for k in ("assembly_s","force_s","constraint_s","solve_s","recovery_s")),
                       "search_wall_s": measured["geometry_build_and_audit_s"]
                           + state["partition_cumulative_s"]
                           + sum(
                               s["constraints_s"]
                               + sum(sum(q[k] for k in ("assembly_s","force_s","constraint_s","solve_s","recovery_s")) for q in s["q6_q8"])
                               for s in path["states"][:state["step"]+1]
                           )}
            prior=(q8["tip_mm"][1],stress8)
        path["first_qualified_state"]=first
    return measured


def main():
    WORK.mkdir(parents=True,exist_ok=True)
    manifest=json.loads(MANIFEST.read_text())
    baseline=json.loads(BASELINE.read_text())
    baseline_cases={tuple(c["radii_ticks"]):c for c in baseline["cases"]}
    executable=ROOT/"target/release/examples/local_fillet_cost"
    build=subprocess.run(["rtk","cargo","build","--release","--locked","--offline","--example","local_fillet_cost"],cwd=ROOT,capture_output=True,text=True,check=True)
    rows=[]
    env=dict(os.environ,OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="1")
    for path in ("from_h1","from_h2"):
      for radii in manifest["mandatory_radii_sixteenths_mm"]:
        case_dir=WORK/("x%d-y%d"%tuple(radii));case_dir.mkdir(parents=True,exist_ok=True)
        raw=case_dir/f"{path}.json";resources=case_dir/f"{path}-resources.txt"
        start=time.monotonic()
        cmd=["/usr/bin/time","-f","%e %M","-o",str(resources),str(executable),"--path",str(radii[0]),str(radii[1]),path,str(raw)]
        try:
            run=subprocess.run(cmd,cwd=ROOT,env=env,capture_output=True,text=True,timeout=1800)
            status="completed" if run.returncode==0 else "failed";diagnostic=(run.stderr+run.stdout)[-2500:]
        except subprocess.TimeoutExpired as exc:
            status="timeout_1800s";diagnostic=str(exc)[-2500:]
        wall=time.monotonic()-start
        row={"radii_ticks":radii,"path":path,"command":cmd,"wall_s":wall,"status":status,"diagnostic":diagnostic}
        if resources.exists():
            vals=resources.read_text().split()
            if len(vals)==2:row["process_wall_s"]=float(vals[0]);row["peak_rss_kib"]=int(vals[1])
        if raw.exists():
            measured=json.loads(raw.read_text());reference=baseline_cases[tuple(radii)]
            reference_run=next(r for r in reference["runs"] if r["method"]=="tetra" and r["resolution"]==1.25 and r["status"]=="completed")
            ref={"status":"stable","point_tip_y_mm":reference_run["result"]["point_tip_y_mm"],"von_mises_mpa":reference_run["probe"]["von_mises_mpa"]} if reference["reference_status"]=="stable" else {"status":"indeterminate"}
            row["reference"]=ref;row["measurements"]=assess(measured,radii,ref,manifest) if ref["status"]=="stable" else measured
        rows.append(row)
        temp=OUTPUT.with_suffix(".tmp")
        report={"schema":"j4d_local_fillet_cost_v1","manifest":manifest,"baseline_report":str(BASELINE.relative_to(ROOT)),
            "baseline_summary_revision":baseline["summary_revision"],"preserved_initial_trial_count":sum(len(c["runs"]) for c in baseline["cases"]),
            "commands":["rtk cargo build --release --locked --offline --example local_fillet_cost","rtk python3 tests/references/local_fillet_cost.py"],"radius_runs":rows}
        temp.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n");temp.replace(OUTPUT)
        print(radii,path,status,round(wall,2),row.get("peak_rss_kib"),flush=True)


def reassess_existing():
    manifest=json.loads(MANIFEST.read_text())
    baseline=json.loads(BASELINE.read_text())
    existing=json.loads(OUTPUT.read_text()) if OUTPUT.exists() else {"radius_runs":[]}
    baseline_cases={tuple(c["radii_ticks"]):c for c in baseline["cases"]}
    prior={(tuple(r["radii_ticks"]),r["path"]):r for r in existing["radius_runs"]}
    rows=[]
    for path in ("from_h1","from_h2"):
      for radii_list in manifest["mandatory_radii_sixteenths_mm"]:
        radii=tuple(radii_list);folder=WORK/("x%d-y%d"%radii);raw=folder/f"{path}.json"
        resource=folder/f"{path}-resources.txt";old=prior.get((radii,path));row=dict(old or {})
        row.update(radii_ticks=list(radii),path=path)
        if old is None:
            row.update(status="not_started_after_prior_h2_timeout" if path=="from_h2" else "missing_raw_output",
                       diagnostic="The serial runner stopped progressing through this path while the first h2 run remained timed out.")
        if raw.exists():
            measured=json.loads(raw.read_text());reference=baseline_cases[radii]
            ref_run=next(r for r in reference["runs"] if r["method"]=="tetra" and r["resolution"]==1.25 and r["status"]=="completed")
            ref={"status":"stable","point_tip_y_mm":ref_run["result"]["point_tip_y_mm"],"von_mises_mpa":ref_run["probe"]["von_mises_mpa"]} if reference["reference_status"]=="stable" else {"status":"indeterminate"}
            row["reference"]=ref
            row["measurements"]=assess(measured,radii,ref,manifest) if ref["status"]=="stable" else measured
            row["raw_output"]=str(raw.relative_to(ROOT))
            row["status"]="completed" if old is None or old.get("status")=="completed" else old["status"]
        if resource.exists():
            values=resource.read_text().split()
            if len(values)==2:
                row["process_wall_s"]=float(values[0]);row["peak_rss_kib"]=int(values[1]);row["resource_record"]="measured"
            else:
                row["resource_record"]="empty_time_output"
        else:
            row["resource_record"]="missing"
        rows.append(row)
    running=[]
    ticks=os.sysconf("SC_CLK_TCK")
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():continue
        try:
            cmd=(entry/"cmdline").read_bytes().replace(b"\0",b" ").decode(errors="replace").strip()
            stat=(entry/"stat").read_text().rsplit(")",1)[1].split();state=stat[0];start_ticks=int(stat[19])
            if "local_fillet_cost.py" in cmd or "local_fillet_cost --path" in cmd:
                uptime=float(Path("/proc/uptime").read_text().split()[0])
                running.append({"pid":int(entry.name),"state":state,"elapsed_s":round(uptime-start_ticks/ticks,2),"command":cmd[:500]})
        except (FileNotFoundError,PermissionError,ProcessLookupError,ValueError):pass
    report={"schema":"j4d_local_fillet_cost_v1","manifest":manifest,"baseline_report":str(BASELINE.relative_to(ROOT)),
        "baseline_summary_revision":baseline["summary_revision"],"preserved_initial_trial_count":sum(len(c["runs"]) for c in baseline["cases"]),
        "commands":["rtk cargo build --release --locked --offline --example local_fillet_cost","rtk python3 tests/references/local_fillet_cost.py"],
        "restoration_source":"raw JSON outputs and per-path /usr/bin/time files; no solves rerun",
        "background_attempts_at_restore":running,"radius_runs":rows}
    temp=OUTPUT.with_suffix(".tmp");temp.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n");temp.replace(OUTPUT)


if __name__=="__main__":
    if os.environ.get("J4D_REASSESS_ONLY")=="1":
        reassess_existing()
    else:
        main()
