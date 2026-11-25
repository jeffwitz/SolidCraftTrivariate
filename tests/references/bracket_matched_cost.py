"""Time three field-qualified chains and verify their full replay coefficients."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

from bracket_cost_fields import read, save
from bracket_field_distance import TetraField
from bracket_serial_cost import (ROOT, MANIFEST, BASELINE, PROTO, active_campaign_processes,
                                 loadavg, physical_qualification, tetra_trial)
from refine_tetra_fields import measure


def monitored(cmd, env, folder, seconds):
    args = ["/usr/bin/taskset", "-c", "2,3", *map(str,cmd)]
    log = folder/"monitored.log"
    before = loadavg()
    running = active_campaign_processes()
    r = measure(args, env, log, cwd=ROOT, seconds=seconds)
    return dict(**r, wall_s=r["elapsed_s"], command=args, loadavg=[before,loadavg()],
                campaign_processes_before=running, diagnostic=log.read_text()[-1000:])


def admissible(comparison, name):
    r = comparison["levels"][-1]
    d = r["distances"][name+"__tetra_h1"]
    return (comparison["all_points_located_once"] and comparison["recovered_energies"]
            and all(d[k]["stable"] and d[k]["mean"] <= lim for k,lim in [("energy",.03),("l2",.005)]))


def verify(method, case, step, dest, fine, tiers):
    if method=="local":
        def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
        a=digest(dest/"field.bin");b=digest(tiers/case/f"s{step}.bin")
        return dict(qualified=a==b, field_sha256=a, reference_sha256=b, geometry_in_checkpoint=True)
    if method=="immersed":
        with np.load(dest/"field.npz") as a,np.load(fine/case/"immersed.npz") as b:
            geom=all(np.array_equal(a[k],b[k]) for k in ("h_mm","radii_mm"))
            u,v=a["displacement_mm"],b["displacement_mm"]
    else:
        a=TetraField(dest/"mesh");b=TetraField(ROOT/"build/bracket-cost"/case/"tetra/h1.5-qNone")
        geom=np.array_equal(a.x,b.x)
        u,v=a.u,b.u
    delta=float(np.linalg.norm(u-v)/np.linalg.norm(v)) if u.shape==v.shape else None
    return dict(qualified=bool(geom and delta is not None and delta<=1e-8),coefficients=int(u.size),
                geometry_equal=geom,relative_displacement_coefficients=delta)


def run(work,output):
    result_dir=ROOT/"docs/results"
    serial=read(result_dir/"bracket-cost-fields-20261125.json")
    candidate=read(result_dir/"bracket-candidate-fields-20261125.json")
    if not(serial["finished"] and candidate["finished"]):raise ValueError("field qualification is incomplete")
    fields={tuple(r["radii_ticks"]):r for r in serial["cases"]}
    native={tuple(r["radii_ticks"]):r for r in candidate["cases"]}
    trials={tuple(c["radii_ticks"]):c for c in read(BASELINE)["cases"]}
    ladder=read(result_dir/"bracket-field-tiers-20261125.json")["states"]
    fine=ROOT/"build/j4e";tiers=ROOT/"build/j4f"
    exe=ROOT/"target/release/examples/local_fillet_cost"
    env=dict(os.environ,OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="2",NUMEXPR_NUM_THREADS="1")
    protocol=read(PROTO);manifest=read(MANIFEST)
    work.mkdir(parents=True,exist_ok=False)
    report=dict(schema="field_matched_cost_v1",base_revision="05f3b7b",schedule=protocol["balanced_order_by_radius_and_repetition"],
                cpu_affinity=[2,3],sampled_session_rss_interval_s=.15,rows=[],finished=False)
    for row in report["schedule"]:
        radii=tuple(row["radii"]);case="x%d-y%d"%radii
        rec=dict(radii_ticks=list(radii),repeat=row["repeat"],trials=[])
        report["rows"].append(rec)
        frow, nrow = fields[radii],native[radii]
        step=nrow.get("result",{}).get("earliest_qualified")
        for method in row["methods"]:
            dest=work/f"{case}-rep{row['repeat']}"/method
            dest.parent.mkdir(exist_ok=True)
            eligible=(step is not None if method=="local" else frow["status"]=="completed" and
                      admissible(frow["result"],"tetra_h15" if method=="tetra" else "immersed_h125"))
            if not eligible:
                rec["trials"].append(dict(method=method,status="missing_field_qualification"));continue
            tier=None
            if method=="local":
                state=next(s for s in ladder if s["radii_ticks"]==list(radii) and s["step"]==step)
                tier=dict(j4d_q8=state["q6_q8"][1],step=step)
                dest.mkdir()
                v=monitored([exe,"--field-tier",*radii,step,8,dest/"result.json",dest/"field.bin"],env,dest,1800)
                if v["status"]=="completed":v["result"]=read(dest/"result.json")
                v["selected_step"]=step
            elif method=="tetra":
                v=tetra_trial(radii,dest,env,1800,runner=monitored)
            else:
                dest.mkdir()
                cmd=[sys.executable,ROOT/"tests/references/immersed_bracket.py","--solver","multigrid",
                     "--h","1.25","--order","8","--radius-x",radii[0]/16,"--radius-y",radii[1]/16,
                     "--field-output",dest/"field.npz"]
                v=monitored(cmd,env,dest,1800)
                if v["status"]=="completed":v["result"]=read(dest/"monitored.log")
            v["end_to_end_with_probe_s"]=v.get("wall_s",0.)+v.get("physical_audit_s",0.)
            if v["status"]=="completed":
                v["observables"]=physical_qualification(method,v,radii,trials[radii],manifest,tier,dest)
                start=time.monotonic()
                v["field_replay"]=verify(method,case,step,dest,fine,tiers)
                v["coefficient_verification_s"]=time.monotonic()-start
                v["qualified"]=v["observables"]["status"]=="qualified" and v["field_replay"]["qualified"]
            v["method"]=method
            rec["trials"].append(v);save(output,report)
            print(case,row["repeat"],method,v["status"],v.get("qualified"),round(v["end_to_end_with_probe_s"],2),flush=True)
    report["finished"]=True
    save(output,report)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args()
    run(a.directory.resolve(),a.output)
