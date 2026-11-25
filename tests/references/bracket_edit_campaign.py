"""Run and compare the three independent edited-bracket mechanics chains."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
Q = (6,8)
H_TET = (3,2,1.5,1.25)
H_BG = (5.,2.5,1.25,1.)
VOLUME = 34080. - math.pi*(360.+8.*(3.5**2+4.**2))


def read(path):
    return json.loads(path.read_text())


def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    tmp.replace(path)


def native_series(directory):
    directory.mkdir(parents=True,exist_ok=True)
    runs=[]
    for level in range(4):
        for q in Q:
            out=directory/f"h{level}-q{q}.json"
            start=time.monotonic()
            peak=0
            with out.open("w") as stdout, out.with_suffix(".log").open("w") as stderr:
                proc=subprocess.Popen([str(ROOT/"target/release/examples/bracket_convergence"),
                                       str(level),str(q),"56","64"],stdout=stdout,stderr=stderr)
                outcome="completed"
                while proc.poll() is None:
                    try:
                        lines=Path(f"/proc/{proc.pid}/status").read_text().splitlines()
                        peak=max(peak,next(int(s.split()[1]) for s in lines if s.startswith("VmHWM:")))
                    except (OSError,StopIteration):
                        pass
                    if time.monotonic()-start>1800 or peak>6*1024**2:
                        outcome="timeout" if time.monotonic()-start>1800 else "rss_budget_exceeded"
                        proc.kill()
                        break
                    time.sleep(.1)
                rc=proc.wait()
            if rc and outcome=="completed":
                outcome="execution_error"
            case=dict(level=level,q=q,exit_code=rc,status=outcome,wall_s=time.monotonic()-start,max_rss_kib=peak)
            if not rc:
                case["result"]=read(out)
            else:
                case["diagnostic"]=out.with_suffix(".log").read_text()[-1500:]
            runs.append(case)
            save(directory/"runs.json",runs)
            print(level,q,outcome,flush=True)


def summarize(directory):
    tet=[]
    for h in H_TET:
        path=directory/f"r14-tet/h{h}"
        result=read(path/"result.json")
        s=read(path/"probe/result.json")
        assert result["radius_x_mm"]==3.5 and result["radius_y_mm"]==4.
        assert abs(result["occ_volume_mm3"]-VOLUME)<1e-7
        rss_upper=result["solver_peak_rss_kib"] + result["mesh_process_peak_rss_kib"]
        tet.append(dict(size_mm=h,displacement=result,stress=s,combined_rss_upper_kib=rss_upper,
                        geometry_solve_s=result["geometry_mesh_deck_s"]+result["solve_s"]))
    ref=tet[-1]
    results=[]
    prev=None
    def assess(tip,stress,volume,force,reaction):
        return dict(tip_relative=abs(tip/ref["displacement"]["point_tip_y_mm"]-1),
                    stress_relative=max(abs(x/ref["stress"]["von_mises_mpa"]-1) for x in stress),
                    volume_relative=abs(volume/VOLUME-1),
                    force_relative=math.dist(force,[0,-1000,0])/1000,
                    balance_relative=math.sqrt(sum((a+b)**2 for a,b in zip(force,reaction)))/1000)
    def passed(a):
        return (a["tip_relative"]<=.01 and a["stress_relative"]<=.05 and a["volume_relative"]<=1e-8
                and a["force_relative"]<=1e-6 and a["balance_relative"]<=1e-6)
    for t in tet:
        r,s=t["displacement"],t["stress"]
        qoi=assess(r["point_tip_y_mm"],[s["von_mises_mpa"]],r["occ_volume_mm3"],[0,-1000.,0],r["reaction_n"])
        change=abs(s["von_mises_mpa"]/prev-1) if prev else None
        ok=passed(qoi) and change is not None and change<=.03 and t["combined_rss_upper_kib"]<=6*1024**2 and t["geometry_solve_s"]<=1800
        results.append(dict(method="tetrahedral",resolution=t["size_mm"],checks=qoi,stress_mesh_change=change,qualified=ok))
        prev=s["von_mises_mpa"]
    assert results[-1]["qualified"],"final tetra reference did not stabilize"
    immersed=read(directory/"r14-immersed/result.json")
    native=read(directory/"r14-native/runs.json")
    for name,cases in (("immersed",immersed["runs"]),("native",native)):
        prev={}
        by_resolution={}
        for complete in cases:
            if complete["status"]=="completed":
                row=complete["result"]
                by_resolution.setdefault(row.get("level",row.get("h_mm")), {})[row["gauss_order"]]=row
        for case in cases:
            if case["status"]!="completed":
                results.append(dict(method=name,case=case,qualified=False))
                continue
            r=case["result"]
            assert r["radius_x_mm"]==3.5 and r["radius_y_mm"]==4.
            stresses=(list(r["fillet_limits_mpa"].values()) if name=="native" else [r["von_mises_mpa"]])
            qoi=assess(r["point_tip_y_mm"],stresses,r["volume_mm3"],r["force_n"],r["reaction_n"])
            q=r["gauss_order"]
            change=max(abs(a/b-1) for a,b in zip(stresses,prev[q])) if q in prev else None
            prev[q]=stresses
            resolution=r.get("level",r.get("h_mm"))
            peer=by_resolution[resolution].get(6 if q==8 else 8)
            quadrature_change=None
            if peer is not None:
                other=list(peer["fillet_limits_mpa"].values()) if name=="native" else [peer["von_mises_mpa"]]
                quadrature_change=max([abs(r["point_tip_y_mm"]/peer["point_tip_y_mm"]-1),*(abs(a/b-1) for a,b in zip(stresses,other))])
            stable_mesh=(name=="native" or (change is not None and change<=.03))
            ok=passed(qoi) and stable_mesh and r["relative_residual"]<=1e-8
            ok=ok and case["wall_s"]<=1800 and case["max_rss_kib"]<=6*1024**2
            results.append(dict(method=name,resolution=r.get("level",r.get("h_mm")),order=q,checks=qoi,
                                quadrature_qoi_change=quadrature_change,stress_mesh_change=change,qualified=ok))
    return {"parameters":{"radius_x_mm":3.5,"radius_y_mm":4.},"ideal_volume_mm3":VOLUME,
            "checks":results,"tetrahedral":tet,"immersed_series":immersed,"native_series":native}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("action",choices=["native","summarize"])
    p.add_argument("directory",type=Path)
    args=p.parse_args()
    if args.action=="native":
        native_series(args.directory)
    else:
        print(json.dumps(summarize(args.directory),indent=2,allow_nan=False))
