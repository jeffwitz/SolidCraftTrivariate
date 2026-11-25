"""Run the preregistered bracket cost corpus and preserve each completed trial."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys

from refine_tetra_fields import measure

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
MANIFEST=HERE/"bracket_cost_corpus.json"


def save(path,data):
    tmp=path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    tmp.replace(path)


def run_trial(case,method,h,q,folder,freecad):
    folder.mkdir(parents=True,exist_ok=True)
    result_file=folder/"trial.json"
    if result_file.exists():
        return json.loads(result_file.read_text())
    rx,ry=case
    env=dict(os.environ,OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="1")
    ccx=freecad/"bin/ccx"
    by_radii=["--radius-x",str(rx/16),"--radius-y",str(ry/16)]
    if method=="tetra":
        env.update(LD_LIBRARY_PATH=str(freecad/"lib"),PYTHONPATH=str(freecad/"lib"))
        cmd=[freecad/"bin/python",HERE/"calculix_bracket.py","--ccx",ccx,"--size",str(h),
             "--store-displacements",*by_radii,folder]
    elif method=="native":
        cmd=[ROOT/"target/release/examples/bracket_convergence",str(h),str(q),str(rx),str(ry)]
    else:
        cmd=[sys.executable,HERE/"immersed_bracket.py","--solver","multigrid","--h",str(h),"--order",str(q),*by_radii]
    elapsed=folder/"resources.txt"
    cmd=["/usr/bin/time","-f","%e %M","-o",str(elapsed),*map(str,cmd)]
    row=measure(cmd,env,folder/"stdout.log",cwd=ROOT)
    row.update(method=method,resolution=h,q=q,radii_ticks=case)
    if elapsed.exists():
        fields=elapsed.read_text().split()
        if fields and fields[-1].isdigit():
            row["largest_process_hwm_kib"]=int(fields[-1])
    if row["status"]=="completed":
        if method=="tetra":
            row["result"]=json.loads((folder/"result.json").read_text())
            probe_env=dict(os.environ,OPENBLAS_NUM_THREADS="1")
            probe=measure([sys.executable,str(HERE/"fillet_probe.py"),str(folder),str(folder/"probe"),"--ccx",str(ccx)],
                          probe_env,folder/"probe.log",seconds=max(0.,1800-row["elapsed_s"]))
            row["postprocessing"]=probe
            if probe["status"]=="completed":
                row["probe"]=json.loads((folder/"probe/result.json").read_text())
            else:
                row["status"]=probe["status"]
        else:
            row["result"]=json.loads((folder/"stdout.log").read_text())
    else:
        row["diagnostic"]=(folder/"stdout.log").read_text()[-1200:]
    save(result_file,row)
    print(case,method,h,q,row["status"],round(row["elapsed_s"],2),flush=True)
    return row


def summarize(root,m):
    cases=[]
    for case in m["mandatory_radii_sixteenths_mm"]:
        d=root/("x%d-y%d"%tuple(case))
        runs=[json.loads(p.read_text()) for p in sorted(d.glob("*/*/trial.json"))]
        table={}
        for method in ("tetra","native","immersed"):
            order=(None,) if method=="tetra" else tuple(m[method+"_orders"])
            levels=m["tetra_sizes_mm"] if method=="tetra" else m["native_levels" if method=="native" else "immersed_sizes_mm"]
            table[method]={(h,q):next((r for r in runs if r["method"]==method and r["resolution"]==h and r["q"]==q),None)
                           for h in levels for q in order}
        summary={"radii_ticks":case,"runs":runs,"qualification":{}}
        cases.append(summary)
        ref=table["tetra"].get((1.25,None))
        prevref=table["tetra"].get((1.5,None))
        if not ref or not prevref or any(x["status"]!="completed" for x in [ref,prevref]):
            summary["reference_status"]="incomplete"
            continue
        def obs(r):
            v=r["result"]
            stress=v.get("fillet_limits_mpa",{})
            sig={}
            if stress:
                for key,value in stress.items():
                    sig.setdefault(key.split("/")[0],[]).append(value)
            else:
                sig["physical_probe"]=[r["probe"]["von_mises_mpa"] if r["method"]=="tetra" else v["von_mises_mpa"]]
            return v["point_tip_y_mm"],sig
        def change(a,b):
            return abs(a/b-1)
        def side_change(current,previous):
            if current.keys()!=previous.keys():
                raise ValueError("Fillet trace supports changed")
            return max(change(a,b) for key in current for a in current[key] for b in previous[key])
        rt,rs=obs(ref)
        old,olds=obs(prevref)
        limits=m["qoi"]
        expected=34080-math.pi*(360+8*sum((a/16)**2 for a in case))
        reference_ok=(change(rt,old)<=limits["successive_tip_relative_change"] and
                      side_change(rs,olds)<=limits["successive_fillet_relative_change"] and
                      (math.dist(ref["result"]["reaction_n"],[0.,1000.,0.])/1000.)<=limits["reaction_relative"])
        summary["reference_status"]="stable" if reference_ok else "indeterminate"
        for method,ladder in table.items():
            first=None
            cumulative=0.
            previous={}
            for (h,q),r in ladder.items():
                if r is None:continue
                elapsed=r["elapsed_s"]+r.get("postprocessing",{}).get("elapsed_s",0.)
                cumulative+=elapsed
                if r["status"]!="completed":continue
                t,s=obs(r);v=r["result"]
                force=v.get("reaction_n")
                vol=v.get("volume_mm3",v.get("occ_volume_mm3"))
                reaction=math.dist(force,[0.,1000.,0.])/1000.
                p=previous.get(q)
                stable=p is not None and change(t,p[0])<=limits["successive_tip_relative_change"] and side_change(s,p[1])<=limits["successive_fillet_relative_change"]
                previous[q]=(t,s)
                peer=ladder.get((h,6)) if q==8 else None
                quad=None
                if peer and peer["status"]=="completed":
                    pt,ps=obs(peer)
                    quad=max(change(t,pt),side_change(s,ps))
                checks=dict(tip_relative=change(t,rt),stress_relative=max(change(x,rs["physical_probe"][0]) for side in s.values() for x in side),
                            volume_relative=change(vol,expected),reaction_relative=reaction,mesh_stable=stable,quadrature_relative=quad)
                meet=(reference_ok and stable and checks["tip_relative"]<=limits["tip_relative_to_finest_tetra"]
                      and checks["stress_relative"]<=limits["fillet_von_mises_relative_to_finest_tetra"]
                      and checks["volume_relative"]<=limits["volume_relative"] and reaction<=limits["reaction_relative"]
                      and (method=="tetra" or (quad is not None and quad<=limits["quadrature_qoi_relative"])))
                r.update(checks=checks,qualified=bool(meet),physical_and_postprocess_s=elapsed)
                if meet and first is None:
                    first=dict(resolution=h,q=q,level_wall_s=elapsed,search_wall_s=cumulative,
                               tree_peak_rss_kib=max(r["tree_peak_rss_kib"],r.get("postprocessing",{}).get("tree_peak_rss_kib",0)))
            summary["qualification"][method]=first
    return dict(manifest=m,summary_revision=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),cases=cases)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("method",choices=["native","immersed","tetra","summarize"])
    p.add_argument("directory",type=Path)
    p.add_argument("--freecad",type=Path,default=Path("/home/jeff/Code/FreeCAD/.pixi/envs/default"))
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    m=json.loads(MANIFEST.read_text())
    root=a.directory.resolve()
    if a.method=="summarize":
        save(a.output,summarize(root,m))
        return
    name=a.method
    levels=m["native_levels" if name=="native" else name+"_sizes_mm"]
    for case in m["mandatory_radii_sixteenths_mm"]:
        for h in levels:
            for q in m.get(name+"_orders",[None]):
                folder=root/("x%d-y%d"%tuple(case))/name/f"h{h}-q{q}"
                run_trial(case,name,h,q,folder,a.freecad.resolve())


if __name__=="__main__":
    main()
