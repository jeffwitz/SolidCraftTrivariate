"""Test J4g decisions against newly solved C3D10 meshes at 0.9 mm."""
import argparse
import os
from pathlib import Path
import sys
import time

from adaptive_field_distance import AdaptiveField
from bracket_cost_fields import ROOT, compare, read, save
from bracket_serial_cost import MANIFEST, metrics
from bracket_field_distance import ImmersedField,TetraField
from iterative_tetra_fields import defect_corrected
from refine_tetra_fields import HERE,measure


def prepare(folder,rx,ry,prefix):
    folder.mkdir(parents=True,exist_ok=True)
    solver=prefix/"bin/ccx"
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="2",
             LD_LIBRARY_PATH=str(prefix/"lib"),PYTHONPATH=str(prefix/"lib"))
    start=time.monotonic()
    mesh=measure([str(prefix/"bin/python"),str(HERE/"calculix_bracket.py"),"--ccx",str(solver),
        "--deck-only","--store-displacements","--radius-x",str(rx/16),"--radius-y",str(ry/16),
        "--size","0.9",str(folder/"mesh")],env,folder/"mesh.log",seconds=1800)
    out=dict(radii_ticks=[rx,ry],size_mm=.9,mesh=mesh,qualified=False)
    if mesh["status"]=="completed":
        out["solve"]=defect_corrected(folder/"mesh",folder/"fine",solver,env,max(0,1800-(time.monotonic()-start)),5)
        out["mesh_and_solve_s"]=time.monotonic()-start
        if out["solve"]["status"]=="completed":
            out["probe_run"]=measure([sys.executable,str(HERE/"fillet_probe.py"),str(folder/"fine"),str(folder/"probe"),
                    "--ccx",str(solver)],dict(os.environ,OPENBLAS_NUM_THREADS="1"),folder/"probe.log")
            if out["probe_run"]["status"]=="completed":
                out["probe"]=read(folder/"probe/result.json")
                out["qualified"]=bool(out["solve"]["result"]["reaction_relative_error"]<=1e-5 and
                                      out["mesh_and_solve_s"]<=1800 and out["probe"]["distance_mm"]<.001)
    return out


def assess(folder,rx,ry):
    case=f"x{rx}-y{ry}";coeffs=ROOT/"build/j4e"/case
    blocks=read(coeffs/"native.package.json")["analysis_space"]["blocks"]
    tier=ROOT/"build/j4f"/case
    fields={"native_local":AdaptiveField(tier/"s2.bin"),"immersed_h125":ImmersedField(coeffs/"immersed.npz")}
    scales={"native_local":read(tier/"s2-q8.json")["q8"]["strain_energy_n_mm"],
            "immersed_h125":read(coeffs/"immersed.json")["strain_energy_n_mm"]}
    for key,path in (("tetra_h15",ROOT/"build/bracket-cost"/case/"tetra/h1.5-qNone"),
                     ("tetra_h1",coeffs/"h1"),("tetra_h09",folder/"fine")):
        fields[key]=TetraField(path);ref=read(path/"result.json")
        scales[key]=ref.get("strain_energy_n_mm",-500*ref["mean_tip_y_mm"])
    return compare(blocks,fields,scales,11,reference_keys=["tetra_h09"])


RADII = [(40,80),(50,74),(80,40),(51,74)]
FIELD_NAMES = {"local":"native_local","tetra":"tetra_h15","immersed":"immersed_h125"}


def measured_assessment(root,rx,ry,prefix):
    case=root/f"x{rx}-y{ry}"
    prep=read(case/"preparation.json")
    if prep.get("solve",{}).get("status")!="completed":
        return
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
    child=[sys.executable,str(Path(__file__).resolve()),"assess",str(root),str(prefix),str(rx),str(ry)]
    proc=measure(child,env,case/"field.log")
    save(case/"field-process.json",proc)
    previous=ROOT/f"build/j4e/x{rx}-y{ry}/h1"
    old=measure([sys.executable,str(HERE/"fillet_probe.py"),str(previous),str(case/"probe-h1"),
                 "--ccx",str(prefix/"bin/ccx")],env,case/"probe-h1.log")
    save(case/"probe-h1-process.json",old)


def aggregate(root,output,original_root=None):
    replay=read(ROOT/"docs/results/bracket-field-matched-cost-20261125.json")
    qoi=read(MANIFEST)["qoi"]
    rows=[]
    for rx,ry in RADII:
        folder=root/f"x{rx}-y{ry}"
        source="initial_1800s"
        if original_root is not None:
            if (folder/"preparation.json").exists():source="continued_3600s_cumulative"
            else:folder=original_root/f"x{rx}-y{ry}"
        row=dict(radii_ticks=[rx,ry],reference_status="pending",matched_trials=[],reference_source=source)
        rows.append(row)
        if not (folder/"preparation.json").exists():continue
        prep=read(folder/"preparation.json")
        row.update(preparation=prep,reference_status="qualified" if prep["qualified"] else "unqualified")
        if source=="continued_3600s_cumulative":
            children=prep.get("solve",{}).get("passes",[])
            bound=prep["continued_process"]["tree_peak_rss_kib"]+max(
                (p["tree_peak_rss_kib"] for p in children if not p.get("reused")),default=0)
            row["continuation_tree_rss_upper_bound_kib"]=bound
            row["rss_accounting"]="sum of separately sampled parent-session and largest non-overlapping child-session peaks"
            if bound>6*1024**2:
                row["reference_status"]="unqualified"
        if not (folder/"field-process.json").exists():continue
        row["field_process"]=read(folder/"field-process.json")
        if row["field_process"]["status"]!="completed":continue
        comp=read(folder/"comparison.json")
        row["field_comparison"]=comp
        if not (folder/"probe-h1-process.json").exists():continue
        old_proc=read(folder/"probe-h1-process.json");row["previous_probe_process"]=old_proc
        if old_proc["status"]!="completed" or "probe" not in prep:continue
        old_tip=read(ROOT/f"build/j4e/x{rx}-y{ry}/h1/result.json")["point_tip_y_mm"]
        old_probe=read(folder/"probe-h1/result.json");row["previous_reference_probe"]=old_probe
        new_tip=prep["solve"]["result"]["point_tip_y_mm"];new_vm=prep["probe"]["von_mises_mpa"]
        row["reference_drift_relative"]=dict(tip=abs(old_tip/new_tip-1),fillet=abs(old_probe["von_mises_mpa"]/new_vm-1))
        for batch in replay["rows"]:
            if batch["radii_ticks"]!=[rx,ry]:continue
            for trial in batch["trials"]:
                method=trial["method"]
                m=metrics(method,trial)
                dist=comp["levels"][-1]["distances"][FIELD_NAMES[method]+"__tetra_h09"]
                tip=abs(m["tip_mm"]/new_tip-1);vm=abs(m["stress_mpa"]/new_vm-1)
                within=(trial["qualified"] and row["reference_status"]=="qualified" and comp["recovered_energies"] and
                        comp["all_points_located_once"] and tip<=qoi["tip_relative_to_finest_tetra"]
                        and vm<=qoi["fillet_von_mises_relative_to_finest_tetra"] and
                        all(dist[k]["stable"] and dist[k]["mean"]<=v for k,v in (("energy",.03),("l2",.005))))
                row["matched_trials"].append(dict(repeat=batch["repeat"],method=method,tip_relative=tip,
                    fillet_relative=vm,field=dist,qualified=bool(within)))
    done=all(r["reference_status"]!="pending" and (r["reference_status"]=="unqualified" or
            "previous_probe_process" in r or r.get("field_process",{}).get("status") in
            ("timeout","rss_budget_exceeded","execution_error")) for r in rows)
    save(output,dict(schema="bracket_finer_reference_extended_v1" if original_root else "bracket_finer_reference_v1",base_revision="40cfd98" if original_root else "4d5ae4c",qualified_reference_count=
         sum(r["reference_status"]=="qualified" for r in rows),matched_trials_qualified=sum(
         t["qualified"] for r in rows for t in r["matched_trials"]),mandatory_trials=24,finished=done,cases=rows))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=["prepare","assess","lane","measured","aggregate","aggregate-extended"])
    p.add_argument("folder",type=Path);p.add_argument("prefix",type=Path)
    p.add_argument("radii",type=int,nargs="*")
    args=p.parse_args()
    if args.mode in ("aggregate","aggregate-extended"):
        aggregate(args.folder,args.prefix,ROOT/"build/j4j" if args.mode=="aggregate-extended" else None);sys.exit(0)
    if not args.radii:raise ValueError("expected radius pairs")
    if len(args.radii)%2:raise ValueError("expected radius pairs")
    for rx,ry in zip(args.radii[::2],args.radii[1::2]):
        case=args.folder/f"x{rx}-y{ry}"
        if args.mode in ("prepare","lane"):
            prep=prepare(case,rx,ry,args.prefix)
            save(case/"preparation.json",prep);print(rx,ry,"reference",prep["qualified"],flush=True)
        if args.mode=="measured":measured_assessment(args.folder,rx,ry,args.prefix)
        if args.mode=="assess":save(case/"comparison.json",assess(case,rx,ry))
