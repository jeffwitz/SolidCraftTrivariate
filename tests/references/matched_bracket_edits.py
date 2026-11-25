"""Compare transactional, locally refined radius edits with the J4f full fields."""
import argparse
import json
import os
from pathlib import Path
import numpy as np

from bracket_checkpoint_geometry import checkpoint, verify_geometry
from refine_tetra_fields import measure

ROOT = Path(__file__).resolve().parents[2]
SEQUENCE = [(50,74),(51,74),(50,74),(81,74),(40,80),(80,40)]


def plan_for(names):
    pending = {}
    for name in names:
        root, *children = name.split("/")
        xyz = (0, 0, 0)
        for depth, bits in enumerate(children):
            pending.setdefault(depth,set()).add((root,xyz))
            xyz = tuple(2*v+int(b) for v,b in zip(xyz,bits))
    return dict(schema="adaptive_refinement_plan_v1",max_cells=3000,max_depth=7,
                steps=[dict(refine=[dict(root=r,level=d,index=xyz) for r,xyz in sorted(values)])
                       for d,values in sorted(pending.items())])


def events(folder, tiers, cold):
    rows = []
    for i, (rx,ry) in enumerate(SEQUENCE):
        path = folder / f"event-{i}.json"
        if not path.exists():
            rows.append(dict(radii_ticks=[rx,ry],status="missing",passed=False))
            continue
        row = json.loads(path.read_text())
        rows.append(row)
        if row["status"] == "refused":
            row["passed"] = (rx == 81 and row["prior_state_preserved"] and row["cache_preserved"])
            continue
        case = f"x{rx}-y{ry}"
        field = folder/f"event-{i}.bin"
        ha,a = checkpoint(field)
        hb,b = checkpoint(tiers/case/"s2.bin")
        root = json.loads((ROOT/"build/j4e"/case/"native.package.json").read_text())
        ids,bad = verify_geometry(root,field.with_suffix(".geometry.bin"))
        same = ha==hb and a.keys()==b.keys()==ids and bad is None
        diff = np.array([a[k][:,1:]-b[k][:,1:] for k in b]) if same else None
        error = float(np.linalg.norm(diff)/np.linalg.norm([v[:,1:] for v in b.values()])) if same else None
        same_weights = same and all(np.array_equal(a[k][:,0],b[k][:,0]) for k in b)
        m=row["metrics"]
        ref = cold[(rx,ry)]["q6_q8"][1]
        fillet = row["fillet_mpa"]
        stress_error = max(abs(v/ref["fillet_limits_mpa"][k]-1) for k,v in fillet.items())
        tip_error = abs(m["tip_displacement_mm"][1]/ref["tip_mm"][1]-1)
        energy_error = abs(m["energy_n_mm"]/ref["strain_energy_n_mm"]-1)
        row.update(fillet_replay_error=stress_error,tip_replay_error=tip_error,
                   energy_replay_error=energy_error,field_coefficients_relative=error,
                   identical_weights=same_weights,exact_geometry=same,leaves=len(a))
        row["passed"] = bool(rx!=81 and same_weights and len(a)==2528 and
                             m["cold_relative_residual"]<=1e-10 and (m["warm_relative_residual"] or 0)<=1e-10 and
                             m["stiffness_action_relative_error"]<=1e-11 and
                             max(stress_error,tip_error,energy_error,error)<=1e-8 and
                             abs(m["load_area_mm2"]/320-1)<1e-8 and
                             (not i or m["reused_contributions"]>0))
    return rows


def run(work,output):
    plans=work/"plans"
    plans.mkdir(parents=True,exist_ok=True)
    tiers=ROOT/"build/j4f"
    cold={tuple(s["radii_ticks"]):s for s in json.loads(
        (ROOT/"docs/results/bracket-field-tiers-20261125.json").read_text())["states"] if s["step"]==2}
    for (rx,ry) in cold:
        case=f"x{rx}-y{ry}"
        _,field=checkpoint(tiers/case/"s2.bin")
        (plans/(case+".json")).write_text(json.dumps(plan_for(field)))
    out = dict(schema="matched_bracket_edit_verification_v1",cache_cap_mib=384,
               max_cells=3000,cpu_affinity=[2,3],trials=[],finished=False)
    env=dict(os.environ,OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="1")
    for repeat in (1,2):
        folder=work/f"repeat-{repeat}"
        folder.mkdir(exist_ok=False)
        cmd=["/usr/bin/taskset","-c","2,3",str(ROOT/"target/release/examples/bracket_edit_cost"),
             "--verified-plan",str(plans),str(folder)]
        proc=measure(cmd,env,folder/"run.log",seconds=1800)
        rows=events(folder,tiers,cold)
        out["trials"].append(dict(repeat=repeat,process=proc,command=cmd,events=rows,
                                  passed=proc["status"]=="completed" and all(r["passed"] for r in rows)))
        output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
        print(repeat,proc["status"],[(r["radii_ticks"],r["passed"]) for r in rows],flush=True)
    out["finished"]=True
    out["qualified"]=all(t["passed"] for t in out["trials"])
    output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    return out["qualified"]


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    raise SystemExit(int(not run(a.work.resolve(),a.output)))
