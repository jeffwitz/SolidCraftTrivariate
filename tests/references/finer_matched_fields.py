"""Test J4g decisions against newly solved C3D10 meshes at 0.9 mm."""
import argparse
import os
from pathlib import Path
import sys
import time

from adaptive_field_distance import AdaptiveField
from bracket_cost_fields import ROOT, compare, read, save
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


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=["prepare","assess","lane"])
    p.add_argument("folder",type=Path);p.add_argument("prefix",type=Path)
    p.add_argument("radii",type=int,nargs="+")
    args=p.parse_args()
    if len(args.radii)%2:raise ValueError("expected radius pairs")
    for rx,ry in zip(args.radii[::2],args.radii[1::2]):
        case=args.folder/f"x{rx}-y{ry}"
        if args.mode in ("prepare","lane"):
            prep=prepare(case,rx,ry,args.prefix)
            save(case/"preparation.json",prep);print(rx,ry,"reference",prep["qualified"],flush=True)
        if args.mode=="assess":save(case/"comparison.json",assess(case,rx,ry))
