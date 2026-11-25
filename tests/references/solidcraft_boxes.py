"""Compare actual SolidCraft box commands with the native certified cartesian program."""
import FreeCAD as App
import Part

import argparse
import copy
import hashlib
from itertools import product
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from orthogonal_program import native_measures, voxel
from solidcraft_bridge import bounds

HERE = Path(__file__).resolve().parent


def script_for(recipe):
    commands = []
    sequence = [("base","new",recipe["base_mm"])]
    sequence += [(f["id"], "join" if f["operation"]=="add" else "cut", f["bounds_mm"])
                 for f in recipe["features"]]
    for index, (name, operation, box) in enumerate(sequence):
        params = dict(zip(("length","width","height"), [b-a for a,b in box]))
        if index % 2:
            params["center"] = [(a+b)/2 for a,b in box]
        else:
            params["corner"] = [a for a,b in box]
        params.update(name=name, operation=operation)
        commands.append({"command":"solid.box","params":params})
    return {"commands":commands}


def process(*args):
    begin=time.monotonic()
    cmd=subprocess.run(list(map(str,args)),capture_output=True,text=True)
    return cmd,time.monotonic()-begin


def shapes_report(shapes, recipe, volume, areas, expected):
    all_bounds=[recipe["base_mm"],*(f["bounds_mm"] for f in recipe["features"])]
    axes=[]
    for d in range(3):
        planes=sorted({n for box in all_bounds for n in box[d]})
        axes.append([planes[0]-.5,*[(a+b)/2 for a,b in zip(planes,planes[1:])],planes[-1]+.5])
    mismatches=0
    for xyz in product(*axes):
        idx=tuple(np.floor(np.array(xyz)-expected["origin"]).astype(int))
        ref=bool(expected["solid"][idx])
        point=App.Vector(*xyz)
        mismatches+=int(any(s.isInside(point,1e-8,False)!=ref for s in shapes))
    verror=max(abs(s.Volume/volume-1) for s in shapes)
    aerror=max(abs(s.Area/sum(areas.values())-1) for s in shapes)
    boxes=np.asarray([bounds(s) for s in shapes])
    return dict(single_valid=[s.isValid() and len(s.Solids)==1 and len(s.Shells)==expected["cavities"]+1 for s in shapes],
                classified_points=int(np.prod(list(map(len,axes)))),classification_mismatches=mismatches,
                volume_relative_error=verror,area_relative_error=aerror,max_bbox_difference_mm=float(np.ptp(boxes,axis=0).max()))


def run(cli,host,work):
    work.mkdir(parents=True,exist_ok=True)
    rows=[]
    baseline=work/"last-admitted.json"
    items=json.loads((HERE/"orthogonal_corpus.json").read_text())
    for item in items:
        folder=work/item["id"];folder.mkdir(exist_ok=True)
        recipe=item["recipe"]; src=script_for(recipe)
        script=folder/"script.json";script.write_text(json.dumps(src))
        editable=folder/"editable.sc";hoststep=folder/"host.step"
        hp,host_s=process(host,"run",script,"--param-check","strict","--save",editable,"--out",hoststep,
                          "--json-out",folder/"host-run.json","--results")
        before=baseline.read_bytes() if baseline.exists() else None
        cp,native_s=process(cli,"import-solidcraft",script,baseline)
        native=json.loads(cp.stdout)
        row=dict(id=item["id"],script=src,host_accepted=hp.returncode==0,host_note=hp.stderr.strip(),
                 native=native,expected=item["admitted"],host_elapsed_s=host_s,native_elapsed_s=native_s)
        if not item["admitted"]:
            unchanged=before==(baseline.read_bytes() if baseline.exists() else None)
            row.update(previous_unchanged=unchanged,qualified=bool(cp.returncode and unchanged))
            if not hp.returncode:
                row["host_solid_count"]=len(Part.read(str(hoststep)).Solids)
        elif hp.returncode or cp.returncode:
            row["qualified"]=False
        else:
            reopened=folder/"reopened.step";native_step=folder/"native.step"
            process_result,_=process(host,"run","--in",editable,"--out",reopened)
            process_result.check_returncode()
            export,_=process(cli,"export-step",baseline,native_step)
            export.check_returncode()
            expected=voxel(recipe)
            vol,areas=native_measures(json.loads(baseline.read_text()))
            evidence=shapes_report([Part.read(str(s)) for s in (native_step,hoststep,reopened)],recipe,vol,areas,expected)
            row.update(evidence,volume_mm3=vol,exact_design_areas=areas==expected["areas"],
                       saved_document_sha256=hashlib.sha256(editable.read_bytes()).hexdigest())
            row["qualified"]=bool(all(evidence["single_valid"]) and evidence["classification_mismatches"]==0
                                  and evidence["volume_relative_error"]<1e-9 and evidence["area_relative_error"]<1e-9
                                  and evidence["max_bbox_difference_mm"]<1e-6 and row["exact_design_areas"] and vol==expected["volume"])
        print(item["id"],row["qualified"],flush=True);rows.append(row)
    base=script_for(items[1]["recipe"]);refusals=[]
    for name,extra in (("explicit_unresolved_target",{"targets":["unresolved"]}),("new_body",{"operation":"new"}),
                       ("unresolved_expression",{"length":"width_var + 1"}),("nonintegral_plane",{"center":[1.25,2,3]}),
                       ("unrecognized_option",{"ignored_geometry":2}),("unsupported_intersection",{"operation":"intersect"})):
        script=copy.deepcopy(base)
        script["commands"][-1]["params"].update(extra)
        src=work/(name+".json");src.write_text(json.dumps(script))
        before=baseline.read_bytes()
        cp,_=process(cli,"import-solidcraft",src,baseline)
        refusals.append(dict(id=name,diagnostic=json.loads(cp.stdout),qualified=bool(cp.returncode and before==baseline.read_bytes())))
    return dict(schema="solidcraft_box_program_parity_v1",solidcraft=process(host,"--version")[0].stdout.strip(),
                occ=Part.OCC_VERSION, cases=rows, additional_refusals=refusals)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","host","work","output"):
        parser.add_argument(k,type=Path)
    a=parser.parse_args()
    report=run(a.cli.resolve(),a.host.resolve(),a.work)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    if not all(r["qualified"] for r in report["cases"]+report["additional_refusals"]):
        raise SystemExit(1)
