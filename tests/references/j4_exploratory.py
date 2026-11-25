"""Replay the prospectively listed operations without changing the corpus denominator."""
import argparse
from decimal import Decimal
import json
from pathlib import Path
import subprocess
import time

from kirsch_campaign import MANIFEST


def ticks(x):
    value=2*Decimal(str(x))
    return int(value) if value==value.to_integral_value() else float(value)


def recipe(case):
    op=case["operation"]
    if op=="square_minus_disk":
        if case["height"]!=8 or case["square_half_width"]!=10:
            return None,"outside fixed plate dimensions"
        if any(ticks(x)!=int(ticks(x)) for x in [*case["center"],case["radius"]]):
            values=[Decimal(str(x))*1000 for x in [*case["center"],case["radius"]]]
            if not all(x==x.to_integral_value() for x in values):
                return None,"coordinates outside available micrometre precision"
            return {"family":"exact_bore_v1","center_um":[int(values[0]),int(values[1])],
                    "radius_um":int(values[2])},None
        return {"family":"movable_bore_v1","center_half_mm":list(map(ticks,case["center"])),
                "radius_half_mm":ticks(case["radius"])},None
    if op=="square_minus_blind_cylinder":
        if case["height"]!=8 or case["square_half_width"]!=10 or case["radius"]!=4 or case["center"]!=[0,0]:
            return None,"outside fixed pocket dimensions"
        return {"family":"pocket_depth_v1","depth_half_mm":ticks(case["depth"])},None
    if op=="square_minus_disks":
        if case["height"]!=8:
            return None,"outside row height"
        return {"family":"row_bores_v1","half_widths_half_mm":[ticks(case["square_half_width"])]*2,
                "centers_half_mm":[list(map(ticks,c)) for c in case["centers"]],
                "radii_half_mm":list(map(ticks,case["radii"]))},None
    return None,"unsupported operation"


def run(args,case):
    start=time.perf_counter()
    request,reason=recipe(case)
    output={"id":case["id"],"recipe":request}
    if request is not None:
        path=args.directory/(case["id"]+".recipe.json")
        path.write_text(json.dumps(request))
        target=args.directory/(case["id"]+".model.json")
        try:
            proc=subprocess.run([args.cli,"build",path,target],capture_output=True,text=True,timeout=60)
            result=json.loads(proc.stdout)
            if proc.returncode==0 and result.get("status")=="valid":
                audited=subprocess.run([args.cli,"audit",target],capture_output=True,text=True,timeout=60)
                reloaded=json.loads(audited.stdout)
                output["status"]="qualified" if audited.returncode==0 and reloaded==result else "execution_error"
                output["certificate"]=result
            else:
                output["status"]="refused"
                reason=result.get("reason","missing diagnostic")
        except subprocess.TimeoutExpired:
            output["status"]="timeout"
        except (OSError,ValueError) as error:
            output["status"]="execution_error"
            reason=str(error)
    else:
        output["status"]="refused"
    if reason:
        output["reason"]=reason
    output["elapsed_s"]=time.perf_counter()-start
    return output


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path)
    p.add_argument("--cli",default="target/debug/solidcraft-trivariate")
    args=p.parse_args()
    args.directory.mkdir(parents=True,exist_ok=True)
    manifest=json.loads(MANIFEST.read_text())
    result={"version":manifest["version"]}
    for group in ("exploratory","invalid_controls"):
        result[group]=[run(args,c) for c in manifest[group]]
    print(json.dumps(result,indent=2))
