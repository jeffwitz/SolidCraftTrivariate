"""Collect measured full-field comparisons for the four preregistered brackets."""
import argparse
import os
from pathlib import Path
import sys
import time

from bracket_cost_fields import read, save
from refine_tetra_fields import measure

HERE = Path(__file__).resolve().parent
CASES = read(HERE/"bracket_cost_corpus.json")["mandatory_radii_sixteenths_mm"]


def collect(kind, inputs, tiers, serial, summary, output, wait):
    rows = []
    deadline = time.monotonic()+wait
    for rx, ry in CASES:
        folder = (inputs if kind=="serial" else tiers)/f"x{rx}-y{ry}"
        result, resources = folder/"comparison.json", folder/"comparison-resources.json"
        if not resources.exists() or not result.exists():
            while True:
                prepared = read(inputs/"preparation.json")
                ready = any(r["radii_ticks"]==[rx,ry] and r.get("fine_qualified") for r in prepared)
                ready &= (inputs/f"x{rx}-y{ry}/native.package.json").exists() if kind=="tiers" else True
                if kind=="tiers":
                    ready &= sum(s["radii_ticks"]==[rx,ry] for s in read(summary)["states"])==4
                if ready or time.monotonic()>=deadline:break
                time.sleep(15)
            if not ready:
                rows.append(dict(radii_ticks=[rx,ry],status="inputs_unavailable"))
                continue
            if kind=="serial":
                cmd=[sys.executable,HERE/"bracket_cost_fields.py",inputs,result,rx,ry,"--serial-report",serial]
            else:
                cmd=[sys.executable,HERE/"bracket_field_candidates.py",inputs,tiers,summary,result,rx,ry]
            actual = measure(list(map(str,cmd)),dict(os.environ),folder/"compare.log")
            save(resources,actual)
        actual = read(resources)
        record = dict(radii_ticks=[rx,ry],status=actual["status"],resources=actual)
        if actual["status"]=="completed":record["result"]=read(result)
        rows.append(record)
        save(output,dict(schema="bracket_physical_field_review_v1",series=kind,mandatory=CASES,cases=rows,finished=False))
        print(kind,rx,ry,actual,flush=True)
    save(output,dict(schema="bracket_physical_field_review_v1",series=kind,mandatory=CASES,cases=rows,finished=True))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("kind",choices=("serial","tiers"))
    for name in ("inputs","tiers","serial","summary","output"):p.add_argument(name,type=Path)
    p.add_argument("--wait-inputs-seconds",type=int,default=0)
    a=p.parse_args()
    collect(a.kind,a.inputs,a.tiers,a.serial,a.summary,a.output,a.wait_inputs_seconds)
