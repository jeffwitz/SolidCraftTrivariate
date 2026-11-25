"""Replay interval edits cold and check exact axial crossing times and surviving selections."""
import argparse
from collections import defaultdict
import copy
from fractions import Fraction
import json
import math
from pathlib import Path
import random
import subprocess
import time

HERE=Path(__file__).resolve().parent


def write(path,data):path.write_text(json.dumps(data))
def read(path):return json.loads(path.read_text())


def shape_class(interval,height):
    return {0:"enclosed",1:"open_start",2:"open_end",3:"through"}[
        int(interval[0]==0)+2*int(interval[1]==height)]


def events(old,id,next_interval):
    limits={"start":0,"end":old["height_mm"]}
    for b in old["bores"]:
        if b["id"]==id:
            previous=b["z_interval_mm"]
        else:
            for end,z in zip(("lower","upper"),b["z_interval_mm"]):
                limits[f"bore/{b['id']}/{end}"]=z
    expected=[]
    for end,a,b in zip(("lower","upper"),previous,next_interval):
        if a!=b:
            for tag,z in limits.items():
                t=Fraction(z-a,b-a)
                if 0<=t<=1:expected.append((t,f"bore/{id}/{end}",tag,z))
    return sorted(expected)


def selector_map(model):
    faces=defaultdict(set)
    for f in model["design_faces"]:
        faces[f["feature"]].add((f["block"],f["side"]))
    return faces


def slab_keys(r):
    cuts=sorted({0,r["height_mm"],*[z for b in r["bores"] for z in b["z_interval_mm"]]})
    return {(a,b,tuple(sorted(f["id"] for f in r["bores"] if f["z_interval_mm"][0] <= (a+b)/2 < f["z_interval_mm"][1])))
            for a,b in zip(cuts,cuts[1:])}


def run(cli,work):
    work.mkdir(exist_ok=True,parents=True)
    mixed=next(x["recipe"] for x in read(HERE/"interval_bores_corpus.json") if x["id"]=="mixed")
    rng=random.Random(2026112522);cases=[];refusals=[]
    recipe,old,cold,edits,out=(work/(k+".json") for k in ("recipe","old","cold","edits","out"))
    def call(*args):
        p=subprocess.run([cli,*args],text=True,capture_output=True)
        return p.returncode,json.loads(p.stdout)
    for count in range(1,5):
        r=copy.deepcopy(mixed);r["bores"]=r["bores"][:count];h=r["height_mm"]
        r["bores"][0]["z_interval_mm"]=[2,h-2]
        write(recipe,r);assert call("build",recipe,old)[0]==0
        for index in range(24):
            key=0 if index<6 else rng.randrange(count)
            bore=r["bores"][key];id=bore["id"];prev=bore["z_interval_mm"]
            target=([0,h-2],[0,h],[2,h-2],[2,h-2],[3,h-2],[3,h-1])[index] if index<6 else sorted(rng.sample(range(h+1),2))
            step=dict(bore=id,z_interval_mm=target,require_features=["start","end",f"bore/{id}/wall"])
            write(edits,dict(schema="interval_bore_edits_v1",steps=[step]))
            old_model=read(old)
            before=time.monotonic();code,result=call("edit-bore-intervals",old,edits,out);warm_time=time.monotonic()-before
            new=copy.deepcopy(r);new["bores"][key]["z_interval_mm"]=target;write(recipe,new)
            before=time.monotonic();assert call("build",recipe,cold)[0]==0;cold_time=time.monotonic()-before
            record=dict(bores=count,index=index,id=id,old=prev,new=target,
                        cli_edit_s=warm_time,cli_cold_s=cold_time,diagnostic=result)
            if code==0:
                rep=result["steps"][0];candidate=read(out)
                oldmap,newmap=selector_map(old_model),selector_map(candidate)
                expected_events=events(r,id,target)
                observed=[(Fraction(*e["time"]),e["moving_cap"],e["fixed_plane"],e["height_mm"]) for e in rep["events"]]
                assert all(math.gcd(*e["time"])==1 and e["time"][1]>0 for e in rep["events"])
                expected_hits=slab_keys(r)&slab_keys(new)
                ranges={(lo,hi) for lo,hi,_ in expected_hits}
                def span(leaf):return tuple(map(int,leaf["id"].rsplit("_z",1)[1].split("_")))
                block_hits=sum(span(b) in ranges for b in candidate["blocks"])
                checks=dict(cold_model_identical=out.read_bytes()==cold.read_bytes(),
                    event_times=observed==expected_events,classes=(rep["old_class"],rep["new_class"])==(shape_class(prev,h),shape_class(target,h)),
                    volume=rep["volume_change_pi_mm3"]==bore["radius_mm"]**2*((prev[1]-prev[0])-(target[1]-target[0])),
                    lost=set(rep["lost_features"])==oldmap.keys()-newmap.keys(),
                    created=set(rep["created_features"])==newmap.keys()-oldmap.keys(),
                    remapped=set(rep["remapped_features"])=={k for k in oldmap.keys()&newmap.keys() if oldmap[k]!=newmap[k]},
                    reused_layers=rep["reused_layers"]==len(expected_hits),
                    rebuilt_layers=rep["built_layers"]+rep["reused_layers"]==len(slab_keys(new)),
                    reused_blocks=rep["reused_blocks"]==block_hits,
                    rebuilt_blocks=rep["built_blocks"]+block_hits==len(candidate["blocks"]))
                record.update(checks=checks,passed=all(checks.values()))
            else:record["passed"]=False
            cases.append(record)
            old.write_bytes(cold.read_bytes());r=new
        print("group",count,sum(x["passed"] for x in cases if x["bores"]==count),flush=True)
        id=r["bores"][0]["id"]
        first=dict(bore=id,z_interval_mm=[2,h-2],require_features=["start"])
        for name,last in [("lost_cap",dict(bore=id,z_interval_mm=[0,h],require_features=[f"bore/{id}/lower"])),
                          ("unknown_bore",dict(bore="missing",z_interval_mm=[0,h],require_features=[])),
                          ("reversed",dict(bore=id,z_interval_mm=[h,0],require_features=[]))]:
            write(edits,dict(schema="interval_bore_edits_v1",steps=[first,last]))
            marker=old.read_bytes();out.write_bytes(marker)
            code,msg=call("edit-bore-intervals",old,edits,out)
            refusals.append(dict(bores=count,case=name,diagnostic=msg,passed=code!=0 and out.read_bytes()==marker))
    return dict(schema="interval_bore_edit_oracle_v1",seed=2026112522,cases=cases,refusals=refusals,
                qualified=all(x["passed"] for x in [*cases,*refusals]))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("cli","work","output"):parser.add_argument(name,type=Path)
    args=parser.parse_args();d=run(args.cli.resolve(),args.work)
    args.output.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    raise SystemExit(int(not d["qualified"]))
