"""Reopen real editable box documents and replay their archived SolidCraft operation stream."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess


def run(cli,host,work):
    source=json.loads((Path(__file__).parent/"orthogonal_corpus.json").read_text())
    rows=[]
    for item in source:
        if not item["admitted"]:
            continue
        root=work/item["id"]
        script=root/"script.json"; archived=root/"canonical.json"
        subprocess.run([host,"script",root/"editable.sc","--out",archived],check=True,capture_output=True)
        outputs=[]
        for path,label in ((script,"direct"),(archived,"archive")):
            dst=root/(label+".native.json")
            subprocess.run([cli,"import-solidcraft",path,dst],check=True,capture_output=True)
            outputs.append(dst.read_bytes())
        document=json.loads(archived.read_text())
        if outputs[0]!=outputs[1]:
            raise ValueError(f"Archive/direct model mismatch: {item['id']}")
        rows.append(dict(id=item["id"],steps=len(document["commands"]),identical_bytes=len(outputs[0]),
                         has_centered=any(c["params"]["feature"].get("centred") for c in document["commands"]),
                         model_sha256=hashlib.sha256(outputs[0]).hexdigest()))
    pocket=work/"pocket"
    base=json.loads((pocket/"canonical.json").read_text())
    script=work/"mutated-timeline.json";target=work/"protected.native.json"
    before=(pocket/"archive.native.json").read_bytes();target.write_bytes(before)
    mutations={"suppression":{"suppressed":True},"repeated_id":{"id":1},"expression":{"width":"d1 + 1"},
               "nonintegral_corner":{"corner":[1.125,3,8]},"duplicate_parameter":{"param_names":["d1","d5","d6"]},
               "new_body":{"operation":"new_body"},"extra_geometry":{"taper":45},"embedded_unit":{"length":"4 in"}}
    refusals=[]
    for name,change in mutations.items():
        d=copy.deepcopy(base);d["commands"][1]["params"]["feature"].update(change)
        script.write_text(json.dumps(d))
        p=subprocess.run([cli,"import-solidcraft",script,target],capture_output=True,text=True)
        if not p.returncode or target.read_bytes()!=before:raise ValueError(name)
        refusals.append(dict(id=name,diagnostic=json.loads(p.stdout),previous_unchanged=True))
    direct=json.loads((pocket/"script.json").read_text());direct["commands"][1]["params"]["length"]=6
    edited=copy.deepcopy(base);edited["commands"][1]["params"]["feature"]["length"]="6"
    edit_outputs=[]
    for label,doc in (("archive",edited),("direct",direct)):
        path=work/(label+"-edited.json");path.write_text(json.dumps(doc))
        dest=work/(label+"-edited.native.json")
        p=subprocess.run([cli,"import-solidcraft",path,dest],capture_output=True,text=True,check=True)
        edit_outputs.append((dest.read_bytes(),json.loads(p.stdout)))
    if edit_outputs[0]!=edit_outputs[1] or edit_outputs[0][0]==before:
        raise ValueError("Edited length did not reproduce direct operation")
    return dict(schema="solidcraft_box_timeline_parity_v1",cases=rows,refusals=refusals,
                length_edit=dict(from_mm=4,to_mm=6,changed_model=True,identical_direct_replay=True,audit=edit_outputs[0][1]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("cli","host","work","output"):p.add_argument(k,type=Path)
    a=p.parse_args()
    a.output.write_text(json.dumps(run(a.cli.resolve(),a.host.resolve(),a.work),indent=2,allow_nan=False)+"\n")
