"""Measure all available h2 local fields against the same two C3D10 references."""
import argparse
from pathlib import Path

from adaptive_field_distance import AdaptiveField
from bracket_cost_fields import compare, read, save
from bracket_field_distance import TetraField


def run(fine_inputs, tiers, summary, radius, output):
    rx, ry = radius
    tag = f"x{rx}-y{ry}"
    source = Path(__file__).resolve().parents[2]/"build/bracket-cost"/tag/"tetra"
    prepared = next(r for r in read(fine_inputs/"preparation.json") if r["radii_ticks"]==radius)
    if not prepared.get("fine_qualified"):raise ValueError("h1 reference not qualified")
    inputs, path = fine_inputs/tag, tiers/tag
    states = [r for r in read(summary)["states"] if r["radii_ticks"]==radius]
    if len(states)!=4:raise ValueError("incomplete radius ladder")
    blocks = read(inputs/"native.package.json")["analysis_space"]["blocks"]
    fields, scales = {},{}
    accepted = []
    for row in states:
        if row["status"]!="completed":continue
        step=row["step"]
        key=f"native_s{step}"
        fields[key]=AdaptiveField(path/f"s{step}.bin")
        scales[key]=row["q6_q8"][1]["strain_energy_n_mm"]
        if row["qualified"]:accepted.append(key)
    for key, directory in (("tetra_h125",source/"h1.25-qNone"),("tetra_h1",inputs/"h1")):
        fields[key]=TetraField(directory)
        d=read(directory/"result.json")
        scales[key]=d.get("strain_energy_n_mm",-500.*d["mean_tip_y_mm"])
    comparison=compare(blocks,fields,scales,11)
    final=comparison["levels"][-1]
    thresholds = dict(energy=.03, l2=.005)
    qualified=[]
    verdicts=[]
    for row in states:
        key=f"native_s{row['step']}"
        status=dict(step=row["step"],completed=row["status"]=="completed",qoi=key in accepted,qualified=False)
        if key in fields:
            distances=final["distances"][key+"__tetra_h1"]
            status["distances"]=distances
            status["reference_difference"]=final["distances"]["tetra_h125__tetra_h1"]
            status["qualified"]=bool(status["qoi"] and comparison["all_points_located_once"] and
                comparison["recovered_energies"] and all(distances[k]["stable"] and distances[k]["mean"]<=v for k,v in thresholds.items()))
            if status["qualified"]:qualified.append(row["step"])
        verdicts.append(status)
    save(output,dict(radii_ticks=radius,thresholds=thresholds,steps=verdicts,
                     earliest_qualified=min(qualified,default=None),comparison=comparison))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs",type=Path);p.add_argument("tiers",type=Path);p.add_argument("summary",type=Path);p.add_argument("output",type=Path)
    p.add_argument("rx",type=int);p.add_argument("ry",type=int)
    a=p.parse_args()
    run(a.inputs,a.tiers,a.summary,[a.rx,a.ry],a.output)
