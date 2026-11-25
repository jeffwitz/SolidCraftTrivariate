"""Plot J4 candidate costs and the distances to the independent h=0.9 mm reference."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(costs, references, output):
    costs=json.loads(costs.read_text())
    references=json.loads(references.read_text())
    cases=[c for c in references["cases"] if c["reference_status"]=="qualified"]
    methods=[("local","Native, local refinement","#117567","o"),
             ("immersed","Immersed, h = 1.25 mm","#b45023","s"),
             ("tetra","C3D10, h = 1.5 mm","#315a96","D")]
    fig,(left,right)=plt.subplots(1,2,figsize=(10.4,4.8),layout="constrained")
    axes=[left,right]
    for ax in axes:
        ax.grid(axis="y",color="#d8dfe1",linewidth=.6)
        ax.set_axisbelow(True)
        ax.spines[["top","right"]].set_visible(False)
        ax.set_xticks(np.arange(len(cases)),[" / ".join(f"{v/16:g}" for v in c["radii_ticks"]) for c in cases])
        ax.set_xlabel("Radii of the two bores (mm)")
    for j,(method,label,color,marker) in enumerate(methods):
        shift=(j-1)*.19
        for i,case in enumerate(cases):
            trials=[v for row in costs["rows"] if row["radii_ticks"]==case["radii_ticks"]
                    for v in row["trials"] if v["method"]==method]
            times=sorted(v["end_to_end_with_probe_s"] for v in trials)
            assert len(times)==2
            x=i+shift
            left.vlines(x,*times,colors=color,lw=2)
            left.scatter([x]*2,times,c=color,marker=marker,s=32,zorder=3,label=label if i==0 else None)
            q=next(v for v in case["matched_trials"] if v["method"]==method)["field"]["energy"]
            right.errorbar(x,100*q["mean"],yerr=100*q["half_width_95"],fmt=marker,color=color,capsize=3)
    right.axhline(3,lw=1,ls="--",color="#545454")
    right.text(.98,3.04,"3% threshold",transform=right.get_yaxis_transform(),ha="right",fontsize=9)
    right.set_ylim(0,3.4)
    left.set_ylim(bottom=0)
    left.set_ylabel("Compute and probe time (s)")
    right.set_ylabel("Relative energy distance (%)")
    left.set_title("Two measurements per state and method",loc="left",fontsize=12)
    right.set_title("Independent C3D10 reference, h = 0.9 mm",loc="left",fontsize=12)
    fig.legend(*left.get_legend_handles_labels(),loc="outside lower center",ncol=3,frameon=False)
    fig.savefig(output,dpi=180,facecolor="white")
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("costs","references","output"):p.add_argument(k,type=Path)
    a=p.parse_args();plot(a.costs,a.references,a.output)
