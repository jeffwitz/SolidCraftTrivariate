"""Plot the source NURBS and its exact knot boundary for the campaign example."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from nurbs_sketch_oracle import curve


def draw(recipe, output):
    lo, hi = recipe["lower"], recipe["upper"]
    t = np.unique(lo["knots"])
    fig, ax = plt.subplots(figsize=(7.8, 4.3))
    for i, (a,b) in enumerate(zip(t,t[1:])):
        grid = np.linspace(a,b,161)
        lower = np.asarray([curve(lo,u) for u in grid])
        upper = np.asarray([curve(hi,u) for u in grid])
        ax.fill_between(lower[:,0],lower[:,1],upper[:,1],color=("#d2e7ee","#f7dec4")[i%2],
                        label=f"Travée {i+1}")
        for line in (lower,upper):
            ax.plot(*line.T,color="#23455c",lw=2)
        ax.plot([lower[0,0],upper[0,0]],[lower[0,1],upper[0,1]],color="#617888",lw=1)
    ax.plot([lower[-1,0],upper[-1,0]],[lower[-1,1],upper[-1,1]],color="#617888",lw=1)
    for spec in (lo,hi):
        points=np.asarray(spec["controls"])
        xy=points[:,:2]/points[:,2,None]
        ax.plot(*xy.T,"o--",ms=3.7,lw=.8,color="#596675",label="Polygones de contrôle" if spec is lo else None)
    ax.set(aspect="equal",xlabel="x (mm)",ylabel="y (mm)")
    ax.set_title("Section de l’extrusion à deux bords NURBS cubiques",loc="left",pad=14)
    ax.spines[["right","top"]].set_visible(False)
    ax.legend(loc="upper center",bbox_to_anchor=(.5,-.17),frameon=False,ncol=3)
    fig.tight_layout()
    fig.savefig(output,dpi=170,bbox_inches="tight")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("recipe",type=Path); p.add_argument("output",type=Path)
    args=p.parse_args()
    draw(json.loads(args.recipe.read_text()),args.output)
