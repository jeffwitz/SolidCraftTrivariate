"""Plot the actual seam layout and measured error per free vector coefficient."""

from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/"docs/results/curved-pocket-fields"
MODEL = ROOT/"docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"
DEST = ROOT/"docs/rtd/reference/curved-pocket-fields.svg"


def arc(block, side):
    ids = {"u0":[0,1,2],"u1":[6,7,8],"v0":[0,3,6],"v1":[2,5,8]}[side]
    net = np.array([block["controls"][2*k] for k in ids])
    pts = net[:,:2]/net[:,3:4]
    t = np.linspace(0,1,90).reshape(-1,1)
    return (1-t)**2*pts[0] + 2*t*(1-t)*pts[1] + t**2*pts[2]


def main():
    model = json.loads(MODEL.read_text())
    gluing = json.loads((DATA/"symbolic-gluing.json").read_text())
    states = json.loads((DATA/"seam-selection.json").read_text())
    blocks = {b["id"]:b for b in model["blocks"]}
    colors = {"C0":"#454b50","radial":"#147052","footprint":"#af4829","C1":"#474189"}
    with plt.rc_context({"svg.hashsalt":"solidcraft-curved-pocket","svg.fonttype":"none","font.size":9}):
        fig,axs = plt.subplots(1,3,figsize=(11,3.75),layout="constrained")
        geometry = axs[0]
        for face in model["design_faces"]:
            if face["block"].endswith("below") and face["feature"]=="outer/body":
                xy = arc(blocks[face["block"]],face["side"])
                geometry.plot(*xy.T,color=colors["C0"])
        for edge in gluing:
            name,side = edge["first"]
            if not name.endswith("below"):
                continue
            xy = arc(blocks[name],side)
            geometry.plot(*xy.T,color=colors["radial" if edge["as_g1"] else "footprint"],lw=2)
        geometry.set_aspect("equal")
        geometry.set(xlabel="x (mm)",ylabel="y (mm)",title="Seams below the pocket floor")
        geometry.legend(handles=[Line2D([],[],color=colors[k],label=s) for k,s in
                                 (("radial","Radial"),("footprint","Footprint"))],loc="lower center",
                        bbox_to_anchor=(0.5,-0.50),ncol=2,frameon=False)
        for p,ax in zip((3,4),axs[1:]):
            for mode,label in (("C0","C⁰"),("radial","Radial C¹"),("footprint","Footprint C¹"),("C1","All C¹")):
                rows=[r for r in states if r["selection"]==mode and r["p"]==p]
                ax.loglog([r["free_vector_dofs"] for r in rows],
                          [r["relative_energy_norm"] for r in rows],"o-",color=colors[mode],label=label)
            ax.set(title=f"Degree {p} · 4 and 8 intervals",xlabel="Free vector coefficients",ylabel="Relative energy-norm error")
            ax.grid(True,which="both",alpha=.22)
        axs[1].legend(fontsize=8)
        fig.savefig(DEST,metadata={"Date":None})
        plt.close(fig)
    DEST.write_text("\n".join(s.rstrip() for s in DEST.read_text().splitlines())+"\n")


if __name__=="__main__":
    main()
