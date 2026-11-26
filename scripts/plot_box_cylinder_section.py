"""Plot the physical images of the native parameter lines at half height."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def section(block, u, v):
    if block["degrees"] != [2, 1, 1]:
        raise ValueError("Expected quadratic-linear-linear maps")
    net = np.asarray(block["controls"], dtype=float).reshape(3, 2, 2, 4)
    u, v = np.broadcast_arrays(u, v)
    bu = np.stack(((1-u)**2, 2*u*(1-u), u**2), axis=-1)
    bv = np.stack((1-v, v), axis=-1)
    h = np.einsum("...i,...j,ijkc,k->...c", bu, bv, net, [0.5, 0.5])
    xyz = h[..., :3] / h[..., 3:4]
    if not np.allclose(xyz[..., 2], 6):
        raise ValueError("The w=1/2 section must lie on z=6 mm")
    return xyz[..., :2].T


def main():
    directory = Path(__file__).resolve().parents[1] / "docs/site/box-cylinder"
    colors = ("#14776e", "#b3562c", "#243536")
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 5.0))
    t = np.linspace(0, 1, 161)
    for ax, name, title in zip(axes, ["result", "classic"], ["Boolean difference", "Classical sketch partition"]):
        for block in json.loads((directory / f"{name}.model.json").read_bytes())["blocks"]:
            for a in np.arange(1, 8) / 8:
                ax.plot(*section(block, a, t), color=colors[0], lw=.85)
                ax.plot(*section(block, t, a), color=colors[1], lw=.85)
            for end in [0, 1]:
                ax.plot(*section(block, end, t), color=colors[2], lw=1.15)
                ax.plot(*section(block, t, end), color=colors[2], lw=1.15)
        ax.set(aspect="equal", title=title, xlabel="x (mm)", xticks=[0,10,20,30,40], yticks=[0,10,20,30,40])
        ax.spines[["right", "top"]].set_visible(False)
    axes[0].set_ylabel("y (mm)")
    fig.subplots_adjust(left=.08, right=.98, top=.9, bottom=.21, wspace=.2)
    fig.legend([Line2D([], [], color=c) for c in colors],
               ["Constant u", "Constant v", "Patch boundaries"], loc="lower center", bbox_to_anchor=(.5, .01), ncols=3, frameon=False)
    fig.savefig(directory / "parameter-lines.png", dpi=145, facecolor="white", bbox_inches="tight", pad_inches=.12)
    plt.close(fig)


if __name__ == "__main__":
    main()
