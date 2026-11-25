"""Draw sampled exterior traces of two models exported by the native builder."""
import argparse
import itertools
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def basis(degree, t):
    return np.stack([math.comb(degree, i) * t**i * (1-t)**(degree-i)
                     for i in range(degree+1)], axis=-1)


def sample_face(block, side):
    degrees = block["degrees"]
    grid = np.linspace(0, 1, 10)
    u, v = np.meshgrid(grid, grid, indexing="ij")
    free = iter((u, v))
    axis = "uvw".index(side[0])
    coords = [np.full_like(u, int(side[1])) if i == axis else next(free) for i in range(3)]
    weights = [basis(n, t) for n, t in zip(degrees, coords)]
    net = np.asarray(block["controls"]).reshape(*(n+1 for n in degrees), 4)
    homogeneous = np.einsum("abi,abj,abk,ijkh->abh", *weights, net)
    return homogeneous[..., :3] / homogeneous[..., 3:]


def panel(ax, model, walls_only):
    blocks = {block["id"]: block for block in model["blocks"]}
    bore_names = sorted(face["feature"] for face in model["design_faces"] if face["feature"].startswith("bore/"))
    first_bore = bore_names[0]
    for face in model["design_faces"]:
        name = face["feature"]
        is_bore = name.startswith("bore/")
        if walls_only and not is_bore:
            continue
        x = sample_face(blocks[face["block"]], face["side"])
        color = ("#c46430" if name == first_bore else "#15798a") if is_bore else "#93b4c4"
        ax.plot_surface(*x.transpose(2,0,1), color=color, linewidth=0, antialiased=True,
                        alpha=0.96 if is_bore else 0.78, shade=True)
        for edge in (x[0], x[-1], x[:,0], x[:,-1]):
            ax.plot(*edge.T, color="#263f4d", linewidth=0.45, alpha=0.85)
    vertices = np.asarray(model["outer"]["vertices"])
    bounds = [list(pair) for pair in zip(vertices.min(axis=0), vertices.max(axis=0))]
    bounds.append([0, model["height_mm"]])
    if walls_only:
        for p, q in itertools.combinations(itertools.product(*bounds), 2):
            if sum(a != b for a, b in zip(p, q)) == 1:
                ax.plot(*zip(p,q), color="#74919d", linewidth=0.8, alpha=0.65)
    ax.set_box_aspect([hi-lo for lo,hi in bounds])
    for setter, pair in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), bounds):
        setter(*pair)
    ax.view_init(24, 52 if not walls_only else -57)
    ax.set_axis_off()
    return len(blocks), len(model["topology"]["interfaces"])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("bracket", "cross", "output"):
        p.add_argument(name, type=Path)
    args = p.parse_args()
    fig = plt.figure(figsize=(11, 5.4), facecolor="white")
    titles = ["Sketch-generated bracket", "Two intersecting bores"]
    for i, file in enumerate((args.bracket, args.cross)):
        ax = fig.add_subplot(1, 2, i+1, projection="3d")
        n, links = panel(ax, json.loads(file.read_text()), bool(i))
        ax.set_title(titles[i], fontsize=16, loc="center", weight="bold", color="#213d50")
        ax.text2D(.5,.03,f"{n} blocks  ·  {links} interfaces",transform=ax.transAxes,
                  fontsize=12,ha="center",color="#213d50")
    fig.text(.5,.045,"Left: outer boundary and chart edges. Right: the two bore walls inside the box.",
             ha="center",fontsize=11, color="#314b59")
    fig.subplots_adjust(left=.015,right=.985,top=.91,bottom=.105,wspace=.02)
    fig.savefig(args.output, dpi=180)
