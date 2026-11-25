"""Draw section boundaries from the six models written by primitive_sweep_cli.py."""
import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def sample_edge(controls):
    t = np.linspace(0, 1, 101)
    p = len(controls) - 1
    basis = np.column_stack([math.comb(p, i) * (1 - t)**(p-i) * t**i for i in range(p+1)])
    values = basis @ controls
    return values[:, :2] / values[:, 3, None]


def section_boundary(block):
    du, dv, dw = block["degrees"]
    net = np.asarray(block["controls"]).reshape(du+1, dv+1, dw+1, 4)[:, :, 0, :]
    edges = (net[:, 0], net[-1, :], net[::-1, -1], net[0, ::-1])
    return np.concatenate([sample_edge(edge) for edge in edges])


def draw(models, output):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, (name, title) in zip(axes, (("rectangle", "Rectangle"), ("circle", "Disk"), ("ellipse", "Filled ellipse"))):
        def blocks(operation):
            return json.loads((models / f"{name}-{operation}.model.json").read_text())["blocks"]

        extrusion = blocks("extrude")
        for block in extrusion:
            points = section_boundary(block)
            ax.fill(*points.T, color="#d7eeea", edgecolor="#26736f", linewidth=1.6)
        ax.set_aspect("equal", adjustable="datalim")
        ax.margins(0.14)
        ax.set_axis_off()
        ax.set_title(title)
        n = len(extrusion)
        ax.text(0.5, -0.04, f"Extrusion · {n} block{'s' if n != 1 else ''}\n"
                f"Off-axis revolution · {len(blocks('revolve'))} blocks",
                transform=ax.transAxes, ha="center", va="top", fontsize=10)

    fig.subplots_adjust(left=.03, right=.97, top=.86, bottom=.22, wspace=.08)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=132, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("models", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    draw(args.models, args.output)
