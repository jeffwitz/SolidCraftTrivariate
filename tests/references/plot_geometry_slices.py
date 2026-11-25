"""Plot one-sided midheight chart conditioning from the qualified VTK exports."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import LogNorm
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def patches(path):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    nodes = vtk_to_numpy(grid.GetPoints().GetData()).reshape(-1, 7, 7, 7, 3)
    kappa = vtk_to_numpy(grid.GetPointData().GetArray("condition_number")).reshape(-1, 7, 7, 7)
    polygons, colors, outlines = [], [], []
    for xyz, values in zip(nodes[:, :, :, 3, :2], kappa[:, :, :, 3]):
        corners = [xyz[:-1, :-1], xyz[1:, :-1], xyz[1:, 1:], xyz[:-1, 1:]]
        polygons.extend(np.stack(corners, axis=-2).reshape(-1, 4, 2))
        colors.extend(((values[:-1, :-1] + values[1:, :-1] +
                        values[1:, 1:] + values[:-1, 1:]) / 4).ravel())
        outlines.extend((xyz[0], xyz[-1], xyz[:, 0], xyz[:, -1]))
    return polygons, colors, outlines


def plot(folder, output):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5), layout="constrained")
    norm = LogNorm(1, 10000)
    cases = [("two_tunnels", "Two polygonal tunnels", "z = 2.5 mm"),
             ("inclined_bores", "Inclined circular bore", "z = 0.5 mm")]
    for ax, (name, title, plane) in zip(axes, cases):
        poly, color, edges = patches(folder / f"{name}-6.vtu")
        field = PolyCollection(poly, array=np.array(color), cmap="viridis", norm=norm,
                               edgecolors="face", linewidths=.2)
        ax.add_collection(field)
        ax.add_collection(LineCollection(edges, colors="#1b2033", linewidths=.42))
        ax.autoscale()
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(f"{title}\n{plane}", fontsize=12)
    fig.colorbar(field, ax=axes, location="bottom", shrink=.65, aspect=45,
                 label="Jacobian condition number κ₂(J) · mean of the four vertices")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vtu_folder", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    plot(args.vtu_folder, args.output)
