"""Plot intersecting bore walls and the exact octant construction section."""
import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def bernstein(degree, values):
    return np.stack([math.comb(degree, i) * values**i * (1-values)**(degree-i)
                     for i in range(degree+1)], axis=-1)


def face_points(block, side, samples=101):
    degrees = block["degrees"]
    fixed_axis = "uvw".index(side[0])
    coords = []
    free = iter(np.meshgrid(np.linspace(0, 1, samples), np.linspace(0, 1, samples), indexing="ij"))
    for axis in range(3):
        coords.append(np.full((samples, samples), int(side[1]), dtype=float) if axis == fixed_axis else next(free))
    weights = [bernstein(degree, value) for degree, value in zip(degrees, coords)]
    net = np.asarray(block["controls"], dtype=float).reshape(*(degree+1 for degree in degrees), 4)
    homogeneous = np.einsum("abi,abj,abk,ijkh->abh", *weights, net)
    return homogeneous[..., :3] / homogeneous[..., 3:]


def draw(model_path: Path, output: Path, clip_half_length=7.0):
    model = json.loads(model_path.read_text())
    blocks = {block["id"]: block for block in model["blocks"]}
    center = np.asarray(model["partition_plan"]["axis_intersection_mm"], dtype=float)
    radius = float(model["partition_plan"]["radius_mm"])
    guard = float(model["partition_plan"]["guard_mm"])
    half_extents = (np.asarray(model["outer"]["vertices"]).max(axis=0)
                    - np.asarray(model["outer"]["vertices"]).min(axis=0)) / 2
    half_height = float(np.ptp(np.asarray(model["outer"]["vertices"])[:, 1])) / 2
    bore_colors = {bore["id"]: color for bore, color in zip(model["bores"], ("#247e99", "#c8692e"))}
    fig = plt.figure(figsize=(12, 6.2), facecolor="white")
    left = fig.add_subplot(1, 2, 1, projection="3d")
    right = fig.add_subplot(1, 2, 2)
    for face in model["design_faces"]:
        feature = face["feature"]
        if not feature.startswith("bore/") or not feature.endswith("/wall"):
            continue
        bore_id = feature.split("/")[1]
        xyz = face_points(blocks[face["block"]], face["side"]) - center
        bore = next(item for item in model["bores"] if item["id"] == bore_id)
        axis = "xyz".index(bore["axis"])
        keep = np.abs(xyz[..., axis]) <= clip_half_length
        xyz[~keep] = np.nan
        left.plot_surface(*xyz.transpose(2, 0, 1), color=bore_colors[bore_id],
                          linewidth=0, antialiased=True, alpha=.78, shade=True)
    left.scatter([0], [0], [0], color="#18252d", s=16, depthshade=False)
    left.set(xlim=(-clip_half_length, clip_half_length), ylim=(-clip_half_length, clip_half_length),
             zlim=(-clip_half_length, clip_half_length), xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)")
    left.set_box_aspect((1, 1, 1))
    left.set_title("Bore walls near their intersection", pad=12)
    left.view_init(24, -57)

    q = np.linspace(0, 1, 240)
    circle_x = radius * (1-q*q) / (1+q*q)
    circle_y = radius * (2*q) / (1+q*q)
    outer = np.array([[0, 0], [guard, 0], [guard, half_height], [0, half_height]])
    material = np.vstack(([radius, 0], outer[1:], [0, half_height], [0, radius],
                          np.column_stack((circle_x, circle_y))[::-1]))
    right.fill(material[:, 0], material[:, 1], color="#d9ecef", zorder=0)
    right.plot([0, guard, guard, 0, 0], [0, 0, half_height, half_height, 0], color="#244655", lw=1.5)
    right.plot(circle_x, circle_y, color="#173742", lw=2.2)
    separator_start = np.array([3*radius/5, 4*radius/5])
    separator_end = np.array([guard, half_height])
    right.plot(*zip(separator_start, separator_end), color="#d76d33", lw=2.4)
    right.annotate("Separator", xy=(2.45, 8.2), xytext=(1.78, 7.3), color="#c85c27",
                   arrowprops={"arrowstyle": "->", "color": "#c85c27"})
    right.text(.38, 12, "Material", color="#183c4a")
    right.text(.38, .72, "Void", color="#183c4a")
    right.set(xlim=(0, guard+0.1), ylim=(0, half_height+0.5), xlabel="min(|X|, |Z|) (mm)", ylabel="|Y| (mm)")
    right.set_title("Construction section in one octant", pad=12)
    right.set_aspect("auto")
    fig.suptitle(f"Two intersecting bores, radius {radius:g} mm", fontsize=16, weight="bold")
    fig.subplots_adjust(left=.04, right=.98, top=.88, bottom=.15, wspace=.25)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=150, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    draw(args.model, args.output)
