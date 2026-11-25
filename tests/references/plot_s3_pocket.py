"""Plot S2/S3 cross-sections from their stored rational trivariate maps.

The figure is a sampled visualization, not a geometric certificate. Every
section patch is evaluated from the model's homogeneous Bezier control net.
The script requires an axis-separable Cartesian coordinate, as in the current
extruded S2/S3 block maps; it refuses a section it cannot derive this way.

Example:
    rtk python3 tests/references/plot_s3_pocket.py \
        --s2 docs/site/data/models/g2-cross-plate.model.json \
        --s3 <S3-model.json> \
        --output docs/figures/s3-face-pocket-sections.png
"""

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.path import Path as MplPath


COORDINATE_INDEX = {"x": 0, "y": 1, "z": 2}


def bernstein(degree, parameters):
    """Return Bernstein basis rows at the requested parameter values."""
    parameters = np.asarray(parameters, dtype=float)
    return np.stack(
        [
            math.comb(degree, i)
            * parameters**i
            * (1.0 - parameters) ** (degree - i)
            for i in range(degree + 1)
        ],
        axis=-1,
    )


def read_model(path):
    model = json.loads(path.read_text())
    if not isinstance(model.get("blocks"), list) or not model["blocks"]:
        raise ValueError(f"{path}: model has no blocks")
    for block in model["blocks"]:
        expected = int(np.prod(np.asarray(block["degrees"]) + 1))
        if len(block["controls"]) != expected:
            raise ValueError(f"{path}: malformed control net in {block.get('id')}")
    return model


def block_net(block):
    degrees = block["degrees"]
    shape = tuple(degree + 1 for degree in degrees) + (4,)
    net = np.asarray(block["controls"], dtype=float).reshape(shape)
    if np.any(net[..., 3] <= 0):
        raise ValueError(f"Non-positive rational weight in {block.get('id')}")
    return degrees, net


def independent_coordinate_axis(net, coordinate, parameter_axis, tolerance=2e-10):
    """Check that a physical coordinate depends on only one parameter axis."""
    weights = net[..., 3]
    physical = net[..., coordinate] / weights
    moved = np.moveaxis(physical, parameter_axis, 0)
    levels = moved[:, 0, 0]
    scale = max(1.0, float(np.max(np.abs(physical))))
    for i in range(len(levels)):
        if float(np.ptp(moved[i])) > tolerance * scale:
            return None

    moved_weights = np.moveaxis(weights, parameter_axis, 0)
    factors = moved_weights / moved_weights[0]
    for i in range(len(levels)):
        if not np.allclose(factors[i], factors[i, 0, 0], rtol=0, atol=2e-10):
            return None
    return levels


def parameter_at_coordinate(levels, target):
    """Invert a monotone rational Bezier coordinate by bisection."""
    delta = np.diff(levels)
    increasing = np.all(delta >= -1e-12) and np.any(delta > 1e-12)
    decreasing = np.all(delta <= 1e-12) and np.any(delta < -1e-12)
    if not (increasing or decreasing):
        return None
    low_value, high_value = min(levels[0], levels[-1]), max(levels[0], levels[-1])
    if target < low_value - 1e-9 or target > high_value + 1e-9:
        return None
    if abs(target - levels[0]) <= 1e-10:
        return 0.0
    if abs(target - levels[-1]) <= 1e-10:
        return 1.0

    lo, hi = 0.0, 1.0
    for _ in range(64):
        mid = (lo + hi) / 2
        basis = bernstein(len(levels) - 1, [mid])[0]
        value = float(basis @ levels)
        if (value < target) == increasing:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def section_patch(block, plane_axis, coordinate, samples):
    """Return a sampled boundary polygon for one block's planar section."""
    degrees, net = block_net(block)
    control_coordinate = net[..., plane_axis] / net[..., 3]
    if coordinate < float(np.min(control_coordinate)) - 1e-9 or coordinate > float(np.max(control_coordinate)) + 1e-9:
        return None
    candidates = []
    separable_axis_found = False
    for parameter_axis in range(3):
        levels = independent_coordinate_axis(net, plane_axis, parameter_axis)
        if levels is None:
            continue
        separable_axis_found = True
        parameter = parameter_at_coordinate(levels, coordinate)
        if parameter is not None:
            candidates.append((parameter_axis, parameter))
    if not candidates:
        if separable_axis_found:
            return None
        raise ValueError(
            f"Block {block.get('id')} intersects the section control hull but its "
            "coordinate is not separable along one Bezier parameter axis"
        )
    parameter_axis, parameter = candidates[0]
    basis = bernstein(degrees[parameter_axis], [parameter])[0]
    restricted = np.tensordot(basis, np.moveaxis(net, parameter_axis, 0), axes=(0, 0))

    remaining_degrees = [degrees[axis] for axis in range(3) if axis != parameter_axis]
    grid = np.linspace(0.0, 1.0, samples)
    basis_a = bernstein(remaining_degrees[0], grid)
    basis_b = bernstein(remaining_degrees[1], grid)
    homogeneous = np.einsum("ai,bj,ijc->abc", basis_a, basis_b, restricted)
    cartesian = homogeneous[..., :3] / homogeneous[..., 3:4]
    plotted_axes = [axis for axis in range(3) if axis != plane_axis]
    projected = cartesian[..., plotted_axes]

    boundary = np.concatenate(
        (
            projected[:, 0],
            projected[-1, 1:],
            projected[-2::-1, -1],
            projected[0, -2:0:-1],
        ),
        axis=0,
    )
    if len(boundary) < 4 or abs(polygon_area(boundary)) < 1e-14:
        return None
    if polygon_area(boundary) < 0:
        boundary = boundary[::-1]
    return boundary


def polygon_area(points):
    return 0.5 * float(
        np.sum(
            points[:, 0] * np.roll(points[:, 1], -1)
            - points[:, 1] * np.roll(points[:, 0], -1)
        )
    )


def collect_sections(model, plane_axis, coordinate, samples):
    polygons = []
    for block in model["blocks"]:
        polygon = section_patch(block, plane_axis, coordinate, samples)
        if polygon is not None:
            polygons.append(polygon)
    if not polygons:
        raise ValueError(
            f"No section patches found at {('x', 'y', 'z')[plane_axis]}={coordinate:g}"
        )
    return polygons


def sample_model_bounds(model, axis):
    recipe = model.get("source_recipe", {})
    bounds = recipe.get("bounds_xyz_mm") if isinstance(recipe, dict) else None
    if bounds and len(bounds) == 3:
        return [float(v) for v in bounds[axis]]
    values = []
    for block in model["blocks"]:
        degrees, net = block_net(block)
        bases = [bernstein(degree, np.linspace(0, 1, 9)) for degree in degrees]
        points = np.einsum("ai,bj,ck,ijkd->abcd", *bases, net)
        values.extend((points[..., axis] / points[..., 3]).ravel())
    return [float(min(values)), float(max(values))]


def rational_value(value):
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return float(value[0]) / float(value[1])
    return float(value)


def pocket_prism(model):
    """Return the source rectangle and extrusion vertices when available."""
    operation = model.get("source_recipe", {}).get("operation", {})
    primitive = operation.get("loop", {}).get("primitive", {})
    if primitive.get("kind") != "rectangle":
        return None
    local_min = [rational_value(v) for v in primitive["min"]]
    local_max = [rational_value(v) for v in primitive["max"]]
    frame = operation["frame"]
    origin = np.asarray(frame["origin_mm"], dtype=float)
    u = np.asarray(frame["u"], dtype=float)
    v = np.asarray(frame["v"], dtype=float)
    corners = np.asarray(
        [
            origin + u * x + v * y
            for x, y in (
                (local_min[0], local_min[1]),
                (local_max[0], local_min[1]),
                (local_max[0], local_max[1]),
                (local_min[0], local_max[1]),
            )
        ]
    )
    direction = np.asarray(operation["direction"], dtype=float)
    end = corners + direction * float(operation["depth_mm"])
    vertices = np.concatenate((corners, end), axis=0)
    return corners, vertices


def polygon_mask(xs, ys, polygon):
    mask = np.zeros((len(ys), len(xs)), dtype=bool)
    if polygon is None or len(polygon) < 3:
        return mask
    path = MplPath(np.asarray(polygon), closed=True)
    lo_x, hi_x = np.min(polygon[:, 0]), np.max(polygon[:, 0])
    lo_y, hi_y = np.min(polygon[:, 1]), np.max(polygon[:, 1])
    i0 = max(0, int(np.searchsorted(xs, lo_x, side="left")) - 1)
    i1 = min(len(xs), int(np.searchsorted(xs, hi_x, side="right")) + 1)
    j0 = max(0, int(np.searchsorted(ys, lo_y, side="left")) - 1)
    j1 = min(len(ys), int(np.searchsorted(ys, hi_y, side="right")) + 1)
    if i0 < i1 and j0 < j1:
        xx, yy = np.meshgrid(xs[i0:i1], ys[j0:j1])
        mask[j0:j1, i0:i1] = path.contains_points(
            np.column_stack((xx.ravel(), yy.ravel()))
        ).reshape(yy.shape)
    return mask


def prism_section_polygon(vertices, axis, coordinate, plotted_axes):
    """Intersect a rectangular extrusion with a plane using its twelve edges."""
    if vertices is None:
        return None
    edges = [(i, (i + 1) % 4) for i in range(4)]
    edges += [(i + 4, (i + 1) % 4 + 4) for i in range(4)]
    edges += [(i, i + 4) for i in range(4)]
    points = []
    for i, j in edges:
        a, b = vertices[i], vertices[j]
        da, db = a[axis] - coordinate, b[axis] - coordinate
        if abs(da) <= 1e-10:
            points.append(a)
        if abs(db) <= 1e-10:
            points.append(b)
        if da * db < 0:
            points.append(a + (b - a) * (da / (da - db)))
    if len(points) < 3:
        return None
    unique = np.unique(np.round(np.asarray(points), decimals=10), axis=0)
    if len(unique) < 3:
        return None
    projected = unique[:, plotted_axes]
    center = np.mean(projected, axis=0)
    angles = np.arctan2(projected[:, 1] - center[1], projected[:, 0] - center[0])
    return projected[np.argsort(angles)]


def section_bounds(polygons):
    all_points = np.concatenate(polygons)
    return (
        float(np.min(all_points[:, 0])),
        float(np.max(all_points[:, 0])),
        float(np.min(all_points[:, 1])),
        float(np.max(all_points[:, 1])),
    )


def rasterize(polygons, bounds, resolution):
    x0, x1, y0, y1 = bounds
    xs = np.linspace(x0, x1, resolution)
    ys = np.linspace(y0, y1, resolution)
    mask = np.zeros((resolution, resolution), dtype=bool)
    for polygon in polygons:
        path = MplPath(polygon, closed=True)
        lo_x, hi_x = np.min(polygon[:, 0]), np.max(polygon[:, 0])
        lo_y, hi_y = np.min(polygon[:, 1]), np.max(polygon[:, 1])
        i0 = max(0, int(np.searchsorted(xs, lo_x, side="left")) - 1)
        i1 = min(resolution, int(np.searchsorted(xs, hi_x, side="right")) + 1)
        j0 = max(0, int(np.searchsorted(ys, lo_y, side="left")) - 1)
        j1 = min(resolution, int(np.searchsorted(ys, hi_y, side="right")) + 1)
        if i0 >= i1 or j0 >= j1:
            continue
        xx, yy = np.meshgrid(xs[i0:i1], ys[j0:j1])
        points = np.column_stack((xx.ravel(), yy.ravel()))
        mask[j0:j1, i0:i1] |= path.contains_points(points).reshape(yy.shape)
    return xs, ys, mask


def zoom_to_removed(ax, xs, ys, removed):
    rows, columns = np.nonzero(removed)
    if not len(rows):
        return False
    x0, x1 = xs[columns.min()], xs[columns.max()]
    y0, y1 = ys[rows.min()], ys[rows.max()]
    dx, dy = max(x1 - x0, 1e-3), max(y1 - y0, 1e-3)
    mx, my = 0.18 * dx, 0.18 * dy
    ax.set_xlim(x0 - mx, x1 + mx)
    ax.set_ylim(y0 - my, y1 + my)
    return True


def draw_material(ax, xs, ys, before, after, *, show_before=False, removal_clip=None):
    removed = before & ~after
    if removal_clip is not None:
        removed &= removal_clip
    classes = np.zeros(before.shape, dtype=np.uint8)
    classes[after] = 1
    classes[removed] = 2
    if show_before:
        classes = np.where(before, 1, 0).astype(np.uint8)
        cmap = ListedColormap(["#ffffff", "#a9c9d2"])
    else:
        cmap = ListedColormap(["#ffffff", "#2f7f78", "#dc8a4c"])
    ax.imshow(
        classes,
        origin="lower",
        interpolation="nearest",
        extent=(xs[0], xs[-1], ys[0], ys[-1]),
        cmap=cmap,
        vmin=0,
        vmax=cmap.N - 1,
        aspect="equal",
    )
    return removed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s2", type=Path, required=True, help="Accepted S2 model JSON")
    parser.add_argument("--s3", type=Path, required=True, help="Accepted S3 model JSON")
    parser.add_argument("--output", type=Path, required=True, help="PNG output path")
    parser.add_argument("--z", type=float, default=3.0, help="XY section height in mm")
    parser.add_argument("--y", type=float, help="XZ section ordinate in mm; defaults to domain midpoint")
    parser.add_argument("--samples", type=int, default=49, help="Bezier section boundary samples per direction")
    parser.add_argument("--raster-resolution", type=int, default=1400, help="Raster pixels per axis")
    parser.add_argument("--zoom-margin", type=float, default=4.0, help="XY sketch zoom margin in mm")
    args = parser.parse_args()
    if args.samples < 5 or args.raster_resolution < 200:
        parser.error("--samples must be at least 5 and --raster-resolution at least 200")

    s2, s3 = read_model(args.s2), read_model(args.s3)
    prism = pocket_prism(s3)
    y_bounds = sample_model_bounds(s2, 1)
    if args.y is not None:
        y_cut = args.y
    elif prism is not None:
        y_cut = float(np.mean(prism[1][:, 1]))
    else:
        y_cut = sum(y_bounds) / 2
    specs = (
        ("z", args.z, [0, 1], "x (mm)", "y (mm)"),
        ("y", y_cut, [0, 2], "x (mm)", "z (mm)"),
    )
    sections = []
    for axis_name, coordinate, _, xlabel, ylabel in specs:
        axis = COORDINATE_INDEX[axis_name]
        before = collect_sections(s2, axis, coordinate, args.samples)
        after = collect_sections(s3, axis, coordinate, args.samples)
        bounds = section_bounds(before + after)
        xs, ys, mask2 = rasterize(before, bounds, args.raster_resolution)
        _, _, mask3 = rasterize(after, bounds, args.raster_resolution)
        sections.append((axis_name, coordinate, xlabel, ylabel, xs, ys, mask2, mask3))

    fig = plt.figure(figsize=(14.8, 7.0), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=(1.15, 0.8))
    axes = (
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[1, :]),
    )
    xy, xz = sections
    draw_material(axes[0], xy[4], xy[5], xy[6], xy[7], show_before=True)
    axes[0].set_title(f"Before — S2, z = {args.z:g} mm")
    axes[0].set_xlabel(xy[2])
    axes[0].set_ylabel(xy[3])

    removed_xy = xy[6] & ~xy[7]
    if prism is not None:
        sketch_xy = prism[0][:, [0, 1]]
        removed_xy &= polygon_mask(xy[4], xy[5], sketch_xy)
    draw_material(axes[1], xy[4], xy[5], xy[6], xy[7], removal_clip=removed_xy)
    axes[1].set_title("After — removed material shown in orange")
    axes[1].set_xlabel(xy[2])
    axes[1].set_ylabel(xy[3])
    if prism is not None:
        sketch_loop = np.vstack((sketch_xy, sketch_xy[0]))
        for ax in axes[:2]:
            ax.plot(*sketch_loop.T, color="#a8421f", ls="--", lw=1.2, zorder=5)
            ax.set_xlim(np.min(sketch_xy[:, 0]) - args.zoom_margin,
                        np.max(sketch_xy[:, 0]) + args.zoom_margin)
            ax.set_ylim(np.min(sketch_xy[:, 1]) - args.zoom_margin,
                        np.max(sketch_xy[:, 1]) + args.zoom_margin)
    elif not zoom_to_removed(axes[1], xy[4], xy[5], removed_xy):
        axes[1].text(0.5, 0.5, "No material removal detected at this coordinate", transform=axes[1].transAxes,
                     ha="center", va="center")

    removed_xz = xz[6] & ~xz[7]
    prism_xz = prism_section_polygon(
        prism[1] if prism is not None else None, 1, y_cut, [0, 2]
    )
    if prism_xz is not None:
        removed_xz &= polygon_mask(xz[4], xz[5], prism_xz)
    draw_material(axes[2], xz[4], xz[5], xz[6], xz[7], removal_clip=removed_xz)
    axes[2].set_title(f"XZ cross-section — y = {y_cut:g} mm")
    axes[2].set_xlabel(xz[2])
    axes[2].set_ylabel(xz[3])
    zoom_to_removed(axes[2], xz[4], xz[5], removed_xz)
    if not removed_xz.any() and prism_xz is not None:
        axes[2].set_xlim(np.min(prism_xz[:, 0]) - args.zoom_margin,
                         np.max(prism_xz[:, 0]) + args.zoom_margin)
        axes[2].set_ylim(np.min(prism_xz[:, 1]) - args.zoom_margin,
                         np.max(prism_xz[:, 1]) + args.zoom_margin)
    fig.legend(
        handles=[
            Patch(facecolor="#2f7f78", label="Remaining material"),
            Patch(facecolor="#dc8a4c", label="Removed material"),
        ],
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.035),
        ncol=2,
    )
    fig.suptitle("Comparison of S2 and S3 trivariate maps", fontsize=15)
    fig.text(
        0.5,
        0.005,
        f"Sections evaluated on rational maps; boundaries sampled at {args.samples} points per direction. "
        f"Display mask: {args.raster_resolution} × {args.raster_resolution}.",
        ha="center",
        fontsize=8,
        color="#52606b",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.output}")
    print(f"S2 blocks: {len(s2['blocks'])}; S3 blocks: {len(s3['blocks'])}")
    print(f"Sections: z={args.z:g} mm; y={y_cut:g} mm")
    print(f"XY removed pixels: {int(removed_xy.sum())}; XZ removed pixels: {int(removed_xz.sum())}")


if __name__ == "__main__":
    main()
