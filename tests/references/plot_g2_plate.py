"""Plot the plate partition and signed wall curvature from stored volume controls."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import BPoly


def top_traces(block, height, exterior):
    shape = tuple(n + 1 for n in block["degrees"]) + (4,)
    net = np.asarray(block["controls"], dtype=float).reshape(shape)
    for axis in range(3):
        for end in (0, 1):
            cap = np.take(net, -1 if end else 0, axis=axis)
            if not np.allclose(cap[..., 2] / cap[..., 3], height, atol=1e-9, rtol=0):
                continue
            others = [d for d in range(3) if d != axis]
            for dim, other in enumerate(others):
                for side in (0, 1):
                    edge = np.take(cap, -1 if side else 0, axis=dim)
                    wall = (block["id"], "uvw"[other] + str(side)) in exterior
                    yield edge, wall


def sample(edge):
    curve = BPoly(edge[:, None, :], [0, 1])
    t = np.linspace(0, 1, 161)
    h, dh, ddh = (curve(t, nu=n) for n in range(3))
    w = h[:, -1:]
    p = h[:, :2] / w
    v = (dh[:, :2] - p * dh[:, -1:]) / w
    a = (ddh[:, :2] - p * ddh[:, -1:] - 2 * v * dh[:, -1:]) / w
    k = (v[:, 0] * a[:, 1] - v[:, 1] * a[:, 0]) / np.linalg.norm(v, axis=1)**3
    return p, k


def chain(curves):
    ordered = [curves.pop(0)]
    while curves:
        last = ordered[-1][0][-1]
        for i, (p, k) in enumerate(curves):
            if np.linalg.norm(last - p[0]) < 1e-8:
                ordered.append(curves.pop(i))
                break
            if np.linalg.norm(last - p[-1]) < 1e-8:
                curves.pop(i)
                ordered.append((p[::-1], -k[::-1]))
                break
        else:
            raise ValueError("Curved exterior traces do not form one cycle")
    if np.linalg.norm(ordered[0][0][0] - ordered[-1][0][-1]) > 1e-8:
        raise ValueError("Open curved wall")
    loop = np.concatenate([p for p, _ in ordered])
    area2 = np.sum(loop[:-1, 0]*loop[1:, 1] - loop[1:, 0]*loop[:-1, 1])
    return ordered if area2 > 0 else [(p[::-1], -k[::-1]) for p,k in reversed(ordered)]


def layer_section(block, z):
    """Evaluate a horizontal slice of a straight extrusion with matching end nets."""
    shape = tuple(n + 1 for n in block["degrees"]) + (4,)
    net = np.asarray(block["controls"], dtype=float).reshape(shape)
    for axis in range(3):
        a = net.take(0, axis=axis)
        b = net.take(-1, axis=axis)
        za = a[..., 2] / a[..., 3]
        zb = b[..., 2] / b[..., 3]
        lo, hi = za.flat[0], zb.flat[0]
        if not min(lo, hi) < z < max(lo, hi):
            continue
        if not np.allclose(za, lo) or not np.allclose(zb, hi):
            continue
        if not np.allclose(a[..., :2] / a[..., 3:], b[..., :2] / b[..., 3:]):
            raise ValueError("Cross-section figure requires straight axial extrusion")
        samples = []
        for dim, end, reverse in ((0, 0, False), (1, -1, False),
                                   (0, -1, True), (1, 0, True)):
            points, _ = sample(np.take(a, end, axis=dim))
            samples.append(points[::-1] if reverse else points)
        return np.concatenate(samples)
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("figure", type=Path)
    parser.add_argument("--slice-z", type=float, help="physical z coordinate inside the transverse cut")
    args = parser.parse_args()
    model = json.loads(args.model.read_text())
    blocks = model["blocks"]
    height = max(c[2] / c[3] for b in blocks for c in b["controls"])
    exterior = {(f["block"], f["side"]) for f in model["design_faces"]}
    sliced = args.slice_z is not None
    ratios = [1.5, 1.5, 1] if sliced else [1.6, 1]
    fig, axes = plt.subplots(len(ratios), 1, figsize=(8, 10 if sliced else 7),
                             layout="constrained", gridspec_kw={"height_ratios": ratios})
    plan, curvature = axes[0], axes[-1]
    if sliced:
        section = axes[1]
        for block in blocks:
            polygon = layer_section(block, args.slice_z)
            if polygon is not None:
                section.fill(*polygon.T, facecolor="#cce5e5", edgecolor="#57888b", lw=.65)
        section.set(aspect="equal", xlabel="x (mm)", ylabel="y (mm)",
                    title=f"Material cross-section at z = {args.slice_z:g} mm")

    traces = []
    for b in blocks:
        for edge, wall in top_traces(b, height, exterior):
            p, k = sample(edge)
            curved = wall and np.max(np.abs(k)) > 1e-8
            if curved:
                traces.append((p, k))
            plan.plot(*p.T, color="#176a72" if curved else "#243845" if wall else "#adc3cc",
                      lw=2.4 if wall else 0.8, zorder=2 if wall else 1)
    plan.set(aspect="equal", xlabel="x (mm)", ylabel="y (mm)",
             title="Top-face partition")
    crossings = model.get("sketch_plan", {}).get("crossings", [])
    cut_levels = sorted({tuple(c["y_mm"]) for c in crossings})
    for val in cut_levels:
        level = val[0] / val[1] if isinstance(val, (list, tuple)) else val
        plan.axhline(level, ls="--", lw=.8, color="#a8421f", alpha=.65)
    crossing_xy = np.array([[a / b for a, b in c["xy_mm"]] for c in crossings])
    labelled = set()
    s0 = 0
    for i, (p, k) in enumerate(chain(traces)):
        s = s0 + np.r_[0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
        curvature.plot(s, k, color="#176a72", lw=2)
        induced = len(crossing_xy) and np.min(np.linalg.norm(crossing_xy - p[0], axis=1)) < 1e-8
        label, color, marker = (("Sketch intersection", "#a8421f", "D") if induced
                                else ("Profile node", "#176a72", "o"))
        curvature.scatter(s[0], k[0], s=23, color=color, marker=marker, zorder=3,
                          label=label if label not in labelled else None)
        labelled.add(label)
        plan.scatter(*p[0], s=19, color=color, marker=marker, zorder=3)
        s0 = s[-1]

    curvature.set(xlabel="Approximate arc length (mm)", ylabel="Signed curvature (mm⁻¹)",
                  title="Signed curvature along one NURBS boundary period")
    curvature.grid(axis="y", alpha=0.25)
    curvature.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.26), ncol=2)
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure, dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
