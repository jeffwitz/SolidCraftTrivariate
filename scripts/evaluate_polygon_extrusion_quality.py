#!/usr/bin/env python3
"""Independent corner-sampled quality metrics for trilinear polygonal extrusions."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys


def _norm(a):
    return math.sqrt(sum(x * x for x in a))


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _singular_values_2(a, b):
    """Singular values of the 3x2 matrix with columns a and b."""
    aa, ab, bb = _dot(a, a), _dot(a, b), _dot(b, b)
    disc = math.sqrt(max(0.0, (aa - bb) ** 2 + 4.0 * ab * ab))
    return math.sqrt(max(0.0, (aa + bb + disc) / 2.0)), math.sqrt(max(0.0, (aa + bb - disc) / 2.0))


def _singular_values_3(columns):
    """Singular values of a 3x3 Jacobian using Jacobi eigenvalue sweeps."""
    gram = [[_dot(columns[i], columns[j]) for j in range(3)] for i in range(3)]
    for _ in range(32):
        p, q = max(((0, 1), (0, 2), (1, 2)), key=lambda ij: abs(gram[ij[0]][ij[1]]))
        if abs(gram[p][q]) <= 1e-15 * max(1.0, max(abs(gram[i][i]) for i in range(3))):
            break
        tau = (gram[q][q] - gram[p][p]) / (2.0 * gram[p][q])
        t = math.copysign(1.0, tau) / (abs(tau) + math.sqrt(1.0 + tau * tau))
        c = 1.0 / math.sqrt(1.0 + t * t)
        s = t * c
        app, aqq, apq = gram[p][p], gram[q][q], gram[p][q]
        gram[p][p] = app - t * apq
        gram[q][q] = aqq + t * apq
        gram[p][q] = gram[q][p] = 0.0
        for k in range(3):
            if k not in (p, q):
                gkp, gkq = gram[k][p], gram[k][q]
                gram[k][p] = gram[p][k] = c * gkp - s * gkq
                gram[k][q] = gram[q][k] = s * gkp + c * gkq
    return tuple(sorted((math.sqrt(max(0.0, gram[i][i])) for i in range(3)), reverse=True))


def block_metrics(block):
    controls = block["controls"]
    if len(controls) != 8:
        raise ValueError(f"{block.get('id', '<block>')}: expected 8 trilinear controls")
    if block.get("degrees") != [1, 1, 1]:
        raise ValueError(f"{block.get('id', '<block>')}: only degree [1, 1, 1] is supported")
    weights = {c[3] for c in controls}
    if len(weights) != 1 or next(iter(weights)) <= 0:
        raise ValueError(f"{block.get('id', '<block>')}: expected one positive common control weight")
    # Model controls use homogeneous (wx, wy, wz, w), in u-v-w bit order.
    xyz = [tuple(float(c[k]) / float(c[3]) for k in range(3)) for c in controls]

    def jacobian(u, v, w):
        cols = [[0.0, 0.0, 0.0] for _ in range(3)]
        for iu in range(2):
            for iv in range(2):
                for iw in range(2):
                    p = xyz[(iu * 2 + iv) * 2 + iw]
                    factors = (v if iv else 1 - v, w if iw else 1 - w)
                    signs = (1.0 if iu else -1.0, 1.0 if iv else -1.0, 1.0 if iw else -1.0)
                    weights = (signs[0] * factors[0] * factors[1],
                               signs[1] * (u if iu else 1 - u) * factors[1],
                               signs[2] * (u if iu else 1 - u) * factors[0])
                    for axis in range(3):
                        for k in range(3):
                            cols[axis][k] += weights[axis] * p[k]
        return tuple(tuple(c) for c in cols)

    best_planar = (0.0, None)
    min_sine = (math.inf, None)
    best_thickness = (0.0, None)
    best_3d = (0.0, None)
    min_det = (math.inf, None)
    grid = tuple(i / 4 for i in range(5))
    for u in grid:
        for v in grid:
            for w in grid:
                ju, jv, jw = jacobian(u, v, w)
                loc = {"block": block["id"], "param": [u, v, w]}
                hi, lo = _singular_values_2(ju, jv)
                # Use area/sigma_max to avoid subtracting nearly equal eigenvalues.
                area = _norm(_cross(ju, jv))
                lo = area / hi if hi else 0.0
                kp = hi / lo if lo else math.inf
                if kp > best_planar[0]:
                    best_planar = (kp, {**loc, "sigma_max_mm": hi, "sigma_min_mm": lo})
                if u in (0, 1) and v in (0, 1) and w in (0, 1):
                    denom = _norm(ju) * _norm(jv)
                    sine = area / denom if denom else 0.0
                    if sine < min_sine[0]:
                        min_sine = (sine, {**loc, "du_mm": ju, "dv_mm": jv})
                normal = tuple(x / area for x in _cross(ju, jv)) if area else (0.0, 0.0, 0.0)
                thickness = abs(_dot(jw, normal))
                at = hi / thickness if thickness else math.inf
                if at > best_thickness[0]:
                    best_thickness = (at, {**loc, "planar_sigma_max_mm": hi,
                                           "normal_thickness_per_parameter_mm": thickness,
                                           "thickness_tangent_3d_mm": _norm(jw)})
                sv = _singular_values_3((ju, jv, jw))
                k3 = sv[0] / sv[-1] if sv[-1] else math.inf
                if k3 > best_3d[0]:
                    best_3d = (k3, {**loc, "singular_values_mm": sv})
                det = _dot(ju, _cross(jv, jw))
                if det < min_det[0]:
                    min_det = (det, {**loc, "determinant_mm3": det})
    return {"planar_spectral_kappa": best_planar,
            "minimum_corner_sine": min_sine,
            "thickness_direction_aspect": best_thickness,
            "max_3d_spectral_kappa": best_3d,
            "minimum_sampled_jacobian_determinant_mm3": min_det}


def summarize_model(model):
    """Return per-model extrema and deterministic block/parameter witnesses."""
    extrema = [block_metrics(block) for block in model["blocks"]]
    specs = (
        ("planar_spectral_kappa", True),
        ("minimum_corner_sine", False),
        ("thickness_direction_aspect", True),
        ("max_3d_spectral_kappa", True),
        ("minimum_sampled_jacobian_determinant_mm3", False),
    )
    out = {"block_count": len(extrema)}
    for name, maximize in specs:
        row = (max if maximize else min)(extrema, key=lambda item: item[name][0])
        out[name] = row[name][0]
        out[name + "_witness"] = row[name][1]
    return out


def evaluate_manifest(manifest_path, data_root, baseline_cases_path):
    manifest = json.loads(Path(manifest_path).read_text())
    entries = manifest if isinstance(manifest, list) else manifest.get("experiments", [])
    # The archived evaluation freezes the 45 baseline case names. The current
    # gallery also contains later polygonal examples, which are intentionally excluded.
    paired = json.loads(Path(baseline_cases_path).read_text())
    ids = {"bench_poly_" + case["id"] for case in paired["paired_cases"]}
    by_id = {e.get("id"): e for e in entries}
    missing = sorted(ids - by_id.keys())
    if missing:
        raise ValueError(f"gallery manifest misses frozen cases: {', '.join(missing)}")
    selected = [by_id[identifier] for identifier in sorted(ids)]
    rows = []
    for entry in selected:
        rel = entry["model"].split("?", 1)[0]
        model_path = Path(data_root) / Path(rel).name if Path(rel).parent == Path("data") else Path(data_root) / rel.removeprefix("data/")
        model = json.loads(model_path.read_text())
        rows.append({"id": entry["id"], "model": str(model_path),
                     "model_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
                     **summarize_model(model)})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="docs/site/experiments.json")
    parser.add_argument("--models", default="docs/site/data")
    parser.add_argument("--baseline-cases", default="docs/results/sketch-planner-polygon-evaluation/paired.json")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = evaluate_manifest(args.manifest, args.models, args.baseline_cases)
    root = Path.cwd()
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, text=True,
                              capture_output=True, check=False).stdout.strip() or None
    digest = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    out = {"method": "independent trilinear Jacobian sampling; 5x5x5 grid for spectral and aspect metrics, corners for corner sine",
           "sample_points_per_block": 125, "polygonal_extrusion_experiment_count": len(rows),
           "python_version": sys.version.split()[0], "source_revision": revision,
           "manifest_sha256": digest(args.manifest),
           "frozen_case_manifest_sha256": digest(args.baseline_cases),
           "evaluator_sha256": digest(__file__),
           "experiment_results": rows}
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}: {len(rows)} polygonal extrusion experiments")


if __name__ == "__main__":
    main()
