#!/usr/bin/env python3
"""Independent planar and OCC oracle for polygonal sketch layer partitions."""
from __future__ import annotations

import FreeCAD as App
import Part

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
import platform
from pathlib import Path
import random
import subprocess
import time

import numpy as np
from shapely.geometry import Polygon, LineString, Point
from shapely.ops import unary_union
from coaxial_program import compile_mass_oracle

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CORPUS = HERE / "sketch_partition_corpus.json"
SEED = 2026112525


def loop_points(loop):
    return [tuple(map(int, p)) for p in loop["vertices"]]


def polygon_moments(points):
    """Return positive exact area and first moments for an integer polygon."""
    cross = [points[i][0] * points[(i + 1) % len(points)][1]
             - points[(i + 1) % len(points)][0] * points[i][1]
             for i in range(len(points))]
    twice = sum(cross)
    if twice == 0:
        raise ValueError("degenerate loop in frozen corpus")
    sign = 1 if twice > 0 else -1
    area = F(sign * twice, 2)
    mx = F(sign * sum((points[i][0] + points[(i + 1) % len(points)][0]) * c
                      for i, c in enumerate(cross)), 6)
    my = F(sign * sum((points[i][1] + points[(i + 1) % len(points)][1]) * c
                      for i, c in enumerate(cross)), 6)
    return area, mx, my


def exact_reference(recipe):
    """Exact volume and first moments, independently from the block model."""
    h = int(recipe["height_mm"])
    shift = recipe.get("top_shift_mm", [0, 0])
    outer_area, outer_mx, outer_my = polygon_moments(loop_points(recipe["outer"]))
    cuts = recipe["cuts"]
    cut_metrics = [(c, polygon_moments(loop_points(c))) for c in cuts]
    levels = sorted({0, h, *(z for c in cuts for z in c["z_interval_mm"])})
    volume = mx = my = mz = F(0)
    layers = []
    for z0, z1 in zip(levels, levels[1:]):
        mid = F(z0 + z1, 2)
        area, base_mx, base_my = outer_area, outer_mx, outer_my
        active = []
        for cut, (a, cx, cy) in cut_metrics:
            lo, hi = cut["z_interval_mm"]
            if lo < mid < hi:
                area -= a
                base_mx -= cx
                base_my -= cy
                active.append(cut["id"])
        if area <= 0:
            raise ValueError("non-positive material area in frozen corpus")
        dz = F(z1 - z0)
        iz = F(z1 * z1 - z0 * z0, 2)
        volume += area * dz
        mx += base_mx * dz + area * F(shift[0], h) * iz
        my += base_my * dz + area * F(shift[1], h) * iz
        mz += area * iz
        layers.append({"z_interval_mm": [z0, z1], "active_cuts": active,
                       "area_mm2": [area.numerator, area.denominator]})
    return {"volume_mm3": volume, "centroid_mm": [mx / volume, my / volume, mz / volume],
            "outer_area_mm2": outer_area, "layers": layers}


def shapely_region(recipe, z):
    outer = Polygon(loop_points(recipe["outer"]))
    holes = [Polygon(loop_points(c)) for c in recipe["cuts"]
             if c["z_interval_mm"][0] < z < c["z_interval_mm"][1]]
    return outer.difference(unary_union(holes)) if holes else outer


def expected_material(recipe, x, y, z):
    h = recipe["height_mm"]
    if not 0 < z < h:
        return False
    pt = Point(x - recipe.get("top_shift_mm", [0, 0])[0] * z / h,
               y - recipe.get("top_shift_mm", [0, 0])[1] * z / h)
    if not Polygon(loop_points(recipe["outer"])).contains(pt):
        return False
    for cut in recipe["cuts"]:
        lo, hi = cut["z_interval_mm"]
        region = Polygon(loop_points(cut))
        if lo < z < hi and region.covers(pt):
            return False
        if z in (lo, hi) and region.covers(pt):
            return False
    return True


def reference_state(recipe, x, y, z, tolerance=1e-9):
    """Return IN/OUT/ON from independent polygon and layer predicates."""
    h = recipe["height_mm"]
    if z < -tolerance or z > h + tolerance:
        return "OUT"
    dx, dy = recipe.get("top_shift_mm", [0, 0])
    point = Point(x - dx * z / h, y - dy * z / h)
    outer = Polygon(loop_points(recipe["outer"]))
    if outer.boundary.distance(point) <= tolerance:
        return "ON"
    if not outer.contains(point):
        return "OUT"
    if abs(z) <= tolerance or abs(z - h) <= tolerance:
        for cut in recipe["cuts"]:
            lo, hi = cut["z_interval_mm"]
            region = Polygon(loop_points(cut))
            if (abs(z - lo) <= tolerance or abs(z - hi) <= tolerance) and region.covers(point):
                if region.boundary.distance(point) <= tolerance:
                    return "ON"
                if (lo == 0 and abs(z) <= tolerance) or (hi == h and abs(z - h) <= tolerance):
                    return "OUT"
                return "ON"
        return "ON"
    for cut in recipe["cuts"]:
        lo, hi = cut["z_interval_mm"]
        region = Polygon(loop_points(cut))
        if lo - tolerance <= z <= hi + tolerance:
            if region.boundary.distance(point) <= tolerance:
                return "ON"
            if lo < z < hi and region.contains(point):
                return "OUT"
            if (abs(z - lo) <= tolerance or abs(z - hi) <= tolerance) and region.covers(point):
                return "ON"
    return "IN"


def canonical_loop(points):
    points = tuple(tuple(p) for p in points)
    variants = []
    for seq in (points, tuple(reversed(points))):
        variants.extend(seq[i:] + seq[:i] for i in range(len(seq)))
    return min(variants)


def physical_signature(recipe):
    cuts = sorted((canonical_loop(loop_points(c)), tuple(c["z_interval_mm"])) for c in recipe["cuts"])
    return (canonical_loop(loop_points(recipe["outer"])), recipe["height_mm"],
            tuple(recipe.get("top_shift_mm", [0, 0])), tuple(cuts))


def independent_cad(recipe):
    """Build a separate OCC prism-minus-prisms solid from the sketch loops."""
    h = recipe["height_mm"]
    dx, dy = recipe.get("top_shift_mm", [0, 0])
    def face(points, z):
        pts = [App.Vector(x + dx * z / h, y + dy * z / h, z) for x, y in points]
        wire = Part.makePolygon([*pts, pts[0]])
        return Part.Face(wire)
    shape = face(loop_points(recipe["outer"]), 0).extrude(App.Vector(dx, dy, h))
    for cut in recipe["cuts"]:
        lo, hi = cut["z_interval_mm"]
        cutter = face(loop_points(cut), lo).extrude(App.Vector(dx * (hi - lo) / h,
                                                                dy * (hi - lo) / h,
                                                                hi - lo))
        shape = shape.cut(cutter)
    return shape.removeSplitter()


def random_recipes(base):
    """Deterministic edits: reversals, cyclic starts, labels, guides and rejects."""
    rng = random.Random(SEED)
    # The seeded corpus was frozen before the dangling-guide refusal. Each
    # randomized case removes its guides, so retain this source rectangle.
    valid_base = [item for item in base if item["accepted"] or item["id"] == "guide_dangling_rejected"]
    for i in range(24):
        source = valid_base[i % len(valid_base)]
        r = json.loads(json.dumps(source["recipe"]))
        r["construction"] = []
        r["height_mm"] = rng.randint(6, 16)
        dx, dy = r.get("top_shift_mm", [0, 0])
        if i % 2:
            r["top_shift_mm"] = [dx + rng.randint(-3, 3), dy + rng.randint(-3, 3)]
        if i % 3 == 0:
            for loop in [r["outer"], *r["cuts"]]:
                pts = loop["vertices"]
                loop["vertices"] = pts[1:] + pts[:1]
        if i % 3 == 1:
            for loop in [r["outer"], *r["cuts"]]:
                loop["vertices"] = list(reversed(loop["vertices"]))
        for c in r["cuts"]:
            lo = rng.randint(0, r["height_mm"] - 2)
            hi = rng.randint(lo + 1, r["height_mm"])
            c["z_interval_mm"] = [lo, hi]
        expected = True
        if i % 6 == 5:
            expected = False
            r["height_mm"] = 0
        yield {"id": f"seed-{i:02d}", "accepted": expected, "recipe": r}


def basis(n, t):
    return [math.comb(n, i) * t ** i * (1 - t) ** (n - i) for i in range(n + 1)]


def basis_d(n, t):
    if n == 0:
        return [0.]
    low = basis(n - 1, t)
    return [n * ((low[i - 1] if i else 0.) - (low[i] if i < n else 0.))
            for i in range(n + 1)]


def jacobian_at(block, uvw):
    deg = block["degrees"]
    net = np.asarray(block["controls"], dtype=np.float64).reshape(deg[0] + 1, deg[1] + 1,
                                                                   deg[2] + 1, 4)
    bs = [basis(n, t) for n, t in zip(deg, uvw)]
    ds = [basis_d(n, t) for n, t in zip(deg, uvw)]
    h = np.zeros(4); dh = np.zeros((3, 4))
    for i in range(deg[0] + 1):
        for j in range(deg[1] + 1):
            for k in range(deg[2] + 1):
                p = net[i, j, k]
                h += bs[0][i] * bs[1][j] * bs[2][k] * p
                dh[0] += ds[0][i] * bs[1][j] * bs[2][k] * p
                dh[1] += bs[0][i] * ds[1][j] * bs[2][k] * p
                dh[2] += bs[0][i] * bs[1][j] * ds[2][k] * p
    if h[3] <= 0:
        raise ValueError("nonpositive evaluated homogeneous weight")
    jac = (dh[:, :3] * h[3] - dh[:, 3, None] * h[:3]) / h[3] ** 2
    return jac


def corner_data(block):
    deg = block["degrees"]
    shape = tuple(d + 1 for d in deg)
    net = np.asarray(block["controls"], dtype=object).reshape(*shape, 4)
    out = {}
    for i in (0, deg[0]):
        for j in (0, deg[1]):
            for k in (0, deg[2]):
                p = net[i, j, k]
                out[(i > 0, j > 0, k > 0)] = tuple(F(int(a), int(p[3])) for a in p[:3])
    return out


def incidence(model):
    vertices, edges, faces = set(), set(), Counter()
    det_min = math.inf
    per_block = []
    for block in model["blocks"]:
        c = corner_data(block)
        verts = list(c)
        projective = {bits: c[bits] for bits in verts}
        vertices.update(projective.values())
        for axis in range(3):
            for bits in verts:
                if bits[axis]:
                    continue
                other = list(bits); other[axis] = True
                edges.add(tuple(sorted((projective[bits], projective[tuple(other)]))))
        for axis in range(3):
            for side in (False, True):
                q = tuple(sorted(projective[b] for b in verts if b[axis] == side))
                faces[q] += 1
        dets = [float(np.linalg.det(jacobian_at(block, p))) for p in itertools.product((0., .5, 1.), repeat=3)]
        local_min = min(dets)
        det_min = min(det_min, local_min)
        per_block.append({"id": block["id"], "determinant_sample_min_mm3": local_min})
    return {"vertices": len(vertices), "edges": len(edges), "faces": len(faces),
            "cells": len(model["blocks"]), "euler": len(vertices) - len(edges) + len(faces) - len(model["blocks"]),
            "max_face_multiplicity": max(faces.values(), default=0),
            "shared_faces": sum(v == 2 for v in faces.values()),
            "exterior_faces": sum(v == 1 for v in faces.values()),
            "determinant_sample_min_mm3": det_min, "blocks": per_block}


def block_faces(block):
    """Yield logical faces as (axis,side, four XYZ points) from corner controls."""
    c = corner_data(block)
    for axis in range(3):
        for side in (False, True):
            face = [bits for bits in c if bits[axis] == side]
            # Order remaining two axes around the quad, independent of axis.
            other = [a for a in range(3) if a != axis]
            ordered = []
            for first, second in ((False, False), (True, False), (True, True), (False, True)):
                bits = [False, False, False]
                bits[axis] = side; bits[other[0]] = first; bits[other[1]] = second
                ordered.append(tuple(bits))
            yield axis, side, [tuple(float(q) for q in c[b]) for b in ordered]


def coverage(model, recipe, reference):
    """Compare each rational block footprint to an independent Shapely section."""
    h = recipe["height_mm"]
    shift = recipe.get("top_shift_mm", [0, 0])
    levels = sorted({0, h, *(z for cut in recipe["cuts"] for z in cut["z_interval_mm"])})
    result = []
    max_gap = max_overlap = max_area_error = 0.
    by_layer = {}
    for block in model["blocks"]:
        c = corner_data(block)
        zs = sorted({float(q[2]) for q in c.values()})
        if len(zs) != 2:
            raise AssertionError(f"block {block['id']} does not span a single layer")
        layer = (round(zs[0], 10), round(zs[1], 10))
        base = []
        for bits in ((False, False, False), (True, False, False), (True, True, False), (False, True, False)):
            x, y, z = map(float, c[bits]); base.append((x - shift[0] * z / h, y - shift[1] * z / h))
        poly = Polygon(base)
        if not poly.is_valid or poly.area <= 0:
            raise AssertionError(f"degenerate projected cell {block['id']}")
        by_layer.setdefault(layer, []).append((block["id"], poly))
    for z0, z1 in zip(levels, levels[1:]):
        rows = by_layer.get((float(z0), float(z1)), [])
        if not rows:
            raise AssertionError(f"missing layer {(z0, z1)}")
        zmid = (z0 + z1) / 2
        expected = shapely_region(recipe, zmid)
        polys = [p for _, p in rows]
        union = unary_union(polys)
        gap = expected.symmetric_difference(union).area
        overlap = sum(polys[i].intersection(polys[j]).area
                      for i in range(len(polys)) for j in range(i + 1, len(polys)))
        area_error = abs(sum(p.area for p in polys) - expected.area)
        max_gap = max(max_gap, gap); max_overlap = max(max_overlap, overlap); max_area_error = max(max_area_error, area_error)
        result.append({"z_interval_mm": [z0, z1], "cell_count": len(rows),
                       "reference_area_mm2": expected.area, "block_area_sum_mm2": sum(p.area for p in polys),
                       "symmetric_difference_mm2": gap, "pairwise_overlap_mm2": overlap})
    triangle_baseline = 0
    for layer in result:
        z0, z1 = layer["z_interval_mm"]
        zmid = (z0 + z1) / 2
        active = [c for c in recipe["cuts"] if c["z_interval_mm"][0] < zmid < c["z_interval_mm"][1]]
        triangle_baseline += 3 * (len(recipe["outer"]["vertices"])
                                  + sum(len(c["vertices"]) for c in active)
                                  + 2 * len(active) - 2)
    return {"layers": result, "triangle_prism_cell_baseline": triangle_baseline,
            "partition_block_count": len(model["blocks"]),
            "max_symmetric_difference_mm2": max_gap,
            "max_pairwise_overlap_mm2": max_overlap, "max_area_error_mm2": max_area_error}


def cad_classification(recipe, points):
    rng = random.Random(SEED)
    outer = np.asarray(loop_points(recipe["outer"]), dtype=float)
    low = outer.min(axis=0) - 1; high = outer.max(axis=0) + 1
    h = recipe["height_mm"]
    probes = [(rng.uniform(low[0], high[0]), rng.uniform(low[1], high[1]), rng.uniform(-.1, h + .1))
              for _ in range(256)]
    for cut in recipe["cuts"]:
        poly = np.asarray(loop_points(cut), dtype=float)
        center = poly.mean(axis=0)
        for lo, hi in (cut["z_interval_mm"],):
            for z in (lo - 1e-4, lo + 1e-4, hi - 1e-4, hi + 1e-4,
                      (lo + hi) / 2):
                for p in (*poly, center):
                    probes.append((float(p[0]), float(p[1]), z))
    shift = recipe.get("top_shift_mm", [0, 0])
    bad = []
    boundary = Part.makeCompound(points.Shells)
    for x, y, z in probes:
        world_xyz = [x + shift[0] * z / h, y + shift[1] * z / h, z]
        expected = reference_state(recipe, *world_xyz)
        world = App.Vector(*world_xyz)
        distance = Part.Vertex(world).distToShape(boundary)[0]
        observed = ("ON" if distance <= 1e-7 else
                    ("IN" if points.isInside(world, 1e-8, False) else "OUT"))
        if observed != expected:
            bad.append({"xyz_mm": world_xyz, "expected_state": expected,
                        "observed_state": observed, "distance_to_step_mm": distance})
    return len(probes), bad


def verify_provenance(model, recipe):
    expected = {"start", "end", f"outer/{recipe['outer']['id']}"}
    for cut in recipe["cuts"]:
        expected.add(f"cut/{cut['id']}/wall")
        lo, hi = cut["z_interval_mm"]
        if lo > 0:
            expected.add(f"cut/{cut['id']}/lower")
        if hi < recipe["height_mm"]:
            expected.add(f"cut/{cut['id']}/upper")
    tags = Counter(face["feature"] for face in model["design_faces"])
    unknown = sorted(set(tags) - expected)
    missing = sorted(expected - set(tags))
    return {"expected_features": sorted(expected), "counts": dict(tags), "missing": missing,
            "unexpected": unknown, "qualified": not missing and not unknown}


def run(cli: Path, work: Path):
    start = time.monotonic()
    cli_sha = hashlib.sha256(cli.read_bytes()).hexdigest()
    work.mkdir(parents=True, exist_ok=True)
    mass_helper = compile_mass_oracle(work)
    frozen = json.loads(CORPUS.read_text())
    cases = [*frozen, *random_recipes(frozen)]
    model_path, recipe_path = work / "model.json", work / "recipe.json"
    step_rows = []; report_rows = []; digests = {}
    physical_signatures = {}
    for item in cases:
        recipe = item["recipe"]; name = item["id"]
        recipe_path.write_text(json.dumps(recipe, allow_nan=False))
        previous = model_path.read_bytes() if model_path.exists() else None
        build_started = time.monotonic()
        proc = subprocess.run([str(cli), "build", str(recipe_path), str(model_path)],
                              capture_output=True, text=True)
        build_elapsed = time.monotonic() - build_started
        try:
            native = json.loads(proc.stdout)
        except json.JSONDecodeError:
            native = {"stdout": proc.stdout, "stderr": proc.stderr, "returncode": proc.returncode}
        row = {"id": name, "expected_accepted": item["accepted"], "recipe": recipe,
               "native_build": native, "build_audit_s": build_elapsed,
               "partition_strategy": recipe.get("partition_strategy", "automatic"),
               "passed": (proc.returncode == 0) == item["accepted"]}
        if proc.returncode:
            row["previous_model_preserved"] = (model_path.read_bytes() if model_path.exists() else None) == previous
            row["passed"] &= row["previous_model_preserved"]
            report_rows.append(row)
            continue
        if not item["accepted"]:
            row["passed"] = False
            report_rows.append(row)
            continue
        model = json.loads(model_path.read_text())
        ref = exact_reference(recipe)
        expected_pair = [ref["volume_mm3"].numerator, ref["volume_mm3"].denominator]
        row["volume_exact_expected_mm3"] = expected_pair
        exact = native.get("volume_exact_mm3")
        exact_ok = exact is not None and F(int(exact[0]), int(exact[1])) == ref["volume_mm3"]
        metrics = incidence(model)
        cover = coverage(model, recipe, ref)
        provenance = verify_provenance(model, recipe)
        digest = hashlib.sha256(json.dumps([model[k] for k in ("blocks", "topology", "design_faces")],
                                           sort_keys=True).encode()).hexdigest()
        digests[name] = digest
        row.update(native_volume_exact=exact, volume_exact_match=exact_ok, block_count=len(model["blocks"]),
                   topology_oracle=metrics, cross_section_coverage=cover,
                   face_provenance=provenance, partition_sha256=digest)
        row["passed"] &= bool(native.get("method") == "certified_sketch_layer_partition" and exact_ok)
        row["passed"] &= metrics["max_face_multiplicity"] <= 2 and metrics["determinant_sample_min_mm3"] > 0
        expected_chi = 1 + sum(0 < c["z_interval_mm"][0] and c["z_interval_mm"][1] < recipe["height_mm"]
                               for c in recipe["cuts"]) - sum(c["z_interval_mm"] == [0, recipe["height_mm"]]
                                                               for c in recipe["cuts"])
        row["expected_euler"] = expected_chi
        row["topology_match"] = (metrics["euler"] == expected_chi
                                 and metrics["shared_faces"] == native.get("shared_faces")
                                 and metrics["exterior_faces"] == native.get("exterior_sides"))
        row["passed"] &= row["topology_match"]
        row["passed"] &= cover["max_symmetric_difference_mm2"] < 1e-8
        row["passed"] &= cover["max_pairwise_overlap_mm2"] < 1e-8
        row["passed"] &= cover["max_area_error_mm2"] < 1e-8 and provenance["qualified"]
        if "equivalent_to" in item:
            row["same_partition_as_reference"] = digest == digests[item["equivalent_to"]]
            row["passed"] &= row["same_partition_as_reference"]
        signature = physical_signature(recipe)
        physical_signatures[name] = signature
        if "equivalent_geometry_to" in item:
            row["same_physical_region"] = signature == physical_signatures[item["equivalent_geometry_to"]]
            row["passed"] &= row["same_physical_region"]
        elif not name.startswith("seed-"):
            step = work / f"{name}.step"
            export_started = time.monotonic()
            exp = subprocess.run([str(cli), "export-step", str(model_path), str(step)],
                                 capture_output=True, text=True)
            row["step_export_s"] = time.monotonic() - export_started
            try:
                row["step_export"] = json.loads(exp.stdout)
            except json.JSONDecodeError:
                row["step_export"] = {"stdout": exp.stdout, "stderr": exp.stderr}
            row["passed"] &= exp.returncode == 0
            if exp.returncode == 0:
                native_shape = Part.read(str(step))
                reference_shape = independent_cad(recipe)
                expected_volume = float(ref["volume_mm3"])
                expected_center = np.asarray([float(x) for x in ref["centroid_mm"]])
                shapes = [reference_shape, native_shape]
                masses = []
                for i, shape in enumerate(shapes):
                    brep = work / f"{name}-{i}.brep"
                    shape.exportBrep(str(brep))
                    measured = json.loads(subprocess.check_output([str(mass_helper), str(brep), "--centroid"], text=True))[-1]
                    masses.append({"volume_mm3": measured["volume_mm3"],
                                   "centroid_mm": measured["centroid_mm"],
                                   "valid": bool(shape.isValid()), "solids": len(shape.Solids),
                                   "shells": len(shape.Shells)})
                nprobe, bad = cad_classification(recipe, native_shape)
                ev = max(abs(m["volume_mm3"] / expected_volume - 1) for m in masses)
                ec = max(float(np.linalg.norm(np.asarray(m["centroid_mm"]) - expected_center)) for m in masses)
                row.update(independent_reference_mass=masses[0], step_mass=masses[1],
                           classification_points=nprobe, classification_mismatches=len(bad),
                           classification_witnesses=bad[:8], volume_relative_error=ev,
                           centroid_error_mm=ec, solid_count_match=len(native_shape.Solids) == 1)
                row["passed"] &= not bad and ev < 1e-8 and ec < 1e-7
                row["passed"] &= all(m["valid"] for m in masses) and masses[1]["solids"] == 1
        report_rows.append(row)
        print(name, row["passed"], flush=True)
    if hashlib.sha256(cli.read_bytes()).hexdigest() != cli_sha:
        raise RuntimeError("CLI executable changed during oracle run")
    strategy_pairs = []
    by_id = {r["id"]: r for r in report_rows}
    for auto_name, tri_name in (("plate_blind_pocket", "plate_blind_pocket_triangulate"),
                                ("c_section_two_holes", "c_section_two_holes_triangulate")):
        auto, tri = by_id[auto_name], by_id[tri_name]
        ac, tc = auto["block_count"], tri["block_count"]
        strategy_pairs.append({"case": auto_name, "automatic_blocks": ac, "triangulation_blocks": tc,
                               "block_ratio_automatic_over_triangulation": ac / tc,
                               "automatic_build_audit_s": auto["build_audit_s"],
                               "triangulation_build_audit_s": tri["build_audit_s"],
                               "same_exact_volume": auto["volume_exact_expected_mm3"] == tri["volume_exact_expected_mm3"],
                               "automatic_min_determinant_sample_mm3": auto["topology_oracle"]["determinant_sample_min_mm3"],
                               "triangulation_min_determinant_sample_mm3": tri["topology_oracle"]["determinant_sample_min_mm3"],
                               "both_cover_reference": (auto["cross_section_coverage"]["max_symmetric_difference_mm2"] < 1e-8
                                                         and tri["cross_section_coverage"]["max_symmetric_difference_mm2"] < 1e-8)})
    return {"schema": "sketch_layer_partition_independent_oracle_v1",
            "provenance": {"cli": str(cli), "cli_sha256": cli_sha,
                           "checkout_revision": subprocess.check_output(["rtk", "git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                           "python": platform.python_version(), "occ": Part.OCC_VERSION,
                           "shapely": __import__("shapely").__version__, "seed": SEED,
                           "elapsed_s": time.monotonic() - start},
            "protocol": {"reference": "exact shoelace moments and layer integrals; independent polygonal OCC difference",
                         "coverage": "Shapely union and pairwise intersections of de-sheared block footprints",
                         "classification": "IN/OUT/ON against exact polygon predicates and OCC distance to the imported STEP boundary",
                         "limits": ["polygonal loops and straight finite guides only", "sampled Jacobian is diagnostic, not a continuous proof"]},
            "legacy_boolean_boundary_diagnostic": json.loads(
                (ROOT / "tests/regressions/sketch_partition_boolean_boundary.json").read_text()),
            "case_count": len(report_rows), "qualified_cases": sum(bool(r["passed"]) for r in report_rows),
            "strategy_comparisons": strategy_pairs,
            "cases": report_rows, "qualified": all(r["passed"] for r in report_rows)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.cli.resolve(), args.work)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    raise SystemExit(int(not report["qualified"]))
