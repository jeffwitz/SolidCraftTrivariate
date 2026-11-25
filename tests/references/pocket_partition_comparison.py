#!/usr/bin/env python3
"""Independent section, incidence, and Euler checks for sketch partitions."""
from __future__ import annotations

import argparse
from collections import Counter
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "tests/references/sketch_partition_corpus.json"
CASE_IDS = (
    "plate_blind_pocket",
    "rectangle_internal_rectangle",
    "rectangle_offset_pocket",
    "c_section_two_holes",
    "l_bracket_polygon_holes_and_pocket",
)


def signed_twice_area(points):
    return sum(points[i][0] * points[(i + 1) % len(points)][1]
               - points[(i + 1) % len(points)][0] * points[i][1]
               for i in range(len(points)))


def polygon_area(points):
    return abs(F(signed_twice_area(points), 2))


def orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(p, a, b):
    return orient(a, b, p) == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) \
        and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def proper_cross(a, b, c, d):
    return (orient(a, b, c) * orient(a, b, d) < 0
            and orient(c, d, a) * orient(c, d, b) < 0)


def point_in_polygon(p, points):
    inside = False
    for a, b in zip(points, points[1:] + points[:1]):
        if on_segment(p, a, b):
            return 0
        if (a[1] > p[1]) != (b[1] > p[1]):
            x = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            if x > p[0]:
                inside = not inside
    return 1 if inside else -1


def polygon_centroid_sample(points):
    return (sum((p[0] for p in points), F(0)) / len(points),
            sum((p[1] for p in points), F(0)) / len(points))


def polygons_overlap(a, b):
    for poly in (a, b):
        turns = [orient(poly[i], poly[(i + 1) % len(poly)], poly[(i + 2) % len(poly)])
                 for i in range(len(poly))]
        if any(turn <= 0 for turn in turns):
            raise AssertionError("SAT oracle requires strictly convex, CCW quads")
    for poly in (a, b):
        for p, q in zip(poly, poly[1:] + poly[:1]):
            normal = (q[1] - p[1], p[0] - q[0])
            pa = [x * normal[0] + y * normal[1] for x, y in a]
            pb = [x * normal[0] + y * normal[1] for x, y in b]
            if max(pa) <= min(pb) or max(pb) <= min(pa):
                return False
    return True


def exact_material_area(recipe, z):
    area = polygon_area([tuple(map(F, p)) for p in recipe["outer"]["vertices"]])
    for cut in recipe.get("cuts", []):
        z0, z1 = map(F, cut["z_interval_mm"])
        if z0 < z < z1:
            area -= polygon_area([tuple(map(F, p)) for p in cut["vertices"]])
    return area


def model_corner_xy(block):
    if block["degrees"] != [1, 1, 1] or len(block["controls"]) != 8:
        raise AssertionError(f"{block['id']}: oracle expects trilinear blocks")
    corners = {}
    for i in range(2):
        for j in range(2):
            for k in range(2):
                index = ((i * 2) + j) * 2 + k
                wx, wy, wz, weight = map(F, block["controls"][index])
                if weight == 0:
                    raise AssertionError(f"{block['id']}: zero corner weight")
                corners[(i, j, k)] = (wx / weight, wy / weight, wz / weight)
    return corners


def jacobian_at(block, uvw):
    """Evaluate the rational trilinear Jacobian from serialized controls."""
    u, v, w = uvw
    params = (u, v, w)
    h = np.zeros(4)
    dh = np.zeros((3, 4))
    for i in range(2):
        for j in range(2):
            for k in range(2):
                bits = (i, j, k)
                basis = np.prod([params[a] if bits[a] else 1.0 - params[a]
                                 for a in range(3)])
                index = ((i * 2) + j) * 2 + k
                control = np.asarray(block["controls"][index], dtype=float)
                h += basis * control
                for axis in range(3):
                    derivative = (1.0 if bits[axis] else -1.0)
                    for other in range(3):
                        if other != axis:
                            derivative *= params[other] if bits[other] else 1.0 - params[other]
                    dh[axis] += derivative * control
    if h[3] <= 0:
        raise AssertionError(f"{block['id']}: non-positive sampled denominator")
    return np.asarray([(dh[a, :3] * h[3] - h[:3] * dh[a, 3]) / h[3] ** 2
                       for a in range(3)])


def max_condition_sample(model):
    samples = list(__import__("itertools").product((0.0, 0.5, 1.0), repeat=3))
    worst = 0.0
    for block in model["blocks"]:
        for uvw in samples:
            singular = np.linalg.svd(jacobian_at(block, uvw), compute_uv=False)
            if singular[-1] <= 0:
                return float("inf")
            worst = max(worst, float(singular[0] / singular[-1]))
    return worst


def layer_cells(model, recipe, z0, z1):
    rows = []
    height = F(recipe["height_mm"])
    shift = tuple(map(F, recipe.get("top_shift_mm", [0, 0])))
    for block in model["blocks"]:
        corners = model_corner_xy(block)
        zs = {p[2] for p in corners.values()}
        if zs != {F(z0), F(z1)}:
            continue
        for i in range(2):
            for j in range(2):
                bottom, top = corners[(i, j, 0)], corners[(i, j, 1)]
                dz = top[2] - bottom[2]
                if (top[0] - bottom[0] != shift[0] * dz / height
                        or top[1] - bottom[1] != shift[1] * dz / height):
                    raise AssertionError(f"{block['id']}: top-face shift does not match recipe")
        # The lower z face is ordered around its four rectangle corners.
        points = []
        for q in ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)):
            x, y, z = corners[q]
            points.append((x - shift[0] * z / height,
                           y - shift[1] * z / height))
        area = F(signed_twice_area(points), 2)
        if area == 0:
            raise AssertionError(f"{block['id']}: collapsed section")
        if area < 0:
            points.reverse()
            area = -area
        rows.append((block["id"], points, area))
    return rows


def planar_euler(rows):
    vertices = set()
    edges = Counter()
    polygons = []
    area_sum = F(0)
    for _, points, area in rows:
        vertices.update(points)
        area_sum += area
        polygons.append(points)
        for a, b in zip(points, points[1:] + points[:1]):
            edges[tuple(sorted((a, b)))] += 1
    if any(n not in (1, 2) for n in edges.values()):
        raise AssertionError("non-manifold planar edge incidence")
    for a, b in edges:
        if any(p not in (a, b) and on_segment(p, a, b) for p in vertices):
            raise AssertionError("nonconforming planar T-junction")
    overlap_pairs = sum(polygons_overlap(polygons[i], polygons[j])
                        for i in range(len(polygons)) for j in range(i + 1, len(polygons)))
    return {"cells": len(rows), "vertices": len(vertices), "edges": len(edges),
            "euler": len(vertices) - len(edges) + len(rows),
            "boundary_edges": sum(n == 1 for n in edges.values()),
            "internal_edges": sum(n == 2 for n in edges.values()),
            "area_exact_mm2": [area_sum.numerator, area_sum.denominator],
            "area_sum_mm2": float(area_sum), "pairwise_overlap_pairs": overlap_pairs,
            "polygons": polygons}


def expected_section(recipe, z):
    outer = [tuple(map(F, p)) for p in recipe["outer"]["vertices"]]
    holes = [[tuple(map(F, p)) for p in c["vertices"]]
             for c in recipe.get("cuts", [])
             if c["z_interval_mm"][0] < z < c["z_interval_mm"][1]]
    return outer, holes


def audit_layer(model, recipe, z0, z1):
    rows = layer_cells(model, recipe, z0, z1)
    metrics = planar_euler(rows)
    zmid = (F(z0) + F(z1)) / 2
    expected_outer, expected_holes = expected_section(recipe, zmid)
    polygons = metrics.pop("polygons")
    # Containment, non-overlap, and exact area conservation together imply
    # coverage up to cell boundaries; all predicates use exact Fractions.
    contained = True
    for points in polygons:
        probes = list(points) + [((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
                                 for a, b in zip(points, points[1:] + points[:1])]
        if any(point_in_polygon(p, expected_outer) < 0 for p in probes):
            contained = False
        if any(point_in_polygon(p, hole) == 1 for hole in expected_holes for p in probes):
            contained = False
        if any(proper_cross(a, b, c, d)
               for a, b in zip(points, points[1:] + points[:1])
               for boundary in [expected_outer, *expected_holes]
               for c, d in zip(boundary, boundary[1:] + boundary[:1])):
            contained = False
        if any(point_in_polygon(polygon_centroid_sample(hole), points) == 1
               for hole in expected_holes):
            contained = False
    overlap_pairs = sum(polygons_overlap(polygons[i], polygons[j])
                        for i in range(len(polygons)) for j in range(i + 1, len(polygons)))
    exact_area = exact_material_area(recipe, zmid)
    area_sum = F(*metrics["area_exact_mm2"])
    gap = abs(area_sum - exact_area)
    metrics.update({"z_interval_mm": [float(z0), float(z1)],
                    "expected_area_exact_mm2": [exact_area.numerator, exact_area.denominator],
                    "area_gap_exact_mm2": [gap.numerator, gap.denominator],
                    "contained_in_expected_material": contained,
                    "pairwise_overlap_pairs": overlap_pairs,
                    "expected_euler_characteristic": 1 - len([
                        c for c in recipe.get("cuts", [])
                        if c["z_interval_mm"][0] < zmid < c["z_interval_mm"][1]])})
    return metrics


def exact_volume(recipe):
    levels = sorted({F(0), F(recipe["height_mm"]),
                     *(F(z) for c in recipe.get("cuts", []) for z in c["z_interval_mm"])})
    return sum((exact_material_area(recipe, (a + b) / 2) * (b - a)
                for a, b in zip(levels, levels[1:])), F(0))


def compare(cli, work, output):
    start = time.monotonic()
    cases = {x["id"]: x for x in json.loads(CORPUS.read_text())}
    offset = copy.deepcopy(cases["rectangle_internal_rectangle"])
    offset["id"] = "rectangle_offset_pocket"
    offset["recipe"]["cuts"][0]["vertices"] = [[7, 1], [13, 1], [13, 5], [7, 5]]
    cases[offset["id"]] = offset
    work.mkdir(parents=True, exist_ok=True)
    cli_sha = hashlib.sha256(cli.read_bytes()).hexdigest()
    results = []
    for name in CASE_IDS:
        item = cases[name]
        if not item["accepted"]:
            raise AssertionError(f"frozen case {name} is not accepted")
        recipe = item["recipe"]
        rp, mp = work / f"{name}.recipe.json", work / f"{name}.model.json"
        rp.write_text(json.dumps(recipe, separators=(",", ":"), allow_nan=False) + "\n")
        proc = subprocess.run([str(cli), "build", str(rp), str(mp)],
                              cwd=ROOT, text=True, capture_output=True)
        if proc.returncode:
            raise RuntimeError(f"CLI build failed for {name}: {proc.stdout}\n{proc.stderr}")
        report = json.loads(proc.stdout)
        model = json.loads(mp.read_text())
        levels = sorted({F(0), F(recipe["height_mm"]),
                         *(F(z) for c in recipe.get("cuts", []) for z in c["z_interval_mm"])})
        layers = [audit_layer(model, recipe, a, b) for a, b in zip(levels, levels[1:])]
        expected_v = exact_volume(recipe)
        actual_v = report.get("volume_exact_mm3")
        volume_ok = actual_v is not None and F(*map(int, actual_v)) == expected_v
        results.append({"id": name, "block_count": len(model["blocks"]),
                        "max_condition_2_sample_27_per_block": max_condition_sample(model),
                        "expected_volume_exact_mm3": [expected_v.numerator, expected_v.denominator],
                        "cli_volume_exact_mm3": actual_v, "volume_match": volume_ok,
                        "layers": layers,
                        "passed": volume_ok and all(
                            x["euler"] == x["expected_euler_characteristic"]
                            and x["area_gap_exact_mm2"] == [0, 1]
                            and x["contained_in_expected_material"]
                            and x["pairwise_overlap_pairs"] == 0 for x in layers)})
        if name == "plate_blind_pocket":
            low = next(x for x in layers if x["z_interval_mm"] == [0, 14])
            crown = next(x for x in layers if x["z_interval_mm"] == [14, 20])
            results[-1]["nine_block_lower_bound"] = {
                "conforming_quad_mesh_assumption": True,
                "lower_slab": {"boundary_vertices_min": 4,
                               "required_interior_reentrant_vertices": 4,
                               "euler_characteristic": 1,
                               "quad_count_bound": "F = V_boundary/2 + V_interior - 1 >= 5",
                               "observed_quad_count": low["cells"],
                               "observed_vertices": low["vertices"],
                               "observed_edges": low["edges"],
                               "observed_euler": low["euler"]},
                "upper_annulus": {"required_boundary_vertices": 8,
                                  "euler_characteristic": 0,
                                  "quad_count_bound": "F = V_boundary/2 + V_interior >= 4",
                                  "observed_quad_count": crown["cells"],
                                  "observed_vertices": crown["vertices"],
                                  "observed_edges": crown["edges"],
                                  "observed_euler": crown["euler"]},
                "observed_total_blocks": len(model["blocks"]),
                "attains_bound": (low["cells"] == 5 and crown["cells"] == 4
                                  and len(model["blocks"]) == 9)}
            results[-1]["passed"] &= results[-1]["nine_block_lower_bound"]["attains_bound"]
        if name == "c_section_two_holes":
            layer = layers[0]
            results[-1]["quad_euler_lower_bound"] = {
                "conforming_quad_mesh_assumption": True,
                "connected_components": 1,
                "holes": 2,
                "boundary_vertices_min": 16,
                "interior_vertices_min": 0,
                "euler_characteristic": -1,
                "formula": "F = V_boundary/2 + V_interior + holes - 1 >= 9",
                "observed_quad_count": layer["cells"],
                "observed_vertices": layer["vertices"],
                "observed_edges": layer["edges"],
                "observed_euler": layer["euler"],
                "expected_volume_mm3": 1568,
                "attains_bound": (layer["cells"] == 9 and len(model["blocks"]) == 9
                                  and expected_v == 1568)}
            results[-1]["passed"] &= results[-1]["quad_euler_lower_bound"]["attains_bound"]
    # A one-millimetre translation of the lower C-section hole tests whether
    # stable block identifiers correspond to actually reusable controls.
    c_recipe = copy.deepcopy(cases["c_section_two_holes"]["recipe"])
    for point in c_recipe["cuts"][0]["vertices"]:
        point[0] += 1
    moved_recipe_path = work / "c_section_lower_hole_x_plus_1.recipe.json"
    moved_model_path = work / "c_section_lower_hole_x_plus_1.model.json"
    moved_recipe_path.write_text(json.dumps(c_recipe, separators=(",", ":"), allow_nan=False) + "\n")
    moved_proc = subprocess.run([str(cli), "build", str(moved_recipe_path), str(moved_model_path)],
                                cwd=ROOT, text=True, capture_output=True)
    if moved_proc.returncode:
        raise RuntimeError("translated lower C-hole build failed: " + moved_proc.stdout + moved_proc.stderr)
    moved_build_report = json.loads(moved_proc.stdout)
    moved_model = json.loads(moved_model_path.read_text())
    base_model = json.loads((work / "c_section_two_holes.model.json").read_text())
    base_blocks = {b["id"]: b for b in base_model["blocks"]}
    moved_blocks = {b["id"]: b for b in moved_model["blocks"]}
    common_ids = set(base_blocks) & set(moved_blocks)
    unchanged_ids = {bid for bid in common_ids
                     if base_blocks[bid]["degrees"] == moved_blocks[bid]["degrees"]
                     and base_blocks[bid]["controls"] == moved_blocks[bid]["controls"]}
    def interfaces(model):
        rows = set()
        for edge in model["topology"]["interfaces"]:
            a, b = tuple(edge["first"]), tuple(edge["second"])
            orientation = edge["orientation"]
            if a <= b:
                key = (a, b, orientation["flip_first"],
                       orientation["flip_second"], orientation["swap"])
            else:
                key = (b, a, orientation["flip_second"],
                       orientation["flip_first"], orientation["swap"])
            rows.add(key)
        return rows
    base_interfaces, moved_interfaces = interfaces(base_model), interfaces(moved_model)
    unchanged_interface_rows = {edge for edge in base_interfaces & moved_interfaces
                                 if edge[0][0] in unchanged_ids and edge[1][0] in unchanged_ids}
    old_lower_hole_corners = {tuple(map(F, point))
                              for point in cases["c_section_two_holes"]["recipe"]["cuts"][0]["vertices"]}
    expected_unchanged_ids = set()
    for bid, block in base_blocks.items():
        corners = model_corner_xy(block)
        if not any(point[:2] in old_lower_hole_corners for point in corners.values()):
            expected_unchanged_ids.add(bid)
    locality_qualified = (unchanged_ids == expected_unchanged_ids
                          and len(expected_unchanged_ids) == 5
                          and len(unchanged_interface_rows) == 5)
    moved_z0, moved_z1 = 0, int(c_recipe["height_mm"])
    moved_layer = audit_layer(moved_model, c_recipe, moved_z0, moved_z1)
    moved_volume = exact_volume(c_recipe)
    base_c = next(x for x in results if x["id"] == "c_section_two_holes")
    results.append({"id": "c_section_lower_hole_x_plus_1", "block_count": len(moved_blocks),
                    "max_condition_2_sample_27_per_block": max_condition_sample(moved_model),
                    "expected_volume_exact_mm3": [moved_volume.numerator, moved_volume.denominator],
                    "cli_volume_exact_mm3": moved_build_report.get("volume_exact_mm3"),
                    "layers": [moved_layer],
                    "preservation_vs_original": {
                        "translation_mm": [1, 0],
                        "common_block_ids": len(common_ids),
                        "unchanged_id_and_degree_and_control_blocks": len(unchanged_ids),
                        "unchanged_block_ids": sorted(unchanged_ids),
                        "expected_unchanged_block_ids": sorted(expected_unchanged_ids),
                        "changed_control_ids": sorted(common_ids - unchanged_ids),
                        "baseline_interface_count": len(base_interfaces),
                        "translated_interface_count": len(moved_interfaces),
                        "exact_shared_interfaces_between_unchanged_blocks": len(unchanged_interface_rows),
                        "locality_qualified": locality_qualified},
                    "baseline_blocks": base_c["block_count"],
                    "volume_match": (moved_build_report.get("volume_exact_mm3") ==
                                     [moved_volume.numerator, moved_volume.denominator]),
                    "passed": (len(moved_blocks) == 9 and moved_layer["euler"] == -1
                               and moved_layer["area_gap_exact_mm2"] == [0, 1]
                               and moved_layer["contained_in_expected_material"]
                               and moved_layer["pairwise_overlap_pairs"] == 0
                               and locality_qualified
                               and moved_build_report.get("volume_exact_mm3") ==
                               [moved_volume.numerator, moved_volume.denominator])})
    if hashlib.sha256(cli.read_bytes()).hexdigest() != cli_sha:
        raise RuntimeError("CLI executable changed during oracle run")
    report = {"schema": "pocket_partition_comparison_v1",
              "provenance": {"cli": str(cli), "cli_sha256": cli_sha,
                             "python": __import__("platform").python_version(),
                             "elapsed_s": time.monotonic() - start},
              "cases": results,
              "limits": ["polygonal orthogonal sections only",
                         "section coverage is established from exact corner polygons, containment, non-overlap and area conservation",
                         "Euler checks prove a combinatorial lower bound only under conforming quadrilateral cells and required corner vertices",
                         "Jacobian samples are intentionally not used as continuous certificates"]}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.cli.resolve(), args.work, args.output)
    for row in result["cases"]:
        print(row["id"], f"blocks={row['block_count']}", "PASS" if row["passed"] else "FAIL")
    raise SystemExit(int(not all(row["passed"] for row in result["cases"])))
