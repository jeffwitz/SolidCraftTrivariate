#!/usr/bin/env python3
"""Independent OCC and analytic oracle for intersecting orthogonal bores."""
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
from coaxial_program import compile_mass_oracle

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CORPUS = HERE / "steinmetz_partition_corpus.json"
SEED = 2026112526
BOUNDARY_TOL_MM = 1e-7


def dimensions(recipe):
    vertices = recipe["outer"]["vertices"]
    x = [p[0] for p in vertices]; y = [p[1] for p in vertices]
    return max(x) - min(x), max(y) - min(y), recipe["height_mm"], min(x), min(y)


def bore_terms(recipe):
    lx, ly, lz, x0, y0 = dimensions(recipe)
    if len(recipe["bores"]) != 2:
        raise ValueError("frozen analytic cases require exactly two bores")
    first_bore, second_bore = recipe["bores"]
    ax, az = ("xyz".index(b["axis"]) for b in (first_bore, second_bore))
    if ax == az:
        raise ValueError("frozen analytic cases require perpendicular bore axes")
    normal = 3 - ax - az
    rx, rz = first_bore["radius_mm"], second_bore["radius_mm"]
    if rx != rz:
        raise ValueError("frozen analytic cases require equal radii")
    r = rx
    crossing = first_bore["center_mm"][normal] == second_bore["center_mm"][normal]
    if not crossing:
        separation = abs(first_bore["center_mm"][normal] - second_bore["center_mm"][normal])
        if separation < rx + rz:
            raise ValueError("non-frozen offset-cylinder intersection needs a separate quadrature oracle")
    intersection = F(16 * r**3, 3) if crossing else F(0)
    lengths = [lx, ly, lz]
    starts = [x0, y0, 0]
    box_center = np.asarray([starts[i] + lengths[i] / 2 for i in range(3)], dtype=float)
    def cylinder_center(bore):
        center = np.asarray(bore["center_mm"], dtype=float).copy()
        axis = "xyz".index(bore["axis"])
        center[axis] = box_center[axis]
        return center
    c0, c1 = cylinder_center(first_bore), cylinder_center(second_bore)
    v0, v1 = math.pi * rx**2 * lengths[ax], math.pi * rz**2 * lengths[az]
    overlap_center = np.zeros(3, dtype=float)
    overlap_center[ax] = second_bore["center_mm"][ax]
    overlap_center[az] = first_bore["center_mm"][az]
    overlap_center[normal] = first_bore["center_mm"][normal]
    box_volume = lx * ly * lz
    # The Steinmetz overlap is a rational volume; cylinder terms carry pi.
    removed_exact_constant = intersection
    rational_volume = F(box_volume) + removed_exact_constant
    pi_coefficient = -(rx**2 * lengths[ax] + rz**2 * lengths[az])
    first = box_volume * box_center
    m_x = v0 * c0
    m_z = v1 * c1
    m_i = float(intersection) * overlap_center
    # At a crossing the overlap is symmetric around the common normal point.
    # For disjoint bores, the intersection term is exactly zero.
    centroid = (first - m_x - m_z + m_i) / (float(rational_volume) + math.pi * pi_coefficient)
    return {"box_volume_mm3": F(box_volume), "intersection_mm3": intersection,
            "volume_rational_mm3": rational_volume, "pi_coefficient_mm3": pi_coefficient,
            "volume_mm3": float(rational_volume) + math.pi * pi_coefficient,
            "centroid_mm": centroid, "crossing": crossing,
            "twice_rational_part_mm3": 2 * rational_volume,
            "four_pi_numerator": 4 * pi_coefficient,
            "cylinder_lengths_mm": [lengths[ax], lengths[az]],
            "axes": [first_bore["axis"], second_bore["axis"]],
            "overlap_center_mm": overlap_center}


def reference_shape(recipe):
    lx, ly, lz, x0, y0 = dimensions(recipe)
    outline = [App.Vector(x, y, 0) for x, y in recipe["outer"]["vertices"]]
    solid = Part.Face(Part.makePolygon([*outline, outline[0]])).extrude(App.Vector(0, 0, lz))
    for bore in recipe["bores"]:
        axis = bore["axis"]
        center = bore["center_mm"]
        if axis == "x":
            start, direction, length = App.Vector(x0, center[1], center[2]), App.Vector(1, 0, 0), lx
        elif axis == "z":
            start, direction, length = App.Vector(center[0], center[1], 0), App.Vector(0, 0, 1), lz
        else:
            start, direction, length = App.Vector(center[0], y0, center[2]), App.Vector(0, 1, 0), ly
        tool = Part.makeCylinder(bore["radius_mm"], length, start, direction)
        solid = solid.cut(tool)
    return solid.removeSplitter()


def expected_state(recipe, xyz, tol=1e-8):
    x, y, z = xyz
    lx, ly, lz, x0, y0 = dimensions(recipe)
    if x < x0 - tol or x > x0 + lx + tol or y < y0 - tol or y > y0 + ly + tol or z < -tol or z > lz + tol:
        return "OUT"
    radial = []
    for b in recipe["bores"]:
        c = b["center_mm"]; axis = b["axis"]
        q = [i for i in range(3) if "xyz"[i] != axis]
        distance = math.hypot(xyz[q[0]] - c[q[0]], xyz[q[1]] - c[q[1]])
        radial.append((distance, b["radius_mm"]))
    # A nominal cylinder surface hidden inside another bore is not part of
    # the boundary of the material complement.
    if any(d < r - tol for d, r in radial):
        return "OUT"
    if any(abs(d - r) <= tol for d, r in radial):
        return "ON"
    if min(abs(x - x0), abs(x - (x0 + lx)), abs(y - y0), abs(y - (y0 + ly)), abs(z), abs(z - lz)) <= tol:
        return "ON"
    return "IN"


def probe_points(recipe, count=320):
    lx, ly, lz, x0, y0 = dimensions(recipe)
    rng = random.Random(SEED)
    points = [(rng.uniform(x0, x0 + lx), rng.uniform(y0, y0 + ly), rng.uniform(0, lz))
              for _ in range(count)]
    bx, bz = recipe["bores"]
    ax, az = ("xyz".index(b["axis"]) for b in (bx, bz))
    normal = 3 - ax - az
    r = bx["radius_mm"]
    cross = [0., 0., 0.]
    cross[ax] = bz["center_mm"][ax]
    cross[az] = bx["center_mm"][az]
    cross[normal] = bx["center_mm"][normal]
    # Surface intersections and their one-sided perturbations.
    for sx, sz in itertools.product((-1, 1), repeat=2):
        p = cross.copy(); p[ax] += sx * r; p[az] += sz * r
        points.append(tuple(p))
        for delta in (-1e-4, 1e-4):
            plus_a = p.copy(); plus_a[ax] += delta * sx
            plus_b = p.copy(); plus_b[az] += delta * sz
            points.extend([tuple(plus_a), tuple(plus_b)])
    for b in recipe["bores"]:
        c = b["center_mm"]; axis = "xyz".index(b["axis"]); perp = [i for i in range(3) if i != axis]
        for a in (-1, 0, 1):
            for sign in (-1, 1):
                for delta in (-1e-4, 0., 1e-4):
                    q = list(c)
                    q[axis] = (x0 + lx * (a + 1) / 2) if axis == 0 else (y0 + ly * (a + 1) / 2) if axis == 1 else lz * (a + 1) / 2
                    q[perp[0]] += sign * (b["radius_mm"] + delta)
                    points.append(tuple(q))
    # One cavity sample in each tube away from the intersection and a sample
    # in material just outside both tubes.
    for bore in (bx, bz):
        axis = "xyz".index(bore["axis"])
        q = list(bore["center_mm"])
        q[axis] = (x0 + lx * .2) if axis == 0 else (y0 + ly * .2) if axis == 1 else lz * .2
        points.append(tuple(q))
    outside_both = cross.copy(); outside_both[ax] += r + .01; outside_both[az] += r + .01
    points.append(tuple(outside_both))
    return points


def bernstein(n, t):
    return [math.comb(n, i) * t**i * (1 - t)**(n - i) for i in range(n + 1)]


def bernstein_d(n, t):
    if n == 0:
        return [0.]
    low = bernstein(n - 1, t)
    return [n * ((low[i - 1] if i else 0.) - (low[i] if i < n else 0.))
            for i in range(n + 1)]


def evaluate(block, uvw):
    degree = block["degrees"]
    net = np.asarray(block["controls"], dtype=float).reshape(*(d + 1 for d in degree), 4)
    b = [bernstein(n, t) for n, t in zip(degree, uvw)]
    db = [bernstein_d(n, t) for n, t in zip(degree, uvw)]
    h = np.zeros(4); dh = np.zeros((3, 4))
    for i, j, k in itertools.product(*(range(d + 1) for d in degree)):
        q = net[i, j, k]
        h += b[0][i] * b[1][j] * b[2][k] * q
        dh[0] += db[0][i] * b[1][j] * b[2][k] * q
        dh[1] += b[0][i] * db[1][j] * b[2][k] * q
        dh[2] += b[0][i] * b[1][j] * db[2][k] * q
    jac = (dh[:, :3] * h[3] - dh[:, 3, None] * h[:3]) / h[3]**2
    return float(np.linalg.det(jac))


def corner_incidence(model):
    vertices = set(); edges = set(); faces = Counter()
    det_min = math.inf
    for b in model["blocks"]:
        d = b["degrees"]
        net = np.asarray(b["controls"], dtype=object).reshape(d[0] + 1, d[1] + 1, d[2] + 1, 4)
        corner = {}
        for i, j, k in itertools.product((0, d[0]), (0, d[1]), (0, d[2])):
            p = net[i, j, k]
            corner[(i > 0, j > 0, k > 0)] = tuple(F(int(p[a]), int(p[3])) for a in range(3))
        vertices.update(corner.values())
        for axis in range(3):
            for bits in itertools.product((False, True), repeat=3):
                if bits[axis]:
                    continue
                other = list(bits); other[axis] = True
                edges.add(tuple(sorted((corner[bits], corner[tuple(other)]))))
            for upper in (False, True):
                faces[tuple(sorted(p for bits, p in corner.items() if bits[axis] == upper))] += 1
        det_min = min(det_min, *(evaluate(b, xyz) for xyz in itertools.product((0., .5, 1.), repeat=3)))
    return {"vertices": len(vertices), "edges": len(edges), "faces": len(faces),
            "cells": len(model["blocks"]), "euler": len(vertices) - len(edges) + len(faces) - len(model["blocks"]),
            "max_face_multiplicity": max(faces.values(), default=0),
            "shared_faces": sum(n == 2 for n in faces.values()),
            "exterior_faces": sum(n == 1 for n in faces.values()),
            "determinant_sample_min_mm3": det_min}


def exact_volume_matches(native, expected_rational, expected_pi):
    """Compare exact rational components while accepting unreduced pairs."""
    try:
        if "rational_part_mm3" in native:
            a, b = native["rational_part_mm3"]
            rational = F(int(a), int(b))
        else:
            rational = F(int(native["twice_rational_part_mm3"]), 2)
        p, q = native["pi_coefficient_mm3"]
        return rational == expected_rational and F(int(p), int(q)) == expected_pi
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return False


def brep_boundary_topology(shape):
    """Count edge uses and vertex links on the imported STEP boundary."""
    edges = list(shape.Edges)
    vertices = list(shape.Vertexes)
    edge_faces = [[] for _ in edges]
    face_vertices = [set() for _ in vertices]
    vertex_edges = [[] for _ in vertices]
    for fi, face in enumerate(shape.Faces):
        for edge in face.Edges:
            ei = next((i for i, known in enumerate(edges) if edge.isSame(known)), None)
            if ei is not None:
                edge_faces[ei].append(fi)
        for vertex in face.Vertexes:
            vi = next((i for i, known in enumerate(vertices) if vertex.isSame(known)), None)
            if vi is not None:
                face_vertices[vi].add(fi)
    for ei, edge in enumerate(edges):
        for vertex in edge.Vertexes:
            vi = next((i for i, known in enumerate(vertices) if vertex.isSame(known)), None)
            if vi is not None:
                vertex_edges[vi].append(ei)
    bad_edge_uses = [i for i, users in enumerate(edge_faces) if len(users) != 2]
    bad_links = []
    for vi, incident_faces in enumerate(face_vertices):
        if not incident_faces:
            bad_links.append(vi); continue
        adjacency = {fi: [] for fi in incident_faces}
        for ei in vertex_edges[vi]:
            users = edge_faces[ei]
            if len(users) == 2 and users[0] in adjacency and users[1] in adjacency:
                adjacency[users[0]].append(users[1]); adjacency[users[1]].append(users[0])
        pending = [next(iter(incident_faces))]; reached = set()
        while pending:
            fi = pending.pop()
            if fi not in reached:
                reached.add(fi); pending.extend(adjacency[fi])
        if reached != incident_faces or any(len(adjacency[fi]) != 2 for fi in incident_faces):
            bad_links.append(vi)
    return {"vertex_count": len(vertices), "edge_count": len(edges), "face_count": len(shape.Faces),
            "euler_characteristic": len(vertices) - len(edges) + len(shape.Faces),
            "closed_shell_count": sum(bool(s.isClosed()) for s in shape.Shells),
            "edge_use_histogram": dict(Counter(len(u) for u in edge_faces)),
            "nonmanifold_edge_count": len(bad_edge_uses),
            "nonmanifold_vertex_link_count": len(bad_links),
            "nonmanifold_vertex_indices": bad_links[:8]}


def compare_model_and_step(cli, recipe_path, model_path, work, recipe, expected, mass_helper):
    started = time.monotonic()
    step = work / f"{recipe_path.stem}.step"
    exp = subprocess.run([str(cli), "export-step", str(model_path), str(step)], text=True, capture_output=True)
    export_result = json.loads(exp.stdout)
    if exp.returncode:
        return {"export": export_result, "qualified": False}
    native = Part.read(str(step))
    reference = reference_shape(recipe)
    shapes = [reference, native]
    native_boundary = Part.makeCompound(list(native.Shells))
    boundary_topology = brep_boundary_topology(native)
    masses = []
    for i, shape in enumerate(shapes):
        path = work / f"{recipe_path.stem}-{i}.brep"
        shape.exportBrep(str(path))
        measure = json.loads(subprocess.check_output([str(mass_helper), str(path), "--centroid"], text=True))[-1]
        masses.append({"volume_mm3": measure["volume_mm3"], "centroid_mm": measure["centroid_mm"],
                       "valid": bool(shape.isValid()), "solids": len(shape.Solids), "shells": len(shape.Shells),
                       "faces": len(shape.Faces), "edges": len(shape.Edges), "vertices": len(shape.Vertexes),
                       "closed_shells": [bool(s.isClosed()) for s in shape.Shells]})
    expected_volume = expected["volume_mm3"]
    expected_center = expected["centroid_mm"]
    volume_rel = max(abs(m["volume_mm3"] / expected_volume - 1) for m in masses)
    centroid_err = max(float(np.linalg.norm(np.asarray(m["centroid_mm"]) - expected_center)) for m in masses)
    differences = {"reference_minus_step_mm3": float(reference.cut(native).Volume),
                   "step_minus_reference_mm3": float(native.cut(reference).Volume)}
    probes = probe_points(recipe)
    mismatch = []; states = Counter()
    for xyz in probes:
        expected_class = expected_state(recipe, xyz)
        point = App.Vector(*xyz)
        inside = bool(native.isInside(point, 1e-8, False))
        dist = None
        if expected_class == "ON":
            # Distance to a solid is zero for every interior point; classify
            # ON against its explicit boundary shell compound instead.
            dist = Part.Vertex(point).distToShape(native_boundary)[0]
            observed_state = "ON" if dist <= BOUNDARY_TOL_MM else ("IN" if inside else "OUT")
        else:
            observed_state = "IN" if inside else "OUT"
        states[expected_class] += 1
        if observed_state != expected_class:
            if dist is None:
                dist = Part.Vertex(point).distToShape(native_boundary)[0]
            mismatch.append({"xyz_mm": list(xyz), "expected": expected_class, "observed": observed_state,
                             "distance_to_step_mm": float(dist)})
    row = {"export": export_result, "occ_masses": masses,
           "analytic_volume_mm3": expected_volume, "analytic_centroid_mm": expected_center.tolist(),
           "max_volume_relative_error": volume_rel, "max_centroid_error_mm": centroid_err,
           "reference_step_boolean_difference_mm3": differences,
           "classification_count": len(probes), "classification_states": dict(states),
           "classification_mismatches": len(mismatch), "first_mismatches": mismatch[:8],
           "native_shell_count": len(native.Shells), "native_single_valid_solid": bool(native.isValid() and len(native.Solids) == 1),
           "native_all_shells_closed": all(s.isClosed() for s in native.Shells),
           "boundary_brep_topology": boundary_topology,
           "step_validation_s": time.monotonic() - started}
    row["qualified"] = (volume_rel < 1e-8 and centroid_err < 1e-7 and not mismatch
                        and max(abs(v) for v in differences.values()) < 1e-7
                        and row["native_single_valid_solid"] and row["native_all_shells_closed"])
    return row


def preserve_initial_refusal(old_cli, work, crossing_recipe, baseline_recipe):
    if old_cli is None:
        return {"status": "not_run", "reason": "no baseline CLI supplied"}
    folder = work / "initial_refusal"; folder.mkdir(parents=True, exist_ok=True)
    model = folder / "previous.json"; base_recipe = folder / "base-recipe.json"
    candidate = folder / "crossing-recipe.json"
    base_recipe.write_text(json.dumps(baseline_recipe))
    p0 = subprocess.run([str(old_cli), "build", str(base_recipe), str(model)], text=True, capture_output=True)
    if p0.returncode:
        return {"status": "baseline_setup_failed", "stdout": p0.stdout, "stderr": p0.stderr}
    before = model.read_bytes()
    candidate.write_text(json.dumps(crossing_recipe))
    p1 = subprocess.run([str(old_cli), "build", str(candidate), str(model)], text=True, capture_output=True)
    try:
        result = json.loads(p1.stdout)
    except json.JSONDecodeError:
        result = {"stdout": p1.stdout, "stderr": p1.stderr}
    preserved = model.read_bytes() == before
    return {"baseline_cli_sha256": hashlib.sha256(old_cli.read_bytes()).hexdigest(),
            "old_crossing_returncode": p1.returncode, "old_crossing_diagnostic": result,
            "previous_model_preserved": preserved,
            "recorded_refusal": p1.returncode != 0 and preserved}


def run(cli: Path, work: Path, baseline_cli: Path | None = None):
    started = time.monotonic(); cli_sha = hashlib.sha256(cli.read_bytes()).hexdigest()
    work.mkdir(parents=True, exist_ok=True); mass_helper = compile_mass_oracle(work)
    cases = json.loads(CORPUS.read_text()); rows = []; model = work / "model.json"
    recipe_path = work / "recipe.json"; prior = None
    baseline_refusal = preserve_initial_refusal(baseline_cli, work, cases[0]["recipe"], cases[1]["recipe"])
    regression = json.loads((HERE.parent / "regressions/steinmetz_hidden_bore_surface.json").read_text())
    regression_recipe = next(x["recipe"] for x in cases if x["id"] == regression["case"])
    regression_actual = expected_state(regression_recipe, regression["point_mm"])
    regression_check = {"fixture": regression["case"], "point_mm": regression["point_mm"],
                        "expected": regression["expected_state"], "observed": regression_actual,
                        "passed": regression_actual == regression["expected_state"]}
    for item in cases:
        recipe = item["recipe"]; recipe_path.write_text(json.dumps(recipe, allow_nan=False))
        previous = model.read_bytes() if model.exists() else None
        t0 = time.monotonic()
        proc = subprocess.run([str(cli), "build", str(recipe_path), str(model)], text=True, capture_output=True)
        build_s = time.monotonic() - t0
        try:
            diag = json.loads(proc.stdout)
        except json.JSONDecodeError:
            diag = {"stdout": proc.stdout, "stderr": proc.stderr}
        row = {"id": item["id"], "recipe": recipe, "expected_accepted": item["accepted"],
               "build_audit": diag, "build_audit_s": build_s,
               "qualified": (proc.returncode == 0) == item["accepted"]}
        if proc.returncode:
            row["previous_model_preserved"] = (model.read_bytes() if model.exists() else None) == previous
            row["qualified"] = not item["accepted"] and row["previous_model_preserved"]
            rows.append(row); continue
        if not item["accepted"]:
            row["qualified"] = False; rows.append(row); continue
        data = json.loads(model.read_text())
        measure = bore_terms(recipe)
        exact = diag.get("volume_exact")
        expected_exact = {"rational_part_mm3": [measure["volume_rational_mm3"].numerator,
                                                   measure["volume_rational_mm3"].denominator],
                          "pi_coefficient_mm3": [measure["pi_coefficient_mm3"], 1]}
        exact_match = exact_volume_matches(exact, measure["volume_rational_mm3"],
                                           F(measure["pi_coefficient_mm3"], 1))
        incidence = corner_incidence(data)
        feature_tags = Counter(f["feature"] for f in data["design_faces"])
        expected_wall_tags = [f"bore/{b['id']}/wall" for b in recipe["bores"]]
        tags_ok = all(feature_tags[tag] > 0 for tag in expected_wall_tags)
        bore_boxes = data.get("partition_plan", {}).get("feature_guards", [])
        row.update(volume_formula={"rational_part_mm3": [measure["volume_rational_mm3"].numerator,
                                                          measure["volume_rational_mm3"].denominator],
                                   "pi_coefficient_mm3": [measure["pi_coefficient_mm3"], 1]},
                   expected_volume_mm3=measure["volume_mm3"], expected_centroid_mm=measure["centroid_mm"].tolist(),
                   volume_exact_expected=expected_exact, volume_exact_native=exact, volume_exact_match=exact_match,
                   topology_oracle=incidence, bore_feature_tags={k: v for k, v in feature_tags.items() if k.startswith("bore/")},
                   feature_guards=bore_boxes)
        row["qualified"] &= bool(diag.get("method") == "certified_orthogonal_sketch_features" and exact_match)
        row["qualified"] &= incidence["max_face_multiplicity"] == 2 and incidence["shared_faces"] == diag.get("shared_faces")
        expected_euler = -2 if measure["crossing"] else -1
        row["expected_euler"] = expected_euler
        row["qualified"] &= incidence["exterior_faces"] == diag.get("exterior_sides") and incidence["euler"] == expected_euler
        row["qualified"] &= incidence["determinant_sample_min_mm3"] > 0 and tags_ok
        row["step"] = compare_model_and_step(cli, recipe_path, model, work, recipe, measure, mass_helper)
        row["qualified"] &= row["step"]["qualified"]
        row["qualified"] = bool(row["qualified"])
        rows.append(row); print(item["id"], row["qualified"], flush=True)
    if hashlib.sha256(cli.read_bytes()).hexdigest() != cli_sha:
        raise RuntimeError("CLI changed during the oracle run")
    return {"schema": "steinmetz_intersecting_bore_independent_oracle_v1",
            "provenance": {"cli": str(cli), "cli_sha256": cli_sha,
                           "baseline_cli": str(baseline_cli) if baseline_cli else None,
                           "checkout_revision": subprocess.check_output(["rtk", "git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                           "python": platform.python_version(), "occ": Part.OCC_VERSION,
                           "seed": SEED, "elapsed_s": time.monotonic() - started},
            "protocol": {"reference": "exact box/cylinder volumes plus 16r^3/3 Steinmetz overlap; independent OCC booleans",
                         "classification": "IN/OUT/ON with analytic cylinder predicates and OCC distance to STEP B-rep",
            "limits": ["equal-radius X/Z intersections and separated axes have closed forms in this oracle",
                                    "sampled determinant values are diagnostic; continuous positivity is audited natively"]},
            "initial_refusal": baseline_refusal, "classification_regression": regression_check,
            "case_count": len(rows), "cases": rows,
            "qualified": baseline_refusal.get("recorded_refusal", True) and regression_check["passed"]
            and all(r["qualified"] for r in rows)}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cli", type=Path, required=True)
    p.add_argument("--baseline-cli", type=Path)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args(); report = run(a.cli.resolve(), a.work, a.baseline_cli.resolve() if a.baseline_cli else None)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    raise SystemExit(int(not report["qualified"]))
