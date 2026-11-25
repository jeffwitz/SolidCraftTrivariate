#!/usr/bin/env python3
"""Independent OCC check of Rust-exported STEP boundaries for J1 and G7.

FreeCAD is deliberately imported before all other modules so its OCC runtime is
initialized in the supported order. The script reads STEP with OCC and only
queries the imported shape; it never calls an external repair or healing tool.
"""

import FreeCAD as App
import Part

import argparse
import json
import math
from pathlib import Path
import subprocess
import time

import j1_geometry


ROOT = Path(__file__).resolve().parents[2]
J1_MANIFEST = ROOT / "tests/references/j1_geometry_corpus.json"
WORK = ROOT / "build/step-boundary"
SAMPLE = (0.0, 0.21132486540518713, 0.5, 0.7886751345948129, 1.0)
DISTANCE_TOL_MM = 1e-7
VOLUME_REL_TOL = 1e-8
VOLUME_STRICT_DIAGNOSTIC = 1e-10
G7_VOLUME = 8.0 * (400.0 - 16.0 * math.pi / 3.0 - 2.0 * math.sqrt(3.0))
G7_TOPOLOGY = (8, 0)  # inter-block faces, boundary Euler characteristic


def source_topology(blocks):
    """Count source boundary faces, edges, vertices and interfaces independently.

    Coincident block faces are recognized by a symmetric 3x3 physical sample
    grid. Boundary edge and vertex counts use the corner incidence of the
    resulting exterior face patches, not the STEP importer's topology.
    """
    face_uses = {}
    descriptors = []
    for bi, block in enumerate(blocks):
        for axis in range(3):
            varying = [index for index in range(3) if index != axis]
            for side in (0, 1):
                sample_points = []
                for a in (0.0, 0.5, 1.0):
                    for b in (0.0, 0.5, 1.0):
                        q = [0.0, 0.0, 0.0]
                        q[axis] = float(side)
                        q[varying[0]], q[varying[1]] = a, b
                        point, _ = j1_geometry.evaluate(block, q)
                        sample_points.append(tuple(round(x, 8) for x in point))
                key = tuple(sorted(sample_points))
                descriptor = (bi, axis, side)
                face_uses.setdefault(key, []).append(descriptor)
                descriptors.append((key, descriptor))

    shared = sum(len(uses) == 2 for uses in face_uses.values())
    exterior = [uses[0] for uses in face_uses.values() if len(uses) == 1]
    nonmanifold = sum(len(uses) > 2 for uses in face_uses.values())
    vertices, edges = set(), set()
    for bi, axis, side in exterior:
        block = blocks[bi]
        varying = [index for index in range(3) if index != axis]
        corners = []
        for a, b in ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)):
            q = [0.0, 0.0, 0.0]
            q[axis] = float(side)
            q[varying[0]], q[varying[1]] = a, b
            point, _ = j1_geometry.evaluate(block, q)
            corners.append(tuple(round(x, 8) for x in point))
            vertices.add(corners[-1])
        for i in range(4):
            edges.add(tuple(sorted((corners[i], corners[(i + 1) % 4]))))
    return {
        "interfaces": shared,
        "exterior_faces": len(exterior),
        "nonmanifold_interfaces": nonmanifold,
        "source_faces": len(exterior),
        "source_edges": len(edges),
        "source_vertices": len(vertices),
        "source_boundary_euler": len(vertices) - len(edges) + len(exterior),
    }


def boundary_parameter_samples(block, axis, side):
    varying = [index for index in range(3) if index != axis]
    for a in SAMPLE:
        for b in SAMPLE:
            q = [0.0, 0.0, 0.0]
            q[axis] = float(side)
            q[varying[0]], q[varying[1]] = a, b
            yield q


def native_samples(model, interfaces):
    blocks = model["blocks"]
    face_use_count = {}
    # The same face signatures used by source_topology identify which block
    # sides are internal. Keep the full signature->uses table to omit those
    # sides from the physical boundary sample set.
    for bi, block in enumerate(blocks):
        for axis in range(3):
            varying = [index for index in range(3) if index != axis]
            for side in (0, 1):
                points = []
                for a in (0.0, 0.5, 1.0):
                    for b in (0.0, 0.5, 1.0):
                        q = [0.0, 0.0, 0.0]
                        q[axis] = float(side)
                        q[varying[0]], q[varying[1]] = a, b
                        point, _ = j1_geometry.evaluate(block, q)
                        points.append(tuple(round(x, 8) for x in point))
                face_use_count.setdefault(tuple(sorted(points)), []).append((bi, axis, side))
    boundary = []
    walls = []
    contours = []
    for uses in face_use_count.values():
        if len(uses) != 1:
            continue
        bi, axis, side = uses[0]
        block = blocks[bi]
        for q in boundary_parameter_samples(block, axis, side):
            point, _ = j1_geometry.evaluate(block, q)
            boundary.append(point)
    by_id = {block["id"]: block for block in blocks}
    design_features = {"circle_r4", "circle_r2", "fillet_r6", "inside_fillet",
                       "hole_x", "hole_y", "bore_left", "bore_right"}
    for design in model.get("design_faces", []):
        if design.get("feature") not in design_features:
            continue
        block = by_id.get(design.get("block"))
        if block is None:
            continue
        side = design.get("side", "")
        if len(side) != 2 or side[0] not in "uvw" or side[1] not in "01":
            continue
        axis = "uvw".index(side[0])
        fixed = float(side[1])
        varying = [index for index in range(3) if index != axis]
        for a in SAMPLE:
            for b in SAMPLE:
                q = [0.0, 0.0, 0.0]
                q[axis] = fixed
                q[varying[0]], q[varying[1]] = a, b
                point, _ = j1_geometry.evaluate(block, q)
                walls.append(point)
        if 2 in varying:
            for a in SAMPLE:
                for w in (0.0, 1.0):
                    q = [0.0, 0.0, 0.0]
                    q[axis] = fixed
                    q[varying[0]] = a
                    q[2] = w
                    point, _ = j1_geometry.evaluate(block, q)
                    contours.append(point)
    return boundary, walls, contours


def point_distance(shape, point):
    vertex = Part.Vertex(App.Vector(point[0], point[1], point[2]))
    return float(shape.distToShape(vertex)[0])


def edge_face_usage(shape):
    usage = []
    for edge in shape.Edges:
        count = 0
        for face in shape.Faces:
            if any(edge.isSame(face_edge) for face_edge in face.Edges):
                count += 1
        usage.append(count)
    return {
        "edge_count": len(usage),
        "minimum_face_usage": min(usage, default=0),
        "maximum_face_usage": max(usage, default=0),
        "edges_used_twice": sum(count == 2 for count in usage),
        "edges_not_used_twice": [index for index, count in enumerate(usage) if count != 2],
    }


def read_step(path):
    shape = Part.read(str(path))
    if shape is None:
        raise RuntimeError(f"OCC returned no shape for {path}")
    return shape


def expected_cases(manifest):
    cases = []
    for case in manifest["cases"]:
        cases.append({"id": case["id"], "recipe": case["recipe"],
                      "expected_blocks": case["expected_blocks"],
                      "analytic_volume": j1_geometry.ANALYTIC_VOLUME[case["id"]],
                      "topology": j1_geometry.TOPOLOGY[case["id"]]})
    cases.append({"id": "overlap_bore", "recipe": {"family": "overlap_bore_v1"},
                  "expected_blocks": 8, "analytic_volume": G7_VOLUME,
                  "topology": G7_TOPOLOGY})
    return cases


def run_case(case, cli, work):
    cid = case["id"]
    recipe_path = work / "recipes" / f"{cid}.json"
    model_path = work / "models" / f"{cid}.json"
    step_path = work / "step" / f"{cid}.step"
    for path in (recipe_path, model_path, step_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    recipe_path.write_text(json.dumps(case["recipe"], indent=2) + "\n")
    start = time.perf_counter_ns()
    build = subprocess.run([cli, "build", str(recipe_path), str(model_path)],
                           cwd=ROOT, capture_output=True, text=True, check=False)
    build_ms = (time.perf_counter_ns() - start) / 1e6
    try:
        build_status = json.loads(build.stdout)
    except json.JSONDecodeError:
        build_status = {"unparsed_stdout": build.stdout.strip()}
    if build.returncode != 0 or not model_path.exists():
        return {"status": "build_refused", "build_returncode": build.returncode,
                "build_wall_ms": build_ms, "build_status": build_status,
                "stderr": build.stderr.strip()}

    start = time.perf_counter_ns()
    export = subprocess.run([cli, "export-step", str(model_path), str(step_path)],
                            cwd=ROOT, capture_output=True, text=True, check=False)
    export_ms = (time.perf_counter_ns() - start) / 1e6
    try:
        export_status = json.loads(export.stdout)
    except json.JSONDecodeError:
        export_status = {"unparsed_stdout": export.stdout.strip()}
    if export.returncode != 0 or not step_path.exists():
        return {"status": "step_export_failed", "build_wall_ms": build_ms,
                "build_status": build_status, "export_wall_ms": export_ms,
                "export_returncode": export.returncode, "export_status": export_status,
                "stderr": export.stderr.strip()}

    model = json.loads(model_path.read_text())
    source = source_topology(model["blocks"])
    expected_interfaces, expected_euler = case["topology"]
    shape = read_step(step_path)
    solids = list(shape.Solids)
    shell_state = [{"closed": bool(shell.isClosed()), "faces": len(shell.Faces),
                    "edges": len(shell.Edges)}
                   for solid in solids for shell in solid.Shells]
    usage = edge_face_usage(shape)
    boundary_points, wall_points, contour_points = native_samples(model, source["interfaces"])
    face_compound = Part.makeCompound(list(shape.Faces))
    edge_compound = Part.makeCompound(list(shape.Edges))
    max_boundary_distance = max((point_distance(face_compound, p) for p in boundary_points), default=None)
    max_wall_distance = max((point_distance(face_compound, p) for p in wall_points), default=None)
    max_contour_distance = max((point_distance(edge_compound, p) for p in contour_points), default=None)
    volume = float(shape.Volume)
    relative_volume_error = abs(volume - case["analytic_volume"]) / case["analytic_volume"]
    source_faces_match = len(shape.Faces) == source["source_faces"]
    source_edges_match = len(shape.Edges) == source["source_edges"]
    interfaces_match = source["interfaces"] == expected_interfaces
    euler_match = source["source_boundary_euler"] == expected_euler
    qualified = (
        len(model["blocks"]) == case["expected_blocks"]
        and len(solids) == 1
        and bool(shape.isValid())
        and bool(solids[0].isClosed())
        and len(shell_state) == 1
        and shell_state[0]["closed"]
        and source_faces_match and source_edges_match
        and interfaces_match and euler_match
        and not usage["edges_not_used_twice"]
        and relative_volume_error <= VOLUME_REL_TOL
        and max_boundary_distance <= DISTANCE_TOL_MM
        and (max_wall_distance is None or max_wall_distance <= DISTANCE_TOL_MM)
        and (max_contour_distance is None or max_contour_distance <= DISTANCE_TOL_MM)
    )
    return {
        "status": "measured",
        "qualified_by_oracle": qualified,
        "recipe": case["recipe"],
        "build_wall_ms": build_ms,
        "build_status": build_status,
        "export_wall_ms": export_ms,
        "export_status": export_status,
        "step_size_bytes": step_path.stat().st_size,
        "occ_shape_type": shape.ShapeType,
        "occ_valid": bool(shape.isValid()),
        "occ_solid_count": len(solids),
        "occ_closed_solid": bool(solids and solids[0].isClosed()),
        "occ_shells": shell_state,
        "occ_faces": len(shape.Faces),
        "occ_edges": len(shape.Edges),
        **usage,
        **source,
        "expected_blocks": case["expected_blocks"],
        "expected_interfaces": expected_interfaces,
        "expected_boundary_euler": expected_euler,
        "source_faces_match_occ": source_faces_match,
        "source_edges_match_occ": source_edges_match,
        "interface_count_matches": interfaces_match,
        "boundary_euler_matches": euler_match,
        "volume_mm3": volume,
        "volume_analytic_mm3": case["analytic_volume"],
        "volume_relative_error": relative_volume_error,
        "volume_meets_strict_1e-10_diagnostic": relative_volume_error <= VOLUME_STRICT_DIAGNOSTIC,
        "native_boundary_sample_count": len(boundary_points),
        "max_native_boundary_distance_mm": max_boundary_distance,
        "native_wall_sample_count": len(wall_points),
        "max_native_wall_distance_mm": max_wall_distance,
        "native_contour_sample_count": len(contour_points),
        "max_native_contour_distance_mm": max_contour_distance,
    }


def atomic_refusal(cli, work, source_model_path):
    """A rejected export must preserve an existing target byte-for-byte."""
    model = json.loads(source_model_path.read_text())
    model["design_faces"][0]["feature"] = "tampered_feature"
    model_path = work / "atomic" / "invalid_model.json"
    target = work / "atomic" / "must_survive.step"
    model_path.parent.mkdir(parents=True, exist_ok=True)
    model_path.write_text(json.dumps(model, indent=2) + "\n")
    sentinel = "preexisting STEP target: keep this byte sequence\n"
    target.write_text(sentinel)
    result = subprocess.run([cli, "export-step", str(model_path), str(target)],
                            cwd=ROOT, capture_output=True, text=True, check=False)
    try:
        status = json.loads(result.stdout)
    except json.JSONDecodeError:
        status = {"unparsed_stdout": result.stdout.strip()}
    survived = target.read_text() == sentinel
    return {
        "status": "measured",
        "export_returncode": result.returncode,
        "export_status": status,
        "target_preserved_exactly": survived,
        "qualified": result.returncode != 0 and survived,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate")
    parser.add_argument("--workdir", type=Path, default=WORK)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--require-qualified", action="store_true")
    args = parser.parse_args()
    cli = args.cli if args.cli.is_absolute() else ROOT / args.cli
    work = args.workdir if args.workdir.is_absolute() else ROOT / args.workdir
    work.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(J1_MANIFEST.read_text())
    cases = expected_cases(manifest)
    try:
        occ_version = getattr(Part, "OCC_VERSION", "unknown")
    except Exception:
        occ_version = "unknown"
    report = {"j1_manifest": manifest["version"], "j1_base": manifest["frozen_base"],
              "occ_version": str(occ_version), "freecad_version": App.Version(),
              "healing_or_repair_called": False,
              "volume_relative_tolerance": VOLUME_REL_TOL,
              "volume_strict_diagnostic_tolerance": VOLUME_STRICT_DIAGNOSTIC,
              "native_distance_tolerance_mm": DISTANCE_TOL_MM,
              "cases": {}}
    for case in cases:
        try:
            result = run_case(case, str(cli), work)
        except Exception as error:
            result = {"status": "oracle_error", "error": f"{type(error).__name__}: {error}"}
        report["cases"][case["id"]] = result
        print(case["id"], result.get("status"), result.get("qualified_by_oracle"), flush=True)
    report["coverage"] = {
        "measured": sum(case.get("status") == "measured" for case in report["cases"].values()),
        "qualified": sum(case.get("qualified_by_oracle", False) for case in report["cases"].values()),
        "required": len(cases),
        "failures": {key: value for key, value in report["cases"].items()
                      if value.get("status") != "measured" or not value.get("qualified_by_oracle")},
    }
    overlap_model = work / "models" / "overlap_bore.json"
    report["atomic_refusal"] = atomic_refusal(str(cli), work, overlap_model) if overlap_model.exists() else {
        "status": "unavailable", "reason": "overlap bore source model missing"
    }
    report_path = args.report or work / "report.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.require_qualified and report["coverage"]["qualified"] != len(cases):
        raise SystemExit(1)
    if args.require_qualified and not report["atomic_refusal"].get("qualified", False):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
