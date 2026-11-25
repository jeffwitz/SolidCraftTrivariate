#!/usr/bin/env python3
"""Independent section, volume, wall and coverage oracle for G10.

The CLI model is treated only as a materialized candidate. The reference
volume, crescent classification and inverse checks are calculated from the
analytic two-circle geometry; quadrature over the emitted rational blocks is
used as an independent measurement.
"""

import argparse
import json
import math
from pathlib import Path
import subprocess
import time

import j1_geometry


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_WORK = ROOT / "build/crescent-pocket"
DEPTH_TICKS = (0, 1, 8, 15, 16)
QUADRATURE_ORDERS = (4, 8, 12, 16, 24, 32)
SAMPLE = (0.0, 0.21132486540518713, 0.5, 0.7886751345948129, 1.0)
SQRT3 = math.sqrt(3.0)
CRESCENT_AREA = 4.0 * math.pi / 3.0 + 2.0 * SQRT3
ANALYTIC_VOLUME = lambda depth: 8.0 * (400.0 - 4.0 * math.pi) - depth * CRESCENT_AREA
VOLUME_REL_TOL = 1e-10
FACE_AREA_REL_TOL = 1e-10
GEOMETRY_TOL = 2e-10


def recipe(depth_ticks):
    return {"family": "crescent_pocket_v1", "depth_half_mm": depth_ticks}


def det2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def norm3(a):
    return math.sqrt(sum(value * value for value in a))


def cross3(a, b):
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]


def source_topology(blocks):
    """Reconstruct face incidence, exterior edge/vertex counts and Euler data."""
    face_uses = {}
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
                face_uses.setdefault(tuple(sorted(points)), []).append((bi, axis, side))
    interfaces = sum(len(uses) == 2 for uses in face_uses.values())
    exterior = [uses[0] for uses in face_uses.values() if len(uses) == 1]
    nonmanifold = sum(len(uses) > 2 for uses in face_uses.values())
    vertices, edges = set(), set()
    for bi, axis, side in exterior:
        block = blocks[bi]
        varying = [index for index in range(3) if index != axis]
        corners = []
        for a, b in ((0., 0.), (1., 0.), (1., 1.), (0., 1.)):
            q = [0.0, 0.0, 0.0]
            q[axis] = float(side)
            q[varying[0]], q[varying[1]] = a, b
            point, _ = j1_geometry.evaluate(block, q)
            corner = tuple(round(x, 8) for x in point)
            vertices.add(corner)
            corners.append(corner)
        for i in range(4):
            edges.add(tuple(sorted((corners[i], corners[(i + 1) % 4]))))
    return {
        "interfaces": interfaces,
        "exterior_faces": len(exterior),
        "nonmanifold_faces": nonmanifold,
        "source_edges": len(edges),
        "source_vertices": len(vertices),
        "boundary_euler": len(vertices) - len(edges) + len(exterior),
    }


def face_axes(side):
    if len(side) != 2 or side[0] not in "uvw" or side[1] not in "01":
        raise ValueError(f"invalid design face side {side!r}")
    axis = "uvw".index(side[0])
    fixed = float(side[1])
    return axis, fixed, [index for index in range(3) if index != axis]


def face_area(block, side, order):
    axis, fixed, varying = face_axes(side)
    nodes, weights = j1_geometry.gauss(order)
    result = 0.0
    for i, a in enumerate(nodes):
        for j, b in enumerate(nodes):
            q = [0.0, 0.0, 0.0]
            q[axis] = fixed
            q[varying[0]], q[varying[1]] = a, b
            _, jac = j1_geometry.evaluate(block, q)
            result += weights[i] * weights[j] * norm3(cross3(jac[varying[0]], jac[varying[1]]))
    return result


def wall_layer(block_id):
    if "_lower_" in block_id or block_id.startswith("g10_crescent_"):
        return "lower"
    if "_upper_" in block_id:
        return "upper"
    return "unknown"


def design_face_areas(model, order):
    blocks = {block["id"]: block for block in model["blocks"]}
    result = {}
    counts = {}
    errors = []
    for face in model.get("design_faces", []):
        key = (wall_layer(face["block"]), face["feature"])
        counts[key] = counts.get(key, 0) + 1
        block = blocks.get(face["block"])
        if block is None:
            errors.append(f"missing design-face block {face['block']}")
            continue
        try:
            area = face_area(block, face["side"], order)
        except (KeyError, ValueError, ZeroDivisionError) as error:
            errors.append(f"{face.get('block')}:{face.get('side')}: {error}")
            continue
        result[key] = result.get(key, 0.0) + area
    return result, counts, errors


def analytic_material(point, depth, block_id=None):
    x, y, z = point
    if x < -10.0 - GEOMETRY_TOL or x > 10.0 + GEOMETRY_TOL:
        return False, "outside_square_x"
    if y < -10.0 - GEOMETRY_TOL or y > 10.0 + GEOMETRY_TOL:
        return False, "outside_square_y"
    if z < -GEOMETRY_TOL or z > 8.0 + GEOMETRY_TOL:
        return False, "outside_plate_z"
    right_d2 = (x - 1.0) ** 2 + y * y
    left_d2 = (x + 1.0) ** 2 + y * y
    crescent_block = block_id is not None and block_id.startswith("g10_crescent_")
    if right_d2 < 4.0 - GEOMETRY_TOL:
        return False, "inside_right_through_bore"
    if crescent_block:
        if left_d2 > 4.0 + GEOMETRY_TOL:
            return False, "crescent_patch_outside_left_disk"
    else:
        if left_d2 < 4.0 - GEOMETRY_TOL:
            return False, "ring_patch_inside_pocket"
    split_z = 8.0 - depth
    if z > split_z + GEOMETRY_TOL and left_d2 < 4.0 - GEOMETRY_TOL:
        return False, "inside_upper_pocket_void"
    return True, "material"


def classify_materialized(blocks, depth):
    grid = (0.1127016653792583, 0.5, 0.8872983346207417)
    failures = []
    positive_determinants = []
    sample_count = 0
    for block in blocks:
        for u in grid:
            for v in grid:
                for w in grid:
                    point, jac = j1_geometry.evaluate(block, (u, v, w))
                    determinant = j1_geometry.determinant(jac)
                    positive_determinants.append(determinant)
                    sample_count += 1
                    valid, reason = analytic_material(point, depth, block["id"])
                    if not valid:
                        failures.append({"block": block["id"], "uvw": [u, v, w],
                                         "point_mm": point, "reason": reason})
    return {
        "sample_count": sample_count,
        "classification_failures": failures,
        "minimum_sampled_determinant_mm3": min(positive_determinants, default=None),
        "nonpositive_determinant_samples": sum(value <= 0 for value in positive_determinants),
    }


def solve2(a, b):
    determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    if abs(determinant) < 1e-18:
        return None
    return [
        (b[0] * a[1][1] - a[0][1] * b[1]) / determinant,
        (a[0][0] * b[1] - b[0] * a[1][0]) / determinant,
    ]


def inverse_crescent_patch(block, target):
    """Newton inverse of the independently sampled planar rational patch."""
    starts = ((0.5, 0.5), (0.2, 0.2), (0.2, 0.8), (0.8, 0.2), (0.8, 0.8))
    for start in starts:
        u, v = start
        for _ in range(40):
            point, jac = j1_geometry.evaluate(block, (u, v, 0.5))
            residual = [point[0] - target[0], point[1] - target[1]]
            error = math.hypot(*residual)
            if error < 2e-10:
                if -1e-8 <= u <= 1.0 + 1e-8 and -1e-8 <= v <= 1.0 + 1e-8:
                    return [min(1.0, max(0.0, u)), min(1.0, max(0.0, v))]
                break
            matrix = [[jac[0][0], jac[1][0]], [jac[0][1], jac[1][1]]]
            delta = solve2(matrix, residual)
            if delta is None:
                break
            accepted = False
            scale = 1.0
            for _ in range(12):
                nu, nv = u - scale * delta[0], v - scale * delta[1]
                if not (-0.5 <= nu <= 1.5 and -0.5 <= nv <= 1.5):
                    scale *= 0.5
                    continue
                next_point, _ = j1_geometry.evaluate(block, (nu, nv, 0.5))
                next_error = math.hypot(next_point[0] - target[0], next_point[1] - target[1])
                if next_error < error:
                    u, v, accepted = nu, nv, True
                    break
                scale *= 0.5
            if not accepted:
                break
    return None


def crescent_points():
    points = []
    nx, ny = 29, 29
    for i in range(1, nx):
        x = -3.0 + 4.0 * i / nx
        for j in range(1, ny):
            y = -2.0 + 4.0 * j / ny
            left_d2 = (x + 1.0) ** 2 + y * y
            right_d2 = (x - 1.0) ** 2 + y * y
            if left_d2 < 4.0 - 0.04 and right_d2 > 4.0 + 0.04 and abs(y) > 1e-6:
                points.append((x, y))
    return points


def inverse_crescent_coverage(blocks, depth):
    crescents = [block for block in blocks if block["id"].startswith("g10_crescent_")]
    points = crescent_points()
    if depth == 8.0:
        return {
            "status": "through_bores_remove_crescent",
            "analytic_samples": len(points),
            "crescent_blocks": len(crescents),
            "unexpected_crescent_material": bool(crescents),
        }
    failures = []
    hit_histogram = {}
    for x, y in points:
        expected_side = "upper" if y > 0 else "lower"
        side_blocks = [block for block in crescents if f"crescent_{expected_side}" in block["id"]]
        hits = []
        for block in side_blocks:
            params = inverse_crescent_patch(block, (x, y))
            if params is not None:
                hits.append((block["id"], params))
        hit_histogram[str(len(hits))] = hit_histogram.get(str(len(hits)), 0) + 1
        if len(hits) != 1:
            failures.append({"point_mm": [x, y], "expected_side": expected_side,
                             "hit_count": len(hits),
                             "reason": "crescent patches do not cover the point exactly once"})
    return {
        "status": "measured",
        "analytic_samples": len(points),
        "crescent_blocks": len(crescents),
        "inverse_hit_histogram": hit_histogram,
        "uncovered_samples": failures,
    }


def geometry_counts(model, depth_ticks):
    depth = depth_ticks * 0.5
    blocks = model.get("blocks", [])
    crescent_count = sum(block["id"].startswith("g10_crescent_") for block in blocks)
    expected_blocks = 8 if depth_ticks == 16 else 14 if depth_ticks == 0 else 22
    expected_interfaces = 8 if depth_ticks == 16 else 20 if depth_ticks == 0 else 36
    topology = source_topology(blocks)
    design_faces = model.get("design_faces", [])
    face_keys = [(face.get("block"), face.get("side")) for face in design_faces]
    block_ids = [block.get("id") for block in blocks]
    floors = [face for face in design_faces if face.get("feature") == "pocket_floor"]
    pocket_walls = [face for face in design_faces if face.get("feature") == "pocket_wall"]
    right_walls = [face for face in design_faces if face.get("feature") in ("bore_right", "right_bore_wall")]
    expected_floor_count = 6 if 0 < depth_ticks < 16 else 0
    return {
        "status": "measured",
        "depth_mm": depth,
        "depth_half_mm": depth_ticks,
        "blocks": len(blocks),
        "expected_blocks": expected_blocks,
        "block_count_matches": len(blocks) == expected_blocks,
        "crescent_blocks": crescent_count,
        "expected_interfaces": expected_interfaces,
        "source_interfaces": topology["interfaces"],
        "source_exterior_faces": topology["exterior_faces"],
        "source_nonmanifold_faces": topology["nonmanifold_faces"],
        "source_boundary_euler": topology["boundary_euler"],
        "interfaces_match": topology["interfaces"] == expected_interfaces,
        "euler_matches": topology["boundary_euler"] == 0,
        "design_face_count": len(design_faces),
        "design_face_identity_unique": len(face_keys) == len(set(face_keys)),
        "floor_face_count": len(floors),
        "expected_floor_face_count": expected_floor_count,
        "floor_selection_matches_depth": len(floors) == expected_floor_count,
        "pocket_wall_face_count": len(pocket_walls),
        "right_bore_wall_face_count": len(right_walls),
        "block_ids": block_ids,
    }


def layer_name(block_id):
    if "g10_crescent_" in block_id or "g10_ring_lower_" in block_id:
        return "lower"
    if "g10_ring_upper_" in block_id:
        return "upper"
    return "unknown"


def surface_measures(model, order):
    blocks = {block["id"]: block for block in model.get("blocks", [])}
    measures, counts = {}, {}
    diagnostics = []
    for face in model.get("design_faces", []):
        feature = face.get("feature", "")
        key = (layer_name(face.get("block", "")), feature)
        counts[key] = counts.get(key, 0) + 1
        block = blocks.get(face.get("block"))
        if block is None:
            diagnostics.append(f"missing block for design face {face}")
            continue
        measures[key] = measures.get(key, 0.0) + face_area(block, face["side"], order)
    return measures, counts, diagnostics


def volume_orders(blocks):
    return {str(order): j1_geometry.volume(blocks, order) for order in QUADRATURE_ORDERS}


def classify_and_measure(model, depth_ticks):
    blocks = model.get("blocks", [])
    depth = depth_ticks * 0.5
    volume_ref = ANALYTIC_VOLUME(depth)
    volume_by_order = volume_orders(blocks)
    surface_by_order = {}
    surface_counts = {}
    surface_diagnostics = []
    for order in QUADRATURE_ORDERS:
        areas, counts, errors = surface_measures(model, order)
        surface_by_order[str(order)] = {f"{layer}:{feature}": value
                                        for (layer, feature), value in sorted(areas.items())}
        if order == QUADRATURE_ORDERS[-1]:
            surface_counts = {f"{layer}:{feature}": value
                              for (layer, feature), value in sorted(counts.items())}
            surface_diagnostics = errors

    area_order16 = {tuple(key.split(":", 1)): value
                    for key, value in surface_by_order[str(QUADRATURE_ORDERS[-1])].items()}
    lower_right = area_order16.get(("lower", "bore_right"), 0.0) + area_order16.get(("lower", "right_bore_wall"), 0.0)
    upper_right = area_order16.get(("upper", "bore_right"), 0.0) + area_order16.get(("upper", "right_bore_wall"), 0.0)
    upper_left = area_order16.get(("upper", "pocket_wall"), 0.0)
    floor = sum(value for (layer, feature), value in area_order16.items()
                if feature == "pocket_floor")
    expected_lower_right = 4.0 * math.pi * (8.0 - depth)
    expected_upper_right = 8.0 * math.pi / 3.0 * depth
    expected_upper_left = 8.0 * math.pi / 3.0 * depth
    expected_floor = CRESCENT_AREA if 0.0 < depth < 8.0 else 0.0

    classification = classify_materialized(blocks, depth)
    inverse = inverse_crescent_coverage(blocks, depth)
    topology = geometry_counts(model, depth_ticks)
    final_order, previous_order = QUADRATURE_ORDERS[-1], QUADRATURE_ORDERS[-2]
    volume_error = abs(volume_by_order[str(final_order)] - volume_ref)
    rel_volume_error = volume_error / volume_ref
    floor_error = abs(floor - expected_floor)
    rel_floor_error = floor_error / expected_floor if expected_floor else floor_error
    relative_wall_errors = {
        "lower_right_bore": abs(lower_right - expected_lower_right) / max(expected_lower_right, 1.0),
        "upper_right_bore": abs(upper_right - expected_upper_right) / max(expected_upper_right, 1.0),
        "upper_pocket_wall": abs(upper_left - expected_upper_left) / max(expected_upper_left, 1.0),
    }
    wall_quadrature_delta = {}
    for name, (layer, feature) in {
        "lower_right_bore": ("lower", "bore_right"),
        "upper_right_bore": ("upper", "bore_right"),
        "upper_pocket_wall": ("upper", "pocket_wall"),
    }.items():
        values = [surface_by_order[str(order)].get(f"{layer}:{feature}", 0.0)
                  for order in (previous_order, final_order)]
        wall_quadrature_delta[name] = abs(values[1] - values[0])
    volume_quad_delta = abs(volume_by_order[str(final_order)] -
                            volume_by_order[str(previous_order)])
    floor_quad_delta = abs(
        sum(value for key, value in surface_by_order[str(final_order)].items()
            if key.endswith(":pocket_floor"))
        - sum(value for key, value in surface_by_order[str(previous_order)].items()
              if key.endswith(":pocket_floor"))
    )
    surface_quad_converged = max(
        [floor_quad_delta, *wall_quadrature_delta.values()], default=0.0
    ) <= FACE_AREA_REL_TOL * max(CRESCENT_AREA, 1.0)
    volume_quad_converged = volume_quad_delta <= VOLUME_REL_TOL * volume_ref
    material_failures = classification["classification_failures"]
    inverse_failures = inverse.get("uncovered_samples", [])
    qualified = (
        topology["block_count_matches"]
        and topology["interfaces_match"] and topology["euler_matches"]
        and topology["source_nonmanifold_faces"] == 0
        and topology["design_face_identity_unique"]
        and topology["floor_selection_matches_depth"]
        and rel_volume_error <= VOLUME_REL_TOL
        and volume_quad_converged
        and surface_quad_converged
        and (not expected_floor or rel_floor_error <= FACE_AREA_REL_TOL)
        and max(relative_wall_errors.values(), default=0.0) <= FACE_AREA_REL_TOL
        and not surface_diagnostics
        and not material_failures
        and classification["nonpositive_determinant_samples"] == 0
        and not inverse_failures
    )
    return {
        **topology,
        "qualified_by_oracle": qualified,
        "volume_by_gauss_order_mm3": volume_by_order,
        "volume_analytic_mm3": volume_ref,
        "final_quadrature_order": final_order,
        "volume_error_final_order_mm3": volume_error,
        "volume_relative_error_final_order": rel_volume_error,
        "volume_quadrature_delta_previous_final_mm3": volume_quad_delta,
        "volume_quadrature_converged": volume_quad_converged,
        "crescent_area_analytic_mm2": CRESCENT_AREA,
        "floor_area_by_gauss_order_mm2": {
            order: sum(value for key, value in data.items() if key.endswith(":pocket_floor"))
            for order, data in surface_by_order.items()
        },
        "floor_area_final_order_mm2": floor,
        "floor_area_expected_mm2": expected_floor,
        "floor_area_error_mm2": floor_error,
        "floor_area_relative_error": rel_floor_error,
        "wall_area_by_layer_feature_final_order_mm2": surface_by_order[str(final_order)],
        "wall_face_counts_by_layer_feature": surface_counts,
        "expected_wall_areas_mm2": {
            "lower_right_bore": expected_lower_right,
            "upper_right_bore": expected_upper_right,
            "upper_pocket_wall": expected_upper_left,
        },
        "measured_wall_areas_mm2": {
            "lower_right_bore": lower_right,
            "upper_right_bore": upper_right,
            "upper_pocket_wall": upper_left,
        },
        "wall_area_relative_errors": relative_wall_errors,
        "wall_quadrature_deltas_mm2": wall_quadrature_delta,
        "floor_quadrature_delta_previous_final_mm2": floor_quad_delta,
        "surface_quadrature_converged": surface_quad_converged,
        "surface_quadrature_orders": surface_by_order,
        "surface_measure_diagnostics": surface_diagnostics,
        "materialized_sample_classification": classification,
        "inverse_crescent_coverage": inverse,
    }


def run_case(cli, work, depth_ticks):
    depth = depth_ticks * 0.5
    name = f"depth_{depth_ticks:02d}"
    recipe_path = work / "recipes" / f"{name}.json"
    model_path = work / "models" / f"{name}.json"
    recipe_path.parent.mkdir(parents=True, exist_ok=True)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    recipe_path.write_text(json.dumps(recipe(depth_ticks), indent=2) + "\n")
    start = time.perf_counter_ns()
    process = subprocess.run([cli, "build", str(recipe_path), str(model_path)],
                             cwd=ROOT, capture_output=True, text=True, check=False)
    elapsed_ms = (time.perf_counter_ns() - start) / 1e6
    try:
        status = json.loads(process.stdout)
    except json.JSONDecodeError:
        status = {"unparsed_stdout": process.stdout.strip()}
    if process.returncode != 0 or not model_path.exists():
        return {"status": "build_refused", "depth_mm": depth,
                "depth_half_mm": depth_ticks, "returncode": process.returncode,
                "build_wall_ms": elapsed_ms, "build_status": status,
                "stderr": process.stderr.strip()}
    model = json.loads(model_path.read_text())
    result = classify_and_measure(model, depth_ticks)
    return {
        "status": "measured",
        "build_wall_ms": elapsed_ms,
        "build_status": status,
        "schema": model.get("schema"),
        **result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate")
    parser.add_argument("--workdir", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--require-qualified", action="store_true")
    args = parser.parse_args()
    cli = args.cli if args.cli.is_absolute() else ROOT / args.cli
    work = args.workdir if args.workdir.is_absolute() else ROOT / args.workdir
    work.mkdir(parents=True, exist_ok=True)
    report = {
        "family": "crescent_pocket_v1",
        "depth_ticks_half_mm": list(DEPTH_TICKS),
        "depths_mm": [tick * 0.5 for tick in DEPTH_TICKS],
        "analytic_crescent_area_mm2": CRESCENT_AREA,
        "analytic_volume": "8*(400-4*pi)-depth*(4*pi/3+2*sqrt(3))",
        "volume_relative_tolerance": VOLUME_REL_TOL,
        "face_area_relative_tolerance": FACE_AREA_REL_TOL,
        "cases": {},
    }
    for depth_ticks in DEPTH_TICKS:
        key = f"depth_{depth_ticks:02d}"
        try:
            case = run_case(str(cli), work, depth_ticks)
        except Exception as error:
            case = {"status": "oracle_error", "error": f"{type(error).__name__}: {error}"}
        report["cases"][key] = case
        print(key, case.get("status"), case.get("qualified_by_oracle"), flush=True)
    report["coverage"] = {
        "measured": sum(case.get("status") == "measured" for case in report["cases"].values()),
        "qualified": sum(case.get("qualified_by_oracle", False) for case in report["cases"].values()),
        "required": len(DEPTH_TICKS),
        "failures": {key: case for key, case in report["cases"].items()
                      if case.get("status") != "measured" or not case.get("qualified_by_oracle")},
    }
    output = args.report or work / "oracle.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.require_qualified and report["coverage"]["qualified"] != len(DEPTH_TICKS):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
