#!/usr/bin/env python3
"""Independent physical and OCC oracle for the translating-bore family.

The physical reference is the analytic difference of a plate, a blind
rectangular pocket, and a circular cylinder clipped to the plate, with the
bore/pocket overlap restored. Optional CLI models and STEP files are treated
as candidates only. FreeCAD is imported first so the intended OCC runtime is
initialized before the rest of the oracle.
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
DEFAULT_WORK = ROOT / "build/travelling-bore"
RADIUS_UM = 2000
CENTERS_UM = (
    3000, 1000, 0, -80, -400, -800, -2000, -3200, -3600, -4000,
    -3600, -3200, -2000, -800, -400, -80, 0,
    3700, 8000, 8080, 8400, 8800, 10000, 11200, 11600, 11920, 12000, 13000,
)
BOUNDED_EXTRA_CENTERS_UM = (
    -100, 11900,
    -4001, -3999, -1, 1, 7999, 8001, 11999, 12001,
    -2001, -1999, 9999, 10001,
    -5000, 14000,
)
BOUNDED_CENTERS_UM = CENTERS_UM + BOUNDED_EXTRA_CENTERS_UM
TANGENCIES_UM = (-4000, 0, 8000, 12000)
LINE_FEATURES = (
    ("plate_x_minus_10", -10000),
    ("pocket_x_minus_6", -6000),
    ("pocket_x_minus_2", -2000),
    ("plate_x_plus_10", 10000),
)
GAUSS_ORDERS = (8, 12, 16, 24, 32)
VOLUME_REL_TOL = 1e-9
AREA_REL_TOL = 1e-8
RADIUS_TOL_MM = 1e-12
DISTANCE_TOL_MM = 1e-7
BOUNDED_FIDELITY_MM = 1e-8
MODEL_SAMPLE_NODES = (0.1127016653792583, 0.5, 0.8872983346207417)
SECTION_GRID = tuple(-9.0 + 3.0 * i for i in range(7))


def circle_primitive(x_mm):
    """Primitive of a disk's vertical chord width, in square mm."""
    radius = 2.0
    x = max(-radius, min(radius, x_mm))
    return x * math.sqrt(max(0.0, radius * radius - x * x)) + radius * radius * math.asin(x / radius)


def disk_strip_area(center_mm, left_mm, right_mm):
    left = max(-2.0, min(2.0, left_mm - center_mm))
    right = max(-2.0, min(2.0, right_mm - center_mm))
    return max(0.0, circle_primitive(right) - circle_primitive(left))


def chord_width(relative_x_mm):
    if abs(relative_x_mm) >= 2.0:
        return 0.0
    return 2.0 * math.sqrt(4.0 - relative_x_mm * relative_x_mm)


def disk_strip_derivative(center_mm, left_mm, right_mm):
    return chord_width(left_mm - center_mm) - chord_width(right_mm - center_mm)


def analytic_volume_mm3(center_mm):
    return (3072.0 - 8.0 * disk_strip_area(center_mm, -10.0, 10.0)
            + 4.0 * disk_strip_area(center_mm, -6.0, -2.0))


def circle_arc_length_in_strip(center_mm, left_mm, right_mm):
    lower = max(-1.0, min(1.0, (left_mm - center_mm) / 2.0))
    upper = max(-1.0, min(1.0, (right_mm - center_mm) / 2.0))
    if upper <= lower:
        return 0.0
    return 4.0 * (math.acos(lower) - math.acos(upper))


def analytic_metrics(center_um):
    center = center_um / 1000.0
    plate_disk = disk_strip_area(center, -10.0, 10.0)
    pocket_overlap = disk_strip_area(center, -6.0, -2.0)
    d_plate_disk = disk_strip_derivative(center, -10.0, 10.0)
    d_pocket_overlap = disk_strip_derivative(center, -6.0, -2.0)
    arc_plate = circle_arc_length_in_strip(center, -10.0, 10.0)
    arc_pocket = circle_arc_length_in_strip(center, -6.0, -2.0)
    derivative = -8.0 * d_plate_disk + 4.0 * d_pocket_overlap
    derivative_step = 1e-5
    derivative_fd = (
        analytic_volume_mm3(center + derivative_step)
        - analytic_volume_mm3(center - derivative_step)
    ) / (2.0 * derivative_step)
    event_mask = tuple(abs(center_um - x) == RADIUS_UM for _, x in LINE_FEATURES)
    if center_um + RADIUS_UM < -6000 or center_um - RADIUS_UM > -2000:
        bore_pocket = "disjoint"
    elif center_um + RADIUS_UM == -6000 or center_um - RADIUS_UM == -2000:
        bore_pocket = "tangent"
    else:
        bore_pocket = "overlapping"
    if center_um + RADIUS_UM < 10000:
        bore_plate = "interior"
    elif center_um + RADIUS_UM == 10000 or center_um - RADIUS_UM == 10000:
        bore_plate = "tangent"
    elif center_um - RADIUS_UM < 10000:
        bore_plate = "open_to_exterior"
    else:
        bore_plate = "disjoint"
    areas = {
        "bore_wall_below_pocket": 4.0 * arc_plate,
        "bore_wall_above_pocket": 4.0 * max(0.0, arc_plate - arc_pocket),
        "pocket_x_min": 4.0 * (8.0 - chord_width(-6.0 - center)),
        "pocket_x_max": 4.0 * (8.0 - chord_width(-2.0 - center)),
        "pocket_y_min": 16.0,
        "pocket_y_max": 16.0,
        "pocket_floor": max(0.0, 32.0 - pocket_overlap),
        "bottom": max(0.0, 400.0 - plate_disk),
        "top": max(0.0, 368.0 - plate_disk + pocket_overlap),
        "x_min": 160.0,
        "x_max": 160.0 - 8.0 * chord_width(10.0 - center),
        "y_min": 160.0,
        "y_max": 160.0,
    }
    return {
        "center_um": center_um,
        "center_mm": center,
        "volume_mm3": analytic_volume_mm3(center),
        "dvolume_dcenter_mm2": derivative,
        "dvolume_centered_fd_mm2": derivative_fd,
        "dvolume_fd_abs_error_mm2": abs(derivative_fd - derivative),
        "areas_by_feature_mm2": areas,
        "total_boundary_area_mm2": sum(areas.values()),
        "bore_pocket_relation": bore_pocket,
        "bore_plate_relation": bore_plate,
        "line_relations": {
            name: "tangent" if tangent else ("crosses" if abs(center_um - x) < RADIUS_UM else "misses")
            for (name, x), tangent in zip(LINE_FEATURES, event_mask)
        },
        "physical_volume_gradient_exists": True,
        "representation_jet_available": not any(event_mask),
        "tangency_events": [name for (name, _), tangent in zip(LINE_FEATURES, event_mask) if tangent],
        "exact_intersections_representable": rational_intersections(center_um),
    }


def rational_intersections(center_um):
    """Whether every transverse arrangement crossing has rational y in mm."""
    for _, line_um in LINE_FEATURES:
        delta = abs(center_um - line_um)
        if delta < RADIUS_UM:
            radicand = RADIUS_UM * RADIUS_UM - delta * delta
            root = math.isqrt(radicand)
            if root * root != radicand:
                return False
    return True


def reference_cad(center_um):
    """Build the independent OCC boolean reference, without shape healing."""
    center = center_um / 1000.0
    plate = Part.makeBox(20.0, 20.0, 8.0, App.Vector(-10.0, -10.0, 0.0))
    bore = Part.makeCylinder(2.0, 8.0, App.Vector(center, 0.0, 0.0), App.Vector(0.0, 0.0, 1.0))
    pocket = Part.makeBox(4.0, 8.0, 4.0, App.Vector(-6.0, -4.0, 4.0))
    return plate.cut(bore).cut(pocket)


def gauss(order):
    nodes, weights = j1_geometry.gauss(order)
    return nodes, weights


def bernstein(degree, nodes):
    import numpy as np
    from math import comb
    return np.array([[comb(degree, i) * x**i * (1.0 - x)**(degree - i)
                      for i in range(degree + 1)] for x in nodes])


def bernstein_derivative(degree, nodes):
    import numpy as np
    if degree == 1:
        return np.tile(np.array([-1.0, 1.0]), (len(nodes), 1))
    low = bernstein(degree - 1, nodes)
    return np.array([
        [degree * ((low[row, i - 1] if i else 0.0)
                   - (low[row, i] if i < degree else 0.0))
         for i in range(degree + 1)] for row in range(len(nodes))
    ])


def evaluate_grid(block, axes):
    """Vectorized rational map and Jacobian on three tensor-product axes."""
    import numpy as np
    degrees = block["degrees"]
    dims = [degree + 1 for degree in degrees]
    controls = np.asarray(block["controls"], dtype=np.float64).reshape((*dims, 4))
    basis = [bernstein(degree, axis) for degree, axis in zip(degrees, axes)]
    dbasis = [bernstein_derivative(degree, axis) for degree, axis in zip(degrees, axes)]
    h = np.einsum("ui,vj,wk,ijkc->uvwc", basis[0], basis[1], basis[2], controls, optimize=True)
    dh = [
        np.einsum("ui,vj,wk,ijkc->uvwc", dbasis[0], basis[1], basis[2], controls, optimize=True),
        np.einsum("ui,vj,wk,ijkc->uvwc", basis[0], dbasis[1], basis[2], controls, optimize=True),
        np.einsum("ui,vj,wk,ijkc->uvwc", basis[0], basis[1], dbasis[2], controls, optimize=True),
    ]
    point = h[..., :3] / h[..., 3, None]
    jac = np.empty((*h.shape[:3], 3, 3), dtype=np.float64)
    for axis in range(3):
        jac[..., axis, :] = (dh[axis][..., :3] - point * dh[axis][..., 3, None]) / h[..., 3, None]
    return point, jac


def block_volume(block, order):
    import numpy as np
    nodes, weights = gauss(order)
    axes = [nodes, nodes, nodes]
    _, jac = evaluate_grid(block, axes)
    determinant = np.linalg.det(np.swapaxes(jac, -1, -2))
    return float(np.einsum("i,j,k,ijk->", weights, weights, weights, determinant, optimize=True))


def model_volume(blocks, order):
    return sum(block_volume(block, order) for block in blocks)


def face_area(block, side, order, z_range=None):
    import numpy as np
    axis_name, side_value = side[0], int(side[1])
    fixed_axis = "uvw".index(axis_name)
    varying = [axis for axis in range(3) if axis != fixed_axis]
    nodes, weights = gauss(order)
    axes = [nodes, nodes, nodes]
    axes[fixed_axis] = [float(side_value)]
    points, jac = evaluate_grid(block, axes)
    tangent_a, tangent_b = jac[..., varying[0], :], jac[..., varying[1], :]
    cross = np.cross(tangent_a, tangent_b)
    density = np.squeeze(np.linalg.norm(cross, axis=-1), axis=fixed_axis)
    measure_weights = np.outer(weights, weights)
    if z_range is not None:
        z = np.squeeze(points[..., 2], axis=fixed_axis)
        mask = (z >= z_range[0] - 1e-10) & (z <= z_range[1] + 1e-10)
        density = np.where(mask, density, 0.0)
    return float(np.sum(density * measure_weights))


def face_feature_areas(model, order):
    blocks = {block["id"]: block for block in model.get("blocks", [])}
    areas, counts, errors = {}, {}, []
    for face in model.get("design_faces", []):
        block = blocks.get(face.get("block"))
        if block is None:
            errors.append(f"unknown block for {face}")
            continue
        feature, side = face.get("feature", ""), face.get("side", "")
        try:
            if feature == "bore_wall":
                for key, zrange in (("bore_wall_below_pocket", (0.0, 4.0)),
                                    ("bore_wall_above_pocket", (4.0, 8.0))):
                    value = face_area(block, side, order, zrange)
                    areas[key] = areas.get(key, 0.0) + value
                    counts[key] = counts.get(key, 0) + int(value > 1e-12)
            else:
                value = face_area(block, side, order)
                key = feature
                areas[key] = areas.get(key, 0.0) + value
                counts[key] = counts.get(key, 0) + 1
        except (KeyError, ValueError, IndexError, ZeroDivisionError) as error:
            errors.append(f"{face}: {error}")
    return areas, counts, errors


def design_face_accuracy(model, center_um):
    errors = {}
    tests = (0.0, 0.21132486540518713, 0.5, 0.7886751345948129, 1.0)
    blocks = {block["id"]: block for block in model.get("blocks", [])}
    for face in model.get("design_faces", []):
        block = blocks.get(face.get("block"))
        if block is None:
            continue
        side = face.get("side", "")
        if len(side) != 2 or side[0] not in "uvw" or side[1] not in "01":
            errors.setdefault("invalid_design_side", []).append(face)
            continue
        axis = "uvw".index(side[0])
        varying = [i for i in range(3) if i != axis]
        feature = face.get("feature", "")
        local = []
        for a in tests:
            for b in tests:
                q = [0.0, 0.0, 0.0]
                q[axis] = float(int(side[1]))
                q[varying[0]], q[varying[1]] = a, b
                point, _ = j1_geometry.evaluate(block, q)
                x, y, z = point
                if feature == "bore_wall":
                    local.append(abs(math.hypot(x - center_um / 1000.0, y) - 2.0))
                elif feature in ("x_min", "x_max", "pocket_x_min", "pocket_x_max"):
                    target = {"x_min": -10.0, "x_max": 10.0,
                              "pocket_x_min": -6.0, "pocket_x_max": -2.0}[feature]
                    local.append(abs(x - target))
                elif feature in ("y_min", "y_max", "pocket_y_min", "pocket_y_max"):
                    target = {"y_min": -10.0, "y_max": 10.0,
                              "pocket_y_min": -4.0, "pocket_y_max": 4.0}[feature]
                    local.append(abs(y - target))
                elif feature in ("bottom", "top", "pocket_floor"):
                    target = {"bottom": 0.0, "top": 8.0, "pocket_floor": 4.0}[feature]
                    local.append(abs(z - target))
        if local:
            errors.setdefault(feature, []).append(max(local))
    return {feature: max(values) for feature, values in errors.items()}


def analytic_material(point):
    x, y, z = point
    if x < -10.0 or x > 10.0 or y < -10.0 or y > 10.0 or z < 0.0 or z > 8.0:
        return False
    if (x - analytic_material.center) ** 2 + y * y < 4.0:
        return False
    if z > 4.0 and -6.0 < x < -2.0 and -4.0 < y < 4.0:
        return False
    return True


def set_material_center(center_um):
    analytic_material.center = center_um / 1000.0


def materialized_samples(blocks, center_um, fidelity_band_mm=0.0):
    set_material_center(center_um)
    failures, boundary_band, determinants = [], [], []
    for block in blocks:
        for u in MODEL_SAMPLE_NODES:
            for v in MODEL_SAMPLE_NODES:
                for w in MODEL_SAMPLE_NODES:
                    point, jac = j1_geometry.evaluate(block, (u, v, w))
                    determinant = j1_geometry.determinant(jac)
                    determinants.append(determinant)
                    if not analytic_material(point):
                        record = {"block": block["id"], "uvw": [u, v, w],
                                  "point_mm": point, "determinant_mm3": determinant,
                                  "ideal_boundary_clearance_mm": ideal_boundary_clearance(point, center_um)}
                        if fidelity_band_mm and record["ideal_boundary_clearance_mm"] <= fidelity_band_mm:
                            boundary_band.append(record)
                        else:
                            failures.append(record)
    return {
        "samples": len(determinants),
        "classification_failures": failures,
        "outside_ideal_within_fidelity_band": boundary_band,
        "fidelity_band_mm": fidelity_band_mm,
        "nonpositive_sampled_determinants": sum(value <= 0.0 for value in determinants),
        "minimum_sampled_determinant_mm3": min(determinants, default=None),
    }


def solve3(matrix, vector):
    import numpy as np
    try:
        return np.linalg.solve(np.asarray(matrix, dtype=float), np.asarray(vector, dtype=float)).tolist()
    except np.linalg.LinAlgError:
        return None


def inverse_block(block, target):
    """Damped Newton inverse in block parameters, using the independent map."""
    starts = ((0.5, 0.5, 0.5), (0.2, 0.2, 0.2), (0.8, 0.2, 0.8),
              (0.2, 0.8, 0.8), (0.8, 0.8, 0.2))
    target = tuple(float(value) for value in target)
    for start in starts:
        uvw = list(start)
        for _ in range(40):
            point, jac = j1_geometry.evaluate(block, uvw)
            residual = [point[i] - target[i] for i in range(3)]
            norm = math.sqrt(sum(value * value for value in residual))
            if norm <= DISTANCE_TOL_MM:
                if all(-1e-8 <= value <= 1.0 + 1e-8 for value in uvw):
                    return [min(1.0, max(0.0, value)) for value in uvw]
                break
            delta = solve3([[jac[axis][component] for axis in range(3)]
                            for component in range(3)], residual)
            if delta is None:
                break
            scale, accepted = 1.0, False
            for _ in range(14):
                candidate = [uvw[i] - scale * delta[i] for i in range(3)]
                if any(value < -0.4 or value > 1.4 for value in candidate):
                    scale *= 0.5
                    continue
                next_point, _ = j1_geometry.evaluate(block, candidate)
                next_norm = math.sqrt(sum((next_point[i] - target[i]) ** 2 for i in range(3)))
                if next_norm < norm:
                    uvw, accepted = candidate, True
                    break
                scale *= 0.5
            if not accepted:
                break
    return None


def physical_bounds(block):
    controls = block["controls"]
    points = [[p[i] / p[3] for i in range(3)] for p in controls]
    return [(min(p[i] for p in points), max(p[i] for p in points)) for i in range(3)]


def ideal_boundary_clearance(point, center_um):
    x, y, z = point
    center = center_um / 1000.0
    return min(
        abs(math.hypot(x - center, y) - 2.0),
        abs(x + 10.0), abs(x - 10.0), abs(y + 10.0), abs(y - 10.0),
        abs(z), abs(z - 8.0),
        abs(x + 6.0), abs(x + 2.0), abs(y + 4.0), abs(y - 4.0), abs(z - 4.0),
    )


def section_samples(blocks, center_um, fidelity_band_mm=0.0):
    """Compare analytic sections at z=2,6 and inverse-map material samples."""
    set_material_center(center_um)
    candidates = [(block, physical_bounds(block)) for block in blocks]
    records, failures = [], []
    for z in (2.0, 6.0):
        count, total_count, missed = 0, 0, []
        for x in SECTION_GRID:
            for y in SECTION_GRID:
                point = (x + 0.137, y + 0.271, z)
                if not analytic_material(point):
                    continue
                total_count += 1
                if fidelity_band_mm and ideal_boundary_clearance(point, center_um) <= fidelity_band_mm:
                    continue
                count += 1
                found = False
                for block, bounds in candidates:
                    if not all(bounds[axis][0] - 1e-8 <= point[axis] <= bounds[axis][1] + 1e-8
                               for axis in range(3)):
                        continue
                    if inverse_block(block, point) is not None:
                        found = True
                        break
                if not found:
                    missed.append(point)
        records.append({"z_mm": z, "analytic_material_points": total_count,
                        "tested_material_points": count,
                        "excluded_near_ideal_boundary_points": total_count - count,
                        "uncovered_inverse_points": missed})
        failures.extend(missed)
    return {"sections": records, "uncovered_inverse_points": failures,
            "excluded_near_ideal_boundary_count": sum(
                row["excluded_near_ideal_boundary_points"] for row in records),
            "fidelity_band_mm": fidelity_band_mm}


def sampled_step_sections(shape, center_um):
    """Classify deterministic section samples with OCC and the analytic set."""
    set_material_center(center_um)
    mismatches = []
    checked = 0
    for z in (2.0, 6.0):
        for x in SECTION_GRID:
            for y in SECTION_GRID:
                point = (x + 0.137, y + 0.271, z)
                expected = analytic_material(point)
                actual = bool(shape.isInside(App.Vector(*point), 1e-8, True))
                checked += 1
                if expected != actual:
                    mismatches.append({"point_mm": point, "expected_material": expected,
                                       "occ_material": actual})
    return {"checked": checked, "mismatches": mismatches}


def topology_signature(model):
    return j1_geometry.topology(model.get("blocks", []))


def physical_continuity(event_um):
    exact = analytic_metrics(event_um)
    left = analytic_metrics(event_um - 1)
    right = analytic_metrics(event_um + 1)
    step = 0.001
    left_secant = (exact["volume_mm3"] - left["volume_mm3"]) / step
    right_secant = (right["volume_mm3"] - exact["volume_mm3"]) / step
    continuity_steps = []
    for h in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
        v_left = analytic_volume_mm3(event_um / 1000.0 - h)
        v_right = analytic_volume_mm3(event_um / 1000.0 + h)
        continuity_steps.append({
            "h_mm": h,
            "left_volume_delta_mm3": abs(v_left - exact["volume_mm3"]),
            "right_volume_delta_mm3": abs(v_right - exact["volume_mm3"]),
            "left_secant_mm2": (exact["volume_mm3"] - v_left) / h,
            "right_secant_mm2": (v_right - exact["volume_mm3"]) / h,
        })
    smallest = continuity_steps[-1]
    return {
        "event_um": event_um,
        "left_volume_mm3": left["volume_mm3"],
        "event_volume_mm3": exact["volume_mm3"],
        "right_volume_mm3": right["volume_mm3"],
        "volume_continuity_checks": continuity_steps,
        "volume_continuous": max(smallest["left_volume_delta_mm3"],
                                  smallest["right_volume_delta_mm3"]) < 1e-8,
        "left_secant_mm2": left_secant,
        "event_gradient_mm2": exact["dvolume_dcenter_mm2"],
        "right_secant_mm2": right_secant,
        "secant_gradient_jump_mm2": abs(left_secant - right_secant),
        "smallest_step_secant_gradient_jump_mm2": abs(
            smallest["left_secant_mm2"] - smallest["right_secant_mm2"]
        ),
        "physical_gradient_exists_at_event": exact["physical_volume_gradient_exists"],
        "representation_jet_at_event": exact["representation_jet_available"],
        "neighboring_representation_jets": [left["representation_jet_available"],
                                              right["representation_jet_available"]],
    }


def run_build(cli, work, center_um, family="travelling_bore_v1"):
    name = f"center_{center_um:+07d}um"
    recipe_path = work / "recipes" / f"{name}.json"
    model_path = work / "models" / f"{name}.json"
    step_path = work / "step" / f"{name}.step"
    recipe_path.parent.mkdir(parents=True, exist_ok=True)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    step_path.parent.mkdir(parents=True, exist_ok=True)
    recipe_path.write_text(json.dumps({"family": family, "center_x_um": center_um}, indent=2) + "\n")
    start = time.perf_counter_ns()
    result = subprocess.run([str(cli), "build", str(recipe_path), str(model_path)],
                            cwd=ROOT, capture_output=True, text=True, check=False)
    build_ms = (time.perf_counter_ns() - start) / 1e6
    try:
        status = json.loads(result.stdout)
    except json.JSONDecodeError:
        status = {"unparsed_stdout": result.stdout.strip()}
    if result.returncode != 0 or not model_path.exists():
        return {"build_status": "refused", "returncode": result.returncode,
                "build_wall_ms": build_ms, "status": status, "stderr": result.stderr.strip()}
    model = json.loads(model_path.read_text())
    export_start = time.perf_counter_ns()
    exported = subprocess.run([str(cli), "export-step", str(model_path), str(step_path)],
                              cwd=ROOT, capture_output=True, text=True, check=False)
    export_ms = (time.perf_counter_ns() - export_start) / 1e6
    step_result = inspect_model(model, center_um, build_ms, status)
    step_result["step_export_wall_ms"] = export_ms
    step_result["step_export_status"] = exported.stdout.strip()
    step_result["step_export_returncode"] = exported.returncode
    if exported.returncode == 0 and step_path.exists():
        try:
            shape = Part.read(str(step_path))
            cad = reference_cad(center_um)
            step_result["step_boundary"] = inspect_shape(shape, analytic_metrics(center_um), center_um)
            step_result["step_vs_cad_volume_error_mm3"] = abs(shape.Volume - cad.Volume)
            step_result["step_vs_cad_area_error_mm2"] = abs(shape.Area - cad.Area)
            step_result["step_sections"] = sampled_step_sections(shape, center_um)
            boundary = step_result["step_boundary"]
            step_result["step_boundary"]["status"] = "read_and_inspected"
            step_result["step_boundary"]["narrow_edge_category"] = (
                "near_tangency_or_sub_tolerance_feature"
                if min(abs(center_um - event) for event in TANGENCIES_UM) <= 1
                else None)
            if family == "travelling_bore_v1":
                step_result["qualified_by_oracle"] = bool(step_result.get("qualified_by_oracle")) and (
                    boundary["solid_count"] == 1 and boundary["all_shells_closed"]
                    and boundary["brep_valid"]
                    and boundary["volume_relative_error"] <= VOLUME_REL_TOL
                    and step_result["step_sections"]["mismatches"] == [])
        except Exception as error:
            step_result["step_boundary"] = {
                "status": "occ_read_or_inspection_failed",
                "error": f"{type(error).__name__}: {error}",
                "narrow_edge_category": "near_tangency_or_sub_tolerance_feature"
                if min(abs(center_um - event) for event in TANGENCIES_UM) <= 1 else "other_occ_failure",
            }
    else:
        step_result["step_boundary"] = {
            "status": "export_refused", "stderr": exported.stderr.strip(),
            "narrow_edge_category": "near_tangency_or_sub_tolerance_feature"
            if min(abs(center_um - event) for event in TANGENCIES_UM) <= 1 else "other_export_failure",
        }
    return step_result


def inspect_shape(shape, metrics, center_um):
    cad = reference_cad(center_um)
    volume_error = abs(shape.Volume - metrics["volume_mm3"])
    volume_relative_error = volume_error / max(metrics["volume_mm3"], 1.0)
    area_error = abs(shape.Area - metrics["total_boundary_area_mm2"])
    shells = [shell for solid in shape.Solids for shell in solid.Shells]
    return {
        "solid_count": len(shape.Solids),
        "shell_count": len(shells),
        "all_shells_closed": bool(shells) and all(bool(shell.isClosed()) for shell in shells),
        "brep_valid": bool(shape.isValid()),
        "brep_closed": bool(shape.isClosed()),
        "face_count": len(shape.Faces),
        "edge_count": len(shape.Edges),
        "volume_mm3": shape.Volume,
        "volume_reference_mm3": metrics["volume_mm3"],
        "volume_error_mm3": volume_error,
        "volume_relative_error": volume_relative_error,
        "area_mm2": shape.Area,
        "area_reference_by_feature_sum_mm2": metrics["total_boundary_area_mm2"],
        "area_error_mm2": area_error,
        "occ_reference_volume_error_mm3": abs(cad.Volume - metrics["volume_mm3"]),
        "occ_reference_area_error_mm2": abs(cad.Area - metrics["total_boundary_area_mm2"]),
    }


def inspect_model(model, center_um, build_ms, build_status):
    metrics = analytic_metrics(center_um)
    blocks = model.get("blocks", [])
    areas_orders = {}
    area_counts = {}
    diagnostics = []
    for order in GAUSS_ORDERS:
        feature_areas, counts, errors = face_feature_areas(model, order)
        areas_orders[str(order)] = feature_areas
        area_counts = counts
        if order == GAUSS_ORDERS[-1]:
            diagnostics = errors
    volume_orders = {str(order): model_volume(blocks, order) for order in GAUSS_ORDERS}
    final_volume = volume_orders[str(GAUSS_ORDERS[-1])]
    volume_error = abs(final_volume - metrics["volume_mm3"])
    volume_relative_error = volume_error / max(metrics["volume_mm3"], 1.0)
    final_areas = areas_orders[str(GAUSS_ORDERS[-1])]
    area_errors = {
        feature: abs(final_areas.get(feature, 0.0) - expected)
        / max(expected, 1.0)
        for feature, expected in metrics["areas_by_feature_mm2"].items()
    }
    qdelta = abs(volume_orders[str(GAUSS_ORDERS[-1])] - volume_orders[str(GAUSS_ORDERS[-2])])
    volume_converged = qdelta <= VOLUME_REL_TOL * max(metrics["volume_mm3"], 1.0)
    area_qdelta = {
        feature: abs(areas_orders[str(GAUSS_ORDERS[-1])].get(feature, 0.0)
                     - areas_orders[str(GAUSS_ORDERS[-2])].get(feature, 0.0))
        for feature in metrics["areas_by_feature_mm2"]
    }
    area_converged = max(area_qdelta.values(), default=0.0) \
        <= AREA_REL_TOL * max(metrics["total_boundary_area_mm2"], 1.0)
    radial_errors = design_face_accuracy(model, center_um)
    bounded = model.get("schema") == "travelling_bore_bounded_v1"
    samples = materialized_samples(blocks, center_um,
                                   BOUNDED_FIDELITY_MM if bounded else 0.0)
    inverse = section_samples(
        blocks, center_um,
        BOUNDED_FIDELITY_MM if model.get("schema") == "travelling_bore_bounded_v1" else 0.0)
    try:
        cad = reference_cad(center_um)
        cad_info = inspect_shape(cad, metrics, center_um)
        cad_info["status"] = "reference_built_and_inspected"
    except Exception as error:
        cad_info = {"status": "reference_occ_failed",
                    "error": f"{type(error).__name__}: {error}"}
    topology = topology_signature(model)
    incidence = model.get("topology", {})
    tags = sorted(face.get("feature", "") for face in model.get("design_faces", []))
    report = {
        "build_status": "built",
        "build_wall_ms": build_ms,
        "native_status": build_status,
        "schema": model.get("schema"),
        "block_count": len(blocks),
        "source_topology": topology,
        "native_topology": incidence,
        "analytic_reference": metrics,
        "volume_by_gauss_order_mm3": volume_orders,
        "volume_gauss_previous_final_delta_mm3": qdelta,
        "volume_quadrature_converged": volume_converged,
        "volume_error_final_mm3": volume_error,
        "volume_relative_error_final": volume_relative_error,
        "surface_area_by_feature_and_order_mm2": areas_orders,
        "surface_feature_counts": area_counts,
        "surface_relative_errors_final": area_errors,
        "surface_quadrature_diagnostics": diagnostics,
        "surface_area_gauss_previous_final_delta_mm2": area_qdelta,
        "surface_quadrature_converged": area_converged,
        "design_face_maximum_radial_or_plane_error_mm": radial_errors,
        "materialized_interior_samples": samples,
        "analytic_section_inverse_samples": inverse,
        "occ_constructed_reference": cad_info,
        "design_face_tags": sorted(set(tags)),
        "signature": {
            "physical": [metrics["bore_pocket_relation"], metrics["bore_plate_relation"]],
            "representation_jet_available": metrics["representation_jet_available"],
            "selection": sorted((face.get("block"), face.get("side"), face.get("feature"))
                                for face in model.get("design_faces", [])),
        },
    }
    if model.get("schema") == "travelling_bore_bounded_v1":
        native = build_status if isinstance(build_status, dict) else {}
        report["representation_fidelity"] = {
            "representation": native.get("representation"),
            "fidelity_target": native.get("fidelity_target"),
            "fidelity_budget_mm": native.get("fidelity_budget_mm"),
            "radial_error_upper_bound_mm": native.get("radial_error_upper_bound_mm"),
            "hausdorff_error_upper_bound_mm": native.get("hausdorff_error_upper_bound_mm"),
            "observed_sampled_radial_error_mm": max(radial_errors.get("bore_wall", [0.0])
                                                       if isinstance(radial_errors.get("bore_wall"), list)
                                                       else [radial_errors.get("bore_wall", 0.0)]),
            "observed_sampled_plane_error_mm": max(
                (value for feature, value in radial_errors.items() if feature != "bore_wall"),
                default=0.0),
        }
    return report


def invalid_recipe_cases(cli, work, family):
    bad = [
        ("below_minimum", {"family": family, "center_x_um": -5001}),
        ("above_maximum", {"family": family, "center_x_um": 14001}),
        ("float_coordinate", {"family": family, "center_x_um": 8.0}),
        ("boolean_coordinate", {"family": family, "center_x_um": True}),
        ("string_coordinate", {"family": family, "center_x_um": "8000"}),
        ("missing_coordinate", {"family": family}),
    ]
    if family == "travelling_bore_v1":
        bad.extend([
            ("unsupported_irrational_left", {"family": family, "center_x_um": -100}),
            ("unsupported_irrational_right", {"family": family, "center_x_um": 11900}),
        ])
    records = {}
    state_path = work / "session" / "current.model.json"
    if not state_path.exists():
        transition_state(cli, work, 3000, family)
    for name, recipe in bad:
        recipe_path = work / "invalid" / f"{name}.json"
        recipe_path.parent.mkdir(parents=True, exist_ok=True)
        recipe_path.write_text(json.dumps(recipe, indent=2) + "\n")
        before = state_path.read_bytes() if state_path.exists() else None
        process = subprocess.run([str(cli), "build", str(recipe_path), str(state_path)],
                                 cwd=ROOT, capture_output=True, text=True, check=False)
        after = state_path.read_bytes() if state_path.exists() else None
        unchanged = before is not None and before == after
        if not unchanged and before is not None:
            state_path.write_bytes(before)
        try:
            status = json.loads(process.stdout)
        except json.JSONDecodeError:
            status = {"unparsed_stdout": process.stdout.strip()}
        records[name] = {
            "returncode": process.returncode,
            "status": status,
            "stderr": process.stderr.strip(),
            "target_unchanged_after_refusal": unchanged,
            "analytically_representable": (rational_intersections(recipe["center_x_um"])
                                             if isinstance(recipe.get("center_x_um"), int)
                                             and not isinstance(recipe.get("center_x_um"), bool)
                                             else False),
        }
    return records


def analytic_only_case(center_um):
    metrics = analytic_metrics(center_um)
    shape = reference_cad(center_um)
    occ = inspect_shape(shape, metrics, center_um)
    return {
        "analytic_reference": metrics,
        "occ_constructed_reference": occ,
        "qualified_analytic_occ_reference": (
            occ["solid_count"] == 1 and occ["all_shells_closed"] and occ["brep_valid"]
            and occ["volume_relative_error"] <= VOLUME_REL_TOL
            and occ["area_error_mm2"] / max(metrics["total_boundary_area_mm2"], 1.0) <= AREA_REL_TOL
        ),
    }


def transition_state(cli, work, center_um, family="travelling_bore_v1"):
    """Exercise repeated replacement and preserve the last state on refusal."""
    path = work / "session" / "current.model.json"
    recipe_path = work / "session" / "next-recipe.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    recipe_path.write_text(json.dumps({"family": family, "center_x_um": center_um}, indent=2) + "\n")
    before = path.read_bytes() if path.exists() else None
    process = subprocess.run([str(cli), "build", str(recipe_path), str(path)],
                             cwd=ROOT, capture_output=True, text=True, check=False)
    after = path.read_bytes() if path.exists() else None
    try:
        status = json.loads(process.stdout)
    except json.JSONDecodeError:
        status = {"unparsed_stdout": process.stdout.strip()}
    accepted = process.returncode == 0 and after is not None and status.get("status") in ("valid", "built")
    return {
        "center_um": center_um,
        "accepted": accepted,
        "returncode": process.returncode,
        "status": status,
        "stderr": process.stderr.strip(),
        "last_state_preserved_on_refusal": accepted or (before is not None and before == after),
        "previous_state_existed": before is not None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", choices=("travelling_bore_v1", "travelling_bore_bounded_v1"),
                        default="travelling_bore_v1")
    parser.add_argument("--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate")
    parser.add_argument("--workdir", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--analytic-only", action="store_true",
                         help="validate the independent analytic and OCC references without calling the CLI")
    parser.add_argument("--extra-centers-um", type=int, nargs="*", default=[])
    parser.add_argument("--centers-um", type=int, nargs="+",
                        help="override the frozen center sequence for a small diagnostic run")
    parser.add_argument("--require-qualified", action="store_true")
    args = parser.parse_args()
    cli = args.cli if args.cli.is_absolute() else ROOT / args.cli
    work = args.workdir if args.workdir.is_absolute() else ROOT / args.workdir
    work.mkdir(parents=True, exist_ok=True)
    frozen_centers = CENTERS_UM if args.family == "travelling_bore_v1" else BOUNDED_CENTERS_UM
    centers = list(args.centers_um) if args.centers_um else list(frozen_centers) + list(args.extra_centers_um)
    report = {
        "family": args.family,
        "contract_base_revision": "fba87ab" if args.family == "travelling_bore_v1" else "5d7bb40",
        "recipe_key": "center_x_um",
        "radius_mm": 2.0,
        "plate_mm": {"x": [-10.0, 10.0], "y": [-10.0, 10.0], "z": [0.0, 8.0]},
        "pocket_mm": {"x": [-6.0, -2.0], "y": [-4.0, 4.0], "z": [4.0, 8.0]},
        "required_centers_um": list(frozen_centers),
        "quadrature_orders": list(GAUSS_ORDERS),
        "volume_relative_tolerance": VOLUME_REL_TOL,
        "area_relative_tolerance": AREA_REL_TOL,
        "representation": "exact_rational_arcs" if args.family == "travelling_bore_v1"
            else "quantized_rational_arcs",
        "fidelity_budget_mm": None if args.family == "travelling_bore_v1" else 1e-8,
        "rational_intersection_rule": "sqrt(2000^2-delta_um^2) must be an integer micrometre",
        "fidelity_target": None if args.family == "travelling_bore_v1"
            else "ideal_clipped_bore_arcs",
        "continuity_at_exact_events": {str(event): physical_continuity(event)
                                        for event in TANGENCIES_UM},
        "cases": {},
    }
    state_transitions = []
    previous_signature = None
    for sequence_index, center_um in enumerate(centers, 1):
        key = f"step_{sequence_index:02d}_center_{center_um:+07d}um"
        try:
            if args.analytic_only:
                case = {"status": "analytic_only", **analytic_only_case(center_um)}
            else:
                case = run_build(cli, work, center_um, args.family)
                transition = transition_state(cli, work, center_um, args.family)
                case["state_transition"] = transition
                reference = analytic_metrics(center_um)
                case.setdefault("analytic_reference", reference)
                expected_admission = not reference["tangency_events"] and (
                    args.family == "travelling_bore_bounded_v1"
                    or reference["exact_intersections_representable"]
                )
                case["expected_native_admission"] = expected_admission
                case["native_outcome_matches_expected_admission"] = (
                    case.get("build_status") == ("built" if expected_admission else "refused")
                )
                if case.get("build_status") == "built":
                    metrics = case["analytic_reference"]
                    samples = case["materialized_interior_samples"]
                    inverse = case["analytic_section_inverse_samples"]
                    cad = case["occ_constructed_reference"]
                    common_qualification = (
                        case["volume_relative_error_final"] <= VOLUME_REL_TOL
                        and case["volume_quadrature_converged"]
                        and not case["surface_quadrature_diagnostics"]
                        and not samples["classification_failures"]
                        and samples["nonpositive_sampled_determinants"] == 0
                        and not inverse["uncovered_inverse_points"]
                        and expected_admission
                    )
                    if args.family == "travelling_bore_bounded_v1":
                        fidelity = case.get("representation_fidelity", {})
                        case["qualified_by_oracle"] = common_qualification and (
                            fidelity.get("representation") == "quantized_rational_arcs"
                            and fidelity.get("fidelity_target") == "ideal_clipped_bore_arcs"
                            and fidelity.get("fidelity_budget_mm") == BOUNDED_FIDELITY_MM
                            and isinstance(fidelity.get("radial_error_upper_bound_mm"), (int, float))
                            and 0.0 <= fidelity["radial_error_upper_bound_mm"] <= BOUNDED_FIDELITY_MM
                            and isinstance(fidelity.get("hausdorff_error_upper_bound_mm"), (int, float))
                            and 0.0 <= fidelity["hausdorff_error_upper_bound_mm"] <= BOUNDED_FIDELITY_MM
                            and fidelity.get("observed_sampled_radial_error_mm", math.inf)
                            <= BOUNDED_FIDELITY_MM
                            and fidelity.get("observed_sampled_plane_error_mm", math.inf) <= 1e-12
                        )
                        case["qualification_scope"] = (
                            "native fidelity bounds and independent analytic/geometric diagnostics; "
                            "STEP/OCC status is reported separately and is not a substitute for certification"
                        )
                    else:
                        case["qualified_by_oracle"] = common_qualification and (
                        case["surface_quadrature_converged"]
                        and max(case["surface_relative_errors_final"].values(), default=0.0) <= AREA_REL_TOL
                        and max(case["design_face_maximum_radial_or_plane_error_mm"].values(), default=0.0)
                        <= RADIUS_TOL_MM
                        and cad["solid_count"] == 1
                        and cad["all_shells_closed"]
                        and cad["brep_valid"]
                        and cad["volume_relative_error"] <= VOLUME_REL_TOL
                        and cad["area_error_mm2"] / max(metrics["total_boundary_area_mm2"], 1.0)
                        <= AREA_REL_TOL
                        )
                    case["signature_change_from_previous_accepted"] = (
                        None if previous_signature is None
                        else case["signature"] != previous_signature
                    )
                    previous_signature = case["signature"]
                    if center_um in TANGENCIES_UM:
                        case["event_expectation"] = "preserve_prior_state_and_report_contact"
            report["cases"][key] = case
            report["cases"][key]["sequence_index"] = sequence_index
        except Exception as error:
            report["cases"][key] = {"status": "oracle_error", "error": f"{type(error).__name__}: {error}"}
        case = report["cases"][key]
        case["sequence_index"] = sequence_index
        if not args.analytic_only:
            state_transitions.append(case.get("state_transition", {}))
        print(key, case.get("build_status", case.get("status")),
              case.get("qualified_by_oracle"), flush=True)
    if not args.analytic_only:
        report["invalid_inputs"] = invalid_recipe_cases(cli, work, args.family)
        report["sequential_state_transitions"] = state_transitions
    report["coverage"] = {
        "cases": len(report["cases"]),
        "built": sum(case.get("build_status") == "built" for case in report["cases"].values()),
        "refused": sum(case.get("build_status") == "refused" for case in report["cases"].values()),
        "qualified": sum(case.get("qualified_by_oracle", False) for case in report["cases"].values()),
        "expected_admissions": sum(
            not analytic_metrics(center)["tangency_events"]
            and (args.family == "travelling_bore_bounded_v1"
                 or analytic_metrics(center)["exact_intersections_representable"])
            for center in centers
        ),
        "expected_event_or_unrepresentable_refusals": sum(
            bool(analytic_metrics(center)["tangency_events"])
            or (args.family == "travelling_bore_v1"
                and not analytic_metrics(center)["exact_intersections_representable"])
            for center in centers
        ),
        "expected_cases": len(centers),
        "failures": {
            key: case for key, case in report["cases"].items()
            if case.get("status") == "oracle_error"
            or (case.get("build_status") in ("built", "refused")
                and not case.get("native_outcome_matches_expected_admission", True))
            or (case.get("build_status") == "built" and not case.get("qualified_by_oracle", False))
            or (case.get("build_status") == "refused"
                and case.get("expected_native_admission") is False
                and not case.get("state_transition", {}).get("last_state_preserved_on_refusal", False))
        },
    }
    output = args.report or work / "oracle.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.require_qualified and report["coverage"]["qualified"] != report["coverage"]["expected_admissions"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
