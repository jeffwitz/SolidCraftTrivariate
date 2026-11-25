#!/usr/bin/env python3
"""Independent exact-topology and STEP oracle for compact FieldCraft brackets.

The model checks use Python Fractions and do not call the Rust audit. The STEP
checks delegate geometry, mass properties, and material probes to the existing
FreeCAD/OCCT replay oracle.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import importlib.util
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import sys


CASES = ("x40-y80", "x50-y74", "x51-y74", "x64-y64", "x80-y40")
ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT / "build/fieldcraft-replay"
HOST = Path("/home/jeff/Code/SolidCraftTrivariate/build/solidcraft-bridge")
REPLAY = ROOT / "tests/references/fieldcraft_kernel_replay.py"


def transformed_control(control: list[int], transform: dict) -> tuple[F, F, F, F]:
    matrix = transform.get("matrix", [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    translation = transform.get("translation", [0, 0, 0])
    denominator = F(transform.get("denominator", 1))
    xyz = [sum(F(matrix[a][j]) * control[j] for j in range(3))
           + F(denominator * translation[a]) * control[3]
           for a in range(3)]
    return tuple(xyz + [denominator * control[3]])


def face_axes(side: str) -> tuple[int, tuple[int, int]]:
    fixed = "uvw".index(side[0])
    return fixed, tuple(axis for axis in range(3) if axis != fixed)


def face_controls(block: dict, side: str, transform: dict) -> tuple[tuple[int, int], list[list[tuple[F, ...]]]]:
    degrees = block["degrees"]
    dims = [degree + 1 for degree in degrees]
    fixed_axis, free_axes = face_axes(side)
    fixed_index = 0 if side[1] == "0" else dims[fixed_axis] - 1
    rows = []
    for a in range(dims[free_axes[0]]):
        row = []
        for b in range(dims[free_axes[1]]):
            ix = [0, 0, 0]
            ix[fixed_axis] = fixed_index
            ix[free_axes[0]] = a
            ix[free_axes[1]] = b
            flat = (ix[0] * dims[1] + ix[1]) * dims[2] + ix[2]
            row.append(transformed_control(block["controls"][flat], transform))
        rows.append(row)
    return tuple(degrees[a] for a in free_axes), rows


def orient_rows(rows: list[list[tuple[F, ...]]], orientation: dict,
                first_shape: tuple[int, int]) -> list[list[tuple[F, ...]]]:
    """Reindex the second trace onto first-face coordinates."""
    result = []
    for a in range(first_shape[0]):
        row = []
        for b in range(first_shape[1]):
            x = first_shape[0] - 1 - a if orientation["flip_first"] else a
            y = first_shape[1] - 1 - b if orientation["flip_second"] else b
            source_a, source_b = (y, x) if orientation["swap"] else (x, y)
            row.append(rows[source_a][source_b])
        result.append(row)
    return result


def canonical_homogeneous(rows: list[list[tuple[F, ...]]]) -> list[list[tuple[F, ...]]]:
    scale = rows[0][0][3]
    return [[tuple(value / scale for value in control) for control in row] for row in rows]


def elevate_curve(controls: list[tuple[F, ...]]) -> list[tuple[F, ...]]:
    degree = len(controls) - 1
    elevated = [controls[0]]
    for index in range(1, degree + 1):
        left = F(index, degree + 1)
        right = 1 - left
        elevated.append(tuple(left * a + right * b
                              for a, b in zip(controls[index - 1], controls[index])))
    elevated.append(controls[-1])
    return elevated


def elevate_face(rows: list[list[tuple[F, ...]]], target: tuple[int, int]) -> list[list[tuple[F, ...]]]:
    result = rows
    while len(result) - 1 < target[0]:
        columns = list(zip(*result))
        elevated_columns = [elevate_curve(list(column)) for column in columns]
        result = [list(row) for row in zip(*elevated_columns)]
    while len(result[0]) - 1 < target[1]:
        result = [elevate_curve(row) for row in result]
    return result


def exact_interface_audit(model: dict) -> dict:
    blocks = {block["id"]: block for block in model["blocks"]}
    transform = model.get("transform", {})
    counts = Counter()
    for interface in model["topology"]["interfaces"]:
        first_id, first_side = interface["first"]
        second_id, second_side = interface["second"]
        first_degree, first = face_controls(blocks[first_id], first_side, transform)
        second_degree, second = face_controls(blocks[second_id], second_side, transform)
        orientation = interface["orientation"]
        if orientation["swap"]:
            second_degree = tuple(reversed(second_degree))
        target_degree = tuple(max(a, b) for a, b in zip(first_degree, second_degree))
        first = elevate_face(first, target_degree)
        second_target = tuple(reversed(target_degree)) if orientation["swap"] else target_degree
        second = elevate_face(second, second_target)
        second = orient_rows(second, orientation,
                             (target_degree[0] + 1, target_degree[1] + 1))
        if canonical_homogeneous(first) != canonical_homogeneous(second):
            raise AssertionError((first_id, first_side, second_id, second_side,
                                  "exact homogeneous trace mismatch"))
        counts["exact_interfaces"] += 1
        if "fillet" in first_id or "fillet" in second_id:
            counts["elbow_interfaces"] += 1

    exterior = model["topology"]["exterior"]
    coplanar = {"x0": [], "y0": []}
    for block_id, side in exterior:
        _, controls = face_controls(blocks[block_id], side, transform)
        flat = [control for row in controls for control in row]
        for axis, key in ((0, "x0"), (1, "y0")):
            if all(control[axis] == 0 for control in flat):
                coplanar[key].append((block_id, side))
    fillet_blocks = sorted(block["id"] for block in model["blocks"]
                           if "fillet" in block["id"])
    return {
        "blocks": len(blocks),
        "interfaces": len(model["topology"]["interfaces"]),
        "exterior_sides": len(exterior),
        "exact_trace_interfaces": counts["exact_interfaces"],
        "elbow_trace_interfaces": counts["elbow_interfaces"],
        "x0_coplanar_exterior_sides": coplanar["x0"],
        "y0_coplanar_exterior_sides": coplanar["y0"],
        "fillet_blocks": fillet_blocks,
    }


def load_replay_module():
    existing = globals().get("_REPLAY_MODULE")
    if existing is not None:
        return existing
    spec = importlib.util.spec_from_file_location("fieldcraft_kernel_replay", REPLAY)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    # FreeCAD's Python bootstrap removes a `json` name from __main__ in some
    # builds, so restore this script's standard-library binding after import.
    globals()["json"] = __import__("json")
    globals()["contextlib"] = __import__("contextlib")
    globals()["io"] = __import__("io")
    globals()["_REPLAY_MODULE"] = module
    return module


def tangency_plane_probes(replay, reference, candidate, source_bores: dict) -> tuple[list[dict], int]:
    """Probe material at z=cz±r, beyond the other radial coordinate, mid-wall."""
    checks = []
    failures = 0
    for axis in ("x", "y"):
        bore = source_bores[axis]
        center = bore["center_mm"]
        radius = float(bore["radius_mm"])
        for z_sign in (-1.0, 1.0):
            z = float(center[2]) + z_sign * radius
            for radial_sign in (-1.0, 1.0):
                offset = radial_sign * (radius + 0.02)
                if axis == "x":
                    point = (4.0, float(center[1]) + offset, z)
                else:
                    point = (float(center[0]) + offset, 4.0, z)
                vector = replay.App.Vector(*point)
                source_inside = bool(reference.isInside(vector, 1e-7, True))
                candidate_inside = bool(candidate.isInside(vector, 1e-7, True))
                passed = source_inside and candidate_inside
                failures += int(not passed)
                checks.append({
                    "axis": axis,
                    "point_mm": list(point),
                    "z_tangency_sign": int(z_sign),
                    "radial_sign": int(radial_sign),
                    "outside_radius_by_mm": 0.02,
                    "expected_inside": True,
                    "reference_inside": source_inside,
                    "candidate_inside": candidate_inside,
                    "pass": passed,
                })
    return checks, failures


def audit_case(case: str, compact_dir: Path, old_dir: Path, host_dir: Path,
               previous_model_dir: Path, occ_volume: Path, full: bool) -> dict:
    model_path = compact_dir / f"{case}.model.json"
    candidate_path = compact_dir / f"{case}.step"
    host_path = host_dir / case / "host.step"
    if not all(path.is_file() for path in (model_path, candidate_path, host_path)):
        raise FileNotFoundError(f"missing model or STEP input for {case}")
    model = json.loads(model_path.read_text())
    topology = exact_interface_audit(model)
    previous_path = previous_model_dir / f"{case}.model.json"
    if not previous_path.is_file():
        previous_path = previous_model_dir / f"fieldcraft_{case.replace('-', '_')}.model.json"
    previous = json.loads(previous_path.read_text())
    previous_topology = exact_interface_audit(previous)
    if len(topology["x0_coplanar_exterior_sides"]) != 5:
        raise AssertionError(f"{case}: expected five exterior map sides on x=0")
    if len(topology["y0_coplanar_exterior_sides"]) != 5:
        raise AssertionError(f"{case}: expected five exterior map sides on y=0")

    replay = load_replay_module()
    script = old_dir / "inputs" / f"{case}.script.json"
    recipe = old_dir / f"{case}.recipe.json"
    with contextlib.redirect_stdout(io.StringIO()):
        record = replay.run(host_path, candidate_path, occ_volume,
                            compact_dir / f"{case}-oracle", fast=not full,
                            source_script=script, candidate_recipe=recipe)
    src_bores = replay.source_holes(script)
    source_shape = replay.read_shape(host_path)
    candidate_shape = replay.read_shape(candidate_path)
    tangencies, tangent_failures = tangency_plane_probes(replay, source_shape, candidate_shape, src_bores)
    record["tangency_plane_probes"] = tangencies
    record["tangency_plane_probe_failures"] = tangent_failures
    record["qualification_pass"] = record["qualification_pass"] and tangent_failures == 0
    # The former 48-block graph has 24 exterior map sides on each wing plane.
    topology["previous_model_exterior_side_count"] = previous_topology["exterior_sides"]
    topology["previous_x0_coplanar_exterior_sides"] = len(previous_topology["x0_coplanar_exterior_sides"])
    topology["previous_y0_coplanar_exterior_sides"] = len(previous_topology["y0_coplanar_exterior_sides"])
    if topology["previous_x0_coplanar_exterior_sides"] != 24 or topology["previous_y0_coplanar_exterior_sides"] != 24:
        raise AssertionError(f"{case}: expected 24 historical exterior map sides on each wing plane")
    return {"case": case, "topology": topology, "step_oracle": record}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compact-dir", type=Path, default=ROOT / "build/fieldcraft-compact")
    parser.add_argument("--old-dir", type=Path, default=OLD)
    parser.add_argument("--previous-model-dir", type=Path)
    parser.add_argument("--host-dir", type=Path, default=HOST)
    parser.add_argument("--occ-volume", type=Path, default=HOST / "occ-volume")
    parser.add_argument("--case", choices=CASES, action="append")
    parser.add_argument("--fast", action="store_true", help="use one coarse classification grid")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not args.occ_volume.is_file():
        # The audited binary is built next to the host exports in current workspaces.
        fallback = HOST.parent / "solidcraft-bridge/occ-volume"
        if not fallback.is_file():
            raise FileNotFoundError("occ-volume helper not found; pass --occ-volume")
        args.occ_volume = fallback
    if args.previous_model_dir is None:
        args.previous_model_dir = args.old_dir
        if (args.compact_dir / "previous").is_dir():
            args.previous_model_dir = args.compact_dir / "previous"
    results = [audit_case(case, args.compact_dir, args.old_dir, args.host_dir,
                          args.previous_model_dir,
                          args.occ_volume, not args.fast)
               for case in (args.case or CASES)]
    payload = {"status": "pass" if all(row["step_oracle"]["qualification_pass"]
                                       for row in results) else "fail",
               "cases": results}
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(encoded, end="")
    if payload["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
