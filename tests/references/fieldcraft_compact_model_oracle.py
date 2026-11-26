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
import hashlib
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import subprocess


CASES = ("x40-y80", "x50-y74", "x51-y74", "x64-y64", "x80-y40")
ROOT = Path(__file__).resolve().parents[2]
CURRENT = ROOT / "docs/site/data"
SOURCE_SCRIPTS = ROOT / "tests/references/fieldcraft_scripts"
HOSTS = ROOT / "build/solidcraft-bridge"
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


def audit_case(case: str, cli: Path, work: Path, current_dir: Path,
               source_dir: Path, host_dir: Path, occ_volume: Path,
               full: bool) -> dict:
    stem = f"fieldcraft_{case.replace('-', '_')}"
    case_dir = current_dir / stem
    model_path = case_dir / f"{stem}.model.json"
    recipe_path = case_dir / f"{stem}.recipe.json"
    candidate_path = case_dir / f"{stem}.step"
    host_path = host_dir / case / "host.step"
    source_script = source_dir / f"{case}.script.json"
    if not all(path.is_file() for path in (model_path, recipe_path, candidate_path,
                                           host_path)):
        raise FileNotFoundError(
            f"missing current recipe/model/STEP or independent source input for {case}"
        )
    model = json.loads(model_path.read_text())
    recipe = json.loads(recipe_path.read_text())
    work.mkdir(parents=True, exist_ok=True)
    rebuilt_model = work / f"{case}.model.json"
    built = subprocess.run([str(cli), "build", str(recipe_path), str(rebuilt_model)],
                           capture_output=True, text=True)
    if built.returncode:
        raise RuntimeError(f"{case}: current recipe failed to build: {built.stderr or built.stdout}")
    build_summary = json.loads(built.stdout)
    audited = subprocess.run([str(cli), "audit", str(rebuilt_model)],
                             capture_output=True, text=True)
    if audited.returncode:
        raise RuntimeError(f"{case}: rebuilt model failed audit: {audited.stderr or audited.stdout}")
    audit_summary = json.loads(audited.stdout)
    if audit_summary.get("status") != "valid":
        raise RuntimeError(f"{case}: rebuilt model audit status is {audit_summary.get('status')!r}")
    if json.loads(rebuilt_model.read_text()) != model:
        raise AssertionError(f"{case}: current recipe rebuild differs from published model")
    topology = exact_interface_audit(model)
    topology["model_path"] = model_path.relative_to(ROOT).as_posix()
    topology["recipe_path"] = recipe_path.relative_to(ROOT).as_posix()
    topology["candidate_step_path"] = candidate_path.relative_to(ROOT).as_posix()
    topology["independent_host_step_path"] = f"<host-dir>/{case}/host.step"
    if len(topology["x0_coplanar_exterior_sides"]) != 5:
        raise AssertionError(f"{case}: expected five exterior map sides on x=0")
    if len(topology["y0_coplanar_exterior_sides"]) != 5:
        raise AssertionError(f"{case}: expected five exterior map sides on y=0")

    replay = load_replay_module()
    with contextlib.redirect_stdout(io.StringIO()):
        record = replay.run(host_path, candidate_path, occ_volume,
                            work / f"{case}-oracle", fast=not full,
                            source_script=source_script if source_script.is_file() else None,
                            candidate_recipe=recipe_path)
    record["reference"] = f"<host-dir>/{case}/host.step"
    record["candidate"] = candidate_path.relative_to(ROOT).as_posix()
    src_bores = (replay.source_holes(source_script) if source_script.is_file()
                 else replay.candidate_holes(recipe_path))
    source_shape = replay.read_shape(host_path)
    candidate_shape = replay.read_shape(candidate_path)
    tangencies, tangent_failures = tangency_plane_probes(replay, source_shape, candidate_shape, src_bores)
    record["tangency_plane_probes"] = tangencies
    record["tangency_plane_probe_failures"] = tangent_failures
    record["qualification_pass"] = record["qualification_pass"] and tangent_failures == 0
    record["probe_coordinate_source"] = (
        source_script.relative_to(ROOT).as_posix() if source_script.is_file()
        else recipe_path.relative_to(ROOT).as_posix()
    )
    return {"case": case,
            "rebuild": {"recipe": recipe_path.relative_to(ROOT).as_posix(),
                        "model_matches_published": True,
                        "build": build_summary, "audit": audit_summary},
            "topology": topology, "step_oracle": record}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate")
    parser.add_argument("--work", type=Path, default=ROOT / "build/fieldcraft-compact-oracle")
    parser.add_argument("--current-dir", type=Path, default=CURRENT)
    parser.add_argument("--source-dir", type=Path, default=SOURCE_SCRIPTS)
    parser.add_argument("--host-dir", type=Path, default=HOSTS,
                        help="directory containing <case>/host.step from independent SolidCore exports")
    parser.add_argument("--occ-volume", type=Path, default=HOSTS / "occ-volume")
    parser.add_argument("--case", choices=CASES, action="append")
    parser.add_argument("--fast", action="store_true", help="use one coarse classification grid")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "docs/results/fieldcraft-compact-current-oracle.json")
    args = parser.parse_args()
    for name in ("work", "current_dir", "source_dir", "host_dir", "occ_volume", "output"):
        value = getattr(args, name, None)
        if value is not None and not value.is_absolute():
            setattr(args, name, ROOT / value)
    cli = args.cli.resolve()
    if not cli.is_file():
        raise FileNotFoundError(f"CLI binary not found: {cli}")
    if not args.occ_volume.is_file():
        raise FileNotFoundError("occ-volume helper not found; pass --occ-volume")
    cli_hash = hashlib.sha256(cli.read_bytes()).hexdigest()
    results = [audit_case(case, cli, args.work, args.current_dir, args.source_dir,
                          args.host_dir, args.occ_volume, not args.fast)
               for case in (args.case or CASES)]
    if hashlib.sha256(cli.read_bytes()).hexdigest() != cli_hash:
        raise RuntimeError("CLI binary changed during oracle run")
    replay = load_replay_module()
    payload = {"status": "pass" if all(row["step_oracle"]["qualification_pass"]
                                       for row in results) else "fail",
               "checkout_revision": subprocess.check_output(
                   ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
               ).strip(),
               "cli_path": cli.relative_to(ROOT).as_posix() if cli.is_relative_to(ROOT) else cli.name,
               "cli_sha256": cli_hash,
               "freecad_version": ".".join(replay.App.Version()[:3]),
               "cases_expected": len(CASES),
               "cases": results}
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded)
    print(encoded, end="")
    if payload["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
