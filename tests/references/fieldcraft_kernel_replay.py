#!/usr/bin/env python3
"""Compare an exported trivariate bracket against the locked FieldCraft STEP.

This uses FreeCAD's OpenCascade reader and Boolean operations as an independent
kernel oracle. It does not heal or modify either input shape.
"""

import FreeCAD as App  # noqa: F401  # Initialize the supported OCC runtime first.
import Part

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess


def read_shape(path):
    shape = Part.read(str(path))
    if shape.isNull():
        raise RuntimeError(f"could not read STEP: {path}")
    return shape


def measure(shape):
    return {
        "valid": bool(shape.isValid()),
        "solids": len(shape.Solids),
        "shells": len(shape.Shells),
        "faces": len(shape.Faces),
        "edges": len(shape.Edges),
        "vertices": len(shape.Vertexes),
        "volume_mm3": float(shape.Volume),
        "area_mm2": float(shape.Area),
        "bounds_mm": [
            float(shape.BoundBox.XMin), float(shape.BoundBox.XMax),
            float(shape.BoundBox.YMin), float(shape.BoundBox.YMax),
            float(shape.BoundBox.ZMin), float(shape.BoundBox.ZMax),
        ],
        "cylinders": sorted(
            [float(face.Surface.Radius) for face in shape.Faces
             if hasattr(face.Surface, "Radius")]
        ),
    }


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_holes(path):
    data = json.loads(path.read_text())
    holes = {}
    for command in data["commands"]:
        if command["command"] != "solid.hole":
            continue
        params = command["params"]
        x, y, z = params["position"]
        if abs(x - 8.0) < 1e-9 and y > 8.0:
            axis, center = "x", (4.0, y, z)
        elif abs(y - 8.0) < 1e-9 and x > 8.0:
            axis, center = "y", (x, 4.0, z)
        else:
            raise ValueError(f"unexpected FieldCraft bore location: {params['position']}")
        holes[axis] = {"center_mm": center, "radius_mm": float(params["diameter"]) / 2}
    if set(holes) != {"x", "y"}:
        raise ValueError(f"expected one x and one y bore, got {sorted(holes)}")
    return holes


def candidate_holes(path):
    data = json.loads(path.read_text())
    source = data.get("source", data)
    transform = data.get("transform", {})
    if transform.get("matrix", [[1, 0, 0], [0, 1, 0], [0, 0, 1]]) != [[1, 0, 0], [0, 1, 0], [0, 0, 1]] or transform.get("translation", [0, 0, 0]) != [0, 0, 0]:
        raise ValueError("oracle currently supports identity-affine FieldCraft replay recipes")
    denominator = float(transform.get("denominator", 1))
    holes = {
        bore["axis"]: {
            "center_mm": [float(value) / denominator for value in bore["center_mm"]],
            "radius_mm": float(bore["radius_mm"]) / denominator,
        }
        for bore in source["bores"]
    }
    if set(holes) != {"x", "y"}:
        raise ValueError(f"expected one x and one y bore, got {sorted(holes)}")
    return holes


def bore_probes(reference, candidate, source_bores, candidate_bores):
    checks = []
    mismatch_count = 0
    for axis in ("x", "y"):
        source = source_bores[axis]
        target = candidate_bores[axis]
        if (abs(float(source["radius_mm"]) - float(target["radius_mm"])) > 1e-12
                or any(abs(float(a) - float(b)) > 1e-12
                       for a, b in zip(source["center_mm"], target["center_mm"]))):
            mismatch_count += 1
        for direction in (0, 1):
            for sign in (-1.0, 1.0):
                for offset in (-0.02, 0.02):
                    radius = float(source["radius_mm"]) + offset
                    cx, cy, cz = map(float, source["center_mm"])
                    if axis == "x":
                        point = App.Vector(cx, cy + (sign * radius if direction == 0 else 0),
                                           cz + (sign * radius if direction == 1 else 0))
                    else:
                        point = App.Vector(cx + (sign * radius if direction == 0 else 0), cy,
                                           cz + (sign * radius if direction == 1 else 0))
                    expected_inside = offset > 0
                    source_inside = bool(reference.isInside(point, 1e-7, True))
                    candidate_inside = bool(candidate.isInside(point, 1e-7, True))
                    failed = source_inside != expected_inside or candidate_inside != expected_inside
                    mismatch_count += int(failed)
                    checks.append({
                        "axis": axis, "direction": "y" if axis == "x" and direction == 0 else
                        "z" if direction == 1 else "x", "sign": int(sign),
                        "radial_offset_mm": offset, "point_mm": [point.x, point.y, point.z],
                        "expected_inside": expected_inside, "reference_inside": source_inside,
                        "candidate_inside": candidate_inside, "pass": not failed,
                    })
    return checks, mismatch_count


def run(reference_path, candidate_path, occ_volume=None, artifact_dir=None, fast=False,
        source_script=None, candidate_recipe=None):
    reference = read_shape(reference_path)
    candidate = read_shape(candidate_path)
    grid = []
    for nx, ny, nz in (((8, 8, 6),) if fast else ((8, 8, 6), (12, 12, 8))):
        mismatch = 0
        tested = 0
        for ix in range(nx):
            x = 60.0 * (ix + 0.5) / nx
            for iy in range(ny):
                y = 50.0 * (iy + 0.5) / ny
                for iz in range(nz):
                    z = 40.0 * (iz + 0.5) / nz
                    point = App.Vector(x, y, z)
                    if reference.isInside(point, 1e-7, True) != candidate.isInside(point, 1e-7, True):
                        mismatch += 1
                    tested += 1
        grid.append({"shape": [nx, ny, nz], "tested": tested, "mismatches": mismatch})
    result = {
        "reference": str(reference_path),
        "candidate": str(candidate_path),
        "reference_shape": measure(reference),
        "candidate_shape": measure(candidate),
        "shape_difference_method": "uniform interior point classification; volume and area from adaptive OCCT Gauss-Kronrod integration",
        "bbox_max_abs_delta_mm": max(abs(a - b) for a, b in zip(
            [reference.BoundBox.XMin, reference.BoundBox.XMax, reference.BoundBox.YMin,
             reference.BoundBox.YMax, reference.BoundBox.ZMin, reference.BoundBox.ZMax],
            [candidate.BoundBox.XMin, candidate.BoundBox.XMax, candidate.BoundBox.YMin,
             candidate.BoundBox.YMax, candidate.BoundBox.ZMin, candidate.BoundBox.ZMax],
        )),
        "interior_classification_grids": grid,
        "input_sha256": {"reference_step": sha256(reference_path),
                         "candidate_step": sha256(candidate_path)},
    }
    if source_script is not None or candidate_recipe is not None:
        if source_script is None or candidate_recipe is None:
            raise ValueError("both source_script and candidate_recipe are required")
        if source_script is not None:
            src_bores = source_holes(source_script)
        elif candidate_recipe is not None:
            # If the historical authoring scripts are unavailable, use the
            # tracked current recipe for probe coordinates while retaining the
            # independent STEP's membership and mass-property checks.
            src_bores = candidate_holes(candidate_recipe)
        else:
            raise ValueError("source_script or candidate_recipe is required for bore checks")
        dst_bores = candidate_holes(candidate_recipe)
        source_radii = sorted([float(src_bores["x"]["radius_mm"]),
                               float(src_bores["y"]["radius_mm"]), 6.0])
        imported_radii = result["reference_shape"]["cylinders"]
        if len(imported_radii) != 3 or any(abs(a - b) > 1e-8 for a, b in zip(source_radii, imported_radii)):
            result["source_cylinder_radius_match"] = False
        else:
            result["source_cylinder_radius_match"] = True
        probes, failures = bore_probes(reference, candidate, src_bores, dst_bores)
        result["bore_parameters_mm"] = {"reference": src_bores, "candidate": dst_bores}
        result["radial_probes"] = probes
        result["radial_probe_failures"] = failures
        result["input_sha256"]["source_script"] = sha256(source_script)
        result["input_sha256"]["candidate_recipe"] = sha256(candidate_recipe)
        model_path = candidate_path.with_suffix(".model.json")
        audit_path = candidate_path.with_suffix(".audit.json")
        if model_path.exists():
            result["input_sha256"]["candidate_model"] = sha256(model_path)
        if audit_path.exists():
            result["input_sha256"]["candidate_audit"] = sha256(audit_path)
        rx, ry = src_bores["x"]["radius_mm"], src_bores["y"]["radius_mm"]
        expected = 34080.0 - math.pi * (360.0 + 8.0 * (rx * rx + ry * ry))
        result["analytic_volume_mm3"] = expected
    if occ_volume is not None:
        if artifact_dir is None:
            raise ValueError("artifact_dir is required with occ_volume")
        artifact_dir.mkdir(parents=True, exist_ok=True)
        adaptive = {}
        for name, shape in (("reference", reference), ("candidate", candidate)):
            brep_path = artifact_dir / f"{name}.brep"
            shape.exportBrep(str(brep_path))
            completed = subprocess.run(
                [str(occ_volume), str(brep_path)], check=True, text=True,
                capture_output=True,
            )
            adaptive[name] = json.loads(completed.stdout)
        result["adaptive_occ_mass_properties"] = adaptive
    checks = [
        result["reference_shape"]["valid"] and result["reference_shape"]["solids"] == 1,
        result["candidate_shape"]["valid"] and result["candidate_shape"]["solids"] == 1,
        result["interior_classification_grids"][-1]["mismatches"] == 0,
        result.get("radial_probe_failures", 0) == 0,
        result.get("source_cylinder_radius_match", True),
    ]
    if "adaptive_occ_mass_properties" in result:
        src_mass = result["adaptive_occ_mass_properties"]["reference"][-1]
        dst_mass = result["adaptive_occ_mass_properties"]["candidate"][-1]
        checks.extend([
            abs(src_mass["volume_mm3"] - dst_mass["volume_mm3"]) <= 1e-8 * src_mass["volume_mm3"],
            abs(src_mass["surface_mm2"] - dst_mass["surface_mm2"]) <= 1e-8 * src_mass["surface_mm2"],
            abs(src_mass["volume_mm3"] - result["adaptive_occ_mass_properties"]["reference"][0]["volume_mm3"]) <= 1e-8 * src_mass["volume_mm3"],
            abs(dst_mass["volume_mm3"] - result["adaptive_occ_mass_properties"]["candidate"][0]["volume_mm3"]) <= 1e-8 * dst_mass["volume_mm3"],
        ])
        if "analytic_volume_mm3" in result:
            checks.extend([
                abs(src_mass["volume_mm3"] - result["analytic_volume_mm3"]) <= 1e-8 * src_mass["volume_mm3"],
                abs(dst_mass["volume_mm3"] - result["analytic_volume_mm3"]) <= 1e-8 * dst_mass["volume_mm3"],
            ])
    checks.append(result["bbox_max_abs_delta_mm"] <= 1e-9)
    result["qualification_pass"] = all(checks)
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--occ-volume", type=Path)
    parser.add_argument("--fast", action="store_true", help="use one coarse classification grid")
    parser.add_argument("--source-script", type=Path)
    parser.add_argument("--candidate-recipe", type=Path)
    args = parser.parse_args()
    result = run(args.reference, args.candidate, args.occ_volume,
                 args.output.parent / args.output.stem if args.output else None,
                 args.fast, args.source_script, args.candidate_recipe)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if not result["qualification_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
