#!/usr/bin/env python3
"""Solve a SolidCraft trivariate job and package its sampled IGA field for FieldCraft display.

The native SolidCraft JSON result is authoritative. The `.fieldcraft` file is a display carrier:
it contains a linear hex mesh sampled chart by chart and nodal fields sampled from that IGA result.
Chart-local samples retain their source ordering and one-sided stress values. FieldCraft does not
run a mechanics solver in this workflow.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def run(args: list[str], *, cwd: Path | None = None) -> str:
    result = subprocess.run(args, cwd=cwd, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(args)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result.stdout.strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def numbers(array: ET.Element) -> list[float]:
    return [float(v) for v in (array.text or "").split()]


def array(piece: ET.Element, section: str, name: str) -> list[float]:
    node = piece.find(f'./{section}/DataArray[@Name="{name}"]')
    if node is None:
        raise ValueError(f"VTU piece is missing {section}/{name}")
    return numbers(node)


def named_array(piece: ET.Element, name: str) -> list[float]:
    node = piece.find(f'.//DataArray[@Name="{name}"]')
    if node is None:
        raise ValueError(f"FieldCraft export is missing array {name}")
    return numbers(node)


def extract_vtu(path: Path) -> tuple[list[dict], list[list[int]]]:
    root = ET.parse(path).getroot()
    pieces = root.findall("./UnstructuredGrid/Piece")
    if not pieces:
        raise ValueError("native field VTU contains no chart pieces")
    all_nodes: list[dict] = []
    cells: list[list[int]] = []

    for piece in pieces:
        xyz = array(piece, "Points", "position_mm")
        disp = array(piece, "PointData", "displacement_mm")
        stress = array(piece, "PointData", "stress_mpa")
        mises = array(piece, "PointData", "von_mises_mpa")
        npoints = int(piece.attrib["NumberOfPoints"])
        if len(xyz) != 3 * npoints or len(disp) != 3 * npoints or len(stress) != 9 * npoints or len(mises) != npoints:
            raise ValueError("native field VTU point-array lengths do not match")
        global_offset = len(all_nodes)
        for i in range(npoints):
            point = xyz[3 * i : 3 * i + 3]
            values = disp[3 * i : 3 * i + 3] + stress[9 * i : 9 * i + 9] + [mises[i]]
            if not all(math.isfinite(v) for v in point + values):
                raise ValueError("native field VTU contains a nonfinite sample")
            # Coincident interface points remain chart-local so C0 one-sided stresses are intact.
            all_nodes.append({
                "label": global_offset + i + 1,
                "xyz_mm": point,
                "displacement_mm": disp[3 * i : 3 * i + 3],
                "stress_mpa": [stress[9 * i + j] for j in (0, 4, 8, 1, 2, 5)],
                "von_mises_mpa": mises[i],
            })

        conn = array(piece, "Cells", "connectivity")
        offsets = [int(v) for v in array(piece, "Cells", "offsets")]
        types = [int(v) for v in array(piece, "Cells", "types")]
        if any(t != 12 for t in types) or len(offsets) != len(types) or (offsets and offsets[-1] != len(conn)):
            raise ValueError("native field VTU is not linear hexahedra")
        previous = 0
        for end in offsets:
            raw = conn[previous:end]
            if len(raw) != 8 or any(v != int(v) or int(v) < 0 or int(v) >= npoints for v in raw):
                raise ValueError("native VTU has invalid hex connectivity")
            cells.append([global_offset + int(v) for v in raw])
            previous = end

    if not cells or not all_nodes:
        raise ValueError("native field VTU has no volume cells or points")
    return all_nodes, cells


def write_inp(path: Path, nodes: list[dict], cells: list[list[int]], material: dict) -> None:
    lines = [
        "*Heading",
        "SolidCraft Trivariate IGA visualization carrier; no FieldCraft solve",
        "*Node",
    ]
    for node in nodes:
        xyz = node["xyz_mm"]
        lines.append(f"{node['label']}, {xyz[0]:.17g}, {xyz[1]:.17g}, {xyz[2]:.17g}")
    lines += ["*Element, type=C3D8, elset=Solid",]
    for i, cell in enumerate(cells, start=1):
        labels = ", ".join(str(v + 1) for v in cell)
        lines.append(f"{i}, {labels}")
    lines += ["*Nset, nset=AllNodes"]
    for start in range(1, len(nodes) + 1, 16):
        lines.append(", ".join(str(v) for v in range(start, min(start + 16, len(nodes) + 1))))
    lines += [
        "*Material, name=IGA-Display-Material",
        "*Elastic",
        f"{float(material['young_mpa']):.17g}, {float(material['poisson']):.17g}",
        "*Solid Section, elset=Solid, material=IGA-Display-Material",
        ",",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def compare_roundtrip(native_path: Path, fieldcraft_path: Path) -> dict:
    native_pieces = ET.parse(native_path).getroot().findall("./UnstructuredGrid/Piece")
    display_pieces = ET.parse(fieldcraft_path).getroot().findall("./UnstructuredGrid/Piece")
    if len(display_pieces) != 1 or not native_pieces:
        raise ValueError("unexpected native or FieldCraft VTU piece count")
    native = {
        "xyz": [],
        "displacement": [],
        "stress": [],
        "von_mises": [],
    }
    for piece in native_pieces:
        native["xyz"].extend(array(piece, "Points", "position_mm"))
        native["displacement"].extend(array(piece, "PointData", "displacement_mm"))
        native["stress"].extend(array(piece, "PointData", "stress_mpa"))
        native["von_mises"].extend(array(piece, "PointData", "von_mises_mpa"))
    shown = display_pieces[0]
    points = named_array(shown, "Points")
    displacement = named_array(shown, "displacement")
    stress = named_array(shown, "stress")
    mises = named_array(shown, "von_mises")
    has_result = named_array(shown, "has_result")
    labels = [int(v) for v in named_array(shown, "node_label")]
    count = len(labels)
    if count != len(native["xyz"]) // 3 or len(has_result) != count or any(v != 1 for v in has_result):
        raise ValueError("FieldCraft result node count or validity mask differs from the native IGA samples")
    expected = {
        "xyz_mm": native["xyz"],
        "displacement_mm": native["displacement"],
        "stress_mpa": [v for i in range(count) for v in (
            native["stress"][9 * i], native["stress"][9 * i + 4], native["stress"][9 * i + 8],
            native["stress"][9 * i + 1], native["stress"][9 * i + 5], native["stress"][9 * i + 2],
        )],
        "von_mises_mpa": native["von_mises"],
    }
    actual = {"xyz_mm": points, "displacement_mm": displacement, "stress_mpa": stress, "von_mises_mpa": mises}
    scale = {"xyz_mm": 100.0, "displacement_mm": 1.0, "stress_mpa": 1000.0, "von_mises_mpa": 1000.0}
    maximum_relative = {}
    for key, values in expected.items():
        if len(values) != len(actual[key]):
            raise ValueError(f"FieldCraft roundtrip {key} length differs from native IGA")
        error = max((abs(a - b) / max(scale[key], abs(a), abs(b)) for a, b in zip(values, actual[key])), default=0.0)
        maximum_relative[key] = error
        if error > 5e-13:
            raise ValueError(f"FieldCraft roundtrip {key} differs from native IGA (scaled max error {error:.3e})")
    if labels != list(range(1, count + 1)):
        raise ValueError("FieldCraft changed the source chart-sample order")
    return {"points": count, "max_scaled_error": maximum_relative, "checks": ["positions", "displacement", "stress_tensor", "von_mises", "chart_sample_order"]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--job", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--case-name", required=True)
    parser.add_argument("--cli", type=Path, required=True, help="native solidcraft-trivariate executable")
    parser.add_argument("--fieldcraft-cli", type=Path, required=True)
    parser.add_argument("--archive-helper", type=Path, required=True, help="pinned FieldCraft v0.8.1 archive adapter")
    parser.add_argument("--subdivisions", type=int, default=4, choices=range(1, 9))
    args = parser.parse_args()
    if not args.recipe.is_file() or not args.job.is_file() or not args.cli.is_file() or not args.fieldcraft_cli.is_file() or not args.archive_helper.is_file():
        parser.error("recipe, job and all three executables must exist")
    if not args.case_name or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for c in args.case_name):
        parser.error("case-name must contain only letters, digits, hyphen and underscore")

    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    case = args.case_name
    model_path = out / f"{case}.model.json"
    result_path = out / f"{case}.result.json"
    vtu_path = out / f"{case}.field.vtu"
    archive_path = out / f"{case}.fieldcraft"
    manifest_path = out / f"{case}.manifest.json"
    work = out / f".{case}.fieldcraft-work"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir()

    job = json.loads(args.job.read_text(encoding="utf-8"))
    run([str(args.cli), "build", str(args.recipe.resolve()), str(model_path)])
    run([str(args.cli), "audit", str(model_path)])
    solved_job = work / "solve.job.json"
    solved_job.write_text(json.dumps(job, indent=2) + "\n")
    run([str(args.cli), "solve-job", str(model_path), str(solved_job), str(result_path)])
    run([str(args.cli), "export-field-vtu", str(result_path), str(args.subdivisions), str(vtu_path)])
    native = json.loads(result_path.read_text(encoding="utf-8"))
    nodes, cells = extract_vtu(vtu_path)
    sample_path = work / "sample.json"
    archive_nodes = [{
        "label": node["label"],
        "xyz_m": [v / 1000.0 for v in node["xyz_mm"]],
        "displacement_m": [v / 1000.0 for v in node["displacement_mm"]],
        "stress_pa": [v * 1e6 for v in node["stress_mpa"]],
        "von_mises_pa": node["von_mises_mpa"] * 1e6,
    } for node in nodes]
    sample_path.write_text(json.dumps({
        "schema": "solidcraft_iga_fieldcraft_sample_v1",
        "nodes": archive_nodes,
        "applied_total_n": native["applied_resultant_n"],
        "reaction_total_n": native["support_reaction_n"],
    }, separators=(",", ":")) + "\n")

    inp_path = work / "sample.inp"
    write_inp(inp_path, nodes, cells, native["job"]["material"])
    base_path = work / "base.fieldcraft"
    display_script = work / "base.commands.json"
    steps = [
        {"command": "deck.read", "params": {"format": "abaqus", "path": str(inp_path), "units": "mm-t-s-N-MPa"}},
        {"command": "model.step.create", "params": {"name": "IGA", "procedure": {"type": "static"}}},
        {"command": "file.save", "params": {"path": str(base_path), "results": False, "geometry": "embed"}},
    ]
    display_script.write_text(json.dumps(steps, indent=2) + "\n")
    run([str(args.fieldcraft_cli), "run", str(display_script)])
    run([str(args.archive_helper), str(sample_path), str(base_path), str(archive_path), "IGA"])

    render_script = work / "render.commands.json"
    render_script.write_text(json.dumps([
        {"command": "file.open", "params": {"path": str(archive_path)}},
        {"command": "results.contour", "params": {"step": "IGA", "field": "von_mises", "bands": 12}},
        {"command": "results.deformed", "params": {"step": "IGA", "show": True, "scale": "auto"}},
    ], indent=2) + "\n")
    preview_path = out / f"{case}.fieldcraft-preview.png"
    run([str(args.fieldcraft_cli), "snapshot", str(render_script), "--out", str(preview_path), "--view", "iso", "--size", "1280x800"])

    export_script = work / "export.commands.json"
    roundtrip_path = out / f"{case}.fieldcraft-export.vtu"
    export_script.write_text(json.dumps([
        {"command": "file.open", "params": {"path": str(archive_path)}},
        {"command": "results.export.vtu", "params": {"step": "IGA", "path": str(roundtrip_path), "encoding": "ascii"}},
    ], indent=2) + "\n")
    run([str(args.fieldcraft_cli), "run", str(export_script)])
    roundtrip = compare_roundtrip(vtu_path, roundtrip_path)

    manifest = {
        "schema": "solidcraft_fieldcraft_iga_display_v1",
        "case": case,
        "provenance": "Native SolidCraft Trivariate IGA solve; FieldCraft is used only to store and render the sampled result.",
        "host_solver_ran": False,
        "analysis": {
            "native_schema": native["schema"],
            "formulation": native["job"]["formulation"],
            "quadrature_order": native["job"]["order"],
            "native_relative_residual": native["relative_residual"],
            "native_strain_energy_n_mm": native["strain_energy_n_mm"],
            "applied_total_n": native["applied_resultant_n"],
            "reaction_total_n": native["support_reaction_n"],
            "nodal_reactions": "not available in the native result; omitted from FieldCraft result arrays",
            "sample_subdivisions_per_chart": args.subdivisions,
            "chart_local_visualization_nodes": len(nodes),
            "visualization_hex8_cells": len(cells),
            "fieldcraft_export_roundtrip": roundtrip,
        },
        "unit_mapping": {"geometry": "mm to m", "displacement": "mm to m", "stress": "MPa to Pa"},
        "inputs": {"recipe_sha256": sha256(args.recipe), "job_sha256": sha256(args.job), "cli_sha256": sha256(args.cli)},
        "outputs": {p.name: sha256(p) for p in (model_path, result_path, vtu_path, archive_path, preview_path, roundtrip_path)},
        "fieldcraft_cli": str(args.fieldcraft_cli.resolve()),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    shutil.rmtree(work)
    print(json.dumps({"status": "complete", "archive": str(archive_path), "native_result": str(result_path), "field_vtu": str(vtu_path), "preview": str(preview_path), "manifest": str(manifest_path)}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"fieldcraft_iga_bridge: {error}", file=sys.stderr)
        raise SystemExit(1)
