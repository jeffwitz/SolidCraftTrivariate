"""Independently read a native IGA case in FieldCraft and compare its exported field."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET


def rows(piece, path):
    node = piece.find(path)
    if node is None or node.get("format", "ascii").lower() != "ascii":
        raise ValueError(f"missing ASCII field {path}")
    width = int(node.get("NumberOfComponents", "1"))
    values = list(map(float, node.text.split()))
    if len(values) % width or not all(map(math.isfinite, values)):
        raise ValueError(f"invalid values for {path}")
    return [tuple(values[i:i + width]) for i in range(0, len(values), width)]


def samples(path, native):
    data = []
    for piece in ET.parse(path).findall(".//Piece"):
        points = rows(piece, "Points/DataArray")
        fields = [rows(piece, f"PointData/DataArray[@Name='{name}']") for name in
                  (("displacement_mm", "von_mises_mpa", "stress_mpa") if native else ("displacement", "von_mises", "stress"))]
        if native:
            fields[2] = [tuple(t[i] for i in (0, 4, 8, 1, 5, 2)) for t in fields[2]]
        data.extend(zip(points, *fields, strict=True))
    if not data:
        raise ValueError(f"empty field: {path}")
    return data


def check(cli, archive, native, out):
    out.mkdir(parents=True, exist_ok=True)
    opened = [{"command": "file.open", "params": {"path": str(archive.resolve())}}]
    reexport = out.resolve() / "fieldcraft-reexport.vtu"
    script = out / "verify.json"
    script.write_text(json.dumps(opened + [{"command": "results.export.vtu",
        "params": {"path": str(reexport), "encoding": "ascii"}}]))
    subprocess.run([cli, "run", script], check=True)
    expected, actual = samples(native, True), samples(reexport, False)
    # Node labels follow the native sampling order. Duplicated interface samples carry
    # distinct one-sided stresses and must remain separate in the visualization.
    errors = [0.0] * 4
    for lhs, rhs in zip(expected, actual, strict=True):
        for i, (a, b) in enumerate(zip(lhs, rhs, strict=True)):
            for v, w in zip(a, b, strict=True):
                errors[i] = max(errors[i], abs(v - w))
                if not math.isclose(v, w, abs_tol=2e-7, rel_tol=2e-6):
                    raise ValueError(f"field {i}: FieldCraft {w}, native {v}")
    if not any(sum(v*v for v in u) > 1e-12 and s[0] > 1e-5 for _, u, s, _ in expected):
        raise ValueError("expected a loaded mechanical case")
    script.write_text(json.dumps(opened))
    subprocess.run([cli, "snapshot", script, "--out", out/"displacement.png",
                    "--size", "840x560", "--field", "displacement"],
                   env={**os.environ, "FIELDCRAFT_GPU": "cpu"}, check=True)
    print(json.dumps({"samples": len(expected), "max_position_error_mm": errors[0],
                      "max_displacement_error_mm": errors[1], "max_von_mises_error_mpa": errors[2],
                      "max_stress_component_error_mpa": errors[3]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("native_vtu", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--fieldcraft-cli", default="fieldcraft-cli")
    args = parser.parse_args()
    check(args.fieldcraft_cli, args.archive, args.native_vtu, args.output_dir)
