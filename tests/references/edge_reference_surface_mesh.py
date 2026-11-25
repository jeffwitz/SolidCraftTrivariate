#!/usr/bin/env python3
"""Export deterministic display meshes for the six edge-treatment STEP files.

Run with ``freecadcmd tests/references/edge_reference_surface_mesh.py``. The
vertices and triangle indices come directly from FreeCAD/OCCT's shape
tessellator; coordinates are millimetres and indices are zero-based.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import FreeCAD as App
import MeshPart
import Part


ROOT = Path(__file__).resolve().parents[2]
STEP_DIR = ROOT / "docs/site/data/steps"
OUTPUT_DIR = ROOT / "docs/site/data/surfaces"
INPUTS = (
    "edges-annulus-fillet",
    "edges-annulus-chamfer",
    "edges-bracket-fillet",
    "edges-bracket-chamfer",
    "edges-S3-fillet",
    "edges-S3-chamfer",
)
LINEAR_DEFLECTION_MM = 0.08
ANGULAR_DEFLECTION_RAD = 0.18


def mesh_step(stem: str) -> tuple[int, int, int]:
    source = STEP_DIR / f"{stem}.step"
    target = OUTPUT_DIR / f"{stem}.json"
    if not source.is_file():
        raise FileNotFoundError(source)

    shape = Part.read(str(source))
    if shape.isNull() or not shape.isValid() or len(shape.Solids) != 1:
        raise ValueError(f"{source.name}: expected one valid solid")
    mesh = MeshPart.meshFromShape(
        Shape=shape,
        LinearDeflection=LINEAR_DEFLECTION_MM,
        AngularDeflection=ANGULAR_DEFLECTION_RAD,
        Relative=False,
    )
    points = mesh.Points
    facets = mesh.Facets
    vertices = [round(float(c), 9) for point in points for c in (point.x, point.y, point.z)]
    triangles = [int(index) for facet in facets for index in facet.PointIndices]
    if not vertices or not triangles or len(triangles) % 3:
        raise ValueError(f"{source.name}: empty or malformed OCCT tessellation")
    vertex_count = len(vertices) // 3
    if min(triangles) < 0 or max(triangles) >= vertex_count:
        raise ValueError(f"{source.name}: tessellation index outside vertex array")

    edge_counts: Counter[tuple[int, int]] = Counter()
    mesh_volume = 0.0
    for offset in range(0, len(triangles), 3):
        a, b, c = triangles[offset:offset + 3]
        pa, pb, pc = (vertices[3 * i:3 * i + 3] for i in (a, b, c))
        cross = (
            pb[1] * pc[2] - pb[2] * pc[1],
            pb[2] * pc[0] - pb[0] * pc[2],
            pb[0] * pc[1] - pb[1] * pc[0],
        )
        mesh_volume += sum(pa[i] * cross[i] for i in range(3)) / 6.0
        for first, second in ((a, b), (b, c), (c, a)):
            edge_counts[tuple(sorted((first, second)))] += 1
    if set(edge_counts.values()) != {2}:
        raise ValueError(f"{source.name}: tessellated surface is not a closed 2-manifold")
    relative_volume_error = abs(mesh_volume - shape.Volume) / shape.Volume
    if mesh_volume <= 0 or relative_volume_error > 0.005:
        raise ValueError(
            f"{source.name}: mesh volume {mesh_volume} differs from STEP "
            f"volume {shape.Volume} by {relative_volume_error:.3%}"
        )

    data = {"vertices": vertices, "triangles": triangles}
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(data, separators=(",", ":")) + "\n", encoding="utf-8")
    print(
        f"{target.relative_to(ROOT)}: {vertex_count} vertices, "
        f"{len(triangles)//3} triangles, volume error {relative_volume_error:.3%}"
    )
    return vertex_count, len(triangles) // 3, target.stat().st_size


print(f"FreeCAD {'.'.join(App.Version()[:3])}; OCCT tessellation")
for item in INPUTS:
    mesh_step(item)
