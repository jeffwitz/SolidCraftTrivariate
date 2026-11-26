"""Read native FieldCraft STEP exports and probe the four manufactured bore mouths."""

import json
from math import cos, pi, sin, sqrt
import os
from pathlib import Path

import FreeCAD as App
import Part

root = Path(__file__).resolve().parents[2]
source = json.loads(
    (
        root / "docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json"
    ).read_text()
)
scale = source["transform"]["denominator"]
directory = Path(os.environ["SCT_STEP_DIR"])
states = os.environ.get("SCT_RIM_STATES", "sharp,chamfer,fillet").split(",")
report = {}
for state in states:
    shape = Part.read(str(directory / (state + ".step")))
    assert len(shape.Solids) == len(shape.Shells) == 1
    solid = shape.Solids[0]
    assert solid.isValid() and solid.isClosed()
    probes = 0
    mismatches = []
    for bore in source["source_model"]["bores"]:
        axis = "xyz".index(bore["axis"])
        center = [c / scale for c in bore["center_mm"]]
        radius = bore["radius_mm"] / scale
        for near in (True, False):
            for depth in (0.0625, 0.125, 0.1875):
                r = 0.25
                delta = {
                    "sharp": 0,
                    "chamfer": r - depth,
                    "fillet": r - sqrt(r * r - (r - depth) ** 2),
                }[state]
                for i in range(12):
                    theta = (i + 0.3) * 2 * pi / 12
                    for eps in (-0.001, 0.001):
                        rad = radius + delta + eps
                        p = center.copy()
                        p[axis] = depth if near else 8.0 - depth
                        p[(axis + 1) % 3] += rad * cos(theta)
                        p[(axis + 2) % 3] += rad * sin(theta)
                        probes += 1
                        if solid.isInside(App.Vector(*p), 1e-7, False) != (eps > 0):
                            mismatches.append(p)
    report[state] = dict(
        solids=1,
        shells=1,
        faces=len(solid.Faces),
        volume_mm3=solid.Volume,
        probes=probes,
        mismatches=mismatches,
        valid=True,
        closed=True,
    )
(directory / "step-probes.json").write_text(json.dumps(report, indent=2) + "\n")
