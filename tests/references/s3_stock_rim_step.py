"""Re-import the three native S3 rim STEP exports with FreeCADCmd."""

import json
from math import sqrt
import os
from pathlib import Path

import FreeCAD as App
import Part

directory = Path(os.environ["SCT_STEP_DIR"])
report = {}
for state in ("sharp", "chamfer", "fillet"):
    shape = Part.read(str(directory / (state + ".step")))
    assert len(shape.Solids) == len(shape.Shells) == 1
    solid = shape.Solids[0]
    assert solid.isValid() and solid.isClosed()
    checks = 0
    for depth in (0.0625, 0.125, 0.1875):
        r = 0.25
        contact = {
            "sharp": 0,
            "chamfer": r - depth,
            "fillet": r - sqrt(r * r - (r - depth) ** 2),
        }[state]
        for z in (depth, 12 - depth):
            for delta in (-0.001, 0.001):
                d = contact + delta
                for x, y in (
                    (d, 13),
                    (100 - d, 35),
                    (4, d),
                    (92, 50 - d),
                    (d, d),
                    (100 - d, 50 - d),
                ):
                    assert solid.isInside(App.Vector(x, y, z), 1e-7, False) == (
                        delta > 0
                    ), (state, x, y, z)
                    checks += 1
    report[state] = dict(
        solids=1,
        shells=1,
        faces=len(solid.Faces),
        volume_mm3=solid.Volume,
        classifications=checks,
        valid=True,
        closed=True,
    )
(directory / "step-probes.json").write_text(json.dumps(report, indent=2) + "\n")
