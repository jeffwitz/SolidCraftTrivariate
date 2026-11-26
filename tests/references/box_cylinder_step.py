"""Run with freecadcmd after exporting the gallery STEP; write the CAD audit."""
import json
import math
import os
import subprocess
import tempfile
from pathlib import Path

import FreeCAD as App
import Part

directory = Path(__file__).resolve().parents[2] / "docs/site/box-cylinder"
shape = Part.read(str(directory / "result.step"))
assert len(shape.Solids) == len(shape.Shells) == 1
solid = shape.Solids[0]
assert solid.isClosed() and solid.isValid()
reference = 40 * 40 * 12 - math.pi * 5**2 * 12
with tempfile.TemporaryDirectory() as temp:
    brep = Path(temp) / "result.brep"
    solid.exportBrep(str(brep))
    rows = json.loads(subprocess.check_output([os.environ["OCC_VOLUME_BIN"], str(brep)]))
    adaptive = [{key: row[key] for key in ("tolerance", "volume_mm3", "estimated_error")} for row in rows]
assert all(abs(row["volume_mm3"] - reference) < 2e-7 for row in adaptive)
checked = 0
for z in (-.01, .01, 6, 11.99, 12.01):
    for radius in (4.99, 5.01, 9, 19, 21, 29):
        for i in range(32):
            theta = 2 * math.pi * (i + .25) / 32
            x, y = 20 + radius * math.cos(theta), 20 + radius * math.sin(theta)
            expected = 0 < x < 40 and 0 < y < 40 and 0 < z < 12 and radius > 5
            assert solid.isInside(App.Vector(x, y, z), 1e-7, False) == expected, (x, y, z)
            checked += 1
result = {"valid": True, "closed": True, "solids": 1, "shells": 1,
          "faces": len(solid.Faces), "default_volume_mm3": solid.Volume, "adaptive_mass": adaptive,
          "analytic_volume_mm3": reference, "point_classifications_matched": checked,
          "occt": Part.OCC_VERSION}
(directory / "step-audit.json").write_text(json.dumps(result, indent=2) + "\n")
