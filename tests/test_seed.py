import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from trivariate.audit import audit
from trivariate.build import build


class SeedTests(unittest.TestCase):
    def test_exact_reference(self):
        recipe = json.loads(Path("tests/references/affine_pair.json").read_text())
        result = audit(build(recipe))
        self.assertEqual((result["volume_mm3"], result["min_det_mm3"]), (5, 2))

    def test_invalid_weight_orientation_and_gap(self):
        for fault in ("weight", "orientation", "gap"):
            with self.subTest(fault=fault):
                model = build({"family": "affine_pair_v1", "widths_mm": [1, 2]})
                if fault == "weight":
                    model["blocks"][0]["controls"][0][3] = 0
                elif fault == "orientation":
                    model["interface"][2] = "flip_v"
                else:
                    model["blocks"][1]["controls"][0][0] += 1
                self.assertEqual(audit(model)["status"], "invalid")

    def test_rejected_rebuild_preserves_saved_model(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "active.json"
            cli = [sys.executable, "-m", "trivariate"]
            subprocess.run(cli + ["build", "tests/references/affine_pair.json", output], check=True)
            previous = output.read_bytes()
            rejected = subprocess.run(cli + ["build", "tests/counterexamples/zero_width.json", output])
            self.assertNotEqual(rejected.returncode, 0)
            self.assertEqual(output.read_bytes(), previous)
            subprocess.run(cli + ["audit", output], check=True)
