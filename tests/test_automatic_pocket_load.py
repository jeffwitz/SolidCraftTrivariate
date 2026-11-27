"""Physical oracle for the single-solve analysis documented in the manual."""
import json
import math
from collections import Counter
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"
JOB = ROOT / "tests/references/automatic_pocket_tension.job.json"


class AutomaticPocketLoad(unittest.TestCase):
    def test_tension_on_top_and_floor_and_atomic_refusal(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = Path(tmp) / "result.json"
            cli = ["cargo", "run", "--quiet", "--locked", "--offline", "--",
                   "solve-job", str(MODEL)]
            executed = subprocess.run([*cli, str(JOB), str(result)], cwd=ROOT,
                                      capture_output=True, text=True)
            self.assertEqual(executed.returncode, 0, executed.stderr)
            data = result.read_bytes()
            row = json.loads(data)
            # E=10 and nu=0 admit sigma_zz=1 on both upward faces. The
            # integrated constant energy density is V/(2E).
            self.assertTrue(math.isclose(row["volume_mm3"], 48976, rel_tol=1e-8))
            self.assertTrue(math.isclose(row["energy_n_mm"], 2448.8, rel_tol=1e-8))
            for key, expected in (("force_n", 2500), ("reaction_n", -2500)):
                for value, target in zip(row[key], (0, 0, expected), strict=True):
                    self.assertAlmostEqual(value, target, delta=2500e-8)
            self.assertLess(row["relative_residual"], 1e-9)
            self.assertEqual(row["continuity_selection"]["qualified_radial_join_count"], 8)
            selected = Counter((d["role"], d["field_order"])
                               for d in row["interface_decisions"])
            self.assertEqual(selected, {("common_parent_subdivision", 1): 4,
                                        ("distinct_support_join", 1): 8,
                                        ("whole_face_join", 0): 4})
            job = json.loads(JOB.read_bytes())
            job["clamp_feature"] = "missing-support"
            invalid = Path(tmp) / "invalid.json"
            invalid.write_text(json.dumps(job))
            failure = subprocess.run([*cli, str(invalid), str(result)], cwd=ROOT,
                                     capture_output=True)
            self.assertNotEqual(failure.returncode, 0)
            self.assertEqual(result.read_bytes(), data)
