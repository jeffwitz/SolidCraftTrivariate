"""Adoption failures for the fixed two-hole bracket."""

import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from trivariate.audit import audit
from trivariate.build import build


class BracketAdoptionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bracket = build({"family":"filleted_bracket_v1"})

    def test_dimensions_and_design_faces(self):
        result = audit(self.bracket)
        self.assertEqual(result["status"],"valid")
        self.assertEqual(result["shared_faces"],55)
        self.assertEqual(len({f["feature"] for f in self.bracket["design_faces"]}),11)

    def test_removed_included_and_translated_blocks(self):
        for defect in ("missing","contained","hole_center"):
            model = copy.deepcopy(self.bracket)
            if defect=="missing":
                model["blocks"].pop()
            elif defect=="contained":
                block = next(b for b in model["blocks"] if b["id"].startswith("arm_x_50_60"))
                for p in block["controls"]:
                    p[0] -= 18*p[3]
            else:
                for block in model["blocks"]:
                    if block["id"].startswith("hole_y"):
                        for p in block["controls"]:
                            p[1] += p[3]
            with self.subTest(defect=defect):
                self.assertNotEqual(audit(model)["status"],"valid")

    def test_topology_and_selection_mutations_are_rejected(self):
        model = copy.deepcopy(self.bracket)
        model["design_faces"][0]["feature"]="obsolete_wall"
        self.assertEqual(audit(model)["reason"],"design_faces")
        model = copy.deepcopy(self.bracket)
        model["topology"]["interfaces"].pop()
        self.assertEqual(audit(model)["reason"],"incidence_mismatch")

    def test_refusal_preserves_existing_file(self):
        with tempfile.TemporaryDirectory() as temp:
            recipe = Path(temp)/"recipe.json"
            output = Path(temp)/"model.json"
            recipe.write_text(json.dumps({"family":"filleted_bracket_v1"}))
            command = [sys.executable,"-m","trivariate","build",str(recipe),str(output)]
            subprocess.run(command,check=True,capture_output=True)
            before = output.read_bytes()
            recipe.write_text(json.dumps({"family":"filleted_bracket_v1","diameter":0}))
            self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0)
            self.assertEqual(before,output.read_bytes())
