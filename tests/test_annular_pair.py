"""Independent analytic and adversarial checks for the annular pair audit."""

import copy
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from trivariate.annular import build_annular_pair
from trivariate.annular_audit import audit
from trivariate.spline import RationalBezierMap


TOL = 1e-12
ROOT = Path(__file__).resolve().parents[1]


def analytic_annular_point(u, v, z0, height, w):
    """Closed-form quarter-circle map, independent of Bernstein controls."""
    q = math.sqrt(0.5)
    a = (1.0 - u) ** 2 + 2.0 * q * u * (1.0 - u)
    b = 2.0 * q * u * (1.0 - u) + u**2
    denominator = (1.0 - u) ** 2 + 2.0 * q * u * (1.0 - u) + u**2
    radius = 2.0 - v
    return (radius * a / denominator, radius * b / denominator, z0 + height * w)


class AnnularPairAuditTests(unittest.TestCase):
    def test_valid_pair_and_analytic_point_samples(self):
        model = build_annular_pair({"family": "annular_pair_v1", "heights_mm": [1, 2]})
        result = audit(model)
        self.assertEqual(result["status"], "valid")
        self.assertEqual(result["shared_faces"], 1)
        self.assertEqual(result["height_total_mm"], 3)

        for block, z0, height in zip(model["blocks"], (0, 1), (1, 2)):
            mapping = RationalBezierMap(block["degrees"], block["controls"])
            for u, v, w in itertools.product((0.0, 0.25, 0.5, 0.75, 1.0), repeat=3):
                actual, _ = mapping.evaluate(u, v, w)
                expected = analytic_annular_point(u, v, z0, height, w)
                for a, e in zip(actual, expected):
                    self.assertAlmostEqual(a, e, delta=TOL)
                self.assertAlmostEqual(actual[0] ** 2 + actual[1] ** 2,
                                       (2.0 - v) ** 2, delta=TOL)

    def test_audit_rejects_altered_coefficients_and_incidence(self):
        original = build_annular_pair({"family": "annular_pair_v1", "heights_mm": [1, 2]})
        changed_control = copy.deepcopy(original)
        changed_control["blocks"][0]["controls"][0][0] += 1e-6
        self.assertEqual(audit(changed_control)["status"], "invalid")

        changed_trace = copy.deepcopy(original)
        changed_trace["blocks"][1]["controls"][0][2] += 1e-6
        self.assertEqual(audit(changed_trace)["reason"], "coefficient_mismatch")

        changed_incidence = copy.deepcopy(original)
        changed_incidence["interface"] = ["lower.v1", "upper.v0", "identity_uv"]
        self.assertEqual(audit(changed_incidence)["reason"], "interface_incidence")

    def test_invalid_dimensions_and_unknown_schema(self):
        for heights in ([0, 2], [1, 1001], [True, 2], [1], (1, 2)):
            model = build_annular_pair({"family": "annular_pair_v1", "heights_mm": [1, 2]})
            model["heights_mm"] = heights
            self.assertEqual(audit(model)["reason"], "invalid_heights")

        self.assertEqual(audit({"schema": "future_family"})["status"], "indeterminate")

    def test_strict_binary_coefficients_degrees_and_large_integer_refusal(self):
        original = build_annular_pair({"family": "annular_pair_v1", "heights_mm": [1, 2]})
        changed = copy.deepcopy(original)
        changed["blocks"][0]["controls"][0][0] = math.nextafter(2.0, math.inf)
        self.assertEqual(audit(changed)["reason"], "coefficient_mismatch")
        for bad_degree in (True, 2.0):
            changed = copy.deepcopy(original)
            changed["blocks"][0]["degrees"][0] = bad_degree
            self.assertEqual(audit(changed)["reason"], "degree_mismatch")
        changed = copy.deepcopy(original)
        changed["blocks"][0]["controls"][0][0] = 10**10000
        self.assertEqual(audit(changed)["reason"], "coefficient_mismatch")

    def test_ulp_regression_fixture_is_rejected(self):
        fixture = json.loads((ROOT / "tests/regressions/annular_trace_ulp.json").read_text())
        model = build_annular_pair(fixture["recipe"])
        point = model["blocks"][fixture["block"]]["controls"][fixture["control"]]
        point[fixture["component"]] = fixture["replacement"]
        self.assertEqual(audit(model)["status"], fixture["expected_status"])

    def run_cli(self, *arguments):
        return subprocess.run(
            [sys.executable, "-m", "trivariate", *map(str, arguments)],
            cwd=ROOT, capture_output=True, text=True, check=False,
        )

    def test_cli_build_roundtrip_modified_heights_and_atomic_refusal(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            recipe = directory / "recipe.json"
            output = directory / "model.json"
            recipe.write_text((ROOT / "tests/references/annular_pair.json").read_text())
            built = self.run_cli("build", recipe, output)
            self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
            self.assertEqual(json.loads(built.stdout)["status"], "valid")
            saved = json.loads(output.read_text())
            self.assertEqual(audit(saved)["status"], "valid")

            recipe.write_text(json.dumps({"family": "annular_pair_v1", "heights_mm": [4, 5]}))
            changed = self.run_cli("build", recipe, output)
            self.assertEqual(changed.returncode, 0, changed.stdout + changed.stderr)
            self.assertEqual(json.loads(output.read_text())["heights_mm"], [4, 5])

            output.write_text("preserve-me\n")
            recipe.write_text((ROOT / "tests/counterexamples/annular_zero_height.json").read_text())
            rejected = self.run_cli("build", recipe, output)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertEqual(output.read_text(), "preserve-me\n")


if __name__ == "__main__":
    unittest.main()
