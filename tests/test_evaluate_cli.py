import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class EvaluateCliTests(unittest.TestCase):
    def run_cli(self, *arguments):
        return subprocess.run(
            [sys.executable, "-m", "trivariate", *map(str, arguments)],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

    def write_mapping(self, degrees, controls):
        temporary = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
        with temporary:
            json.dump({"degrees": degrees, "controls": controls}, temporary)
        self.addCleanup(Path(temporary.name).unlink, missing_ok=True)
        return Path(temporary.name)

    def test_evaluate_affine_map_emits_position_and_jacobian_columns(self):
        controls = [
            [float(i), float(j), float(k), 1.0]
            for i in range(2)
            for j in range(2)
            for k in range(2)
        ]
        mapping = self.write_mapping([1, 1, 1], controls)

        result = self.run_cli("evaluate", mapping, 0.2, 0.3, 0.4)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["status"], "evaluated")
        for actual, expected in zip(payload["position"], [0.2, 0.3, 0.4]):
            self.assertAlmostEqual(actual, expected, places=14)
        expected_columns = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        for actual_column, expected_column in zip(
            payload["jacobian_columns"], expected_columns
        ):
            for actual, expected in zip(actual_column, expected_column):
                self.assertAlmostEqual(actual, expected, places=14)

    def test_evaluate_rejects_inconsistent_control_count_with_json_error(self):
        mapping = self.write_mapping([1, 1, 1], [[0, 0, 0, 1]])

        result = self.run_cli("evaluate", mapping, 0.5, 0.5, 0.5)

        self.assertNotEqual(result.returncode, 0)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["status"], "error")
        self.assertTrue(payload["reason"])

    def test_evaluate_reports_unsupported_degree_family_as_json_rejection(self):
        controls = [
            [float(i), float(j), float(k), 1.0]
            for i in range(4)
            for j in range(2)
            for k in range(2)
        ]
        mapping = self.write_mapping([3, 1, 1], controls)

        result = self.run_cli("evaluate", mapping, 0.5, 0.5, 0.5)

        self.assertNotEqual(result.returncode, 0)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["status"], "unsupported")
        self.assertTrue(payload["reason"])

    def test_evaluate_rejects_nonnumeric_parameter_as_json_error(self):
        controls = [
            [float(i), float(j), float(k), 1.0]
            for i in range(2)
            for j in range(2)
            for k in range(2)
        ]
        mapping = self.write_mapping([1, 1, 1], controls)

        result = self.run_cli("evaluate", mapping, "not-a-number", 0.5, 0.5)

        self.assertNotEqual(result.returncode, 0)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["status"], "error")
        self.assertTrue(payload["reason"])

    def test_capabilities_lists_evaluator_command_and_supported_degrees(self):
        result = self.run_cli("capabilities")

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        payload = json.loads(result.stdout)
        self.assertIn("evaluate", payload["commands"])
        self.assertEqual(
            payload["evaluators"]["rational_bezier_tensor"]["degrees"],
            [[1, 1, 1], [2, 1, 1]],
        )


if __name__ == "__main__":
    unittest.main()
