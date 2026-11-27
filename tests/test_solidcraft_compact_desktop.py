import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location(
    "solidcraft_compact_desktop", ROOT / "scripts/solidcraft_compact_desktop.py"
)
desktop = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(desktop)


class FieldCraftDiscoveryTests(unittest.TestCase):
    def _executable(self, path):
        path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        path.chmod(0o755)
        return path

    def test_cad_only_when_external_application_is_absent(self):
        with tempfile.TemporaryDirectory() as root:
            empty = Path(root)
            app = empty / "no-fieldcraft"
            cli = empty / "no-fieldcraft-cli"
            found = desktop.detect_fieldcraft(str(app), str(cli))
        self.assertEqual(found, {"app": None, "cli": None, "available": False})

    def test_explicit_supported_install_is_detected(self):
        with tempfile.TemporaryDirectory() as root:
            app = self._executable(Path(root) / "fieldcraft")
            cli = self._executable(Path(root) / "fieldcraft-cli")
            found = desktop.detect_fieldcraft(str(app), str(cli))
        self.assertTrue(found["available"])
        self.assertEqual(found["app"].name, "fieldcraft")
        self.assertEqual(found["cli"].name, "fieldcraft-cli")

    def test_invalid_override_does_not_accidentally_select_path_install(self):
        with tempfile.TemporaryDirectory() as root:
            app = self._executable(Path(root) / "fieldcraft")
            cli = self._executable(Path(root) / "fieldcraft-cli")
            env = {**os.environ, "PATH": root}
            old = os.environ.copy()
            try:
                os.environ.clear()
                os.environ.update(env)
                found = desktop.detect_fieldcraft(str(Path(root) / "missing"), str(cli))
            finally:
                os.environ.clear()
                os.environ.update(old)
        self.assertIsNone(found["app"])
        self.assertTrue(found["cli"])
        self.assertFalse(found["available"])

    def test_companion_request_paths_are_leaf_names_only(self):
        self.assertEqual(desktop._safe_leaf("job-1.json", "job"), "job-1.json")
        for bad in ("../job.json", "dir/job.json", ".private", "", "job name.json"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                desktop._safe_leaf(bad, "job")

    def test_qualified_inputs_seed_once_without_replacing_local_edits(self):
        with tempfile.TemporaryDirectory() as root:
            companion = Path(root)
            self.assertTrue(desktop.seed_analysis_inputs(companion))
            job = companion / "fieldcraft_bracket_hex_pa5.json"
            job.write_text('{"user_edit":true}\n', encoding="utf-8")
            self.assertTrue(desktop.seed_analysis_inputs(companion))
            self.assertEqual(json.loads(job.read_text(encoding="utf-8")), {"user_edit": True})

    def test_worker_consumes_snapshotted_request_and_opens_archive(self):
        with tempfile.TemporaryDirectory() as root:
            base = Path(root)
            companion = base / "companion"
            companion.mkdir()
            request = companion / desktop.REQUEST_NAME
            request.write_text(json.dumps({
                "schema": "solidcraft_fieldcraft_analysis_request_v1",
                "recipe_source": "qualified.recipe.json",
                "job_source": "job.json",
                "recipe": {"family": "test_family_v1"},
                "job": {"step": "Static"},
                "case_name": "test-case",
            }), encoding="utf-8")
            bridge = base / "bridge.py"
            bridge.write_text("# fake bridge", encoding="utf-8")
            output = companion / "results"
            app = self._executable(base / "fieldcraft")
            cli = self._executable(base / "fieldcraft-cli")

            def fake_run(command, **kwargs):
                out_dir = Path(command[command.index("--out-dir") + 1])
                case_name = command[command.index("--case-name") + 1]
                out_dir.mkdir(parents=True, exist_ok=True)
                (out_dir / f"{case_name}.fieldcraft").write_bytes(b"verified archive")

            with patch.object(desktop.subprocess, "run", side_effect=fake_run) as run, \
                 patch.object(desktop.subprocess, "Popen") as launch:
                desktop._run_fieldcraft_request(
                    request, companion,
                    {"app": app, "cli": cli}, base / "solidcraft-trivariate",
                    base / "fieldcraft-iga-archive-bridge", bridge,
                )
                first = json.loads((companion / desktop.STATUS_NAME).read_text(encoding="utf-8"))
                request.write_text(json.dumps({
                    "schema": "solidcraft_fieldcraft_analysis_request_v1",
                    "recipe_source": "qualified.recipe.json",
                    "job_source": "job.json",
                    "recipe": {"family": "test_family_v1"},
                    "job": {"step": "Static"},
                    "case_name": "test-case",
                }), encoding="utf-8")
                desktop._run_fieldcraft_request(
                    request, companion,
                    {"app": app, "cli": cli}, base / "solidcraft-trivariate",
                    base / "fieldcraft-iga-archive-bridge", bridge,
                )
            self.assertFalse(request.exists())
            self.assertEqual(run.call_count, 2)
            self.assertEqual(launch.call_count, 2)
            status = json.loads((companion / desktop.STATUS_NAME).read_text(encoding="utf-8"))
            self.assertEqual(status["state"], "complete")
            self.assertEqual(status["result_carrier"], "sampled IGA visualization; not a FieldCraft solve")
            self.assertNotEqual(first["archive"], status["archive"])
            self.assertTrue(Path(first["archive"]).is_file())
            self.assertNotEqual(first["recipe_sha256"], "")
            self.assertEqual(len(status["job_sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
