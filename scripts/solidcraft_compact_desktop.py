"""Open a document in the isolated SolidCraft trivariate desktop profile."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import hashlib

from solidcraft_compact_profile import seed


ROOT = Path(__file__).resolve().parents[1]
REQUEST_NAME = "fieldcraft-analysis-request.json"
STATUS_NAME = "fieldcraft-analysis-status.json"


def executable(candidate):
    """Return a usable executable path from a path override or a PATH lookup."""
    if not candidate:
        return None
    path = Path(candidate).expanduser()
    resolved = shutil.which(str(path)) if not path.is_absolute() and path.parent == Path(".") else str(path)
    if not resolved:
        return None
    resolved = Path(resolved).resolve()
    return resolved if resolved.is_file() and os.access(resolved, os.X_OK) else None


def detect_fieldcraft(app_override=None, cli_override=None):
    """Detect an installed FieldCraft desktop and CLI without installing anything."""
    app_choice = app_override or os.environ.get("FIELDCRAFT_BIN")
    cli_choice = cli_override or os.environ.get("FIELDCRAFT_CLI")
    app = executable(app_choice) if app_choice else executable("fieldcraft")
    cli = executable(cli_choice) if cli_choice else executable("fieldcraft-cli")
    return {"app": app, "cli": cli, "available": bool(app and cli)}


def _write_status(companion, state):
    companion.mkdir(parents=True, exist_ok=True)
    path = companion / STATUS_NAME
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=companion,
                                     prefix="fieldcraft-status-", delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(json.dumps(state, indent=2) + "\n")
    temporary.replace(path)


def seed_analysis_inputs(companion):
    """Copy the qualified example inputs once; never replace user-edited companions."""
    defaults = (
        (ROOT / "docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.recipe.json",
         companion / "fieldcraft_x64_y64.recipe.json"),
        (ROOT / "plugins/solidcraft-trivariate/examples/fieldcraft_bracket_hex_pa5.json",
         companion / "fieldcraft_bracket_hex_pa5.json"),
    )
    for source, destination in defaults:
        if source.is_file() and not destination.exists():
            shutil.copyfile(source, destination)
    return all(destination.is_file() for _, destination in defaults)


def _safe_leaf(value, label):
    if (not isinstance(value, str) or not value or Path(value).name != value
            or value in {".", ".."} or value.startswith(".")):
        raise ValueError(f"{label} must be a simple filename in the approved companion folder")
    if any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-" for c in value):
        raise ValueError(f"{label} contains unsupported characters")
    return value


def _run_fieldcraft_request(request_path, companion, fieldcraft, trivariate_cli, archive_helper, bridge,
                           open_viewer=True):
    """Run one queued request with fixed executables and companion-confined inputs."""
    processing = request_path.with_suffix(".processing")
    try:
        request_path.replace(processing)
    except FileNotFoundError:
        return
    case_name = ""
    request = None
    recipe_hash = None
    job_hash = None
    output = None
    try:
        request = json.loads(processing.read_text(encoding="utf-8"))
        if not isinstance(request, dict) or request.get("schema") != "solidcraft_fieldcraft_analysis_request_v1":
            raise ValueError("unsupported analysis request schema")
        recipe_source = _safe_leaf(request.get("recipe_source"), "recipe_source")
        job_source = _safe_leaf(request.get("job_source"), "job_source")
        case_name = _safe_leaf(request.get("case_name"), "case_name")
        recipe_data = request.get("recipe")
        job_data = request.get("job")
        if not isinstance(recipe_data, dict) or not isinstance(job_data, dict):
            raise ValueError("request must include snapshotted recipe and job JSON objects")
        recipe_bytes = (json.dumps(recipe_data, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
        job_bytes = (json.dumps(job_data, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
        if len(recipe_bytes) > 1_048_576 or len(job_bytes) > 1_048_576:
            raise ValueError("recipe or mechanical job snapshot exceeds 1 MiB")
        results_root = companion / "results"
        results_root.mkdir(parents=True, exist_ok=True)
        output = Path(tempfile.mkdtemp(prefix=f"{case_name}-", dir=results_root))
        archive = output / f"{case_name}.fieldcraft"
        with tempfile.TemporaryDirectory(prefix="solidcraft-fieldcraft-") as request_dir:
            recipe = Path(request_dir) / "recipe.json"
            job = Path(request_dir) / "job.json"
            recipe.write_bytes(recipe_bytes)
            job.write_bytes(job_bytes)
            recipe_hash = hashlib.sha256(recipe_bytes).hexdigest()
            job_hash = hashlib.sha256(job_bytes).hexdigest()
            _write_status(companion, {
                "available": True, "state": "running", "case_name": case_name,
                "recipe_source": recipe_source, "job_source": job_source,
                "recipe_sha256": recipe_hash, "job_sha256": job_hash,
                "solver": "SolidCraft native IGA", "viewer": "FieldCraft",
                "result_carrier": "FieldCraft archive; visualization only, not a FieldCraft solve"
            })
            command = [
                sys.executable, str(bridge), "--recipe", str(recipe), "--job", str(job),
                "--out-dir", str(output), "--cli", str(trivariate_cli),
                "--fieldcraft-cli", str(fieldcraft["cli"]),
                "--archive-helper", str(archive_helper),
                "--case-name", case_name,
            ]
            subprocess.run(command, cwd=ROOT, check=True)
        if not archive.is_file():
            raise RuntimeError("bridge finished without the expected FieldCraft visualization archive")
        _write_status(companion, {
            "available": True, "state": "complete", "case_name": case_name,
            "archive": str(archive), "solver": "SolidCraft native IGA",
            "output_dir": str(output),
            "recipe_source": recipe_source, "job_source": job_source,
            "recipe_sha256": recipe_hash, "job_sha256": job_hash,
            "viewer": "FieldCraft", "result_carrier": "sampled IGA visualization; not a FieldCraft solve"
        })
        if open_viewer:
            subprocess.Popen([str(fieldcraft["app"]), str(archive)], cwd=output)
    except (OSError, ValueError, TypeError, KeyError,
            subprocess.CalledProcessError, RuntimeError) as error:
        _write_status(companion, {
            "available": True, "state": "failed", "case_name": case_name or None,
            "message": str(error), "solver": "SolidCraft native IGA",
            "recipe_source": request.get("recipe_source") if isinstance(request, dict) else None,
            "job_source": request.get("job_source") if isinstance(request, dict) else None,
            "recipe_sha256": recipe_hash, "job_sha256": job_hash,
            "output_dir": str(output) if output else None,
            "viewer": "FieldCraft", "result_carrier": "not produced"
        })
    finally:
        processing.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", type=Path, help="Installable trivariate plugin ZIP")
    parser.add_argument("document", type=Path, nargs="?",
                        default=ROOT / "plugins/solidcraft-trivariate/tests/ribbon-bracket.sc")
    parser.add_argument("--profile", type=Path, default=ROOT / "build/compact-solidcraft")
    parser.add_argument("--prepare-only", action="store_true", help="Install and verify without opening a window")
    parser.add_argument("--fieldcraft", help="FieldCraft desktop executable (auto-detected when omitted)")
    parser.add_argument("--fieldcraft-cli", help="FieldCraft CLI executable (auto-detected when omitted)")
    parser.add_argument("--trivariate-cli", help="SolidCraft Trivariate CLI used by the analysis bridge")
    parser.add_argument("--archive-helper", help="Prebuilt FieldCraft archive codec helper")
    parser.add_argument("--process-once", action="store_true",
                        help="Process one queued FieldCraft analysis without opening SolidCraft")
    parser.add_argument("--external-consumer", action="store_true",
                        help="Keep the bridge available for a separately managed consumer")
    parser.add_argument("--control", type=int, metavar="PORT",
                        help="Enable SolidCraft's localhost control interface on PORT")
    args = parser.parse_args()
    profile = args.profile.resolve()
    document = args.document.resolve(strict=True)
    package = args.package.resolve(strict=True)
    settings = profile / "settings"
    seed(settings)
    companion = profile / "companion"
    companion.mkdir(exist_ok=True)
    demo_inputs_ready = seed_analysis_inputs(companion)
    fieldcraft = detect_fieldcraft(args.fieldcraft, args.fieldcraft_cli)
    trivariate_cli = executable(args.trivariate_cli or os.environ.get("SOLIDCRAFT_TRIVARIATE_CLI"))
    if not trivariate_cli:
        for candidate in ("solidcraft-trivariate", ROOT / "target/release/solidcraft-trivariate",
                          ROOT / "target/debug/solidcraft-trivariate"):
            trivariate_cli = executable(candidate)
            if trivariate_cli:
                break
    archive_helper = executable(args.archive_helper or os.environ.get("FIELDCRAFT_IGA_ARCHIVE_HELPER"))
    if not archive_helper:
        for candidate in ("fieldcraft-iga-archive-bridge",
                          ROOT / "target/release/fieldcraft-iga-archive-bridge"):
            archive_helper = executable(candidate)
            if archive_helper:
                break
    bridge = ROOT / "scripts/fieldcraft_iga_bridge.py"
    bridge_ready = bool(fieldcraft["available"] and trivariate_cli and archive_helper
                         and bridge.is_file() and demo_inputs_ready)
    consumer_active = not args.prepare_only or args.process_once or args.external_consumer
    state = {
        "available": bridge_ready and consumer_active,
        "state": "ready" if bridge_ready and consumer_active else "unavailable",
        "message": ("FieldCraft IGA solve and result display are ready." if bridge_ready and consumer_active else
                    "FieldCraft analysis integration is unavailable; SolidCraft CAD features remain available."),
        "fieldcraft_app": str(fieldcraft["app"]) if fieldcraft["app"] else None,
        "fieldcraft_cli": str(fieldcraft["cli"]) if fieldcraft["cli"] else None,
        "solidcraft_trivariate_cli": str(trivariate_cli) if trivariate_cli else None,
        "fieldcraft_archive_helper": str(archive_helper) if archive_helper else None,
    }
    _write_status(companion, state)
    env = {**os.environ, "XDG_DATA_HOME": str(profile / "data"), "SOLIDCRAFT_SETTINGS_DIR": str(settings)}
    if args.control is not None:
        if not 1 <= args.control <= 65535:
            parser.error("--control must be a TCP port in 1..65535")
        env.setdefault("SOLIDCRAFT_CONTROL_DIR", str(profile / "control"))

    def command(name, params):
        raw = subprocess.check_output(["solidcraft-cli", "--keep-plugin-approvals", "exec",
                                       name, json.dumps(params)], env=env)
        return json.loads(raw)

    review = command("plugin.install", {"path": str(package), "dry_run": True})
    if review["id"] != "solidcraft-trivariate":
        parser.error("expected the SolidCraft Trivariate plugin package")
    grants = command("plugin.grants", {})["grants"]
    if not any(g["id"] == review["id"] and g["status"] == "current" and
               g["package_hash"] == review["approval"] and "file" in g["command_groups"] and
               any(f["name"] == "companion" and f["folder"] == str(companion) and f["read"] and f["write"]
                   for f in g["folders"]) for g in grants):
        command("plugin.install", {"path": str(package), "approval": review["approval"],
                                   "archive": review["archive"], "replace": True,
                                   "approve": {"command_groups": ["file"], "folders": [
                                       {"name": "companion", "path": str(companion), "read": True, "write": True}]}})
    inventory = command("plugin.list", {})["plugins"]
    if not any(p["manifest"]["id"] == review["id"] and not p["disabled"] for p in inventory):
        parser.error("the trivariate plugin is disabled in this profile")
    if args.process_once:
        if not bridge_ready:
            parser.error("FieldCraft bridge executables are unavailable")
        request_path = companion / REQUEST_NAME
        if not request_path.is_file():
            parser.error("no pending FieldCraft analysis request is present")
        _run_fieldcraft_request(request_path, companion, fieldcraft,
                                trivariate_cli, archive_helper, bridge, open_viewer=False)
        worker_status = json.loads((companion / STATUS_NAME).read_text(encoding="utf-8"))
        worker_status["available"] = False
        if worker_status.get("state") == "complete":
            worker_status["message"] = "One-shot analysis finished; generated FieldCraft results are available."
        else:
            worker_status["message"] = "One-shot analysis failed; inspect the recorded error and output directory."
        _write_status(companion, worker_status)
        if worker_status.get("state") != "complete":
            parser.error(worker_status["message"])
        return
    if not args.prepare_only:
        app = shutil.which("solidcraft")
        if not app:
            parser.error("could not find the SolidCraft desktop executable")
        launch = [app, "--new-instance", "--view", "model"]
        if args.control is not None:
            launch.extend(["--control", str(args.control)])
        launch.append(str(document))
        process = subprocess.Popen(launch, env=env)
        request_path = companion / REQUEST_NAME
        try:
            while process.poll() is None:
                if bridge_ready and request_path.exists():
                    _run_fieldcraft_request(request_path, companion, fieldcraft,
                                            trivariate_cli, archive_helper, bridge)
                time.sleep(0.4)
        except KeyboardInterrupt:
            process.terminate()
            process.wait()
            _write_status(companion, {
                **state, "available": False, "state": "stopped",
                "message": "The SolidCraft desktop bridge has stopped; CAD features are unchanged.",
            })
            raise
        _write_status(companion, {
            **state, "available": False, "state": "stopped",
            "message": "The SolidCraft desktop bridge has stopped; CAD features are unchanged.",
        })


if __name__ == "__main__":
    main()
