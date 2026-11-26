#!/usr/bin/env python3
"""Run one frozen phase of the external polygon partition review in FreeCADCmd."""

from __future__ import annotations

from collections import Counter
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import FreeCAD as App
import Part


ROOT = Path(__file__).resolve().parents[3]
INPUTS = ROOT / "docs/results/planar-external-review/inputs"
REFS = ROOT / "tests/references"
sys.path.insert(0, str(REFS))
import sketch_partition_oracle as oracle  # noqa: E402

METRICS_PATH = Path(os.environ["SCT_REVIEW_METRICS_SCRIPT"])
spec = importlib.util.spec_from_file_location("external_planar_quality", METRICS_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot load independent metric script: {METRICS_PATH}")
metrics = importlib.util.module_from_spec(spec)
spec.loader.exec_module(metrics)

PHASE = os.environ["SCT_REVIEW_PHASE"]
CLI = Path(os.environ["SCT_REVIEW_CLI"]).resolve()
CHECKOUT = Path(os.environ["SCT_REVIEW_CHECKOUT"]).resolve()
WORK = Path(os.environ["SCT_REVIEW_WORK"]).resolve()
OUTPUT = Path(os.environ["SCT_REVIEW_OUTPUT"]).resolve()
STRATEGY = os.environ.get("SCT_REVIEW_STRATEGY", "")
SUPPLIED_MODEL = Path(os.environ.get("SCT_REVIEW_SUPPLIED_MODEL", ""))
MANIFEST = json.loads((ROOT / "docs/results/planar-external-review/inputs.sha256.json").read_text())
INPUT_SHA = {Path(row["path"]).name: row["sha256"] for row in MANIFEST["inputs"]}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_safe(value):
    if isinstance(value, float):
        if value != value:
            return "NaN"
        if value == float("inf"):
            return "Infinity"
        if value == float("-inf"):
            return "-Infinity"
        return value
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def read_json(path: Path):
    return json.loads(path.read_text())


def run_command(args):
    started = time.perf_counter()
    proc = subprocess.run([str(arg) for arg in args], cwd=ROOT, capture_output=True, text=True)
    return proc, time.perf_counter() - started


def expected_euler(recipe):
    return (1 + sum(0 < cut["z_interval_mm"][0] and cut["z_interval_mm"][1] < recipe["height_mm"]
                    for cut in recipe["cuts"])
            - sum(cut["z_interval_mm"] == [0, recipe["height_mm"]] for cut in recipe["cuts"]))


def audit_case(name: str, source: Path):
    source_bytes = source.read_bytes()
    source_hash = hashlib.sha256(source_bytes).hexdigest()
    if source_hash != INPUT_SHA[name]:
        raise RuntimeError(f"input SHA mismatch for {name}")
    recipe = json.loads(source_bytes)
    build_recipe = json.loads(source_bytes)
    if STRATEGY:
        build_recipe["partition_strategy"] = STRATEGY

    case_work = WORK / name.removesuffix(".recipe.json")
    case_work.mkdir(parents=True, exist_ok=True)
    build_recipe_path = case_work / "build.recipe.json"
    build_recipe_path.write_text(json.dumps(build_recipe, separators=(",", ":"), allow_nan=False) + "\n")
    model_path = case_work / "model.json"
    sentinel = b"external-review-sentinel\n"
    model_path.write_bytes(sentinel)

    build_proc, build_s = run_command([CLI, "build", build_recipe_path, model_path])
    try:
        build_report = json.loads(build_proc.stdout)
    except json.JSONDecodeError:
        build_report = {"stdout": build_proc.stdout[-1000:], "stderr": build_proc.stderr[-1000:]}
    row = {
        "id": name.removesuffix(".recipe.json"),
        "source": "inputs/" + name,
        "input_sha256": source_hash,
        "build_recipe_sha256": digest(build_recipe_path),
        "strategy": STRATEGY or "automatic",
        "build_seconds": build_s,
        "build_returncode": build_proc.returncode,
        "build_report": build_report,
    }
    if build_proc.returncode:
        row["refusal_preserved_previous_target"] = model_path.read_bytes() == sentinel
        row["passed"] = False
        row["failure"] = {"stage": "build", "stderr": build_proc.stderr[-2000:], "stdout": build_proc.stdout[-2000:]}
        return row

    model = read_json(model_path)
    audit_proc, audit_s = run_command([CLI, "audit", model_path])
    try:
        audit_report = json.loads(audit_proc.stdout)
    except json.JSONDecodeError:
        audit_report = {"stdout": audit_proc.stdout[-1000:], "stderr": audit_proc.stderr[-1000:]}
    row.update({
        "block_count": len(model["blocks"]),
        "native_model_sha256": digest(model_path),
        "shared_faces": build_report.get("shared_faces"),
        "exterior_sides": build_report.get("exterior_sides"),
        "determinant_lower_bound_mm3": build_report.get("determinant_lower_bound_mm3"),
        "volume_exact_mm3": build_report.get("volume_exact_mm3"),
        "audit_seconds": audit_s,
        "audit_returncode": audit_proc.returncode,
        "audit_report": audit_report,
    })

    oracle_started = time.perf_counter()
    exact = oracle.exact_reference(recipe)
    topo = oracle.incidence(model)
    coverage = oracle.coverage(model, recipe, exact)
    provenance = oracle.verify_provenance(model, recipe)
    sample_metrics = metrics.summarize_model(model)
    step_path = case_work / "model.step"
    export_proc, export_s = run_command([CLI, "export-step", model_path, step_path])
    independent_shape = oracle.independent_cad(recipe)
    step_shape = Part.read(str(step_path)) if export_proc.returncode == 0 else None
    if step_shape is not None:
        probe_count, mismatches = oracle.cad_classification(recipe, step_shape)
        cad_rows = [independent_shape, step_shape]
        cad = []
        for shape in cad_rows:
            center = shape.Solids[0].CenterOfMass if shape.Solids else shape.CenterOfMass
            cad.append({
                "valid": bool(shape.isValid()),
                "solid_count": len(shape.Solids),
                "shell_count": len(shape.Shells),
                "volume_mm3": shape.Volume,
                "centroid_mm": [center.x, center.y, center.z],
            })
    else:
        probe_count, mismatches, cad = 0, [], []
    oracle_s = time.perf_counter() - oracle_started - export_s

    native_volume = audit_report.get("volume_exact_mm3", build_report.get("volume_exact_mm3"))
    volume_match = (native_volume is not None
                    and oracle.F(int(native_volume[0]), int(native_volume[1])) == exact["volume_mm3"])
    topology_summary = {k: v for k, v in topo.items() if k != "blocks"}
    expected_chi = expected_euler(recipe)
    topology_match = (topo["euler"] == expected_chi
                      and topo["shared_faces"] == build_report.get("shared_faces")
                      and topo["exterior_faces"] == build_report.get("exterior_sides"))
    expected_volume_float = float(exact["volume_mm3"])
    volume_relative_errors = [abs(row_['volume_mm3'] / expected_volume_float - 1.0) for row_ in cad]
    centroids = [[float(x) for x in exact["centroid_mm"]]]
    centroid_errors = []
    for cad_row in cad:
        centroid_errors.append(sum((cad_row["centroid_mm"][axis] - centroids[0][axis]) ** 2
                                   for axis in range(3)) ** 0.5)

    row.update({
        "independent_oracle_seconds": oracle_s,
        "independent_exact_volume_mm3": [exact["volume_mm3"].numerator, exact["volume_mm3"].denominator],
        "independent_centroid_mm": [[x.numerator, x.denominator] for x in exact["centroid_mm"]],
        "independent_layers": exact["layers"],
        "volume_exact_match": volume_match,
        "independent_topology": topology_summary,
        "expected_euler": expected_chi,
        "topology_match": topology_match,
        "cross_section_coverage": coverage,
        "face_provenance": provenance,
        "step_export_seconds": export_s,
        "step_sha256": digest(step_path) if step_path.is_file() else None,
        "step_export_returncode": export_proc.returncode,
        "step_export_report": json.loads(export_proc.stdout) if export_proc.stdout.strip().startswith("{") else export_proc.stdout[-500:],
        "independent_cad": cad,
        "classification_probe_count": probe_count,
        "classification_mismatches": mismatches,
        "cad_max_volume_relative_error": max(volume_relative_errors, default=float("inf")),
        "cad_max_centroid_error_mm": max(centroid_errors, default=float("inf")),
        "planar_metrics": sample_metrics,
    })
    passed = (
        audit_proc.returncode == 0
        and audit_report.get("status") == "valid"
        and volume_match
        and topo["max_face_multiplicity"] <= 2
        and topo["determinant_sample_min_mm3"] > 0
        and topology_match
        and coverage["max_symmetric_difference_mm2"] < 1e-8
        and coverage["max_pairwise_overlap_mm2"] < 1e-8
        and coverage["max_area_error_mm2"] < 1e-8
        and provenance["qualified"]
        and export_proc.returncode == 0
        and all(shape["valid"] and shape["solid_count"] == 1 for shape in cad)
        and max(volume_relative_errors, default=float("inf")) < 1e-8
        and max(centroid_errors, default=float("inf")) < 1e-7
        and not mismatches
    )
    row["passed"] = bool(passed)
    if not passed:
        row["failure"] = {
            "stage": "independent_audit",
            "audit_returncode": audit_proc.returncode,
            "audit_status": audit_report.get("status"),
            "volume_exact_match": volume_match,
            "topology_match": topology_match,
            "coverage": coverage,
            "face_provenance": provenance,
            "step_export_returncode": export_proc.returncode,
            "classification_mismatches": mismatches,
            "cad_max_volume_relative_error": max(volume_relative_errors, default=float("inf")),
            "cad_max_centroid_error_mm": max(centroid_errors, default=float("inf")),
        }
    if PHASE == "baseline" and row["id"] == "mixed_sizes":
        row["supplied_model_sha256"] = digest(SUPPLIED_MODEL)
        row["supplied_model_byte_identical"] = model_path.read_bytes() == SUPPLIED_MODEL.read_bytes()
        row["supplied_model_audit_status"] = audit_report.get("status")
    return json_safe(row)


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    cli_sha = digest(CLI)
    rows = []
    names = sorted(INPUT_SHA)
    for index, name in enumerate(names, 1):
        row = audit_case(name, INPUTS / name)
        rows.append(row)
        partial = {
            "schema": "heldout_planar_external_review_phase_v1",
            "phase": PHASE,
            "strategy": STRATEGY or "automatic",
            "checkout_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=CHECKOUT, text=True).strip(),
            "cli_sha256": cli_sha,
            "metrics_script_sha256": digest(METRICS_PATH),
            "independent_oracle_script_sha256": digest(REFS / "sketch_partition_oracle.py"),
            "input_manifest_sha256": digest(ROOT / "docs/results/planar-external-review/inputs.sha256.json"),
            "case_count": len(rows),
            "completed_count": index,
            "passed_count": sum(bool(item["passed"]) for item in rows),
            "failed_count": sum(not bool(item["passed"]) for item in rows),
            "cases": rows,
        }
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(json.dumps(json_safe(partial), indent=2, sort_keys=True, allow_nan=False) + "\n")
        print(f"{PHASE} {index}/{len(names)} {row['id']} passed={row['passed']} blocks={row.get('block_count')} build={row['build_seconds']:.3f}s audit={row.get('audit_seconds',0):.3f}s", flush=True)
    if digest(CLI) != cli_sha:
        raise RuntimeError("CLI binary changed during evaluation")
    print(f"{PHASE} complete: {sum(bool(item['passed']) for item in rows)}/{len(rows)} independently audited", flush=True)
    raise SystemExit(int(any(not item["passed"] for item in rows)))


main()
